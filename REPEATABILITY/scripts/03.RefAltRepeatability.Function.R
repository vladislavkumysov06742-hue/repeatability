# === Required libraries ===
library(Biostrings)
library(stringdist)
library(IRanges)
library(dplyr)
library(tidyr)

analyze_mtDNA_repeats <- function(fasta_path,
                                  pos,
                                  ref_nuc,
                                  alt_nuc,
                                  max_flank_left = 20,
                                  max_flank_right = 20,
                                  min_length = 5,
                                  mismatch_percent = 0.2,
                                  major_arc_start,
                                  major_arc_end,
                                  output_path = "../data/2_derived/",
                                  output_prefix = "repeatability_results") {
  # Ensure output path exists and create subfolder for pos:Ref>Alt
  if(!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
  folder_name <- paste0("pos", pos, "_", ref_nuc, "to", alt_nuc)
  full_output_dir <- file.path(output_path, folder_name)
  if (!dir.exists(full_output_dir)) dir.create(full_output_dir)
  
  # Load mtDNA sequence
  mtDNA <- readDNAStringSet(fasta_path)
  if (length(mtDNA) != 1) stop("Fasta must contain exactly one sequence.")
  seq <- mtDNA[[1]]
  
  # Helper: Extract motif sequence with asymmetric flanks
  get_motif_around_position <- function(seq, pos, left_flank, right_flank) {
    start_pos <- max(pos - left_flank, 1)
    end_pos <- min(pos + right_flank, length(seq))
    subseq(seq, start = start_pos, end = end_pos)
  }
  
  # Helper: Find approximate repeats allowing max mismatch (hamming distance)
  find_approximate_repeats <- function(seq, motif, max_mismatch) {
    motif_str <- as.character(motif)
    motif_length <- nchar(motif_str)
    seq_length <- length(seq)
    all_kmers <- substring(as.character(seq), 1:(seq_length - motif_length + 1),
                           motif_length:(seq_length))
    distances <- stringdist(motif_str, all_kmers, method = "hamming")
    matches_pos <- which(distances <= max_mismatch)
    matched_kmers <- all_kmers[matches_pos]
    data.frame(
      repeat.seq = matched_kmers,
      repeat.start = matches_pos,
      repeat.end = matches_pos + motif_length - 1,
      repeat.hamming.distance = distances[matches_pos],
      stringsAsFactors = FALSE
    )
  }
  
  motif_max_length <- max_flank_left + max_flank_right + 1
  
  # Main repeat detection function
  get_repeatability_table <- function(seq, pos, ref_nuc, alt_nuc,
                                      min_length, max_flank_left, max_flank_right, mismatch_percent) {
    all_results <- data.frame()
    for (allele in c(ref_nuc, alt_nuc)) {
      allele_label <- ifelse(allele == ref_nuc, "Ref", "Alt")
      seq_allele <- seq
      if (allele != ref_nuc) {
        seq_chars <- unlist(strsplit(as.character(seq), split = ""))
        seq_chars[pos] <- allele
        seq_allele <- DNAString(paste(seq_chars, collapse = ""))
      }
      for (motif_len in min_length:motif_max_length) {
        for (left_flank in 0:(motif_len - 1)) {
          right_flank <- motif_len - left_flank - 1
          if (left_flank <= max_flank_left && right_flank <= max_flank_right) {
            motif_seq <- get_motif_around_position(seq_allele, pos, left_flank, right_flank)
            motif_string <- as.character(motif_seq)
            this_length <- nchar(motif_string)
            max_mismatch_allowed <- floor(mismatch_percent * this_length)
            repeats_df <- find_approximate_repeats(seq_allele, motif_seq, max_mismatch_allowed)
            if (nrow(repeats_df) > 0) {
              repeats_df$motif.seq <- motif_string
              repeats_df$motif.length <- this_length
              repeats_df$motif.start <- max(pos - left_flank, 1)
              repeats_df$motif.end <- min(pos + right_flank, length(seq_allele))
              repeats_df$nuc <- allele
              repeats_df$RefAlt <- allele_label
              repeats_df$pos <- pos
              repeats_df <- repeats_df[, c(
                "pos", "nuc", "RefAlt", "motif.seq", "motif.length", "motif.start", "motif.end",
                "repeat.seq", "repeat.start", "repeat.end", "repeat.hamming.distance"
              )]
              all_results <- rbind(all_results, repeats_df)
            }
          }
        }
      }
    }
    rownames(all_results) <- NULL
    return(all_results)
  }
  
  # Filter 1: remove repeats overlapping motif region (self-overlaps)
  filter_self_overlaps <- function(df) {
    df %>% filter(motif.end < repeat.start | repeat.end < motif.start)
  }
  
  # Filter 2: remove nested repeats fully contained inside longer repeats per allele
  remove_nested_repeats <- function(repeats_df) {
    filtered_repeats <- data.frame()
    for (allele in unique(repeats_df$nuc)) {
      allele_repeats <- repeats_df %>% filter(nuc == allele)
      allele_repeats <- allele_repeats %>% arrange(desc(EffectiveLength))
      kept_intervals <- IRanges()
      for (i in seq_len(nrow(allele_repeats))) {
        current_start <- allele_repeats$repeat.start[i]
        current_end <- allele_repeats$repeat.end[i]
        current_range <- IRanges(start = current_start, end = current_end)
        is_nested <- if (length(kept_intervals) == 0) FALSE else
          any(start(current_range) >= start(kept_intervals) & end(current_range) <= end(kept_intervals))
        if (!is_nested) {
          filtered_repeats <- rbind(filtered_repeats, allele_repeats[i, ])
          kept_intervals <- c(kept_intervals, current_range)
        }
      }
    }
    rownames(filtered_repeats) <- NULL
    return(filtered_repeats)
  }
  
  # --- Execute pipeline ---
  all_repeats <- get_repeatability_table(seq, pos, ref_nuc, alt_nuc,
                                         min_length, max_flank_left, max_flank_right, mismatch_percent)
  filtered_repeats <- filter_self_overlaps(all_repeats)
  filtered_repeats$EffectiveLength <- filtered_repeats$motif.length - filtered_repeats$repeat.hamming.distance
  non_nested_repeats <- remove_nested_repeats(filtered_repeats)
  # Correct factor ordering to order Ref before Alt in sorting
  non_nested_repeats$RefAlt <- factor(non_nested_repeats$RefAlt, levels = c("Ref", "Alt"))
  
  # Sort by descending EffectiveLength, then RefAlt with Ref repeats first
  all_repeats_sorted <- non_nested_repeats %>%
    arrange(desc(EffectiveLength), RefAlt)
  
  # --- Top5 summary within major arc region ---
  major_arc_repeats <- non_nested_repeats %>%
    filter(repeat.start >= major_arc_start, repeat.end <= major_arc_end)
  
  top5_lengths <- major_arc_repeats %>%
    group_by(RefAlt, EffectiveLength) %>%
    summarise(Count = n(), .groups = "drop") %>%
    filter(EffectiveLength > 0) %>%
    arrange(desc(EffectiveLength)) %>%
    group_by(RefAlt) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    filter(rank <= 5)
  
  efflen_top <- unique(sort(top5_lengths$EffectiveLength, decreasing = TRUE))[1:5]
  
  summary_table <- tibble(EffectiveLength = efflen_top) %>%
    left_join(top5_lengths %>% filter(RefAlt == "Ref") %>%
                select(EffectiveLength, Ref_Count = Count), by = "EffectiveLength") %>%
    left_join(top5_lengths %>% filter(RefAlt == "Alt") %>%
                select(EffectiveLength, Alt_Count = Count), by = "EffectiveLength") %>%
    mutate(
      Ref_Count = ifelse(is.na(Ref_Count), 0, Ref_Count),
      Alt_Count = ifelse(is.na(Alt_Count), 0, Alt_Count)
    )
  
  # --- Write outputs ---
  write.csv(all_repeats_sorted,
            file = file.path(full_output_dir, paste0(output_prefix, "_all_repeats.csv")),
            row.names = FALSE)
  
  write.csv(summary_table,
            file = file.path(full_output_dir, paste0(output_prefix, "_major_arc_summary_top5.csv")),
            row.names = FALSE)
  
  message("Output files saved to folder: ", full_output_dir)
  
  invisible(list(
    all_repeats = all_repeats_sorted,
    summary_top5_major_arc = summary_table,
    output_folder = full_output_dir
  ))
}

# === Example usage ===
# 8473 T>C 
 analyze_mtDNA_repeats(
   fasta_path      = "../data/1_raw/Homo_sapients.mtDNA.fasta",
   pos             = 8473,
   ref_nuc         = "T",
   alt_nuc         = "C",
   max_flank_left  = 20,
   max_flank_right = 20,
   min_length      = 5,
   mismatch_percent= 0.2,
   major_arc_start = 5798,
   major_arc_end   = 16568,
   output_path     = "../data/2_derived/",
   output_prefix   = "my_mtDNA_repeat"
 )

 
#  8251 G>A
 analyze_mtDNA_repeats(
   fasta_path      = "../data/1_raw/Homo_sapients.mtDNA.fasta",
   pos             = 8251,
   ref_nuc         = "G",
   alt_nuc         = "A",
   max_flank_left  = 20,
   max_flank_right = 20,
   min_length      = 5,
   mismatch_percent= 0.2,
   major_arc_start = 5798,
   major_arc_end   = 16568,
   output_path     = "../data/2_derived/",
   output_prefix   = "my_mtDNA_repeat"
 )
 
 
 
 #  8472 C>T
 analyze_mtDNA_repeats(
   fasta_path      = "../data/1_raw/Homo_sapients.mtDNA.fasta",
   pos             = 8472,
   ref_nuc         = "C",
   alt_nuc         = "T",
   max_flank_left  = 20,
   max_flank_right = 20,
   min_length      = 5,
   mismatch_percent= 0.2,
   major_arc_start = 5798,
   major_arc_end   = 16568,
   output_path     = "../data/2_derived/",
   output_prefix   = "my_mtDNA_repeat"
 )
 
 #  12705 C>T
 analyze_mtDNA_repeats(
   fasta_path      = "../data/1_raw/Homo_sapients.mtDNA.fasta",
   pos             = 12705,
   ref_nuc         = "C",
   alt_nuc         = "T",
   max_flank_left  = 20,
   max_flank_right = 20,
   min_length      = 5,
   mismatch_percent= 0.2,
   major_arc_start = 5798,
   major_arc_end   = 16568,
   output_path     = "../data/2_derived/",
   output_prefix   = "my_mtDNA_repeat"
 )
 
#  14798 T>C
 
 analyze_mtDNA_repeats(
   fasta_path      = "../data/1_raw/Homo_sapients.mtDNA.fasta",
   pos             = 14798,
   ref_nuc         = "T",
   alt_nuc         = "C",
   max_flank_left  = 20,
   max_flank_right = 20,
   min_length      = 5,
   mismatch_percent= 0.2,
   major_arc_start = 5798,
   major_arc_end   = 16568,
   output_path     = "../data/2_derived/",
   output_prefix   = "my_mtDNA_repeat"
 )
 
 # 16223 C>T
 
 analyze_mtDNA_repeats(
   fasta_path      = "../data/1_raw/Homo_sapients.mtDNA.fasta",
   pos             = 16223,
   ref_nuc         = "C",
   alt_nuc         = "T",
   max_flank_left  = 20,
   max_flank_right = 20,
   min_length      = 5,
   mismatch_percent= 0.2,
   major_arc_start = 5798,
   major_arc_end   = 16568,
   output_path     = "../data/2_derived/",
   output_prefix   = "my_mtDNA_repeat"
 )
 
 