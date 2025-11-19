# Load required libraries at the top of your R script, not in the function
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
                                  output_prefix = "repeatability_results") {
  # 1. Read sequence
  mtDNA <- readDNAStringSet(fasta_path)
  if (length(mtDNA) != 1) stop("Fasta should contain exactly one sequence.")
  seq <- mtDNA[[1]]
  seq_len <- length(seq)
  # Ensure that the nucleotide at `pos` equals the provided reference allele
  current_nuc <- as.character(subseq(seq, start = pos, width = 1))
  if (current_nuc != ref_nuc) {
    seq_chars <- unlist(strsplit(as.character(seq), split = ""))
    seq_chars[pos] <- ref_nuc
    seq <- DNAString(paste(seq_chars, collapse = ""))
  }
  
  # Motif extraction helper
  get_motif_around_position <- function(seq, pos, left_flank, right_flank) {
    start_pos <- max(pos - left_flank, 1)
    end_pos <- min(pos + right_flank, length(seq))
    subseq(seq, start = start_pos, end = end_pos)
  }
  
  # Repeat finder helper
  find_approximate_repeats <- function(seq, motif, max_mismatch) {
    motif_str <- as.character(motif)
    motif_length <- nchar(motif_str)
    seq_length <- length(seq)
    all_kmers <- substring(as.character(seq), 1:(seq_length - motif_length + 1), motif_length:(seq_length))
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
  
  # Main repeat detection - max motif length now set by flanks
  motif_max_length <- max_flank_left + max_flank_right + 1
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
  
  # Self-overlap filter
  filter_self_overlaps <- function(df) {
    df %>% filter(motif.end < repeat.start | repeat.end < motif.start)
  }
  
  # Remove nested
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
        if (length(kept_intervals) == 0) {
          is_nested <- FALSE
        } else {
          is_nested <- any(start(current_range) >= start(kept_intervals) &
                             end(current_range) <= end(kept_intervals))
        }
        if (!is_nested) {
          filtered_repeats <- rbind(filtered_repeats, allele_repeats[i, ])
          kept_intervals <- c(kept_intervals, current_range)
        }
      }
    }
    rownames(filtered_repeats) <- NULL
    return(filtered_repeats)
  }
  
  # ---- Run Analysis ----
  all_repeats <- get_repeatability_table(seq, pos, ref_nuc, alt_nuc,
                                         min_length, max_flank_left, max_flank_right, mismatch_percent)
  filtered_repeats <- filter_self_overlaps(all_repeats)
  filtered_repeats$EffectiveLength <- filtered_repeats$motif.length - filtered_repeats$repeat.hamming.distance
  # remove nested within each allele
  non_nested_repeats <- remove_nested_repeats(filtered_repeats)
  # Explicitly ensure RefAlt column is factor
  non_nested_repeats$RefAlt <- factor(non_nested_repeats$RefAlt, levels = c("Ref","Alt"))
  
  long_repeats <- non_nested_repeats %>% filter(EffectiveLength > 10)
  # Also make sure long_repeats has the RefAlt column and it's correct
  long_repeats$RefAlt <- factor(long_repeats$RefAlt, levels = c("Ref","Alt"))
  
  # ---- Create top5 EffectiveLength summary, wide format as per your image ----
  top5_lengths <- non_nested_repeats %>%
    group_by(RefAlt, EffectiveLength) %>%
    summarise(Count = n(), .groups = "drop") %>%
    filter(EffectiveLength > 0) %>%
    arrange(desc(EffectiveLength)) %>%
    group_by(RefAlt) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    filter(rank <= 5)
  # Pivot wider: EffectiveLength in descending order, one row per EffectiveLength, columns Ref_Count, Alt_Count
  efflen_top <- unique(sort(top5_lengths$EffectiveLength, decreasing = TRUE))[1:5]
  summary_table <- tibble(EffectiveLength = efflen_top) %>%
    left_join(top5_lengths %>% filter(RefAlt == "Ref") %>% select(EffectiveLength, Ref_Count = Count), by = "EffectiveLength") %>%
    left_join(top5_lengths %>% filter(RefAlt == "Alt") %>% select(EffectiveLength, Alt_Count = Count), by = "EffectiveLength") %>%
    mutate(Ref_Count = ifelse(is.na(Ref_Count), 0, Ref_Count),
           Alt_Count = ifelse(is.na(Alt_Count), 0, Alt_Count))
  
  # ---- Export ----
  write.csv(non_nested_repeats,
            file = paste0(output_prefix, "_all_repeats.csv"),
            row.names = FALSE)
  write.csv(long_repeats,
            file = paste0(output_prefix, "_long_repeats.csv"),
            row.names = FALSE)
  write.csv(summary_table,
            file = paste0(output_prefix, "_repeat_summary_top5.csv"),
            row.names = FALSE)

  # ---- Greedy-style per-allele outputs (01KP.<pos>.<nuc>.txt) with effective.length ----
  if (nrow(non_nested_repeats) > 0) {
    non_nested_repeats$effective.length <- non_nested_repeats$motif.length - non_nested_repeats$repeat.hamming.distance
    # sort descending by effective.length
    non_nested_repeats <- non_nested_repeats[order(-non_nested_repeats$effective.length), ]
    # write one file per allele (nuc)
    alleles <- unique(non_nested_repeats$nuc)
    for (a in alleles) {
      subdf <- non_nested_repeats[non_nested_repeats$nuc == a, ]
      if (nrow(subdf) > 0) {
        outfn <- paste0("01KP.", pos, ".", a, ".txt")
        outpath <- file.path(dirname(fasta_path), outfn)
        # prepare column order similar to Greedy function
        cols <- c("pos", "nuc", "motif.seq", "motif.length", "motif.start", "motif.end",
                  "repeat.seq", "repeat.start", "repeat.end", "repeat.hamming.distance", "effective.length")
        # ensure columns exist
        cols <- cols[cols %in% colnames(subdf)]
        write.table(subdf[, cols, drop = FALSE], file = outpath, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    }
  }
  
  invisible(list(
    all_repeats = non_nested_repeats,
    long_repeats = long_repeats,
    summary_top5 = summary_table
  ))
}

##### RUN IT:


# Required packages: Biostrings, stringdist, IRanges, dplyr, tidyr
# If not installed, run:
# install.packages(c("stringdist", "dplyr", "tidyr"))
# BiocManager::install(c("Biostrings", "IRanges"))

## Example / manual run (kept as a commented block so that sourcing this
## file won't execute a full analysis). To run interactively or from the
## command line, copy the call below and adjust paths as needed.
# analyze_mtDNA_repeats(
#   fasta_path = "../data/1_raw/Homo_sapients.mtDNA.fasta",
#   pos = 8473,
#   ref_nuc = "T",
#   alt_nuc = "C",
#   max_flank_left = 20,
#   max_flank_right = 20,
#   min_length = 5,
#   mismatch_percent = 0.2,
#   output_prefix = "my_mtDNA_repeat"
# )
