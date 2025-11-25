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
    # === Required libraries ===
    # Подключаем библиотеки, которые используются в функции ниже.
    # Biostrings: структуры для работы с ДНК-последовательностями (DNAString, DNAStringSet)
    # stringdist: вычисление расстояний между строками (в т.ч. Hamming)
    # IRanges: удобный тип для хранения и сравнения геномных интервалов
    # dplyr/tidyr: манипуляции с таблицами (data.frame / tibble)
    library(Biostrings)
    library(stringdist)
    library(IRanges)
    library(dplyr)
    library(tidyr)

    ## Функция: analyze_mtDNA_repeats
    ## Описание: запускает весь пайплайн поиска повторов вокруг заданной позиции
    ## Параметры:
    ## - fasta_path: путь к FASTA-файлу с mtDNA (ожидается одна последовательность)
    ## - pos: позиция варианта в геноме (1-based)
    ## - ref_nuc, alt_nuc: буквы нуклеотидов (например, "T", "C")
    ## - max_flank_left/right: максимальный размер левого/правого фланка при формировании мотивов
    ## - min_length: минимальная длина мотивов для поиска
    ## - mismatch_percent: доля допустимых несовпадений (например 0.2 = 20%)
    ## - major_arc_start/major_arc_end: границы 'major arc' (фильтрация итогов)
    ## - output_path: директория, в которую будут записаны результаты (создается, если нет)
    ## - output_prefix: префикс для имён выходных файлов
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
      # --- Настройка директорий ---
      # Если родительская output-папка не существует, создаём её рекурсивно
      if(!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
      # Создаём подпапку, уникальную для позиции и пары Ref->Alt
      folder_name <- paste0("pos", pos, "_", ref_nuc, "to", alt_nuc)
      full_output_dir <- file.path(output_path, folder_name)
      if (!dir.exists(full_output_dir)) dir.create(full_output_dir)
  
      # --- Чтение последовательности ---
      # readDNAStringSet возвращает объект DNAStringSet (может содержать несколько последовательностей)
      mtDNA <- readDNAStringSet(fasta_path)
      # Проверяем, что в FASTA ровно одна запись — иначе алгоритм ожидал единичную последовательность
      if (length(mtDNA) != 1) stop("Fasta must contain exactly one sequence.")
      seq <- mtDNA[[1]]  # извлекаем сам объект DNAString
  
      # --- Вспомогательная функция: извлечение мотива вокруг позиции ---
      # seq: объект DNAString
      # pos: позиция (1-based)
      # left_flank/right_flank: сколько нуклеотидов взять влево/вправо
      get_motif_around_position <- function(seq, pos, left_flank, right_flank) {
        # Учитываем границы последовательности: start не меньше 1, end не больше длины
        start_pos <- max(pos - left_flank, 1)
        end_pos <- min(pos + right_flank, length(seq))
        # subseq возвращает объект DNAString — подстроку исходной последовательности
        subseq(seq, start = start_pos, end = end_pos)
      }
  
      # --- Вспомогательная функция: поиск приближённых повторов ---
      # seq: объект DNAString (вся последовательность)
      # motif: объект DNAString с искомым мотивом
      # max_mismatch: максимально допустимый Hamming distance (целое число)
      find_approximate_repeats <- function(seq, motif, max_mismatch) {
        # Преобразуем motif в строку и узнаём его длину
        motif_str <- as.character(motif)
        motif_length <- nchar(motif_str)
        seq_length <- length(seq) # длина полной последовательности
        # Получаем все k-mers указанной длины в последовательности.
        # substring(str, starts, ends) возвращает вектор строк — все кандидаты для сравнения
        all_kmers <- substring(as.character(seq), 1:(seq_length - motif_length + 1),
                               motif_length:(seq_length))
        # Вычисляем расстояния Хэмминга от motif_str до каждой k-mer
        distances <- stringdist(motif_str, all_kmers, method = "hamming")
        # Выбираем позиции, где расстояние <= max_mismatch
        matches_pos <- which(distances <= max_mismatch)
        matched_kmers <- all_kmers[matches_pos]
        # Возвращаем data.frame с информацией о каждом найденном повторе
        data.frame(
          repeat.seq = matched_kmers,
          repeat.start = matches_pos,
          repeat.end = matches_pos + motif_length - 1,
          repeat.hamming.distance = distances[matches_pos],
          stringsAsFactors = FALSE
        )
      }
  
      # Максимальная длина мотива, которую мы будем рассматривать — определяется фланками
      motif_max_length <- max_flank_left + max_flank_right + 1
  
      # --- Основная функция: перебор всех мотивов и сбор таблицы повторов ---
      # Проходит по двум аллелям: Ref и Alt; для каждого аллеля строит модифицированную последовательность
      # и ищет все подходящие повторы для всех допустимых длин мотивов и фланков
      get_repeatability_table <- function(seq, pos, ref_nuc, alt_nuc,
                                          min_length, max_flank_left, max_flank_right, mismatch_percent) {
        all_results <- data.frame()  # сюда будем накапливать найденные повторы
        # Проходим по двум нуклеотидам: референсному и альтернативному
        for (allele in c(ref_nuc, alt_nuc)) {
          allele_label <- ifelse(allele == ref_nuc, "Ref", "Alt")
          seq_allele <- seq
          # Если рассматриваемая аллель не равна референсу — заменяем нуклеотид в позиции
          if (allele != ref_nuc) {
            seq_chars <- unlist(strsplit(as.character(seq), split = ""))
            seq_chars[pos] <- allele
            seq_allele <- DNAString(paste(seq_chars, collapse = ""))
          }
          # Для каждой длины мотива (min_length..motif_max_length) перебираем
          for (motif_len in min_length:motif_max_length) {
            # left_flank принимает значения от 0 до motif_len-1, а right_flank вычисляется как остаток
            for (left_flank in 0:(motif_len - 1)) {
              right_flank <- motif_len - left_flank - 1
              # Убеждаемся, что фланки не превышают заданных лимитов
              if (left_flank <= max_flank_left && right_flank <= max_flank_right) {
                # Берём мотив (DNAString) с ассиметричными фланками
                motif_seq <- get_motif_around_position(seq_allele, pos, left_flank, right_flank)
                motif_string <- as.character(motif_seq)
                this_length <- nchar(motif_string)
                # Максимально допустимое количество различий для данного мотива
                max_mismatch_allowed <- floor(mismatch_percent * this_length)
                # Ищем подходящие повторы в последовательности с учётом максимального числа несовпадений
                repeats_df <- find_approximate_repeats(seq_allele, motif_seq, max_mismatch_allowed)
                if (nrow(repeats_df) > 0) {
                  # Добавляем метаинформацию о мотиве и аллели
                  repeats_df$motif.seq <- motif_string
                  repeats_df$motif.length <- this_length
                  repeats_df$motif.start <- max(pos - left_flank, 1)
                  repeats_df$motif.end <- min(pos + right_flank, length(seq_allele))
                  repeats_df$nuc <- allele
                  repeats_df$RefAlt <- allele_label
                  repeats_df$pos <- pos
                  # Упорядочиваем колонки для удобства
                  repeats_df <- repeats_df[, c(
                    "pos", "nuc", "RefAlt", "motif.seq", "motif.length", "motif.start", "motif.end",
                    "repeat.seq", "repeat.start", "repeat.end", "repeat.hamming.distance"
                  )]
                  # Склеиваем найденные строки к общему результату
                  all_results <- rbind(all_results, repeats_df)
                }
              }
            }
          }
        }
        rownames(all_results) <- NULL
        return(all_results)
      }
  
      # --- Фильтрация self-overlap ---
      # Убираем те строки, где мотив (область вокруг позиции) и найденный повтор пересекаются
      # (т.е. мотив сам находит себя — такие совпадения не интересуют для нашей метрики)
      filter_self_overlaps <- function(df) {
        df %>% filter(motif.end < repeat.start | repeat.end < motif.start)
      }
  
      # --- Удаление вложенных (nested) повторов ---
      # Для каждой аллели оставляем только непересекающиеся интервалы, начиная с самых длинных
      remove_nested_repeats <- function(repeats_df) {
        filtered_repeats <- data.frame()
        for (allele in unique(repeats_df$nuc)) {
          allele_repeats <- repeats_df %>% filter(nuc == allele)
          # Сортируем по EffectiveLength (будет добавлено позже), от большего к меньшему
          allele_repeats <- allele_repeats %>% arrange(desc(EffectiveLength))
          kept_intervals <- IRanges()
          for (i in seq_len(nrow(allele_repeats))) {
            current_start <- allele_repeats$repeat.start[i]
            current_end <- allele_repeats$repeat.end[i]
            current_range <- IRanges(start = current_start, end = current_end)
            # is_nested = TRUE, если текущий интервал полностью находится внутри одного из уже
            # сохранённых (kept_intervals). В этом случае мы его пропустим.
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
  
      # ---- Выполнение пайплайна анализа ----
      # 1) Собираем все найденные повторы для Ref и Alt
      all_repeats <- get_repeatability_table(seq, pos, ref_nuc, alt_nuc,
                                             min_length, max_flank_left, max_flank_right, mismatch_percent)
      # 2) Убираем self-overlaps
      filtered_repeats <- filter_self_overlaps(all_repeats)
      # 3) Вычисляем EffectiveLength = длина мотива - Hamming distance (приближённая "эффективная" длина совпадения)
      filtered_repeats$EffectiveLength <- filtered_repeats$motif.length - filtered_repeats$repeat.hamming.distance
      # 4) Удаляем nested-повторы
      non_nested_repeats <- remove_nested_repeats(filtered_repeats)
      # Гарантируем корректный порядок уровней фактора для RefAlt (Ref до Alt)
      non_nested_repeats$RefAlt <- factor(non_nested_repeats$RefAlt, levels = c("Ref", "Alt"))
  
      # Сортируем итоговую таблицу: сначала по убыванию EffectiveLength, затем Ref/Alt
      all_repeats_sorted <- non_nested_repeats %>%
        arrange(desc(EffectiveLength), RefAlt)
  
      # --- Формируем Top-5 summary внутри major arc ---
      # Ограничиваем повторы по координатам major_arc_start..major_arc_end
      major_arc_repeats <- non_nested_repeats %>%
        filter(repeat.start >= major_arc_start, repeat.end <= major_arc_end)
  
      # Группируем по RefAlt и EffectiveLength и считаем количество
      top5_lengths <- major_arc_repeats %>%
        group_by(RefAlt, EffectiveLength) %>%
        summarise(Count = n(), .groups = "drop") %>%
        filter(EffectiveLength > 0) %>%
        arrange(desc(EffectiveLength)) %>%
        group_by(RefAlt) %>%
        mutate(rank = row_number()) %>%
        ungroup() %>%
        filter(rank <= 5)
  
      # Берём 5 уникальных наиболее длинных EffectiveLength (если их меньше — NA появится и далее заменится 0)
      efflen_top <- unique(sort(top5_lengths$EffectiveLength, decreasing = TRUE))[1:5]
  
      # Собираем итоговую таблицу с колонками Ref_Count и Alt_Count для тех EffectiveLength
      summary_table <- tibble(EffectiveLength = efflen_top) %>%
        left_join(top5_lengths %>% filter(RefAlt == "Ref") %>%
                    select(EffectiveLength, Ref_Count = Count), by = "EffectiveLength") %>%
        left_join(top5_lengths %>% filter(RefAlt == "Alt") %>%
                    select(EffectiveLength, Alt_Count = Count), by = "EffectiveLength") %>%
        mutate(
          Ref_Count = ifelse(is.na(Ref_Count), 0, Ref_Count),
          Alt_Count = ifelse(is.na(Alt_Count), 0, Alt_Count)
        )
  
      # --- Запись результатов на диск ---
      write.csv(all_repeats_sorted,
                file = file.path(full_output_dir, paste0(output_prefix, "_all_repeats.csv")),
                row.names = FALSE)
  
      write.csv(summary_table,
                file = file.path(full_output_dir, paste0(output_prefix, "_major_arc_summary_top5.csv")),
                row.names = FALSE)
  
      message("Output files saved to folder: ", full_output_dir)
  
      # Возвращаем список с основными результатами; invisible чтобы результат не печатался автоматически при source()
      invisible(list(
        all_repeats = all_repeats_sorted,
        summary_top5_major_arc = summary_table,
        output_folder = full_output_dir
      ))
    }

    # === Примеры использования ===
    # Ниже приведены примеры вызова функции для нескольких позиций. Комментарии поясняют аргументы.
    # Если хотите, можно раскомментировать любой вызов и выполнить скрипт напрямую через Rscript.
    # Пример: 8473 T>C
    ## analyze_mtDNA_repeats(
    ##   fasta_path      = "../data/1_raw/Homo_sapients.mtDNA.fasta",
    ##   pos             = 8473,
    ##   ref_nuc         = "T",
    ##   alt_nuc         = "C",
    ##   max_flank_left  = 20,
    ##   max_flank_right = 20,
    ##   min_length      = 5,
    ##   mismatch_percent= 0.2,
    ##   major_arc_start = 5798,
    ##   major_arc_end   = 16568,
    ##   output_path     = "../data/2_derived/",
    ##   output_prefix   = "my_mtDNA_repeat"
    ## )

