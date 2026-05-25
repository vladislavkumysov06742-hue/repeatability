"""
KP function to Get All Repeats for All Motivs Overlapping The Position with a Given Nucleotide
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Tuple, Optional
from numba import jit, njit
import warnings

warnings.filterwarnings("ignore")


@njit
def hamming_distance_numba(s1: str, s2: str) -> int:
    """Вычисляет расстояние Хэмминга между двумя строками (оптимизировано с Numba)"""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def hamming_distance_vectorized(seq_str: str, motif: str) -> np.ndarray:
    """Векторизованное вычисление расстояний Хэмминга для всех k-mer"""
    motif_len = len(motif)
    seq_len = len(seq_str)

    if motif_len > seq_len:
        return np.array([])

    kmers = np.array(
        [seq_str[i : i + motif_len] for i in range(seq_len - motif_len + 1)]
    )

    motif_arr = np.frombuffer(motif.encode(), dtype="uint8")

    distances = np.zeros(len(kmers), dtype=int)
    for i, kmer in enumerate(kmers):
        kmer_arr = np.frombuffer(kmer.encode(), dtype="uint8")
        distances[i] = np.sum(motif_arr != kmer_arr)

    return distances


def find_approximate_repeats_vectorized(
    seq_str: str, motif: str, max_mismatch: int
) -> pd.DataFrame:
    """Векторизованный поиск приблизительных повторов"""
    motif_len = len(motif)
    seq_len = len(seq_str)

    if motif_len > seq_len:
        return pd.DataFrame()

    starts = np.arange(seq_len - motif_len + 1)

    distances = np.zeros(len(starts), dtype=int)
    motif_arr = np.frombuffer(motif.encode(), dtype="uint8")

    for i, start in enumerate(starts):
        kmer = seq_str[start : start + motif_len]
        kmer_arr = np.frombuffer(kmer.encode(), dtype="uint8")
        distances[i] = np.sum(motif_arr != kmer_arr)

    mask = distances <= max_mismatch
    filtered_starts = starts[mask]
    filtered_distances = distances[mask]

    if len(filtered_starts) == 0:
        return pd.DataFrame()

    results = []
    for start, dist in zip(filtered_starts, filtered_distances):
        results.append(
            {
                "repeat.seq": seq_str[start : start + motif_len],
                "repeat.start": start + 1,  # 1-based
                "repeat.end": start + motif_len,  # 1-based
                "repeat.hamming.distance": dist,
            }
        )

    return pd.DataFrame(results)


def retrieve_all_repeats_covering_given_position_and_nucleotide_in_sequence(
    seq: str,
    pos: int,
    fixed_nuc: str,
    output_folder: str,
    max_flank_left: int = 20,
    max_flank_right: int = 20,
    min_length: int = 5,
    max_length: int = 41,
    mismatch_fraction: float = 0.2,
) -> pd.DataFrame:
    """
    Основная функция для поиска всех повторов (оптимизированная)
    """
    seq_list = list(seq)
    if seq_list[pos - 1] != fixed_nuc:
        seq_list[pos - 1] = fixed_nuc
    seq_modified = "".join(seq_list)

    seq_len = len(seq_modified)

    def get_motif_around_position(
        position: int, left_flank: int, right_flank: int
    ) -> Tuple[str, int, int]:
        """Извлекает мотив вокруг позиции с заданными флангами"""
        start_pos = max(position - left_flank - 1, 0)
        end_pos = min(position + right_flank, seq_len)
        motif = seq_modified[start_pos:end_pos]
        return motif, start_pos + 1, end_pos

    all_results = []

    for motif_length in range(min_length, max_length + 1):
        max_start = seq_len - motif_length + 1

        for left_flank in range(motif_length):
            right_flank = motif_length - left_flank - 1

            if left_flank <= max_flank_left and right_flank <= max_flank_right:
                motif_seq, motif_start, motif_end = get_motif_around_position(
                    pos, left_flank, right_flank
                )

                this_length = len(motif_seq)
                if this_length == 0:
                    continue

                max_mismatch_allowed = int(mismatch_fraction * this_length)

                repeats_df = find_approximate_repeats_vectorized(
                    seq_modified, motif_seq, max_mismatch_allowed
                )

                if len(repeats_df) > 0:
                    repeats_df["motif.seq"] = motif_seq
                    repeats_df["motif.length"] = this_length
                    repeats_df["motif.start"] = motif_start
                    repeats_df["motif.end"] = motif_end
                    repeats_df["nuc"] = fixed_nuc
                    repeats_df["effective.length"] = (
                        repeats_df["motif.length"]
                        - repeats_df["repeat.hamming.distance"]
                    )
                    repeats_df["pos"] = pos

                    all_results.append(repeats_df)

    if all_results:
        results_df = pd.concat(all_results, ignore_index=True)

        results_df = results_df.sort_values("effective.length", ascending=False)

        column_order = [
            "pos",
            "nuc",
            "motif.seq",
            "motif.length",
            "motif.start",
            "motif.end",
            "repeat.seq",
            "repeat.start",
            "repeat.end",
            "repeat.hamming.distance",
            "effective.length",
        ]
        results_df = results_df[column_order]

        output_filename = f"01KP.{pos}.{fixed_nuc}.txt"
        output_path = Path(output_folder) / output_filename
        output_path.parent.mkdir(parents=True, exist_ok=True)

        results_df.to_csv(output_path, sep="\t", index=False)
        print(f"Сохранено: {output_path}, записей: {len(results_df)}")
    else:
        print("Повторы не найдены")
        results_df = pd.DataFrame()

    return results_df


class RepeatFinderOptimized:
    """Класс для оптимизированного поиска повторов с кешированием"""

    def __init__(self, seq: str):
        self.seq = seq
        self.seq_len = len(seq)
        self._kmers_cache = {}

    def get_kmers_for_length(self, length: int) -> np.ndarray:
        """Получает все k-mer заданной длины (с кешированием)"""
        if length not in self._kmers_cache:
            kmers = np.array(
                [self.seq[i : i + length] for i in range(self.seq_len - length + 1)]
            )
            self._kmers_cache[length] = kmers
        return self._kmers_cache[length]

    def find_repeats_fast(self, pos: int, fixed_nuc: str, **kwargs) -> pd.DataFrame:
        """Быстрый поиск повторов с использованием кеширования"""
        seq_list = list(self.seq)
        if seq_list[pos - 1] != fixed_nuc:
            seq_list[pos - 1] = fixed_nuc
        seq_modified = "".join(seq_list)

        finder = RepeatFinderOptimized(seq_modified)

        return finder._internal_find(pos, fixed_nuc, **kwargs)

    def _internal_find(self, pos: int, fixed_nuc: str, **kwargs) -> pd.DataFrame:
        """Внутренний метод поиска"""
        pass
