import argparse
from pathlib import Path
from .reference import read_reference
from .clean import clean_file
from .plot import run_plot
from .stats import run_stats, run_comparison
from .pipeline import run_all
from concurrent.futures import ProcessPoolExecutor, as_completed

def main():
    parser = argparse.ArgumentParser(description='Анализ повторяемости митохондриальной ДНК')
    parser.add_argument('mode', choices=['clean','plot','stats','all','compare'])
    parser.add_argument('--input_dir', type=str, help='Директория с входными файлами (для clean/plot/stats/all)')
    parser.add_argument('--output_dir', type=str, help='Директория для результатов')
    parser.add_argument('--workers', type=int, default=8, help='Количество параллельных процессов')
    parser.add_argument('--reference', type=str, help='Путь к референсному геному в формате FASTA')
    parser.add_argument('--window_size', type=int, default=0, help='Размер окна для оконного анализа (0 = отключено)')
    parser.add_argument('--snv_file', type=str, help='CSV файл с однонуклеотидными заменами (колонки: position, ref_allele, alt_allele)')
    
    # Аргументы для сравнения двух геномов (режим compare)
    parser.add_argument('--input_dir1', type=str, help='Директория с очищенными файлами первого генома')
    parser.add_argument('--reference1', type=str, help='Референс первого генома')
    parser.add_argument('--input_dir2', type=str, help='Директория с очищенными файлами второго генома')
    parser.add_argument('--reference2', type=str, help='Референс второго генома')
    parser.add_argument('--name1', type=str, default='Genome1', help='Название первого генома для подписей на графиках')
    parser.add_argument('--name2', type=str, default='Genome2', help='Название второго генома для подписей на графиках')
    parser.add_argument('--compare_only', action='store_true', help='Только сравнительные графики (без пересчёта индивидуальных)')
    
    args = parser.parse_args()

    if args.mode == 'clean':
        input_dir = Path(args.input_dir or "01KP")
        output_dir = Path(args.output_dir or input_dir / "cleaned")
        if not args.reference:
            print("Ошибка: для режима clean необходимо указать --reference")
            return
        ref_seq = read_reference(args.reference)
        MT_LEN = len(ref_seq)
        print(f"Режим очистки: {input_dir} -> {output_dir}")
        txt_files = list(input_dir.glob("*.txt"))
        if not txt_files:
            print("Нет файлов с расширением .txt.")
            return
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futures = {ex.submit(clean_file, f, output_dir, MT_LEN): f for f in txt_files}
            for future in as_completed(futures):
                f = futures[future]
                try:
                    future.result()
                except Exception as e:
                    print(f"   Ошибка при обработке {f.name}: {e}")
        print("Очистка завершена.")

    elif args.mode == 'plot':
        input_dir = Path(args.input_dir or "mut")
        output_dir = Path(args.output_dir or input_dir / "plots")
        if not args.reference:
            print("Ошибка: для режима plot необходимо указать --reference")
            return
        ref_seq = read_reference(args.reference)
        MT_LEN = len(ref_seq)
        run_plot(input_dir, output_dir, MT_LEN)

    elif args.mode == 'stats':
        if not args.input_dir:
            print("Ошибка: укажите --input_dir")
            return
        input_dir = Path(args.input_dir)
        output_dir = Path(args.output_dir or input_dir / "stats")
        run_stats(
            input_dir=input_dir,
            output_dir=output_dir,
            ref_path=args.reference,
            workers=args.workers,
            window_size=args.window_size,
            snv_file=args.snv_file,
            input_dir2=None,
            ref_path2=None,
            name1=args.name1,
            name2=args.name2
        )

    elif args.mode == 'all':
        if not args.input_dir:
            print("Ошибка: укажите --input_dir")
            return
        input_dir = Path(args.input_dir)
        output_dir = Path(args.output_dir or input_dir / "results")
        run_all(
            input_dir=input_dir,
            output_dir=output_dir,
            ref_path=args.reference,
            workers=args.workers,
            window_size=args.window_size,
            snv_file=args.snv_file
        )

    elif args.mode == 'compare':
        if not (args.input_dir1 and args.reference1 and args.input_dir2 and args.reference2):
            print("Ошибка: для режима compare укажите --input_dir1, --reference1, --input_dir2, --reference2")
            return
        output_dir = Path(args.output_dir or "comparison_results")
        run_comparison(
            input_dir1=args.input_dir1,
            ref_path1=args.reference1,
            input_dir2=args.input_dir2,
            ref_path2=args.reference2,
            output_dir=output_dir,
            workers=args.workers,
            window_size=args.window_size,
            snv_file=args.snv_file,
            name1=args.name1,
            name2=args.name2,
            compare_only=args.compare_only   # передаём флаг
        )

if __name__ == "__main__":
    main()
