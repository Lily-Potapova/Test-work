#!/usr/bin/env python3
import argparse
import logging
import sys
import os
from pysam import FastaFile

# Настройка логгирования
def setup_logging(log_file):
    """Настройка логгирования с выводом в файл и консоль."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),  # Запись логов в файл
            logging.StreamHandler(sys.stdout)  # Вывод логов в консоль
        ]
    )
    return logging.getLogger(__name__)

def parse_args():
    """Парсинг аргументов командной строки."""
    parser = argparse.ArgumentParser(description="Определение референсного и альтернативного аллелей для SNP.")
    parser.add_argument('-i', '--input', required=True, help="Путь к входному файлу (формат: #CHROM POS ID allele1 allele2).")
    parser.add_argument('-o', '--output', required=True, help="Путь к выходному файлу (формат: #CHROM POS ID REF ALT).")
    parser.add_argument('-r', '--ref', required=True, help="Путь к папке с референсными геномами (FASTA-файлы).")
    parser.add_argument('--log', default="/app/out/log.txt", help="Путь к файлу для записи логов.")
    return parser.parse_args()

def check_file(file_path, description):
    """Проверка наличия файла."""
    if not os.path.exists(file_path):
        logger.error(f"Файл {description} не найден: {file_path}")
        sys.exit(1)

def process_snps(input_file, output_file, ref_dir):
    """Обработка SNP и определение REF/ALT аллелей."""
    # Счётчики для статистики
    mismatched_alleles_count = 0  # Референсный аллель не совпадает с allele1 или allele2
    out_of_bounds_count = 0       # Координата выходит за пределы длины последовательности
    total_snps_processed = 0      # Общее количество обработанных SNP

    try:
        # Открываем входной и выходной файлы
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            logger.info(f"Обработка входного файла: {input_file}")
            logger.info(f"Результат будет сохранен в: {output_file}")

            # Читаем заголовок входного файла
            header = infile.readline().strip()
            if not header.startswith("#CHROM"):
                logger.error("Неправильный формат входного файла. Ожидается заголовок #CHROM POS ID allele1 allele2.")
                sys.exit(1)

            # Записываем заголовок в выходной файл
            outfile.write("#CHROM\tPOS\tID\tREF\tALT\n")

            # Обрабатываем каждую строку
            for line in infile:
                total_snps_processed += 1
                chrom, pos, snp_id, allele1, allele2 = line.strip().split()
                pos = int(pos)

                # Определяем путь к FASTA-файлу для текущей хромосомы
                if chrom.startswith("chr"):
                    ref_file = os.path.join(ref_dir, f"{chrom}.fa")
                else:
                    ref_file = os.path.join(ref_dir, f"chr{chrom}.fa")

                if not os.path.exists(ref_file):
                    logger.warning(f"Файл референсного генома для хромосомы {chrom} не найден: {ref_file}")
                    continue  # Пропустить строку, если файл не найден

                # Открываем референсный геном
                fasta = FastaFile(ref_file)
                try:
                    # Получаем длину последовательности
                    seq_length = fasta.get_reference_length(chrom)

                    # Проверяем, что координата не выходит за пределы длины последовательности
                    if pos > seq_length:
                        logger.warning(f"SNP {snp_id} ({chrom}:{pos}): Координата выходит за пределы длины последовательности ({seq_length}).")
                        out_of_bounds_count += 1
                        continue  # Пропустить строку, если координата некорректна

                    # Получаем референсный аллель
                    ref_allele = fasta.fetch(chrom, pos - 1, pos).upper()  # pysam использует 0-based координаты

                    # Определяем REF и ALT
                    if ref_allele == allele1:
                        ref, alt = allele1, allele2
                    elif ref_allele == allele2:
                        ref, alt = allele2, allele1
                    else:
                        logger.warning(f"SNP {snp_id} ({chrom}:{pos}): Референсный аллель ({ref_allele}) не совпадает с allele1 ({allele1}) или allele2 ({allele2}).")
                        mismatched_alleles_count += 1
                        continue  # Пропустить строку, если референсный аллель не совпадает

                    # Записываем результат
                    outfile.write(f"{chrom}\t{pos}\t{snp_id}\t{ref}\t{alt}\n")

                except Exception as e:
                    logger.error(f"Ошибка при обработке SNP {snp_id} ({chrom}:{pos}): {e}")
                    continue  # Пропустить строку в случае ошибки

            # Выводим статистику
            logger.info(f"Обработка завершена. Обработано SNP: {total_snps_processed}")
            logger.info(f"Количество SNP с несовпадающими аллелями: {mismatched_alleles_count}")
            logger.info(f"Количество SNP с координатами за пределами последовательности: {out_of_bounds_count}")

    except Exception as e:
        logger.error(f"Ошибка при обработке SNP: {e}")
        sys.exit(1)

def main():
    args = parse_args()

    # Создаем папку для логов, если её нет
    log_dir = os.path.dirname(args.log)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    # Настройка логгирования
    global logger
    logger = setup_logging(args.log)

    # Проверка наличия файлов
    check_file(args.input, "входной файл")
    check_file(args.ref, "папка с референсными геномами")

    # Обработка SNP
    process_snps(args.input, args.output, args.ref)

if __name__ == "__main__":
    main()