#!/usr/bin/env python3
import argparse
import logging
import sys
import os
from pysam import FastaFile

# Настройка логгирования
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_args():
    """Парсинг аргументов командной строки."""
    parser = argparse.ArgumentParser(description="Определение референсного и альтернативного аллелей для SNP.")
    parser.add_argument('-i', '--input', required=True, help="Путь к входному файлу (формат: #CHROM POS ID allele1 allele2).")
    parser.add_argument('-o', '--output', required=True, help="Путь к выходному файлу (формат: #CHROM POS ID REF ALT).")
    parser.add_argument('-r', '--ref', required=True, help="Путь к папке с референсными геномами (FASTA-файлы).")
    return parser.parse_args()

def check_file(file_path, description):
    """Проверка наличия файла."""
    if not os.path.exists(file_path):
        logger.error(f"Файл {description} не найден: {file_path}")
        sys.exit(1)

def process_snps(input_file, output_file, ref_dir):
    """Обработка SNP и определение REF/ALT аллелей."""
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
                chrom, pos, snp_id, allele1, allele2 = line.strip().split()
                pos = int(pos)

                # Определяем путь к FASTA-файлу для текущей хромосомы
                if chrom.startswith("chr"):
                    ref_file = os.path.join(ref_dir, f"{chrom}.fa")
                else:
                    ref_file = os.path.join(ref_dir, f"chr{chrom}.fa")

                if not os.path.exists(ref_file):
                    logger.warning(f"Файл референсного генома для хромосомы {chrom} не найден: {ref_file}")
                    ref, alt = allele1, allele2  # Если файл не найден, оставляем как есть
                else:
                    # Открываем референсный геном
                    fasta = FastaFile(ref_file)
                    try:
                        # Получаем длину последовательности
                        seq_length = fasta.get_reference_length(chrom)

                        # Проверяем, что координата не выходит за пределы длины последовательности
                        if pos > seq_length:
                            logger.warning(f"SNP {snp_id} ({chrom}:{pos}): Координата выходит за пределы длины последовательности ({seq_length}).")
                            ref, alt = allele1, allele2  # Если координата некорректна, оставляем как есть
                        else:
                            # Получаем референсный аллель
                            ref_allele = fasta.fetch(chrom, pos - 1, pos).upper()  # pysam использует 0-based координаты

                            # Определяем REF и ALT
                            if ref_allele == allele1:
                                ref, alt = allele1, allele2
                            elif ref_allele == allele2:
                                ref, alt = allele2, allele1
                            else:
                                logger.warning(f"SNP {snp_id} ({chrom}:{pos}): Референсный аллель ({ref_allele}) не совпадает с allele1 ({allele1}) или allele2 ({allele2}).")
                                ref, alt = allele1, allele2  # Если референсный аллель не совпадает, оставляем как есть
                    except Exception as e:
                        logger.error(f"Ошибка при обработке SNP {snp_id} ({chrom}:{pos}): {e}")
                        ref, alt = allele1, allele2  # В случае ошибки оставляем как есть

                # Записываем результат
                outfile.write(f"{chrom}\t{pos}\t{snp_id}\t{ref}\t{alt}\n")

            logger.info("Обработка завершена.")

    except Exception as e:
        logger.error(f"Ошибка при обработке SNP: {e}")
        sys.exit(1)
    """Обработка SNP и определение REF/ALT аллелей."""
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
                chrom, pos, snp_id, allele1, allele2 = line.strip().split()
                pos = int(pos)

                # Определяем путь к FASTA-файлу для текущей хромосомы
                if chrom.startswith("chr"):
                    ref_file = os.path.join(ref_dir, f"{chrom}.fa")
                else:
                    ref_file = os.path.join(ref_dir, f"chr{chrom}.fa")

                if not os.path.exists(ref_file):
                    logger.warning(f"Файл референсного генома для хромосомы {chrom} не найден: {ref_file}")
                    ref, alt = allele1, allele2  # Если файл не найден, оставляем как есть
                else:
                    # Открываем референсный геном
                    fasta = FastaFile(ref_file)
                    if pos > seq_length:
                        logger.warning(f"SNP {snp_id} ({chrom}:{pos}): Координата выходит за пределы длины последовательности ({seq_length}).")
                        ref, alt = allele1, allele2  # Если координата некорректна, оставляем как есть
                    else:
                        ref_allele = fasta.fetch(chrom, pos - 1, pos).upper()  # pysam использует 0-based координаты

                    # Определяем REF и ALT
                    if ref_allele == allele1:
                        ref, alt = allele1, allele2
                    elif ref_allele == allele2:
                        ref, alt = allele2, allele1
                    else:
                        logger.warning(f"SNP {snp_id} ({chrom}:{pos}): Референсный аллель ({ref_allele}) не совпадает с allele1 ({allele1}) или allele2 ({allele2}).")
                        ref, alt = allele1, allele2  # Если референсный аллель не совпадает, оставляем как есть

                # Записываем результат
                outfile.write(f"{chrom}\t{pos}\t{snp_id}\t{ref}\t{alt}\n")

            logger.info("Обработка завершена.")

    except Exception as e:
        logger.error(f"Ошибка при обработке SNP: {e}")
        sys.exit(1)

def main():
    args = parse_args()

    # Проверка наличия файлов
    check_file(args.input, "входной файл")
    check_file(args.ref, "папка с референсными геномами")

    # Обработка SNP
    process_snps(args.input, args.output, args.ref)

if __name__ == "__main__":
    main()