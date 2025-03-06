#!/usr/bin/env python3

import argparse
import os
import sys
import logging
from datetime import datetime
import pysam

# Настройка логирования
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

def parse_args():
    """Парсинг аргументов командной строки."""
    parser = argparse.ArgumentParser(
        description="Преобразует файл из формата '#CHROM POS ID allele1 allele2' в формат '#CHROM POS ID REF ALT'.",
        epilog="Пример использования: python script.py -i input.tsv -o output.tsv -r /path/to/reference/"
    )
    parser.add_argument(
        "-i", "--input", required=True, 
        help="Путь к входному файлу в формате '#CHROM POS ID allele1 allele2'."
    )
    parser.add_argument(
        "-o", "--output", required=True, 
        help="Путь к выходному файлу в формате '#CHROM POS ID REF ALT'."
    )
    parser.add_argument(
        "-r", "--reference", required=True, 
        help="Путь к директории с референсным геномом (chr[1-22,M,X,Y].fa)."
    )
    return parser.parse_args()

def check_header(header_line):
    """Проверяет корректность заголовка входного файла."""
    expected_columns = ["#CHROM", "POS", "ID", "allele1", "allele2"]
    actual_columns = header_line.strip().split('\t')
    if actual_columns != expected_columns:
        logger.error(f"Ошибка в заголовке файла. Ожидалось: {expected_columns}, получено: {actual_columns}")
        sys.exit(1)

def is_valid_nucleotide(allele):
    """Проверяет, является ли аллель допустимым нуклеотидом."""
    return allele.upper() in {'A', 'T', 'C', 'G', 'N'}

def process_file(input_path, output_path, reference_dir):
    """Обрабатывает входной файл и генерирует выходной."""
    logger.info(f"Начало обработки файла: {input_path}")
    
    with open(input_path, 'r', newline='') as infile, \
         open(output_path, 'w', newline='') as outfile:

        header = infile.readline()
        check_header(header)
        
        outfile.write("##fileformat=VCFv4.2\n")
        outfile.write("#CHROM\tPOS\tID\tREF\tALT\n")

        fasta_cache = {}

        for line_num, line in enumerate(infile, start=1):
            line = line.strip()
            if not line:
                continue

            try:
                chrom, pos_str, variant_id, allele1, allele2 = line.split('\t')
                
                try:
                    pos = int(pos_str)
                    if pos <= 0:
                        raise ValueError
                except ValueError:
                    logger.warning(f"Некорректная позиция в строке {line_num}: {pos_str}. Строка пропущена.")
                    continue

                if not (is_valid_nucleotide(allele1) and is_valid_nucleotide(allele2)):
                    logger.warning(f"Недопустимые аллели в строке {line_num}: {allele1}, {allele2}. Строка пропущена.")
                    continue

                if chrom not in fasta_cache:
                    fasta_path = os.path.join(reference_dir, f"chr{chrom}.fa")
                    if not os.path.exists(fasta_path):
                        logger.warning(f"Файл референсного генома для хромосомы {chrom} не найден. Строка {line_num} пропущена.")
                        continue
                    try:
                        fasta_cache[chrom] = pysam.FastaFile(fasta_path)
                    except Exception as e:
                        logger.error(f"Ошибка при открытии файла {fasta_path}: {e}")
                        del fasta_cache[chrom]
                        continue

                try:
                    ref_allele = fasta_cache[chrom].fetch(chrom, pos-1, pos).upper()
                except Exception as e:
                    logger.error(f"Ошибка при чтении позиции {chrom}:{pos}: {e}")
                    continue

                if allele1 == ref_allele:
                    ref, alt = allele1, allele2
                elif allele2 == ref_allele:
                    ref, alt = allele2, allele1
                else:
                    logger.warning(f"Ни один аллель не совпадает с REF ({ref_allele}) в строке {line_num}. Строка пропущена.")
                    continue

                outfile.write(f"{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\n")

            except Exception as e:
                logger.error(f"Ошибка обработки строки {line_num}: {e}. Строка пропущена.")
                continue

        for fa in fasta_cache.values():
            fa.close()

    logger.info(f"Обработка завершена. Результат сохранен в: {output_path}")

def main():
    args = parse_args()

    if not os.path.isfile(args.input):
        logger.error(f"Входной файл не найден: {args.input}")
        sys.exit(1)

    if not os.path.isdir(args.reference):
        logger.error(f"Директория с референсным геномом не существует: {args.reference}")
        sys.exit(1)

    process_file(args.input, args.output, args.reference)

if __name__ == "__main__":
    main()
