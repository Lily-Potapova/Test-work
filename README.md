# Test-work
Техзадание 3

# Bioinformatics Docker Image

Этот Docker-образ содержит набор специализированных программ для работы с биоинформатическими данными.

## Специализированные программы

В Docker-образе установлены следующие программы:

- **samtools** (версия 1.21) — утилиты для работы с SAM/BAM/CRAM файлами.
- **htslib** (версия 1.21) — библиотека для работы с высокопроизводительными sequencing данными.
- **bcftools** (версия 1.21) — утилиты для работы с VCF/BCF файлами.
- **vcftools** (версия 0.1.16) — инструменты для работы с VCF файлами.

## Сборка Docker-образа

Чтобы собрать Docker-образ, выполните следующую команду в терминале:

```bash
docker build -t my-bioinformatics-image -f Dockerfile.ubuntu .
```

## Запуск Docker-контейнера

Чтобы запустить контейнер в интерактивном режиме, выполните следующую команду:

```bash
docker run -it --rm my-bioinformatics-image
```
