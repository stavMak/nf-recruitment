#!/bin/bash
# Generate samplesheet.csv by scanning the Novogene sample folders.

# Folder that directly contains the sample subfolders (LUC_338, LUC_339, ...)
RAW="/media/penbio24/sata4/20260810_Lucas_shotgun_metagenomes_novogene/rawdata I/X204SC26072518-Z01-F002_01/01.RawData"

# Write the CSV header
echo "sample,fastq_1,fastq_2" > samplesheet.csv

# Loop over each sample folder inside $RAW
for dir in "$RAW"/*/; do
    sample=$(basename "$dir")               # folder name, e.g. LUC_338
    r1=$(ls "$dir"*_1.fq.gz 2>/dev/null)    # the R1 file
    r2=$(ls "$dir"*_2.fq.gz 2>/dev/null)    # the R2 file
    if [ -n "$r1" ] && [ -n "$r2" ]; then
        echo "${sample},${r1},${r2}" >> samplesheet.csv
    else
        echo "WARNING: no R1/R2 found for $sample" >&2
    fi
done

echo "Done. Check samplesheet.csv"
