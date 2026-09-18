#!/bin/bash
# Generate samplesheet.csv by scanning multiple Novogene delivery folders.
#
# Handles samples sequenced across MULTIPLE LANES: writes one CSV row per
# R1/R2 lane pair, repeating the sample name. The pipeline groups rows by
# sample and merges the lanes before QC/trimming.

RAW_FOLDERS=(
        "/media/penbio24/sata4/20260810_Lucas_shotgun_metagenomes_novogene/rawdata II/X204SC26072518-Z01-F003_001/01.RawData"
)

echo "sample,fastq_1,fastq_2" > samplesheet.csv

for RAW in "${RAW_FOLDERS[@]}"; do
    for dir in "$RAW"/*/; do
        sample=$(basename "$dir")
        found_any=0
        # One R1 per lane; find its R2 mate by swapping the _1/_2 suffix
        for r1 in "$dir"*_1.fq.gz; do
            [ -e "$r1" ] || continue           # skip if the glob matched nothing
            r2="${r1%_1.fq.gz}_2.fq.gz"         # same name, _2 instead of _1
            if [ -e "$r2" ]; then
                echo "${sample},${r1},${r2}" >> samplesheet.csv
                found_any=1
            else
                echo "WARNING: no R2 mate for $r1" >&2
            fi
        done
        if [ "$found_any" -eq 0 ]; then
            echo "WARNING: no R1/R2 found for $sample in $RAW" >&2
        fi
    done
done

echo "Done. Check samplesheet.csv"