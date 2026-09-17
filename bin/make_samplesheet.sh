#!/bin/bash
# Generate samplesheet.csv by scanning multiple Novogene delivery folders.

# List every 01.RawData folder you want included.
# Add a new line here each time you verify and extract a new delivery.
RAW_FOLDERS=(
<<<<<<< HEAD
    .......
    )
=======
  ..........

)
>>>>>>> e624f53 (change input paths. rename samplesheet to create a new one with the rest samples)

# Write the CSV header
echo "sample,fastq_1,fastq_2" > samplesheet.csv

# Loop over every delivery folder, then every sample subfolder inside it
for RAW in "${RAW_FOLDERS[@]}"; do
    for dir in "$RAW"/*/; do
        sample=$(basename "$dir")
        r1=$(ls "$dir"*_1.fq.gz 2>/dev/null)
        r2=$(ls "$dir"*_2.fq.gz 2>/dev/null)
        if [ -n "$r1" ] && [ -n "$r2" ]; then
            echo "${sample},${r1},${r2}" >> samplesheet.csv
        else
            echo "WARNING: no R1/R2 found for $sample in $RAW" >&2
        fi
    done
done

echo "Done. Check samplesheet.csv"
