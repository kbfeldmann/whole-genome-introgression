#!/bin/bash

# MODIFIED FROM ERIK R FUNK
# Merge multiple bam files from same individual using the short ID prefix from
# previous steps. This step uses an underscore to separate wildcard
# Output file will then be re-indexed

samples="/data5/K_Feldmann_data/sample_names_merge.txt"
outdir="/data5/K_Feldmann_data/sorted_bam_files"

while read -r ID; do
echo "merging: "$ID
samtools merge $outdir/"$ID"_sorted_RGadded_dupmarked_merged.bam $outdir/"$ID"_L1_*.bam $outdir/"$ID"_L2_*.bam
samtools index $outdir/"$ID"_sorted_RGadded_dupmarked_merged.bam

done<"$samples"
