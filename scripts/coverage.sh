#!/bin/bash

# Set the region of interest (chromosome and positions)
REGION="4:144349551-145199151"
WINDOW_SIZE=1600

# Input and output directories
INPUT_DIR="PUT_YOU_DIRECTORY_HERE"
OUTPUT_DIR="PUT_YOU_DIRECTORY_HERE"
MAPPABILITY_FILE="PUT_YOU_DIRECTORY_HERE/crg_mappability_filtered.bed"

# Create output directory if not exists
mkdir -p "$OUTPUT_DIR"

# Create region BED file
printf "4\t144349551\t145199151\n" > "$OUTPUT_DIR/region.bed"

# Generate windows
bedtools makewindows -b "$OUTPUT_DIR/region.bed" -w $WINDOW_SIZE > "$OUTPUT_DIR/windows.bed"

# Create the combined file with the header
printf "Sample\t" > "$OUTPUT_DIR/combined_coverage.txt"
awk '{print $2}' "$OUTPUT_DIR/windows.bed" | tr '\n' '\t' | sed 's/\t$/\n/' >> "$OUTPUT_DIR/combined_coverage.txt"

# Loop through all BAM files in the input directory
for BAM in "$INPUT_DIR"/*.bam; do
    # Get the sample name
    SAMPLE=$(basename "$BAM" .bam)
    
    # Filter BAM for the region of interest
    samtools view -b "$BAM" "$REGION" > "$OUTPUT_DIR/${SAMPLE}_filtered.bam"
    
    # Calculate coverage over the windows
    bedtools coverage -a "$OUTPUT_DIR/windows.bed" -b "$OUTPUT_DIR/${SAMPLE}_filtered.bam" -mean > "$OUTPUT_DIR/${SAMPLE}_coverage.txt"
    
    # Filter by mappability
    bedtools intersect -a "$OUTPUT_DIR/${SAMPLE}_coverage.txt" -b "$MAPPABILITY_FILE" -wa -wb \
  | awk '!seen[$1,$2,$3]++' \
  > "$OUTPUT_DIR/${SAMPLE}_filtered_coverage.txt"

    #bedtools intersect -a "$OUTPUT_DIR/${SAMPLE}_coverage.txt" -b "$MAPPABILITY_FILE" -wa -wb > "$OUTPUT_DIR/${SAMPLE}_filtered_coverage.txt"
    
    # Sort files by the second column (window start position)
    sort -k3,3 "$OUTPUT_DIR/windows.bed" > "$OUTPUT_DIR/windows_sorted.bed"
    sort -k3,3 "$OUTPUT_DIR/${SAMPLE}_filtered_coverage.txt" > "$OUTPUT_DIR/${SAMPLE}_filtered_coverage_sorted.txt"
    
    # Join windows and coverage, ensuring matched order
    join -1 2 -2 2 "$OUTPUT_DIR/windows_sorted.bed" "$OUTPUT_DIR/${SAMPLE}_filtered_coverage_sorted.txt" | \
    awk '{print $6}' | tr '\n' '\t' | sed 's/\t$//' > "$OUTPUT_DIR/${SAMPLE}_matched_coverage.txt"
    
    # Add sample and matched coverage to combined file
    printf "${SAMPLE}\t" >> "$OUTPUT_DIR/combined_coverage.txt"
    cat "$OUTPUT_DIR/${SAMPLE}_matched_coverage.txt" >> "$OUTPUT_DIR/combined_coverage.txt"
    printf "\n" >> "$OUTPUT_DIR/combined_coverage.txt"
done

echo "Processing complete. Combined coverage file saved to $OUTPUT_DIR/combined_coverage.txt"
