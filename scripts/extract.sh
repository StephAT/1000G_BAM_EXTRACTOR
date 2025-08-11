#!/bin/bash

#Directory containing Bam files
BAM_DIR="PUT_YOU_DIRECTORY_HERE/cram_files"

#output directory for extracted regions
OUTPUT_DIR="PUT_YOU_DIRECTORY_HERE/Bam"
mkdir -p "$OUTPUT_DIR"

# Specify the region to extract (e.g., "chr1:1000000-2000000")
REGION="4:143349551-145050000"

# Path to the reference genome file
REFERENCE="PUT_YOU_DIRECTORY_HEREscripts/hg19.fa"

# Loop through all BAM files in the directory
for BAM_FILE in "$BAM_DIR"/*.cram; do
    # Check if the file is actually a file
    if [ ! -f "$BAM_FILE" ]; then
        continue
    fi

    # Define the output BAM file path for the extracted region
    BASENAME=$(basename "$BAM_FILE")
    OUTPUT_BAM="$OUTPUT_DIR/${BASENAME%.bam}_extracted.bam"

    # Extract the region with reference genome specified
    echo "Extracting region $REGION from $BAM_FILE..."
    samtools view -b -T "$REFERENCE" "$BAM_FILE" "$REGION" -o "$OUTPUT_BAM"

    # Check if the extraction was successful
    if [ $? -eq 0 ]; then
        echo "Extraction completed successfully. Extracted region saved to $OUTPUT_BAM."

        # Optionally, index the extracted BAM file
        echo "Indexing $OUTPUT_BAM..."
        samtools index "$OUTPUT_BAM"
        echo "Indexing completed."

        # Delete the original BAM file to save space
        echo "Deleting original BAM file: $BAM_FILE..."
        rm "$BAM_FILE"
        echo "Deletion completed."
    else
        echo "Failed to extract region from $BAM_FILE."
    fi
done

echo "All BAM files processed."

