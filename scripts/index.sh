#!/bin/bash

# Set the directory containing CRAM files
CRAM_DIR="PUT_YOU_DIRECTORY_HERE"  # Replace with the directory containing your CRAM files

# Loop through all .cram files in the directory
for cram_file in "${CRAM_DIR}"/*.bam; do
  # Check if the file exists
  if [[ -f "${cram_file}" ]]; then
    echo "Indexing ${cram_file}..."
    
    # Use samtools to create the index
    samtools index "${cram_file}"
    
    # Check if the index (.crai) was successfully created
    if [[ -f "${cram_file}.bai" ]]; then
      echo "Successfully indexed ${cram_file}"
    else
      echo "Failed to index ${cram_file}" >&2
    fi
  else
    echo "No CRAM files found in the directory." >&2
    break
  fi
done

echo "All CRAM files processed."
