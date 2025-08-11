# Genomics Coverage Analysis Project

This project provides a complete pipeline for downloading, processing, and analyzing genomic data from the 1000 Genomes Project, specifically focusing on coverage analysis across African populations.

## Overview

The pipeline downloads CRAM files from the 1000 Genomes Project Phase 3 data, extracts specific genomic regions, and calculates coverage statistics across defined windows. It's designed to analyze coverage patterns in African populations with mappability filtering.

## Project Structure

```
project/
├── Download.py              # Downloads CRAM files from 1000 Genomes
├── scripts/
│   ├── extract.sh          # Extracts specific genomic regions from CRAM files
│   ├── index.sh           # Creates indices for BAM/CRAM files
│   ├── coverage.sh        # Calculates coverage statistics
│   └── fixed_windows.bed  # Predefined genomic windows
└── README.md
```

## Prerequisites

### Software Requirements
- Python 3.x with libraries:
  - `requests`
  - `beautifulsoup4`
  - `tqdm`
- SAMtools
- BEDtools
- Access to a reference genome file (hg19.fa)

### Hardware Requirements
- Sufficient storage space (CRAM files can be several GB each)
- Internet connection for downloading
- Adequate computational resources for parallel processing

## Installation

1. Clone or download this project
2. Install Python dependencies:
   ```bash
   pip install requests beautifulsoup4 tqdm
   ```
3. Install SAMtools and BEDtools through your package manager or from source
4. Download the hg19 reference genome

## Configuration

Before running the scripts, update the following paths in each script:

### In `scripts/extract.sh`:
```bash
BAM_DIR="PUT_YOU_DIRECTORY_HERE/cram_files"
OUTPUT_DIR="PUT_YOU_DIRECTORY_HERE/Bam"
REFERENCE="PUT_YOU_DIRECTORY_HERE/scripts/hg19.fa"
```

### In `scripts/coverage.sh`:
```bash
INPUT_DIR="PUT_YOU_DIRECTORY_HERE"
OUTPUT_DIR="PUT_YOU_DIRECTORY_HERE"
MAPPABILITY_FILE="PUT_YOU_DIRECTORY_HERE/crg_mappability_filtered.bed"
```

### In `scripts/index.sh`:
```bash
CRAM_DIR="PUT_YOU_DIRECTORY_HERE"
```

## Usage

### Step 1: Download CRAM Files

The download script fetches CRAM files from specific African populations:

```bash
python Download.py
```

**Populations included:**
- LWK, CHA, PAR, WAS, ACB, ASW, FUL, JOL, MAN
- WOF, GWD, MSL, MOS, ESN, YRI, BAN, SBA

**Features:**
- Downloads up to 60 files per population
- Supports resume functionality for interrupted downloads
- Automatically indexes CRAM files using samtools
- Uses 3 concurrent threads for faster processing

### Step 2: Extract Genomic Regions

Extract specific regions from the downloaded CRAM files:

```bash
bash scripts/extract.sh
```

**Default region:** chromosome 4:143349551-145050000

This script:
- Converts CRAM to BAM format for the specified region
- Creates indices for the extracted BAM files
- Removes original CRAM files to save space

### Step 3: Index BAM Files (if needed)

If you need to create additional indices:

```bash
bash scripts/index.sh
```

### Step 4: Calculate Coverage Statistics

Analyze coverage across predefined windows:

```bash
bash scripts/coverage.sh
```

**Default analysis:**
- **Region:** chromosome 4:144349551-145199151
- **Window size:** 1,600 bp
- **Output:** Combined coverage matrix for all samples

This script:
- Creates genomic windows using BEDtools
- Calculates mean coverage per window
- Filters by mappability regions
- Generates a combined coverage file with all samples

## Output Files

### From Download.py:
- `cram_files/`: Directory containing downloaded CRAM files and indices

### From extract.sh:
- `Bam/`: Directory containing extracted BAM files for the specified region
- `*.bam.bai`: Index files for each extracted BAM file

### From coverage.sh:
- `combined_coverage.txt`: Matrix with samples as rows and window positions as columns
- `*_coverage.txt`: Individual coverage files per sample
- `*_filtered_coverage.txt`: Mappability-filtered coverage files
- `windows.bed`: BED file defining the analysis windows

## Key Features

### Population Focus
The project specifically targets African populations from the 1000 Genomes Project, making it suitable for population genetics studies in these groups.

### Mappability Filtering
Coverage analysis includes mappability filtering to exclude regions where reads cannot be uniquely mapped, improving analysis quality.

### Scalable Processing
- Parallel downloading with concurrent threads
- Batch processing of multiple samples
- Resume capability for interrupted downloads

### Quality Control
- File integrity checks during download
- Automatic indexing of genomic files
- Error handling and logging

## Customization

### Changing Target Regions
Modify the `REGION` variable in `extract.sh` and `coverage.sh` to analyze different genomic regions.

### Adjusting Window Size
Change the `WINDOW_SIZE` variable in `coverage.sh` to use different window sizes for coverage analysis.

### Population Selection
Modify the `POPULATION_CODES` list in `Download.py` to include different populations.

### Sample Limits
Adjust the sample limit (currently 60 per population) in the `Download.py` script as needed.

## Troubleshooting

### Common Issues

1. **Download failures:** Check internet connection and 1000 Genomes FTP availability
2. **Permission errors:** Ensure write permissions in output directories
3. **Memory issues:** Reduce concurrent threads in `Download.py`
4. **Missing reference:** Ensure hg19.fa is properly downloaded and indexed

### Log Files
Check the console output for detailed error messages and progress information.

## Data Sources

- **1000 Genomes Project Phase 3:** ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/
- **Reference Genome:** hg19 (GRCh37)
- **Mappability tracks:** CRG mappability (user-provided)

## Citation

If you use this pipeline in your research, please cite:
- The 1000 Genomes Project Consortium
- SAMtools and BEDtools publications
- Any relevant population genetics studies

## License

This project is provided as-is for research purposes. Please ensure compliance with 1000 Genomes Project data usage policies.

## Support

For issues or questions, please check:
1. Software documentation (SAMtools, BEDtools)
2. 1000 Genomes Project documentation
3. Error logs and console output
