import requests
from bs4 import BeautifulSoup
import concurrent.futures
from tqdm import tqdm
import os
import subprocess
import logging
from collections import defaultdict

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Base URL
base_url = "https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/"

# List of population codes to filter (preserve order)
POPULATION_CODES = [
    "LWK", "CHA", "PAR", "WAS", "ACB", "ASW", "FUL", "JOL", "MAN", 
    "WOF", "GWD", "MSL", "MOS", "ESN", "YRI", "BAN", "SBA"
]

def get_subdirectories(url):
    """Get all subdirectories from the URL in reverse order"""
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    return sorted([a['href'] for a in soup.find_all('a', href=True) if a['href'].endswith('/')], reverse=True)

def get_cram_files(url):
    """Get all CRAM files from a directory in reverse order"""
    response = requests.get(url)
    if response.status_code != 200:
        return []
    soup = BeautifulSoup(response.text, 'html.parser')
    return sorted([a['href'] for a in soup.find_all('a', href=True) if a['href'].endswith('.cram')], reverse=True)

def extract_population(filename):
    """Extract the first matching population code from the filename"""
    for code in POPULATION_CODES:
        if code in filename:
            return code
    return None

def check_file_integrity(file_path, expected_size):
    """Check if file exists and is completely downloaded"""
    if not os.path.exists(file_path):
        return False, 0
    actual_size = os.path.getsize(file_path)
    return actual_size == expected_size, actual_size

def resume_download(url, file_path, existing_size, total_size):
    """Resume downloading a partially downloaded file"""
    headers = {'Range': f'bytes={existing_size}-'}
    response = requests.get(url, stream=True, headers=headers)
    progress_bar = tqdm(total=total_size, initial=existing_size, unit='iB', unit_scale=True, desc=f"Resuming {os.path.basename(file_path)}")
    with open(file_path, 'ab') as f:
        for chunk in response.iter_content(chunk_size=1024):
            if chunk:
                progress_bar.update(len(chunk))
                f.write(chunk)
    progress_bar.close()
    return os.path.getsize(file_path) == total_size

def index_cram_file(cram_file):
    """Index a CRAM file using samtools"""
    index_file = f"{cram_file}.crai"
    if os.path.exists(index_file):
        logging.info(f"Index already exists for {cram_file}, skipping download.")
        return True
    try:
        logging.info(f"Indexing {cram_file}")
        cmd = ['samtools', 'index', cram_file]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            logging.error(f"Failed to index {cram_file}: {result.stderr.decode()}")
            return False
        logging.info(f"Successfully indexed {cram_file}")
        return True
    except Exception as e:
        logging.error(f"Error indexing {cram_file}: {str(e)}")
        return False

def process_cram_file(cram_url, output_dir):
    """Download and index a CRAM file"""
    cram_file = os.path.join(output_dir, os.path.basename(cram_url))
    index_file = f"{cram_file}.crai"
    if os.path.exists(index_file):
        logging.info(f"Skipping {cram_file}, index file already exists.")
        return True
    response = requests.head(cram_url)
    total_size = int(response.headers.get('content-length', 0))
    file_exists, existing_size = check_file_integrity(cram_file, total_size)
    if file_exists and existing_size == total_size:
        logging.info(f"File {cram_file} already exists and is complete.")
    elif existing_size > 0 and existing_size < total_size:
        logging.info(f"Resuming download of {cram_file}")
        if not resume_download(cram_url, cram_file, existing_size, total_size):
            logging.error(f"Failed to resume download of {cram_file}")
            return None
    else:
        logging.info(f"Starting new download of {cram_file}")
        response = requests.get(cram_url, stream=True)
        progress_bar = tqdm(total=total_size, unit='iB', unit_scale=True, desc=os.path.basename(cram_file))
        with open(cram_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=1024):
                if chunk:
                    progress_bar.update(len(chunk))
                    f.write(chunk)
        progress_bar.close()
        if os.path.getsize(cram_file) != total_size:
            logging.error(f"Download verification failed for {cram_file}")
            return None
    return index_cram_file(cram_file)

def main():
    output_dir = "cram_files"
    os.makedirs(output_dir, exist_ok=True)
    
    # Collect all CRAM URLs in reverse order
    all_cram_urls = []
    subdirs = get_subdirectories(base_url)
    
    for subdir in subdirs:
        alignment_url = f"{base_url}{subdir}alignment/"
        cram_files = get_cram_files(alignment_url)
        for cram_file in cram_files:
            cram_url = f"{alignment_url}{cram_file}"
            all_cram_urls.append(cram_url)
    
    # Group URLs by population
    population_urls = defaultdict(list)
    for url in all_cram_urls:
        population = extract_population(os.path.basename(url))
        if population:
            population_urls[population].append(url)
    
    # Select up to 60 URLs per population (starting from the bottom)
    selected_urls = []
    for code in POPULATION_CODES:
        urls = population_urls.get(code, [])[:60]
        selected_urls.extend(urls)
    
    # Process selected URLs
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        futures = []
        for url in selected_urls:
            futures.append(executor.submit(process_cram_file, url, output_dir))
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                logging.error(f"An error occurred: {str(e)}")

if __name__ == "__main__":
    main()
