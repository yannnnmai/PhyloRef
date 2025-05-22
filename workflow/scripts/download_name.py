#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    This script downloads mitochondrial genomes for a large list of species names.
    It automatically splits the list into chunks (default=1000 species per chunk)
    to mitigate potential issues with very large inputs. Each chunk is processed
    with multi-threading and includes a timeout+retry mechanism.

    Results are saved in:
        <output_directory>/outputs/
            ├─ final_download.log          (unified log file)
            ├─ final_download_summary.csv  (unified CSV summary, sorted by Listed Species Name)
            └─ gb_files/                  (all .gb files)

    The final_download.log ends with a Summary section in the following format:

    ----------------------------------------
    Summary:

    Total species in input list: X
    Total species with sequences found: Y (Total sequences downloaded: Z)

    Species skipped due to duplication: W
      <duplicated species list>

    Species with no sequences found: N
      <no-sequence species list>
    ----------------------------------------

Usage Example:
    python download_mtgenomes_name.py \
        -i /path/to/large_species_list.txt \
        -o /path/to/desired_output \
        -e your_email@example.com \
        -k your_ncbi_api_key \
        -c animals \
        -s 1000 \ # Chunk size (default: 1000)
        -m 5 \ # Maximum records to fetch for each species (default: 5)
        -t 8 # Number of threads (recommended value; not intended to be changed)

"""

import os
import sys
import re
import csv
import time
import signal
import argparse
import datetime
import platform
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
USE_ALARM = platform.system() != "Windows" 
try:
    from Bio import Entrez
except ImportError:
    sys.exit("Biopython is required. Please install it via: pip install biopython")

DEFAULT_TIMEOUT = 3600     
DEFAULT_MAX_RETRIES = 3    
DEFAULT_CHUNK_SIZE = 1000  

ACCESSION_RE = re.compile(
    r"^(?:[A-Z]{1,2}_[0-9]{6}|[A-Z]{1,2}[0-9]{5,6})(?:\.[0-9]+)?$", re.I
)

def is_accession(token: str) -> bool:
    token = token.strip()
    return bool(ACCESSION_RE.match(token))

class TimeoutException(Exception):
    pass

def handler(signum, frame):
    raise TimeoutException()

if USE_ALARM:
    signal.signal(signal.SIGALRM, handler)   

def robust_fetch(db, id, rettype, retmode):
    attempts = 0
    while attempts < 3:
        try:
            with Entrez.efetch(db=db, id=id, rettype=rettype, retmode=retmode) as fetch_handle:
                data = fetch_handle.read()
            return data
        except Exception as e:
            print(f"Failed to fetch {id} (attempt {attempts+1}/3): {str(e)}", flush=True)
            attempts += 1
            time.sleep(5 * attempts)  # exponential backoff
    return None

def process_sequence(seq_id, output_dir):
    seq_record = robust_fetch("nucleotide", seq_id, "gb", "text")
    if seq_record is None:
        return None, None, None

    ncbi_species_name = None
    LOCUS_number = None
    for line in seq_record.split('\n'):
        if line.startswith("  ORGANISM"):
            ncbi_species_name = ' '.join(line.split()[1:])
        if line.startswith("LOCUS"):
            LOCUS_number = line.split()[1]

    if LOCUS_number:
        safe_species_name = ncbi_species_name.replace(' ', '_') if ncbi_species_name else "Unknown"
        file_path = os.path.join(output_dir, f"{safe_species_name}|{LOCUS_number}.gb")
        with open(file_path, "w") as f:
            f.write(seq_record)
        return ncbi_species_name, LOCUS_number, file_path
    return None, None, None


def download_single_accession(acc_id,
                              output_dir,
                              log_fp,
                              csv_data,
                              ncbi_processed,
                              skipped_species):
    species_name, version, _ = process_sequence(acc_id, output_dir)
    if species_name:
        if species_name not in ncbi_processed:
            ncbi_processed[species_name] = f"Accession:{acc_id}"
        else:
            skipped_species.append(f"{acc_id}(dup-species:{species_name})")

        csv_data.append([
            acc_id,               
            species_name,
            1,                    
            version,
            "Accession provided"
        ])
        return True
    else:
        csv_data.append([
            acc_id, 'NA', 'NA', acc_id, 'Failed to download accession'
        ])
        return False
    
def download_mt_genomes_for_chunk(category,
                                  species_names,
                                  email,
                                  api_key,
                                  output_dir,
                                  log_fp,
                                  max_records=10,
                                  max_threads=8,
                                  ncbi_processed=None,
                                  skipped_species=None,
                                  no_sequence_species=None,
                                  csv_data=None):
    Entrez.email = email
    Entrez.api_key = api_key

    if ncbi_processed is None:
        ncbi_processed = {}
    if skipped_species is None:
        skipped_species = []
    if no_sequence_species is None:
        no_sequence_species = []
    if csv_data is None:
        csv_data = []

    for listed_species in species_names:
        search_query = f"{listed_species}[orgn] AND {category}[filter] AND mitochondrion[filter] AND complete genome[title]"
        try:
            search_handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=max_records)
            search_results = Entrez.read(search_handle)
            search_handle.close()
        except Exception as e:
            log_fp.write(f"ERROR searching {listed_species}: {e}\n")
            no_sequence_species.append(listed_species)
            csv_data.append([listed_species, 'NA', 'NA', 'NA', 'Search failed'])
            continue

        if not search_results['IdList']:
            no_sequence_species.append(listed_species)
            csv_data.append([listed_species, 'NA', 'NA', 'NA', 'No sequences found'])
            continue

        downloaded_data = []
        with ThreadPoolExecutor(max_workers=max_threads) as executor:
            future_to_seq_id = {
                executor.submit(process_sequence, seq_id, output_dir): seq_id
                for seq_id in search_results['IdList']
            }
            for future in as_completed(future_to_seq_id):
                species_name, version, file_path = future.result()
                if species_name:
                    downloaded_data.append((species_name, version, file_path))

        if downloaded_data:
            ncbi_species_name, version, _ = downloaded_data[0]
            accession_list = [data[1] for data in downloaded_data]
            if ncbi_species_name not in ncbi_processed:
                ncbi_processed[ncbi_species_name] = listed_species
                note = ("Consistent name" if listed_species == ncbi_species_name
                        else f"Synonym for {ncbi_species_name}")
                csv_data.append([
                    listed_species,
                    ncbi_species_name,
                    len(downloaded_data),
                    ", ".join(accession_list),
                    note
                ])
            else:
                skipped_species.append(listed_species)
                csv_data.append([
                    listed_species,
                    ncbi_species_name,
                    '0',
                    'NA',
                    'Duplicate, skipped download'
                ])
        else:
            no_sequence_species.append(listed_species)
            csv_data.append([listed_species, 'NA', 'NA', 'NA', 'No valid record downloaded'])

    return ncbi_processed, skipped_species, no_sequence_species, csv_data

def run_chunk_with_timeout_and_retry(species_chunk,
                                     category,
                                     email,
                                     api_key,
                                     output_dir,
                                     log_fp,
                                     max_records,
                                     max_threads,
                                     ncbi_processed,
                                     skipped_species,
                                     no_sequence_species,
                                     csv_data):
    remaining_species = species_chunk.copy()
    for attempt in range(1, DEFAULT_MAX_RETRIES + 2):  
        if not remaining_species:
            break

        log_fp.write(f"[RETRY] Chunk retry attempt {attempt}, {len(remaining_species)} species remaining.\n")
        log_fp.flush()
        print(f"[INFO] Attempt {attempt} for current chunk with {len(remaining_species)} species...")

        if USE_ALARM:
            signal.alarm(DEFAULT_TIMEOUT)            

        try:
            round_csv_data = []
            round_no_seq = []
            round_skipped = []

            ncbi_processed, round_skipped, round_no_seq, round_csv_data = download_mt_genomes_for_chunk(
                category=category,
                species_names=remaining_species,
                email=email,
                api_key=api_key,
                output_dir=output_dir,
                log_fp=log_fp,
                max_records=max_records,
                max_threads=max_threads,
                ncbi_processed=ncbi_processed,
                skipped_species=skipped_species,
                no_sequence_species=no_sequence_species,
                csv_data=[]
            )

            if USE_ALARM:
                signal.alarm(0)                          

            csv_data.extend(round_csv_data)
            skipped_species.extend(round_skipped)
            no_sequence_species.extend(round_no_seq)

            failed_species = [row[0] for row in round_csv_data if row[4] in (
                'Search failed', 'No sequences found', 'No valid record downloaded')]
            remaining_species = failed_species

        except TimeoutException:
            msg = (f"[TIMEOUT] Chunk attempt {attempt} timed out (> {DEFAULT_TIMEOUT}s). "
                   f"Retrying with remaining species...")
            print(msg)
            log_fp.write(msg + "\n")
            log_fp.flush()
            continue

        except Exception as e:
            msg = (f"[ERROR] Chunk attempt {attempt} failed unexpectedly: {e}. Retrying...")
            print(msg)
            log_fp.write(msg + "\n")
            log_fp.flush()
            continue

    if remaining_species:
        for sp in remaining_species:
            no_sequence_species.append(sp)
            csv_data.append([sp, 'NA', 'NA', 'NA', 'Chunk failed after max retries'])

        log_fp.write(f"[ERROR] {len(remaining_species)} species failed after all retries.\n")
        log_fp.flush()

    return ncbi_processed, skipped_species, no_sequence_species, csv_data

def load_species_list(filename):
    species_names = []
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                species_names.append(line)
    return species_names

def chunked_species_list(species_list, chunk_size):
    for i in range(0, len(species_list), chunk_size):
        yield species_list[i:i+chunk_size]

def main():
    parser = argparse.ArgumentParser(
        description="Automatically split a large species list into chunks, "
                    "download mitochondrial genomes with multi-threading, "
                    "and produce a single log and CSV under <output_dir>/outputs/."
    )
    parser.add_argument("-i", "--input_file", required=True,
                        help="Path to the large species list (one species name per line).")
    parser.add_argument("-o", "--output_directory", required=True,
                        help="Directory where 'outputs/' folder will be created.")
    parser.add_argument("-e", "--email", required=True,
                        help="Your NCBI Entrez email.")
    parser.add_argument("-k", "--api_key", default=None,
                        help="Your NCBI Entrez API key (optional but recommended).")
    parser.add_argument("-c", "--category", default="animals",
                        help="NCBI category filter (e.g., animals, plants, fungi...). Default=animals.")
    parser.add_argument("-s", "--chunk_size", type=int, default=DEFAULT_CHUNK_SIZE,
                        help=f"Number of species names per chunk. Default={DEFAULT_CHUNK_SIZE}.")
    parser.add_argument("-m", "--max_records", type=int, default=5,
                        help="Maximum records to fetch for each species. Default=5.")
    parser.add_argument("-t", "--max_threads", type=int, default=8,
                        help="Number of threads for parallel downloading. Default=8.")

    args = parser.parse_args()

    all_entries = load_species_list(args.input_file)
    if not all_entries:
        print("[INFO] No valid entries found in the input file. Exiting.")
        sys.exit(0)

    accession_list   = [x for x in all_entries if is_accession(x)]     ### NEW ###
    species_list_raw = [x for x in all_entries if not is_accession(x)] ### NEW ###

    total_species_entries   = len(species_list_raw)
    total_accession_entries = len(accession_list)

    outputs_dir = os.path.join(args.output_directory, "outputs")
    os.makedirs(outputs_dir, exist_ok=True)

    gb_dir = os.path.join(outputs_dir, "gb_files")
    os.makedirs(gb_dir, exist_ok=True)

    log_path = os.path.join(outputs_dir, "final_download.log")
    csv_path = os.path.join(outputs_dir, "final_download_summary.csv")

    start_time = datetime.datetime.now()

    with open(log_path, 'w') as log_fp:
        log_fp.write(f"Script started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        log_fp.write(f"[INFO] Loaded {total_species_entries} species-name entry(ies) "
                     f"+ {total_accession_entries} accession ID(s).\n")        
        log_fp.flush()

        ncbi_processed = {}
        skipped_species = []
        no_sequence_species = []
        csv_data = []  

        ###################################################################
        if accession_list:
            log_fp.write(f"[INFO] Downloading {len(accession_list)} accession ID(s)…\n")
            log_fp.flush()

            Entrez.email = args.email
            Entrez.api_key = args.api_key

            for acc in accession_list:
                ok = download_single_accession(acc,
                                               output_dir=gb_dir,
                                               log_fp=log_fp,
                                               csv_data=csv_data,
                                               ncbi_processed=ncbi_processed,
                                               skipped_species=skipped_species)

                if not ok:
                    no_sequence_species.append(acc)

        if species_list_raw:
            chunks = list(chunked_species_list(species_list_raw, args.chunk_size))
            log_fp.write(f"[INFO] Species names split into {len(chunks)} chunk(s), "
                         f"chunk_size={args.chunk_size}.\n")
            log_fp.flush()

            for idx, chunk_species in enumerate(chunks, start=1):
                msg_chunk = f"[INFO] Processing chunk #{idx} with {len(chunk_species)} species."
                print(msg_chunk)
                log_fp.write(msg_chunk + "\n")
                log_fp.flush()

                ncbi_processed, skipped_species, no_sequence_species, csv_data = run_chunk_with_timeout_and_retry(
                    species_chunk=chunk_species,
                    category=args.category,
                    email=args.email,
                    api_key=args.api_key,
                    output_dir=gb_dir,
                    log_fp=log_fp,
                    max_records=args.max_records,
                    max_threads=args.max_threads,
                    ncbi_processed=ncbi_processed,
                    skipped_species=skipped_species,
                    no_sequence_species=no_sequence_species,
                    csv_data=csv_data
                )

        gb_files = [f for f in os.listdir(gb_dir) if f.endswith(".gb")]
        total_downloaded_sequences = len(gb_files)
        end_time = datetime.datetime.now()
        duration = end_time - start_time

        log_fp.write("\n----------------------------------------\n")
        log_fp.write("Summary:\n\n")
        log_fp.write(f"Total species-name entries in input list: {total_species_entries}\n")
        log_fp.write(f"Total accession-ID entries in input list: {total_accession_entries}\n")
        log_fp.write(
            f"Unique NCBI species with sequences found: {len(ncbi_processed)}\n"
            f"Total sequences downloaded: {total_downloaded_sequences}\n\n"
        )

        log_fp.write(f"Species skipped due to duplication: {len(skipped_species)}\n")
        for sp in skipped_species:
            log_fp.write(f"  {sp}\n")
        log_fp.write("\n")

        failed_notes = {
            "Search failed",
            "No sequences found",
            "No valid record downloaded",
            "Chunk failed after max retries",
            "Failed to download accession"
        }
        failed_species_count = len({row[0] for row in csv_data if row[4] in failed_notes})
        log_fp.write(f"Entries with no sequences found / failed: {failed_species_count}\n")
        log_fp.write("----------------------------------------\n\n")
        log_fp.write(
            f"Script finished at {end_time.strftime('%Y-%m-%d %H:%M:%S')}, "
            f"duration: {duration}\n"
        )
        log_fp.flush()


    with open(csv_path, 'w', newline='') as out_csv:
        writer = csv.writer(out_csv)
        writer.writerow([
            'Listed Species Name (or Accession)',
            'NCBI Species Name',
            'Number of Sequences Downloaded',
            'VERSION',
            'Notes'
        ])
        csv_data_sorted = sorted(csv_data, key=lambda row: row[0])
        writer.writerows(csv_data_sorted)

    failed_species_path = Path(outputs_dir) / "failed_species.txt"
    failed_species_set = {row[0] for row in csv_data if row[4] in failed_notes}
    with failed_species_path.open("w") as f:
        for sp in sorted(failed_species_set):
            f.write(f"{sp}\n")

    print("[INFO] Failed-species list saved to:", failed_species_path)
    print("[INFO] All downloads attempted. Check the log for details.")
    print("[INFO] Final log:", log_path)
    print("[INFO] Final CSV:", csv_path)


if __name__ == "__main__":
    main()