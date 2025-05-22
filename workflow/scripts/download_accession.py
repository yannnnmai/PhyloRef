#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Description:
    This script downloads sequences from NCBI by a large list of accession
    numbers. It splits the list into chunks (default=1000) to handle large
    input more reliably. Then, each chunk is processed with multi-threading
    and includes a timeout + retry mechanism.

    Results are saved in:
        <output_directory>/outputs/
            ├─ final_download.log          (unified log file)
            ├─ final_download_summary.csv  (unified CSV summary)
            └─ gb_files/                  (all .gb files)

Usage Example:
    python download_accession.py \
        -i /path/to/large_accession_list.txt \
        -o /path/to/output_directory \
        -e your_email@example.com \
        -k your_ncbi_api_key \
        -s 1000 \  # Chunk size (default: 1000)
        -t 8  # Number of threads (recommended value; not intended to be changed)
"""

import os
import re
import sys
import csv
import time
import argparse
import datetime
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

try:
    from Bio import Entrez
except ImportError:
    sys.exit("Biopython is required. Please install it via: pip install biopython")

MAX_RETRIES = 3              
BACKOFF = [5, 10, 20]        

def load_accession_ids(input_file):
    with open(input_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]


def split_accession_ids(accession_ids, chunk_size):
    for i in range(0, len(accession_ids), chunk_size):
        yield accession_ids[i:i+chunk_size]


def robust_fetch(db, id, rettype, retmode):
    for attempt in range(MAX_RETRIES):
        try:
            with Entrez.efetch(db=db, id=id,
                               rettype=rettype, retmode=retmode) as h:
                return h.read()
        except Exception as e:
            wait = BACKOFF[min(attempt, len(BACKOFF) - 1)]
            print(f"[WARN] {id} attempt {attempt + 1}/{MAX_RETRIES} failed: {e}",
                  flush=True)
            time.sleep(wait)
    return None


def process_sequence(seq_id, gb_folder):
    data = robust_fetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
    if data is None:
        return None, None

    species_name = None
    locus_number = None
    for line in data.split('\n'):
        if line.startswith("  ORGANISM"):
            species_name = ' '.join(line.split()[1:])
        elif line.startswith("LOCUS"):
            locus_number = line.split()[1]

    if not locus_number:
        return None, None

    safe_species_name = re.sub(r"[^A-Za-z0-9._-]", "_",
                               species_name) if species_name else "Unknown"
    file_name = f"{safe_species_name}|{seq_id}.gb"   
    file_path = os.path.join(gb_folder, file_name)

    try:
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(data)
        return species_name, seq_id
    except Exception as e:
        print(f"[ERROR] Cannot write {seq_id} → {file_path}: {e}", flush=True)
        return None, None

def download_sequences_in_chunk(chunk_accessions,
                                email,
                                api_key,
                                gb_folder,
                                max_threads=8):

    Entrez.email = email
    Entrez.api_key = api_key
    Entrez.sleep_between_tries = 0.2

    if not os.path.exists(gb_folder):
        os.makedirs(gb_folder, exist_ok=True)

    success_dict = defaultdict(set)
    fail_list = []

    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        future_to_seq_id = {executor.submit(process_sequence, seq_id, gb_folder): seq_id
                            for seq_id in chunk_accessions}
        for future in as_completed(future_to_seq_id):
            seq_id = future_to_seq_id[future]
            try:
                sp_name, acc_id = future.result()
                if sp_name:
                    success_dict[sp_name].add(acc_id)
                else:
                    fail_list.append(seq_id)
            except Exception as e:
                print(f"[ERROR] Accession {seq_id} crashed: {e}", flush=True)
                fail_list.append(seq_id)
    return success_dict, fail_list


def main():
    parser = argparse.ArgumentParser(
        description="Split a large accession list into chunks and download them all into one folder under outputs/."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Path to the file containing large accession numbers (one per line).")
    parser.add_argument("-o", "--output", required=True,
                        help="Directory where results will be saved inside 'outputs/'.")
    parser.add_argument("-e", "--email", required=True,
                        help="NCBI Entrez email.")
    parser.add_argument("-k", "--api_key", default=None,
                        help="NCBI Entrez API key (recommended).")
    parser.add_argument("-s", "--chunk_size", type=int, default=1000,
                        help="Number of accessions per chunk. Default=1000.")
    parser.add_argument("-t", "--threads", type=int, default=8,
                        help="Number of threads for parallel downloading. Default=8 (recommended ≤8).")
    
    args = parser.parse_args()

    all_accessions = load_accession_ids(args.input)
    if not all_accessions:
        print("[INFO] No valid accession found in the input file. Exiting.")
        sys.exit(0)

    outputs_dir = os.path.join(args.output, "outputs")
    os.makedirs(outputs_dir, exist_ok=True)

    gb_folder = os.path.join(outputs_dir, "gb_files")
    os.makedirs(gb_folder, exist_ok=True)

    log_path = os.path.join(outputs_dir, "final_download.log")
    csv_path = os.path.join(outputs_dir, "final_download_summary.csv")

    start_time = datetime.datetime.now()
    with open(log_path, 'w') as log_fp:
        log_fp.write(f"Script started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        log_fp.flush()

        chunks = list(split_accession_ids(all_accessions, args.chunk_size))
        msg = f"[INFO] Loaded {len(all_accessions)} accession(s). Split into {len(chunks)} chunk(s)."
        print(msg)
        log_fp.write(msg + "\n")
        log_fp.flush()

        global_species_dict = defaultdict(set)
        global_fail_list = []

        for idx, chunk_list in enumerate(chunks, start=1):
            chunk_msg = f"[INFO] Processing chunk {idx}/{len(chunks)} with {len(chunk_list)} accession(s)."
            print(chunk_msg)
            log_fp.write(chunk_msg + "\n")
            log_fp.flush()

            s_dict, f_list = download_sequences_in_chunk(
                chunk_accessions=chunk_list,
                email=args.email,
                api_key=args.api_key,
                gb_folder=gb_folder,
                max_threads=args.threads
            )

            log_fp.write(f"[INFO] Chunk {idx} finished: {len(f_list)} failed.\n")
            log_fp.flush()

            for sp_name, acc_set in s_dict.items():
                global_species_dict[sp_name].update(acc_set)
            global_fail_list.extend(f_list)

        msg_fail = f"[INFO] Download completed. {len(global_fail_list)} accession(s) failed in total."
        print(msg_fail)
        log_fp.write(msg_fail + "\n")

        end_time = datetime.datetime.now()
        duration = end_time - start_time
        finish_msg = (f"\nScript finished at {end_time.strftime('%Y-%m-%d %H:%M:%S')}, "
                      f"duration: {duration}\n")
        print(finish_msg)
        log_fp.write(finish_msg)
        log_fp.flush()

    with open(csv_path, 'w', newline='') as csv_fp:
        writer = csv.writer(csv_fp)
        writer.writerow(["Species Name", "Sequence Count", "Accession IDs"])

        for sp_name in sorted(global_species_dict.keys()):
            acc_list = sorted(global_species_dict[sp_name])
            writer.writerow([sp_name, len(acc_list), ", ".join(acc_list)])

    if global_fail_list:
        fail_txt = os.path.join(outputs_dir, "failed_accessions.txt")
        with open(fail_txt, 'w') as f_txt:
            for acc in global_fail_list:
                f_txt.write(acc + "\n")

    print("[INFO] All results are saved in:", outputs_dir)
    sys.exit(1 if global_fail_list else 0)

if __name__ == "__main__":
    main()