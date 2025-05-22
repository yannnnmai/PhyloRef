#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import subprocess
import sys
from pathlib import Path
from datetime import datetime
import os
from Bio import SeqIO
from io import StringIO
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

def check_dependencies(mafft_path: Path, fasttree_path: Path) -> None:
    """
    Check whether the specified MAFFT and FastTree executables are available and executable.

    Parameters
    ----------
    mafft_path : Path
        Path to the MAFFT executable.
    fasttree_path : Path
        Path to the FastTree executable.

    Raises
    ------
    FileNotFoundError
        If either MAFFT or FastTree is not found at the specified path
        or does not have execute permissions.
    """
    if not mafft_path.is_file() or not os.access(mafft_path, os.X_OK):
        raise FileNotFoundError(f"MAFFT executable not found or not executable: {mafft_path}")
    if not fasttree_path.is_file() or not os.access(fasttree_path, os.X_OK):
        raise FileNotFoundError(f"FastTree executable not found or not executable: {fasttree_path}")
    
def run_alignment(input_file: Path, alignment_dir: Path, mafft_path: Path, threads: int) -> Path:
    output_path = alignment_dir / (input_file.stem + '_aligned.fa')
    log_path = alignment_dir / (input_file.stem + '_alignment.log')

    start_time = datetime.now()
    try:
        mafft_command = [
            str(mafft_path), 
            "--thread", str(threads), 
            str(input_file)
        ]
        result = subprocess.run(mafft_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        with open(output_path, 'w') as aligned_file:
            aligned_file.write(result.stdout)

        end_time = datetime.now()
        total_time = end_time - start_time

        with open(log_path, 'w') as log_file:
            log_file.write(f"Start Time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            log_file.write(f"End Time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            log_file.write(f"Total running time: {total_time}\n")
            if result.stderr:
                log_file.write("\n=== STDERR ===\n")
                log_file.write(result.stderr)

        return output_path

    except subprocess.CalledProcessError as e:
        error_message = f"MAFFT alignment failed: {e.stderr.decode().strip()}"
        raise RuntimeError(error_message) from e
    
def build_fasttree_tree(input_file: Path, tree_dir: Path, fasttree_path: Path) -> Path:
    output_file = tree_dir / (input_file.stem + '_fasttree.nwk')
    log_path = tree_dir / (input_file.stem + '_fasttree.log')

    start_time = datetime.now()
    try:
        command = [
            str(fasttree_path),
            "-fastest", 
            "-nt", 
            "-log", str(log_path),
            str(input_file)
        ]
        with open(output_file, 'w') as out_f:
            result = subprocess.run(command, check=True, stdout=out_f, stderr=subprocess.PIPE, text=True)

        end_time = datetime.now()
        total_time = end_time - start_time

        with open(log_path, 'w') as log_file: 
            log_file.write(f"Start Time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            log_file.write(f"End Time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            log_file.write(f"Total running time: {total_time}\n")
            
            if result.stderr:
                log_file.write("\n=== STDERR ===\n")
                log_file.write(result.stderr)

        return output_file

    except subprocess.CalledProcessError as e:
        error_message = f"FastTree tree construction failed: {e.stderr.decode().strip()}"
        raise RuntimeError(error_message) from e

def parse_arguments() -> tuple:
    FASTA_EXTS = {".fa", ".fasta", ".fas", ".fna"}

    parser = argparse.ArgumentParser(
        description="Perform alignment and tree inference with MAFFT and FastTree.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        nargs="+",
        help="Path(s) to input FASTA file(s) or directory. Multiple inputs supported.",
    )
    parser.add_argument(
        "-o", "--output", required=True, type=Path, help="Path to the output directory."
    )
    parser.add_argument(
        "-d",
        "--mode",
        choices=["align", "tree", "both"],
        default="both",
        help="Operation mode: 'align', 'tree', or 'both'. Default is 'both'.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=30,
        help="Number of threads to use for MAFFT alignment. Default is 30.",
    )
    parser.add_argument(
        "-mf",
        "--mafft",
        type=Path,
        default=Path("/usr/bin/mafft"),
        help="Path to the MAFFT executable. If in $PATH, just use 'mafft'.",
    )
    parser.add_argument(
        "-ft",
        "--fasttree",
        type=Path,
        default=Path("~/users/zhoutao/conda/bin/fasttree").expanduser(),
        help="Path to the FastTree executable. If in $PATH, just use 'FastTree'.",
    )
    args = parser.parse_args()

    expanded_inputs = []
    for p in args.input:
        if p.is_file():
            expanded_inputs.append(p)
        elif p.is_dir():
            for f in p.iterdir():
                if f.suffix.lower() in FASTA_EXTS and f.is_file():
                    expanded_inputs.append(f)
        else:
            parser.error(f"Path does not exist: {p}")

    if not expanded_inputs:
        parser.error("No FASTA files found. Please check the -i / --input argument.")

    args.input = list(dict.fromkeys(expanded_inputs))

    outputs_dir = args.output / "outputs"
    alignment_dir = outputs_dir / "alignment"
    tree_dir = outputs_dir / "fasttree"
    alignment_dir.mkdir(parents=True, exist_ok=True)
    tree_dir.mkdir(parents=True, exist_ok=True)

    if args.threads < 1:
        parser.error("Number of threads must be at least 1.")

    return args, alignment_dir, tree_dir

def process_file(input_file: Path, args: argparse.Namespace, alignment_dir: Path, tree_dir: Path) -> dict:
    result = {
        'input_file': str(input_file),
        'aligned_output': None,
        'tree_output': None,
        'error': None
    }

    try:
        aligned_file = None
        # 1) If mode is 'align' or 'both', perform MAFFT alignment.
        if args.mode in ('align', 'both'):
            aligned_file = run_alignment(input_file, alignment_dir, args.mafft, args.threads)
            result['aligned_output'] = str(aligned_file)

        # 2) If mode is 'tree' or 'both', perform FastTree phylogenetic inference.
        if args.mode in ('tree', 'both'):
            #    - If mode is 'both', use the aligned file from MAFFT as input for FastTree.
            #    - If mode is 'tree', use the original input_file directly for tree construction.   
            tree_input = aligned_file if args.mode == 'both' else input_file
            fasttree_file = build_fasttree_tree(tree_input, tree_dir, args.fasttree)
            result['tree_output'] = str(fasttree_file)

    except Exception as e:
        result['error'] = str(e)

    return result

def write_log(log_path: Path, start_time: datetime, end_time: datetime, results: list, errors: list) -> None:
    total_time = end_time - start_time
    with log_path.open('w') as log_file:
        log_file.write(f"Script Start Time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        for res in results:
            log_file.write(f"File: {res['input_file']}\n")
            if res['aligned_output']:
                log_file.write(f"Aligned Output: {res['aligned_output']}\n")
            if res['tree_output']:
                log_file.write(f"FastTree Output: {res['tree_output']}\n")
            log_file.write("\n")

        if errors:
            log_file.write("Errors:\n")
            for err in errors:
                log_file.write(f"File: {err['input_file']}\n")
                log_file.write(f"Error: {err['error']}\n\n")

        log_file.write(f"Script End Time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        log_file.write(f"Total Running Time: {total_time}\n")

def main():
    args, alignment_dir, tree_dir = parse_arguments()

    outputs_dir = args.output / "outputs"
    outputs_dir.mkdir(parents=True, exist_ok=True)

    log_path = outputs_dir / "process.log"

    if args.mode == "align":
        tree_dir.rmdir() 
    elif args.mode == "tree":
        alignment_dir.rmdir()  

    start_time = datetime.now()

    results = []
    errors = []

    try:
        check_dependencies(args.mafft, args.fasttree)

        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            future_to_file = {
                executor.submit(process_file, input_file, args, alignment_dir, tree_dir): input_file
                for input_file in args.input
            }

            for future in tqdm(as_completed(future_to_file), total=len(future_to_file), desc="处理文件", unit="file"):
                res = future.result()
                if res['error']:
                    errors.append(res)
                else:
                    results.append(res)

        end_time = datetime.now()

        write_log(log_path, start_time, end_time, results, errors)
        print(f"\nProcessing complete ✅")
        print(f"Results saved to: {outputs_dir}")
        print(f"Log saved to: {log_path}\n")

    except Exception as e:
        end_time = datetime.now()
        write_log(log_path, start_time, end_time, results, errors + [{'input_file': 'N/A', 'error': str(e)}])
        print(f"An error occurred during processing. See log file for details: {log_path}")
        sys.exit(1)

if __name__ == "__main__":
    main()