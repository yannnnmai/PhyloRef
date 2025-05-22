#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Mitochondrial Gene Extraction Script
====================================

This script extracts specified mitochondrial gene sequences from GenBank (.gb or .genbank) files
located in a given directory. It automatically outputs raw gene sequences, concatenated sequences,
and a log of missing genes.

Usage Example:
    python extract_genes.py -i /path/to/genbank_files -o /path/to/outputs -g 12S COX1 ND2

Arguments:
    -i, --input      Path to the input directory containing .gb or .genbank files.
    -o, --outputs    Path to the output directory. A subdirectory named 'outputs/' will be created automatically.
    -g, --genes      List of genes to extract (case-insensitive). If not specified, all supported genes will be extracted.
                    Supported genes include: 12S, 16S, ND1, ND2, COX1, COX2, ATP8, ATP6, COX3, ND3, ND4L, ND4, ND5, ND6, CYTB
"""

import os
import re
import logging
import datetime
from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError
from pathlib import Path
from collections import deque
import argparse
import csv

GENE_PATTERNS = {
    "12S": r"12S ribos?o?m[ae]l\s+rna|12S rRNA|s-rRNA|rrn12|small subunit ribosomal rna|"
           r"small ribosomal rna subunit rna|12S ribosomal rna subunit|12S ribosormal rna|"
           r"12S small subunit|ssu rRNA|rrnS|12S|l2S ribosomal RNA",
    "16S": r"16S ribos?o?m[ae]l\s+rna|16S rRNA|l-rRNA|rrn16|large subunit ribosomal rna|"
           r"large ribosomal rna subunit rna|16S ribosomal rna subunit|16S ribosormal rna|"
           r"16S large subunit|lsu rRNA|rrnL|16S|l6S ribosomal RNA",
    "ND1": r"nd[-_]?1|nad1|nadh1|NADH dehydrogenase subunit 1|nd1|NADH hydrogenase subunit 1",
    "ND2": r"nd[-_]?2|nad2|nadh2|NADH dehydrogenase subunit 2|nd2|NADH hydrogenase subunit 2",
    "COX1": r"co[-_]?1|cox1|coi|cytochrome c oxidase (subunit )?i|COX I|COXI|cytochrome oxidase I|"
            r"COX1|cytochrome oxidase subunit 1|cytochrome oxidase subunit I|CO I|cytochrome c oxidase subunit 1",
    "COX2": r"co[-_]?2|cox2|coii|cytochrome c oxidase (subunit )?ii|COX II|COXII|cytochrome oxidase II|"
            r"COX2|cytochrome oxidase subunit 2|cytochrome oxidase subunit II|CO II|cytochrome c oxidase subunit 2",
    "ATP8": r"atp[-_]?8|atpase8|ATP synthase (subunit )?8|ATPase subunit 8|atp8|ATP synthase F0 subunit 8|"
            r"ATPase 8|ATPase subuint 8",
    "ATP6": r"atp[-_]?6|atpase6|ATP synthase (subunit )?6|ATPase subunit 6|atp6|ATP synthase F0 subunit 6|"
            r"ATPase 6|ATPase subuint 6",
    "COX3": r"co[-_]?3|cox3|coiii|cytochrome c oxidase (subunit )?iii|COX III|COXIII|cytochrome oxidase III|"
            r"COX3|cytochrome oxidase subunit 3|cytochrome oxidase subunit III|CO III|cytochrome c oxidase subunit 3",
    "ND3": r"nd[-_]?3|nad3|nadh3|NADH dehydrogenase subunit 3|nd3|NADH hydrogenase subunit 3",
    "ND4L": r"nd[-_]?4l|nad4l|nadh4l|NADH dehydrogenase subunit 4L|nd4l|NADH hydrogenase subunit 4L",
    "ND4": r"nd[-_]?4|nad4|nadh4|NADH dehydrogenase subunit 4|nd4|NADH hydrogenase subunit 4",
    "ND5": r"nd[-_]?5|nad5|nadh5|NADH dehydrogenase subunit 5|nd5|NADH hydrogenase subunit 5",
    "ND6": r"nd[-_]?6|nad6|nadh6|NADH dehydrogenase subunit 6|nd6|NADH hydrogenase subunit 6",
    "CYTB": r"cytb|cytochrome b|cytochrome b-245 heavy chain|cytb-245|cytochrome b subunit|cytb"
}

class CustomFormatter(logging.Formatter):
    def format(self, record):
        level_to_prefix = {
            logging.DEBUG: "DEBUG - ",
            logging.INFO: "",
            logging.WARNING: "WARNING - ",
            logging.ERROR: "ERROR - "
        }
        prefix = level_to_prefix.get(record.levelno, "")
        if record.levelno == logging.INFO and ("Run started" in record.msg or "Run ended" in record.msg):
            prefix = ""
        return f"{prefix}{record.msg}"


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Extract specified mitochondrial gene sequences from GenBank files and save the results to the output directory. "
            "Supported genes include: 12S, 16S, ND1, ND2, COX1, COX2, ATP8, ATP6, "
            "COX3, ND3, ND4L, ND4, ND5, ND6, CYTB."
        )
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to the input directory containing .gb or .genbank files."
    )
    parser.add_argument(
        "-o", "--outputs", required=True,
        help="Path to the output directory. A subdirectory named 'outputs/' will be created automatically."
    )
    parser.add_argument(
        "-g", "--genes", nargs="+", default=None,
        help="List of genes to extract (case-insensitive). If not specified, all supported genes will be extracted."
    )
    return parser.parse_args()

def build_taxonomy_maps(folder_path):
 
    genus_to_family = {}
    family_to_order = {}

    for file_name in os.listdir(folder_path):
        if file_name.endswith(".gb") or file_name.endswith(".genbank"):
            file_path = os.path.join(folder_path, file_name)
            with open(file_path, 'r') as file:
                for record in SeqIO.parse(file, 'genbank'):
                    taxonomy = record.annotations.get('taxonomy', [])
                    genus = family = order = 'NA'
                    for taxon in taxonomy:
                        if taxon.endswith('idae'):
                            family = taxon
                        elif taxon.endswith('formes'):
                            order = taxon
                    if 'organism' in record.annotations:
                        parts = record.annotations['organism'].split()
                        if len(parts) > 1:
                            genus = parts[0]

                    if genus != 'NA' and family != 'NA':
                        genus_to_family[genus] = family
                    if family != 'NA' and order != 'NA':
                        family_to_order[family] = order

    return genus_to_family, family_to_order

def safe_str(seq, logger, context=""):

    try:
        return str(seq)
    except UndefinedSequenceError:
        logger.warning(f"Undefined sequence encountered in {context}. Skipping this record.")
        return None

def extract_and_concatenate_genes(seq_record, gene_patterns, selected_genes, output_dir,
                                  base_filename, concat_file, taxonomy_data, logger, full_seq_str):
    all_genes_order = list(gene_patterns.keys())
    sorted_genes = sorted(selected_genes, key=lambda gene: all_genes_order.index(gene))

    concat_sequences = []
    missing_genes = []
    is_any_gene_missing = False

    taxonomy_suffix = f"|{taxonomy_data.get('order', 'o__NA')}|{taxonomy_data.get('family', 'f__NA')}|{taxonomy_data.get('genus', 'g__NA')}"

    for gene in sorted_genes:
        pattern = gene_patterns.get(gene)
        gene_found = False
        sequence = None

        for feature in seq_record.features:
            qualifiers = (feature.qualifiers.get('product', []) +
                          feature.qualifiers.get('gene', []) +
                          feature.qualifiers.get('note', []))
            if any(re.search(pattern, qualifier, re.IGNORECASE) for qualifier in qualifiers):
                try:
                    extracted = feature.extract(seq_record.seq)
                    gene_str = safe_str(extracted, logger, f"Feature extraction for gene {gene} in {base_filename}")
                except Exception as e:
                    logger.warning(f"Error extracting gene {gene} from {base_filename}: {e}")
                    gene_str = None
                if gene_str is not None:
                    sequence = gene_str
                    gene_found = True
                    break

        if not gene_found:
            missing_genes.append(gene)
            is_any_gene_missing = True  
            sequence = full_seq_str

        concat_sequences.append(sequence)

        with open(output_dir / f"{gene}.fa", 'a') as gene_file:
            gene_file.write(f">{base_filename}|{gene}{taxonomy_suffix}\n{sequence}\n")

    if is_any_gene_missing:
        concat_file.write(
            f">{base_filename}|Concat|{taxonomy_data.get('order', 'o__NA')}|"
            f"{taxonomy_data.get('family', 'f__NA')}|{taxonomy_data.get('genus', 'g__NA')}\n"
            f"{full_seq_str}\n"
        )
    else:
        concat_file.write(
            f">{base_filename}|Concat|{taxonomy_data.get('order', 'o__NA')}|"
            f"{taxonomy_data.get('family', 'f__NA')}|{taxonomy_data.get('genus', 'g__NA')}\n"
            f"{''.join(concat_sequences)}\n"
        )

    return missing_genes

def fill_missing_taxonomy(folder_path, genus_to_family, family_to_order, output_csv, log_file):
    headers = ['version', 'order', 'family', 'genus', 'species']
    all_taxonomy_data = []
    missing_info = {'order': [], 'family': [], 'genus': []}

    with open(output_csv, 'w', newline='', encoding='utf-8') as csvfile, open(log_file, 'w', encoding='utf-8') as logf:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()

        for file_name in os.listdir(folder_path):
            if file_name.endswith(".gb") or file_name.endswith(".genbank"):
                file_path = os.path.join(folder_path, file_name)
                with open(file_path, 'r') as file:
                    for record in SeqIO.parse(file, 'genbank'):
                        version = record.annotations['accessions'][0] if 'accessions' in record.annotations else 'NA'
                        taxonomy = record.annotations.get('taxonomy', [])
                        genus = species = family = order = 'NA'

                        for taxon in taxonomy:
                            if taxon.endswith('idae'):
                                family = f"f__{taxon}"
                            elif taxon.endswith('formes'):
                                order = f"o__{taxon}"

                        if 'organism' in record.annotations:
                            parts = record.annotations['organism'].split()
                            if len(parts) > 1:
                                genus = f"g__{parts[0]}"
                                species = f"s__{' '.join(parts[:2])}"

                        if family == 'f__NA' and genus in genus_to_family:
                            family = f"f__{genus_to_family[genus]}"

                        if order == 'o__NA':
                            family_name = family.replace("f__", "")
                            if family_name in family_to_order:
                                order = f"o__{family_to_order[family_name]}"

                        row = {
                            'version': version,
                            'order': order,
                            'family': family,
                            'genus': genus,
                            'species': species
                        }
                        writer.writerow(row)
                        all_taxonomy_data.append(row)

                        for key, value in {'order': order, 'family': family, 'genus': genus}.items():
                            if value.endswith('NA'):
                                missing_info[key].append(file_name)

        for key, files in missing_info.items():
            if files:
                logf.write(f"Missing {key} in files:\n")
                for f in files:
                    logf.write(f"WARNING - {f}\n")
                logf.write("\n")

    return all_taxonomy_data

def main():
    args = parse_arguments()

    all_available_genes = list(GENE_PATTERNS.keys())
    if args.genes:
        selected_genes = [g.upper() for g in args.genes if g.upper() in all_available_genes]
        if not selected_genes:
            selected_genes = all_available_genes
    else:
        selected_genes = all_available_genes

    input_dir = Path(args.input)
    if not input_dir.is_dir():
        print("Error: The input path does not exist or is not a directory.")
        return

    output_base = Path(args.outputs)
    output_dir = output_base / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger('gene_extraction')
    logger.setLevel(logging.DEBUG)
    log_path = output_dir / 'gene_extraction.log'
    file_handler = logging.FileHandler(log_path, mode='w', encoding='utf-8')
    file_handler.setFormatter(CustomFormatter())
    logger.addHandler(file_handler)

    start_time = datetime.datetime.now()

    genus_to_family, family_to_order = build_taxonomy_maps(str(input_dir))

    complete_fa_path = output_dir / "Complete.fa"
    concat_fa_path = output_dir / "Concat.fa"
    with open(complete_fa_path, 'w') as all_seq_file, open(concat_fa_path, 'w') as concat_file:
        records = []
        for genbank_file in input_dir.glob('*.gb'):
            seq_record = SeqIO.read(genbank_file, "genbank")
            records.append((seq_record, genbank_file.stem))
        for genbank_file in input_dir.glob('*.genbank'):
            seq_record = SeqIO.read(genbank_file, "genbank")
            records.append((seq_record, genbank_file.stem))

        records.sort(key=lambda x: x[0].annotations.get('organism', ''))

        debug_logs = deque()
        total_files_processed = 0
        total_missing_genes = 0

        for seq_record, base_filename in records:
            full_seq_str = safe_str(seq_record.seq, logger, f"Complete sequence for {base_filename}")
            if full_seq_str is None:
                logger.warning(f"Skipping file {base_filename}因为其完整序列未定义。")
                continue

            taxonomy_data = {'order': 'o__NA', 'family': 'f__NA', 'genus': 'g__NA'}
            for taxon in seq_record.annotations.get('taxonomy', []):
                if taxon.endswith('idae'):
                    taxonomy_data['family'] = f"f__{taxon}"
                elif taxon.endswith('formes'):
                    taxonomy_data['order'] = f"o__{taxon}"
            if 'organism' in seq_record.annotations:
                parts = seq_record.annotations['organism'].split()
                if len(parts) > 1:
                    taxonomy_data['genus'] = f"g__{parts[0]}"

            raw_family = taxonomy_data['family'].replace("f__", "")
            raw_genus = taxonomy_data['genus'].replace("g__", "")
            if raw_family == 'NA' and raw_genus in genus_to_family:
                taxonomy_data['family'] = f"f__{genus_to_family[raw_genus]}"
            raw_family = taxonomy_data['family'].replace("f__", "")
            if taxonomy_data['order'] == 'o__NA' and raw_family in family_to_order:
                taxonomy_data['order'] = f"o__{family_to_order[raw_family]}"

            taxonomy_suffix = (f"|Complete|{taxonomy_data['order']}|"
                               f"{taxonomy_data['family']}|{taxonomy_data['genus']}")
            all_seq_file.write(
                f">{base_filename}{taxonomy_suffix}\n{full_seq_str}\n"
            )

            missing_genes = extract_and_concatenate_genes(
                seq_record, GENE_PATTERNS, selected_genes, output_dir,
                base_filename, concat_file, taxonomy_data, logger, full_seq_str
            )

            if missing_genes:
                debug_msg = f"Missing genes in {base_filename}: {', '.join(missing_genes)}"
                debug_logs.append(debug_msg)
                for mg in missing_genes:
                    debug_logs.append(
                        f"Gene {mg} in {base_filename} is missing and has been replaced with full length sequence."
                    )
                total_missing_genes += len(missing_genes)

            total_files_processed += 1

        end_time = datetime.datetime.now()
        run_time = end_time - start_time

        logger.info("========================================================")
        logger.info("      Mitochondrial Gene Extraction Report")
        logger.info("========================================================")
        logger.info(f"Start Time           : {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"Total Sequences      : {total_files_processed}")
        logger.info(f"Total Missing Genes  : {total_missing_genes}")
        logger.info("--------------------------------------------------------")
        logger.info("Missing Gene Report:")
        logger.info("--------------------------------------------------------")

        if total_missing_genes == 0:
            logger.info("\nNo missing genes found. All genes extracted successfully.")
        else:
            parsed_missing = {}  
            for msg in debug_logs:
                if msg.startswith("Missing genes in "):
                    base_f = msg.split("Missing genes in ")[1].split(": ")[0]
                    genes_str = msg.split(": ")[1]
                    genes = [g.strip() for g in genes_str.split(",")]
                    parsed_missing[base_f] = genes
                elif msg.startswith("Gene "):
                    pass

            for base_f, genes in parsed_missing.items():
                logger.info(f"\nFile: {base_f}")
                for g in genes:
                    logger.info(f"    - {g} : Missing, full-length sequence used.")

        logger.info(f"\n--------------------------------------------------------")
        logger.info(f"End Time             : {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"Total Run Time       : {run_time}")
        logger.info("========================================================\n")

    taxonomy_output_csv_path = output_dir / 'taxonomy_output.csv'
    missing_taxonomy_log_path = output_dir / 'missing_taxonomy.log'
    fill_missing_taxonomy(
        str(input_dir), genus_to_family, family_to_order,
        str(taxonomy_output_csv_path), str(missing_taxonomy_log_path)
    )

    print(f"[INFO] Processing complete.")
    print(f"[INFO] Log saved to: {log_path}")
    
if __name__ == "__main__":
    main()