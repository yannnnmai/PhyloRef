#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import argparse
import logging
import shutil
from pathlib import Path
from Bio import SeqIO
from ete3 import Tree  # type: ignore

from datetime import datetime
def get_version_suffix():
    return "V" + datetime.now().strftime("%y%m%d")

deleted_accessions = set()
clean_accessions = set()

def parse_species_from_header(header: str) -> str:
    parts = header.split("|")
    if len(parts) < 1:
        return ""
    return parts[0].strip()

def extract_accession(header: str):
    parts = header.split("|")
    if len(parts) < 2:
        return None
    return re.sub(r"\.\d+$", "", parts[1].strip())

def find_representative(seq_list):
    nc_candidates = [x for x in seq_list if 'NC_' in x]
    if nc_candidates:
        return sorted(nc_candidates)[0]
    else:
        return sorted(seq_list)[0]

def build_similar_annotation_for_acc(acc, acc_species, groups_for_acc, all_groups, acc_to_species):
    if not groups_for_acc:
        return ""  

    all_reps = set()
    for gidx in groups_for_acc:
        group = all_groups[gidx] 
        species_dict = {}
        for a in group:
            sp = acc_to_species.get(a, "")
            species_dict.setdefault(sp, []).append(a)

        reps_for_species = {}
        for sp, acclist in species_dict.items():
            reps_for_species[sp] = find_representative(acclist)

        reps_others = []
        for sp, repacc in reps_for_species.items():
            if sp != acc_species:
                reps_others.append(repacc)

        if reps_others:
            all_reps.update(reps_others)

    if not all_reps:
        return "similar_to=None"
    else:
        reps_sorted = sorted(all_reps)
        return "similar_to=" + ",".join(reps_sorted)

def load_similar_file(filepath: str):
    groups = []
    with open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            group = sorted(set(line.split()))
            groups.append(group)
    return groups

def merge_similar_groups(base_groups, new_groups):
    base_groups.extend(new_groups)
    return base_groups

def load_similar_all(files):
    final_groups = []
    for f in files:
        g = load_similar_file(f)
        final_groups = merge_similar_groups(final_groups, g)
    return final_groups

def load_anomalies_file(filepath: str):
    result = set()
    with open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if re.match(r"^(Green|Blue|Red) anomalies for", line) or line == "None":
                continue
            parts = line.split("|")
            if len(parts) >= 2:
                acc = re.sub(r"\.\d+$", "", parts[1].strip())
                result.add(acc)
            else:
                result.add(line.strip())
    return result

def load_anomalies_all(files):
    all_set = set()
    for f in files:
        subset = load_anomalies_file(f)
        all_set |= subset
    return all_set

def load_species_map_all(files):
    species_map = {}
    for f in files:
        with open(f, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if (not line) or re.match(r"^(Green|Blue|Red) anomalies for", line) or line == "None":
                    continue
                parts = line.split("|")
                if len(parts) >= 2:
                    sp = parts[0].strip()
                    acc = re.sub(r"\.\d+$", "", parts[1].strip())
                    species_map[acc] = sp
    return species_map

def process_fasta(infile: Path, out_cleaned: Path, out_blacklisted: Path,
                  remove_set: set, similar_groups: list, logger: logging.Logger):
    global deleted_accessions, clean_accessions
    total, black_count, modified = 0, 0, 0

    records = list(SeqIO.parse(str(infile), "fasta"))
    acc_to_species = {}
    for rec in records:
        acc = extract_accession(rec.description)
        sp = parse_species_from_header(rec.description)
        if acc:
            acc_to_species[acc] = sp

    acc_to_groups = {}
    for g_idx, group in enumerate(similar_groups):
        for a in group:
            acc_to_groups.setdefault(a, []).append(g_idx)

    with open(infile, "r", encoding="utf-8") as fin, \
         open(out_cleaned, "w", encoding="utf-8") as foutC, \
         open(out_blacklisted, "w", encoding="utf-8") as foutB:

        for rec in records:
            total += 1
            acc = extract_accession(rec.description)
            if not acc:
                foutC.write(rec.format("fasta"))
                continue

            if acc in remove_set:
                black_count += 1
                foutB.write(rec.format("fasta"))
                deleted_accessions.add(acc)
            else:
                sp = acc_to_species.get(acc, "")
                groups_for_acc = acc_to_groups.get(acc, [])

                similar_str = build_similar_annotation_for_acc(
                    acc, sp, groups_for_acc, similar_groups, acc_to_species
                )

                base_header = re.sub(r"\|similar_to=.*", "", rec.description)

                if similar_str:
                    new_header = base_header + "|" + similar_str
                    rec.id = rec.name = rec.description = new_header
                    modified += 1
                else:
                    rec.id = rec.name = rec.description = base_header

                foutC.write(rec.format("fasta"))
                clean_accessions.add(acc)

    remain = total - black_count
    logger.info(f"[FASTA] {infile.name}; Total: {total}; Removed: {black_count}; Retained: {remain}; Modified (similar): {modified}")
    return total, black_count, remain

def process_gb_by_filename(infile: Path, out_cleaned: Path, out_blacklisted: Path,
                           remove_set: set, logger: logging.Logger):
    global deleted_accessions, clean_accessions
    total = 1
    black_count = 0
    base_name = infile.stem
    if "|" in base_name:
        parts = base_name.split("|")
        acc_raw = parts[-1].strip()
        acc_nover = re.sub(r"\.\d+$", "", acc_raw)
        if acc_nover in remove_set:
            black_count += 1
            shutil.copy2(infile, out_blacklisted)
            deleted_accessions.add(acc_nover)
        else:
            shutil.copy2(infile, out_cleaned)
            clean_accessions.add(acc_nover)
    else:
        shutil.copy2(infile, out_cleaned)
    remain = total - black_count
    return total, black_count, remain

def process_nwk_file(infile: Path, out_cleaned_file: Path,
                     remove_set: set, similar_map: dict, species_map: dict,
                     logger: logging.Logger):
    try:
        t = Tree(str(infile), format=1)
    except Exception as e:
        logger.info(f"[SKIPPED NWK] Failed to read file: {infile.name}, Reason: {e}")
        return (0, 0, 0)

    all_leaves = t.get_leaves()
    total_leaves = len(all_leaves)

    keep_names = []
    remove_count = 0
    for leaf in all_leaves:
        acc = extract_accession(leaf.name)
        if acc and acc in remove_set:
            remove_count += 1
        else:
            keep_names.append(leaf.name)

    t.prune(keep_names, preserve_branch_length=True)

    remain = len(t.get_leaves())
    t.write(format=1, outfile=str(out_cleaned_file))
    logger.info(f"[NWK] {infile.name}; Total leaves: {total_leaves}; Removed: {remove_count}; Retained: {remain}")
    return total_leaves, remove_count, remain

def process_single_file(infile: Path,
                        cleaned_dir: Path,
                        blacklisted_dir: Path,
                        remove_set: set,
                        similar_groups: list,
                        species_map: dict,
                        logger: logging.Logger):

    version_tag = get_version_suffix()  
    suffix = infile.suffix.lower()

    # FASTA
    if suffix in {".fa", ".fasta", ".fas", ".faa", ".fna"}:
        out_cleaned_file = cleaned_dir / f"{infile.stem}_cleaned_{version_tag}{infile.suffix}"
        out_blacklisted_file = blacklisted_dir / f"{infile.stem}_blacklisted_{version_tag}{infile.suffix}"
        return process_fasta(infile,
                             out_cleaned_file,
                             out_blacklisted_file,
                             remove_set, similar_groups, logger)

    # GB
    elif suffix in {".gb", ".genbank"}:
        out_cleaned_file = cleaned_dir / f"{infile.stem}_cleaned_{version_tag}{infile.suffix}"
        out_blacklisted_file = blacklisted_dir / f"{infile.stem}_blacklisted_{version_tag}{infile.suffix}"
        return process_gb_by_filename(infile,
                                      out_cleaned_file,
                                      out_blacklisted_file,
                                      remove_set, logger)

    # NWK
    elif suffix == ".nwk":
        out_cleaned_file = cleaned_dir / f"{infile.stem}_cleaned_{version_tag}{infile.suffix}"
        dummy_map = {}  
        return process_nwk_file(infile,
                                out_cleaned_file,
                                remove_set, dummy_map, species_map, logger)
    else:
        logger.info(f"[SKIPPED] Unsupported file type: {infile.name}")
        return 0, 0, 0

def process_folder(folder: Path,
                   base_cleaned: Path,
                   base_blacklisted: Path,
                   remove_set: set,
                   similar_groups: list,
                   species_map: dict,
                   logger: logging.Logger):

    version_tag = get_version_suffix()

    folder_cleaned = base_cleaned / f"{folder.name}_cleaned_{version_tag}"
    folder_cleaned.mkdir(parents=True, exist_ok=True)

    need_blackfolder = False
    relevant_files = []

    for root, dirs, files in os.walk(folder):
        for fname in files:
            p = Path(root) / fname
            suffix = p.suffix.lower()
            if suffix in {".fa", ".fasta", ".fas", ".faa", ".fna", 
                          ".gb", ".genbank", ".nwk"}:
                relevant_files.append(p)
                if suffix in {".fa", ".fasta", ".fas", ".faa", ".fna", ".gb", ".genbank"}:
                    need_blackfolder = True

    if need_blackfolder:
        folder_blacklisted = base_blacklisted / f"{folder.name}_blacklisted_{version_tag}"
        folder_blacklisted.mkdir(parents=True, exist_ok=True)
    else:
        folder_blacklisted = None

    total_count, black_count = 0, 0
    for p in relevant_files:
        suffix = p.suffix.lower()
        if suffix in {".fa", ".fasta", ".fas", ".faa", ".fna"}:
            out_cleaned_file = folder_cleaned / f"{p.stem}_cleaned{p.suffix}"
            if folder_blacklisted:
                out_blacklisted_file = folder_blacklisted / f"{p.stem}_blacklisted{p.suffix}"
            else:
                out_blacklisted_file = base_blacklisted / f"{p.name}_blacklisted{p.suffix}"

            ctot, cblk, _ = process_fasta(
                p, out_cleaned_file, out_blacklisted_file,
                remove_set, similar_groups, logger
            )
        elif suffix in {".gb", ".genbank"}:
            out_cleaned_file = folder_cleaned / p.name
            if folder_blacklisted:
                out_blacklisted_file = folder_blacklisted / p.name
            else:
                out_blacklisted_file = base_blacklisted / p.name

            ctot, cblk, _ = process_gb_by_filename(
                p, out_cleaned_file, out_blacklisted_file,
                remove_set, logger
            )
        elif suffix == ".nwk":
            out_cleaned_file = folder_cleaned / p.name
            dummy_map = {}
            ctot, cblk, _ = process_nwk_file(
                p, out_cleaned_file,
                remove_set, dummy_map, species_map, logger
            )
        else:
            logger.info(f"[SKIPPED] Unsupported file type: {p.name}")
            ctot, cblk, _ = 0, 0, 0

        total_count += ctot
        black_count += cblk

    remain = total_count - black_count
    logger.info(f"[FOLDER] {folder.name}; Total: {total_count}; Removed: {black_count}; Retained: {remain}")
    return total_count, black_count, remain

def parse_args():
    parser = argparse.ArgumentParser(
        description="Remove sequences listed in anomalies.txt from NWK/FA/GB files, and annotate retained sequences listed in similar.txt as 'similar'."
    )
    parser.add_argument(
        "-i", "--input", nargs="+", required=True,
        help="One or more input files or folders. Supports mixed formats (NWK, FA, GB)."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output directory. Will contain outputs/cleaned/ and outputs/blacklisted/ subdirectories."
    )
    parser.add_argument(
        "-a", "--anomalies", nargs="+", required=True,
        help="One or more clustering_anomalies.txt files. All entries will be merged for sequence removal."
    )
    parser.add_argument(
        "-s", "--similar", nargs="+", required=True,
        help="One or more similar.txt files. All entries will be merged to annotate retained sequences."
    )
    return parser.parse_args()

def main():
    global deleted_accessions, clean_accessions
    args = parse_args()

    logger = logging.getLogger("remove_similar")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    out_base = Path(args.output).resolve()
    outputs_dir = out_base / "outputs"
    outputs_dir.mkdir(parents=True, exist_ok=True)

    log_path = outputs_dir / "remove_similar.log"
    fh = logging.FileHandler(log_path, mode='w', encoding='utf-8')
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter("%(message)s")
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    cleaned_dir = outputs_dir / "cleaned"
    cleaned_dir.mkdir(parents=True, exist_ok=True)
    blacklisted_dir = outputs_dir / "blacklisted"
    blacklisted_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=========== remove_with_similar + ETE3 Prune ===========")
    logger.info(f"Input path(s): {args.input}")
    logger.info(f"Output directory: {args.output}")
    logger.info(f"Anomalies file(s): {args.anomalies}")
    logger.info(f"Similar file(s): {args.similar}")
    logger.info("--------------------------------------------------------")

    anomalies_set = load_anomalies_all(args.anomalies)
    species_map = load_species_map_all(args.anomalies)

    similar_groups = load_similar_all(args.similar)
    remove_set = anomalies_set

    logger.info(f"[LOADED] Sequences to remove: {len(remove_set)}; Similarity groups: {len(similar_groups)}")

    total_targets = 0
    for path_str in args.input:
        p = Path(path_str).resolve()
        if p.is_dir():
            process_folder(p, cleaned_dir, blacklisted_dir,
                           remove_set, similar_groups, species_map, logger)
            total_targets += 1
        elif p.is_file():
            process_single_file(
                p, cleaned_dir, blacklisted_dir,
                remove_set, similar_groups, species_map, logger
            )
            total_targets += 1
        else:
            logger.info(f"[SKIPPED] Invalid path: {p}")

    blacklist_path = blacklisted_dir / "blacklist.txt"
    cleanlist_path = cleaned_dir / "cleanlist.txt"
    with open(blacklist_path, "w", encoding="utf-8") as bf:
        for acc in sorted(deleted_accessions):
            bf.write(acc + "\n")
    with open(cleanlist_path, "w", encoding="utf-8") as cf:
        for acc in sorted(clean_accessions):
            cf.write(acc + "\n")

    logger.info(f"Blacklist saved to: {blacklist_path}")
    logger.info(f"Cleaned list saved to: {cleanlist_path}")
    logger.info("======================================================")
    print(f"Processing complete. All results are saved in: {outputs_dir}")

if __name__=="__main__":
    main()