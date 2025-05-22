#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import shutil
import re
import random
from datetime import datetime
from collections import defaultdict
from Bio import SeqIO

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


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Single-pass scanning with priority rules and species-level filtering."
    )

    parser.add_argument(
        '-i', '--input_dir', required=True,
        help='Input directory containing .gb files.'
    )

    parser.add_argument(
        '-o', '--base_output_dir', required=True,
        help='Base output directory. An "outputs" subdirectory will be created inside it.'
    )

    parser.add_argument(
        '-m', '--max-sequences', type=int, default=5,
        help='Maximum number of sequences to retain per species (default: 5).'
    )

    parser.add_argument(
        '-x', action='store_true',
        help='Retain hybrid sequences (do not filter sequences containing " x ").'
    )

    parser.add_argument(
        '-cf', action='store_true',
        help='Retain sequences labeled with "cf." (e.g., "Cyprinus cf. carpio").'
    )

    parser.add_argument(
        '-sp', action='store_true',
        help='Retain sequences labeled with "sp." or "ssp." (e.g., "Cyprinus sp.").'
    )

    parser.add_argument(
        '-s', action='store_true',
        help='Retain subspecies (by default, subspecies are filtered out).'
    )

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument(
        '-c', '--complete-only', action='store_true',
        help='Retain only complete mitochondrial genomes. All partial or fragmentary sequences will be filtered out.'
    )

    group.add_argument(
        '-g', '--gene', nargs='+', metavar=('GENE', 'MINLEN'),
        help=(
            "Target gene name and optional minimum length (in bp).\n"
            "Examples:\n"
            "  -g 12S 200   → retain sequences where 12S ≥ 200 bp\n"
            "  -g 12S       → default minimum length is 1 bp\n"
            "If complete genomes are insufficient, the longest available fragments containing the target gene will be used as fallback."
        )
    )

    args = parser.parse_args()

    args.min_len = 1           
    if args.gene:
        if len(args.gene) not in (1, 2):
            parser.error("The '-g' option requires 1 or 2 arguments: GENE [MINLEN]")

        gene_name = args.gene[0].upper()
        if gene_name not in GENE_PATTERNS:
            parser.error(f"Unknown gene: {gene_name}. Supported genes are: {', '.join(GENE_PATTERNS)}")

        if len(args.gene) == 2:
            try:
                ml = int(args.gene[1])
            except ValueError:
                parser.error("Minimum length must be an integer.")
            if ml < 1:
                parser.error("Minimum length must be ≥ 1.")
            args.min_len = ml
        else:
            args.min_len = 1

        args.gene = gene_name
    else:
        args.min_len = 0          

    return args

def setup_directories(base_output_dir):
    outputs_dir = os.path.join(base_output_dir, 'outputs')
    os.makedirs(outputs_dir, exist_ok=True)
    return outputs_dir

def read_all_gb_files(input_dir):
    return [f for f in os.listdir(input_dir) if f.endswith('.gb')]

def rename_filtered_dir_if_needed(directory_path, base_name, moved_count):
    if directory_path and os.path.exists(directory_path):
        if moved_count > 0:
            parent_dir = os.path.dirname(directory_path)
            new_name = f"{base_name}_{moved_count}"
            new_dir_path = os.path.join(parent_dir, new_name)
            os.rename(directory_path, new_dir_path)
            return new_dir_path
    return directory_path

def is_subspecies(filename):

    name_part = filename.split('|')[0]
    parts = name_part.split('_')
    return (len(parts) > 2)

def get_first_author_simplified(references):

    for ref in references:
        if hasattr(ref, "authors") and ref.authors:
            authors_line = ref.authors.strip()
            if not authors_line:
                continue
            parts = authors_line.split(',')
            if len(parts) > 0:
                return parts[0].strip()
    return None

def collect_sequences_for_species(input_dir, files, filtered_out, args):

    species_sequences = defaultdict(list)
    valid_files = [f for f in files if f.endswith(".gb") and f not in filtered_out]

   
    gene_regex = re.compile(GENE_PATTERNS[args.gene], re.IGNORECASE) if args.gene else None
    

    for filename in valid_files:
        file_path = os.path.join(input_dir, filename)
        try:
            with open(file_path, "r") as gb_file:
                for record in SeqIO.parse(gb_file, "genbank"):
                    species_name = record.annotations.get("organism", "Unknown_species")

                    refs = record.annotations.get("references", [])
                    first_author = get_first_author_simplified(refs)

                    geo_loc = None
                    for feature in record.features:
                        if feature.type == "source":
                            if "geo_loc_name" in feature.qualifiers:
                                geo_loc = feature.qualifiers["geo_loc_name"][0]
                                break

                    acc_list = record.annotations.get("accessions", [])
                    is_nc = any(acc.startswith("NC_") for acc in acc_list)

                    authors_set = set()
                    if first_author:
                        authors_set.add(first_author)

                    definition = record.description.lower()
                    is_full_length = ("complete genome" in definition) or ("complete mitochondrial genome" in definition)
                    seq_len = len(record.seq)

                    target_gene_len = 0
                    if gene_regex:
                        for feature in record.features:
                            if feature.type in ("gene", "CDS", "rRNA"):
                                text = " ".join(feature.qualifiers.get("gene", []) +
                                                feature.qualifiers.get("product", []) +
                                                feature.qualifiers.get("note", []))
                                if re.search(gene_regex, text):
                                    target_gene_len = max(target_gene_len, len(feature.location))
                    if target_gene_len < args.min_len:
                        target_gene_len = 0          

                    species_sequences[species_name].append(
                        (file_path, authors_set, geo_loc, is_nc, seq_len, is_full_length, target_gene_len)
                    )
        except:
            pass

    return species_sequences

def parse_country(geo_loc):

    if not geo_loc:
        return None
    parts = geo_loc.split(':')
    country = parts[0].strip()
    return country

def old_logic_filter(sequences, max_sequences):

    selected_files = []
    operation_notes = []

    from collections import defaultdict
    authors_count = defaultdict(int)
    geo_count = defaultdict(int)

    nc_sequences = [s for s in sequences if s[3]]   # s[3] = is_nc
    non_nc_sequences = [s for s in sequences if not s[3]]

    if nc_sequences:
        selected_files.extend(nc_sequences)
        operation_notes.append(f"NC sequences retained: {len(nc_sequences)}")

    if len(selected_files) >= max_sequences:
        operation_notes.append("Maximum number of NC sequences reached. Skipping further filtering.")
        return [sf[0] for sf in selected_files[:max_sequences]], operation_notes

    remaining_slots = max_sequences - len(selected_files)

    country_map = defaultdict(list)
    for item in non_nc_sequences:
        (f, a, g, i, sl, fl, *_) = item
        c = parse_country(g)
        if c is not None:
            country_map[c].append(item)

    all_countries = list(country_map.keys())
    random.shuffle(all_countries)
    if len(all_countries) > remaining_slots:
        chosen_countries = all_countries[:remaining_slots]
        operation_notes.append(f"Number of countries exceeds {remaining_slots}. Randomly retained {len(chosen_countries)} countries.")
    else:
        chosen_countries = all_countries

    geo_selected = []
    for c in chosen_countries:
        items = country_map[c]
        if not items:
            continue
        random.shuffle(items)
        picked = items[0]
        geo_selected.append(picked)
        geo_count[c] += 1

    if geo_selected:
        loc_count = len(geo_selected)
        detail_str = ", ".join([f"{k}: {v} sequences" for k, v in geo_count.items()])
        operation_notes.append(f"{loc_count} sequences selected after country-level deduplication ({detail_str})")

    selected_files.extend(geo_selected)

    if len(selected_files) >= max_sequences:
        return [sf[0] for sf in selected_files[:max_sequences]], operation_notes

    remaining_slots = max_sequences - len(selected_files)
    already_paths = set(x[0] for x in selected_files)
    leftover_for_author = [s for s in non_nc_sequences if s[0] not in already_paths]

    author_files_map = defaultdict(list)
    for (fpath, authors_set, geo_loc, is_nc, sl, fl, *_) in leftover_for_author:
        for au in authors_set:
            author_files_map[au].append((fpath, authors_set, geo_loc, is_nc, sl, fl))

    all_authors = list(author_files_map.keys())
    random.shuffle(all_authors)
    if len(all_authors) > remaining_slots:
        chosen_authors = all_authors[:remaining_slots]
        operation_notes.append(f"Number of authors exceeds {remaining_slots}. Randomly retained {len(chosen_authors)} authors.")
    else:
        chosen_authors = all_authors

    author_selected_files = []
    for au in chosen_authors:
        items = author_files_map[au]
        random.shuffle(items)
        picked = items[0]
        author_selected_files.append(picked)
        authors_count[au] += 1

    if author_selected_files:
        selected_num = len(author_selected_files)
        author_detail_str = ", ".join([f"{k}:{v}sequences" for k,v in authors_count.items()])
        operation_notes.append(f"{selected_num} sequences selected after author-level deduplication ({author_detail_str})")

    selected_files.extend(author_selected_files)

    if len(selected_files) >= max_sequences:
        return [sf[0] for sf in selected_files[:max_sequences]], operation_notes

    remaining_slots = max_sequences - len(selected_files)
    already_paths = set(x[0] for x in selected_files)
    leftover = [s for s in sequences if s[0] not in already_paths]
    if leftover:
        random.shuffle(leftover)
        needed = remaining_slots
        picked = leftover[:needed]
        selected_files.extend(picked)
        operation_notes.append(f"Randomly supplemented {min(needed, len(picked))} sequences.")

    return [sf[0] for sf in selected_files[:max_sequences]], operation_notes


def filter_sequences_for_species(species, sequences, max_sequences, args):

    total_count = len(sequences)
    nc_cnt = sum(1 for s in sequences if s[3])

    if total_count <= max_sequences:
        notes = [f"Total sequences for this species: {total_count} (≤ {max_sequences}). All retained."]
        if nc_cnt:
            notes.append(f"NC sequences retained: {nc_cnt}")
        return [x[0] for x in sequences], notes

# Separate full-length genomes from fragments:
    full_with_gene    = [s for s in sequences if s[5] and s[6] > 0]
    full_without_gene = [s for s in sequences if s[5] and s[6] == 0]
    partial_list      = [s for s in sequences if not s[5]]

    # ① Full sequences containing the target gene are sorted by gene length (descending)
    full_with_gene.sort(key=lambda x: x[6], reverse=True)
    full_length_list = full_with_gene + full_without_gene

    # ② Filtering is applied only to fragment sequences
    if args.gene:
        partial_list = [s for s in partial_list if s[6] > 0]

    # ③ If no sequences remain after filtering → discard the entire species
    if not full_length_list and not partial_list:
        return [], [f"No fragments contain {args.gene} ≥ {args.min_len} bp. Species discarded."]

    # ④ Then merge the retained full sequences and filtered fragments
    sequences = full_length_list + partial_list

    # If no sequences remain after merging → discard the species entirely
    if not sequences:
        msg = f"Only fragments present, but none contain {args.gene} ≥ {args.min_len} bp. Species discarded."
        return [], [msg]

    if len(full_length_list) > max_sequences:
        sel_files, notes = old_logic_filter(full_length_list, max_sequences)
        notes.append("Number of complete genomes exceeds max_sequences. Fragments will not be considered.")
        return sel_files, notes

    kept_full, notes_full = old_logic_filter(full_length_list, max_sequences)
    kept_full_count = len(kept_full)
    remain = max_sequences - kept_full_count

    if remain <= 0:
        return kept_full, notes_full

    partial_sorted = sorted(partial_list, key=lambda x: x[6], reverse=True)
    chosen_partials = partial_sorted[:remain]

    notes_partial = []
    if partial_list:
        notes_partial.append(
            f"{kept_full_count} complete genomes retained; {len(partial_list)} fragments found; "
            f"{len(chosen_partials)} fragments selected based on target gene length for supplementation."
        )
    final_files = kept_full + [x[0] for x in chosen_partials]
    final_notes = notes_full + notes_partial
    return final_files, final_notes

def run_filter_sequences(species_sequences, max_sequences, outputs_dir, args):
    species_count = len(species_sequences)
    details, output_files = [], []

    for species, seqs in sorted(species_sequences.items()):
        selected_files, notes = filter_sequences_for_species(
            species, seqs, max_sequences, args
        )
        if not any("NC sequences retained" in n for n in notes):
            notes.append("No NC sequences found.")

        output_files.extend(selected_files)
        details.append((species, len(seqs), notes))

    final_dir = os.path.join(outputs_dir, f"final_sequences_gb")
    os.makedirs(final_dir, exist_ok=True)

    copied_set = set()
    for fp in output_files:
        if fp not in copied_set:
            shutil.copy(fp, final_dir)
            copied_set.add(fp)

    if args.gene:
        part_dir = os.path.join(outputs_dir, "filtered_out", "7.unselected_partials_gene")
        os.makedirs(part_dir, exist_ok=True)

        selected_set = set(output_files)
        for seqs in species_sequences.values():
            for fp, *_ , is_full, _tglen in seqs:      
                if (not is_full) and fp not in selected_set:
                    shutil.copy(fp, part_dir)

        cnt = len(os.listdir(part_dir))
        if cnt == 0:
            os.rmdir(part_dir)
        else:
            new_dir = f"{part_dir}_{cnt}"
            suffix = 1
            while os.path.exists(new_dir):
                new_dir = f"{part_dir}_{cnt}_{suffix}"
                suffix += 1
            os.rename(part_dir, new_dir)

    return {
        "species_count": species_count,
        "final_count": len(output_files),
        "details": details,
        "final_dir": final_dir,
    }


def write_final_log(
    log_path,
    args,
    all_files_count,
    step0_info,
    step1_info,
    step2_info,
    step3_info,
    step4_info,
    step5_info,
    step6_info,
    step7_info
):
    start_time = step1_info.get('start_time').strftime("%Y-%m-%d %H:%M:%S")
    end_time   = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    cf_status = "retained" if args.cf else "filtered"
    sp_status = "retained" if args.sp else "filtered"
    hybrid_status = "retained" if args.x else "filtered"
    subspecies_status = "retained" if args.s else "filtered"
    complete_only_status = "yes" if args.complete_only else "no"

    if args.gene:
        gene_info = f"{args.gene} (≥ {args.min_len} bp)"
    else:
        gene_info = "not specified"


    with open(log_path, 'w') as log:
        log.write("=====================================\n")
        log.write(" Species-Level Filtering Log\n")
        log.write("=====================================\n")
        log.write(f"Start time: {start_time}\n")
        log.write(f"End time: {end_time}\n\n")

        log.write("Parameters:\n")
        log.write(f"  Input directory: {args.input_dir}\n")
        log.write(f"  Output directory: {args.base_output_dir}\n")
        log.write(f"  Max sequences per species: {args.max_sequences}\n")
        log.write(f"  Retain subspecies: {subspecies_status}\n")
        log.write(f"  Retain 'cf.': {cf_status}\n")
        log.write(f"  Retain 'sp.': {sp_status}\n")
        log.write(f"  Retain hybrids: {hybrid_status}\n")
        log.write(f"  Complete genomes only (-c): {complete_only_status}\n")
        log.write(f"  Gene filtering (-g): {gene_info}\n")
        log.write(f"  Total input files: {all_files_count}\n\n")

        # Step 0
        log.write("------------------------------------------------------------\n")
        log.write("Step 0: Empty Sequences\n")
        log.write("------------------------------------------------------------\n")
        moved = step0_info['moved']
        log.write(f"  Files moved/filtered: {moved}\n\n")
        if moved == 0:
            log.write("  No files were moved.\n\n")
        else:
            for fn in step0_info['details']:
                log.write(f"  {fn}\n")
            log.write("\n")


        # Step 1
        log.write("------------------------------------------------------------\n")
        log.write("Step 1: Incomplete Mitochondrial Genomes\n")
        log.write("------------------------------------------------------------\n")
        moved = step1_info['moved']
        log.write(f"  Files moved/filtered: {moved}\n\n")

        if moved == 0:
            log.write("  No files were moved.\n")
            if args.complete_only:
                # In complete-only mode, but no incomplete fragments were found
                log.write("  (-c mode: all sequences are complete mitochondrial genomes)\n")
            else:
                # In gene-based mode, incomplete fragments are not moved
                log.write("  (-g mode: this step only moves incomplete fragments in -c mode)\n")
            log.write("\n")
        else:
            for fn in step1_info['details']:
                log.write(f"  {fn}\n")
            log.write("\n")


        # Step 2
        log.write("------------------------------------------------------------\n")
        log.write("Step 2: NC Duplicates\n")
        log.write("------------------------------------------------------------\n")
        moved = step2_info['moved']
        log.write(f"  Files moved/filtered: {moved}\n\n")

        if moved == 0:
            log.write("  No files were moved.\n\n")
        else:
            for fn in step2_info['details']:
                log.write(f"  {fn}\n")
            log.write("\n")


        # Step 3
        log.write("------------------------------------------------------------\n")
        log.write("Step 3: Hybrid Sequences\n")
        log.write("------------------------------------------------------------\n")
        moved = step3_info['moved']
        log.write(f"  Files moved/filtered: {moved}\n\n")

        if moved == 0:
            log.write("  No files were moved.\n\n")
        else:
            for fn in step3_info['details']:
                log.write(f"  {fn}\n")
            log.write("\n")


        # Step 4
        log.write("------------------------------------------------------------\n")
        log.write("Step 4: 'cf.' Sequences\n")
        log.write("------------------------------------------------------------\n")
        moved = step4_info['moved']
        log.write(f"  Files moved/filtered: {moved}\n\n")

        if moved == 0:
            log.write("  No files were moved.\n\n")
        else:
            for fn in step4_info['details']:
                log.write(f"  {fn}\n")
            log.write("\n")


        # Step 5
        log.write("------------------------------------------------------------\n")
        log.write("Step 5: 'sp.' / 'ssp.' Sequences\n")
        log.write("------------------------------------------------------------\n")
        moved = step5_info['moved']
        log.write(f"  Files moved/filtered: {moved}\n\n")

        if moved == 0:
            log.write("  No files were moved.\n\n")
        else:
            for fn in step5_info['details']:
                log.write(f"  {fn}\n")
            log.write("\n")


        # Step 6
        log.write("------------------------------------------------------------\n")
        log.write("Step 6: Subspecies\n")
        log.write("------------------------------------------------------------\n")
        moved = step6_info['moved']
        log.write(f"  Files moved/filtered: {moved}\n\n")

        if moved == 0:
            log.write("  No files were moved.\n\n")
        else:
            for fn in step6_info['details']:
                log.write(f"  {fn}\n")
            log.write("\n")


        # Step 7
        log.write("------------------------------------------------------------\n")
        log.write("Step 7: Species-Level Filtering\n")
        log.write("------------------------------------------------------------\n")
        log.write(f"  Total species: {step7_info['species_count']}\n")
        log.write(f"  Final sequences retained: {step7_info['final_count']}\n")
        log.write(f"  Output directory: {step7_info['final_dir']}\n\n")

        for species, seq_count, notes in step7_info['details']:
            log.write(f"  Species: {species}\n")
            log.write(f"    Original sequence count: {seq_count}\n")
            if notes:
                log.write("    Filtering steps:\n")
                for n in notes:
                    log.write(f"      * {n}\n")
            else:
                log.write("    No additional filtering applied.\n")
            log.write("\n")


        log.write("------------------------------------------------------------\n")
        log.write("Summary:\n")
        log.write(f"- Total input files: {all_files_count}\n")
        log.write(f"- Step 0 (Empty Sequences):     {step0_info['moved']} files moved\n")
        log.write(f"- Step 1 (Incomplete Mitochondrial Genomes):  {step1_info['moved']} files moved\n")
        log.write(f"- Step 2 (NC Duplicates):       {step2_info['moved']} files moved\n")
        log.write(f"- Step 3 (Hybrids):             {step3_info['moved']} files moved\n")
        log.write(f"- Step 4 (cf. Sequences):       {step4_info['moved']} files moved\n")
        log.write(f"- Step 5 (sp./ssp. Sequences):  {step5_info['moved']} files moved\n")
        log.write(f"- Step 6 (Subspecies):          {step6_info['moved']} files moved\n")

        # -------- Step 7 Fragment Handling --------
        if args.complete_only:
            log.write("- Step 7 Fragment Handling: All fragments were filtered in Step 1 (-c mode)\n")
        elif args.gene:
            log.write(f"- Step 7 Fragment Handling: Fragments were ranked by {args.gene.upper()} length. "
                    "Unselected fragments saved to: 7.unselected_partials_gene/\n")
        else:
            log.write("- Step 7 Fragment Handling: No fragment processing performed\n")
        # ------------------------------------------

        log.write(f"- Final sequences retained: {step7_info['final_count']}\n")
        log.write(f"- Final output directory: {step7_info['final_dir']}\n")
        log.write("------------------------------------------------------------\n")
        log.write("[End of Report]\n")


def single_pass_priority_filter(args, input_dir, all_files, outputs_dir):
    """
    Core Workflow:

    1) For each file, evaluate the following flags:
    - is_incomplete, is_hybrid, is_subspecies, is_cf, is_sp, is_nc_duplicate

    2) Apply filtering in the following priority order:
    incomplete > nc_duplicates > hybrid > cf > sp > subspecies > keep

    3) Move filtered files to their corresponding directories (renamed_xxx)

    4) Keep retained files untouched, then perform species-level filtering;
    copy final selections to 'final_sequences_gb/'

    5) Summarize and return the results from all 7 filtering steps for logging
    """
    file_info_list = []
    for fn in all_files:
        filepath = os.path.join(input_dir, fn)
        info = {
            'filename': fn,
            'filepath': filepath,
            'is_empty_seq': False,
            'is_incomplete': False,
            'is_hybrid': False,
            'is_subspecies': False,
            'is_cf': False,
            'is_sp': False,
            'is_nc_duplicate': False,
            'label': 'keep'
        }
        try:
            with open(filepath, "r") as gbfile:
                record = SeqIO.read(gbfile, "genbank")
                definition = record.description.lower()

                if len(record.seq) == 0 or record.annotations.get("contig"):
                    info['is_empty_seq'] = True
                    info['label'] = '0.empty_sequence'
                    file_info_list.append(info)
                    continue        

                is_full_length = ("complete genome" in definition) or ("complete mitochondrial genome" in definition)

                if (not is_full_length) and args.complete_only:
                    info['is_incomplete'] = True
                    info['label'] = '1.incomplete_genomes'
                else:
                    info['is_incomplete'] = False

                organism = record.annotations.get("organism", "")
                if " x " in organism:
                    info['is_hybrid'] = True

                if " cf." in organism:
                    info['is_cf'] = True

                if (" sp." in organism) or (" ssp." in organism):
                    info['is_sp'] = True

                if is_subspecies(fn):
                    info['is_subspecies'] = True

        except:
            pass
        file_info_list.append(info)

    ref_pattern = re.compile(r"reference sequence (?:is identical to|was derived from) ([\w\.]+)")
    file_contents = {}
    for f in file_info_list:
        try:
            with open(f['filepath'], 'r') as r:
                file_contents[f['filename']] = r.read()
        except:
            file_contents[f['filename']] = ""

    nc_names = set(fn for fn in all_files if 'NC_' in fn)
    identified_dup = {}

    for f in file_info_list:
        if f['filename'] in nc_names:
            txt = file_contents[f['filename']]
            matches = set(ref_pattern.findall(txt))
            if matches:
                matched_files = []
                for m in matches:
                    for o in file_info_list:
                        if o['filename'] not in nc_names:
                            if m in o['filename']:
                                matched_files.append(o['filename'])
                if matched_files:
                    identified_dup[f['filename']] = matched_files

    for nc, orig_list in identified_dup.items():
        for ofn in orig_list:
            for i in file_info_list:
                if i['filename'] == ofn:
                    i['is_nc_duplicate'] = True

    for f in file_info_list:
        if f['label'] != 'keep':
            continue
        if f['is_incomplete']:
            f['label'] = '1.incomplete_genomes'
            continue
        if f['is_nc_duplicate']:
            f['label'] = '2.nc_duplicates'
            continue
        if f['is_hybrid'] and not args.x:
            f['label'] = '3.hybrid'
            continue
        if f['is_cf'] and not args.cf:
            f['label'] = '4.cf'
            continue
        if f['is_sp'] and not args.sp:
            f['label'] = '5.sp'
            continue
        if f['is_subspecies'] and not args.s:
            f['label'] = '6.subspecies'
            continue

    filtered_out_dir = os.path.join(outputs_dir, "filtered_out")
    os.makedirs(filtered_out_dir, exist_ok=True)
    label_to_dir = {
    '0.empty_sequence': os.path.join(filtered_out_dir, '0.empty_sequence'),  
    '1.incomplete_genomes': os.path.join(filtered_out_dir, '1.incomplete_genomes'),
    '2.nc_duplicates': os.path.join(filtered_out_dir, '2.nc_duplicates'),
    '3.hybrid':        os.path.join(filtered_out_dir, '3.hybrid'),
    '4.cf':            os.path.join(filtered_out_dir, '4.cf'),
    '5.sp':            os.path.join(filtered_out_dir, '5.sp'),
    '6.subspecies':    os.path.join(filtered_out_dir, '6.subspecies'),
    'keep':            None
    }


    for d in label_to_dir.values():
        if d:
            os.makedirs(d, exist_ok=True)

    for f in file_info_list:
        lb = f['label']
        if lb == 'keep':
            continue            
        dest = label_to_dir.get(lb)
        if dest is None:        
           continue          
        shutil.copy(f['filepath'], dest)

    def rename_dir(label):
        dd = label_to_dir[label]
        if not dd or not os.path.exists(dd):
            return
        c = len(os.listdir(dd))          
        new_path = f"{dd}_{c}"
        suffix = 1
        while os.path.exists(new_path):
            new_path = f"{dd}_{c}_{suffix}"
            suffix += 1
        os.rename(dd, new_path)
        label_to_dir[label] = new_path

    rename_dir('0.empty_sequence')
    rename_dir('1.incomplete_genomes')
    rename_dir('2.nc_duplicates')
    rename_dir('3.hybrid')
    rename_dir('4.cf')
    rename_dir('5.sp')
    rename_dir('6.subspecies')

    keep_files = [x['filename'] for x in file_info_list if x['label'] == 'keep']
    filtered_out_files = set(x['filename'] for x in file_info_list if x['label'] != 'keep')

    sp_seq = collect_sequences_for_species(input_dir, all_files, filtered_out_files, args)
    step7 = run_filter_sequences(sp_seq, args.max_sequences, outputs_dir, args)

    empty_list = [x['filename'] for x in file_info_list if x['label'] == '0.empty_sequence']
    step0_info = {'moved': len(empty_list), 'details': empty_list}
    # -----------------------------------
    incomplete_list = [x['filename'] for x in file_info_list if x['label'] == '1.incomplete_genomes']
    step1_info = {
        'moved': len(incomplete_list),
        'details': incomplete_list,
        'start_time': datetime.now()
    }
    ncdup_list = [x['filename'] for x in file_info_list if x['label'] == '2.nc_duplicates']
    step2_info = {'moved': len(ncdup_list), 'details': ncdup_list}
    hybrid_list = [x['filename'] for x in file_info_list if x['label'] == '3.hybrid']
    step3_info = {'moved': len(hybrid_list), 'details': hybrid_list}
    cf_list = [x['filename'] for x in file_info_list if x['label'] == '4.cf']
    step4_info = {'moved': len(cf_list), 'details': cf_list}
    sp_list = [x['filename'] for x in file_info_list if x['label'] == '5.sp']
    step5_info = {'moved': len(sp_list), 'details': sp_list}
    sub_list = [x['filename'] for x in file_info_list if x['label'] == '6.subspecies']
    step6_info = {'moved': len(sub_list), 'details': sub_list}
    step7_info = step7

    return step0_info, step1_info, step2_info, step3_info, step4_info, step5_info, step6_info, step7_info


def main():
    args = parse_arguments()
    start_time = datetime.now()

    outputs_dir = setup_directories(args.base_output_dir)
    all_files = read_all_gb_files(args.input_dir)
    all_files_count = len(all_files)

    (step0_info, step1_info, step2_info, step3_info, 
     step4_info, step5_info, step6_info, step7_info) = single_pass_priority_filter(
        args, args.input_dir, all_files, outputs_dir
    )

    step1_info['start_time'] = start_time

    log_path = os.path.join(outputs_dir, "process.log")
    write_final_log(
        log_path, args, all_files_count,
        step0_info,
        step1_info,
        step2_info,
        step3_info,
        step4_info,
        step5_info,
        step6_info,
        step7_info
    )

    print(f"[INFO] Processing complete.")
    print(f"[INFO] Log saved to: {log_path}")


if __name__ == "__main__":
    main()