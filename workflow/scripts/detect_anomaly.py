#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Clustering Anomaly Detection and Tree Visualization Tool

Usage:
  python clustering_anomaly_detector.py -i 1.nwk 2.nwk -o <output_directory> \
       [-og <outgroup1, outgroup2> <outgroup3, outgroup4> ...] \
       [-s <similar_file1> <similar_file2> ...]

Example: 
  python clustering_anomaly_detector.py -i tree1.nwk tree2.nwk \
       -o /path/to/output \
       -og NC_023455,NC_035057 NC_008106,NC_014177 \
       -s similar1.txt similar2.txt
"""

import os
import sys
import argparse
import logging
from datetime import datetime
from collections import defaultdict, Counter
import re
from Bio import Phylo
os.environ['QT_QPA_PLATFORM'] = 'offscreen'
from ete3 import Tree, TreeStyle, NodeStyle, TextFace  # type: ignore


def setup_logging(log_file_path):
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(log_file_path)
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger


def save_multiple_nwk_files(trees, output_dir):
    nwk_dir = os.path.join(output_dir, 'nwk_files')
    os.makedirs(nwk_dir, exist_ok=True)
    for i, tree in enumerate(trees, start=1):
        output_file = os.path.join(nwk_dir, f"tree_{i}.nwk")
        try:
            Phylo.write(tree, output_file, 'newick')
            logging.info(f"Saved tree {i} to {output_file}")
        except Exception as e:
            logging.error(f"Failed to save tree {i} to NWK file: {e}")
            sys.exit(1)


def assign_parents(tree):
    def assign_parent_recursive(clade, parent=None):
        clade.parent = parent
        for child in clade.clades:
            assign_parent_recursive(child, clade)
    assign_parent_recursive(tree.root)


def parse_taxonomic_info(leaf_name):
    parts = leaf_name.split('|')
    if len(parts) < 6:
        logging.warning(f"Invalid leaf node format: {leaf_name}. Taxonomic fields set to 'NA'.")
        return leaf_name, {
            'accession': 'NA',
            'order': 'NA',
            'family': 'NA',
            'genus': 'NA'
        }
    species = parts[0]
    accession = parts[1]
    order = parts[3][2:] if parts[3].startswith("o_") else parts[3]
    family = parts[4][2:] if parts[4].startswith("f_") else parts[4]
    genus = parts[5][2:] if parts[5].startswith("g_") else parts[5]
    return species, {'accession': accession, 'order': order, 'family': family, 'genus': genus}


def find_genus_parent(leaf, taxonomic_info):
    current = leaf
    sp = leaf.name.split('|')[0]
    genus = taxonomic_info.get(sp, {}).get('genus')
    if not genus:
        return None
    while hasattr(current, 'parent') and current.parent:
        parent = current.parent
        child_genera = {taxonomic_info.get(child.name.split('|')[0], {}).get('genus')
                        for child in parent.get_terminals() if child.name}
        if len(child_genera) == 1 and any(
            taxonomic_info.get(child.name.split('|')[0], {}).get('genus') == genus
            for child in parent.get_terminals() if child.name and child != leaf
        ):
            return parent
        current = parent
    return None


def find_family_parent(leaf, taxonomic_info):
    current = leaf
    sp = leaf.name.split('|')[0]
    family = taxonomic_info.get(sp, {}).get('family')
    if not family:
        return None
    while hasattr(current, 'parent') and current.parent:
        parent = current.parent
        child_families = {taxonomic_info.get(child.name.split('|')[0], {}).get('family')
                          for child in parent.get_terminals() if child.name}
        if len(child_families) == 1 and any(
            taxonomic_info.get(child.name.split('|')[0], {}).get('family') == family
            for child in parent.get_terminals() if child.name and child != leaf
        ):
            return parent
        current = parent
    return None


def find_order_parent(leaf, taxonomic_info):
    current = leaf
    sp = leaf.name.split('|')[0]
    order = taxonomic_info.get(sp, {}).get('order')
    if not order:
        return None
    while hasattr(current, 'parent') and current.parent:
        parent = current.parent
        child_orders = {taxonomic_info.get(child.name.split('|')[0], {}).get('order')
                        for child in parent.get_terminals() if child.name}
        if len(child_orders) == 1 and any(
            taxonomic_info.get(child.name.split('|')[0], {}).get('order') == order
            for child in parent.get_terminals() if child.name and child != leaf
        ):
            return parent
        current = parent
    return None


def check_taxonomic_level(leaf1, leaf2, taxonomic_info, level='genus'):
    if level == 'genus':
        p1 = find_genus_parent(leaf1, taxonomic_info)
        p2 = find_genus_parent(leaf2, taxonomic_info)
    elif level == 'family':
        p1 = find_family_parent(leaf1, taxonomic_info)
        p2 = find_family_parent(leaf2, taxonomic_info)
    elif level == 'order':
        p1 = find_order_parent(leaf1, taxonomic_info)
        p2 = find_order_parent(leaf2, taxonomic_info)
    else:
        return None
    if p1 and p2:
        return [leaf1.name, leaf2.name]
    elif p1 and not p2:
        return [leaf2.name]
    elif not p1 and p2:
        return [leaf1.name]
    else:
        return None


def analyze_clustering(tree, outgroup_names):
    leaves_by_species = defaultdict(list)
    taxonomic_info = {}
    family_counts = defaultdict(int)
    anomalies = defaultdict(list)
    classified = set()

    for leaf in tree.get_terminals():
        if leaf.name and leaf.name not in outgroup_names:
            sp, tax_info = parse_taxonomic_info(leaf.name)
            taxonomic_info[sp] = tax_info
            leaves_by_species[sp].append(leaf)
            fam = tax_info['family']
            if fam != 'NA':
                family_counts[fam] += 1

    for sp, lf_list in leaves_by_species.items():
        if len(lf_list) == 1:
            leaf = lf_list[0]
            fam = taxonomic_info[sp]['family']
            odr = taxonomic_info[sp]['order']
            if fam == 'NA' and odr != 'NA':
                if not find_order_parent(leaf, taxonomic_info):
                    anomalies['single'].append(leaf.name)
                    classified.add(leaf.name)
            elif fam != 'NA':
                if family_counts[fam] > 1:
                    if not find_family_parent(leaf, taxonomic_info):
                        anomalies['single'].append(leaf.name)
                        classified.add(leaf.name)
                elif family_counts[fam] == 1 and odr != 'NA':
                    if not find_order_parent(leaf, taxonomic_info):
                        anomalies['single'].append(leaf.name)
                        classified.add(leaf.name)
                else:
                    if not find_family_parent(leaf, taxonomic_info):
                        anomalies['single'].append(leaf.name)
                        classified.add(leaf.name)
        elif len(lf_list) == 2:
            if all(leaf.name not in classified for leaf in lf_list):
                if lf_list[0].parent != lf_list[1].parent:
                    anomalies['two'].extend(leaf.name for leaf in lf_list)
                    classified.update(leaf.name for leaf in lf_list)
        elif len(lf_list) > 2:
            for leaf in lf_list:
                if leaf.name not in classified:
                    parent = leaf.parent
                    sibs = [s for s in parent.get_terminals() if s != leaf and s.name and s.name not in outgroup_names]
                    c = Counter(sib.name.split('|')[0] for sib in sibs if sib.name)
                    total = sum(c.values())
                    sp_count = c.get(sp, 0)
                    if total > 0 and (sp_count / total) < 0.5:
                        anomalies['three_or_more'].append(leaf.name)
                        classified.add(leaf.name)

    blue_seqs = anomalies.get('two', [])
    new_blue = []
    species_in_blue = set(seq.split('|')[0] for seq in blue_seqs)
    for spc in species_in_blue:
        lf_list = [x for x in leaves_by_species[spc] if x.name in blue_seqs]
        if len(lf_list) != 2:
            continue
        leaf1, leaf2 = lf_list
        need_blue = check_taxonomic_level(leaf1, leaf2, taxonomic_info, 'genus')
        if need_blue is None:
            need_blue = check_taxonomic_level(leaf1, leaf2, taxonomic_info, 'family')
            if need_blue is None:
                need_blue = check_taxonomic_level(leaf1, leaf2, taxonomic_info, 'order')
                if need_blue is None:
                    need_blue = [leaf1.name, leaf2.name]
        new_blue.extend(need_blue)
    anomalies['two'] = new_blue
    return anomalies, taxonomic_info


def save_anomalies(anomalies, output_file):
    categories = {'single': 'Green', 'two': 'Blue', 'three_or_more': 'Red'}
    try:
        with open(output_file, 'w') as f:
            for cat, color in categories.items():
                f.write(f"{color} anomalies for {cat} sequences:\n")
                arr = anomalies.get(cat, [])
                if arr:
                    for e in sorted(set(arr)):
                        f.write(f"{e}\n")
                else:
                    f.write("None\n")
                f.write("\n")
        logging.info("Clustering anomalies successfully saved to file.")
    except Exception as e:
        logging.error(f"Failed to save clustering anomalies to file: {e}")
        sys.exit(1)


def load_anomalies(anomalies_file):
    an = {'single': [], 'two': [], 'three_or_more': []}
    try:
        with open(anomalies_file, 'r') as f:
            curcat = None
            for line in f:
                line = line.strip()
                if line.startswith("Green anomalies for single sequences:"):
                    curcat = 'single'
                elif line.startswith("Blue anomalies for two sequences:"):
                    curcat = 'two'
                elif line.startswith("Red anomalies for three_or_more sequences:"):
                    curcat = 'three_or_more'
                elif line and curcat:
                    if line != "None":
                        an[curcat].append(line)
        logging.debug("Anomaly file successfully loaded.")
    except Exception as e:
        logging.error(f"Failed to read anomaly file: {e}")
        sys.exit(1)
    return an


def load_tree_ete3(file_path):
    try:
        tree = Tree(file_path, format=1)
        logging.info("ETE3 tree file successfully loaded.")
        return tree
    except Exception as e:
        logging.error(f"Failed to read Newick file using ETE3: {e}")
        sys.exit(1)


def set_outgroups(tree, outgroup_ids):
    outgroup_nodes = []
    outgroup_names = set()
    try:
        found_ogs_log = []

        for og in outgroup_ids:
            nodes = [node for node in tree.iter_leaves() if og in node.orig_name]
            if nodes:
                outgroup_nodes.append(nodes[0])
                outgroup_names.add(nodes[0].name)
                sp, _ = parse_taxonomic_info(nodes[0].orig_name)
                found_ogs_log.append((og, sp))
            else:
                logging.warning(f"Outgroup ID not found in tree: {og}")

        if found_ogs_log:
            logging.info("Outgroup node(s) found:")
            for og, sp in found_ogs_log:
                logging.info(f"  - {og} ({sp})")

        if outgroup_nodes:
            if len(outgroup_nodes) == 1:
                tree.set_outgroup(outgroup_nodes[0])
            else:
                mrca_node = tree.get_common_ancestor(outgroup_nodes)
                tree.set_outgroup(mrca_node)
            logging.info("Root node successfully set.")
        else:
            logging.warning("No valid outgroup node found.")

    except Exception as e:
        logging.error(f"Error while setting outgroup: {e}")
        sys.exit(1)

    return outgroup_names


def color_tree(tree, anomalies, output_pdf, outgroup_names):
    colors = {'single': 'green', 'two': 'blue', 'three_or_more': 'red'}

    def custom_layout(node):
        node_color = 'black'
        if node.name in outgroup_names:
            node_color = 'black'
        else:
            for cat, nodes in anomalies.items():
                if node.name in nodes:
                    node_color = colors.get(cat, 'black')
                    break
        nstyle = NodeStyle()
        nstyle["fgcolor"] = node_color
        nstyle["size"] = 6
        node.set_style(nstyle)
        if node.is_leaf():
            face = TextFace(" " + node.name, fgcolor=node_color)
            node.add_face(face, column=0, position="branch-right")

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.layout_fn = custom_layout
    ts.scale = 1000
    ts.branch_vertical_margin = 10
    ts.margin_left = 20
    ts.margin_right = 20
    try:
        tree.render(output_pdf, tree_style=ts)
        logging.info(f"Tree coloring complete. PDF saved to: {output_pdf}")
    except Exception as e:
        logging.error(f"Failed to render tree to PDF: {e}")
        sys.exit(1)


def load_user_similar_file(similar_file):
    groups = []
    with open(similar_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                arr = line.split()  
                groups.append(arr)
    return groups


def parse_species(seq_name):
    clean_name = seq_name.split("|similar_to=")[0]
    return clean_name.split('|')[0]


def extract_accession(seq_name):
    parts = seq_name.split('|')
    if len(parts) >= 2:
        return parts[1]
    return seq_name


def apply_similar_to_inplace(node, new_name, anomalies):
    old_name = node.name
    node.name = new_name
    for cat in anomalies:
        if old_name in anomalies[cat]:
            anomalies[cat].remove(old_name)
            anomalies[cat].append(new_name)


def find_representative(seq_list):
    nc_candidates = [x for x in seq_list if 'NC_' in x]
    if nc_candidates:
        return sorted(nc_candidates)[0]
    else:
        return sorted(seq_list)[0]


def annotate_user_defined_similar(similar_file, tree_ete, anomalies):
    groups = load_user_similar_file(similar_file)
    for line_seqs in groups:
        unique_accessions = sorted(set(line_seqs))
        species_group = defaultdict(list)
        for acc in unique_accessions:
            for leaf in tree_ete.iter_leaves():
                leaf_acc = extract_accession(leaf.orig_name)
                if leaf_acc == acc:
                    species = leaf.orig_name.split('|')[0]
                    species_group[species].append(leaf_acc)
                    break
        species_representative = {}
        for sp, accs in species_group.items():
            species_representative[sp] = find_representative(accs)
        for acc in unique_accessions:
            cur_species = None
            for leaf in tree_ete.iter_leaves():
                if extract_accession(leaf.orig_name) == acc:
                    cur_species = leaf.orig_name.split('|')[0]
                    break
            if cur_species is None:
                logging.warning(f"[similar] Accession not found in tree: {acc}. Skipped.")
                continue
            other_reps = [species_representative[sp] for sp in species_representative if sp != cur_species]
            if other_reps:
                sim_str = ",".join(sorted(set(other_reps)))
            else:
                sim_str = "None"
            nds = []
            for leaf in tree_ete.iter_leaves():
                if extract_accession(leaf.orig_name) == acc:
                    nds.append(leaf)
            if not nds:
                logging.warning(f"[similar] Accession {acc} not found in tree. Skipped.")
                continue
            node = nds[0]
            old_name_clean = re.sub(r"\|similar_to=.*", "", node.name)
            new_name = old_name_clean + f"|similar_to={sim_str}"
            apply_similar_to_inplace(node, new_name, anomalies)


def main():
    parser = argparse.ArgumentParser(
        description="Clustering anomaly detection and phylogenetic tree visualization tool."
    )
    parser.add_argument(
        '-i', '--input', nargs='+', required=True,
        help="One or more Newick tree file paths, separated by spaces."
    )
    parser.add_argument(
        '-o', '--output_prefix', required=True,
        help="Output filename prefix. The script will generate files such as *_anomaly.txt and *.pdf."
    )
    parser.add_argument(
        '-og', '--outgroups', nargs='*',
        help=(
            "Outgroup definitions (optional). Three usage patterns are supported:\n"
            "1) If omitted, no outgroup will be set;\n"
            "2) If one group is provided, it will be used for all trees;\n"
            "3) If the number of groups matches the number of input trees, each group will apply to its corresponding tree.\n"
            "Within a group, separate multiple outgroup IDs with commas (no spaces). Use spaces to separate different groups."
        )
    )
    parser.add_argument(
        '-s', nargs='*',
        help="Optional: One or more user-provided similar.txt files. If provided, matching sequences will be annotated with '|similar_to=...'. "
             "If multiple files are given, they should correspond to input trees in order."
    )

    args = parser.parse_args()

    prefix_dir   = os.path.dirname(os.path.abspath(args.output_prefix))
    os.makedirs(prefix_dir, exist_ok=True)

    log_file = f"{args.output_prefix}.log"                   
    setup_logging(log_file)

    start_time = datetime.now()
    logging.info("========== Execution started ==========")
    logging.info(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M')}")

    trees_bp = []
    trees_ete = []
    loaded_files = []

    for file_path in args.input:
        try:
            bp_tree = Phylo.read(file_path, 'newick')
            assign_parents(bp_tree)
            trees_bp.append(bp_tree)

            newick_str = bp_tree.format('newick')
            ete_tree = Tree(newick_str, format=1)
            for leaf in ete_tree.iter_leaves():
                leaf.add_feature("orig_name", leaf.name)
            trees_ete.append(ete_tree)

            loaded_files.append(os.path.basename(file_path))
        except Exception as e:
            logging.error(f"Failed to load file: {file_path}. Error: {e}")
            sys.exit(1)

    logging.info("Tree files successfully loaded:")
    for fname in loaded_files:
        logging.info(f"  - {fname}")

    if args.outgroups:
        outgroups_list = []
        for og in args.outgroups:
            if ', ' in og:
                logging.error(
                    f"Whitespace detected after comma in outgroup string: '{og}'\n"
                    "Please ensure that IDs within the same outgroup group are comma-separated with no spaces (e.g., NC_023455,NC_035057)."
                )
                sys.exit(1)
            splitted = og.split(',')
            splitted = [x.strip() for x in splitted if x.strip()]
            outgroups_list.append(splitted)

        if len(outgroups_list) == 1:
            outgroups_list = outgroups_list * len(args.input)
        elif len(outgroups_list) != len(args.input):
            logging.error("Invalid input: the number of outgroup groups (-og) must be either 1 or match the number of input tree files (-i).")
            sys.exit(1)
    else:
        outgroups_list = [None] * len(args.input)

    if args.s:
        if len(args.s) not in {1, len(args.input)}:
            logging.error("Invalid input: the number of '-s' files must either match the number of '-i' tree files or be exactly one.")
            sys.exit(1)

    total_green = 0
    total_blue = 0
    total_red = 0

    for i, (bp_tree, ete_tree) in enumerate(zip(trees_bp, trees_ete), start=1):
        tree_start_time = datetime.now()
        logging.info("========================================")
        logging.info(f"Tree {i} analysis started (file: {loaded_files[i - 1]})")
        logging.info(f"Start time: {tree_start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        logging.info("----------------------------------------")

        outgroup_names = set()
        if outgroups_list[i-1]:
            outgroup_names = set_outgroups(ete_tree, outgroups_list[i-1])
        else:
            logging.info("No outgroup information provided. Skipping outgroup assignment.")

        logging.info("------ Clustering Anomaly Analysis ------")
        anomalies, _ = analyze_clustering(bp_tree, outgroup_names)
        
        anomalies_file = f"{args.output_prefix}_anomaly.txt"      
        save_anomalies(anomalies, anomalies_file)
        logging.info(f"Anomalies written to: {anomalies_file}")

        total_seq = len([leaf for leaf in bp_tree.get_terminals() if leaf.name not in outgroup_names])
        if outgroup_names:
            logging.info(f"Total sequences: {total_seq}; Outgroup sequences: {len(outgroup_names)}")
        else:
            logging.info(f"Total sequences: {total_seq}; Outgroup sequences: 0")

        single_count = len(anomalies.get('single', []))
        two_count = len(anomalies.get('two', []))
        three_count = len(anomalies.get('three_or_more', []))

        logging.info(f"Labeled green: {single_count}; blue: {two_count}; red: {three_count}")

        total_green += single_count
        total_blue += two_count
        total_red += three_count

        if args.s and len(args.s) > 0:
            if len(args.s) == len(trees_bp):
                similar_file = args.s[i-1]
            else:
                similar_file = args.s[0]
            annotate_user_defined_similar(similar_file, ete_tree, anomalies)

        logging.info("------ PDF Export ------")
        pdf_out = f"{args.output_prefix}.pdf"                   
        color_tree(ete_tree, anomalies, pdf_out, outgroup_names)

        tree_end_time = datetime.now()
        td_seconds = (tree_end_time - tree_start_time).total_seconds()
        logging.info(f"Analysis finished at: {tree_end_time.strftime('%Y-%m-%d %H:%M:%S')}, Duration: {int(td_seconds)}s")
        logging.info("----------------------------------------")

    save_multiple_nwk_files(trees_bp, prefix_dir)

    logging.info("========================================")
    logging.info("Summary of anomalous sequences:")
    logging.info(f"  - Green-labeled: {total_green}")
    logging.info(f"  - Blue-labeled: {total_blue}")
    logging.info(f"  - Red-labeled: {total_red}")
    total_anomaly = total_green + total_blue + total_red
    logging.info(f"  - Total anomalies: {total_anomaly}")
    logging.info("========================================")

    end_time = datetime.now()
    runtime = end_time - start_time
    logging.info("========== Execution finished ==========")
    logging.info(f"End time: {end_time.strftime('%Y-%m-%d %H:%M')}")
    total_seconds = int(runtime.total_seconds())
    if total_seconds < 60:
        show_runtime = f"{total_seconds} seconds"
    else:
        show_runtime = f"{total_seconds // 60} minutes"

    logging.info(f"Total runtime: {show_runtime}")
    logging.info(f"All results saved to: {prefix_dir}")

    final_msg = [
        f"Results saved in: {prefix_dir}",
        "- Clustering anomaly files for each tree (clustering_anomalies_*.txt)",
        "- Colored PDF tree visualizations (colored_tree_*.pdf)",
        f"- Newick files saved in: {os.path.join(prefix_dir, 'nwk_files')}",
        f"\nTotal runtime: {runtime}",
    ]
    print("\n".join(final_msg))

if __name__ == "__main__":
    main()