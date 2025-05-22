import os
import argparse
from collections import defaultdict
from datetime import datetime

def parse_args():
    parser = argparse.ArgumentParser(
        description="Dynamically group sequences by order based on a maximum per file threshold."
    )
    parser.add_argument("-i", "--input_file", type=str, required=True,
                        help="Path to the input FASTA file")
    parser.add_argument("-o", "--output_dir", type=str, required=True,
                        help="Directory where the groups folder will be created to save the grouped sequences")
    parser.add_argument("-og", "--outgroup_ids", type=str, nargs='+', default=[],
                        help="List of outgroup sequence IDs (optional). Example: -og NC_023455 NC_035057")
    parser.add_argument("-m", "--max_per_file", type=int, default=2000,
                        help="Maximum total number of sequences per file (default: 2000)")
    return parser.parse_args()

def get_order(header):
    if "|o__" in header:
        return header.split("|o__")[1].split("|")[0]
    return "No_order"

def is_outgroup(header, outgroup_ids):
    for out_id in outgroup_ids:
        if out_id in header:
            return True
    return False

def read_fasta(file_path, outgroup_ids, outgroup_sequences):
    total = 0
    grouped_sequences = defaultdict(list)
    with open(file_path, "r") as infile:
        current_header = ""
        current_sequence = []
        for line in infile:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    total += 1
                    if is_outgroup(current_header, outgroup_ids):
                        outgroup_sequences.append((current_header, "".join(current_sequence)))
                    else:
                        order = get_order(current_header)
                        grouped_sequences[order].append((current_header, "".join(current_sequence)))
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)
        if current_header:
            total += 1
            if is_outgroup(current_header, outgroup_ids):
                outgroup_sequences.append((current_header, "".join(current_sequence)))
            else:
                order = get_order(current_header)
                grouped_sequences[order].append((current_header, "".join(current_sequence)))
    return total, grouped_sequences

def write_category(file_name, category, outgroup_sequences, output_dir):
    file_path = os.path.join(output_dir, file_name)
    with open(file_path, "w") as outfile:
        for out_header, out_seq in outgroup_sequences:
            outfile.write(f"{out_header}\n{out_seq}\n")
        for order in sorted(category.keys()):
            for header, sequence in category[order]:
                outfile.write(f"{header}\n{sequence}\n")

def classify_orders(grouped_sequences, max_per_file):
    sorted_orders = sorted(grouped_sequences.items(), key=lambda x: len(x[1]), reverse=True)
    categories = []
    current_category = {}
    current_count = 0
    for order, seq_list in sorted_orders:
        order_count = len(seq_list)
        if order_count > max_per_file:
            if current_category:
                categories.append(current_category)
                current_category = {}
                current_count = 0
            categories.append({order: seq_list})
        else:
            if current_count + order_count > max_per_file:
                categories.append(current_category)
                current_category = {order: seq_list}
                current_count = order_count
            else:
                current_category[order] = seq_list
                current_count += order_count
    if current_category:
        categories.append(current_category)
    return categories

def main():
    args = parse_args()
    max_per_file = args.max_per_file

    groups_output_dir = args.output_dir
    os.makedirs(args.output_dir, exist_ok=True)
    log_file = os.path.join(groups_output_dir, "grouping.log")

    outgroup_sequences = []
    total_input, grouped_sequences = read_fasta(args.input_file, args.outgroup_ids, outgroup_sequences)

    outgroup_species = []
    for i, (header, seq) in enumerate(outgroup_sequences):
        original_species = header[1:].split("|")[0]
        sequence_id = header.split("|")[1] if "|" in header else "Unknown"
        if original_species.startswith("Outgroup"):
            new_header = header
            species_name = original_species.split("_", 1)[-1]
        else:
            new_header = f">Outgroup{i+1}_{original_species}" + header[len(">" + original_species):]
            species_name = original_species
        outgroup_sequences[i] = (new_header, seq)
        outgroup_species.append((species_name, sequence_id))
    
    non_outgroup_count = sum(len(seq_list) for seq_list in grouped_sequences.values())
    
    categories = classify_orders(grouped_sequences, max_per_file)

    for i, category in enumerate(categories):
        total_in_cat = sum(len(seq_list) for seq_list in category.values())
        file_name = f"groups_{i+1}_{total_in_cat}.fa"
        write_category(file_name, category, outgroup_sequences, groups_output_dir)

    try:
        with open(log_file, "w") as log:
            log.write("FASTA Grouping Log\n")
            log.write("="*50 + "\n")
            log.write(f"Date and Time: {datetime.now()}\n\n")
            log.write(f"Total sequences in input file: {total_input}\n")
            log.write(f"Total outgroup sequences (matching IDs {', '.join(args.outgroup_ids) if args.outgroup_ids else 'None'}): {len(outgroup_sequences)}\n")
            log.write(f"Total non-outgroup sequences processed into groups: {non_outgroup_count}\n\n")
            
            num_groups = len(categories)
            log.write("Classification Summary:\n")
            log.write("="*50 + "\n")
            log.write(f"We have classified the sequences into {num_groups} major groups based on the number of sequences:\n\n")
            for i, group in enumerate(categories):
                total_in_group = sum(len(seq_list) for seq_list in group.values())
                order_count = len(group)
                log.write(f"Group {i+1}:\n")
                log.write(f"   - This group contains {order_count} order{'s' if order_count != 1 else ''} with a total of {total_in_group} sequences.\n")
                log.write("   Orders in this group:\n")
                for order in sorted(group.keys()):
                    log.write(f"   - {order}: {len(group[order])} sequences\n")
                log.write("\n")
            
            if args.outgroup_ids:
                log.write("Outgroup species:\n")
                for species_name, sequence_id in outgroup_species:
                    log.write(f"  - {species_name} (Sequence ID: {sequence_id})\n")
            else:
                log.write("Outgroup species: None provided.\n")
            
            log.write("\nAdditional Information:\n")
            if args.outgroup_ids:
                log.write("  - Each output file includes all outgroup sequences at the beginning.\n")
            else:
                log.write("  - No outgroup sequences were specified; therefore, no outgroup sequences are included in the output files.\n")
            log.write("="*50 + "\n")
            log.write("Grouping completed successfully.\n")
        print(f"[INFO] Processing complete.")
        print(f"[INFO] Log saved to: {log_file}")

    except Exception as e:
        print(f"Error occurred while generating the log: {e}")
        with open(log_file, "w") as log:
            log.write("FASTA Grouping Log\n")
            log.write("="*50 + "\n")
            log.write(f"Date and Time: {datetime.now()}\n")
            log.write("An error occurred while processing the sequences.\n")
            log.write(f"Error details: {str(e)}\n")
            log.write("="*50 + "\n")

if __name__ == "__main__":
    main()