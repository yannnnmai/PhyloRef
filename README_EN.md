[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15489512.svg)](https://doi.org/10.5281/zenodo.15489512)
> 📌 **Archived release available via Zenodo**  
> Cite this DOI when referencing PhyloRef in publications.

# PhyloRef 🧬 — Build & QC Mitochondrial Reference Libraries with Snakemake

*Semi-automated: download ▶︎ filter ▶︎ extract ▶︎ align ▶︎ phylogeny ▶︎ anomaly detection ▶︎ clean.*

> **Why PhyloRef?**  
> Environmental‐DNA (eDNA) studies rely on clean reference databases. PhyloRef takes raw GenBank records and turns them into a quality-controlled, phylogenetically validated FASTA+taxonomy library—reproducibly and at scale.

---

## 🌟 Key Features

| Step | Rule                    | Script                                       | What it does                                                                          |
|------|-------------------------|----------------------------------------------|----------------------------------------------------------------------------------------|
| 1    | `download_mtgenomes`    | `download_accession.py` / `download_name.py` | Batch-download GenBank files by accession **or** species name                         |
| 2    | `filter_sequences`      | `filter_sequences.py`                        | Retain complete genomes *or* target-gene records; optional hybrid / cf. / sp. filters |
| 3    | `extract_genes`         | `extract_genes.py`                           | Extract one or more genes and build *Concat.fa*                                       |
| 4    | `group_by_order`        | `group_by_order.py`                          | Split Concat by taxonomic **Order** (checkpoint)                                      |
| 5    | `build_tree_single`     | (built-in shell)                             | MAFFT alignment → FastTree phylogeny per group                                        |
| 6    | `detect_anomaly_single` | `detect_anomaly.py`                          | Tree-based anomaly detection, outputs PDF for manual review                           |
| 7    | `clean_database`        | `clean_by_anomaly.py`                        | Remove flagged sequences, write final **cleaned DB**                                  |

*Fully checkpoint-safe & parallelisable.*

---

## 🚀 Quick Start

```bash
# clone & enter project
git clone https://github.com/yannnnmai/PhyloRef.git
cd PhyloRef

# --- Micromamba bootstrap (no root needed) -------------------------------
mkdir -p ~/micromamba_bin
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C ~/micromamba_bin
export PATH=~/micromamba_bin/bin:$PATH          # add to ~/.bashrc for permanence
micromamba --version                            # sanity check
# ------------------------------------------------------------------------

# create & activate environment
micromamba create -n phyloref_env -f environment.yaml
eval "$(micromamba shell hook --shell bash)"
micromamba activate phyloref_env

# run demo until anomaly detection (steps 1–6½)
snakemake -j 8 --until anomaly_done
# └─ Uses example input from resources/ (~100 sequences)
# └─ PDF + TXT reports saved in results/6_anomaly/

# after manual review of anomaly reports
snakemake clean_database -j 8
# final cleaned DB saved to: results/7_final_database/outputs/cleaned/

All outputs go to results/, logs to logs/.

⸻

📥 Supported Input

File	Description
resources/accessions.txt	List of GenBank accession numbers (1 per line)
resources/species.txt	List of Latin species names (e.g. Cyprinus carpio)


⸻

🛠 Full Workflow
	1.	Edit config/config.yaml — set download mode, filters, outgroups, etc.
	2.	Place your accession/species list in resources/.
	3.	Run up to anomaly detection (first 6 steps):

snakemake -j <NCPU> --until anomaly_done


	4.	After reviewing trees and editing reports, clean the database:

snakemake clean_database -j <NCPU>



⸻

🔍 Manual Review (Step 6½)

PhyloRef separates anomaly detection from cleaning to enable manual expert review.

💡 Review Workflow

[Run until anomaly_done]
        ↓
[Check groups_*.pdf]
        ↓
[Edit *_anomaly.txt + similar.txt]
        ↓
[Run clean_database]
        ↓
[Get final cleaned DB]

📊 Anomaly Types

Anomalies are auto-detected and color-coded in the trees:
	•	🟢 Type I (green): Singleton outlier — far from expected clade
	•	🔵 Type II (blue): Split pair — two sequences from the same species cluster apart
	•	🔴 Type III (red): Minority outlier — deviating subset in species with ≥3 sequences

📁 Files You Review

File	Description
results/6_anomaly/*.pdf	Color-coded phylogenetic trees
results/6_anomaly/*.txt	Corresponding anomaly reports with flagged sequence IDs
resources/similar.txt	Sequences indistinguishable from others (but not wrong)

✅ After Review
	•	Remove false positives from *_anomaly.txt
	•	Add indistinguishable but correct sequences to resources/similar.txt
	•	Then run:

snakemake clean_database -j 8


⸻

📂 Output Folders

results/
├── 1_download/outputs/...
├── 2_filter/outputs/...
├── 3_extracted/outputs/...
├── 4_groups/...
├── 5_tree/...
├── 6_anomaly/
│   ├── groups_*.pdf
│   └── groups_*_anomaly.txt
└── 7_final_database/outputs/
    ├── cleaned/             # Final QC-ed sequences
    ├── blacklisted/         # Removed problematic sequences
    └── remove_similar.log   # Log of similarity-based removals


⸻

### 🔧 Main Configuration Options (`config/config.yaml`)

| Key                   | Example                  | Description                                             |
|-----------------------|--------------------------|---------------------------------------------------------|
| `download.mode`       | `accession` / `species_name` | Download by accessions or species names                 |
| `filter.mode`         | `gene` / `complete`      | Filter by gene + length or "complete genome" tag       |
| `filter.gene`         | `12S`                    | Target gene name for filtering                         |
| `extract.genes`       | `["12S", "16S"]`         | Genes to extract into `Concat.fa`                      |
| `grouping.max_per_file` | `2000`                 | Max sequences per group before splitting               |
| `tree.threads`        | `30`                     | Number of threads for MAFFT                            |
| `outgroups`           | `NC_023455, NC_035057`   | Optional outgroup IDs for rooting trees                |


⸻

⚙️ Reproducibility
	•	Environment fully defined in environment.yaml
	•	CI tested: GitHub Actions validates demo config on every push
	•	Version control: Archived on Zenodo

⸻

📖 Citing PhyloRef

Mai Y. et al. PhyloRef: A Snakemake-Based Workflow for Semi-Automated Reference Library Curation via Phylogenetic Anomaly Detection (in prep.).

Zenodo DOI: 10.5281/zenodo.15489512

⸻

🤝 Contributing
	1.	Fork → create feature branch (feat/<name>)
	2.	Use Conventional Commits (feat:, fix:, docs:…)
	3.	Submit Pull Request — CI must pass
	4.	Suggestions, issues, and feature requests welcome!

Planned: pre-commit + black + ruff for linting and style

⸻

🛡 License

MIT License — see LICENSE

⸻

👤 Contact

Maintainer: Yan Mai
College of Fisheries and Life Science, Shanghai Ocean University
📧 Email: yann_maii@foxmail.com

⸻

Happy referencing — may your trees be monophyletic!
