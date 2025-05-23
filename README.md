[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15489512.svg)](https://doi.org/10.5281/zenodo.15489512)
## PhyloRef 🧬 — Build & QC Mitochondrial Reference Libraries with Snakemake

*Semi-automated download ▶︎ filter ▶︎ extract ▶︎ align ▶︎ phylogeny ▶︎ anomaly detection ▶︎ clean.*

> **Why PhyloRef?**
> Environmental‐DNA (eDNA) studies rely on clean reference databases. PhyloRef takes raw GenBank records and turns them into a quality-controlled, phylogenetically validated FASTA+taxonomy library—reproducibly and at scale.

---

### 🌟 Key Features

| Step | Rule                    | Script                                       | What it does                                                                          |
| ---- | ----------------------- | -------------------------------------------- | ------------------------------------------------------------------------------------- |
| 1    | `download_mtgenomes`    | `download_accession.py` / `download_name.py` | Batch-download GenBank files by accession **or** species name                         |
| 2    | `filter_sequences`      | `filter_sequences.py`                        | Retain complete genomes *or* target-gene records; optional hybrid / cf. / sp. filters |
| 3    | `extract_genes`         | `extract_genes.py`                           | Extract one or more genes and build *Concat.fa*                                       |
| 4    | `group_by_order`        | `group_by_order.py`                          | Split Concat by taxonomic **Order** (checkpoint)                                      |
| 5    | `build_tree_single`     | (built-in shell)                             | MAFFT alignment → FastTree phylogeny per group                                        |
| 6    | `detect_anomaly_single` | `detect_anomaly.py`                          | Tree-based anomaly detection, outputs PDF for manual review                           |
| 7    | `clean_database`        | `clean_by_anomaly.py`                        | Remove flagged sequences, write final **cleaned DB**                                  |

*Fully checkpoint-safe & parallelisable.*

---

### 🚀 Quick Start (Demo)

```bash
# clone & enter project
git clone https://github.com/<YOUR_GITHUB>/PhyloRef.git
cd PhyloRef

# create environment (≈1–2 min with micromamba)
micromamba env create -f environment.yaml
micromamba activate phyloref_env

# run minimal test dataset (resources/demo + config/demo_config.yaml)
snakemake \
    --configfile config/demo_config.yaml \
    -j 4 --use-conda --conda-frontend mamba
```

Outputs land in `results/`, logs in `logs/`.
Check `results/anomaly/*.pdf` and mark sequences to blacklist if necessary.

---

### 🛠 Full Workflow

1. Edit `config/config.yaml` — switch *download.mode* to `accession` **or** `species_name`, adjust filter/extract options.
2. Place your accession/species list in `resources/`.
3. Run:

```bash
snakemake -j <CPU> --use-conda --conda-frontend mamba \
          --rerun-incomplete --printshellcmds
```

Large datasets can be resumed any time; Snakemake tracks completed jobs.

---

### 🔍 Manual Review Step (Tree-based anomaly detection)

PhyloRef intentionally separates anomaly detection from final filtering to allow **manual review and expert judgment**.

#### 💡 Workflow Summary:

1. **Run first 6 steps automatically:**

   ```bash
   snakemake -j 8 --forcerun all --use-conda --conda-frontend mamba
   ```

   This will generate anomaly reports at:

   ```
   results/anomaly/groups_*.pdf
   results/anomaly/groups_*_anomaly.txt
   ```

2. **Open each `*.pdf` and visually inspect the tree**:
   Anomalies are categorized into three types:

   * 🟢 **Type I (green)**: Single sequences failing to cluster with their taxonomic group.
   * 🔵 **Type II (blue)**: Two sequences from the same species clustering separately.
   * 🔴 **Type III (red)**: Minority sequences deviating from the main cluster in species with three or more records.

   These are flagged both in `*_anomaly.txt` and visualized via color-coded phylogenetic trees (`*.pdf`) to assist manual review, significantly reducing inspection workload.

3. **Minority deviations within multi-sequence clusters**

   All flagged cases are visually indicated in color-coded phylogenetic trees for manual review.

4. **Edit the corresponding `*_anomaly.txt` file**:

   * Remove any **sequences you judge to be correct** (e.g. type specimens, museum vouchers).
   * Also remove sequences that are **not wrong but indistinguishable** from others.

5. **Add all indistinguishable but non-error sequences to:**

   ```
   resources/similar.txt
   ```

   Format: one sequence ID per line (e.g., `NC_012345`)

6. **Then run the final cleaning step:**

   ```bash
   snakemake clean_database -j 8 --forcerun --use-conda --conda-frontend mamba
   ```

7. Final results are saved to:

   ```
   results/final_database/outputs/cleaned/
   ```

---

### 📌 File Definitions

| File                      | Description                                                         |
| ------------------------- | ------------------------------------------------------------------- |
| `groups_1_63_anomaly.txt` | List of automatically flagged outliers                              |
| `groups_1_63.pdf`         | Tree diagram with color-coded anomaly branches                      |
| `similar.txt`             | Sequences that are indistinguishable from others, but not erroneous |

---

### 📂 Output Folders (default: `results/final_database/outputs/`)

```
final_database/
└── outputs/
    ├── cleaned/                 
    ├── blacklisted/             
    └── remove_similar.log       
```

---

### 🔧 Main Configuration Options

| Key                     | Example                      | Meaning                                         |
| ----------------------- | ---------------------------- | ----------------------------------------------- |
| `download.mode`         | `accession` / `species_name` | Select download strategy                        |
| `filter.mode`           | `gene` / `complete`          | Filter by gene+length or “complete genome” flag |
| `filter.gene`           | `12S`                        | Target gene when `filter.mode: gene`            |
| `extract.genes`         | `["12S", "16S"]`             | Genes to extract into Concat.fa                 |
| `grouping.max_per_file` | `2000`                       | Split point for large orders                    |
| `tree.threads`          | `30`                         | MAFFT CPU threads                               |
| `outgroups`             | `- NC_023455, NC_035057`     | IDs passed to anomaly detector                  |

All parameters are documented inline in `config/config.yaml`.

---

### ⚙️ Reproducibility

* **Environment** fully specified in [`environment.yaml`](environment.yaml).
* **Continuous Integration** via GitHub Actions runs the demo workflow on every push.
* First public release archived on Zenodo → citeable DOI.

---

### 📖 Citing PhyloRef Tool

```text
Mai Y. et al. PhyloRef: A Snakemake-Based Workflow for Semi-Automated Reference Library Curation via Phylogenetic Anomaly Detection (in prep.).
```

DOI will be available after the first GitHub release is archived on Zenodo.

---

### 🤝 Contributing

1. Fork → create feature branch (`feat/<name>`).
2. Commit following *Conventional Commits* (`feat:`, `fix:`, `docs:`…).
3. Open Pull Request; CI must pass.
4. Feel free to open issues for bugs or feature requests.

Code style: *black ∘ ruff ∘ pre-commit* (coming soon).

---

### 🛡 License

This project is licensed under the [MIT License](LICENSE).

---

### 👤 Contact

*Maintainer*: **Yan Mai** · `yann_maii@foxmail.com`
College of Fisheries and Life Science, Shanghai Ocean University, Shanghai, China
---

> *Happy referencing & may your trees be monophyletic!*
