[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15489512.svg)](https://doi.org/10.5281/zenodo.15489512)  
> 📌 **已归档版本可通过 Zenodo 获取**  
> 在引用 PhyloRef 时请使用该 DOI。

# PhyloRef 🧬 — 基于 Snakemake 构建和质控线粒体参考数据库

*半自动流程：下载 ▶︎ 筛选 ▶︎ 提取 ▶︎ 比对 ▶︎ 构树 ▶︎ 异常检测 ▶︎ 清洗。*

> **为什么选择 PhyloRef？**  
> 环境 DNA（eDNA）研究高度依赖干净的参考数据库。PhyloRef 可将原始 GenBank 序列转化为质量受控、系统发育验证的 FASTA+分类库，具备良好的可重复性和扩展性。

---

## 🌟 核心功能

| 步骤 | 规则名                  | 脚本名                                        | 功能说明                                                                                   |
|------|--------------------------|-----------------------------------------------|--------------------------------------------------------------------------------------------|
| 1    | `download_mtgenomes`     | `download_accession.py` / `download_name.py`  | 批量下载 GenBank 文件（按 accession 或物种名）                                             |
| 2    | `filter_sequences`       | `filter_sequences.py`                         | 保留完整线粒体或目标基因记录，支持去除杂交 / cf. / sp. 等标签                              |
| 3    | `extract_genes`          | `extract_genes.py`                            | 提取指定基因并生成 *Concat.fa*                                                              |
| 4    | `group_by_order`         | `group_by_order.py`                           | 按分类 Order 拆分 Concat（为 checkpoint 设置）                                              |
| 5    | `build_tree_single`      | （内置 shell 命令）                           | 每个分组进行 MAFFT 比对 → FastTree 构树                                                    |
| 6    | `detect_anomaly_single`  | `detect_anomaly.py`                           | 基于系统发育树的异常检测，输出 PDF 供人工审查                                              |
| 7    | `clean_database`         | `clean_by_anomaly.py`                         | 移除异常序列，输出最终 **清洗后的数据库**                                                   |

---

## 🚀 快速开始

```bash
# 克隆并进入项目目录
git clone https://github.com/yannnnmai/PhyloRef.git
cd PhyloRef

# --- Micromamba 安装（无需 root 权限） -------------------------------
mkdir -p ~/micromamba_bin
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C ~/micromamba_bin
export PATH=~/micromamba_bin/bin:$PATH          # 可添加至 ~/.bashrc 以永久生效
micromamba --version                            # 检查是否安装成功
# ------------------------------------------------------------------------

# 创建并激活环境
micromamba create -n phyloref_env -f environment.yaml
eval "$(micromamba shell hook --shell bash)"
micromamba activate phyloref_env

# 运行示例数据至异常检测（步骤 1–6½）
snakemake -j 8 --until anomaly_done
# └─ 使用 resources/ 目录下提供的示例输入（约 100 条序列）
# └─ PDF 和 TXT 报告保存在 results/6_anomaly/ 中

# 人工审查异常检测报告后，运行清洗流程
snakemake clean_database -j 8
# 最终清洗后的数据库保存在：results/7_final_database/outputs/cleaned/

所有输出结果保存在 results/，日志信息位于 logs/。

⸻

---

📥 支持的输入格式

| 文件路径                   | 描述                                                 |
|----------------------------|------------------------------------------------------|
| resources/accessions.txt   | GenBank 序列号列表，每行一个                         |
| resources/species.txt      | 拉丁学名的物种列表（例如 Cyprinus carpio）           |

⸻

🛠 完整工作流程

1. 编辑 `config/config.yaml` —— 设置下载方式、筛选模式、目标基因、系统发育树参数等；
2. 将 accession 列表或物种列表放入 `resources/` 文件夹；
3. 运行前六步流程（至异常检测阶段）：snakemake -j <NCPU> --until anomaly_done
4. 审查系统发育树（PDF）和异常报告（TXT），完成人工标注后，运行最终清洗：snakemake clean_database -j <NCPU>

🔍 人工审查（第 6½ 步）

PhyloRef 将“异常检测”与“数据库清洗”分离，以支持手动专家审查。

💡 审查工作流程

[执行 anomaly_done 步骤]
        ↓
[查看 groups_*.pdf]
        ↓
[编辑 *_anomaly.txt 与 similar.txt]
        ↓
[执行 clean_database]
        ↓
[获得最终清洗后的数据库]

📊 异常类型

PhyloRef 自动检测并在系统发育树中以颜色标记三类常见异常：

- 🟢 类型 I（绿色）：单个序列远离其所属的纲、科、属的其他序列 → 孤立异常；
- 🔵 类型 II（蓝色）：同一物种的两个序列分布在树的不同位置 → 分裂对；
- 🔴 类型 III（红色）：某物种中 ≥3 条序列有一部分偏离主簇 → 少数偏离群体。

📁 审查的文件

| 文件路径                     | 描述                                               |
|------------------------------|----------------------------------------------------|
| results/6_anomaly/*.pdf      | 彩色系统发育树，标记潜在异常                      |
| results/6_anomaly/*.txt      | 与 PDF 匹配的异常序列报告，含待清洗的序列 ID      |
| resources/similar.txt        | 可保留的相似序列（与其他序列几乎一致，但并非错误）|

✅ 审查后操作

- 若判断某些异常是误报（如系统发育树构建误差），请将对应序列 ID 从 `*_anomaly.txt` 中删除；
- 若发现一些与其他序列非常相似但无系统发育偏离，请将其添加至 `resources/similar.txt`；
- 然后运行：

```bash
snakemake clean_database -j 8


📂 输出文件结构

results/
├── 1_download/outputs/...         # 下载的 GenBank 文件
├── 2_filter/outputs/...           # 筛选后的序列文件
├── 3_extracted/outputs/...        # 提取的基因序列（Concat.fa）
├── 4_groups/...                   # 按纲（Order）分组的序列
├── 5_tree/...                     # 系统发育树构建结果
├── 6_anomaly/                     # 异常检测输出
│   ├── groups_*.pdf               # 异常标记的系统发育树
│   └── groups_*_anomaly.txt       # 对应的异常序列报告
└── 7_final_database/outputs/      # 最终数据库输出
    ├── cleaned/                   # 清洗后保留的高质量序列
    ├── blacklisted/               # 被删除的异常序列
    └── remove_similar.log         # 标记相似序列的日志记录

⸻

### 🔧 主配置参数说明（`config/config.yaml`）

| 参数键名               | 示例值                        | 描述                                                  |
|------------------------|-------------------------------|-------------------------------------------------------|
| `download.mode`        | `accession` / `species_name`  | 选择使用 Accession 编号或物种学名进行下载            |
| `filter.mode`          | `gene` / `complete`           | 选择按基因 + 序列长度筛选，或仅保留完整基因组        |
| `filter.gene`          | `12S`                         | 筛选时使用的目标基因名称（如 12S、COI）              |
| `extract.genes`        | `["12S", "16S"]`              | 从序列中提取的基因列表，并用于拼接成 Concat.fa       |
| `grouping.max_per_file`| `2000`                        | 每个分组最大序列数，超出将拆分为多个子组             |
| `tree.threads`         | `30`                          | 系统发育树构建中 MAFFT 的线程数                      |
| `outgroups`            | `NC_023455, NC_035057`        | 系统发育树构建时用于定根的外类群序列 ID（可选）      |

⸻

⚙️ 可重复性保障

- 完整环境配置记录于 `environment.yaml`，可通过 micromamba 创建环境
- 持续集成：每次提交自动使用 GitHub Actions 运行测试（基于 demo 配置）
- 数据库版本控制：定期归档至 Zenodo，确保可引用性与长期可用性

⸻

📖 引用方式

Mai Y. 等人，PhyloRef: 基于 Snakemake 的半自动参考数据库系统发育异常检测与清洗流程（论文准备中）

Zenodo DOI: 10.5281/zenodo.15489512

> 请在学术引用中使用上述 DOI。

⸻

🤝 贡献指南

1. Fork 本项目，并创建新的功能分支（例如：`feat/<功能名>`）
2. 遵循 Conventional Commits 风格提交代码（例如：`feat:`、`fix:`、`docs:` 等）
3. 提交 Pull Request，必须通过 CI 测试
4. 欢迎任何问题反馈、建议和功能请求！

🚧 规划中：将集成 pre-commit、black、ruff 以实现统一的代码风格检查和格式化

⸻

🛡 授权许可

MIT License — 详见 LICENSE 文件

⸻

👤 联系方式

维护者：Yan Mai  
单位：上海海洋大学 水产与生命学院  
📧 邮箱：yann_maii@foxmail.com

⸻

🌿 祝你构建的参考数据库系统发育关系清晰、单系演化！
