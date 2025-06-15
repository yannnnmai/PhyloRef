📌 Zenodo 上提供归档版本
在学术论文中引用 PhyloRef 时请使用此 DOI。

PhyloRef 🧬 — 使用 Snakemake 构建与质控线粒体参考库

半自动化流程：下载 ▶︎ 筛选 ▶︎ 提取 ▶︎ 比对 ▶︎ 构树 ▶︎ 异常检测 ▶︎ 清洗.

为什么选择 PhyloRef？
eDNA（环境 DNA）研究依赖于干净、可靠的参考数据库。PhyloRef 可将原始 GenBank 序列转换为高质量、经过系统发育验证的 FASTA+分类数据库，支持可重复和大规模分析。

⸻

🌟 主要功能

步骤	Snakemake 规则	脚本文件	功能说明
1	download_mtgenomes	download_accession.py / download_name.py	通过 accession 或物种名批量下载 GenBank 文件
2	filter_sequences	filter_sequences.py	保留完整线粒体或目标基因序列，支持 hybrid / cf. / sp. 过滤
3	extract_genes	extract_genes.py	提取一个或多个基因，生成 Concat.fa
4	group_by_order	group_by_order.py	按纲目（Order）分组序列，便于并行分析
5	build_tree_single	（内置 shell）	每组进行 MAFFT 比对 → FastTree 构建系统发育树
6	detect_anomaly_single	detect_anomaly.py	基于系统发育树的异常检测，输出 PDF 供人工审核
7	clean_database	clean_by_anomaly.py	移除标记异常的序列，生成最终 清洗数据库

全流程支持断点续跑与并行执行。

⸻

🚀 快速开始

# 克隆并进入项目目录
git clone https://github.com/yannnnmai/PhyloRef.git
cd PhyloRef

# --- 安装 micromamba（无需管理员权限） -----------------------------
mkdir -p ~/micromamba_bin
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C ~/micromamba_bin
export PATH=~/micromamba_bin/bin:$PATH        # 可加入 ~/.bashrc 永久生效
micromamba --version                          # 检查安装结果
# -------------------------------------------------------------------

# 创建并激活环境
micromamba create -n phyloref_env -f environment.yaml
eval "$(micromamba shell hook --shell bash)"
micromamba activate phyloref_env

# 使用示例数据运行前六步（含异常检测）
snakemake -j 8 --until anomaly_done
# └─ 示例输入位于 resources/（约 100 条序列）
# └─ PDF 和 TXT 异常报告保存在 results/6_anomaly/

# 完成人工审核后执行清洗步骤
snakemake clean_database -j 8
# └─ 最终清洗数据库输出至 results/7_final_database/outputs/cleaned/

所有结果保存在 results/，日志位于 logs/。

⸻

📥 支持的输入文件

文件名	描述
resources/accessions.txt	GenBank accession 编号列表（每行一个）
resources/species.txt	拉丁学名物种列表（如 Cyprinus carpio）

⸻

🛠 完整流程
	1.	编辑 config/config.yaml 配置文件：设置下载模式、筛选条件、外群等。
	2.	将 accession 或物种列表放入 resources/ 文件夹。
	3.	执行前六步：

snakemake -j <NCPU> --until anomaly_done

	4.	审核 PDF 树图与 TXT 报告后执行清洗：

snakemake clean_database -j <NCPU>

⸻

🔍 异常检测与人工审核（第 6½ 步）

PhyloRef 将异常检测与清洗分开，以便手动核查可疑序列。

💡 审核流程

[执行 anomaly_done]
      ↓
[查看 groups_*.pdf 树图]
      ↓
[编辑 *_anomaly.txt 与 similar.txt]
      ↓
[执行 clean_database]
      ↓
[获得清洗后的最终数据库]

📊 异常类型说明

系统自动识别并使用颜色标记异常类型：
	•	🟢 Type I：单一异常序列，远离其应属类群
	•	🔵 Type II：同一物种的两个序列未聚类
	•	🔴 Type III：多个序列中部分偏离（当序列数 ≥3）

📁 需要审核的文件

文件	说明
results/6_anomaly/*.pdf	带颜色标记的系统发育树
results/6_anomaly/*.txt	对应的异常序列列表
resources/similar.txt	与其他序列极度相似但无明显错误的序列列表

✅ 审核完成后：
	•	从 *_anomaly.txt 中移除误判的条目
	•	将正确但不可区分的序列添加至 resources/similar.txt
	•	然后执行：

snakemake clean_database -j 8

⸻

📂 输出目录结构

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
    ├── cleaned/             # 最终保留序列
    ├── blacklisted/         # 被移除的序列
    └── remove_similar.log   # 相似序列处理记录

⸻

🔧 主配置选项（config/config.yaml）

配置键	示例值	描述
download.mode	accession / species_name	下载模式：按 accession 或物种名下载
filter.mode	gene / complete	筛选方式：目标基因 + 最长 / 完整线粒体
filter.gene	12S	目标基因名（与 filter.mode: gene 联用）
extract.genes	["12S", "16S"]	要提取用于建树的目标基因
grouping.max_per_file	2000	每组最大序列数，超出自动分组
tree.threads	30	MAFFT 多线程数量
outgroups	NC_023455, NC_035057	可选的外群序列 ID

⸻

⚙️ 可重复性保障
	•	所有依赖已封装于 environment.yaml
	•	GitHub Actions 自动验证示例配置
	•	Zenodo 长期归档，支持版本控制与引用

⸻

📖 引用 PhyloRef

Mai Y. et al. PhyloRef: A Snakemake-Based Workflow for Semi-Automated Reference Library Curation via Phylogenetic Anomaly Detection (in prep.)

Zenodo DOI: 10.5281/zenodo.15489512

⸻

🤝 贡献方式
	1.	Fork 本仓库 → 创建功能分支（如 feat/<功能名>）
	2.	使用 Conventional Commits 格式（如 feat: / fix: / docs:）
	3.	提交 Pull Request，CI 验证需通过
	4.	欢迎提出建议、报告 bug 或新功能请求

未来将集成 pre-commit、black、ruff 进行格式化与风格校验。

⸻

🛡 许可协议

本项目采用 MIT License（详见 LICENSE 文件）

⸻

👤 联系方式

维护者：颜迈（Yan Mai）
单位：上海海洋大学 水产与生命学院
📧 邮箱：yann_maii@foxmail.com

⸻

祝你的系统发育树永远是单系的 🧬🌿
