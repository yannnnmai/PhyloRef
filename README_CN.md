
📌 Zenodo 上提供归档版本
在论文中引用 PhyloRef 时，请使用此 DOI。

PhyloRef 🧬 — 使用 Snakemake 构建与质控线粒体参考库

半自动化流程：下载 ▶︎ 筛选 ▶︎ 基因提取 ▶︎ 比对 ▶︎ 构建系统发育树 ▶︎ 异常检测 ▶︎ 序列清洗

为什么使用 PhyloRef？
环境 DNA（eDNA）研究依赖于干净、可靠的参考数据库。PhyloRef 可将原始 GenBank 序列转化为经过质量控制、系统发育验证的 FASTA+分类标签库，具有高度可重复性和良好的扩展性。

⸻

🌟 核心功能概览

步骤	Snakemake 规则	Python 脚本	功能说明
1	download_mtgenomes	download_accession.py / download_name.py	批量下载 GenBank 文件（可通过 accession 或物种名）
2	filter_sequences	filter_sequences.py	保留完整线粒体或目标基因序列；支持 hybrid / cf. / sp. 筛选
3	extract_genes	extract_genes.py	提取一个或多个基因，构建 Concat.fa
4	group_by_order	group_by_order.py	按纲目（Order）对序列分组（可并行处理）
5	build_tree_single	内置 shell 脚本	使用 MAFFT + FastTree 构建系统发育树
6	detect_anomaly_single	detect_anomaly.py	检测系统发育异常，输出 PDF 报告供人工审核
7	clean_database	clean_by_anomaly.py	根据人工审核结果清除异常序列，输出最终参考库


⸻

🚀 快速上手

# 克隆项目并进入目录
git clone https://github.com/yannnnmai/PhyloRef.git
cd PhyloRef

# --- 安装 micromamba（不需要管理员权限） -----------------------------
mkdir -p ~/micromamba_bin
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C ~/micromamba_bin
export PATH=~/micromamba_bin/bin:$PATH      # 可写入 ~/.bashrc 永久生效
micromamba --version                        # 检查是否安装成功
# ----------------------------------------------------------------------

# 创建并激活环境
micromamba create -n phyloref_env -f environment.yaml
eval "$(micromamba shell hook --shell bash)"
micromamba activate phyloref_env

# 使用测试数据运行至第 6 步（异常检测）
snakemake -j 8 --until anomaly_done
# └─ 测试数据位于 resources/，约 100 条序列
# └─ PDF 和 TXT 异常报告保存在 results/6_anomaly/

# 人工审核并修改后，运行清洗步骤生成最终数据库
snakemake clean_database -j 8
# └─ 清洗后数据库输出到 results/7_final_database/outputs/cleaned/

# 所有中间结果保存在 results/，日志输出至 logs/


⸻

📥 支持的输入文件

文件名	说明
resources/accessions.txt	accession 编号列表（每行一个）
resources/species.txt	拉丁学名物种列表（如 Cyprinus carpio）


⸻

🛠 完整流程（7 步）
	1.	编辑配置文件：config/config.yaml（下载模式、筛选条件、外群设置等）
	2.	放入你的物种列表或 accession 编号文件到 resources/ 文件夹
	3.	运行至第六步（系统发育树和异常检测）：

snakemake -j <NCPU> --until anomaly_done

	4.	人工审核 PDF 树图和 TXT 报告后，执行清洗：

snakemake clean_database -j <NCPU>


⸻

🔍 异常检测与人工审核（第 6½ 步）

PhyloRef 专门将异常检测和最终清洗分离，便于人工干预和判断。

✅ 审核流程：

运行 anomaly_done
      ↓
查看 groups_*.pdf 树图
      ↓
编辑 anomaly.txt + similar.txt
      ↓
执行 clean_database
      ↓
获得清洗后的最终数据库

📊 异常类型说明

类型	颜色	含义
I	🟢	单一序列偏离其应归属的类群（属/科/目）
II	🔵	同种的两个序列在树上不聚为一类
III	🔴	多序列种中有少数偏离主群（≥3 条序列的情况）

📁 你需要审核的文件

文件	说明
results/6_anomaly/*.pdf	树图，带颜色标记异常序列
results/6_anomaly/*.txt	对应的异常序列 ID 列表
resources/similar.txt	无明显错误但与其他序列高度相似的条目


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
    ├── blacklisted/         # 被移除的异常序列
    └── remove_similar.log   # 相似序列处理日志


⸻

🔧 主配置项（config/config.yaml）

配置键	示例值	含义
download.mode	accession / species_name	下载模式：通过 accession 或物种名
filter.mode	gene / complete	按基因+长度或完整基因组标签筛选
filter.gene	12S	筛选基因名
extract.genes	["12S", "16S"]	提取的目标基因
grouping.max_per_file	2000	每组最大序列数，超出将自动分组
tree.threads	30	MAFFT 使用的线程数
outgroups	NC_023455, NC_035057	可选外群序列 ID，用于树的根定向


⸻

⚙️ 可重复性保证
	•	所有依赖环境在 environment.yaml 中完整定义
	•	GitHub Actions 自动测试 demo 流程
	•	版本已在 Zenodo 存档，确保长期可引用性

⸻

📖 引用方法

Mai Y. et al. PhyloRef: A Snakemake-Based Workflow for Semi-Automated Reference Library Curation via Phylogenetic Anomaly Detection (in prep.).
DOI: 10.5281/zenodo.15489512


⸻

🤝 欢迎贡献
	1.	Fork 仓库 → 创建新分支（如 feat/xxx）
	2.	提交前请遵循 Conventional Commits 格式（如：feat:、fix:、docs:）
	3.	发起 Pull Request，自动测试需通过
	4.	欢迎提出建议、报告 bug 或提交新功能！

后续计划集成：pre-commit、black、ruff 等自动格式化工具。

⸻

🛡 项目许可

本项目使用 MIT 协议。

⸻

👤 联系方式

维护者：颜迈（Yan Mai）
单位：上海海洋大学 水产与生命学院
邮箱：📧 yann_maii@foxmail.com

⸻

🌿 祝你的参考库构建顺利，愿你的系统发育树单系一致！
