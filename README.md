# 🧬 AKT1 Inhibitor RNA-seq Project
## STAR-based analysis of immune checkpoint regulation in myeloid cell lines

This repository contains a **classical bulk RNA-seq workflow** for studying the transcriptional response to **PU-001**, a putative **AKT1 inhibitor**, in three MDS-relevant myeloid cell lines:

- **KG-1**
- **Mono-Mac-1**
- **THP-1**

The main biological goal is to test whether **AKT1 inhibition is associated with reduced expression of immune checkpoint genes**, including **CD274 (PD-L1)**, **HAVCR2 (TIM-3)**, **PDCD1**, **PDCD1LG2**, **CTLA4**, **LAG3**, **TIGIT**, **BTLA**, **VSIR**, **CD80**, **CD86**, **CD276**, **VTCN1**, and **IDO1**.

---

## ✨ Project summary

### Biological question
Does treatment with **PU-001** reduce immune checkpoint-associated transcriptional programs in MDS-related myeloid cell lines?

### Experimental design
Each cell line contains **3 biological replicates**, each with a matched **control** and **treated** sample:

- **KG-1**: K1–K6
- **Mono-Mac-1**: M7–M12
- **THP-1**: T13–T18

Treatment conditions:

- **CTRL** — untreated cells
- **PU001** — cells treated with **200 µM PU-001 for 24 h**

### Main hypothesis
**AKT1 inhibition leads to downregulation of checkpoint genes at the transcriptomic level.**

---

## 🧪 Pipeline overview

This repository uses a **STAR-based RNA-seq pipeline**:

1. **FastQC** — raw read quality control
2. **fastp** — adapter and quality trimming
3. **STAR genome index generation**
4. **STAR alignment**
5. **samtools sort + index**
6. **featureCounts** — gene-level counting
7. **DESeq2** — differential expression analysis
8. **Checkpoint-focused summary and heatmap**
9. **MultiQC** — consolidated QC report

---

## 📌 What each file does

### `Snakefile`
Main Snakemake workflow definition. Controls the full execution order from raw FASTQ files to final DESeq2 outputs.

### `config.yaml`
Stores project-wide configuration:
- sample sheet path
- genome FASTA path
- GTF annotation path
- STAR index path
- thread settings
- featureCounts parameters

### `metadata/samples.tsv`
Sample metadata table containing:
- `sample_id`
- `fq1`
- `fq2`
- `cell_line`
- `treatment`
- `replicate`

### `scripts/get_ref.sh`
Downloads and prepares:
- **GRCh38 primary assembly genome FASTA**
- matching **GENCODE annotation GTF**

### `scripts/featurecounts_to_tsv.py`
Converts raw `featureCounts` output into a tidy count matrix with columns matching `sample_id`.

### `scripts/deseq2_pipeline.R`
Runs:
- global DESeq2 model
- per-cell-line DESeq2 models
- PCA
- checkpoint summary table
- checkpoint heatmap

---

## 🚀 How to run the project

## 1️⃣ Create the software environment

Recommended with **mamba**:

```bash
mamba create -n rnaseq -c conda-forge -c bioconda -y \
  python=3.11 \
  snakemake \
  fastqc multiqc fastp \
  star samtools subread \
  r-base=4.3 \
  bioconductor-deseq2 \
  r-data.table r-ggplot2 r-pheatmap \
  bioconductor-rtracklayer
conda activate rnaseq
```

If `mamba` is unavailable, replace it with `conda`.

---

## 2️⃣ Place raw FASTQ files

Put all paired-end FASTQ files into:

```text
fastq/
```

Example filenames:

```text
Unknown_CQ888-001U0001_1.fq.gz
Unknown_CQ888-001U0001_2.fq.gz
...
Unknown_CQ888-001U0018_1.fq.gz
Unknown_CQ888-001U0018_2.fq.gz
```

---

## 3️⃣ Prepare `metadata/samples.tsv`

The file should look like this:

```tsv
sample_id	fq1	fq2	cell_line	treatment	replicate
K1	fastq/Unknown_CQ888-001U0001_1.fq.gz	fastq/Unknown_CQ888-001U0001_2.fq.gz	KG-1	CTRL	1
K2	fastq/Unknown_CQ888-001U0002_1.fq.gz	fastq/Unknown_CQ888-001U0002_2.fq.gz	KG-1	PU001	1
K3	fastq/Unknown_CQ888-001U0003_1.fq.gz	fastq/Unknown_CQ888-001U0003_2.fq.gz	KG-1	CTRL	2
K4	fastq/Unknown_CQ888-001U0004_1.fq.gz	fastq/Unknown_CQ888-001U0004_2.fq.gz	KG-1	PU001	2
...
```

---

## 4️⃣ Download the reference genome and annotation

Run:

```bash
bash scripts/get_ref.sh 49 comprehensive ref
```

This prepares:

- `ref/GRCh38.primary_assembly.genome.fa`
- `ref/gencode.annotation.gtf`

> ⚠️ The FASTA and GTF must match in **assembly** and **release**.

---

## 5️⃣ Review `config.yaml`

Example:

```yaml
samples_tsv: "metadata/samples.tsv"

ref:
  fasta: "ref/GRCh38.primary_assembly.genome.fa"
  gtf:   "ref/gencode.annotation.gtf"
  star_index: "ref/STAR_GRCh38"

threads:
  star: 12
  fastqc: 4
  featureCounts: 8

featurecounts:
  extra: "-p -B -C -t exon -g gene_id"
```

---

## 6️⃣ Run a dry run first

```bash
snakemake -n -p
```

This checks planned jobs and dependencies without executing anything.

---

## 7️⃣ Run the full workflow

```bash
snakemake --cores 12 -p --rerun-incomplete --latency-wait 60 --keep-going
```

---

## 🔬 Workflow details

### 🟢 Rule: `fastqc_raw`
Runs **FastQC** on raw paired-end reads.

Outputs:
- `results/fastqc/raw/<sample>_1_fastqc.html`
- `results/fastqc/raw/<sample>_2_fastqc.html`

### ✂️ Rule: `fastp_trim`
Runs **fastp** trimming.

Outputs:
- `results/trimmed/<sample>_1.fq.gz`
- `results/trimmed/<sample>_2.fq.gz`
- `results/fastp/<sample>.html`
- `results/fastp/<sample>.json`

### 🧱 Rule: `star_index`
Builds the STAR reference index once.

Output:
- `ref/STAR_GRCh38/`

### 🎯 Rule: `star_align`
Runs STAR alignment, then sorts and indexes BAM files with samtools.

Outputs:
- `results/star/<sample>.Aligned.sortedByCoord.out.bam`
- `results/star/<sample>.Aligned.sortedByCoord.out.bam.bai`
- `results/star/<sample>.Log.final.out`

### 🔢 Rule: `featurecounts`
Counts reads at the gene level using **featureCounts**.

Output:
- `results/counts/gene_counts.tsv`

### 📊 Rule: `deseq2`
Runs DESeq2 globally and per cell line.

Outputs:
- `results/deseq2/DE_allCellLines_treatment.tsv`
- `results/deseq2/DE_KG-1_PU001_vs_CTRL.tsv`
- `results/deseq2/DE_Mono-Mac-1_PU001_vs_CTRL.tsv`
- `results/deseq2/DE_THP-1_PU001_vs_CTRL.tsv`
- `results/deseq2/checkpoints_summary.tsv`
- `results/deseq2/figures/PCA_vst.png`
- `results/deseq2/figures/checkpoints_heatmap.png`

### 📎 Rule: `multiqc`
Collects QC and alignment statistics into one HTML report.

Output:
- `results/multiqc/multiqc_report.html`

---

## 🧠 How to interpret the outputs

### `results/deseq2/checkpoints_summary.tsv`
This is the key file for testing the central hypothesis.

It summarizes, for each checkpoint gene and each cell line:
- `log2FoldChange`
- `pvalue`
- `padj`
- support for the expected **downregulation after PU001**

### `results/deseq2/DE_*`
Full DE tables for downstream interpretation, pathway analysis, and figure generation.

### `results/deseq2/figures/checkpoints_heatmap.png`
A checkpoint-focused heatmap based on VST-normalized expression.

### `results/multiqc/multiqc_report.html`
A unified quality-control report.

---

## 🧬 Recommended interpretation strategy

To support the hypothesis, prioritize genes showing:

- **negative log2FoldChange**
- **adjusted p-value (`padj`) < 0.05**
- consistent downregulation across **KG-1**, **Mono-Mac-1**, and **THP-1**

Pay particular attention to:
- **CD274 (PD-L1)**
- **HAVCR2 (TIM-3)**
- **PDCD1**
- **PDCD1LG2**
- **CTLA4**
- **LAG3**
- **TIGIT**
- **BTLA**
- **VSIR**
- **CD80**
- **CD86**
- **CD276**
- **VTCN1**
- **IDO1**

---

## 🤝 Acknowledgments

This repository was developed for an academic bioinformatics project focused on **AKT1 inhibition in MDS-related myeloid cell lines** and its relationship to **immune checkpoint regulation**.

---

## 👤 Author

Maintained by: **Alexey Demaev**  
Affiliation: **ITMO University** 
