# Unified KHV Pangenome & Metagenomic Pipeline

A comprehensive Python pipeline for quality control, pangenome construction, and metagenomic analysis of viral genomes (KHV - Koi Herpesvirus).

This pipeline integrates multiple bioinformatics tools into a unified, modular workflow with support for:
- Quality assessment (QUAST, BUSCO)
- Contamination detection (Kraken2)
- Pangenome construction (PGGB, Minigraph-Cactus)
- Variant calling and analysis
- PanSN format support for reference-based pangenomes

---

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Pipeline Features](#pipeline-features)
4. [Usage Examples](#usage-examples)
5. [Output Structure](#output-structure)
6. [Troubleshooting](#troubleshooting)
7. [Documentation](#documentation)

---

## Installation

### Step 1: Install Miniconda

If you don't have Miniconda installed:

```bash
# Download Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run installer
bash Miniconda3-latest-Linux-x86_64.sh

# Follow the prompts, then restart your terminal
```

### Step 2: Create Conda Environment

Create a new conda environment with all required tools:

```bash
# Create environment named 'BILL' with Python 3.9
conda create -n BILL python=3.9

# Activate the environment
conda activate BILL

# Install bioinformatics tools
conda install -c bioconda -c conda-forge \
    quast \
    busco \
    pggb \
    odgi \
    vg \
    bcftools \
    samtools

# Optional: Install Kraken2 for contamination detection
conda install -c bioconda kraken2

# Optional: Install additional tools
conda install -c bioconda mafft fasttree
```

### Step 3: Verify Installation

Check that all tools are installed:

```bash
# Test the installation
python3 test_portability.py

# Should show "EXCELLENT - Le pipeline est prêt à l'emploi!"
```

### Step 4: Download BUSCO Database (First Time Only)

```bash
# BUSCO will auto-download lineages on first use
# Or manually download:
busco --download viruses_odb10
```

---

## Quick Start

### Method 1: Automatic Mode (Recommended)

Run the entire pipeline automatically:

```bash
# 1. Activate conda environment
conda activate BILL

# 2. Edit configuration in run_all_steps.sh
nano run_all_steps.sh
# Update CLEANED_GENOMES and REFERENCE_FASTA paths

# 3. Run the complete pipeline
bash run_all_steps.sh
```

### Method 2: Interactive Mode

Run with user prompts for each optional step:

```bash
# 1. Activate conda environment
conda activate BILL

# 2. Edit configuration in run_full_pipeline.sh
nano run_full_pipeline.sh
# Update genome paths in configuration section

# 3. Run interactively
bash run_full_pipeline.sh
```

### Method 3: Custom Steps

Run specific pipeline steps:

```bash
# List all available steps
python3 unified_pipeline.py --list-steps

# Run specific steps
python3 unified_pipeline.py \
  --genomes sample1:data/sample1.fasta sample2:data/sample2.fasta \
  --pangenome-reference-fasta data/reference.fasta \
  --steps pangenome-pggb \
  --pggb-threads 8
```

---

## Pipeline Features

### Quality Control & Assessment

| Step | Tool | Description |
|------|------|-------------|
| pre-qc | Custom | Assembly statistics before cleaning |
| post-qc | Custom | Assembly statistics after cleaning |
| qc-merge | Custom | Compare pre/post QC metrics |
| quality-flags | Custom | Flag genomes based on expected characteristics |

### Genome Quality Evaluation

| Step | Tool | Description |
|------|------|-------------|
| quast-pre | QUAST | Structural assessment before cleaning |
| quast-post | QUAST | Structural assessment after cleaning |
| busco-pre | BUSCO | Completeness assessment before cleaning |
| busco-post | BUSCO | Completeness assessment after cleaning |

### Contamination Detection

| Step | Tool | Description |
|------|------|-------------|
| kraken | Kraken2 | Taxonomic classification and contamination detection |

### Pangenome Construction

| Step | Tool | Description |
|------|------|-------------|
| pangenome-pggb | PGGB | Pangenome graph with PanSN format support |
| pangenome-cactus | Minigraph-Cactus | Alternative pangenome construction |

### Reporting

| Step | Description |
|------|-------------|
| summary | Generate global summary report |

---

## Usage Examples

### Example 1: Complete Pipeline with All Steps

```bash
conda activate BILL

python3 unified_pipeline.py \
  --genomes \
    p15_1:clean_fasta/p15_1_clean.fasta \
    p15_6:clean_fasta/p15_6_clean.fasta \
    p90_2:clean_fasta/p90_2_clean.fasta \
    p90_6:clean_fasta/p90_6_clean.fasta \
  --pangenome-reference-fasta clean_fasta/ref/KHV-U_trunc.fasta \
  --all \
  --quast-threads 8 \
  --busco-threads 8 \
  --pangenome-threads 8 \
  --pggb-threads 8
```

### Example 2: Quality Control Only

```bash
conda activate BILL

python3 unified_pipeline.py \
  --genomes sample1:data/sample1.fasta sample2:data/sample2.fasta \
  --steps pre-qc,quast-pre,busco-pre,summary \
  --quast-threads 8 \
  --busco-threads 8 \
  --busco-lineage viruses_odb10
```

### Example 3: Pangenome Construction Only (PGGB)

```bash
conda activate BILL

python3 unified_pipeline.py \
  --genomes \
    sample1:clean_fasta/sample1_clean.fasta \
    sample2:clean_fasta/sample2_clean.fasta \
  --pangenome-reference-fasta clean_fasta/ref/reference.fasta \
  --steps pangenome-pggb \
  --pggb-segment 3000 \
  --pggb-pid 95 \
  --pggb-passes 10 \
  --pggb-threads 8
```

### Example 4: Using run_all_steps.sh Script

The easiest way to run the complete pipeline:

```bash
# 1. Activate conda environment
conda activate BILL

# 2. Edit configuration at the top of run_all_steps.sh
nano run_all_steps.sh

# Update these variables:
#   CLEANED_GENOMES - paths to your genome files
#   REFERENCE_FASTA - path to reference genome
#   THREADS - number of CPU threads to use
#   BUSCO_LINEAGE - appropriate lineage for your organism
#   EXPECTED_SIZE - expected genome size in bp
#   EXPECTED_GC - expected GC content percentage

# 3. Run the script
bash run_all_steps.sh

# Or with a custom conda environment name
bash run_all_steps.sh my_env_name
```

Configuration example in `run_all_steps.sh`:

```bash
# Input genomes (cleaned/assembled)
CLEANED_GENOMES=(
    "sample1:clean_fasta/sample1_clean.fasta"
    "sample2:clean_fasta/sample2_clean.fasta"
    "sample3:clean_fasta/sample3_clean.fasta"
)

# Reference genome for pangenome
REFERENCE_FASTA="clean_fasta/ref/reference.fasta"

# Number of threads (auto-detect or set manually)
THREADS=8

# BUSCO lineage (adjust for your organism)
BUSCO_LINEAGE="viruses_odb10"

# Expected genome characteristics
EXPECTED_SIZE=295146
EXPECTED_GC=59.2
```

---

## Output Structure

The pipeline creates the following output structure:

```
results/
├── reports/
│   ├── pre_decontam_qc_metrics.tsv     # Pre-QC statistics
│   ├── post_decontam_qc_metrics.tsv    # Post-QC statistics
│   ├── qc_compare_pre_post.tsv         # QC comparison
│   ├── quality_flags.tsv               # Quality assessment flags
│   └── global_summary.tsv              # Combined summary report
├── kraken/
│   ├── sample_kraken.out               # Kraken2 classifications
│   ├── sample_kraken_report.txt        # Kraken2 report
│   └── clean/summary_samples.tsv       # Contamination summary
├── quast_pre/                          # QUAST results (pre-cleaning)
│   └── report.html
├── quast_post/                         # QUAST results (post-cleaning)
│   └── report.html
├── busco_pre/                          # BUSCO results (pre-cleaning)
│   └── busco_sample_pre/
├── busco_post/                         # BUSCO results (post-cleaning)
│   └── busco_sample_post/
├── pggb/                               # PGGB pangenome
│   ├── input_genomes_pansn.fa          # Combined input in PanSN format
│   ├── *.gfa                           # Pangenome graph
│   ├── *.vcf                           # Variant calls
│   └── *.png                           # Visualizations
├── minigraph_cactus/                   # Cactus pangenome
│   ├── *.hal                           # HAL alignment
│   └── *.gfa                           # Graph format
└── validation/                         # Validation results
```

---

## Troubleshooting

### Common Issues and Solutions

| Problem | Cause | Solution |
|---------|-------|----------|
| "conda: command not found" | Conda not in PATH | Add to PATH: `export PATH="$HOME/miniconda3/bin:$PATH"` |
| "python3: command not found" | Not in conda environment | Run `conda activate BILL` |
| "pggb: command not found" | Tool not installed | Install: `conda install -c bioconda pggb` |
| QUAST not found | Binary missing | Install: `conda install -c bioconda quast` |
| BUSCO fails | Lineage not downloaded | Download: `busco --download viruses_odb10` |
| PGGB execution fails | Genomes too small/fragmented | Use complete, high-quality assemblies |
| Kraken2 step skipped | Database not configured | Set `KRAKEN2_DB` environment variable |
| Permission denied on scripts | Not executable | Run: `chmod +x run_*.sh` |

### Checking Installation

```bash
# Verify tools are available
conda activate BILL
which python3 pggb quast.py busco bcftools samtools

# Test portability
python3 test_portability.py

# List available pipeline steps
python3 unified_pipeline.py --list-steps
```

### Cache and Re-runs

The pipeline implements intelligent caching:

- **QUAST**: Skips if results exist (force with `--force-quast`)
- **BUSCO**: Skips if results exist (force with `--force-busco`)
- **Kraken2**: Skips if summary exists

To force a complete re-run:

```bash
# Remove results directory
rm -rf results/

# Run pipeline again
bash run_all_steps.sh
```

---

## Documentation

Complete documentation is available:

- **SETUP_GUIDE.md** - Detailed installation and configuration
- **ENVIRONMENT_SETUP.md** - Environment variable configuration
- **TRANSFER_GUIDE.md** - Moving pipeline to a new system
- **PANGENOME_USAGE.md** - Pangenome construction guide
- **PANSN_FORMAT_GUIDE.md** - PanSN naming format explained
- **PIPELINE_WORKFLOW.md** - Workflow and step descriptions
- **PORTABILITY_UPDATE.md** - Portability changes documentation
- **PORTABILITY_CHECK.md** - Portability verification results
- **QUICK_REFERENCE.sh** - Quick help (run: `bash QUICK_REFERENCE.sh`)

### Getting Help

```bash
# View pipeline help
python3 unified_pipeline.py --help

# List all available steps
python3 unified_pipeline.py --list-steps

# Quick reference guide
bash QUICK_REFERENCE.sh

# Test installation
python3 test_portability.py
```

---

## System Requirements

### Minimum Requirements

- **OS**: Linux (tested on Ubuntu/Debian)
- **CPU**: 4+ cores
- **RAM**: 8 GB
- **Storage**: 20 GB free space

### Recommended Requirements

- **CPU**: 8+ cores
- **RAM**: 16-32 GB
- **Storage**: 50-100 GB free space

---

## Citation

If you use this pipeline in your research, please cite the relevant tools:

- **QUAST**: Gurevich et al. (2013)
- **BUSCO**: Simão et al. (2015)
- **PGGB**: Garrison et al. (2023)
- **Kraken2**: Wood et al. (2019)
- **Minigraph-Cactus**: Armstrong et al. (2020)

---

## License

This pipeline is for academic and research use.

---

## Support

For issues or questions:

1. Check the documentation in this repository
2. Run `python3 test_portability.py` to verify installation
3. Review error messages carefully
4. Ensure all dependencies are installed

---

**Last updated**: October 2, 2025
