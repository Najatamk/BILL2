# 🧬 KHV Pangenome Analysis Pipeline - Complete Workflow

## 📋 Overview
This pipeline performs comprehensive pangenome analysis of Koi Herpesvirus (KHV) genomes with focus on thermal adaptation. It includes contamination detection, pangenome construction, variant analysis, and advanced comparative genomics.

## 🚀 Quick Start
```python
# Initialize pipeline
pipeline = PipelineExecutor()
pipeline.initialize()

# Run complete analysis
run_full_pipeline()

# Or execute step-by-step
pipeline.execute_step("validate_input_genomes_rigorously")
```

---

## 📊 Complete Pipeline Workflow

### **PHASE 1: SETUP & VALIDATION** 🔧

#### 1.1 Initial Setup
- **Function**: `validate_paths()`
- **Purpose**: Validate file paths and create directory structure
- **Outputs**: 
  - `results/` directory structure
  - Path validation report

#### 1.2 Dependency Verification
- **Function**: `check_dependencies()`
- **Purpose**: Check for required bioinformatics tools
- **Tools Verified**:
  - Python libraries (BioPython, pandas, matplotlib)
  - QUAST, BUSCO, FastANI
  - Kraken2, BLAST, bcftools
  - Pangenome tools (PGGB, Cactus)
  - Workflow management (Snakemake)

#### 1.3 Genome Validation
- **Function**: `validate_input_genomes_rigorously()`
- **Purpose**: Comprehensive validation of input genomes
- **Checks**:
  - File format (FASTA)
  - Sequence integrity
  - Expected genome size (~295 kb for KHV)
  - GC content (~59.2%)
  - N content analysis

---

### **PHASE 2: QUALITY ASSESSMENT** 📊

#### 2.1 Genome Quality Analysis
- **Function**: `run_quast_analysis()`
- **Purpose**: Comprehensive genome quality assessment
- **Metrics**:
  - Assembly statistics
  - N50, L50 values
  - GC content distribution
  - Contig analysis
- **Output**: `results/quast_analysis/`

#### 2.2 Completeness Assessment
- **Function**: `run_comprehensive_busco_analysis()`
- **Purpose**: Assess genome completeness using viral lineages
- **Analysis**:
  - Core gene completeness
  - Missing genes identification
  - Fragmented genes analysis
- **Output**: `results/busco_analysis/`

#### 2.3 Similarity Analysis
- **Function**: `run_fastani_comparison()`
- **Purpose**: Calculate Average Nucleotide Identity (ANI)
- **Analysis**:
  - Pairwise genome comparison
  - Strain similarity assessment
  - Phylogenetic distance estimation
- **Output**: ANI matrix and plots

---

### **PHASE 3: CONTAMINATION DETECTION & REMOVAL** 🧽

#### 3.1 Comprehensive Contamination Check
- **Function**: `run_comprehensive_contamination_check()`
- **Tool**: Kraken2 + KrakenTools
- **Process**:
  1. Taxonomic classification of all sequences
  2. Identification of viral vs non-viral content
  3. Detection of host, bacterial, and other contamination
  4. Extraction of clean viral sequences

#### 3.2 Alternative Validation Methods
- **Function**: `run_blast_contamination_check()`
- **Tool**: BLAST against NCBI
- **Purpose**: Independent validation of contamination detection
- **Target Organisms**:
  - KHV-specific sequences
  - Other herpesviruses
  - Host sequences (fish)
  - Bacterial contamination

#### 3.3 Sequence Extraction & Cleaning
- **Function**: `extract_viral_sequences()`
- **Tool**: KrakenTools (`extract_kraken_reads.py`)
- **Target TaxIDs**:
  - 10292 (Herpesviridae)
  - 548681 (Cyprinid herpesvirus 3 - KHV)
  - 10376 (Varicellovirus)
  - 28285 (Cytomegalovirus)
- **Output**: `data/{genome_name}_cleaned.fasta`

#### 3.4 Contamination Reporting
- **Function**: `create_contamination_summary_report()`
- **Output**: Detailed contamination report with recommendations

---

### **PHASE 4: PANGENOME CONSTRUCTION** 🧬

#### 4.1 Minigraph-Cactus Pipeline
- **Function**: `run_cactus_pipeline_locally()` or `run_cactus_pipeline_slurm()`
- **Purpose**: Create reference-guided pangenome
- **Process**:
  1. Progressive alignment using Cactus
  2. Graph construction with Minigraph
  3. HAL format generation
  4. VCF variant calling
- **Outputs**:
  - `khv_thermal.hal` (alignment)
  - `khv_thermal.gfa` (graph)
  - `khv_thermal.vcf` (variants)
  - `khv_thermal.paf` (pairwise alignments)

#### 4.2 PGGB Pipeline (Alternative)
- **Function**: `run_pggb_pipeline()`
- **Purpose**: Reference-free pangenome construction
- **Process**:
  1. All-vs-all sequence alignment
  2. Graph induction and optimization
  3. Variant detection
- **Outputs**: Similar format to Cactus

---

### **PHASE 5: VARIANT ANALYSIS** 🔬

#### 5.1 Thermal Adaptation Analysis
- **Function**: `analyze_khv_thermal_adaptation_variants()`
- **Purpose**: Identify variants associated with thermal stress
- **Key Genes Analyzed**:
  - ORF25 (Major capsid protein)
  - ORF36 (Large subunit terminase)
  - ORF56 (DNA polymerase)
  - ORF57 (Single-strand DNA binding)
  - ORF81 (Ribonucleotide reductase)
  - ORF134 (Helicase)
  - ORF136 (Thymidine kinase)

#### 5.2 Variant Hotspot Analysis
- **Function**: `analyze_variant_hotspots()`
- **Purpose**: Identify regions of high variability
- **Process**:
  1. Genome segmentation (5kb windows)
  2. Variant density calculation
  3. Hotspot identification
  4. Functional annotation correlation

#### 5.3 VCF Processing & Analysis
- **Function**: `run_variant_analysis()`
- **Purpose**: Comprehensive variant characterization
- **Analysis**:
  - SNP/Indel classification
  - Variant effect prediction
  - Population genetics metrics
  - Functional impact assessment

---

### **PHASE 6: ADVANCED ANALYSES** 🔬

#### 6.1 Phylogenetic Analysis
- **Function**: `run_phylogenetic_analysis()`
- **Purpose**: Evolutionary relationship reconstruction
- **Methods**:
  - Multiple sequence alignment
  - Phylogenetic tree construction
  - Bootstrap analysis
  - Molecular clock analysis

#### 6.2 Functional Impact Prediction
- **Function**: `FunctionalImpactPredictor`
- **Purpose**: Predict functional consequences of variants
- **Tools**:
  - SnpEff annotation
  - PROVEAN impact prediction
  - Conservation analysis
  - Structural impact assessment

#### 6.3 Synteny Analysis
- **Function**: `SyntenyAnalyzer`
- **Purpose**: Analyze structural variations and gene order
- **Analysis**:
  - Collinearity assessment
  - Inversion detection
  - Translocation identification
  - Gene duplication/loss

#### 6.4 Selection Pressure Analysis
- **Function**: `SelectionPressureAnalyzer`
- **Purpose**: Detect signatures of natural selection
- **Metrics**:
  - dN/dS ratios
  - Tajima's D
  - McDonald-Kreitman test
  - Positive selection detection

---

### **PHASE 7: QUALITY CONTROL** ✅

#### 7.1 Comprehensive Quality Control
- **Function**: `run_comprehensive_quality_control()`
- **Purpose**: Validate all pipeline outputs
- **Validation**:
  - PAF file integrity
  - GFA graph validation
  - VCF format compliance
  - Statistical consistency checks

#### 7.2 Output Validation
- **Function**: `execute_quality_control_pipeline()`
- **Purpose**: Systematic validation of results
- **Checks**:
  - File completeness
  - Format validation
  - Content verification
  - Cross-method consistency

---

### **PHASE 8: VISUALIZATION & REPORTING** 📊

#### 8.1 Statistical Visualizations
- **Function**: `create_statistical_plots()`
- **Purpose**: Generate comprehensive plots
- **Plots**:
  - Variant distribution
  - Phylogenetic trees
  - Synteny plots
  - Selection pressure maps

#### 8.2 Quality Visualizations
- **Function**: `create_quality_visualizations()`
- **Purpose**: QC plots and summaries
- **Outputs**:
  - Assembly quality plots
  - Contamination reports
  - Coverage statistics
  - Comparison matrices

#### 8.3 Comparative Analysis
- **Function**: `create_variant_comparison_plots()`
- **Purpose**: Method comparison and validation
- **Analysis**:
  - Cactus vs PGGB comparison
  - Concordance analysis
  - Performance benchmarking

#### 8.4 Summary Reports
- **Function**: `create_summary_table()`
- **Purpose**: Generate final analysis report
- **Outputs**:
  - `results/reports/pangenome_comparison_summary.csv`
  - Comprehensive analysis summary
  - Method recommendations

---

### **PHASE 9: BENCHMARKING & PERFORMANCE** ⚡

#### 9.1 Performance Monitoring
- **Class**: `PerformanceMonitor`
- **Purpose**: Track resource usage
- **Metrics**:
  - CPU usage
  - Memory consumption
  - Disk I/O
  - Execution time

#### 9.2 Method Comparison
- **Class**: `BenchmarkComparison`
- **Purpose**: Compare pangenome construction methods
- **Analysis**:
  - Speed comparison
  - Resource efficiency
  - Output quality
  - Scalability assessment

---

### **PHASE 10: VALIDATION & REPRODUCIBILITY** 🔄

#### 10.1 Validation Framework
- **Function**: `create_validation_framework()`
- **Purpose**: Ensure result reliability
- **Components**:
  - Automated testing
  - Reproducibility checks
  - Parameter validation
  - Output verification

#### 10.2 Reproducibility Documentation
- **Function**: `create_validation_checklist()`
- **Purpose**: Document analysis parameters
- **Outputs**:
  - Software versions
  - Parameter settings
  - Environment info
  - Execution logs

---

## 📂 Key Output Directories

```
results/
├── benchmarks/           # Performance comparisons
├── minigraph_cactus/     # Cactus pipeline outputs
├── pggb/                 # PGGB pipeline outputs
├── plots/                # All visualization outputs
├── quast_analysis/       # Genome quality reports
├── reports/              # Summary reports
├── validation/           # Validation and testing
└── vcf_comparison/       # Variant analysis results
```

## 🛠️ Required Dependencies

### Core Tools
- Python 3.8+ with BioPython, pandas, matplotlib
- QUAST, BUSCO, FastANI
- Kraken2, KrakenTools, BLAST
- bcftools, samtools, vcftools

### Pangenome Construction
- Cactus + Minigraph
- PGGB (PanGenome Graph Builder)
- odgi, vg, wfmash

### Workflow Management
- Snakemake
- Singularity (for containerized execution)

## 🚀 Execution Options

### 1. Interactive Execution
```python
pipeline = PipelineExecutor()
pipeline.initialize()
pipeline.list_available_steps()
pipeline.execute_step("step_name")
```

### 2. Complete Workflow
```python
run_full_pipeline()
```

### 3. Individual Analyses
```python
run_comprehensive_contamination_check()
run_cactus_pipeline_locally()
run_comprehensive_advanced_analysis()
```

### 4. Cluster Execution
```bash
sbatch results/minigraph_cactus/submit_khv_analysis.sh
```

## 📊 Expected Results

1. **Clean genomes** free of contamination
2. **High-quality pangenome** with comprehensive variant catalog
3. **Thermal adaptation variants** with functional annotations
4. **Phylogenetic relationships** and evolutionary insights
5. **Comprehensive quality reports** and validation
6. **Publication-ready visualizations** and tables

---

## 📞 Support

For questions or issues, refer to:
- `QUICK_START_GUIDE.md` - Step-by-step execution guide
- `results/validation/validation_checklist.md` - Quality control checklist
- Pipeline logs and error messages for troubleshooting

---

**Pipeline Version**: Advanced KHV Pangenome Analysis v2.0
**Last Updated**: October 2025
**Authors**: BILL2 Development Team