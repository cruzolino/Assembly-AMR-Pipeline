# Prokaryotic Genome Assembly & AMR Analysis Pipeline

A bash pipeline for de novo assembly of prokaryotic genomes from paired-end Illumina reads, followed by antimicrobial resistance (AMR) gene detection and phenotypic interpretation.

---

## Overview

The pipeline runs three sequential stages:

```
Raw FASTQ reads
      │
      ▼
┌─────────────────────────────────────┐
│  Stage 1 — QC & Assembly            │
│  FastQC → Trimmomatic → FastQC      │
│              → Unicycler            │
└─────────────────────┬───────────────┘
                      │ assembly.fasta
                      ▼
┌─────────────────────────────────────┐
│  Stage 2 — AMRFinder+               │
│  AMR gene & mutation detection      │
└─────────────────────┬───────────────┘
                      │ amrfinder.tsv
                      ▼
┌─────────────────────────────────────┐
│  Stage 3 — AMR Rules                │
│  Phenotypic resistance prediction   │
└─────────────────────────────────────┘
```

---

## Dependencies

All tools must be available in your `$PATH` before running the pipeline.

| Tool | Version tested | Purpose | Install |
|------|---------------|---------|---------|
| FastQC | ≥ 0.11.9 | Read quality metrics | `conda install -c bioconda fastqc` |
| Trimmomatic | ≥ 0.39 | Adapter trimming & quality filtering | `conda install -c bioconda trimmomatic` |
| Unicycler | ≥ 0.5.0 | De novo assembly | `conda install -c bioconda unicycler` |
| AMRFinder+ | ≥ 3.11 | AMR gene & point mutation detection | `conda install -c bioconda ncbi-amrfinderplus` |
| AMR Rules | ≥ 1.0 | Phenotypic resistance interpretation | `pip install amr-rules` |
| Python 3 | ≥ 3.8 | Required by AMR Rules | system / conda |

### Recommended: Conda environment

```bash
conda create -n amr_pipeline \
    -c bioconda -c conda-forge \
    fastqc trimmomatic unicycler ncbi-amrfinderplus python=3.10

conda activate amr_pipeline
pip install amr-rules
```

### Adapter file

Trimmomatic requires an adapter FASTA file. Common options shipped with Trimmomatic:

- `NexteraPE-PE.fa` — Nextera XT libraries (default in this pipeline)
- `TruSeq3-PE-2.fa` — Illumina TruSeq libraries
- `TruSeq2-PE.fa` — older Illumina GA/HiSeq libraries

Locate your Trimmomatic installation's adapter directory with:

```bash
dirname $(which trimmomatic)/../share/trimmomatic/adapters/
```

---

## Usage

```bash
bash prokaryotic_amr_pipeline.sh \
    -1 <R1.fastq.gz> \
    -2 <R2.fastq.gz> \
    -s <sample_name> \
    -o <output_dir> \
    [options]
```

### Required arguments

| Flag | Description |
|------|-------------|
| `-1` | Forward reads (FASTQ or FASTQ.gz) |
| `-2` | Reverse reads (FASTQ or FASTQ.gz) |
| `-s` | Sample name — used for all output file naming |
| `-o` | Output directory (created if it does not exist) |

### Optional arguments

| Flag | Default | Description |
|------|---------|-------------|
| `-t` | `8` | Number of CPU threads |
| `-g` | _(none)_ | Organism name for AMRFinder+ species-specific detection (see list below) |
| `-a` | `NexteraPE-PE.fa` | Path to Trimmomatic adapter FASTA file |
| `-h` | — | Print help and exit |

### `-g` organism values (AMRFinder+)

Providing an organism enables point mutation detection for known chromosomal resistance mechanisms. Accepted values include:

`Acinetobacter_baumannii`, `Burkholderia_cepacia`, `Burkholderia_pseudomallei`, `Campylobacter`, `Clostridioides_difficile`, `Enterococcus_faecalis`, `Enterococcus_faecium`, `Escherichia`, `Klebsiella_oxytoca`, `Klebsiella_pneumoniae`, `Neisseria_gonorrhoeae`, `Neisseria_meningitidis`, `Pseudomonas_aeruginosa`, `Salmonella`, `Staphylococcus_aureus`, `Staphylococcus_pseudintermedius`, `Streptococcus_agalactiae`, `Streptococcus_pneumoniae`, `Streptococcus_pyogenes`, `Vibrio_cholerae`

For unlisted species, omit `-g` — gene-based detection still runs.

---

## Example

```bash
# E. coli isolate, 16 threads
bash prokaryotic_amr_pipeline.sh \
    -1 data/ECO_001_R1.fastq.gz \
    -2 data/ECO_001_R2.fastq.gz \
    -s ECO_001 \
    -o results/ECO_001 \
    -t 16 \
    -g Escherichia \
    -a /opt/trimmomatic/adapters/TruSeq3-PE-2.fa
```

---

## Output structure

```
<output_dir>/
├── 01_fastqc_raw/
│   ├── <sample>_R1_fastqc.html        # Raw read QC report (forward)
│   └── <sample>_R2_fastqc.html        # Raw read QC report (reverse)
│
├── 02_trimmomatic/
│   ├── <sample>_R1_paired.fastq.gz    # Trimmed forward reads (used downstream)
│   ├── <sample>_R2_paired.fastq.gz    # Trimmed reverse reads (used downstream)
│   ├── <sample>_R1_unpaired.fastq.gz  # Orphaned forward reads
│   └── <sample>_R2_unpaired.fastq.gz  # Orphaned reverse reads
│
├── 03_fastqc_trimmed/
│   ├── <sample>_R1_paired_fastqc.html # Post-trim QC report (forward)
│   └── <sample>_R2_paired_fastqc.html # Post-trim QC report (reverse)
│
├── 04_assembly/
│   ├── <sample>_assembly.fasta        # Final assembled contigs (used downstream)
│   ├── assembly.fasta                 # Original Unicycler output
│   └── assembly.gfa                   # Assembly graph
│
├── 05_amrfinder/
│   └── <sample>_amrfinder.tsv         # AMR gene/mutation report
│
├── 06_amr_rules/
│   └── <sample>_amr_rules.*           # Phenotypic resistance predictions
│
├── logs/
│   ├── <sample>_fastqc_raw.log
│   ├── <sample>_trimmomatic.log
│   ├── <sample>_fastqc_trimmed.log
│   ├── <sample>_unicycler.log
│   ├── <sample>_amrfinder_update.log
│   ├── <sample>_amrfinder.log
│   └── <sample>_amr_rules.log
│
└── <sample>_pipeline_summary.txt      # Run summary (inputs, contig count, AMR hits)
```

---

## Stage details

### Stage 1 — Quality control & assembly

**FastQC (raw)** generates per-base quality scores, GC content, adapter content, and duplication metrics for the raw reads. These reports serve as baseline documentation.

**Trimmomatic** runs in paired-end mode with the following defaults:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ILLUMINACLIP` | `2:30:10` | Adapter clipping (seed mismatches : palindrome threshold : simple clip threshold) |
| `LEADING` | `3` | Remove leading bases below quality 3 |
| `TRAILING` | `3` | Remove trailing bases below quality 3 |
| `SLIDINGWINDOW` | `4:15` | Cut when 4-base window average quality drops below 15 |
| `MINLEN` | `50` | Discard reads shorter than 50 bp after trimming |

These defaults are conservative and suitable for most Illumina short-read data. Adjust them directly in the script header if your library chemistry requires different thresholds.

**FastQC (trimmed)** repeats quality assessment on the paired trimmed reads. Compare the pre- and post-trimming reports (or aggregate them with [MultiQC](https://multiqc.info/)) to confirm adapter removal and quality improvement.

**Unicycler** assembles the trimmed paired reads using a short-read-only mode. It internally orchestrates SPAdes and Pilon to produce a polished, circularised assembly where possible. Only paired reads are passed to the assembler; orphaned reads from Trimmomatic are retained in the output directory but not used.

### Stage 2 — AMRFinder+

AMRFinder+ (NCBI's Antimicrobial Resistance Gene Finder) screens the assembly nucleotide FASTA against the NCBI Reference Gene Database. The pipeline:

- Runs `amrfinder --update` at the start of each run to pull the latest database version. If the update fails (e.g. offline environment), a warning is issued and the existing local database is used.
- Passes `--plus` to extend detection beyond AMR genes to include stress response, virulence, and biocide resistance elements.
- Passes `--organism` when the `-g` flag is provided, enabling point mutation screening for chromosomal resistance determinants (e.g. *gyrA*, *parC* mutations conferring fluoroquinolone resistance in *E. coli*).

The output TSV contains one row per detected element with columns for gene name, sequence name, coordinates, strand, element type, subtype, drug class, and resistance mechanism.

### Stage 3 — AMR Rules

AMR Rules applies a curated rule set to the AMRFinder+ TSV to translate detected genes and mutations into predicted phenotypic resistance profiles (e.g. predicted resistance to specific antibiotic classes per clinical breakpoint logic). The pipeline auto-detects whether AMR Rules is installed as a command-line tool (`amr_rules`) or as a Python module (`python3 -m amr_rules`).

---

## Troubleshooting

**Pipeline exits immediately with "Required tool not found in PATH"**
Activate your conda environment before running: `conda activate amr_pipeline`

**Trimmomatic fails with "Unable to detect quality encoding"**
Your reads may use Phred+33 or Phred+64 encoding. Add `-phred33` or `-phred64` to the Trimmomatic call in the script.

**Trimmomatic fails to find the adapter file**
Pass the full absolute path with `-a /path/to/adapters/NexteraPE-PE.fa`. Locate available adapter files with `find $(conda info --base) -name "*.fa" | grep -i adapter`.

**Unicycler produces 0 contigs or fails**
Check `logs/<sample>_unicycler.log`. Common causes: insufficient read depth (< 20×), very short reads after trimming, or SPAdes running out of memory. For low-depth samples consider reducing SPAdes k-mer sizes.

**AMRFinder+ update fails in an HPC environment**
Run `amrfinder --update` once on a login node with internet access before submitting the pipeline as a job. The database is cached locally and the pipeline will use the cached version.

**AMR Rules not found**
Install with `pip install amr-rules` inside the activated conda environment. Verify with `python3 -c "import amr_rules; print(amr_rules.__version__)"`.

---

## Citation

If you use this pipeline in published work, please cite the underlying tools:

- **FastQC** — Andrews S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **Trimmomatic** — Bolger AM, Lohse M, Usadel B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics*, 30(15):2114–2120.
- **Unicycler** — Wick RR, Judd LM, Gorrie CL, Holt KE. (2017). Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads. *PLOS Computational Biology*, 13(6):e1005595.
- **AMRFinder+** — Feldgarden M et al. (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Scientific Reports*, 11:12728.
- **AMR Rules** — cite per the tool's repository/publication.

---

## License

This pipeline script is released under the MIT License. The tools it wraps are subject to their own respective licenses.
