# Prokaryotic Genome Assembly & AMR Analysis Pipeline

**Version:** 2.0.0  
**Script:** `prokaryotic_amr_pipeline.sh`  
**Language:** Bash (requires Bash ≥ 4.2)

---

## Overview

A four-stage shell pipeline for prokaryotic genome assembly and antimicrobial resistance (AMR) characterisation from bacterial sequencing data. The pipeline supports three sequencing modes — short reads (Illumina), long reads (Nanopore/PacBio), and hybrid (Illumina + Nanopore/PacBio) — and produces a polished assembly, a full AMRFinder+ gene report, and phenotypic resistance predictions from the AMR Rules tool.

```
Reads  ──►  Stage 1: QC & Trimming
            Stage 2: Assembly
            Stage 3: AMR Gene Detection   (AMRFinder+)
            Stage 4: Phenotypic Prediction (AMR Rules)  ──►  Summary report
```

---

## Modes at a glance

| Mode | Input | Assembler | Polisher | Typical use case |
|---|---|---|---|---|
| `short` | Illumina PE | Unicycler | — | Standard WGS on short-read platforms |
| `long` | ONT / PacBio | Flye | Medaka | Nanopore-only sequencing runs |
| `hybrid` | Illumina PE + ONT/PacBio | Unicycler (hybrid) | — | Long reads for scaffolding, Illumina for accuracy |

---

## Dependencies

### Conda environments

AMR Rules **must** be installed in a dedicated conda environment named `amrrules_beta`. The pipeline calls it with `conda run`, so the calling shell does **not** need to have the environment activated.

```bash
conda create -n amrrules_beta python=3.10
conda activate amrrules_beta
pip install amr-rules
conda deactivate
```

The `conda` executable must be reachable in `PATH` (i.e. conda must be initialised in your shell profile).

### Tools by mode

The table below lists every tool that must be available in `PATH` (outside of `amrrules_beta`). Install via conda/mamba or the method recommended by each tool's documentation.

| Tool | Version tested | Mode | Purpose |
|---|---|---|---|
| `fastqc` | ≥ 0.12 | short, hybrid | Raw and trimmed read QC |
| `trimmomatic` | ≥ 0.39 | short, hybrid | Illumina adapter trimming and quality filtering |
| `unicycler` | ≥ 0.5 | short, hybrid | SPAdes-based short-read / hybrid assembly |
| `nanostat` | ≥ 1.6 | long, hybrid | Long-read QC statistics |
| `nanofilt` | ≥ 2.8 | long, hybrid | Long-read quality and length filtering |
| `flye` | ≥ 2.9 | long | Long-read de novo assembly |
| `medaka_consensus` | ≥ 1.8 | long | Nanopore consensus polishing |
| `amrfinder` | ≥ 3.11 | all | AMR gene detection (NCBI AMRFinder+) |
| `minimap2` | ≥ 2.24 | hybrid | Required internally by Unicycler hybrid mode |

> **Tip — suggested conda setup for non-AMR-Rules tools:**
> ```bash
> mamba create -n amr_pipeline -c bioconda -c conda-forge \
>     fastqc trimmomatic unicycler \
>     nanostat nanofilt flye medaka \
>     ncbi-amrfinderplus minimap2
> conda activate amr_pipeline
> ```

---

## Installation

```bash
# Clone or download the script
git clone https://github.com/your-org/prokaryotic-amr-pipeline.git
cd prokaryotic-amr-pipeline

# Make executable
chmod +x prokaryotic_amr_pipeline.sh
```

No compilation is required.

---

## Usage

```
prokaryotic_amr_pipeline.sh -m <short|long|hybrid> [options]
```

### Required flags (all modes)

| Flag | Description |
|---|---|
| `-m` | Assembly mode: `short`, `long`, or `hybrid` |
| `-s` | Sample name — used as a prefix for all output files |
| `-o` | Output directory (created if it does not exist) |

### Read input flags

| Flag | Required for | Description |
|---|---|---|
| `-1` | `short`, `hybrid` | Illumina forward reads (FASTQ / FASTQ.gz) |
| `-2` | `short`, `hybrid` | Illumina reverse reads (FASTQ / FASTQ.gz) |
| `-l` | `long`, `hybrid` | Long reads (FASTQ / FASTQ.gz) |

### Optional flags

| Flag | Default | Description |
|---|---|---|
| `-t` | `8` | Number of CPU threads |
| `-g` | _(none)_ | Organism name for AMRFinder+ (e.g. `Escherichia`, `Klebsiella`, `Salmonella`). Enables organism-specific point mutation detection. |
| `-a` | `NexteraPE-PE.fa` | Trimmomatic adapter FASTA file |
| `-q` | `8` | NanoFilt minimum Phred quality score |
| `-L` | `500` | NanoFilt minimum read length (bp) |
| `-x` | _(none, **required** for `long`/`hybrid`)_ | Estimated genome size for Flye (e.g. `5m`, `4.5m`, `2.8m`) |
| `-r` | `nano-raw` | Flye read type: `nano-raw`, `nano-hq`, `nano-corr`, `pacbio-raw`, `pacbio-corr`, `pacbio-hifi` |
| `-h` | — | Print help and exit |

### Medaka model

The Medaka polishing model (long mode only) is not a command-line flag. Set it via the environment variable `MEDAKA_MODEL` before calling the script. If unset, it defaults to `r941_min_sup_g507`.

```bash
export MEDAKA_MODEL="r1041_e82_400bps_sup_v5.0.0"
bash prokaryotic_amr_pipeline.sh -m long ...
```

Refer to the [Medaka documentation](https://github.com/nanoporetech/medaka#models) for the correct model for your flowcell and basecalling version.

---

## Examples

### Short-read assembly (Illumina WGS)

```bash
bash prokaryotic_amr_pipeline.sh \
    -m short \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -s KPNEUMO_001 \
    -o results/KPNEUMO_001 \
    -t 16 \
    -g Klebsiella
```

### Long-read assembly (Nanopore only)

```bash
bash prokaryotic_amr_pipeline.sh \
    -m long \
    -l sample_ont.fastq.gz \
    -s KPNEUMO_001 \
    -o results/KPNEUMO_001 \
    -t 16 \
    -x 5m \
    -r nano-hq \
    -g Klebsiella
```

### Hybrid assembly (Illumina + Nanopore)

```bash
bash prokaryotic_amr_pipeline.sh \
    -m hybrid \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -l sample_ont.fastq.gz \
    -s KPNEUMO_001 \
    -o results/KPNEUMO_001 \
    -t 16 \
    -x 5m \
    -g Klebsiella
```

### PacBio HiFi assembly

```bash
bash prokaryotic_amr_pipeline.sh \
    -m long \
    -l sample_hifi.fastq.gz \
    -s ECOLI_042 \
    -o results/ECOLI_042 \
    -t 16 \
    -x 5m \
    -r pacbio-hifi \
    -g Escherichia
```

---

## Pipeline stages

### Stage 1 — Quality control

**Short / hybrid mode (Illumina reads)**

1. `FastQC` is run on raw paired reads to generate per-base quality, GC content, and adapter content reports.
2. `Trimmomatic` (PE mode) removes adapter sequences and low-quality bases using a sliding-window filter (window 4, mean Q ≥ 15), hard-clips leading/trailing bases below Q3, and discards reads shorter than 50 bp.
3. `FastQC` is re-run on the trimmed paired output.

**Long / hybrid mode (Nanopore / PacBio reads)**

1. `NanoStat` produces a summary of read N50, mean quality, total bases, and read-length distribution.
2. `NanoFilt` filters reads below the specified minimum Q score (default Q8) and minimum length (default 500 bp), writing a compressed FASTQ.

### Stage 2 — Assembly

**Short mode** — `Unicycler` assembles trimmed Illumina paired reads using its SPAdes-based pipeline. Unicycler selects the best SPAdes k-mer graph, trims overlaps, and attempts to circularise replicons.

**Long mode** — `Flye` assembles NanoFilt-filtered long reads into a draft genome, producing contigs and an assembly graph. The draft is then polished with `Medaka` (`medaka_consensus`), which realigns the original reads to the draft assembly and calls a consensus using a neural network model. The final polished FASTA is used for all downstream steps.

**Hybrid mode** — `Unicycler` is run in hybrid mode, supplying both trimmed Illumina paired reads (`-1/-2`) and NanoFilt-filtered long reads (`-l`). Unicycler uses the long reads to bridge SPAdes assembly graph gaps and resolve repeats, while the Illumina reads correct base-level errors. No separate polishing step is needed.

In all modes, the final assembly is copied to `<outdir>/03_assembly/<sample>_assembly.fasta`.

### Stage 3 — AMR gene detection (AMRFinder+)

`AMRFinder+` is run in nucleotide mode (`--nucleotide`) against the assembly FASTA with the `--plus` flag, which adds stress, virulence, and biocide gene detection on top of the core AMR gene set. If an organism name is provided with `-g`, species-specific point mutation rules are also applied.

The AMRFinder+ database is updated automatically at the start of this stage. If the update fails (e.g. offline HPC node), a warning is printed and the pipeline continues with the existing local database.

Output: a TSV report (`<sample>_amrfinder.tsv`) with one row per detected element, including gene name, sequence name, coverage, identity, element type, and predicted drug class.

### Stage 4 — Phenotypic AMR interpretation (AMR Rules)

`AMR Rules` translates the AMRFinder+ gene-detection output into predicted phenotypic resistance profiles, applying curated rule sets that map genotype to clinical breakpoint categories (S/I/R).

Because AMR Rules is only available inside the `amrrules_beta` conda environment, this stage uses `conda run -n amrrules_beta`, which spawns the tool as a subprocess inside that environment without requiring the calling shell to activate it. This approach is safe with `set -euo pipefail` and compatible with non-interactive schedulers (SLURM, PBS, Snakemake, Nextflow).

Before execution, the pipeline verifies that `amr_rules` (or the `amr_rules` Python module) is importable inside `amrrules_beta`, and exits with a descriptive error if not found.

---

## Output structure

```
<outdir>/
├── 01_qc/
│   ├── <sample>_R1_fastqc.html          # FastQC report — raw R1 (short/hybrid)
│   ├── <sample>_R2_fastqc.html          # FastQC report — raw R2 (short/hybrid)
│   ├── <sample>_R1_paired_fastqc.html   # FastQC report — trimmed R1 (short/hybrid)
│   ├── <sample>_R2_paired_fastqc.html   # FastQC report — trimmed R2 (short/hybrid)
│   └── <sample>_nanostat_raw.txt        # NanoStat summary (long/hybrid)
│
├── 02_trimmed/
│   ├── <sample>_R1_paired.fastq.gz      # Trimmomatic output — R1 paired (short/hybrid)
│   ├── <sample>_R2_paired.fastq.gz      # Trimmomatic output — R2 paired (short/hybrid)
│   ├── <sample>_R1_unpaired.fastq.gz    # Trimmomatic output — R1 singletons
│   ├── <sample>_R2_unpaired.fastq.gz    # Trimmomatic output — R2 singletons
│   └── <sample>_long_filtered.fastq.gz  # NanoFilt output (long/hybrid)
│
├── 03_assembly/
│   ├── unicycler/                        # Full Unicycler output dir (short mode)
│   ├── unicycler_hybrid/                 # Full Unicycler output dir (hybrid mode)
│   ├── flye/                             # Full Flye output dir (long mode)
│   ├── medaka/                           # Full Medaka output dir (long mode)
│   └── <sample>_assembly.fasta          # Final assembly used downstream (all modes)
│
├── 04_amrfinder/
│   └── <sample>_amrfinder.tsv           # AMRFinder+ tabular results
│
├── 05_amr_rules/
│   └── <sample>_amr_rules*              # AMR Rules phenotypic output
│
├── logs/
│   ├── <sample>_fastqc_raw.log
│   ├── <sample>_trimmomatic.log
│   ├── <sample>_fastqc_trimmed.log
│   ├── <sample>_nanostat.log
│   ├── <sample>_nanofilt.log
│   ├── <sample>_flye.log
│   ├── <sample>_medaka.log
│   ├── <sample>_unicycler.log           # or _unicycler_hybrid.log
│   ├── <sample>_amrfinder_update.log
│   ├── <sample>_amrfinder.log
│   └── <sample>_amr_rules.log
│
└── <sample>_pipeline_summary.txt        # Plain-text run summary
```

---

## Troubleshooting

**`conda environment 'amrrules_beta' not found`**  
The `amrrules_beta` environment does not exist. Create it and install `amr-rules` as described in the [Dependencies](#dependencies) section. Confirm with `conda env list`.

**`amr_rules not found inside conda env 'amrrules_beta'`**  
The environment exists but `amr-rules` is not installed inside it. Activate the environment manually (`conda activate amrrules_beta`) and run `pip install amr-rules`.

**`Flye assembly.fasta not found`**  
Flye failed to produce an assembly. Common causes: genome size estimate too far off (`-x`), read depth below ~10×, or corrupted input FASTQ. Inspect `logs/<sample>_flye.log` for details. Flye is more tolerant of uneven depth than most assemblers, but very low coverage (< 5–8×) assemblies are frequently incomplete.

**`Medaka consensus.fasta not found`**  
Medaka completed without writing output, which typically indicates a model mismatch. Set the correct model for your chemistry via `export MEDAKA_MODEL=<model>` before running the pipeline. Run `medaka tools list_models` inside the relevant conda environment to see available models.

**`Unicycler assembly.fasta not found`**  
Unicycler failed. Inspect `logs/<sample>_unicycler.log`. Unicycler requires paired reads with sufficient overlap; very low coverage (< 20×), highly fragmented reads, or missing `minimap2` (in hybrid mode) are the most common causes.

**`AMRFinder+ database update failed`**  
This is a non-fatal warning. The pipeline continues with the currently installed database. To update manually: `amrfinder --update`. On air-gapped systems, follow the [offline update instructions](https://github.com/ncbi/amr/wiki/Installing-AMRFinder#offline-install) from NCBI.

**Pipeline exits immediately with `set -euo pipefail`**  
Any non-zero exit from any command will abort the pipeline. Check the most recent log file in `logs/` for the actual error. Pipe failures (e.g. broken gzip stream in NanoFilt) are also caught.

---

## Notes on assembly mode selection

**Use `short`** when you have Illumina PE reads only. Unicycler typically produces near-complete assemblies for organisms with genomes under 10 Mb, though chromosomal repetitive regions (IS elements, rRNA operons) frequently remain unresolved.

**Use `long`** when you have Nanopore or PacBio reads only and want complete, circularised chromosomes and plasmids. Requires minimum ~20–30× depth for Flye; Medaka polishing requires the same reads used for assembly. Choose the Medaka model matching your flowcell chemistry and basecaller version.

**Use `hybrid`** when you have both data types. Unicycler's hybrid mode uses long reads to span repeats and circularise replicons while using Illumina reads to correct substitution and indel errors introduced by long-read basecalling. This typically yields the highest-quality assemblies and is particularly recommended for AMR studies, where accurate gene sequences are critical for AMRFinder+ point mutation calling.

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **FastQC** — Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **Trimmomatic** — Bolger AM et al. (2014). *Bioinformatics*, 30(15):2114–2120.
- **Unicycler** — Wick RR et al. (2017). *PLOS Computational Biology*, 13(6):e1005595.
- **NanoStat** — De Coster W & Rademakers R. (2023). *Bioinformatics*, 39(5):btad311.
- **NanoFilt** — De Coster W et al. (2018). *Bioinformatics*, 34(15):2666–2669.
- **Flye** — Kolmogorov M et al. (2019). *Nature Methods*, 16:1050–1054.
- **Medaka** — Oxford Nanopore Technologies. https://github.com/nanoporetech/medaka
- **AMRFinder+** — Feldgarden M et al. (2021). *Scientific Reports*, 11:12728.
- **AMR Rules** — Cite per the tool's own documentation / repository.

---

## Licence

This pipeline script is released under the MIT Licence. See `LICENSE` for details.
