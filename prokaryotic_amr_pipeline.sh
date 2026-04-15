#!/usr/bin/env bash
# =============================================================================
# Prokaryotic Genome Assembly & AMR Analysis Pipeline
# Version: 2.0.0
# =============================================================================
# Stages:
#   1. Quality control (FastQC + Trimmomatic / NanoStat + NanoFilt)
#   2. Assembly  — short reads : Unicycler (SPAdes-based)
#                — long reads  : Flye
#                — hybrid      : Unicycler (-1/-2 + -l)
#   3. AMR gene detection   (AMRFinder+)
#   4. Phenotypic AMR interpretation (AMR Rules — run inside conda env
#                                      'amrrules_beta')
#
# Dependencies (short / hybrid):
#   fastqc, trimmomatic, unicycler, amrfinder
# Dependencies (long / hybrid):
#   nanostat, nanofilt, flye, medaka  (long); unicycler also needs minimap2
# AMR Rules (all modes):
#   conda env 'amrrules_beta' containing amr_rules
#
# Usage:
#   bash prokaryotic_amr_pipeline.sh \
#       -m <short|long|hybrid> \
#       [-1 reads_R1.fastq.gz]  \      # required for short / hybrid
#       [-2 reads_R2.fastq.gz]  \      # required for short / hybrid
#       [-l long_reads.fastq.gz]\      # required for long  / hybrid
#       -s <sample_name>        \
#       -o <output_dir>         \
#       [-t <threads>]          \
#       [-g <organism>]         \      # e.g. Escherichia (AMRFinder+)
#       [-q <min_quality>]      \      # NanoFilt min Q (default 8)
#       [-L <min_length>]       \      # NanoFilt min length (default 500)
#       [-x <genome_size>]      \      # Flye estimated genome size, e.g. 5m
#       [-r <flye_read_type>]         # nano-raw|nano-hq|nano-corr|pacbio-raw|
#                                     # pacbio-corr|pacbio-hifi (default nano-raw)
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
MODE=""
R1=""; R2=""; LONG_READS=""
THREADS=8
ORGANISM=""
# Illumina trimming
TRIMMOMATIC_ADAPTERS="NexteraPE-PE.fa"
MIN_LEN=50; LEADING=3; TRAILING=3; SLIDINGWINDOW="4:15"
# Nanopore / long-read filtering
NANO_MIN_Q=8
NANO_MIN_LEN=500
# Flye
GENOME_SIZE=""          # e.g. 5m — required for long / hybrid; fatal if absent
FLYE_READ_TYPE="nano-raw"
# AMR Rules conda environment
AMRRULES_ENV="amrrules_beta"

# ---------------------------------------------------------------------------
# Colour helpers
# ---------------------------------------------------------------------------
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'

log()  { echo -e "${CYAN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $*"; }
ok()   { echo -e "${GREEN}[OK]${NC} $*"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }
die()  { echo -e "${RED}[ERROR]${NC} $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
usage() {
cat <<EOF
${BOLD}Usage:${NC}
  $(basename "$0") -m <short|long|hybrid> [options]

${BOLD}Assembly mode:${NC}
  -m  Assembly mode: short | long | hybrid  (REQUIRED)

${BOLD}Read inputs:${NC}
  -1  Illumina forward reads  (required for short / hybrid)
  -2  Illumina reverse reads  (required for short / hybrid)
  -l  Long reads FASTQ        (required for long  / hybrid)

${BOLD}Required for all modes:${NC}
  -s  Sample name
  -o  Output directory

${BOLD}Optional:${NC}
  -t  Threads                      (default: ${THREADS})
  -g  Organism for AMRFinder+      (e.g. Escherichia, Klebsiella, Salmonella)
  -a  Trimmomatic adapters FASTA   (default: ${TRIMMOMATIC_ADAPTERS})
  -q  NanoFilt minimum Q score     (default: ${NANO_MIN_Q})
  -L  NanoFilt minimum read length (default: ${NANO_MIN_LEN})
  -x  Estimated genome size for Flye, e.g. 5m (REQUIRED for long/hybrid)
  -r  Flye read type               (default: ${FLYE_READ_TYPE})
        nano-raw | nano-hq | nano-corr | pacbio-raw | pacbio-corr | pacbio-hifi
  -h  Show this help

${BOLD}AMR Rules environment:${NC}
  AMR Rules is always executed inside conda env '${AMRRULES_ENV}'.
  Conda must be initialised in the shell PATH (i.e. 'conda' must be findable).

${BOLD}Examples:${NC}
  # Short-read assembly
  $(basename "$0") -m short -1 R1.fq.gz -2 R2.fq.gz -s SAMPLE -o results/

  # Long-read (Nanopore) assembly
  $(basename "$0") -m long -l ont.fq.gz -s SAMPLE -o results/ -x 5m

  # Hybrid assembly
  $(basename "$0") -m hybrid -1 R1.fq.gz -2 R2.fq.gz -l ont.fq.gz \\
                  -s SAMPLE -o results/ -x 5m
EOF
exit 0
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while getopts ":m:1:2:l:s:o:t:g:a:q:L:x:r:h" opt; do
  case $opt in
    m) MODE="$OPTARG" ;;
    1) R1="$OPTARG" ;;
    2) R2="$OPTARG" ;;
    l) LONG_READS="$OPTARG" ;;
    s) SAMPLE="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    g) ORGANISM="$OPTARG" ;;
    a) TRIMMOMATIC_ADAPTERS="$OPTARG" ;;
    q) NANO_MIN_Q="$OPTARG" ;;
    L) NANO_MIN_LEN="$OPTARG" ;;
    x) GENOME_SIZE="$OPTARG" ;;
    r) FLYE_READ_TYPE="$OPTARG" ;;
    h) usage ;;
    :) die "Option -$OPTARG requires an argument." ;;
    \?) die "Unknown option: -$OPTARG" ;;
  esac
done

# ---------------------------------------------------------------------------
# Validate mode and inputs
# ---------------------------------------------------------------------------
[[ -z "${MODE:-}"   ]] && die "Assembly mode (-m) is required. Choose: short | long | hybrid"
[[ -z "${SAMPLE:-}" ]] && die "Sample name (-s) is required."
[[ -z "${OUTDIR:-}" ]] && die "Output directory (-o) is required."

case "$MODE" in
  short)
    [[ -z "$R1" || -z "$R2" ]] && die "Short mode requires -1 and -2 (Illumina paired reads)."
    [[ -f "$R1" ]] || die "File not found: $R1"
    [[ -f "$R2" ]] || die "File not found: $R2"
    ;;
  long)
    [[ -z "$LONG_READS" ]] && die "Long mode requires -l (long reads FASTQ)."
    [[ -f "$LONG_READS" ]] || die "File not found: $LONG_READS"
    [[ -z "$GENOME_SIZE" ]] && die "Long mode requires -x <genome_size> for Flye (e.g. -x 5m)."
    ;;
  hybrid)
    [[ -z "$R1" || -z "$R2" ]] && die "Hybrid mode requires -1 and -2 (Illumina paired reads)."
    [[ -z "$LONG_READS" ]]     && die "Hybrid mode requires -l (long reads FASTQ)."
    [[ -f "$R1" ]]             || die "File not found: $R1"
    [[ -f "$R2" ]]             || die "File not found: $R2"
    [[ -f "$LONG_READS" ]]     || die "File not found: $LONG_READS"
    [[ -z "$GENOME_SIZE" ]]    && die "Hybrid mode requires -x <genome_size> for Flye (e.g. -x 5m)."
    ;;
  *) die "Unknown mode '${MODE}'. Choose: short | long | hybrid" ;;
esac

# ---------------------------------------------------------------------------
# Dependency checks
# ---------------------------------------------------------------------------
check_tool() {
  command -v "$1" &>/dev/null || die "Required tool not found in PATH: $1"
}

log "Checking dependencies for mode: ${MODE}..."

# Tools required in all modes
for tool in amrfinder; do
  check_tool "$tool"
done

case "$MODE" in
  short)
    for tool in fastqc trimmomatic unicycler; do check_tool "$tool"; done ;;
  long)
    for tool in nanostat nanofilt flye medaka_consensus; do check_tool "$tool"; done ;;
  hybrid)
    for tool in fastqc trimmomatic nanostat nanofilt unicycler; do check_tool "$tool"; done ;;
esac

# Verify conda is accessible and the AMR Rules environment exists
check_tool conda
conda env list | grep -qw "${AMRRULES_ENV}" \
  || die "Conda environment '${AMRRULES_ENV}' not found. Create it or adjust -e flag."

ok "All required tools found."

# ---------------------------------------------------------------------------
# Directory structure
# ---------------------------------------------------------------------------
QC_DIR="${OUTDIR}/01_qc"
TRIM_DIR="${OUTDIR}/02_trimmed"
ASSEMBLY_DIR="${OUTDIR}/03_assembly"
AMRFINDER_DIR="${OUTDIR}/04_amrfinder"
AMRRULES_DIR="${OUTDIR}/05_amr_rules"
LOG_DIR="${OUTDIR}/logs"

mkdir -p "$QC_DIR" "$TRIM_DIR" "$ASSEMBLY_DIR" \
         "$AMRFINDER_DIR" "$AMRRULES_DIR" "$LOG_DIR"

log "Output directory : ${OUTDIR}"
log "Sample           : ${SAMPLE}"
log "Assembly mode    : ${MODE}"
log "Threads          : ${THREADS}"
[[ -n "$ORGANISM" ]] && log "Organism         : ${ORGANISM}"

# =============================================================================
# STAGE 1 — Quality Control
# =============================================================================
echo ""
log "${BOLD}========== STAGE 1: QUALITY CONTROL ==========${NC}"

# ── SHORT / HYBRID — Illumina QC ────────────────────────────────────────────
if [[ "$MODE" == "short" || "$MODE" == "hybrid" ]]; then

  log "[Illumina] Running FastQC on raw reads..."
  fastqc \
    --outdir "${QC_DIR}" \
    --threads "$THREADS" \
    "$R1" "$R2" \
    &> "${LOG_DIR}/${SAMPLE}_fastqc_raw.log"
  ok "FastQC (raw) → ${QC_DIR}"

  log "[Illumina] Running Trimmomatic..."
  R1_PAIRED="${TRIM_DIR}/${SAMPLE}_R1_paired.fastq.gz"
  R2_PAIRED="${TRIM_DIR}/${SAMPLE}_R2_paired.fastq.gz"
  R1_UNPAIRED="${TRIM_DIR}/${SAMPLE}_R1_unpaired.fastq.gz"
  R2_UNPAIRED="${TRIM_DIR}/${SAMPLE}_R2_unpaired.fastq.gz"

  if command -v trimmomatic &>/dev/null; then
    TRIM_CMD="trimmomatic"
  else
    TRIM_JAR="$(find /usr -name 'trimmomatic*.jar' 2>/dev/null | head -1)"
    [[ -z "$TRIM_JAR" ]] && die "Trimmomatic wrapper and JAR not found."
    TRIM_CMD="java -jar $TRIM_JAR"
  fi

  $TRIM_CMD PE \
    -threads "$THREADS" \
    "$R1" "$R2" \
    "$R1_PAIRED" "$R1_UNPAIRED" \
    "$R2_PAIRED" "$R2_UNPAIRED" \
    ILLUMINACLIP:"${TRIMMOMATIC_ADAPTERS}":2:30:10 \
    LEADING:"${LEADING}" \
    TRAILING:"${TRAILING}" \
    SLIDINGWINDOW:"${SLIDINGWINDOW}" \
    MINLEN:"${MIN_LEN}" \
    &> "${LOG_DIR}/${SAMPLE}_trimmomatic.log"
  ok "Trimmomatic → ${TRIM_DIR}"

  log "[Illumina] Running FastQC on trimmed reads..."
  fastqc \
    --outdir "${QC_DIR}" \
    --threads "$THREADS" \
    "$R1_PAIRED" "$R2_PAIRED" \
    &> "${LOG_DIR}/${SAMPLE}_fastqc_trimmed.log"
  ok "FastQC (trimmed) → ${QC_DIR}"
fi

# ── LONG / HYBRID — Nanopore / long-read QC ─────────────────────────────────
if [[ "$MODE" == "long" || "$MODE" == "hybrid" ]]; then

  log "[Long reads] Running NanoStat on raw reads..."
  NanoStat \
    --fastq "$LONG_READS" \
    --threads "$THREADS" \
    --outdir "${QC_DIR}" \
    --name "${SAMPLE}_nanostat_raw.txt" \
    &> "${LOG_DIR}/${SAMPLE}_nanostat.log"
  ok "NanoStat → ${QC_DIR}/${SAMPLE}_nanostat_raw.txt"

  log "[Long reads] Filtering with NanoFilt (min Q=${NANO_MIN_Q}, min len=${NANO_MIN_LEN})..."
  LONG_FILTERED="${TRIM_DIR}/${SAMPLE}_long_filtered.fastq.gz"
  gunzip -c "$LONG_READS" \
    | NanoFilt -q "$NANO_MIN_Q" -l "$NANO_MIN_LEN" \
    | gzip > "$LONG_FILTERED" \
    2> "${LOG_DIR}/${SAMPLE}_nanofilt.log"
  ok "NanoFilt → ${LONG_FILTERED}"
fi

# =============================================================================
# STAGE 2 — Assembly
# =============================================================================
echo ""
log "${BOLD}========== STAGE 2: ASSEMBLY (mode: ${MODE}) ==========${NC}"

ASSEMBLY_FASTA=""

# ── SHORT — Unicycler (Illumina only) ───────────────────────────────────────
if [[ "$MODE" == "short" ]]; then

  log "[Short] Running Unicycler (Illumina PE)..."
  unicycler \
    -1 "$R1_PAIRED" \
    -2 "$R2_PAIRED" \
    --out "${ASSEMBLY_DIR}/unicycler" \
    --threads "$THREADS" \
    &> "${LOG_DIR}/${SAMPLE}_unicycler.log"

  RAW_FASTA="${ASSEMBLY_DIR}/unicycler/assembly.fasta"
  [[ -f "$RAW_FASTA" ]] \
    || die "Unicycler assembly.fasta not found. Check: ${LOG_DIR}/${SAMPLE}_unicycler.log"

  ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${SAMPLE}_assembly.fasta"
  cp "$RAW_FASTA" "$ASSEMBLY_FASTA"
  ok "Unicycler assembly → ${ASSEMBLY_FASTA}"
fi

# ── LONG — Flye + Medaka ────────────────────────────────────────────────────
if [[ "$MODE" == "long" ]]; then

  FLYE_OUTDIR="${ASSEMBLY_DIR}/flye"
  log "[Long] Running Flye (read type: ${FLYE_READ_TYPE}, genome size: ${GENOME_SIZE})..."
  flye \
    "--${FLYE_READ_TYPE}" "$LONG_FILTERED" \
    --genome-size "$GENOME_SIZE" \
    --out-dir "$FLYE_OUTDIR" \
    --threads "$THREADS" \
    &> "${LOG_DIR}/${SAMPLE}_flye.log"

  FLYE_FASTA="${FLYE_OUTDIR}/assembly.fasta"
  [[ -f "$FLYE_FASTA" ]] \
    || die "Flye assembly.fasta not found. Check: ${LOG_DIR}/${SAMPLE}_flye.log"
  ok "Flye assembly → ${FLYE_FASTA}"

  # Medaka consensus polishing
  log "[Long] Running Medaka polishing..."
  MEDAKA_OUTDIR="${ASSEMBLY_DIR}/medaka"

  # Resolve medaka model: use environment variable if set, else default
  MEDAKA_MODEL="${MEDAKA_MODEL:-r941_min_sup_g507}"
  log "[Long] Medaka model: ${MEDAKA_MODEL}"

  medaka_consensus \
    -i "$LONG_FILTERED" \
    -d "$FLYE_FASTA" \
    -o "$MEDAKA_OUTDIR" \
    -t "$THREADS" \
    -m "$MEDAKA_MODEL" \
    &> "${LOG_DIR}/${SAMPLE}_medaka.log"

  MEDAKA_FASTA="${MEDAKA_OUTDIR}/consensus.fasta"
  [[ -f "$MEDAKA_FASTA" ]] \
    || die "Medaka consensus.fasta not found. Check: ${LOG_DIR}/${SAMPLE}_medaka.log"

  ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${SAMPLE}_assembly.fasta"
  cp "$MEDAKA_FASTA" "$ASSEMBLY_FASTA"
  ok "Medaka polished assembly → ${ASSEMBLY_FASTA}"
fi

# ── HYBRID — Unicycler (short + long) ───────────────────────────────────────
if [[ "$MODE" == "hybrid" ]]; then

  log "[Hybrid] Running Unicycler (Illumina PE + long reads)..."
  unicycler \
    -1 "$R1_PAIRED" \
    -2 "$R2_PAIRED" \
    -l "$LONG_FILTERED" \
    --out "${ASSEMBLY_DIR}/unicycler_hybrid" \
    --threads "$THREADS" \
    &> "${LOG_DIR}/${SAMPLE}_unicycler_hybrid.log"

  RAW_FASTA="${ASSEMBLY_DIR}/unicycler_hybrid/assembly.fasta"
  [[ -f "$RAW_FASTA" ]] \
    || die "Unicycler hybrid assembly.fasta not found. Check: ${LOG_DIR}/${SAMPLE}_unicycler_hybrid.log"

  ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${SAMPLE}_assembly.fasta"
  cp "$RAW_FASTA" "$ASSEMBLY_FASTA"
  ok "Unicycler hybrid assembly → ${ASSEMBLY_FASTA}"
fi

# Sanity check and basic stats
[[ -f "$ASSEMBLY_FASTA" ]] || die "Assembly FASTA missing after Stage 2 — this should not happen."
CONTIGS=$(grep -c '^>' "$ASSEMBLY_FASTA" || true)
log "Assembly contains ${CONTIGS} contig(s)."

# =============================================================================
# STAGE 3 — AMRFinder+
# =============================================================================
echo ""
log "${BOLD}========== STAGE 3: AMRFinder+ ==========${NC}"

log "Updating AMRFinder+ database..."
amrfinder --update \
  &> "${LOG_DIR}/${SAMPLE}_amrfinder_update.log" \
  || warn "AMRFinder+ database update failed — proceeding with existing database."

AMRFINDER_TSV="${AMRFINDER_DIR}/${SAMPLE}_amrfinder.tsv"

ORG_FLAG=""
[[ -n "$ORGANISM" ]] && ORG_FLAG="--organism ${ORGANISM}"

log "Running AMRFinder+..."
# shellcheck disable=SC2086
amrfinder \
  --nucleotide "$ASSEMBLY_FASTA" \
  --output "$AMRFINDER_TSV" \
  --threads "$THREADS" \
  --plus \
  $ORG_FLAG \
  &> "${LOG_DIR}/${SAMPLE}_amrfinder.log"

[[ -f "$AMRFINDER_TSV" ]] \
  || die "AMRFinder+ output not found. Check: ${LOG_DIR}/${SAMPLE}_amrfinder.log"

AMR_HITS=$(tail -n +2 "$AMRFINDER_TSV" | wc -l)
ok "AMRFinder+ complete — ${AMR_HITS} AMR element(s) detected → ${AMRFINDER_TSV}"

# =============================================================================
# STAGE 4 — AMR Rules  (runs inside conda env 'amrrules_beta')
# =============================================================================
# Strategy: use 'conda run -n <env>' which:
#   • activates the named environment for a single subprocess
#   • works in non-interactive shells (set -euo pipefail safe)
#   • does NOT require 'conda activate' in the calling shell
#   • does NOT require 'source conda.sh' / 'conda init'
# =============================================================================
echo ""
log "${BOLD}========== STAGE 4: AMR Rules (env: ${AMRRULES_ENV}) ==========${NC}"

AMRRULES_OUT="${AMRRULES_DIR}/${SAMPLE}_amr_rules"

# Verify amr_rules is available inside the target environment before running
log "Verifying amr_rules availability in '${AMRRULES_ENV}'..."
conda run -n "${AMRRULES_ENV}" python3 -c "import amr_rules" 2>/dev/null \
  || conda run -n "${AMRRULES_ENV}" bash -c "command -v amr_rules" 2>/dev/null \
  || die "amr_rules not found inside conda env '${AMRRULES_ENV}'. " \
         "Activate the environment manually and run: pip install amr-rules"

log "Running AMR Rules on AMRFinder+ output..."
# Prefer the CLI entry point; fall back to python -m if the binary is absent
if conda run -n "${AMRRULES_ENV}" bash -c "command -v amr_rules" &>/dev/null; then
  conda run -n "${AMRRULES_ENV}" \
    amr_rules \
      --input  "$AMRFINDER_TSV" \
      --output "$AMRRULES_OUT" \
    &> "${LOG_DIR}/${SAMPLE}_amr_rules.log"
else
  conda run -n "${AMRRULES_ENV}" \
    python3 -m amr_rules \
      --input  "$AMRFINDER_TSV" \
      --output "$AMRRULES_OUT" \
    &> "${LOG_DIR}/${SAMPLE}_amr_rules.log"
fi

ok "AMR Rules complete → ${AMRRULES_DIR}"

# =============================================================================
# SUMMARY REPORT
# =============================================================================
echo ""
SUMMARY="${OUTDIR}/${SAMPLE}_pipeline_summary.txt"
{
  echo "======================================================"
  echo "  Prokaryotic AMR Pipeline — Summary"
  echo "  Sample : ${SAMPLE}"
  echo "  Mode   : ${MODE}"
  echo "  Date   : $(date)"
  echo "======================================================"
  echo ""
  echo "Input reads"
  case "$MODE" in
    short)
      echo "  R1 : $R1"
      echo "  R2 : $R2"
      ;;
    long)
      echo "  Long reads : $LONG_READS"
      ;;
    hybrid)
      echo "  R1         : $R1"
      echo "  R2         : $R2"
      echo "  Long reads : $LONG_READS"
      ;;
  esac
  echo ""
  echo "Assembly"
  echo "  Tool    : $(case "$MODE" in short) echo "Unicycler";; long) echo "Flye + Medaka";; hybrid) echo "Unicycler (hybrid)";; esac)"
  echo "  FASTA   : $ASSEMBLY_FASTA"
  echo "  Contigs : $CONTIGS"
  echo ""
  echo "AMRFinder+"
  echo "  Report   : $AMRFINDER_TSV"
  echo "  AMR hits : $AMR_HITS"
  echo ""
  echo "AMR Rules"
  echo "  Conda env : ${AMRRULES_ENV}"
  echo "  Output    : ${AMRRULES_DIR}"
  echo ""
  echo "Logs : $LOG_DIR"
  echo "======================================================"
} | tee "$SUMMARY"

echo ""
ok "${BOLD}Pipeline completed successfully — sample: ${SAMPLE}${NC}"
ok "Summary written to: ${SUMMARY}"
