#!/usr/bin/env bash
# =============================================================================
# Prokaryotic Genome Assembly & AMR Analysis Pipeline
# =============================================================================
# Stages:
#   1. Quality control (FastQC + Trimmomatic) + Assembly (Unicycler)
#   2. AMR gene detection (AMRFinder+)
#   3. Phenotypic AMR interpretation (AMR Rules)
#
# Dependencies:
#   fastqc, trimmomatic, unicycler, amrfinder, python3 (for amr_rules)
#
# Usage:
#   bash prokaryotic_amr_pipeline.sh \
#     -1 reads_R1.fastq.gz \
#     -2 reads_R2.fastq.gz \
#     -s <sample_name> \
#     -o <output_dir> \
#     [-t <threads>] \
#     [-g <organism>]          # e.g. Escherichia, Salmonella (optional, for AMRFinder+)
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Default parameters
# ---------------------------------------------------------------------------
THREADS=8
ORGANISM=""                     # Leave blank if organism flag is not needed
TRIMMOMATIC_ADAPTERS="NexteraPE-PE.fa"   # Adjust to your adapter file
MIN_LEN=50                      # Minimum read length after trimming
LEADING=3
TRAILING=3
SLIDINGWINDOW="4:15"

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
  $(basename "$0") -1 R1.fastq.gz -2 R2.fastq.gz -s SAMPLE -o OUTDIR [options]

${BOLD}Required:${NC}
  -1  Forward reads (FASTQ / FASTQ.gz)
  -2  Reverse reads (FASTQ / FASTQ.gz)
  -s  Sample name (used for file naming)
  -o  Output directory

${BOLD}Optional:${NC}
  -t  Threads (default: ${THREADS})
  -g  Organism for AMRFinder+ (e.g. Escherichia, Salmonella, Klebsiella)
  -a  Trimmomatic adapters FASTA (default: ${TRIMMOMATIC_ADAPTERS})
  -h  Show this help
EOF
exit 0
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
while getopts ":1:2:s:o:t:g:a:h" opt; do
  case $opt in
    1) R1="$OPTARG" ;;
    2) R2="$OPTARG" ;;
    s) SAMPLE="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    g) ORGANISM="$OPTARG" ;;
    a) TRIMMOMATIC_ADAPTERS="$OPTARG" ;;
    h) usage ;;
    :) die "Option -$OPTARG requires an argument." ;;
    \?) die "Unknown option: -$OPTARG" ;;
  esac
done

# ---------------------------------------------------------------------------
# Validate required inputs
# ---------------------------------------------------------------------------
[[ -z "${R1:-}" ]]     && die "Forward reads (-1) required."
[[ -z "${R2:-}" ]]     && die "Reverse reads (-2) required."
[[ -z "${SAMPLE:-}" ]] && die "Sample name (-s) required."
[[ -z "${OUTDIR:-}" ]] && die "Output directory (-o) required."
[[ -f "$R1" ]]         || die "File not found: $R1"
[[ -f "$R2" ]]         || die "File not found: $R2"

# ---------------------------------------------------------------------------
# Check dependencies
# ---------------------------------------------------------------------------
check_tool() {
  command -v "$1" &>/dev/null || die "Required tool not found in PATH: $1"
}
for tool in fastqc trimmomatic unicycler amrfinder python3; do
  check_tool "$tool"
done
ok "All required tools found."

# ---------------------------------------------------------------------------
# Directory structure
# ---------------------------------------------------------------------------
QC_RAW_DIR="${OUTDIR}/01_fastqc_raw"
TRIM_DIR="${OUTDIR}/02_trimmomatic"
QC_TRIM_DIR="${OUTDIR}/03_fastqc_trimmed"
ASSEMBLY_DIR="${OUTDIR}/04_assembly"
AMRFINDER_DIR="${OUTDIR}/05_amrfinder"
AMRRULES_DIR="${OUTDIR}/06_amr_rules"
LOG_DIR="${OUTDIR}/logs"

mkdir -p "$QC_RAW_DIR" "$TRIM_DIR" "$QC_TRIM_DIR" \
         "$ASSEMBLY_DIR" "$AMRFINDER_DIR" "$AMRRULES_DIR" "$LOG_DIR"

log "Output directory: ${OUTDIR}"
log "Sample:           ${SAMPLE}"
log "Threads:          ${THREADS}"
[[ -n "$ORGANISM" ]] && log "Organism:         ${ORGANISM}"

# =============================================================================
# STAGE 1 — Quality Control & Assembly
# =============================================================================
echo ""
log "${BOLD}========== STAGE 1: QC & ASSEMBLY ==========${NC}"

# --- 1a. FastQC on raw reads -------------------------------------------------
log "Running FastQC on raw reads..."
fastqc \
  --outdir "$QC_RAW_DIR" \
  --threads "$THREADS" \
  "$R1" "$R2" \
  &> "${LOG_DIR}/${SAMPLE}_fastqc_raw.log"
ok "FastQC (raw) complete → ${QC_RAW_DIR}"

# --- 1b. Trimmomatic ---------------------------------------------------------
log "Running Trimmomatic..."

R1_PAIRED="${TRIM_DIR}/${SAMPLE}_R1_paired.fastq.gz"
R2_PAIRED="${TRIM_DIR}/${SAMPLE}_R2_paired.fastq.gz"
R1_UNPAIRED="${TRIM_DIR}/${SAMPLE}_R1_unpaired.fastq.gz"
R2_UNPAIRED="${TRIM_DIR}/${SAMPLE}_R2_unpaired.fastq.gz"

# Trimmomatic may be invoked as a JAR or as a wrapper; try wrapper first
if command -v trimmomatic &>/dev/null; then
  TRIM_CMD="trimmomatic"
else
  TRIM_JAR="$(find /usr -name 'trimmomatic*.jar' 2>/dev/null | head -1)"
  [[ -z "$TRIM_JAR" ]] && die "Trimmomatic JAR not found."
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
ok "Trimmomatic complete → ${TRIM_DIR}"

# --- 1c. FastQC on trimmed reads ---------------------------------------------
log "Running FastQC on trimmed reads..."
fastqc \
  --outdir "$QC_TRIM_DIR" \
  --threads "$THREADS" \
  "$R1_PAIRED" "$R2_PAIRED" \
  &> "${LOG_DIR}/${SAMPLE}_fastqc_trimmed.log"
ok "FastQC (trimmed) complete → ${QC_TRIM_DIR}"

# --- 1d. Unicycler assembly --------------------------------------------------
log "Running Unicycler assembly..."
unicycler \
  -1 "$R1_PAIRED" \
  -2 "$R2_PAIRED" \
  --out "${ASSEMBLY_DIR}" \
  --threads "$THREADS" \
  &> "${LOG_DIR}/${SAMPLE}_unicycler.log"

ASSEMBLY_FASTA="${ASSEMBLY_DIR}/assembly.fasta"
[[ -f "$ASSEMBLY_FASTA" ]] || die "Assembly FASTA not found. Check Unicycler log: ${LOG_DIR}/${SAMPLE}_unicycler.log"

# Rename for clarity
cp "$ASSEMBLY_FASTA" "${ASSEMBLY_DIR}/${SAMPLE}_assembly.fasta"
ASSEMBLY_FASTA="${ASSEMBLY_DIR}/${SAMPLE}_assembly.fasta"

ok "Unicycler assembly complete → ${ASSEMBLY_FASTA}"

# Basic assembly stats
CONTIGS=$(grep -c '^>' "$ASSEMBLY_FASTA" || true)
log "Assembly contains ${CONTIGS} contig(s)."

# =============================================================================
# STAGE 2 — AMRFinder+
# =============================================================================
echo ""
log "${BOLD}========== STAGE 2: AMRFinder+ ==========${NC}"

# Update AMRFinder+ database (skip with --no-update if offline)
log "Updating AMRFinder+ database..."
amrfinder --update &> "${LOG_DIR}/${SAMPLE}_amrfinder_update.log" \
  || warn "AMRFinder+ database update failed — proceeding with existing database."

AMRFINDER_TSV="${AMRFINDER_DIR}/${SAMPLE}_amrfinder.tsv"

# Build organism flag if provided
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

[[ -f "$AMRFINDER_TSV" ]] || die "AMRFinder+ output not found. Check: ${LOG_DIR}/${SAMPLE}_amrfinder.log"

AMR_HITS=$(tail -n +2 "$AMRFINDER_TSV" | wc -l)
ok "AMRFinder+ complete — ${AMR_HITS} AMR element(s) detected → ${AMRFINDER_TSV}"

# =============================================================================
# STAGE 3 — AMR Rules
# =============================================================================
echo ""
log "${BOLD}========== STAGE 3: AMR Rules ==========${NC}"

AMRRULES_OUT="${AMRRULES_DIR}/${SAMPLE}_amr_rules"

# Locate amr_rules entry point (adjust if installed differently)
if command -v amr_rules &>/dev/null; then
  AMR_RULES_CMD="amr_rules"
elif python3 -c "import amr_rules" &>/dev/null 2>&1; then
  AMR_RULES_CMD="python3 -m amr_rules"
else
  die "amr_rules not found. Install via: pip install amr-rules"
fi

log "Running AMR Rules on AMRFinder+ output..."

$AMR_RULES_CMD \
  --input "$AMRFINDER_TSV" \
  --output "$AMRRULES_OUT" \
  &> "${LOG_DIR}/${SAMPLE}_amr_rules.log"

ok "AMR Rules complete → ${AMRRULES_DIR}"

# =============================================================================
# SUMMARY REPORT
# =============================================================================
echo ""
SUMMARY="${OUTDIR}/${SAMPLE}_pipeline_summary.txt"
{
  echo "======================================================"
  echo "  Pipeline Summary — ${SAMPLE}"
  echo "  $(date)"
  echo "======================================================"
  echo ""
  echo "Input reads"
  echo "  R1:  $R1"
  echo "  R2:  $R2"
  echo ""
  echo "Assembly"
  echo "  FASTA:    $ASSEMBLY_FASTA"
  echo "  Contigs:  $CONTIGS"
  echo ""
  echo "AMRFinder+"
  echo "  Report:   $AMRFINDER_TSV"
  echo "  AMR hits: $AMR_HITS"
  echo ""
  echo "AMR Rules"
  echo "  Output directory: $AMRRULES_DIR"
  echo ""
  echo "Logs: $LOG_DIR"
  echo "======================================================"
} | tee "$SUMMARY"

echo ""
ok "${BOLD}Pipeline completed successfully for sample: ${SAMPLE}${NC}"
ok "Summary written to: ${SUMMARY}"
