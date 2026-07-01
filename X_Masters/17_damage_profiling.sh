#!/bin/bash
#SBATCH --job-name=17dmgprof
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

### 17. Damage Profiling — mapDamage2 per (sample × gut taxon) ====================================================
# Runs mapDamage2 on every dedup BAM produced by 15_genome_mapping.sh to quantify
# aDNA damage patterns (5' C→T, 3' G→A, read-length distribution) per reference genome.
#
# USAGE:  sbatch --array 4-8 scripts/17_damage_profiling.sh <accesslist> <BAMDIR> <REFDIR> <OUTDIR>
# EXAMPLE:
#   sbatch --array 4-8 scripts/17_damage_profiling.sh \
#       00_accessions/picareiro \
#       15_genome_mapping/picareiro/bam \
#       resources/gut_taxa_genomes \
#       17_damage_profiling/picareiro
#
# INPUTS:
#   - Dedup BAMs:   ${BAMDIR}/${PREFIX}_${GENOME}.dedup.bam   (from 15_genome_mapping)
#   - References:   ${REFDIR}/${GENOME}.fa.gz                  (24 gut taxa)
#
# OUTPUTS per sample:
#   - ${OUTDIR}/${PREFIX}/${GENOME}/                          mapDamage2 full output
#       ├── Fragmisincorporation_plot.pdf                    5'/3' misincorporation plot
#       ├── Length_plot.pdf                                  read-length distribution
#       ├── 5pCtoT_freq.txt                                  5' C→T table (per position)
#       ├── 3pGtoA_freq.txt                                  3' G→A table
#       ├── misincorporation.txt                             full misincorporation matrix
#       ├── lgdistribution.txt                               read length distribution
#       └── Runtime_log.txt                                  run log
#   - ${OUTDIR}/${PREFIX}.damage_summary.tsv                  cross-genome summary
#
# WHY mapDamage2:
#   - DamageProfiler is not available on LISC; mapDamage2 is the original tool
#     it was designed to replace and remains a community standard
#   - Operates per BAM + per reference, matching our BBSplit-split outputs
#
# UDG-HALF LIBRARY PREP — IMPORTANT INTERPRETATION NOTE:
#   ** Picareiro libraries were prepared with UDG-HALF treatment **
#   Partial uracil-DNA-glycosylase digestion removes ~half of deaminated cytosines
#   from internal positions while leaving terminal 1–2 bp largely intact. Expected
#   damage is therefore HALF the magnitude of untreated aDNA, concentrated at the
#   very ends. Adjusted authentication thresholds:
#     - 5' C→T position 1   ≥ 2–3%   (vs ≥5% for non-UDG; ≥5% strong, ≥10% pristine)
#     - 5' C→T position 2   sharp drop-off vs pos 1 (UDG signature)
#     - Mean read length     30–70 bp  (>100 bp suggests modern contamination)
#     - Minimum reads        ≥ 1000    per genome for reliable damage estimation
#
#   mapDamage2 has NO UDG-aware option. Its only related flag is `--rescale`,
#   which rewrites BAM quality scores at damaged positions for downstream variant
#   calling — it is a corrective step, not an estimation mode, and we deliberately
#   do NOT use it here. The reported C→T / G→A frequencies are RAW observations;
#   the UDG-half adjustment is applied at INTERPRETATION time (thresholds above),
#   not in this script.
#
# SKIPPING LOW-COVERAGE BAMs:
#   - MIN_READS env var (default 100) — skip BAMs with fewer reads (mapDamage returns noise)
#   - Override:  MIN_READS=500 sbatch ...
# ===================================================================================================================

set -euo pipefail

### Variables
HOMEDIR=$(pwd)
ACCESS="$1"
DATASET=$(basename "$1")
BAMDIR="${2%/}"
REFDIR="${3%/}"
OUTDIR="${4%/}"
PREFIX=$(awk '{if (NR=='${SLURM_ARRAY_TASK_ID}') print $1}' "$ACCESS")
MIN_READS="${MIN_READS:-100}"

### Load modules
# mapDamage2 (foss-2022a) is used in place of DamageProfiler — the latter is
# not available on LISC. mapDamage2 is functionally equivalent for our purposes:
# 5' C→T / 3' G→A misincorporation tables + read-length distribution.
# It is the original tool DamageProfiler was built to replace and is still the
# default fallback in nf-core/eager when DamageProfiler is unavailable.
#
# The mapDamage module does NOT expose `samtools` on $PATH, so we also load
# SAMtools explicitly for `samtools view -c`, `samtools faidx`, and the
# combined-reference build below. We pin SAMtools/1.16.1-GCC-11.3.0 because it
# shares the GCC-11.3.0 toolchain underlying foss-2022a (mapDamage's toolchain),
# avoiding the Lmod toolchain warning that newer SAMtools builds trigger.
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
set +u
module purge
module load mapDamage/2.2.3-foss-2022a
module load SAMtools/1.16.1-GCC-11.3.0
set -u

### Create output directories
mkdir -p "$HOMEDIR/$OUTDIR/$PREFIX"
mkdir -p "$TMPDIR/$OUTDIR/$PREFIX"
mkdir -p "$TMPDIR/bams"

### Cache decompressed + indexed references to node-local storage
### (shared across array tasks of this job array via /tmp + flock; same pattern as script 15)
CACHE_DIR="/tmp/slurm_cache_${SLURM_ARRAY_JOB_ID}_dmgprof"
CACHE_LOCK="/tmp/slurm_cache_${SLURM_ARRAY_JOB_ID}_dmgprof.lock"
(
  flock -x -w 3600 200
  if [[ ! -f "$CACHE_DIR/.cached" ]]; then
    mkdir -p "$CACHE_DIR"
    echo "Decompressing + indexing references into $CACHE_DIR ..."
    for refgz in "$HOMEDIR/$REFDIR"/*.fa.gz; do
      base=$(basename "$refgz" .fa.gz)
      zcat "$refgz" > "$CACHE_DIR/${base}.fa"
      samtools faidx "$CACHE_DIR/${base}.fa"
    done
    # Combined reference: BBSplit in 15_genome_mapping.sh is called with
    # ref=<all 24 .fa.gz concatenated>, so every per-genome dedup BAM inherits
    # @SQ headers for the UNION of all 24 genomes' contigs. mapDamage refuses
    # to run if any BAM @SQ contig is missing from the reference dictionary,
    # so we build one combined FASTA that mirrors what BBSplit indexed. Each
    # BAM still only contains reads for its own genome's contigs, so the
    # damage tables come out per-genome correctly.
    echo "Building combined reference FASTA ..."
    cat "$CACHE_DIR"/*.fa > "$CACHE_DIR/combined.fa"
    samtools faidx "$CACHE_DIR/combined.fa"
    touch "$CACHE_DIR/.cached"
    echo "  Cached $(ls $CACHE_DIR/*.fa | grep -v combined.fa | wc -l) references (+ combined.fa)"
  fi
) 200>"$CACHE_LOCK"
REF_CACHE="$CACHE_DIR"
COMBINED_FA="$REF_CACHE/combined.fa"

echo "============================================================="
echo "Sample:    $PREFIX"
echo "Dataset:   $DATASET"
echo "BAM dir:   $BAMDIR"
echo "Ref cache: $REF_CACHE ($(ls "$REF_CACHE"/*.fa 2>/dev/null | wc -l) refs)"
echo "Output:    $OUTDIR/$PREFIX"
echo "Min reads: $MIN_READS"
echo "============================================================="

### Discover BAMs for this sample
shopt -s nullglob
BAMS=( "$HOMEDIR/$BAMDIR/${PREFIX}_"*.dedup.bam )
shopt -u nullglob

if [[ ${#BAMS[@]} -eq 0 ]]; then
    echo "ERROR: no dedup BAMs found matching $BAMDIR/${PREFIX}_*.dedup.bam"
    exit 1
fi
echo "Found ${#BAMS[@]} BAM(s) for $PREFIX"

### Stage all BAMs (+ their .bai) for this sample into TMPDIR for fast local I/O
for B in "${BAMS[@]}"; do
    rsync -a "$B" "$TMPDIR/bams/"
    [[ -f "${B}.bai" ]] && rsync -a "${B}.bai" "$TMPDIR/bams/"
done

### Initialise summary table
# CT_drop = CT_5p_pos1 / CT_5p_pos2 ; GA_drop = GA_3p_pos1 / GA_3p_pos2.
# Authentic aDNA (esp. UDG-half) shows a sharp pos1->pos2 cliff (drop >= 5x);
# modern DNA shows ratios near 1 (no cliff). These ratios make the per-genome
# authentication call directly readable from the summary without re-deriving.
SUMMARY="$HOMEDIR/$OUTDIR/${PREFIX}.damage_summary.tsv"
echo -e "Sample\tGenome\tNreads\tMean_ReadLen\tCT_5p_pos1\tCT_5p_pos2\tCT_drop\tGA_3p_pos1\tGA_3p_pos2\tGA_drop" > "$SUMMARY"

### Process each BAM (reading from TMPDIR-local copies)
for BAM_SRC in "${BAMS[@]}"; do
    BAM_BASE=$(basename "$BAM_SRC" .dedup.bam)
    GENOME="${BAM_BASE#${PREFIX}_}"
    BAM="$TMPDIR/bams/$(basename "$BAM_SRC")"
    # Use the combined FASTA for every BAM — see cache block above for why.
    REF_FA="$COMBINED_FA"

    echo ""
    echo "--- $GENOME ---"

    if [[ ! -f "$REF_FA" ]]; then
        echo "  WARN: combined reference $REF_FA not found in cache, skipping"
        continue
    fi

    # Pre-flight: count reads, skip if below threshold
    NREADS=$(samtools view -c "$BAM")
    if [[ "$NREADS" -lt "$MIN_READS" ]]; then
        echo "  SKIP: $NREADS reads < MIN_READS=$MIN_READS"
        continue
    fi
    echo "  Reads: $NREADS"

    GENOME_OUT="$TMPDIR/$OUTDIR/$PREFIX/$GENOME"
    mkdir -p "$GENOME_OUT"

    # mapDamage2 writes all output files directly into -d (no nested dir),
    # so no flatten step needed (unlike DamageProfiler).
    if ! mapDamage \
            -i "$BAM" \
            -r "$REF_FA" \
            -d "$GENOME_OUT" \
            --no-stats; then
        echo "  ERROR: mapDamage failed on $GENOME"
        # Record FAIL row so absent genomes in the summary mean "not attempted",
        # while explicit failures stay visible for debugging.
        echo -e "${PREFIX}\t${GENOME}\t${NREADS}\tFAIL\tFAIL\tFAIL\tFAIL\tFAIL\tFAIL\tFAIL" >> "$SUMMARY"
        continue
    fi
    # --no-stats: skip Bayesian posterior estimation (the slow R/JAGS step).
    # We only need the descriptive 5'/3' frequency tables for authentication.

    # Extract key metrics for summary table from mapDamage2 output files.
    # File formats (mapDamage 2.2):
    #   5pCtoT_freq.txt:  pos<TAB>5pC>T   (header line, then 1-indexed positions)
    #   3pGtoA_freq.txt:  pos<TAB>3pG>A
    #   lgdistribution.txt: Std<TAB>Length<TAB>Occurences  (Std = '+' or '-')
    FIVE="$GENOME_OUT/5pCtoT_freq.txt"
    THREE="$GENOME_OUT/3pGtoA_freq.txt"
    LEN="$GENOME_OUT/lgdistribution.txt"

    CT1="NA"; CT2="NA"; GA1="NA"; GA2="NA"; MEANLEN="NA"

    if [[ -f "$FIVE" ]]; then
        CT1=$(awk -F'\t' 'NR>1 && $1==1 {printf "%.4f",$2}' "$FIVE")
        CT2=$(awk -F'\t' 'NR>1 && $1==2 {printf "%.4f",$2}' "$FIVE")
    fi
    if [[ -f "$THREE" ]]; then
        GA1=$(awk -F'\t' 'NR>1 && $1==1 {printf "%.4f",$2}' "$THREE")
        GA2=$(awk -F'\t' 'NR>1 && $1==2 {printf "%.4f",$2}' "$THREE")
    fi
    if [[ -f "$LEN" ]]; then
        # lgdistribution.txt: skip header + comment lines; sum length*occurrences
        # over both strands (Std='+' and Std='-').
        MEANLEN=$(awk '/^[+-]\t/ {tot+=$2*$3; cnt+=$3} END {if(cnt>0) printf "%.2f", tot/cnt; else print "NA"}' "$LEN")
    fi

    echo "  5' C→T pos1: $CT1 / pos2: $CT2"
    echo "  3' G→A pos1: $GA1 / pos2: $GA2"
    echo "  Mean read length: $MEANLEN"

    # Compute pos1/pos2 drop ratios (sharp cliff = UDG-half / aDNA signature).
    # Guard against pos2=0 (-> "Inf") and missing values (-> "NA").
    drop_ratio() {
        awk -v p1="$1" -v p2="$2" 'BEGIN {
            if (p1=="NA" || p2=="NA") { print "NA"; exit }
            p1+=0; p2+=0
            if (p2==0) { print "Inf"; exit }
            printf "%.2f", p1/p2
        }'
    }
    CT_DROP=$(drop_ratio "$CT1" "$CT2")
    GA_DROP=$(drop_ratio "$GA1" "$GA2")

    echo -e "${PREFIX}\t${GENOME}\t${NREADS}\t${MEANLEN}\t${CT1}\t${CT2}\t${CT_DROP}\t${GA1}\t${GA2}\t${GA_DROP}" >> "$SUMMARY"
done

### Copy results back
rsync -a --no-p "$TMPDIR/$OUTDIR/$PREFIX/" "$HOMEDIR/$OUTDIR/$PREFIX/"
chgrp -R dome "$HOMEDIR/$OUTDIR" 2>/dev/null || true

echo ""
echo "============================================================="
echo "Damage profiling summary for $PREFIX"
echo "============================================================="
column -t "$SUMMARY"

### Cleanup
set +u
module purge
set -u
rm -rf "$TMPDIR"

echo ""
echo "Done! Results in: $OUTDIR/$PREFIX/"
echo "  Per-genome dirs: $OUTDIR/$PREFIX/<GENOME>/"
echo "  Summary table:   $OUTDIR/${PREFIX}.damage_summary.tsv"
