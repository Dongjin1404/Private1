#!/bin/bash
#SBATCH --job-name=20ismA
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

### 20. ismA functional screen — coprostanol-production capacity via DIAMOND blastx ================================
# Screens the merged single-end reads (from 04_bbmerge) for the ismA gene, the enzyme that converts
# cholesterol -> coprostanol in the gut (Kenny et al. 2020, Cell Host & Microbe 28:245-257,
# "Cholesterol Metabolism by Uncultured Human Gut Bacteria Influences Host Cholesterol Level").
#
# WHY A FUNCTIONAL (not taxonomic) SCREEN:
#   The sylph taxonomic profile of the picareiro coprolites contained NONE of the cultured coprostanol
#   producers (Phocaeicola dorei / Eubacterium coprostanoligenes / Sterolibacterium all = 0). ismA is
#   carried overwhelmingly by UNCULTURED gut taxa, so a species-level profile cannot see it and there is
#   nothing to "pre-filter" to. We therefore search ALL reads directly against an ismA PROTEIN reference
#   with a translated (blastx) search, which (a) catches divergent/uncultured ismA clades and (b) is more
#   robust to aDNA C->T deamination than a nucleotide search (many base changes are synonymous at protein
#   level, and --frameshift tolerates indels).
#
# USAGE:  sbatch --array 1-N scripts/20_ismA_screen.sh <accesslist> <MERGED_DIR> <OUTDIR> <ismA_proteins.faa>
#   <accesslist>        00_accessions/<dataset>            (one sample ID per line; array indexes lines)
#   <MERGED_DIR>        04_bbmerge/<dataset>               (holds <PREFIX>.merged.fastq.gz)
#   <OUTDIR>            20_ismA_screen                     (SHARED across datasets; per-sample files go in
#                                                          <OUTDIR>/<dataset>/, overview is <OUTDIR>/...tsv)
#   <ismA_proteins.faa> resources/ismA/ismA_proteins.faa  (protein FASTA of ismA references; see NOTE below)
#
# EXAMPLE (run each dataset as its own array job; they all accumulate into one shared overview):
#   sbatch --array 1-3 scripts/20_ismA_screen.sh 00_accessions/picareiro   04_bbmerge/picareiro   20_ismA_screen resources/ismA/ismA_proteins.faa
#   sbatch --array 1-3 scripts/20_ismA_screen.sh 00_accessions/controls    04_bbmerge/controls    20_ismA_screen resources/ismA/ismA_proteins.faa
#   sbatch --array 1-N scripts/20_ismA_screen.sh 00_accessions/moderngut   04_bbmerge/moderngut   20_ismA_screen resources/ismA/ismA_proteins.faa
#   ... (Wibowo_a, Maixner, modernsed as biological / environmental comparisons)
#
# NOTE ON THE REFERENCE FASTA (must be supplied by the user, ONCE):
#   Obtain the IsmA protein sequences from the Kenny et al. 2020 supplement (they form a few ismA sequence
#   clades) and save them as resources/ismA/ismA_proteins.faa. Include several clade representatives so the
#   database spans the known ismA diversity (reduces false negatives on divergent ancient carriers). The
#   DIAMOND database (.dmnd) is built automatically from this FASTA on first run (lock-guarded).
#
# INPUTS  (per sample):
#   - ${MERGED_DIR}/${PREFIX}.merged.fastq.gz            merged SE reads from 04_bbmerge
#
# OUTPUTS (per sample, in <OUTDIR>/<dataset>/):
#   - ${PREFIX}.ismA.blastx.tsv.gz     full DIAMOND hit table (fmt6 + qlen + full_qseq)
#   - ${PREFIX}.ismA.hits.faa.gz       hit reads as FASTA (for optional nucleotide-level mapDamage follow-up)
#   - ${PREFIX}.ismA_summary.tsv       one-row per-sample summary
# Shared across all datasets (flock-guarded):
#   - <OUTDIR>/ismA_screen_overview.tsv   one row per sample -> the cross-dataset comparison table
#
# INTERPRETATION:
#   - CONTROLS (SED2/SED3/EXTBL, modernsed) must be ~zero. Non-zero hits in a blank => contamination flag.
#   - Compare HitsPerMillion across samples (depth-normalised). moderngut = positive anchor.
#   - MedianHitReadLen should match the sample's aDNA fragment-length regime for ancient samples; a hit-read
#     length distribution that looks modern in an "ancient" sample is a red flag (dump FASTA -> mapDamage).
# ===================================================================================================================

set -euo pipefail

### Variables
HOMEDIR=$(pwd)
ACCESS="$1"
DATASET=$(basename "$1")
MERGED_DIR="${2%/}"
OUTDIR="${3%/}"
ISMA_FAA="$4"
PREFIX=$(awk '{if (NR=='${SLURM_ARRAY_TASK_ID}') print $1}' "$ACCESS")

### Tunable parameters (env-overridable)
MIN_PID="${MIN_PID:-60}"          # min amino-acid % identity of a hit (permissive for divergent/ancient ismA)
MIN_AALEN="${MIN_AALEN:-15}"      # min alignment length in amino acids (short aDNA reads -> keep modest)
EVALUE="${EVALUE:-1e-5}"          # DIAMOND e-value cutoff
SENS="${SENS:---sensitive}"       # sensitivity mode (e.g. --very-sensitive for maximum recall)
FRAMESHIFT="${FRAMESHIFT:-15}"    # frameshift penalty; enables indel/damage-tolerant alignment
DIAMOND_MODULE="${DIAMOND_MODULE:-DIAMOND/2.1.8-GCC-12.3.0}"

if [[ -z "$PREFIX" ]]; then
    echo "ERROR: no sample at line ${SLURM_ARRAY_TASK_ID} of $ACCESS"
    exit 1
fi

### Load modules ---------------------------------------------------------------------------------------------
# Lmod init references unbound LD_LIBRARY_PATH -> aborts under `set -u`; pre-seed it and relax nounset around
# every module call (errexit stays on). Standard LISC pattern (see 16/17/18).
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
set +u
module purge
module load "$DIAMOND_MODULE"
set -u

if ! command -v diamond >/dev/null 2>&1; then
    echo "ERROR: 'diamond' not on PATH after loading $DIAMOND_MODULE"
    echo "       Override with:  DIAMOND_MODULE=DIAMOND/<ver> sbatch ..."
    exit 1
fi
diamond version || true

### Locate + validate inputs ---------------------------------------------------------------------------------
READS="$HOMEDIR/$MERGED_DIR/${PREFIX}.merged.fastq.gz"
if [[ ! -s "$READS" ]]; then
    # fall back to an uncompressed merged file if that is how the dataset was written
    if [[ -s "$HOMEDIR/$MERGED_DIR/${PREFIX}.merged.fastq" ]]; then
        READS="$HOMEDIR/$MERGED_DIR/${PREFIX}.merged.fastq"
    else
        echo "ERROR: merged reads not found: $HOMEDIR/$MERGED_DIR/${PREFIX}.merged.fastq[.gz]"
        exit 1
    fi
fi

if [[ ! -s "$HOMEDIR/$ISMA_FAA" && ! -s "$ISMA_FAA" ]]; then
    echo "ERROR: ismA reference FASTA not found: $ISMA_FAA"
    echo "       Supply the IsmA protein sequences (Kenny et al. 2020) at resources/ismA/ismA_proteins.faa"
    exit 1
fi
[[ -s "$HOMEDIR/$ISMA_FAA" ]] && ISMA_FAA="$HOMEDIR/$ISMA_FAA"

### Build the DIAMOND protein DB once (lock-guarded so parallel array tasks don't collide) -------------------
DBBASE="$HOMEDIR/resources/ismA/ismA_proteins"   # diamond appends .dmnd
DB="${DBBASE}.dmnd"
mkdir -p "$(dirname "$DBBASE")"
(
    flock -x 210
    if [[ ! -s "$DB" || "$ISMA_FAA" -nt "$DB" ]]; then
        echo "Building DIAMOND DB from $ISMA_FAA ..."
        diamond makedb --in "$ISMA_FAA" -d "$DBBASE" -p "${SLURM_CPUS_PER_TASK}"
        chgrp -R dome "$(dirname "$DBBASE")" 2>/dev/null || true
    else
        echo "DIAMOND DB already present: $DB"
    fi
) 210>"$HOMEDIR/resources/ismA/.dbbuild.lock"

### Output + work directories -------------------------------------------------------------------------------
DOUT="$OUTDIR/$DATASET"
mkdir -p "$HOMEDIR/$DOUT"
mkdir -p "$TMPDIR/$PREFIX"
WORK="$TMPDIR/$PREFIX"

echo "============================================================="
echo "Sample:     $PREFIX"
echo "Dataset:    $DATASET"
echo "Reads:      $READS"
echo "ismA DB:    $DB"
echo "Filters:    pident>=${MIN_PID}  aln_aa>=${MIN_AALEN}  e<=${EVALUE}  ${SENS}  frameshift=${FRAMESHIFT}"
echo "============================================================="

### DIAMOND blastx (translated search: reads -> ismA proteins) ----------------------------------------------
# -k 1 keeps the single best hit per read (we are counting reads that carry ismA, not annotating every HSP).
# Output columns: qseqid sseqid pident length evalue bitscore qlen full_qseq
#   qlen      = read length in nt (for the hit-read length distribution / aDNA sanity check)
#   full_qseq = the read sequence (lets us dump hit reads to FASTA without re-reading the fastq)
RAW="$WORK/${PREFIX}.ismA.blastx.tsv"
diamond blastx \
    -d "$DBBASE" \
    -q "$READS" \
    -o "$RAW" \
    -p "${SLURM_CPUS_PER_TASK}" \
    -f 6 qseqid sseqid pident length evalue bitscore qlen full_qseq \
    -e "$EVALUE" $SENS --frameshift "$FRAMESHIFT" -k 1

### Filter hits + build the hit-read FASTA ------------------------------------------------------------------
# Keep rows passing identity + aligned-length thresholds. -k1 gives one row per read, but sort -u on the
# read id defends against any duplicates. Emit: (a) filtered table, (b) hit reads as FASTA.
FILT="$WORK/${PREFIX}.ismA.filtered.tsv"
awk -F'\t' -v pid="$MIN_PID" -v al="$MIN_AALEN" \
    '($3+0)>=pid && ($4+0)>=al' "$RAW" | sort -u -k1,1 > "$FILT"

HITFAA="$WORK/${PREFIX}.ismA.hits.faa"
awk -F'\t' '{print ">"$1"\n"$8}' "$FILT" > "$HITFAA"

### Per-sample metrics --------------------------------------------------------------------------------------
median_col() {  # $1=file $2=1-based column ; prints median of numeric column, or NA
    awk -F'\t' -v c="$2" 'NF{v[n++]=$c+0} END{
        if(n==0){print "NA"; exit}
        asort(v)
        if(n%2==1) printf "%.2f", v[(n+1)/2]
        else       printf "%.2f", (v[n/2]+v[n/2+1])/2
    }' "$1" 2>/dev/null || echo "NA"
}

# gawk asort may be unavailable; fall back to a sort|awk median if the above prints nothing.
median_safe() {  # $1=file $2=column
    local m; m=$(median_col "$1" "$2")
    if [[ -z "$m" || "$m" == "NA" ]]; then
        m=$(cut -f"$2" "$1" | sort -n | awk '{a[c++]=$1} END{
            if(c==0){print "NA"} else if(c%2==1){printf "%.2f",a[int(c/2)]}
            else{printf "%.2f",(a[c/2-1]+a[c/2])/2}}')
    fi
    echo "$m"
}

TOTALREADS=$( { zcat -f "$READS" 2>/dev/null || cat "$READS"; } | awk 'END{print NR/4}')
HITREADS=$(wc -l < "$FILT" | tr -d ' ')
UNIQREFS=$(cut -f2 "$FILT" | sort -u | grep -c . || true)
MEDPID=$(median_safe "$FILT" 3)
MEDLEN=$(median_safe "$FILT" 7)
HPM=$(awk -v h="$HITREADS" -v t="$TOTALREADS" 'BEGIN{ if(t>0) printf "%.3f", h/t*1e6; else print "NA"}')

### Per-sample summary + shared overview --------------------------------------------------------------------
SUMMARY="$WORK/${PREFIX}.ismA_summary.tsv"
echo -e "Sample\tDataset\tTotalReads\tHitReads\tHitsPerMillion\tUniqRefs\tMedianPctID\tMedianHitLen" > "$SUMMARY"
echo -e "${PREFIX}\t${DATASET}\t${TOTALREADS}\t${HITREADS}\t${HPM}\t${UNIQREFS}\t${MEDPID}\t${MEDLEN}" >> "$SUMMARY"

OVERVIEW="$HOMEDIR/$OUTDIR/ismA_screen_overview.tsv"
OVLOCK="$HOMEDIR/$OUTDIR/.ismA_overview.lock"
(
    flock -x 202
    [[ -f "$OVERVIEW" ]] || echo -e "Sample\tDataset\tTotalReads\tHitReads\tHitsPerMillion\tUniqRefs\tMedianPctID\tMedianHitLen" > "$OVERVIEW"
    echo -e "${PREFIX}\t${DATASET}\t${TOTALREADS}\t${HITREADS}\t${HPM}\t${UNIQREFS}\t${MEDPID}\t${MEDLEN}" >> "$OVERVIEW"
) 202>"$OVLOCK"

### Finalize: compress + ship back --------------------------------------------------------------------------
gzip -f "$RAW" "$HITFAA"
chgrp -R dome "$WORK" 2>/dev/null || true
rsync -a --no-p "$WORK/${PREFIX}.ismA.blastx.tsv.gz" "$HOMEDIR/$DOUT/" 2>/dev/null || true
rsync -a --no-p "$WORK/${PREFIX}.ismA.hits.faa.gz"   "$HOMEDIR/$DOUT/" 2>/dev/null || true
rsync -a --no-p "$SUMMARY"                            "$HOMEDIR/$DOUT/" 2>/dev/null || true
chgrp -R dome "$HOMEDIR/$DOUT" 2>/dev/null || true

echo ""
echo "============================================================="
echo "ismA screen done for $PREFIX ($DATASET)"
echo "============================================================="
column -t "$SUMMARY"
echo ""
echo "  Hit table:   $DOUT/${PREFIX}.ismA.blastx.tsv.gz"
echo "  Hit reads:   $DOUT/${PREFIX}.ismA.hits.faa.gz   (for optional mapDamage authentication)"
echo "  Summary:     $DOUT/${PREFIX}.ismA_summary.tsv"
echo "  Overview:    $OUTDIR/ismA_screen_overview.tsv   (shared across datasets)"

### Cleanup
set +u
module purge
set -u
rm -rf "$TMPDIR"
