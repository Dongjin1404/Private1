#!/bin/bash
#SBATCH --job-name=18varcall
#SBATCH --cpus-per-task=8
#SBATCH --mem=48GB
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err

### 18. Variant Calling + SNP alignment — per reference genome (SPAAM Ch.14 "Genotyping") =========================
# Reproduces the SPAAM "Genome Mapping" genotyping step (which uses GATK UnifiedGenotyper ->
# MultiVCFAnalyzer) with a modern bcftools toolchain. For ONE reference genome it gathers EVERY
# sample BAM that mapped to it across ALL datasets (ancient: picareiro, Maixner, Tett_a, Wibowo_a;
# modern: moderngut_genomemapp, modernsed_genomemapp), genotypes each sample, builds a per-sample
# consensus pseudo-genome, and combines them into a multi-sample SNP alignment for phylogenetics.
#
# This is the "reference axis" complement to scripts 15-17 (which array over SAMPLES): here we array
# over REFERENCE GENOMES, because a SNP alignment is one-reference x many-samples.
#
# USAGE:  sbatch --array 1-8 scripts/18_variant_calling.sh <genome_list> <MAPROOT> <REFDIR> <OUTDIR>
# EXAMPLE:
#   sbatch --array 1-8 scripts/18_variant_calling.sh \
#       00_accessions/phylo_genomes \
#       15_genome_mapping \
#       resources/gut_taxa_genomes \
#       18_variant_calling
#
# INPUTS:
#   - Dedup BAMs:  ${MAPROOT}/<DATASET>/bam/<SAMPLE>_<GENOME>.dedup.bam   (from 15_genome_mapping)
#   - Reference:   ${REFDIR}/<GENOME>.fa.gz                                (one of the 24 gut taxa)
#
# OUTPUTS per genome:
#   - <GENOME>/vcf/<DATASET>__<SAMPLE>.vcf.gz          per-sample variant calls (bcftools)
#   - <GENOME>/consensus/<DATASET>__<SAMPLE>.fa        per-sample consensus pseudo-genome (ref length)
#   - <GENOME>/<GENOME>.exclude.bed                    masked regions (dustmasker + cov-outlier + rRNA)
#   - <GENOME>/<GENOME>.wga.fasta                      whole-genome alignment (reference + passing samples)
#   - <GENOME>/<GENOME>.snpAlignment.fasta             SNP-only alignment (-> RAxML/MEGA, SPAAM Ch.15)
#   - <GENOME>/<GENOME>.snpAlignment.positions.txt     SNP index -> contig:pos map
#   - <GENOME>/<GENOME>.merged.vcf.gz                  merged multi-sample VCF (-> VCFtools/IQ-TREE)
#   - <GENOME>/<GENOME>.variant_summary.tsv            per-sample breadth / SNP / pass table
#   - phylo_snp_overview.tsv                           one row per genome (shared; flock-guarded)
#
# ---------------------------------------------------------------------------------------------------------------
# WHY bcftools (not GATK UnifiedGenotyper + MultiVCFAnalyzer):
#   - GATK3 UG + MVA are legacy Java tools; bcftools is current and on LISC.
#   - reproduce MVA's decision semantics EXACTLY at the consensus step (see "SNP FILTER" below).
#
# PLOIDY 1 (haploid) vs SPAAM's diploid+allele-frequency trick:
#   - Bacteria are haploid: one true base per position. SPAAM runs GATK at ploidy 2 only to expose
#     per-allele read fractions, then MVA decides the base by allele frequency. We call at ploidy 1
#     (biologically correct) and apply the SAME allele-fraction rule ourselves. Same outcome.
#
# SNP FILTER (== MultiVCFAnalyzer semantics):
#   - depth >= MINDP (3x)                      -> else position = N (missing)
#   - QUAL >= MINQUAL (30)
#   - major(ALT) allele fraction >= MINAF (0.90) -> call the base
#   - mixed allele fraction in (MIX_LO, MINAF) = (0.10, 0.90) -> N (uncertain: contamination/strain mix/damage)
#
# STRAIN-HETEROGENEITY REPORT (reporting only -- does NOT change any base call):
#   - the Het_SNPs_<lo>_<hi> column counts SNPs whose ALT fraction sits in [HET_LO, HET_HI]
#     (default 0.40-0.80), the mixed-allele band most diagnostic of within-sample strain mixture.
#     These are a SUBSET of the masked "ambiguous" positions; the column simply quantifies how many
#     variable sites were removed specifically because of strain heterogeneity.
#   - override the window: `HET_LO=0.30 HET_HI=0.70 sbatch ...`  (or a single value: `HET_LO=0.50 HET_HI=0.50`).
#
# OTTONI INCLUSION FILTER:
#   - a sample enters the alignment only if >= BREADTH_MIN (70%) of the genome is covered at >= MINDP (3x).
#   - failing samples are reported in the summary but excluded from the WGA/SNP alignment.
#
# UNIFORM MASKING (== SPAAM "exclude problematic regions"; identical for RefSeq + MAG references):
#   - dustmasker        low-complexity / simple repeats
#   - coverage-outlier  positions with merged depth >> genome median (collapsed multi-copy repeats / rRNA;
#                       these are the tall spikes in plot_coverage.R) -> COV_OUTLIER_FACTOR x median
#   - barrnap rRNA      ab-initio rRNA prediction (runs on any FASTA, so RefSeq and MAGs are treated the same)
#   masked positions become N in EVERY sample, so they can never contribute a SNP.
#
# DAMAGE (UDG-half libraries): optional terminal end-trim (TRIM_ENDS bp) via bamUtil `bam trimBam`
#   (sets trimmed-base qualities to 0 so the -Q base-quality filter ignores them). If bamUtil is not
#   available the script WARNS and proceeds untrimmed (the 3x + 0.90 + mixed->N filters still guard
#   against damage-derived false SNPs). For strict SPAAM-faithful rescaling, pre-run mapDamage --rescale.
#
# SNP EXTRACTION: snp-sites is not on LISC, so the SNP alignment is built with a small awk pass
#   (the "custom SNP-FASTA builder"): SNP positions = union of confident-ALT positions across passing
#   samples, then re-verified variable against the masked WGA (>= 2 distinct A/C/G/T). No dependency.
#
# DEPTH-CUTOFF NOTE (borderline genomes): run all genomes at MINDP=3 first. For the lower-coverage
#   genomes (SPIREOTU_02219533, SPIREOTU_02217808, GCA_900553245, GCA_937900695) inspect
#   phylo_snp_overview.tsv; if a genome yields < KEEP_MIN_SAMPLES passing samples or < KEEP_MIN_SNPS
#   SNP sites it is flagged LOW_SIGNAL. You can then re-run just those with `MINDP=2 sbatch ...` and
#   compare, otherwise keep the 3x result.
#
# CONVENTION NOTE on the two quality flags (they are SWAPPED between tools!):
#   bcftools mpileup -q = min MAPPING quality, -Q = min BASE quality
#   samtools  depth   -q = min BASE quality,    -Q = min MAPPING quality
# ===================================================================================================================

set -euo pipefail

### Variables
HOMEDIR=$(pwd)
GENOME_LIST="$1"
MAPROOT="${2%/}"
REFDIR="${3%/}"
OUTDIR="${4%/}"
GENOME=$(awk -v n="${SLURM_ARRAY_TASK_ID}" 'NR==n {print $1}' "$GENOME_LIST")

### Tunable parameters (env-overridable)
MINDP="${MINDP:-3}"                       # minimum depth to call a base (else N)
MINQUAL="${MINQUAL:-30}"                   # minimum site QUAL
MINAF="${MINAF:-0.90}"                     # major-allele fraction for a confident homozygous call
HET_LO="${HET_LO:-0.40}"                   # strain-heterogeneity reporting band: lower ALT-fraction bound
HET_HI="${HET_HI:-0.80}"                   # strain-heterogeneity reporting band: upper ALT-fraction bound
MIX_LO="${MIX_LO:-0.10}"                   # below this ALT fraction => treat as reference (no SNP)
MAPQ_MIN="${MAPQ_MIN:-20}"                 # min mapping quality (matches 15_genome_mapping.sh)
BASEQ_MIN="${BASEQ_MIN:-30}"              # min base quality
BREADTH_MIN="${BREADTH_MIN:-0.70}"        # Ottoni: fraction of genome at >= MINDP to include a sample
TRIM_ENDS="${TRIM_ENDS:-2}"               # terminal bp to trim for UDG-half damage (0 = disable)
COV_OUTLIER_FACTOR="${COV_OUTLIER_FACTOR:-5}"  # mask positions with merged depth > FACTOR x median
MASK_RRNA="${MASK_RRNA:-1}"               # 1 = run barrnap and mask rRNA loci
KEEP_MIN_SAMPLES="${KEEP_MIN_SAMPLES:-4}" # decision threshold for the overview "Decision" column
KEEP_MIN_SNPS="${KEEP_MIN_SNPS:-50}"      # decision threshold for the overview "Decision" column
MPILEUP_EXTRA="${MPILEUP_EXTRA:-}"        # e.g. set to "-B" to disable BAQ

# Integer breadth percentage for self-documenting column labels (tracks BREADTH_MIN, e.g. 0.30 -> 30).
BREADTH_PCT=$(awk -v b="$BREADTH_MIN" 'BEGIN{printf "%d", b*100}')
# Integer percent labels for the strain-heterogeneity band column name (e.g. 0.40,0.80 -> 40,80).
HET_LO_PCT=$(awk -v x="$HET_LO" 'BEGIN{printf "%d", x*100}')
HET_HI_PCT=$(awk -v x="$HET_HI" 'BEGIN{printf "%d", x*100}')

if [[ -z "$GENOME" ]]; then
    echo "ERROR: no genome at line ${SLURM_ARRAY_TASK_ID} of $GENOME_LIST"
    exit 1
fi

### Load modules (LISC-verified builds)
# The four core EasyBuild tools all sit on the GCCcore-12.3.0 toolchain (GCC-12.3.0 / gompi-2023a),
# so they co-load without a GCCcore conflict. dustmasker ships inside BLAST+.
# barrnap and bamUtil are provided ONLY as Conda environments on LISC (and on a different GCCcore),
# so they are NOT module-loaded here -- loading them would swap the GCCcore the core tools need.
# Instead they run inside isolated Conda subshells (see masking + per-sample trim below).
# Lmod's init script references unbound shell vars (LD_LIBRARY_PATH) and `conda activate` trips over
# PS1/etc -- both abort under `set -u`. Pre-seed LD_LIBRARY_PATH and relax nounset around every
# module/conda interaction (errexit stays on, so a genuinely missing module still aborts the job).
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
set +u
module purge
module load SAMtools/1.18-GCC-12.3.0
module load BCFtools/1.18-GCC-12.3.0
module load BEDTools/2.31.0-GCC-12.3.0
module load BLAST+/2.14.1-gompi-2023a          # gompi-2023a == GCCcore-12.3.0; provides dustmasker
set -u

# barrnap (rRNA) -> Conda env. Probe once; disable rRNA masking if the env is unavailable.
BARRNAP_ENV="barrnap-1.10.6"
if [[ "$MASK_RRNA" == "1" ]]; then
    if ! ( set +u; module load Conda >/dev/null 2>&1 && conda activate "$BARRNAP_ENV" >/dev/null 2>&1 \
           && command -v barrnap >/dev/null 2>&1 ); then
        echo "WARN: Conda env '$BARRNAP_ENV' unavailable -> disabling rRNA masking"
        MASK_RRNA=0
    fi
fi

# bamUtil (UDG-half end-trim) -> Conda env. Resolve the 'bam' binary path once; the per-sample loop
# calls it directly. If unavailable, the loop proceeds WITHOUT trimming (filters still guard damage).
BAMUTIL_BIN=""
if [[ "$TRIM_ENDS" -gt 0 ]]; then
    BAMUTIL_BIN=$( ( set +u; module load Conda >/dev/null 2>&1 && conda activate bamUtil-1.0.15 >/dev/null 2>&1 \
                     && command -v bam ) 2>/dev/null || true )
    [[ -n "$BAMUTIL_BIN" ]] || echo "WARN: Conda env 'bamUtil-1.0.15' unavailable -> proceeding WITHOUT end-trim"
fi

### Output + work directories
GOUT="$OUTDIR/$GENOME"
mkdir -p "$HOMEDIR/$GOUT"
mkdir -p "$TMPDIR/$GOUT/vcf" "$TMPDIR/$GOUT/consensus" "$TMPDIR/bams"
WORK="$TMPDIR/$GOUT"

echo "============================================================="
echo "Genome:    $GENOME"
echo "Map root:  $MAPROOT"
echo "Ref dir:   $REFDIR"
echo "Output:    $GOUT"
echo "Params:    MINDP=$MINDP MINQUAL=$MINQUAL MINAF=$MINAF MIX_LO=$MIX_LO"
echo "           MAPQ_MIN=$MAPQ_MIN BASEQ_MIN=$BASEQ_MIN BREADTH_MIN=$BREADTH_MIN"
echo "           TRIM_ENDS=$TRIM_ENDS COV_OUTLIER_FACTOR=$COV_OUTLIER_FACTOR MASK_RRNA=$MASK_RRNA"
echo "           HET_BAND=[${HET_LO}, ${HET_HI}]  (strain-heterogeneity reporting window)"
echo "============================================================="

### -----------------------------------------------------------------------
### Reference prep
### -----------------------------------------------------------------------
GFA="$WORK/${GENOME}.fa"
if [[ ! -f "$HOMEDIR/$REFDIR/${GENOME}.fa.gz" ]]; then
    echo "ERROR: reference $REFDIR/${GENOME}.fa.gz not found"
    exit 1
fi
zcat "$HOMEDIR/$REFDIR/${GENOME}.fa.gz" > "$GFA"
samtools faidx "$GFA"

# Region restriction: BBSplit BAMs carry @SQ headers for ALL 24 genomes' contigs. We restrict
# bcftools/samtools to THIS genome's contigs so the other 23 genomes' header entries are ignored.
# (We do NOT reheader: BAM records store the reference as an integer index into the @SQ list, so
#  dropping @SQ lines would corrupt those indices. Region restriction is the safe approach.)
REGIONS=$(cut -f1 "${GFA}.fai" | paste -sd, -)
GBED="$WORK/${GENOME}.contigs.bed"
awk -v OFS='\t' '{print $1, 0, $2}' "${GFA}.fai" > "$GBED"
GLEN=$(awk '{s+=$2} END {print s+0}' "${GFA}.fai")
echo "Reference: $(wc -l < "${GFA}.fai") contig(s), ${GLEN} bp total"

### -----------------------------------------------------------------------
### Discover sample BAMs for THIS genome across ALL datasets
### -----------------------------------------------------------------------
shopt -s nullglob
BAMS=( "$HOMEDIR/$MAPROOT"/*/bam/*_"${GENOME}".dedup.bam )
shopt -u nullglob
if [[ ${#BAMS[@]} -eq 0 ]]; then
    echo "ERROR: no BAMs matching $MAPROOT/*/bam/*_${GENOME}.dedup.bam"
    exit 1
fi
echo "Found ${#BAMS[@]} sample BAM(s) for $GENOME"

# Stage all BAMs (+ .bai) into TMPDIR for fast local I/O, recording NAME = <DATASET>__<SAMPLE>.
declare -a NAMES STAGED
BAMLIST="$WORK/bamlist.txt"; : > "$BAMLIST"
for B in "${BAMS[@]}"; do
    DSET=$(basename "$(dirname "$(dirname "$B")")")              # 15_genome_mapping/<DSET>/bam -> <DSET>
    SAMPLE=$(basename "$B" "_${GENOME}.dedup.bam")               # strips known suffix (handles underscores)
    NAME="${DSET}__${SAMPLE}"
    DST="$TMPDIR/bams/${NAME}.bam"
    rsync -a "$B" "$DST"
    if [[ -f "${B}.bai" ]]; then rsync -a "${B}.bai" "${DST}.bai"; else samtools index "$DST"; fi
    NAMES+=( "$NAME" )
    STAGED+=( "$DST" )
    echo "$DST" >> "$BAMLIST"
done

### -----------------------------------------------------------------------
### Build the uniform exclude.bed  (dustmasker + coverage-outlier + barrnap rRNA)
### -----------------------------------------------------------------------
echo ""
echo "=== Building exclude.bed (dustmasker + coverage-outlier + rRNA) ==="
EXCLUDE="$WORK/${GENOME}.exclude.bed"
: > "$EXCLUDE"

# (a) dustmasker low-complexity: emit soft-masked FASTA, parse lowercase runs -> BED (unambiguous).
dustmasker -in "$GFA" -outfmt fasta 2>/dev/null \
  | awk 'BEGIN{OFS="\t"}
         /^>/ { sub(/^>/,""); split($0,h," "); contig=h[1]; pos=0; inrun=0; next }
         {
           n=length($0)
           for(i=1;i<=n;i++){
             pos++
             c=substr($0,i,1)
             low=(c ~ /[acgtn]/)
             if(low && !inrun){ inrun=1; start=pos-1 }
             else if(!low && inrun){ print contig, start, pos-1; inrun=0 }
           }
         }
         END{ if(inrun) print contig, start, pos }' >> "$EXCLUDE" || true

# (b) coverage-outlier: merge all sample BAMs, compute per-position merged depth on this genome,
#     mask positions whose depth exceeds COV_OUTLIER_FACTOR x (median of covered positions).
MERGED="$WORK/merged_all.bam"
samtools merge -f -@ "$SLURM_CPUS_PER_TASK" -b "$BAMLIST" "$MERGED"
samtools index "$MERGED"
samtools depth -a -b "$GBED" -q "$BASEQ_MIN" -Q "$MAPQ_MIN" "$MERGED" > "$WORK/merged_depth.txt" || true
MED=$(awk '$3>0{print $3}' "$WORK/merged_depth.txt" | sort -n \
      | awk '{a[c++]=$1} END{if(c==0){print 0} else if(c%2){print a[int(c/2)]} else {printf "%.1f",(a[c/2-1]+a[c/2])/2}}')
if awk -v m="$MED" 'BEGIN{exit !(m>0)}'; then
    THRESH=$(awk -v m="$MED" -v f="$COV_OUTLIER_FACTOR" 'BEGIN{printf "%.4f", m*f}')
    awk -v t="$THRESH" -v OFS='\t' '$3>t{print $1, $2-1, $2}' "$WORK/merged_depth.txt" >> "$EXCLUDE" || true
    echo "  coverage-outlier: median=${MED}x threshold=${THRESH}x"
else
    echo "  coverage-outlier: median depth = 0 (no coverage); skipped"
fi

# (c) barrnap rRNA loci -> BED. Run inside a Conda subshell so the env change does not leak.
if [[ "$MASK_RRNA" == "1" ]]; then
    ( set +u; module load Conda >/dev/null 2>&1; conda activate "$BARRNAP_ENV" >/dev/null 2>&1; \
      barrnap --kingdom bac --quiet "$GFA" 2>/dev/null ) \
      | awk -v OFS='\t' '$0 !~ /^#/ && $3=="rRNA" {print $1, $4-1, $5}' >> "$EXCLUDE" || true
    echo "  rRNA: barrnap done"
fi

# Merge the exclude intervals (sort + collapse).
if [[ -s "$EXCLUDE" ]]; then
    sort -k1,1 -k2,2n "$EXCLUDE" | bedtools merge > "${EXCLUDE}.tmp" && mv "${EXCLUDE}.tmp" "$EXCLUDE"
fi
EXCL_BP=$(awk '{s+=$3-$2} END{print s+0}' "$EXCLUDE")
echo "  exclude.bed: $(wc -l < "$EXCLUDE") interval(s), ${EXCL_BP} bp masked uniformly"

### -----------------------------------------------------------------------
### Per-sample genotyping + consensus
### -----------------------------------------------------------------------
echo ""
echo "=== Per-sample genotyping + consensus (${#STAGED[@]} samples) ==="

SUMMARY="$HOMEDIR/$GOUT/${GENOME}.variant_summary.tsv"
echo -e "Genome\tDataset\tSample\tNreads\tCovBreadth_Pct\tNonN_Pct\tConfident_SNPs\tAmbiguous_Masked\tHet_SNPs_${HET_LO_PCT}_${HET_HI_PCT}\tPass${BREADTH_PCT}" > "$SUMMARY"

# AD-based ALT fraction for a normalised (biallelic) record: AD[0:1] = alt depth, AD[0:0] = ref depth.
AFEXPR='(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1])'

CAND_POS="$WORK/cand_positions.txt"; : > "$CAND_POS"   # union of confident-ALT positions (passing samples)
declare -a PASS_NAMES MERGE_VCFS
NPASS=0

for i in "${!STAGED[@]}"; do
    NAME="${NAMES[$i]}"
    DSET="${NAME%%__*}"
    SAMPLE="${NAME#*__}"
    SRC="${STAGED[$i]}"

    # Damage end-trim (UDG-half). bamUtil sets trimmed-base qualities to 0 -> ignored by -Q filter.
    CALLBAM="$SRC"
    if [[ -n "$BAMUTIL_BIN" ]]; then
        "$BAMUTIL_BIN" trimBam "$SRC" "$WORK/trim.bam" "$TRIM_ENDS" 2>/dev/null \
            && samtools index "$WORK/trim.bam" && CALLBAM="$WORK/trim.bam" || CALLBAM="$SRC"
    fi

    NREADS=$(samtools view -c "$SRC" 2>/dev/null || echo 0)

    # Depth (samtools depth: -q=BASE qual, -Q=MAPPING qual -- swapped vs mpileup!)
    samtools depth -a -b "$GBED" -q "$BASEQ_MIN" -Q "$MAPQ_MIN" "$CALLBAM" > "$WORK/depth.txt" || true
    COVERED=$(awk -v d="$MINDP" '$3>=d{c++} END{print c+0}' "$WORK/depth.txt")
    COVBR=$(awk -v c="$COVERED" -v g="$GLEN" 'BEGIN{printf "%.4f", (g>0)?100*c/g:0}')
    # low-coverage mask = genome MINUS well-covered positions (robust to missing zero-depth rows).
    awk -v d="$MINDP" -v OFS='\t' '$3>=d{print $1, $2-1, $2}' "$WORK/depth.txt" | sort -k1,1 -k2,2n | bedtools merge > "$WORK/covok.bed" || : > "$WORK/covok.bed"
    bedtools subtract -a "$GBED" -b "$WORK/covok.bed" > "$WORK/lowcov.bed" || cp "$GBED" "$WORK/lowcov.bed"

    # Genotype (bcftools mpileup -q=MAPPING qual, -Q=BASE qual), variants only, haploid, split multiallelic.
    bcftools mpileup -f "$GFA" -q "$MAPQ_MIN" -Q "$BASEQ_MIN" -r "$REGIONS" \
        -a FORMAT/AD,FORMAT/DP $MPILEUP_EXTRA -Ou "$CALLBAM" 2>/dev/null \
      | bcftools call -mv --ploidy 1 -Oz -o "$WORK/calls.raw.vcf.gz" 2>/dev/null
    bcftools index -f "$WORK/calls.raw.vcf.gz"
    bcftools norm -f "$GFA" -m- "$WORK/calls.raw.vcf.gz" -Oz -o "$WORK/calls.vcf.gz" 2>/dev/null
    bcftools index -f "$WORK/calls.vcf.gz"

    # Confident ALT SNPs (== call the base) and ambiguous/mixed SNPs (== mask to N).
    bcftools view -v snps -i "FORMAT/DP>=${MINDP} && QUAL>=${MINQUAL} && ${AFEXPR}>=${MINAF}" \
        "$WORK/calls.vcf.gz" -Oz -o "$WORK/confident.vcf.gz" 2>/dev/null
    bcftools index -f "$WORK/confident.vcf.gz"
    NCONF=$(bcftools view -H "$WORK/confident.vcf.gz" 2>/dev/null | wc -l | tr -d ' ')

    bcftools view -v snps -i "FORMAT/DP>=${MINDP} && ${AFEXPR}>${MIX_LO} && ${AFEXPR}<${MINAF}" \
        "$WORK/calls.vcf.gz" 2>/dev/null \
      | bcftools query -f '%CHROM\t%POS\n' 2>/dev/null \
      | awk -v OFS='\t' '{print $1, $2-1, $2}' > "$WORK/ambig.bed" || : > "$WORK/ambig.bed"
    NAMBIG=$(wc -l < "$WORK/ambig.bed" | tr -d ' ')

    # Strain-heterogeneity band [HET_LO, HET_HI] (default 40-80%): SNPs whose ALT fraction lands in the
    # mixed-allele window most diagnostic of within-sample strain mixture. A subset of the ambiguous
    # positions above; reported (not re-filtered) so we can quantify variable sites lost to heterogeneity.
    NHET=$(bcftools view -v snps -i "FORMAT/DP>=${MINDP} && ${AFEXPR}>=${HET_LO} && ${AFEXPR}<=${HET_HI}" \
        "$WORK/calls.vcf.gz" -H 2>/dev/null | wc -l | tr -d ' ')

    # Final mask = exclude (uniform) U low-coverage U ambiguous  -> N in the consensus.
    cat "$EXCLUDE" "$WORK/lowcov.bed" "$WORK/ambig.bed" 2>/dev/null \
      | sort -k1,1 -k2,2n | bedtools merge > "$WORK/mask.bed" || : > "$WORK/mask.bed"

    # Consensus pseudo-genome (apply confident ALTs, mask uncertain -> N). Concatenate contigs into
    # one record so every sample is a single sequence of length GLEN (alignment columns line up).
    bcftools consensus -f "$GFA" -m "$WORK/mask.bed" --mask-with N -H 1 "$WORK/confident.vcf.gz" \
        > "$WORK/cons_multi.fa" 2>/dev/null
    { echo ">$NAME"; awk '!/^>/{printf "%s",$0} END{print ""}' "$WORK/cons_multi.fa"; } > "$WORK/consensus/${NAME}.fa"
    NONN=$(awk -v g="$GLEN" '!/^>/{for(i=1;i<=length($0);i++){c=substr($0,i,1); if(c!="N"&&c!="n") k++}} END{printf "%.4f",(g>0)?100*k/g:0}' "$WORK/consensus/${NAME}.fa")

    # Keep the per-sample calls (rename sample to NAME for a clean merged VCF later).
    echo "$NAME" > "$WORK/sn.txt"
    bcftools reheader -s "$WORK/sn.txt" "$WORK/calls.vcf.gz" -o "$WORK/vcf/${NAME}.vcf.gz"
    bcftools index -f "$WORK/vcf/${NAME}.vcf.gz"

    # Ottoni inclusion decision.
    if awk -v b="$COVBR" -v t="$BREADTH_MIN" 'BEGIN{exit !((b/100)>=t)}'; then
        PASS="yes"; NPASS=$((NPASS+1))
        PASS_NAMES+=( "$NAME" )
        MERGE_VCFS+=( "$WORK/vcf/${NAME}.vcf.gz" )
        bcftools query -f '%CHROM\t%POS\n' "$WORK/confident.vcf.gz" 2>/dev/null >> "$CAND_POS" || true
    else
        PASS="no"
    fi

    printf "%-28s reads=%-8s cov%%@%dx=%-8s nonN%%=%-8s SNP=%-6s ambig=%-6s het=%-6s %s\n" \
        "$NAME" "$NREADS" "$MINDP" "$COVBR" "$NONN" "$NCONF" "$NAMBIG" "$NHET" "$PASS"
    echo -e "${GENOME}\t${DSET}\t${SAMPLE}\t${NREADS}\t${COVBR}\t${NONN}\t${NCONF}\t${NAMBIG}\t${NHET}\t${PASS}" >> "$SUMMARY"

    rm -f "$WORK/trim.bam" "$WORK/trim.bam.bai" "$WORK/calls.raw.vcf.gz"* "$WORK/calls.vcf.gz"* \
          "$WORK/confident.vcf.gz"* "$WORK/cons_multi.fa" "$WORK/depth.txt" \
          "$WORK/covok.bed" "$WORK/lowcov.bed" "$WORK/ambig.bed" "$WORK/mask.bed"
done

echo ""
echo "Passing samples (>= ${BREADTH_MIN} breadth @ ${MINDP}x): $NPASS / ${#STAGED[@]}"

### -----------------------------------------------------------------------
### Assemble WGA + SNP alignment + merged VCF
### -----------------------------------------------------------------------
WGA="$WORK/${GENOME}.wga.fasta"
SNPFA="$WORK/${GENOME}.snpAlignment.fasta"
SNPPOS="$WORK/${GENOME}.snpAlignment.positions.txt"

# Reference baseline tip (SPAAM/MVA includes the reference in the SNP table).
{ echo ">REFERENCE_${GENOME}"; awk '!/^>/{printf "%s",$0} END{print ""}' "$GFA"; } > "$WGA"
for NAME in "${PASS_NAMES[@]:-}"; do
    [[ -n "$NAME" ]] && cat "$WORK/consensus/${NAME}.fa" >> "$WGA"
done

NSNP=0
if [[ "$NPASS" -ge 1 && -s "$CAND_POS" ]]; then
    # Candidate concat indices = confident-ALT positions mapped to the concatenated coordinate.
    sort -u "$CAND_POS" > "$WORK/cand_uniq.txt"
    awk '
        FNR==NR { ord[++n]=$1; len[$1]=$2; next }
        !init   { off=0; for(i=1;i<=n;i++){offset[ord[i]]=off; off+=len[ord[i]]}; init=1 }
        { print (offset[$1]+$2)"\t"$1"\t"$2 }
    ' "${GFA}.fai" "$WORK/cand_uniq.txt" | sort -n -u -k1,1 > "$WORK/cand_idx.txt"

    # Extract variable columns (re-verify >= 2 distinct A/C/G/T against the masked WGA).
    awk -v snpout="$SNPFA" -v posout="$SNPPOS" '
        FNR==NR { ci[++m]=$1; cc[m]=$2; cp[m]=$3; next }
        /^>/    { name=substr($0,2); names[++ns]=name; seq[ns]=""; next }
                { seq[ns]=seq[ns] $0 }
        END {
            kept=0
            for(j=1;j<=m;j++){
                idx=ci[j]; delete seen; distinct=0
                for(s=1;s<=ns;s++){
                    b=toupper(substr(seq[s],idx,1))
                    if(b=="A"||b=="C"||b=="G"||b=="T"){ if(!(b in seen)){seen[b]=1;distinct++} }
                }
                if(distinct>=2){ kept++; kidx[kept]=idx; kcon[kept]=cc[j]; kpos[kept]=cp[j] }
            }
            print "snp_index\tcontig\tref_pos" > posout
            for(k=1;k<=kept;k++) print k"\t"kcon[k]"\t"kpos[k] > posout
            for(s=1;s<=ns;s++){
                printf(">%s\n", names[s]) > snpout
                line=""
                for(k=1;k<=kept;k++) line=line substr(seq[s],kidx[k],1)
                print line > snpout
            }
            print kept > "/dev/stderr"
        }
    ' "$WORK/cand_idx.txt" "$WGA" 2> "$WORK/nsnp.txt"
    NSNP=$(cat "$WORK/nsnp.txt" 2>/dev/null || echo 0)
else
    { echo ">REFERENCE_${GENOME}"; echo ""; } > "$SNPFA"
    echo -e "snp_index\tcontig\tref_pos" > "$SNPPOS"
    echo "  (no passing samples or no candidate SNPs -> empty SNP alignment)"
fi
echo "SNP alignment: ${NSNP} variable site(s) across $((NPASS+1)) sequences (incl. reference)"

# Merged multi-sample VCF (for VCFtools / IQ-TREE routes). Use NPASS (== number of passing samples,
# one MERGE_VCFS entry each) so we never take the length of a possibly-empty array under `set -u`
# (bash < 4.4 throws "unbound variable" for ${#arr[@]} when arr is empty).
if [[ "$NPASS" -ge 2 ]]; then
    bcftools merge -Oz -o "$WORK/${GENOME}.merged.vcf.gz" "${MERGE_VCFS[@]}" 2>/dev/null \
        && bcftools index -f "$WORK/${GENOME}.merged.vcf.gz" \
        && echo "Merged VCF: ${GENOME}.merged.vcf.gz ($NPASS samples)"
elif [[ "$NPASS" -eq 1 ]]; then
    cp "${MERGE_VCFS[0]}" "$WORK/${GENOME}.merged.vcf.gz"
    cp "${MERGE_VCFS[0]}.csi" "$WORK/${GENOME}.merged.vcf.gz.csi" 2>/dev/null || true
fi

# (exclude.bed already lives at $WORK/${GENOME}.exclude.bed == $EXCLUDE; rsync ships it, no copy needed.)

### -----------------------------------------------------------------------
### Per-genome row in the shared overview (flock-guarded: array tasks append concurrently)
### -----------------------------------------------------------------------
OVERVIEW="$HOMEDIR/$OUTDIR/phylo_snp_overview.tsv"
OVLOCK="$HOMEDIR/$OUTDIR/.overview.lock"
if awk -v s="$NPASS" -v n="$NSNP" -v ms="$KEEP_MIN_SAMPLES" -v mn="$KEEP_MIN_SNPS" \
       'BEGIN{exit !(s>=ms && n>=mn)}'; then DECISION="KEEP"; else DECISION="LOW_SIGNAL"; fi
(
    flock -x 201
    [[ -f "$OVERVIEW" ]] || echo -e "Genome\tMINDP\tSamples_Total\tSamples_Pass${BREADTH_PCT}\tRef_Len_bp\tExcluded_bp\tSNP_Sites\tDecision" > "$OVERVIEW"
    echo -e "${GENOME}\t${MINDP}\t${#STAGED[@]}\t${NPASS}\t${GLEN}\t${EXCL_BP}\t${NSNP}\t${DECISION}" >> "$OVERVIEW"
) 201>"$OVLOCK"

### Finalize
chgrp -R dome "$TMPDIR/$GOUT" 2>/dev/null || true
rsync -a --no-p "$TMPDIR/$GOUT/" "$HOMEDIR/$GOUT/"
chgrp -R dome "$HOMEDIR/$GOUT" 2>/dev/null || true

echo ""
echo "============================================================="
echo "Variant-calling summary for $GENOME"
echo "============================================================="
column -t "$SUMMARY" | head -n 40
echo "  ... ($(($(wc -l < "$SUMMARY")-1)) samples total)"
echo ""
echo "Decision: $DECISION  (passing=$NPASS/${#STAGED[@]}, SNP sites=$NSNP)"

### Cleanup
set +u
module purge
rm -rf "$TMPDIR"

echo ""
echo "Done! Results in: $GOUT/"
echo "  SNP alignment: $GOUT/${GENOME}.snpAlignment.fasta   (-> RAxML/MEGA)"
echo "  Full WGA:      $GOUT/${GENOME}.wga.fasta            (-> IQ-TREE / VCFtools route)"
echo "  Merged VCF:    $GOUT/${GENOME}.merged.vcf.gz"
echo "  Summary:       $GOUT/${GENOME}.variant_summary.tsv"
echo "  Overview row:  $OUTDIR/phylo_snp_overview.tsv"