# What `18_variant_calling.sh` Does — Overview

**Script:** `18_variant_calling.sh`
**Pipeline stage:** SPAAM Ch.14 "Genotyping" → produces the SNP alignment used to build the phylogenetic tree (script 19)
**One sentence:** For one reference genome, it takes every sample's mapped reads, decides the most reliable base at every genome position for each sample, and assembles those per-sample sequences into a single SNP alignment for phylogenetics.

---

## 1. Visual scheme of the steps

```mermaid
flowchart TD
    REF["Reference genome\n(e.g. GCF_039521565, 2.94 Mbp)"]
    BAMS["All sample BAMs for this genome\nancient: picareiro, Maixner, Tett\nmodern: gut samples\n(from script 15 mapping)"]

    REF --> MASK
    BAMS --> MASK

    MASK["STEP A — Build ONE uniform exclusion mask\n(applied identically to every sample)\n• repetitive / low-complexity regions (dustmasker)\n• rRNA genes (barrnap)\n• abnormal coverage spikes (collapsed repeats)"]

    MASK --> LOOP

    subgraph LOOP["STEP B — For EACH sample, independently"]
        direction TB
        P1["Count reads at every genome position"]
        P1 --> P2{"Depth >= 3 reads?\n(and not in masked region)"}
        P2 -->|"No / masked"| N1["position = N\n(missing / unknown)"]
        P2 -->|"Yes"| P3{"Do >= 90% of reads\nagree on one base?"}
        P3 -->|"Yes, >=90%"| CALL["call that base\n(A/C/G/T)"]
        P3 -->|"mixed 10-90%"| N2["position = N\n(uncertain: damage / contamination)"]
        N1 & N2 & CALL --> PG["Per-sample pseudo-genome\nFULL reference length\nA/C/G/T where confident, N elsewhere"]
        PG --> BF{"STEP C — Breadth filter\n>= 30% of genome\ncovered at >= 3x?"}
        BF -->|"FAIL"| EXC["excluded from the tree\n(still reported in summary)"]
        BF -->|"PASS"| INC["included"]
    end

    INC --> STACK["STEP D — Stack reference + all passing samples\n(every sequence same length, positions aligned)"]
    STACK --> SNP["STEP E — Keep only SNP columns\n(positions where >= 2 different A/C/G/T\nappear across the sequences,\nreference tip included in the count)"]
    SNP --> OUT["snpAlignment.fasta\n= input to RAxML (script 19)"]
    SNP --> SUM["variant_summary.tsv\nConfident SNPs vs Ambiguous (N) per sample"]
    SNP --> OV["phylo_snp_overview.tsv\nKEEP / LOW_SIGNAL decision per genome"]
```

---

## 2. Step-by-step in plain language

| Step | What happens | Why it matters |
|------|--------------|----------------|
| **A. Uniform mask** | Mark unreliable genome regions (repeats, rRNA, coverage spikes) **once**, then blank them in **every** sample. | Guarantees all samples are judged on the same fair set of positions; stops repeats from creating false SNPs. |
| **B. Per-sample base calling** | For each sample, at each position: require ≥3 reads and ≥90% agreement to call a base; otherwise write **N**. | This is the core variant call. The 90% rule and N-masking protect against ancient-DNA damage and contamination. |
| **C. Breadth filter** | A sample is kept only if **≥30%** of the genome is covered ≥3×. | Samples with too little data are excluded so they can't add noise. |
| **D. Stacking** | The reference and all passing samples are lined up — because all are reference length, positions already correspond. | This is the "alignment" — no gap-insertion software needed (mapping already aligned everything to the reference). |
| **E. SNP extraction** | Keep only the columns where the sequences actually differ (≥2 distinct A/C/G/T), **counting the reference tip together with the passing samples**. | These variable positions carry the phylogenetic signal; invariant positions are dropped to make the tree computation efficient. Because the reference is a genuine tip in the tree, a column where all samples share one base but the reference differs is a real SNP (it defines the reference-vs-samples branch) — so the reference must be included in the ≥2 count, not excluded. |

**Key output:** `snpAlignment.fasta` — one row per sample (+ reference), one column per SNP, with `N` wherever a sample lacked reliable data.

---

## 3. The coverage / low-breadth concern — answered point by point

> *Context: picareiro covers only ~30% of this genome. Does that distort the SNP calls and the tree?*

### Q1. "Is it an alignment of SNPs where regions may not be shared between samples?"
**Yes — and that is expected and handled correctly.** Each sample's sequence spans the **entire** reference length. Where a sample has no reads, that position is **`N` (missing data)**, *not* the reference base. So the alignment honestly records "this sample has no information here" rather than pretending it matches.

### Q2. "Is the 30% covered *consistently* (same region every time)?"
**Largely yes, and not by chance.** The covered ~30% is driven by which parts of this organism's genome are (a) conserved enough for reads to map and (b) not in the uniform mask. Because the **same exclusion mask** and the **same mapping criteria** apply to every sample, the regions that survive tend to be the **same well-behaved regions** across samples — not a random different 30% per sample. You can verify this directly (see §4).

### Q3. "Do the other samples also cover that region?"
**Not necessarily — and here is exactly what happens to a column where only one sample has data.** The retained-column test is "≥2 distinct A/C/G/T across **reference + passing samples**." Because the **reference is a complete genome with no `N`**, a position where exactly one sample has a base (say an ALT `G`) and every other sample is `N` still has two distinct bases (reference `A` + sample `G`) — so it is **KEPT, not dropped**. Such a column is a **singleton / autapomorphy**: it adds to that one sample's **terminal branch length** but **cannot change topology** (a single differing tip cannot group tips).

**Consequence to be aware of:** the mutation is real (that sample genuinely carries it), so it is not an artifact — but the *number* of these private columns scales with **breadth**, so a high-coverage sample accumulates more of them and its terminal branch looks longer purely because more of its genome was visible. This biases **branch-length comparability**, never topology. If branch-length comparison across samples of very different breadth matters, apply a **minimum taxon-occupancy filter** (drop columns with real data in fewer than *k* samples) or restrict to **parsimony-informative** columns (variant shared by ≥2 samples). Both are compatible with the `ASC_LEWIS` correction because they remove columns, not taxa. See §4 for the exact site counts at each threshold. The phylogenetically **informative** columns are, by definition, the ones with overlapping coverage between samples.

### Q4. "Is the tree built only on the covered region, or on the entire genome? (could it push samples artificially apart?)"
**Only on the positions with actual data — the model explicitly ignores `N`.** Two important safeguards:

1. **`N` = missing, not mismatch.** RAxML treats `N` as "no information." A sample is **never** penalised (pushed away on the tree) for positions it didn't cover. Distances are computed **only** over positions where both sequences being compared have real bases.
2. **Ascertainment-bias correction (`ASC_LEWIS`) in script 19.** Because we feed only SNP columns (not the whole genome), branch lengths would otherwise be inflated. The Lewis correction mathematically compensates for the dropped invariant positions, so branch lengths stay biologically meaningful.

**Net effect:** low breadth reduces the **amount** of evidence (wider confidence, i.e. it could lower bootstrap support), but it does **not** systematically push the picareiro samples artificially far from or close to others. The danger from low breadth is *under-powered* support, **not** *biased* topology.

### Q5. "What is the residual risk, then?"
The honest caveats:
- **Private SNPs** (a position covered by only one sample) add to that sample's **branch length** but not to topology. With 30% breadth these are a minority but non-zero.
- **Shared SNPs** (≥2 samples covered) are what actually resolve the tree. The more breadth overlap between samples, the higher the support.
- So the right question is **"how many SNP columns are shared vs. private?"** — quantified next.

---

## 4. How to quantify the overlap (run on the cluster)

**"Threshold" = the minimum number of samples that must have a real base (non-`N`) in a column for that column to survive.** The rule "remove any position where **more than 1 sample is `N`**" is the threshold *"at most 1 `N` allowed"* — i.e. keep a column only if at least `(N_samples − 1)` samples have real data. The script below prints the surviving-site count for that exact rule (highlighted) plus the full range of thresholds for comparison:

```bash
python3 - <<'EOF'
from collections import Counter
fa = "18_variant_calling_0.3/GCF_039521565/GCF_039521565.snpAlignment.fasta"

aln, name = {}, None
for line in open(fa):
    line = line.rstrip()
    if line.startswith(">"): name = line[1:]; aln[name] = []
    else: aln[name].extend(line)

taxa    = list(aln)
samples = [t for t in taxa if not t.startswith("REFERENCE_")]   # exclude ref from occupancy
n       = len(aln[taxa[0]])
ntax    = len(samples)

def col_ok(i):                      # non-N SAMPLES at column i (reference excluded)
    return sum(aln[t][i].upper() not in ("N","-") for t in samples)

def informative(i):                 # >=2 samples share each of >=2 alleles (removes singletons)
    c = Counter(aln[t][i].upper() for t in samples if aln[t][i].upper() in "ACGT")
    return sum(v >= 2 for v in c.values()) >= 2

occ = [col_ok(i) for i in range(n)]
print(f"{n} SNP columns | {ntax} samples (+reference)\n")

# ---- YOUR RULE: at most 1 sample may be N (i.e. >= ntax-1 samples have data) ----
strict = sum(o >= ntax - 1 for o in occ)
print(f">>> YOUR RULE  (<=1 sample N  ==  >= {ntax-1}/{ntax} samples with data): "
      f"{strict} sites  ({100*strict/n:.1f}%)\n")

print("full threshold sweep  (min samples with data -> columns kept):")
for k in [1, 2, 3, 4, max(1,ntax//2), max(1,int(0.75*ntax)), ntax-1, ntax]:
    kept = sum(o >= k for o in occ)
    tag  = "   <-- <=1 N (your rule)" if k == ntax-1 else ("   <-- 0 N (all samples)" if k == ntax else "")
    print(f"  >= {k:>3} samples ({100*k/ntax:5.1f}% occ) : {kept:>8}  ({100*kept/n:4.1f}%){tag}")

pinf = sum(informative(i) for i in range(n))
print(f"\nparsimony-informative columns (no singletons): {pinf}  ({100*pinf/n:.1f}%)")
EOF
```

**Reading the result:**
- The `>>> YOUR RULE` line is the direct answer to "how many sites remain if I drop every column with more than one `N`."
- The **full threshold sweep** shows how the site count grows as you relax the rule (allow more `N`), so you can see the trade-off between completeness and number of sites.
- The **parsimony-informative** line is the strictest "remove private SNPs" count: only columns whose variant is shared by ≥2 samples, i.e. the subset that actually resolves topology.
- Note: with uneven breadth, "≤1 `N`" can be very aggressive — if it leaves too few sites, relax toward `≥ 50%` occupancy and report both.

---

## 5. Bottom line

The script produces a **reference-anchored SNP alignment**: every sample is genotyped across the whole genome, uncovered or unreliable positions become `N` (missing), and only positions that vary between samples are kept for tree-building. Low breadth (picareiro ~30%) means **less data**, not **biased data** — the phylogenetic model ignores missing positions entirely and only compares samples where they both have real bases, with an ascertainment-bias correction keeping branch lengths honest. The covered region is consistent across samples because identical masking and mapping rules are applied to all, and the informative SNP columns are by definition the ones with overlapping coverage. The remaining limitation is statistical power (reflected in bootstrap support), which §4 lets you measure directly.
