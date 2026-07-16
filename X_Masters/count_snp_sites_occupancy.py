#!/usr/bin/env python3
"""
count_snp_sites_occupancy.py — how many SNP sites remain after occupancy filtering?

Answers: "If I remove every column where MORE THAN 1 sample is N, how many sites
remain for phylogeny?"  (your cutoff = at most 1 sample may be N per column,
i.e. a column is kept only if >= N_samples - 1 samples carry a real A/C/G/T base.)

The reference tip (>REFERENCE_...) is EXCLUDED from the occupancy count: it is a
complete genome with no N, so it would otherwise make every column look 1 deeper.
It is still a valid tree tip; we only measure how many SAMPLES have real data.

USAGE:
  # one alignment
  python3 count_snp_sites_occupancy.py 18_variant_calling_0.3/GCF_039521565/GCF_039521565.snpAlignment.fasta

  # every genome under one results dir (auto-discovers *.snpAlignment.fasta)
  python3 count_snp_sites_occupancy.py 18_variant_calling_0.3

  # several dirs / files at once, and also WRITE the filtered alignments
  python3 count_snp_sites_occupancy.py --write 18_variant_calling_0.3 18_variant_calling_0.3_2

OPTIONS:
  --max-n K     keep a column only if AT MOST K samples are N  (default 1 = your cutoff)
  --write       also write <name>.snpAlignment.maxN<K>.fasta next to each input
  --quiet       print only the final summary table (no per-genome detail block)
  --no-per-sample  suppress the per-sample non-N site table (still shown by default)
"""

import argparse
import glob
import os
import sys
from collections import Counter


def read_fasta(path):
    """Return dict name -> list(chars). Streams line by line (alignments can be large)."""
    aln, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                aln[name] = []
            elif name is not None:
                aln[name].extend(line)
    return aln


def analyse(path, max_n, write, quiet, per_sample, write_pi=False):
    aln = read_fasta(path)
    if not aln:
        print(f"[skip] {path}: no sequences")
        return None

    taxa = list(aln)
    samples = [t for t in taxa if not t.startswith("REFERENCE_")]
    if not samples:
        print(f"[skip] {path}: no non-reference samples")
        return None

    n = len(aln[taxa[0]])
    ntax = len(samples)
    if n == 0:
        print(f"[skip] {path}: empty alignment (0 columns)")
        return None

    # occupancy = number of SAMPLES with a real base (non-N/-) at each column
    occ = [
        sum(aln[t][i].upper() not in ("N", "-") for t in samples)
        for i in range(n)
    ]

    # your cutoff: at most `max_n` samples may be N  <=>  occupancy >= ntax - max_n
    keep_thresh = max(1, ntax - max_n)
    keep_idx = [i for i in range(n) if occ[i] >= keep_thresh]
    kept = len(keep_idx)

    # parsimony-informative among samples (variant shared by >=2 samples; no singletons)
    def informative(i):
        c = Counter(
            aln[t][i].upper() for t in samples if aln[t][i].upper() in "ACGT"
        )
        return sum(v >= 2 for v in c.values()) >= 2

    pi_idx = [i for i in range(n) if informative(i)]
    pinf = len(pi_idx)

    genome = os.path.basename(path).replace(".snpAlignment.fasta", "")

    if not quiet:
        print(f"=== {genome} ===")
        print(f"  file           : {path}")
        print(f"  samples        : {ntax}  (+1 reference tip)")
        print(f"  total SNP cols : {n}")
        print(f"  CUTOFF (<= {max_n} sample N  ==  >= {keep_thresh}/{ntax} samples with data):")
        print(f"      sites remaining : {kept}  ({100*kept/n:.1f}% of {n})")
        print(f"  parsimony-informative (>=2 samples share variant): {pinf}  ({100*pinf/n:.1f}%)")
        # small sweep so you can see the trade-off of relaxing the rule
        print("  sweep (min samples with data -> sites):")
        for k in sorted({1, 2, 3, max(1, ntax // 2),
                         max(1, int(0.75 * ntax)), ntax - 1, ntax}):
            s = sum(o >= k for o in occ)
            tag = "  <-- your cutoff" if k == keep_thresh else ""
            print(f"      >= {k:>3}/{ntax} ({100*k/ntax:5.1f}%): {s:>8}  ({100*s/n:4.1f}%){tag}")
        print()

    # --- Per-sample non-N site counts (over the KEPT columns only) ---
    if per_sample and not quiet:
        print(f"  Per-sample non-N sites (over {kept} retained columns):")
        print(f"    {'Sample':<40} {'Non-N':>8}  {'Non-N%':>7}")
        ref_names = [t for t in taxa if t.startswith("REFERENCE_")]
        sample_names = [t for t in taxa if not t.startswith("REFERENCE_")]
        for t in ref_names + sorted(sample_names):
            nn = sum(aln[t][i].upper() in "ACGT" for i in keep_idx) if kept > 0 else 0
            pct = 100 * nn / kept if kept > 0 else 0.0
            tag = "  (reference)" if t.startswith("REFERENCE_") else ""
            print(f"    {t:<40} {nn:>8}  {pct:>6.1f}%{tag}")
        print()

    if write and kept > 0:
        out = path.replace(".snpAlignment.fasta", f".snpAlignment.maxN{max_n}.fasta")
        with open(out, "w") as fh:
            for t in taxa:  # keep reference + all samples as rows
                fh.write(f">{t}\n")
                fh.write("".join(aln[t][i] for i in keep_idx) + "\n")
        print(f"  [written] {out}  ({kept} columns x {len(taxa)} taxa)")

    if write_pi and pinf > 0:
        out = path.replace(".snpAlignment.fasta", ".snpAlignment.PI.fasta")
        with open(out, "w") as fh:
            for t in taxa:  # keep reference + all samples as rows
                fh.write(f">{t}\n")
                fh.write("".join(aln[t][i] for i in pi_idx) + "\n")
        print(f"  [written] {out}  ({pinf} parsimony-informative columns x {len(taxa)} taxa)")

    return {"genome": genome, "samples": ntax, "total": n,
            "kept": kept, "pinf": pinf, "thresh": keep_thresh}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+",
                    help="snpAlignment.fasta file(s) and/or results dir(s) to scan")
    ap.add_argument("--max-n", type=int, default=1,
                    help="max samples allowed to be N per column (default 1 = your cutoff)")
    ap.add_argument("--write", action="store_true",
                    help="also write filtered <name>.snpAlignment.maxN<K>.fasta")
    ap.add_argument("--write-pi", dest="write_pi", action="store_true",
                    help="also write parsimony-informative-only <name>.snpAlignment.PI.fasta")
    ap.add_argument("--quiet", action="store_true",
                    help="print only the final summary table")
    ap.add_argument("--no-per-sample", dest="per_sample", action="store_false",
                    help="suppress the per-sample non-N site table")
    ap.set_defaults(per_sample=True)
    args = ap.parse_args()

    # expand inputs -> concrete FASTA paths
    fastas = []
    for item in args.inputs:
        if os.path.isdir(item):
            fastas += sorted(glob.glob(os.path.join(item, "**", "*.snpAlignment.fasta"),
                                       recursive=True))
        elif os.path.isfile(item):
            fastas.append(item)
        else:
            print(f"[warn] not found: {item}", file=sys.stderr)
    if not fastas:
        sys.exit("No *.snpAlignment.fasta inputs found.")

    rows = [r for f in fastas if (r := analyse(f, args.max_n, args.write, args.quiet, args.per_sample, args.write_pi))]

    # summary table
    print("=" * 78)
    print(f"SUMMARY  (cutoff: <= {args.max_n} sample N per column)")
    print("=" * 78)
    print(f"{'Genome':<22}{'Samples':>8}{'TotalSNP':>10}{'Kept':>10}{'Kept%':>8}{'ParsInf':>9}")
    for r in rows:
        pct = 100 * r["kept"] / r["total"] if r["total"] else 0
        print(f"{r['genome']:<22}{r['samples']:>8}{r['total']:>10}"
              f"{r['kept']:>10}{pct:>7.1f}%{r['pinf']:>9}")


if __name__ == "__main__":
    main()
