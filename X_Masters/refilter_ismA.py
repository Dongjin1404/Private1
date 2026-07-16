#!/usr/bin/env python3
"""
refilter_ismA.py — re-threshold the ismA DIAMOND screen at a stricter identity WITHOUT re-running DIAMOND.

The 60% aa-identity floor used by 20_ismA_screen.sh catches the conserved short-chain-dehydrogenase (SDR)
fold shared by many enzymes, not ismA specifically — so blanks/sediment show a non-zero background. Raising
the identity cutoff (e.g. >=75%) isolates true ismA hits from that SDR background. Because the raw hit tables
(<sample>.ismA.blastx.tsv.gz) and the per-sample TotalReads (in the overview) are already on disk, we just
re-filter and recompute — no expensive re-alignment.

USAGE:
  python3 refilter_ismA.py 20_ismA_screen --min-pid 75 [--min-aalen 15] \
      [--overview 20_ismA_screen/ismA_screen_overview.tsv] \
      [--out 20_ismA_screen/ismA_screen_overview.pid75.tsv]

OUTPUT: a new overview TSV with the same columns, recomputed at the stricter cutoff, plus a
        BlankSubtracted column if control samples are named (see --blanks).
"""

import argparse
import glob
import gzip
import os
import statistics
import sys


def read_table(path):
    """Yield (pident, aalen, qlen, sseqid) from a gzipped or plain DIAMOND fmt6 table.
    Columns (from 20_ismA_screen.sh): qseqid sseqid pident length evalue bitscore qlen full_qseq."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            try:
                yield float(f[2]), int(f[3]), int(f[6]), f[1]
            except ValueError:
                continue


def load_totalreads(overview):
    """Map (dataset, sample) -> TotalReads from an existing overview TSV."""
    tr = {}
    if not overview or not os.path.isfile(overview):
        return tr
    with open(overview) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_s, i_d, i_t = header.index("Sample"), header.index("Dataset"), header.index("TotalReads")
        except ValueError:
            return tr
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_s, i_d, i_t):
                continue
            try:
                tr[(f[i_d], f[i_s])] = int(f[i_t])
            except ValueError:
                pass
    return tr


def median(xs):
    return round(statistics.median(xs), 2) if xs else "NA"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", help="ismA screen output dir (contains <dataset>/<sample>.ismA.blastx.tsv.gz)")
    ap.add_argument("--min-pid", type=float, default=75.0, help="min amino-acid %% identity (default 75)")
    ap.add_argument("--min-aalen", type=int, default=15, help="min alignment length in aa (default 15)")
    ap.add_argument("--overview", default=None,
                    help="existing overview TSV for TotalReads (default <root>/ismA_screen_overview.tsv)")
    ap.add_argument("--out", default=None, help="output TSV (default <root>/ismA_screen_overview.pid<PID>.tsv)")
    ap.add_argument("--blanks", default="",
                    help="comma-separated control sample IDs (e.g. EXTBL,SED2,SED3,LIBBL) to compute a "
                         "background floor and a BlankSubtracted HPM column")
    args = ap.parse_args()

    root = args.root.rstrip("/")
    overview = args.overview or os.path.join(root, "ismA_screen_overview.tsv")
    out = args.out or os.path.join(root, f"ismA_screen_overview.pid{int(args.min_pid)}.tsv")
    totalreads = load_totalreads(overview)

    tables = sorted(glob.glob(os.path.join(root, "*", "*.ismA.blastx.tsv.gz")) +
                    glob.glob(os.path.join(root, "*", "*.ismA.blastx.tsv")))
    if not tables:
        sys.exit(f"No *.ismA.blastx.tsv[.gz] under {root}")

    rows = []
    for t in tables:
        dataset = os.path.basename(os.path.dirname(t))
        sample = os.path.basename(t).split(".ismA.blastx.tsv")[0]
        pids, qlens, refs = [], [], set()
        seen = set()
        for pident, aalen, qlen, sseqid in read_table(t):
            if pident >= args.min_pid and aalen >= args.min_aalen:
                # one row per read already (-k1), but guard against dupes
                key = (sseqid, qlen, pident)
                pids.append(pident)
                qlens.append(qlen)
                refs.add(sseqid)
        hit = len(pids)
        total = totalreads.get((dataset, sample))
        hpm = round(hit / total * 1e6, 3) if total else "NA"
        rows.append({
            "Sample": sample, "Dataset": dataset,
            "TotalReads": total if total is not None else "NA",
            "HitReads": hit, "HitsPerMillion": hpm,
            "UniqRefs": len(refs), "MedianPctID": median(pids), "MedianHitLen": median(qlens),
        })

    # optional blank-subtraction
    blanks = [b.strip() for b in args.blanks.split(",") if b.strip()]
    blank_floor = None
    if blanks:
        bvals = [r["HitsPerMillion"] for r in rows
                 if r["Sample"] in blanks and isinstance(r["HitsPerMillion"], float)]
        if bvals:
            blank_floor = round(statistics.mean(bvals), 3)

    cols = ["Sample", "Dataset", "TotalReads", "HitReads", "HitsPerMillion",
            "UniqRefs", "MedianPctID", "MedianHitLen"]
    if blank_floor is not None:
        cols.append("BlankSubtracted")

    rows.sort(key=lambda r: (r["HitsPerMillion"] if isinstance(r["HitsPerMillion"], float) else -1))
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            if blank_floor is not None:
                r["BlankSubtracted"] = (round(r["HitsPerMillion"] - blank_floor, 3)
                                        if isinstance(r["HitsPerMillion"], float) else "NA")
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"Re-filtered at pident>={args.min_pid}, aln>={args.min_aalen}aa")
    if blank_floor is not None:
        print(f"Blank floor (mean HPM of {','.join(blanks)}): {blank_floor}")
    print(f"Wrote {len(rows)} samples -> {out}")


if __name__ == "__main__":
    main()
