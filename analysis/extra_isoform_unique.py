#!/usr/bin/env python3
"""Isoforms reached: what each protease contributes to the discriminating set.

    # one pass over the reports -> the table
    python3 extra_isoform_unique.py scan 'data/search/*-60min-Phospho.parquet' \\
        --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta --tsv derived/da_iso/isoform_unique_fraction.tsv

    # the inset that sits in Figure 2a's second row
    python3 extra_isoform_unique.py --letter "" --panel-w 300 --panel-h 280 \\
        --out figures/fig2b_unique.svg
"""

import argparse
import collections
import csv
import glob
import itertools
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_svg import Canvas                                  # noqa: E402
from lib_palette import (AXIS, FONT, GRID, INK, INK_MUTED, INK_SECONDARY,
                     assign, display, ink_on)                      # noqa: E402

ORDER = ["GluC", "LysC", "Trypsin"]
FIELDS = ["base", "isoform", "digest", "unique", "shared"]


def scan(paths, fasta, out, precursor_q=0.01):
    """One row per (isoform, digest): unique and shared peptide counts."""
    import fig2b_isoform_strip as ie
    from lib_fasta import read_fasta
    seqs = read_fasta(fasta)
    groups = ie.collect_groups(sorted(paths), precursor_q)

    bybase = collections.defaultdict(lambda: collections.defaultdict(set))
    for g, byd in groups.items():
        accs = [a.strip() for a in g.split(";") if a.strip()]
        bases = {ie.ISO.sub("", a) for a in accs}
        if len(bases) != 1:
            continue
        base = bases.pop()
        for d, v in byd.items():
            bybase[base][d] |= v

    isoforms = collections.defaultdict(list)
    for a in seqs:
        if ie.ISO.search(a):
            isoforms[ie.ISO.sub("", a)].append(a)

    rows = []
    for base, byd in bybase.items():
        if base not in seqs:
            continue
        cseq = seqs[base]
        for iso in isoforms.get(base, ()):
            iseq = seqs.get(iso)
            if not iseq or iseq == cseq:
                continue
            for d, peps in byd.items():
                uniq = shared = 0
                for p in peps:
                    in_i, in_c = p in iseq, p in cseq
                    if in_i and not in_c:
                        uniq += 1
                    elif in_i and in_c:
                        shared += 1
                if uniq or shared:
                    rows.append({"base": base, "isoform": iso, "digest": d,
                                 "unique": uniq, "shared": shared})
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {out}  ({len(rows):,} isoform x digest rows, "
          f"{len({r['isoform'] for r in rows}):,} isoforms)")
    return rows


def load(path):
    """-> {isoform: {digest: [unique, shared]}}."""
    per = collections.defaultdict(lambda: collections.defaultdict(
        lambda: [0, 0]))
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            a = per[r["isoform"]][r["digest"]]
            a[0] += int(r["unique"])
            a[1] += int(r["shared"])
    return per


def accrual(per, order=None, min_unique=1):
    """-> {subset: isoforms reached}, evidence pooled over each subset of proteases."""
    order = order or ORDER
    isos = set(per)
    out = {}
    for n in range(1, len(order) + 1):
        for sub in itertools.combinations(order, n):
            out[sub] = sum(
                1 for i in isos
                if sum(per[i][d][0] for d in sub if d in per[i]) >= min_unique)
    return out


def shares(per, order=None):
    """-> {digest|'union': [per-isoform discriminating share]}."""
    order = order or ORDER
    isos = set(per)
    out = {}
    for d in order:
        out[d] = [per[i][d][0] / (per[i][d][0] + per[i][d][1]) for i in isos
                  if d in per[i] and per[i][d][0] >= 1]
    tot = [(sum(per[i][d][0] for d in per[i]),
            sum(per[i][d][1] for d in per[i])) for i in isos]
    out["union"] = [u / (u + s) for u, s in tot if u >= 1]
    return out


def panel(per, out, font, letter="b", width=1215.0, ts=1.7, order=None,
          min_unique=2, panel_w=None, panel_h=None):
    """A pie of each protease's incremental contribution to the isoforms reached."""
    order = order or ORDER
    reach = accrual(per, order, min_unique)
    colour = assign(list(order) + ["All"])

    seq = sorted(order, key=lambda d: -reach[(d,)])
    segs, prev = [], 0
    for n in range(1, len(seq) + 1):
        sub = tuple(d for d in order if d in seq[:n])
        segs.append((seq[n - 1], reach[sub] - prev))
        prev = reach[sub]
    total = prev

    compact = panel_w is not None and panel_w < 500
    W = panel_w or 1010.0
    ml = mr = 14.0 if compact else 40.0
    mt = 52.0 if compact else 74.0
    pw = W - ml - mr
    H = panel_h or 314.0
    ph = H - mt - (34 if compact else 44)

    c = Canvas(W, H, font, font_scale=1.45 * ts / 1.7, out_w=width)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    if compact:
        c.text(ml + pw / 2, 26, f"{total:,} isoforms reached", 10.4, INK,
               "middle", "600")
        c.text(ml + pw / 2, 40, f"≥ {min_unique} discriminating peptides", 8.2,
               INK_MUTED, "middle")
    else:
        c.text(ml, 30, f"{total:,} isoforms reached — each protease reaches "
               f"some the others cannot", 11.5, INK, "start", "600")
        c.text(ml, 48, f"≥ {min_unique} peptides absent from their canonical "
               f"sequence", 8.8, INK_MUTED, "start")

    # Label radius grows until the slice's arc can hold the number.
    r = min(58.0, ph * 0.30)
    cx, cy = ml + pw / 2, mt + ph * 0.50
    fsz, fs = 9.0, 1.45 * ts / 1.7
    ang = -math.pi / 2
    for dig, v in segs:
        sweep = 2 * math.pi * v / total
        x1, y1 = cx + r * math.cos(ang), cy + r * math.sin(ang)
        x2 = cx + r * math.cos(ang + sweep)
        y2 = cy + r * math.sin(ang + sweep)
        big = 1 if sweep > math.pi else 0
        c.add(f'<path d="M {cx:.1f} {cy:.1f} L {x1:.1f} {y1:.1f} '
              f'A {r:.1f} {r:.1f} 0 {big} 1 {x2:.1f} {y2:.1f} Z" '
              f'fill="{colour[dig]}" fill-opacity="0.9" stroke="#ffffff" '
              f'stroke-width="1.4"/>')
        lab = f"{v:,}"
        need = (0.58 * fsz * fs * len(lab) + 6.0) / max(sweep, 1e-9)
        rl = min(max(need, 0.45 * r), 0.80 * r)
        mid = ang + sweep / 2
        c.text(cx + rl * math.cos(mid), cy + rl * math.sin(mid) + 3.2, lab,
               fsz, ink_on(colour[dig]), "middle", "600")
        ang += sweep

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({W:.0f}x{H:.0f})  threshold >= {min_unique}")
    for dig, v in segs:
        print(f"    {display(dig):<10s} contributes {v:>6,}")
    best1 = max(v for k, v in reach.items() if len(k) == 1)
    print(f"    total {total:,};  best single -> all three: {best1:,} -> "
          f"{total:,} (+{100 * (total - best1) / best1:.0f}%)")
    for k in (1, 2, 3):
        rr = accrual(per, order, k)
        b = max(v for kk, v in rr.items() if len(kk) == 1)
        print(f"    at >= {k}: {b:,} -> {rr[tuple(order)]:,}  "
              f"(+{100 * (rr[tuple(order)] - b) / b:.0f}%)")
    sh = shares(per, order)
    print("    median discriminating share: " + ", ".join(
        f"{display(k) if k in order else k} {statistics.median(v):.1%}"
        for k, v in sh.items()))
    print("    all subsets: " + ", ".join(
        f"{'+'.join(display(d) for d in k)} {v:,}"
        for k, v in sorted(reach.items(), key=lambda kv: (len(kv[0]), kv[1]))))
    return segs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", nargs="?", default="panel",
                    choices=("panel", "scan"))
    ap.add_argument("reports", nargs="*", default=["data/search/*-60min-Phospho.parquet"])
    ap.add_argument("--fasta", default="data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta")
    ap.add_argument("--tsv", default="derived/da_iso/isoform_unique_fraction.tsv")
    ap.add_argument("--precursor-q", type=float, default=0.01)
    ap.add_argument("--min-unique", type=int, default=2,
                    help="discriminating peptides needed, pooled over digests")
    ap.add_argument("--panel-w", type=float, default=None)
    ap.add_argument("--panel-h", type=float, default=None)
    ap.add_argument("--letter", default="b")
    ap.add_argument("--out", default="figures/fig2b_unique.svg")
    ap.add_argument("--font", default=FONT)
    ap.add_argument("--width", type=float, default=1215.0)
    args = ap.parse_args(argv)

    if args.mode == "scan":
        paths = sorted(p for pat in args.reports for p in glob.glob(pat))
        if not paths:
            sys.exit("no reports matched")
        scan(paths, args.fasta, args.tsv, args.precursor_q)
        return 0
    panel(load(args.tsv), args.out, args.font, args.letter, args.width,
          min_unique=args.min_unique, panel_w=args.panel_w,
          panel_h=args.panel_h)
    return 0


if __name__ == "__main__":
    sys.exit(main())
