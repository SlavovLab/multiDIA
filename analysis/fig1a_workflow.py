#!/usr/bin/env python3
"""Figure 1a: the multiDIA workflow and one protein's observed coverage by protease, as SVG.

    python3 analysis/fig1a_workflow.py 'data/search/*-60min-Phospho.parquet' --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta --protein Q96RQ3 --label MCCC1 --acquisition '65 min · 60 windows' --out figures/fig1a.svg
"""

import argparse
import collections
import glob
import os
import sys

from lib_fasta import read_fasta
from lib_svg import Canvas
from lib_palette import (AXIS, FONT, GRID, INK, INK_MUTED, STRIP_FILL, UNION,
                         UNION_LABEL, TEXT_BOOST, assign)

BOX_LINES = [("Post-mortem brain", 10.5), ("6 LBD · 6 control", 9.2),
             ("100 ng", 9.2), ("SP3 clean-up", 9.2)]

# (display name, residues cleaved C-terminal to)
SPECIFICITY = [("Glu-C", "E"), ("Lys-C", "K"), ("Trypsin", "KR")]


def coverage_mask(peptides, n):
    mask = bytearray(n)
    for s, e in peptides:
        for i in range(s, e):
            mask[i] = 1
    return mask


def observed(paths, fasta, acc):
    """-> {protease: [(start, end)]} for the reports' peptides found only in `acc`'s gene."""
    import lib_report as rp
    from lib_fasta import gene_map
    from lib_palette import display
    from prep_counts import PRECURSOR_Q
    seqs = read_fasta(fasta)
    gm = gene_map(fasta)
    seq = seqs[acc]
    others = [s for a, s in seqs.items() if gm.get(a) != gm.get(acc)]
    found = collections.defaultdict(set)
    for r in rp.open_reports(paths):
        for b in r.batches(["peptide", "precursor_q"]):
            for k in range(b["_n"]):
                p, q = b["peptide"][k], b["precursor_q"][k]
                if p and q is not None and q <= PRECURSOR_Q and p in seq:
                    found[display(r.protease)].add(p)
    out = {}
    for name, peps in found.items():
        spans = []
        for p in peps:
            if any(p in s for s in others):
                continue
            i = seq.find(p)
            while i != -1:
                spans.append((i, i + len(p)))
                i = seq.find(p, i + 1)
        out[name] = sorted(spans)
    return out


def runs(mask):
    """-> [(start, end)] stretches where `mask` is set."""
    out, i, n = [], 0, len(mask)
    while i < n:
        while i < n and not mask[i]:
            i += 1
        j = i
        while j < n and mask[j]:
            j += 1
        if j > i:
            out.append((i, j))
        i = j
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reports", nargs="+", help="the search reports, one per protease")
    ap.add_argument("--fasta", required=True, help="FASTA holding --protein")
    ap.add_argument("--protein", required=True, help="accession whose coverage is drawn")
    ap.add_argument("--label", default=None, help="gene name shown with the accession")
    ap.add_argument("--acquisition", required=True, help="sub-label under dia-PASEF")
    ap.add_argument("--letter", default="a", help="panel letter; '' for none")
    ap.add_argument("--out", default="fig1a.svg")
    args = ap.parse_args(argv)

    seqs = read_fasta(args.fasta)
    if args.protein not in seqs:
        sys.exit(f"{args.protein} not in {args.fasta}")
    seq_name, seq = args.protein, seqs[args.protein]
    paths = sorted(p for pat in args.reports for p in glob.glob(pat))
    if not paths:
        sys.exit("no reports matched")
    seen = observed(paths, args.fasta, args.protein)

    names = [nm for nm, _ in SPECIFICITY]
    colour = assign(names + ["All"])
    n = len(seq)
    peptides = {name: seen.get(name, []) for name in names}
    masks = {name: coverage_mask(peptides[name], n) for name in names}
    depth = [sum(masks[k][i] for k in masks) for i in range(n)]

    W, H = 1010, 360
    band2 = 30.0            # offset of the bottom band
    c = Canvas(W, H, FONT, font_scale=1.25, out_w=1215.0)

    if args.letter:
        c.text(22, 32, args.letter, 13, INK, "start", "600")

    # band 1: sample -> three digests -> shared stages -> union
    lanes = [92.0, 136.0, 180.0]
    mid_y = lanes[1]
    # sample box sized to its widest line
    k_type = c.fs * TEXT_BOOST                    # drawn units per font unit
    need = max(len(t) * sz * 0.53 * k_type for t, sz in BOX_LINES)
    sx0 = 40.0
    sx1 = sx0 + max(112.0, need + 16.0)
    c.rect(sx0, 76, sx1 - sx0, 120, fill=STRIP_FILL, stroke=AXIS, sw=1.1, rx=7)
    scx = (sx0 + sx1) / 2
    head, head_sz = BOX_LINES[0]
    c.text(scx, 104, head, head_sz, INK, "middle", "600")
    for k, (line, sz) in enumerate(BOX_LINES[1:]):
        c.text(scx, 128 + k * sz * k_type * 1.2, line, sz, INK_MUTED, "middle")

    px0 = max(196.0, sx1 + 28.0)
    px1 = px0 + 100.0
    stage_x = [356.0, 552.0, 716.0]          # dia-PASEF | directDIA | end
    ux0, ux1 = 752.0, 980.0

    # shared stages, labelled once above the lanes
    c.rect(stage_x[0], 26, stage_x[2] - stage_x[0], 42, fill=STRIP_FILL,
           stroke=AXIS, sw=1.0, rx=4)
    c.line(stage_x[1], 26, stage_x[1], 68, stroke=AXIS, sw=1.0)
    for (a, b), title, sub in ((stage_x[:2], "dia-PASEF", args.acquisition),
                               (stage_x[1:], "directDIA", "Spectronaut · 1% FDR")):
        c.text((a + b) / 2, 45, title, 10, INK, "middle", "600")
        c.text((a + b) / 2, 62, sub, 8.4, INK_MUTED, "middle")

    for i, (name, residues) in enumerate(SPECIFICITY):
        y = lanes[i]
        col = colour[name]
        c.arrow(sx1, mid_y, px0, y, xmid=(sx1 + px0) / 2)
        c.chip(px0, y - 19, px1 - px0, 38, name,
               "after " + ", ".join(residues), colour=col)
        c.line(px1, y, stage_x[2], y, stroke=col, sw=2.2)
        for sx in stage_x:
            c.line(sx, y - 4, sx, y + 4, stroke=col, sw=2.2)
        c.arrow(stage_x[2], y, ux0, mid_y, xmid=740)

    bh = 44.0
    c.rect(ux0, mid_y - bh / 2, ux1 - ux0, bh, fill=UNION, fo=0.10, stroke=UNION,
           sw=1.6, rx=7)
    c.text((ux0 + ux1) / 2, mid_y + 4, "All proteases", 11, INK, "middle", "600")

    # band 2: the protein's observed coverage by protease, and pooled
    c.add(f'<g transform="translate(0 {band2:g})">')
    c.line(44, 196, W - 30, 196, stroke=GRID, sw=1)

    bx0, bx1 = 116.0, W - 74.0
    pitch = (bx1 - bx0) / n

    def X(res):
        return bx0 + res * pitch

    ruler_y = 224.0
    c.text(26, ruler_y - 8,
           f"{args.label} · {seq_name}" if args.label else seq_name,
           9.5, INK, "start", "600")
    tick = 10 if n <= 100 else (20 if n <= 260 else (50 if n <= 600 else 100))
    for r in range(tick, n + 1, tick):
        c.line(X(r), ruler_y - 5, X(r), ruler_y - 1, stroke=AXIS, sw=1)
        c.text(X(r), ruler_y - 8, str(r), 8.2, INK_MUTED, "middle")
    c.text(bx1 + 8, ruler_y - 8, "observed", 8.2, INK_MUTED, "start")

    block_h, gap_h = 12.0, 10.0
    top = ruler_y + 10
    lanes2 = [(name, colour[name], masks[name]) for name in names]
    lanes2.append((UNION_LABEL, UNION, depth))
    for i, (name, col, mask) in enumerate(lanes2):
        by = top + i * (block_h + gap_h)
        c.rect(bx0 - 13, by + 7, 9, 9, fill=col, rx=2)
        c.text(bx0 - 19, by + 15, name, 9.5, INK, "end")
        for s, e in runs(mask):
            c.rect(X(s), by + 6, max(X(e) - X(s), 1.0), 12, fill=col, fo=0.85, rx=2)
        c.text(bx1 + 8, by + 15, f"{sum(1 for v in mask if v) / n:.0%}", 9.5, INK, "start",
               "600" if name == UNION_LABEL else None)

    c.add("</g>")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(c.out())

    print(f"wrote {args.out}  ({W}x{H}px)")
    print(f"  sequence  {seq_name}, {n} aa")
    for name in names:
        print(f"  {name:<8s} {len(peptides[name]):>4d} peptides, "
              f"{sum(masks[name]) / n:6.1%} of residues   {colour[name]}")
    cov = sum(1 for d in depth if d)
    print(f"  {'union':<8s} {cov:>4d} residues covered ({cov / n:.1%}), "
          f"{sum(1 for d in depth if d >= 2)} by >=2 digests")
    return 0


if __name__ == "__main__":
    sys.exit(main())
