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
from lib_svg import Canvas, text_width
from lib_palette import (AXIS, FONT, GRID, INK, INK_MUTED, STRIP_FILL, UNION,
                         UNION_LABEL, TEXT_BOOST, assign)

BOX_HEAD = ["Post-mortem", "brain"]
BOX_LINES = ["6 LBD · 6 control", "100 ng", "SP3 clean-up"]

# font sizes before font_scale: titles, labels, small type
TITLE, LABEL, SMALL = 10.5, 9.6, 8.8

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

    W, fs, m = 1010.0, 1.9, 12.0
    k_type = fs * TEXT_BOOST                      # drawn units per font unit

    def tw(s, size, bold=False):
        """Estimated drawn width of `s` at font size `size`."""
        return text_width(s, size * k_type, bold)

    def mid(size):
        """Baseline offset that centres capitals on a line."""
        return 0.35 * size * k_type

    # band 1: sample -> three digests -> shared stages -> union
    stage_top, stage_h = 12.0, 60.0
    chip_h, pitch = 58.0, 66.0
    lanes = [stage_top + stage_h + 10.0 + chip_h / 2 + i * pitch for i in range(3)]
    mid_y = lanes[1]
    subs = [(nm, "after " + ", ".join(res)) for nm, res in SPECIFICITY]
    sx0 = m
    sx1 = sx0 + max([tw(t, TITLE, True) for t in BOX_HEAD]
                    + [tw(t, SMALL) for t in BOX_LINES]) + 24.0
    px0 = sx1 + 24.0
    px1 = px0 + max(max(tw(nm, 11, True), tw(sub, SMALL)) for nm, sub in subs) + 20.0
    stage_x = [px1 + 10.0]                        # dia-PASEF | directDIA | end
    stages = [("dia-PASEF", args.acquisition), ("directDIA", "Spectronaut · 1% FDR")]
    for title, sub in stages:
        stage_x.append(stage_x[-1] + max(tw(title, TITLE, True), tw(sub, SMALL)) + 28.0)
    ux0, ux1 = stage_x[2] + 24.0, W - m

    # band 2 sits in a translated group; y below is relative to its top
    band2 = lanes[2] + chip_h / 2 + 14.0
    ruler_y, row_h, bar_h = 38.0, 27.0, 15.0
    top = ruler_y + 6.0
    H = band2 + top + 4 * row_h + 14.0
    c = Canvas(W, H, FONT, font_scale=fs, out_w=1215.0)

    if args.letter:
        c.text(22, 42, args.letter, 11.4, INK, "start", "600")

    lh = [TITLE * k_type * 1.12] * len(BOX_HEAD) + [SMALL * k_type * 1.2] * len(BOX_LINES)
    bh = sum(lh) + 22.0
    c.rect(sx0, mid_y - bh / 2, sx1 - sx0, bh, fill=STRIP_FILL, stroke=AXIS, sw=1.1, rx=7)
    scx = (sx0 + sx1) / 2
    y = mid_y - bh / 2 + 11.0
    for k, line in enumerate(BOX_HEAD + BOX_LINES):
        y += lh[k]
        bold = k < len(BOX_HEAD)
        c.text(scx, y - 0.25 * lh[k], line, TITLE if bold else SMALL,
               INK if bold else INK_MUTED, "middle", "600" if bold else None)

    # shared stages, labelled once above the lanes
    c.rect(stage_x[0], stage_top, stage_x[2] - stage_x[0], stage_h, fill=STRIP_FILL,
           stroke=AXIS, sw=1.0, rx=4)
    c.line(stage_x[1], stage_top, stage_x[1], stage_top + stage_h, stroke=AXIS, sw=1.0)
    for (a, b), (title, sub) in zip((stage_x[:2], stage_x[1:]), stages):
        c.text((a + b) / 2, stage_top + 25, title, TITLE, INK, "middle", "600")
        c.text((a + b) / 2, stage_top + 50, sub, SMALL, INK_MUTED, "middle")

    for i, (name, sub) in enumerate(subs):
        y = lanes[i]
        col = colour[name]
        c.arrow(sx1, mid_y, px0, y, xmid=(sx1 + px0) / 2, sw=1.8, head=6.0)
        c.chip(px0, y - chip_h / 2, px1 - px0, chip_h, name, sub, colour=col)
        c.line(px1, y, stage_x[2], y, stroke=col, sw=2.8)
        for sx in stage_x:
            c.line(sx, y - 6, sx, y + 6, stroke=col, sw=2.8)
        c.arrow(stage_x[2], y, ux0, mid_y, xmid=(stage_x[2] + ux0) / 2, sw=1.8, head=6.0)

    c.rect(ux0, mid_y - chip_h / 2, ux1 - ux0, chip_h, fill=UNION, fo=0.10, stroke=UNION,
           sw=1.6, rx=7)
    c.text((ux0 + ux1) / 2, mid_y + mid(TITLE), "All proteases", TITLE, INK, "middle",
           "600")

    # band 2: the protein's observed coverage by protease, and pooled
    c.add(f'<g transform="translate(0 {band2:g})">')
    c.line(m, 0, W - m, 0, stroke=GRID, sw=1)

    head = f"{args.label} · {seq_name}" if args.label else seq_name
    tick = 10 if n <= 100 else (20 if n <= 260 else (50 if n <= 600 else 100))
    lanes2 = [(name, colour[name], masks[name]) for name in names]
    lanes2.append((UNION_LABEL, UNION, depth))
    bx1 = W - m - tw("observed", SMALL) - 20.0
    # the first ruler label clears the header
    p = min(tick / n, 0.5)
    clear = m + tw(head, LABEL, True) + 24.0 + tw(str(tick), SMALL) / 2
    bx0 = max(m + max(tw(nm, LABEL) for nm, _c, _m in lanes2) + 34.0,
              (clear - bx1 * p) / (1 - p))
    step = (bx1 - bx0) / n

    def X(res):
        return bx0 + res * step

    c.text(m, ruler_y - 8, head, LABEL, INK, "start", "600")
    for r in range(tick, n + 1, tick):
        c.line(X(r), ruler_y - 5, X(r), ruler_y - 1, stroke=AXIS, sw=1)
        c.text(X(r), ruler_y - 8, str(r), SMALL, INK_MUTED, "middle")
    c.text(bx1 + 12, ruler_y - 8, "observed", SMALL, INK_MUTED, "start")

    for i, (name, col, mask) in enumerate(lanes2):
        cy = top + (i + 0.5) * row_h
        c.rect(bx0 - 18, cy - 6, 12, 12, fill=col, rx=2)
        c.text(bx0 - 24, cy + mid(LABEL), name, LABEL, INK, "end")
        for s, e in runs(mask):
            c.rect(X(s), cy - bar_h / 2, max(X(e) - X(s), 1.0), bar_h, fill=col, fo=0.85,
                   rx=2)
        c.text(bx1 + 12, cy + mid(LABEL), f"{sum(1 for v in mask if v) / n:.0%}", LABEL,
               INK, "start", "600" if name == UNION_LABEL else None)

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
