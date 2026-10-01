#!/usr/bin/env python3
"""Figure 1a: the multiDIA workflow schematic and in-silico digest, as SVG.

    python3 fig1a_workflow.py --out figures/fig1a.svg
    python3 fig1a_workflow.py --fasta P49189.fasta --window 90
"""

import argparse
import os
import random
import sys

from lib_svg import Canvas, esc                                     # noqa: E402
from lib_palette import (AXIS, DEPTH, FONT, GRID, INK, INK_MUTED, INK_SECONDARY,
                     STRIP_FILL, SURFACE, UNION, UNION_LABEL, TEXT_BOOST, assign,
                     boost_type)

# approximate human proteome amino-acid composition (%)
COMPOSITION = {
    "A": 7.0, "R": 5.6, "N": 3.6, "D": 4.7, "C": 2.3, "Q": 4.7, "E": 7.1,
    "G": 6.6, "H": 2.6, "I": 4.3, "L": 9.9, "K": 5.7, "M": 2.1, "F": 3.6,
    "P": 6.3, "S": 8.3, "T": 5.4, "W": 1.2, "Y": 2.7, "V": 6.0,
}

# (display name, residues cleaved C-terminal to, blocked by a following proline)
SPECIFICITY = [
    ("Glu-C", "E", False),
    ("Lys-C", "K", False),
    ("Trypsin", "KR", True),
]


def illustrative_sequence(length, seed):
    rng = random.Random(seed)
    letters = list(COMPOSITION)
    weights = [COMPOSITION[a] for a in letters]
    return "".join(rng.choices(letters, weights=weights, k=length))


def read_fasta(path, accession=None):
    """-> (name, sequence). Without `accession`, the first entry in the file."""
    want = (accession or "").strip()
    name, seq, cur, hit = None, [], None, False
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if hit and seq:
                    break
                head = line[1:].split()[0]
                parts = [p for p in head.split("|") if p]
                cur = parts[1] if len(parts) >= 3 else parts[0]
                hit = (not want) or cur == want or head == want
                if hit:
                    name, seq = cur, []
            elif line and hit:
                seq.append("".join(ch for ch in line.upper() if ch.isalpha()))
    if not seq:
        sys.exit(f"{path}: no sequence found"
                 + (f" for {want}" if want else ""))
    return name or os.path.basename(path), "".join(seq)


def digest(seq, residues, block_p, missed, lo, hi):
    """Fully/partially cleaved peptides as (start, end) 0-based half-open."""
    sites = [0]
    for i, aa in enumerate(seq[:-1]):
        if aa in residues and not (block_p and seq[i + 1] == "P"):
            sites.append(i + 1)
    sites.append(len(seq))
    out = []
    for i in range(len(sites) - 1):
        for j in range(i + 1, min(i + 2 + missed, len(sites))):
            start, end = sites[i], sites[j]
            if lo <= end - start <= hi:
                out.append((start, end))
    return out


def coverage_mask(peptides, n):
    mask = bytearray(n)
    for s, e in peptides:
        for i in range(s, e):
            mask[i] = 1
    return mask


def longest_run(predicate, n):
    """-> (start, end) of the longest maximal run where predicate(i) holds."""
    best = cur = None
    for i in range(n + 1):
        ok = i < n and predicate(i)
        if ok and cur is None:
            cur = i
        elif not ok and cur is not None:
            if best is None or i - cur > best[1] - best[0]:
                best = (cur, i)
            cur = None
    return best


def union_gains(counts_csv, units=("peptides", "protein_isoform_groups")):
    """-> ['262K peptides  +59%', ...] medians from the count table."""
    import csv as _csv
    import collections as _c
    vals = _c.defaultdict(lambda: _c.defaultdict(list))
    with open(counts_csv, newline="") as fh:
        for r in _csv.DictReader(fh):
            vals[r["unit"].strip()][r["protease"].strip()].append(int(r["count"]))

    def med(v):
        v = sorted(v)
        n = len(v)
        return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2

    def short(x):
        return (f"{x / 1e3:.0f}K" if x >= 1e5 else
                f"{x / 1e3:.1f}K" if x >= 1e3 else f"{x:,.0f}")

    label = {"peptides": "peptides", "precursors": "precursors",
             "protein_isoform_groups": "proteoforms",
             "proteins_canonical": "proteins"}
    out = []
    for u in units:
        if u not in vals or "All" not in vals[u] or "Trypsin" not in vals[u]:
            continue
        a, t = med(vals[u]["All"]), med(vals[u]["Trypsin"])
        out.append(f"{short(a)} {label.get(u, u)}   +{a / t - 1:.0%}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="fig1a.svg")
    ap.add_argument("--fasta", default=None, help="draw a real protein instead")
    ap.add_argument("--protein", default=None,
                    help="accession in a multi-entry --fasta")
    ap.add_argument("--label", default=None,
                    help="gene name shown with the accession")
    ap.add_argument("--length", type=int, default=300,
                    help="illustrative sequence length (default 300)")
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--min-len", type=int, default=7, help="detectable peptide floor")
    ap.add_argument("--max-len", type=int, default=30, help="detectable peptide ceiling")
    ap.add_argument("--missed-cleavages", type=int, default=0)
    ap.add_argument("--detect-prob", type=float, default=1.0,
                    help="fraction of in-silico peptides to keep")
    ap.add_argument("--samples", default="12", help="cohort size shown in the first box")
    ap.add_argument("--tissue", default="Post-mortem brain",
                    help="title of the first box")
    ap.add_argument("--groups", default="6 LBD · 6 control",
                    help="cohort composition")
    ap.add_argument("--protein-amount", default="100 ng",
                    help="protein per sample shown in the first box")
    ap.add_argument("--acquisition", default="22 min gradient · 30 min run",
                    help="sub-label under dia-PASEF")
    ap.add_argument("--prep", default="SP3 clean-up",
                    help="sample-prep line in the first box")
    ap.add_argument("--window", type=int, default=0,
                    help="residues shown; 0 draws the whole protein")
    ap.add_argument("--window-start", type=int, default=None,
                    help="first residue of the window (0-based)")
    ap.add_argument("--window-rule", choices=("cuts", "balance"), default="cuts",
                    help="most cleavage sites, or most balanced coverage")
    ap.add_argument("--min-union", type=float, default=0.90)
    ap.add_argument("--min-unique", type=float, default=0.05)
    ap.add_argument("--gluc-residues", default="E",
                    help="Glu-C specificity: 'E' or 'ED'")
    ap.add_argument("--letter", default="a", help="panel letter; '' for none")
    ap.add_argument("--font", default=FONT)
    ap.add_argument("--counts", default=None, metavar="COUNTS_CSV",
                    help="count table")
    ap.add_argument("--text-scale", type=float, default=1.25,
                    help="enlarge every label without moving the geometry")
    ap.add_argument("--width", type=float, default=1215.0,
                    help="rendered width")
    args = ap.parse_args(argv)

    if args.fasta:
        seq_name, seq = read_fasta(args.fasta, args.protein)
        seq_note = f"in-silico digest of {seq_name} ({len(seq)} aa)"
    else:
        seq = illustrative_sequence(args.length, args.seed)
        seq_name = "illustrative sequence"
        seq_note = (f"illustrative {len(seq)}-residue sequence at human average "
                    f"composition (seed {args.seed})")

    spec = [(n, args.gluc_residues if n == "Glu-C" else r, p)
            for n, r, p in SPECIFICITY]
    colour = assign([n for n, _, _ in spec] + ["All"])

    n = len(seq)
    rng = random.Random(args.seed + 1)
    peptides, masks = {}, {}
    for name, residues, block_p in spec:
        peps = digest(seq, set(residues), block_p,
                      args.missed_cleavages, args.min_len, args.max_len)
        if args.detect_prob < 1.0:
            peps = [p for p in peps if rng.random() < args.detect_prob]
        peptides[name] = peps
        masks[name] = coverage_mask(peps, n)
    depth = [sum(masks[k][i] for k in masks) for i in range(n)]

    # pick the window to show
    win = n if args.window <= 0 else min(args.window, n)
    sites = {nm: set() for nm, _, _ in spec}
    for nm, residues, block_p in spec:
        for i, aa in enumerate(seq[:-1]):
            if aa in set(residues) and not (block_p and seq[i + 1] == "P"):
                sites[nm].add(i + 1)
    if win >= n:
        w0 = 0
    elif args.window_start is not None:
        w0 = max(0, min(args.window_start, n - win))
    elif args.window_rule == "balance":
        # most similar coverage, subject to the --min-union/--min-unique floors
        cand = []
        for start in range(0, max(n - win, 0) + 1):
            cov = {k: sum(masks[k][start:start + win]) / win for k in masks}
            dep = [sum(masks[k][i] for k in masks) for i in range(start, start + win)]
            uni = sum(1 for d in dep if d) / win
            uq = min(sum(1 for i in range(start, start + win)
                         if masks[k][i] and not any(masks[o][i]
                                                    for o in masks if o != k)) / win
                     for k in masks)
            spread = max(cov.values()) - min(cov.values())
            cand.append((uni >= args.min_union and uq >= args.min_unique,
                         -spread, start))
        ok = [c for c in cand if c[0]]
        w0 = max(ok or cand)[2]
        if not ok:
            print("  NOTE: no window meets --min-union/--min-unique; "
                  "took the most balanced anyway")
    else:
        allsites = sorted(s for st in sites.values() for s in st)
        best, w0 = -1, 0
        for start in range(0, max(n - win, 0) + 1):
            k = sum(1 for s in allsites if start < s < start + win)
            if k > best:
                best, w0 = k, start
    w1 = w0 + win

    W, H = 1010, 360
    band2 = 30.0            # offset of the bottom band
    c = Canvas(W, H, args.font, font_scale=args.text_scale, out_w=args.width)
    names = [nm for nm, _, _ in spec]

    if args.letter:
        c.text(22, 32, args.letter, 13, INK, "start", "600")

    # band 1: sample -> three digests -> shared stages -> union
    lanes = [92.0, 136.0, 180.0]
    mid_y = lanes[1]
    # sample box sized to its widest line
    box_lines = [(args.tissue, 10.5), (args.groups, 9.2),
                 (args.protein_amount, 9.2), (args.prep, 9.2)]
    k_type = c.fs * TEXT_BOOST                    # drawn units per font unit
    need = max((len(t) * sz * 0.53 * k_type for t, sz in box_lines if t),
               default=60.0)
    sx0 = 40.0
    sx1 = sx0 + max(112.0, need + 16.0)
    c.rect(sx0, 76, sx1 - sx0, 120, fill=STRIP_FILL, stroke=AXIS, sw=1.1, rx=7)
    scx = (sx0 + sx1) / 2
    c.text(scx, 104, args.tissue, 10.5, INK, "middle", "600")
    for k, line in enumerate([l for l in (args.groups, args.protein_amount,
                                          args.prep) if l]):
        c.text(scx, 128 + k * 9.2 * k_type * 1.2, line, 9.2, INK_MUTED, "middle")

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

    for i, (name, residues, block_p) in enumerate(spec):
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

    # band 2: the in-silico peptides of each protease
    c.add(f'<g transform="translate(0 {band2:g})">')
    c.line(44, 196, W - 30, 196, stroke=GRID, sw=1)

    bx0, bx1 = 116.0, W - 34.0
    pitch = (bx1 - bx0) / win

    def X(res):
        return bx0 + (res - w0) * pitch

    ruler_y = 224.0
    if args.fasta:
        c.text(26, ruler_y - 8,
               f"{args.label} · {seq_name}" if args.label else seq_name,
               9.5, INK, "start", "600")
    tick = 10 if win <= 100 else (20 if win <= 260 else 50)
    first = ((w0 // tick) + 1) * tick
    for r in range(first, w1 + 1, tick):
        c.line(X(r), ruler_y - 5, X(r), ruler_y - 1, stroke=AXIS, sw=1)
        c.text(X(r), ruler_y - 8, str(r), 8.2, INK_MUTED, "middle")

    block_h, gap_h = 12.0, 10.0
    top = ruler_y + 10
    for i, name in enumerate(names):
        by = top + i * (block_h + gap_h)
        col = colour[name]
        c.rect(bx0 - 13, by + 7, 9, 9, fill=col, rx=2)
        c.text(bx0 - 19, by + 15, name, 9.5, INK, "end")
        # a 2px gap keeps consecutive peptides apart
        for s, e in peptides[name]:
            if e <= w0 or s >= w1:
                continue
            a, b = max(s, w0), min(e, w1)
            x, wd = X(a) + 1, X(b) - X(a) - 2
            c.rect(x, by + 6, wd, 12, fill=col, fo=0.85,
                   rx=3 if (s >= w0 and e <= w1) else 0)

    # the union: bars wherever any protease covers
    dy = top + len(names) * (block_h + gap_h)
    c.rect(bx0 - 13, dy + 7, 9, 9, fill=UNION, rx=2)
    c.text(bx0 - 19, dy + 15, UNION_LABEL, 9.5, INK, "end")
    i = w0
    while i < w1:
        while i < w1 and not depth[i]:
            i += 1
        j = i
        while j < w1 and depth[j]:
            j += 1
        if j > i:
            c.rect(X(i) + 1, dy + 6, X(j) - X(i) - 2, 12, fill=UNION, fo=0.85,
                   rx=3)
        i = j

    c.add("</g>")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(c.out())

    print(f"wrote {args.out}  ({W}x{H}px)")
    print(f"  sequence           {seq_name}, {n} aa")
    print(f"  window shown       residues {w0 + 1}-{w1}"
          + ("  (--window-start)" if args.window_start is not None
             else f"  (most balanced)" if args.window_rule == "balance"
             else "  (most cleavage sites)"))
    print(f"  Glu-C specificity  after {', '.join(args.gluc_residues)}  "
          f"<- CHECK: buffer-dependent (E in EPPS/bicarbonate, E+D in phosphate)")
    print(f"  detectable window  {args.min_len}-{args.max_len} aa, "
          f"{args.missed_cleavages} missed cleavage(s)")
    for name in names:
        inwin = [s for s in sites[name] if w0 < s < w1]
        print(f"  {name:<8s} {len(peptides[name]):>4d} peptides, "
              f"{sum(masks[name]) / n:6.1%} of residues, "
              f"{len(inwin):>3d} cut sites in window   {colour[name]}")
    cov = sum(1 for d in depth if d)
    print(f"  {'union':<8s} {cov:>4d} residues covered ({cov / n:.1%}), "
          f"{sum(1 for d in depth if d >= 2)} by >=2 digests")
    return 0


if __name__ == "__main__":
    sys.exit(main())
