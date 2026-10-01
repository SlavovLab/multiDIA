#!/usr/bin/env python3
"""Figure 2b: every proteoform ratio that changes, and which digests measured it.

    python3 extra_changed_proteoforms.py --iso derived/da_iso/isoform_ratios.tsv \\
        --letter b --out figures/fig2b_changed.svg
"""

import argparse
import collections
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_svg import Canvas                                 # noqa: E402
from lib_palette import (AXIS, GRID, INK, INK_MUTED, INK_SECONDARY, assign,
                     display)                                     # noqa: E402

ORDER = ["GluC", "LysC", "Trypsin"]


def changed(path, alpha=0.05):
    """-> (changed rows sorted by p, n_tested)."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out = sorted((r for r in rows if float(r["p"]) < alpha),
                 key=lambda r: float(r["p"]))
    assert all(any(d in r["digests"] for d in ORDER) for r in out), \
        "a row names no digest at all"
    return out, len(rows)


def panel(rows, out, font, letter="", n_sig=None, n_tested=None,
          width=1215.0, ts=1.7):
    """One row per changed ratio, in two columns."""
    colour = assign(ORDER + ["All"])
    W = 1010.0
    ncol = 2
    per = (len(rows) + ncol - 1) // ncol
    mt, row_h, bar_h = 88.0, 15.5, 9.0
    gutter = 34.0
    colw = (W - 60 - 40 - gutter) / ncol
    lab_w, dot_w = 132.0, 62.0
    pw = colw - lab_w - dot_w
    H = mt + per * row_h + 26.0
    lim = max(abs(float(r["log2fc_ratio"])) for r in rows) * 1.08

    c = Canvas(W, H, font, font_scale=1.45 * ts / 1.7, out_w=width)
    if letter:
        c.text(22, 30, letter, 13, INK, "start", "600")
    c.text(60, 30, "Differential isoform / canonical ratios by digest "
           "accessibility", 11.5, INK, "start", "600")
    by = collections.Counter()
    for r in rows:
        for d in ORDER:
            if d in r["digests"]:
                by[d] += 1
    c.text(60, 47, f"LBD \u2212 Control, log2  ·  n = {len(rows)} of "
           f"{n_tested:,} at p < 0.05  ·  measurable by "
           + ", ".join(f"{display(d)} {by[d]}" for d in ORDER),
           8.8, INK_MUTED, "start")

    for ci in range(ncol):
        x0 = 60 + ci * (colw + gutter)
        bx = x0 + lab_w + dot_w
        zero = bx + pw / 2

        def X(v, z=zero):
            return z + v / lim * pw / 2

        for t in (-2, -1, 0, 1, 2):
            if abs(t) > lim:
                continue
            c.line(X(t), mt - 10, X(t), mt + per * row_h,
                   stroke=AXIS if t == 0 else GRID, sw=1)
            c.text(X(t), mt - 15, f"{t:+g}" if t else "0", 8.0, INK_MUTED,
                   "middle")
        for j, d in enumerate(ORDER):
            c.add(f'<g transform="translate({x0 + lab_w + 8 + j * 19} '
                  f'{mt - 14}) rotate(-90)"><text x="0" y="0" font-size="7.4" '
                  f'text-anchor="start" fill="{INK_MUTED}">{display(d)}</text></g>')

        for i, r in enumerate(rows[ci * per:(ci + 1) * per]):
            y = mt + i * row_h
            fc = float(r["log2fc_ratio"])
            digs = [d for d in ORDER if d in r["digests"]]
            c.text(x0, y + row_h * 0.76, r["gene"], 9.0, INK, "start", "600")
            c.text(x0 + 74, y + row_h * 0.76, r["isoform_group"], 7.8,
                   INK_MUTED, "start")
            for j, d in enumerate(ORDER):
                c.add(f'<circle cx="{x0 + lab_w + 8 + j * 19}" '
                      f'cy="{y + row_h * 0.5:.1f}" r="4.2" '
                      f'fill="{colour[d] if d in digs else "#e6e5e0"}"/>')
            a_, b_ = (X(fc), X(0)) if fc < 0 else (X(0), X(fc))
            c.rect(a_, y + (row_h - bar_h) / 2, b_ - a_, bar_h,
                   fill=INK_SECONDARY, fo=0.72, rx=1.5)

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({W:.0f}x{H:.0f} -> {width:.0f} wide)")
    print(f"  {len(rows)} changed of {n_tested:,};  " +
          ", ".join(f"{display(d)} {by[d]}" for d in ORDER))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--iso", default="derived/da_iso/isoform_ratios.tsv")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--letter", default="")
    ap.add_argument("--width", type=float, default=1215.0)
    ap.add_argument("--font", default=None)
    ap.add_argument("--out", default="figures/fig2b_changed.svg")
    args = ap.parse_args(argv)

    from lib_palette import FONT
    rows, n_tested = changed(args.iso, args.alpha)
    panel(rows, args.out, args.font or FONT, args.letter, None, n_tested,
          args.width)
    return 0


if __name__ == "__main__":
    sys.exit(main())
