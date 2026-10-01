#!/usr/bin/env python3
"""Figure 2a: which digests can measure which proteoform ratios.

    python3 extra_digest_upset.py --iso derived/da_iso/isoform_ratios.tsv \\
        --letter a --out figures/fig2a_upset.svg
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


def combos(path, alpha=0.05):
    """-> ([(digests, n_testable, n_changed)], set_sizes, n_tested, n_changed)."""
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    testable = collections.Counter()
    changed = collections.Counter()
    sizes = collections.Counter()
    for r in rows:
        key = frozenset(d for d in ORDER if d in r["digests"])
        if not key:
            continue
        testable[key] += 1
        for d in key:
            sizes[d] += 1
        if float(r["p"]) < alpha:
            changed[key] += 1
    cols = sorted(testable, key=lambda k: (-testable[k], len(k)))
    out = [(k, testable[k], changed[k]) for k in cols]
    assert sum(t for _k, t, _c in out) == sum(testable.values())
    assert sum(c for _k, _t, c in out) == sum(changed.values())
    return out, sizes, sum(testable.values()), sum(changed.values())


def panel(cols, sizes, n_tested, n_changed, out, font, letter="",
          width=1215.0, ts=1.7, title=None, subtitle=None,
          primary="Ratios", secondary="p < 0.05", union_total=None):
    """The UpSet. `secondary=None` draws one bar row instead of two."""
    colour = assign(ORDER + ["All"])
    W = 1010.0
    ml, mr = 250.0, 60.0
    pw = W - ml - mr
    cw = pw / len(cols)
    bar_top, bar_h = 82.0, 150.0
    dot_r, row_h = 5.2, 21.0
    mat_top = bar_top + bar_h + 14
    n_rows = len(ORDER) + (1 if union_total else 0)
    chg_top = mat_top + n_rows * row_h + 12
    chg_h = 44.0 if secondary else 0.0
    H = chg_top + chg_h + 18.0

    c = Canvas(W, H, font, font_scale=1.45 * ts / 1.7, out_w=width)
    if letter:
        c.text(22, 30, letter, 13, INK, "start", "600")
    c.text(60, 30, title or "Digest accessibility of isoform / canonical ratios",
           11.5, INK, "start", "600")
    if subtitle is not False:
        c.text(60, 47, subtitle or
               f"LBD vs Control  ·  n = {n_tested:,} ratios quantified in "
               f"\u2265 3 samples per group", 8.8, INK_MUTED, "start")

    hi = max(t for _k, t, _c in cols) * 1.16
    for i, (key, t, ch) in enumerate(cols):
        x = ml + i * cw + cw / 2
        col = colour[next(iter(key))] if len(key) == 1 else colour["All"]
        bw = cw * 0.52
        y = bar_top + bar_h - t / hi * bar_h
        c.rect(x - bw / 2, y, bw, bar_top + bar_h - y, fill=col, fo=0.88, rx=2)
        c.text(x, y - 7, f"{t:,}", 9.4, INK, "middle", "600")
        for j, d in enumerate(ORDER):
            cy = mat_top + j * row_h + row_h / 2
            on = d in key
            c.add(f'<circle cx="{x:.1f}" cy="{cy:.1f}" r="{dot_r}" '
                  f'fill="{colour[d] if on else "#e6e5e0"}"/>')
        if len(key) > 1:
            ys = [mat_top + ORDER.index(d) * row_h + row_h / 2 for d in key]
            c.line(x, min(ys), x, max(ys), stroke=colour["All"], sw=2.0)
        chi = max(cc for _k, _t, cc in cols) * 1.3
        cy = chg_top + chg_h - (ch / chi * chg_h if chi else 0)
        if ch and secondary:
            c.rect(x - bw / 2, cy, bw, chg_top + chg_h - cy, fill=col, fo=0.88,
                   rx=2)
            c.text(x, cy - 6, f"{ch}", 8.8, INK_SECONDARY, "middle", "600")

    for j, d in enumerate(ORDER):
        cy = mat_top + j * row_h + row_h / 2
        c.text(ml - 14, cy + 3.4, display(d), 9.4, INK, "end", "600")
        c.text(ml - 14 - 62, cy + 3.4, f"{sizes[d]:,}", 8.8, INK_MUTED, "end")
    if union_total:
        cy = mat_top + len(ORDER) * row_h + row_h / 2
        c.line(ml - 14 - 62 - 34, cy - row_h / 2 + 2, ml - 10,
               cy - row_h / 2 + 2, stroke=GRID, sw=1)
        c.text(ml - 14, cy + 3.8, union_total[0], 9.4, colour["All"], "end",
               "600")
        c.text(ml - 14 - 62, cy + 3.8, f"{union_total[1]:,}", 9.4, INK, "end",
               "600")
    c.text(ml - 14 - 62, mat_top - 6, "set size", 8.2, INK_MUTED, "end")
    c.text(ml - 14, bar_top + bar_h * 0.5, primary, 9.6, INK, "end", "600")
    c.line(ml, bar_top + bar_h, ml + pw, bar_top + bar_h, stroke=AXIS, sw=1)
    if secondary:
        c.text(ml - 14, chg_top + chg_h * 0.62, secondary, 9.6, INK, "end",
               "600")
        c.line(ml, chg_top + chg_h, ml + pw, chg_top + chg_h, stroke=AXIS, sw=1)

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({W:.0f}x{H:.0f} -> {width:.0f} wide)")
    for key, t, ch in cols:
        print(f"  {' + '.join(display(d) for d in ORDER if d in key):<24s} "
              f"{t:>5,}" + (f"   {secondary} {ch:>3}" if secondary else ""))
    return cols


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--iso", default="derived/da_iso/isoform_ratios.tsv")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--letter", default="")
    ap.add_argument("--width", type=float, default=1215.0)
    ap.add_argument("--font", default=None)
    ap.add_argument("--out", default="figures/fig2a_upset.svg")
    args = ap.parse_args(argv)
    from lib_palette import FONT
    cols, sizes, n_t, n_c = combos(args.iso, args.alpha)
    panel(cols, sizes, n_t, n_c, args.out, args.font or FONT, args.letter,
          args.width)
    return 0


if __name__ == "__main__":
    sys.exit(main())
