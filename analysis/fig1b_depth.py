#!/usr/bin/env python3
"""Boxplot + jittered points from a tidy `sample,protease,unit,count` table -> SVG.

    python3 fig1b_depth.py derived/counts/counts.csv --units 'precursors,peptides|protein_isoform_groups' --out figures/fig1b.svg
"""

import argparse
import csv
import hashlib
import os
import statistics
import sys
from collections import defaultdict

from lib_palette import (AXIS, FONT, GRID, INK, INK_MUTED, STRIP_FILL, SURFACE,
                     assign, TEXT_BOOST)
from lib_svg import esc

TEXT_SCALE = 2.12

UNIT_LABELS = {
    "precursors": "Precursors",
    "peptides": "Peptides",
    "protein_isoform_groups": "Protein groups",
    "proteins_canonical": "Protein (canonical)",
}


def quantile(sorted_vals, p):
    """Type-7 quantile, the one R and ggplot's boxplot use."""
    n = len(sorted_vals)
    if n == 0:
        return None
    if n == 1:
        return float(sorted_vals[0])
    h = (n - 1) * p
    lo = int(h)
    hi = min(lo + 1, n - 1)
    return sorted_vals[lo] + (h - lo) * (sorted_vals[hi] - sorted_vals[lo])


def box_stats(vals):
    v = sorted(vals)
    q1, q2, q3 = quantile(v, 0.25), quantile(v, 0.5), quantile(v, 0.75)
    iqr = q3 - q1
    lo_fence, hi_fence = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    inside = [x for x in v if lo_fence <= x <= hi_fence]
    return {
        "q1": q1, "median": q2, "q3": q3,
        "lower": min(inside) if inside else q1,
        "upper": max(inside) if inside else q3,
        "outliers": [x for x in v if x < lo_fence or x > hi_fence],
    }


def jitter(sample, unit, width):
    """Deterministic offset in [-width/2, width/2] so the figure is reproducible."""
    h = hashlib.blake2b(f"{sample}|{unit}".encode(), digest_size=4).digest()
    frac = int.from_bytes(h, "big") / 0xFFFFFFFF
    return (frac - 0.5) * width


def nice_ticks(lo, hi, target=6):
    """Round tick steps (1/2/2.5/5 x 10^k) covering [lo, hi]."""
    if hi <= lo:
        hi = lo + 1
    raw = (hi - lo) / max(target - 1, 1)
    mag = 10 ** (len(str(int(abs(raw)))) - 1) if abs(raw) >= 1 else 10 ** -1
    best = None
    for m in (1, 2, 2.5, 5, 10):
        step = m * mag
        n = (hi - lo) / step
        if best is None or abs(n - (target - 1)) < abs(best[1] - (target - 1)):
            best = (step, n)
    step = best[0]
    # stay strictly inside [lo, hi]
    ticks = []
    t = step * (int(lo / step) - (1 if lo < 0 else 0))
    while t <= hi:
        if t >= lo:
            ticks.append(t)
        t += step
    return ticks


def fmt(v):
    return f"{v:,.0f}" if abs(v) >= 1 else f"{v:g}"


def fmt_compact(v, exact=False):
    """Three significant figures with a K/M suffix: 376550 -> '377K'."""
    a = abs(v)
    if a >= 1e6:
        out = f"{v / 1e6:.2f}M"
        return out.replace(".00M", "M") if exact else out
    if a >= 1e5:
        return f"{v / 1e3:.0f}K"
    if a >= 1e3:
        out = f"{v / 1e3:.1f}K"
        # `exact` (axis ticks) drops a trailing .0; data labels keep it
        return out.replace(".0K", "K") if exact else out
    return fmt(v)


def load(path):
    """-> {unit: {protease: {sample: count}}}"""
    data = defaultdict(lambda: defaultdict(dict))
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh)
        need = {"sample", "protease", "unit", "count"}
        if not need.issubset({c.strip() for c in (rdr.fieldnames or [])}):
            sys.exit(f"{path}: expected columns {sorted(need)}, got {rdr.fieldnames}")
        for row in rdr:
            try:
                c = float(row["count"])
            except (TypeError, ValueError):
                continue
            data[row["unit"].strip()][row["protease"].strip()][row["sample"].strip()] = c
    return data


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("counts_csv")
    ap.add_argument("--out", default="fig.svg")
    ap.add_argument("--units", required=True,
                    help="panel order, comma-separated; '|' starts a new y-scale group")
    ap.add_argument("--xticks", action=argparse.BooleanOptionalAction, default=True,
                    help="category labels under each panel")
    ap.add_argument("--letter", default=None, help="panel letter")
    ap.add_argument("--no-outlier-labels", action="store_true")
    ap.add_argument("--width", type=float, default=1215.0, help="total figure width")
    ap.add_argument("--height", type=float, default=422.0, help="total figure height")
    args = ap.parse_args(argv)

    data = load(args.counts_csv)
    if not data:
        sys.exit(f"{args.counts_csv}: no usable rows")

    groups = [[u.strip() for u in g.split(",") if u.strip()] for g in args.units.split("|")]
    groups = [g for g in groups if g]
    units = [u for g in groups for u in g]
    missing = [u for u in units if u not in data]
    if missing:
        sys.exit(f"unit(s) not in {args.counts_csv}: {missing}\navailable: {sorted(data)}")

    # x order is ascending median; colour comes from assign
    first = data[units[0]]
    singles = [p for p in first if p != "All"]
    singles.sort(key=lambda p: statistics.median(first[p].values()) if first[p] else 0)
    cats = singles + (["All"] if "All" in first else [])
    colour = assign(cats)

    def limits(unit_list, head_px, panel_h):
        """Axis range with pixel headroom for the median label."""
        vals = [v for u in unit_list for p in cats for v in data[u].get(p, {}).values()]
        if not vals:
            return 0.0, 1.0
        lo, hi = min(vals), max(vals)
        span = (hi - lo) or max(abs(hi), 1.0)
        bot = lo - span * 0.05
        f = min(head_px / panel_h, 0.45)
        return bot, bot + (hi - bot) / (1.0 - f)

    ts = TEXT_SCALE * TEXT_BOOST
    sz = {"letter": 13.0 * ts, "strip": 11.5 * ts, "ytick": 9.5 * ts,
          "median": 9.5 * ts, "outlier": 9.5 * ts, "xcat": 10.5 * ts,
          "axis": 11.0 * ts}

    ml, mr, mt = 34.0 + 42.0 * ts, 16.0, 14.0
    mb = (16.0 + 50.0 * ts) if args.xticks else 22.0
    strip_h = 18.0 * ts
    gap_in, gap_out = 14.0, 40.0 + 22.0 * ts
    letter_h = 0.0

    n_panels = sum(len(g) for g in groups)
    gaps = gap_in * sum(max(len(g) - 1, 0) for g in groups) \
        + gap_out * max(len(groups) - 1, 0)
    pw = (args.width - ml - mr - gaps) / n_panels
    ph = args.height - mt - letter_h - strip_h - mb
    if pw < 40 or ph < 60:
        sys.exit(f"--width {args.width:g} --height {args.height:g} leave a {pw:.0f}x{ph:.0f}px "
                 f"panel; enlarge the figure")

    layout = []            # (unit, x, group_index, first_in_group)
    x = ml
    for gi, g in enumerate(groups):
        for j, u in enumerate(g):
            layout.append((u, x, gi, j == 0))
            x += pw + (gap_in if j < len(g) - 1 else 0)
        if gi < len(groups) - 1:
            x += gap_out
    W = x + mr
    H = mt + letter_h + strip_h + ph + mb
    # offset + cap height of the median label, plus a little air
    head_px = 9.0 * ts + sz["median"] * 0.78 + 3.0
    glim = [limits(g, head_px, ph) for g in groups]

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.0f}" height="{H:.0f}" '
        f'viewBox="0 0 {W:.0f} {H:.0f}" font-family=\'{FONT}\'>',
        f'<!-- generated by fig1b_depth.py from {esc(os.path.basename(args.counts_csv))}; '
        f'colours {esc(", ".join(f"{c}={colour[c]}" for c in cats))}; '
        f'y-groups {esc(" | ".join(",".join(g) for g in groups))} -->',
        f'<rect width="{W:.0f}" height="{H:.0f}" fill="{SURFACE}"/>',
    ]

    if args.letter:
        out.append(f'<text x="{ml - 46 * ts:.1f}" y="{mt + 12 * ts:.1f}" '
                   f'font-size="{sz["letter"]:.1f}" '
                   f'font-weight="600" fill="{INK}">{esc(args.letter)}</text>')

    for unit, px0, gi, first_in_group in layout:
        py0 = mt + letter_h + strip_h
        lo, hi = glim[gi]

        def Y(v):
            return py0 + ph - (v - lo) / (hi - lo) * ph

        out.append(f'<rect x="{px0:.1f}" y="{mt + letter_h:.1f}" width="{pw:.1f}" '
                   f'height="{strip_h:.1f}" fill="{STRIP_FILL}" stroke="{AXIS}" '
                   f'stroke-width="1"/>')
        out.append(f'<text x="{px0 + pw / 2:.1f}" '
                   f'y="{mt + letter_h + strip_h / 2 + 4 * ts:.1f}" '
                   f'font-size="{sz["strip"]:.1f}" text-anchor="middle" fill="{INK}">'
                   f'{esc(UNIT_LABELS.get(unit, unit))}</text>')

        ticks = nice_ticks(lo, hi)
        for t in ticks:
            y = Y(t)
            out.append(f'<line x1="{px0:.1f}" y1="{y:.2f}" x2="{px0 + pw:.1f}" y2="{y:.2f}" '
                       f'stroke="{GRID}" stroke-width="1"/>')
            if first_in_group:
                out.append(f'<text x="{px0 - 8:.1f}" y="{y + 3.5 * ts:.2f}" '
                           f'font-size="{sz["ytick"]:.1f}" '
                           f'text-anchor="end" fill="{INK_MUTED}" '
                           f'style="font-variant-numeric:tabular-nums">{fmt_compact(t, exact=True)}</text>')
        out.append(f'<line x1="{px0:.1f}" y1="{py0 + ph:.1f}" x2="{px0 + pw:.1f}" '
                   f'y2="{py0 + ph:.1f}" stroke="{AXIS}" stroke-width="1"/>')

        band = pw / len(cats)
        for ci, cat in enumerate(cats):
            cx = px0 + band * (ci + 0.5)
            vals = list(data[unit].get(cat, {}).values())
            if not vals:
                continue
            col = colour[cat]
            st = box_stats(vals)
            bw = min(band * 0.46, 46.0)

            out.append(f'<line x1="{cx:.1f}" y1="{Y(st["lower"]):.2f}" x2="{cx:.1f}" '
                       f'y2="{Y(st["upper"]):.2f}" stroke="{col}" stroke-width="1" '
                       f'stroke-opacity="0.55"/>')
            for v in (st["lower"], st["upper"]):
                out.append(f'<line x1="{cx - bw * 0.28:.1f}" y1="{Y(v):.2f}" '
                           f'x2="{cx + bw * 0.28:.1f}" y2="{Y(v):.2f}" stroke="{col}" '
                           f'stroke-width="1" stroke-opacity="0.55"/>')
            out.append(f'<rect x="{cx - bw / 2:.1f}" y="{Y(st["q3"]):.2f}" width="{bw:.1f}" '
                       f'height="{max(Y(st["q1"]) - Y(st["q3"]), 0.6):.2f}" fill="{col}" '
                       f'fill-opacity="0.10" stroke="{col}" stroke-width="1" '
                       f'stroke-opacity="0.45"/>')
            out.append(f'<line x1="{cx - bw / 2:.1f}" y1="{Y(st["median"]):.2f}" '
                       f'x2="{cx + bw / 2:.1f}" y2="{Y(st["median"]):.2f}" stroke="{col}" '
                       f'stroke-width="2"/>')

            for sample, v in sorted(data[unit][cat].items()):
                jx = cx + jitter(sample, unit, bw * 1.05)
                out.append(f'<circle cx="{jx:.2f}" cy="{Y(v):.2f}" r="4" fill="{col}" '
                           f'fill-opacity="0.85" stroke="{SURFACE}" stroke-width="2"/>')

            top = max(st["upper"], max(vals))
            out.append(f'<text x="{cx:.1f}" y="{Y(top) - 9 * ts:.2f}" '
                       f'font-size="{sz["median"]:.1f}" '
                       f'text-anchor="middle" fill="{INK}" '
                       f'style="font-variant-numeric:tabular-nums">'
                       f'{fmt_compact(st["median"])}</text>')

            if not args.no_outlier_labels and st["outliers"]:
                # push near-coincident outlier labels apart
                marked = sorted(((v, s) for s, v in data[unit][cat].items()
                                 if v in st["outliers"]), reverse=True)
                prev_y = None
                for v, sample in marked:
                    ty = Y(v) + 3.5
                    if prev_y is not None and ty - prev_y < 9.5 * ts:
                        ty = prev_y + 9.5 * ts
                    prev_y = ty
                    out.append(f'<text x="{cx + bw * 0.62:.1f}" y="{ty:.2f}" '
                               f'font-size="{sz["outlier"]:.1f}" fill="{INK_MUTED}">'
                               f'{esc(sample)}</text>')

            if args.xticks:
                out.append(f'<text x="{cx:.1f}" y="{py0 + ph + 16 * ts:.1f}" '
                           f'font-size="{sz["xcat"]:.1f}" text-anchor="middle" '
                           f'fill="{INK}">{esc(cat)}</text>')

    py0 = mt + letter_h + strip_h
    if args.xticks:
        out.append(f'<text x="{ml + (W - ml - mr) / 2:.1f}" y="{py0 + ph + 44 * ts:.1f}" '
                   f'font-size="{sz["axis"]:.1f}" text-anchor="middle" fill="{INK}">'
                   f'Protease</text>')
    yc = py0 + ph / 2
    out.append(f'<g transform="translate({18 + 10 * ts:.1f} {yc:.1f}) rotate(-90)">'
               f'<text x="0" y="0" font-size="{sz["axis"]:.1f}" text-anchor="middle" fill="{INK}">'
               f'Count / sample</text></g>')
    out.append("</svg>")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\n".join(out) + "\n")

    n_samples = len({s for u in units for p in cats for s in data[u].get(p, {})})
    print(f"wrote {args.out}  ({len(units)} panel(s), {len(cats)} groups, "
          f"{n_samples} samples, {W:.0f}x{H:.0f}px)")
    for unit in units:
        for cat in cats:
            vals = list(data[unit].get(cat, {}).values())
            if vals:
                st = box_stats(vals)
                flag = f"  outliers: {[fmt(o) for o in st['outliers']]}" if st["outliers"] else ""
                print(f"  {UNIT_LABELS.get(unit, unit):<22s} {cat:<10s} "
                      f"n={len(vals):<3d} median {fmt(st['median']):>9s}{flag}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
