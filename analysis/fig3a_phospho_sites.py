#!/usr/bin/env python3
"""Fig. 3a: phosphosites as a digest Venn beside per-digest peptide-count bars.

    python3 fig3a_phospho_sites.py --sites derived/mods/sites_phospho_complete.tsv --letter a --out figures/fig3a.svg
"""

import argparse
import collections
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_svg import Canvas                                  # noqa: E402
from lib_palette import (AXIS, FONT, GRID, INK, INK_SECONDARY, UNION, assign,
                         display, ink_on, TEXT_BOOST)                           # noqa: E402

ORDER = ["GluC", "LysC", "Trypsin"]
NCAP = 5                  # peptide-count classes: 1..4, 5+


def pep_shade(k):
    """fill-opacity for a site placed by `k` distinct peptides (k <= NCAP)."""
    return 0.18 + 0.54 * k / NCAP


def text_width(txt, size, fsc):
    """Estimated drawn width of `txt` in canvas units, after TEXT_BOOST."""
    return 0.55 * size * fsc * TEXT_BOOST * len(txt)


# Modified sequence: _[Nterm]X[mod]YZ_ -- N-terminal mods precede the first residue.
NTERM = re.compile(r"^_((?:\[[^\]]+\])+)")
RESIDUE = re.compile(r"([A-Z])((?:\[[^\]]+\])*)")


def parse_mods(modseq):
    """-> [(offset_in_peptide, modification)] with N-terminal mods at offset 0."""
    out = []
    m = NTERM.match(modseq)
    if m:
        for x in re.findall(r"\[([^\]]+)\]", m.group(1)):
            out.append((0, x))
    body = modseq[m.end():] if m else modseq.strip("_")
    body = body.rstrip("_")
    i = 0
    for tok in RESIDUE.finditer(body):
        for x in re.findall(r"\[([^\]]+)\]", tok.group(2) or ""):
            out.append((i, x))
        i += 1
    return out


def load(path):
    out = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            byd = {}
            for chunk in (r["peptides_by_digest"] or "").split("|"):
                if ":" in chunk:
                    d, ps = chunk.split(":", 1)
                    byd[d] = [x for x in ps.split(",") if x]
            out.append({"acc": r["accession"], "res": int(r["residue"]),
                        "mod": r["modification"],
                        "digests": r["digests"].split(","),
                        "by_digest": byd,
                        "peps": sorted({p for v in byd.values() for p in v})})
    return out


def peptide_classes(sel, digest=None):
    """-> Counter {k: sites} by distinct peptides placing each site (all digests if None)."""
    cnt = collections.Counter()
    for r in sel:
        peps = set(r["by_digest"].get(digest, [])) if digest else set(r["peps"])
        if peps:
            cnt[min(len(peps), NCAP)] += 1
    return cnt


BY_DIGEST = ["Trypsin", "LysC", "GluC", None]       # None is All proteases


def log_height(v, decades, ph):
    """Bar height for `v` on a log axis from 1 to 10**decades drawn `ph` tall."""
    import math
    return 0.0 if v <= 1 else min(math.log10(v) / decades, 1.0) * ph


def draw_grouped(c, groups, decades, x0, x_end, ml, mt, ph, bw, gap):
    """Grouped bars on a log axis: one group per digest, one bar per peptide class."""
    y0 = mt + ph
    for e in range(decades + 1):
        y = y0 - log_height(10 ** e, decades, ph) if e else y0
        c.line(ml, y, x_end, y, GRID if e else AXIS, 1)
        c.text(ml - 8, y + 3.5, f"{10 ** e:,}", 9, INK_SECONDARY, "end")
    c.text(x0 + 30, (mt + y0) / 2, "Number of phosphosites", 9.6, INK_SECONDARY,
           "middle", rot=-90)

    slot = (x_end - ml) / len(groups)
    span = NCAP * bw + (NCAP - 1) * gap
    for i, (d, cnt) in enumerate(groups):
        x = ml + slot * i + (slot - span) / 2
        for k in range(1, NCAP + 1):
            v = cnt.get(k, 0)
            h = log_height(v, decades, ph)
            if h:
                c.rect(x, y0 - h, bw, h, UNION, fo=pep_shade(k))
            c.text(x + bw / 2 + 3.5, y0 - h - 4, f"{v:,}", 7.6, INK_SECONDARY,
                   "start", "600", rot=-90)
            x += bw + gap
        c.text(ml + slot * (i + 0.5), y0 + 18, display(d) if d else "All proteases",
               10, INK, "middle")
        print(f"    {display(d) if d else 'All':<8} "
              + " ".join(f"{k}:{cnt.get(k, 0)}" for k in range(1, NCAP + 1)))


def grouped_key_width(fsc):
    """-> the width `grouped_key` draws."""
    labs = [f"{k}+" if k == NCAP else f"{k}" for k in range(1, NCAP + 1)]
    return (text_width("distinct peptides containing the site", 9, fsc) + 10
            + sum(14 + text_width(lab, 9, fsc) + 16 for lab in labs) - 16)


def grouped_key(c, x, y, fsc):
    """The 1-5+ class key on one line, starting at `x`."""
    title = "distinct peptides containing the site"
    kx = x + text_width(title, 9, fsc) + 10
    c.text(kx - 8, y, title, 9, INK_SECONDARY, "end")    # hugs the first box
    for k in range(1, NCAP + 1):
        lab = f"{k}+" if k == NCAP else f"{k}"
        c.rect(kx, y - 9, 10, 10, UNION, fo=pep_shade(k))
        c.text(kx + 14, y, lab, 9, INK_SECONDARY)
        kx += 14 + text_width(lab, 9, fsc) + 16
    return kx


def venn_regions(sel, sets):
    """-> {frozenset of digests: sites placed by exactly those digests}."""
    out = {frozenset(c): 0 for k in range(1, len(sets) + 1)
           for c in __import__("itertools").combinations(sets, k)}
    for r in sel:
        key = frozenset(r["digests"]) & frozenset(sets)
        if key:
            out[key] += 1
    return out


def venn_layout(regions, sets, R):
    """-> (centres, radii): circle areas match the totals, overlaps fitted."""
    import math
    import numpy as np
    from scipy.optimize import brentq, minimize
    tots = [sum(v for k, v in regions.items() if d in k) for d in sets]
    r = [math.sqrt(t / tots[0]) for t in tots]           # largest = 1
    per = math.pi / tots[0]                               # area per site

    def lens(r1, r2, d):
        if d >= r1 + r2:
            return 0.0
        if d <= abs(r1 - r2):
            return math.pi * min(r1, r2) ** 2
        return (r1 * r1 * math.acos((d * d + r1 * r1 - r2 * r2) / (2 * d * r1))
                + r2 * r2 * math.acos((d * d + r2 * r2 - r1 * r1) / (2 * d * r2))
                - 0.5 * math.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2)
                                  * (d + r1 + r2)))

    def dist(i, j):
        want = per * sum(v for k, v in regions.items() if sets[i] in k and sets[j] in k)
        lo, hi = abs(r[i] - r[j]) + 1e-9, r[i] + r[j]
        return hi if want <= 0 else brentq(lambda d: lens(r[i], r[j], d) - want, lo, hi)
    d01, d02, d12 = dist(0, 1), dist(0, 2), dist(1, 2)
    gx = (d02 ** 2 - d12 ** 2 + d01 ** 2) / (2 * d01)
    gy = math.sqrt(max(d02 ** 2 - gx ** 2, 0.0))
    g = np.linspace(-1.1, 1.1 + d01 + r[1], 700)
    X, Y = np.meshgrid(g, g)
    cell = (g[1] - g[0]) ** 2
    keys = list(regions)

    def loss(q):
        cen = [(0.0, 0.0), (q[0], 0.0), (q[1], q[2])]
        m = [(X - cx) ** 2 + (Y - cy) ** 2 < ri * ri for (cx, cy), ri in zip(cen, r)]
        err = 0.0
        for k in keys:
            z = np.ones_like(X, bool)
            for i, d in enumerate(sets):
                z &= m[i] if d in k else ~m[i]
            err += (z.sum() * cell / per - regions[k]) ** 2
        return err
    q = minimize(loss, [d01, gx, gy], method="Nelder-Mead",
                 options={"xatol": 1e-3, "fatol": 0.5}).x
    cen = [(0.0, 0.0), (q[0], 0.0), (q[1], q[2])]
    return [(R * x, R * y) for x, y in cen], [R * ri for ri in r]


def venn_label_points(centres, radii, box=(0.0, 0.0), step=1.0):
    """-> {frozenset of circle indices: (x, y, fits)} the roomiest point per region.

    `fits` says a `box` (w, h) centred there stays inside the region."""
    import itertools

    def inside(x, y):
        return frozenset(i for i, ((cx, cy), r) in enumerate(zip(centres, radii))
                         if (x - cx) ** 2 + (y - cy) ** 2 < r * r)

    def clear(x, y):
        return min(abs(((x - cx) ** 2 + (y - cy) ** 2) ** 0.5 - r)
                   for (cx, cy), r in zip(centres, radii))
    x0 = min(cx - r for (cx, _y), r in zip(centres, radii))
    x1 = max(cx + r for (cx, _y), r in zip(centres, radii))
    y0 = min(cy - r for (_x, cy), r in zip(centres, radii))
    y1 = max(cy + r for (_x, cy), r in zip(centres, radii))
    w, h = box
    corners = [(dx * w / 2, dy * h / 2) for dx in (-1, 0, 1) for dy in (-1, 0, 1)]
    best = {}
    y = y0
    while y <= y1:
        x = x0
        while x <= x1:
            k = inside(x, y)
            if k:
                fits = all(inside(x + dx, y + dy) == k for dx, dy in corners)
                score = (fits, clear(x, y))
                if k not in best or score > best[k][0]:
                    best[k] = (score, (x, y, fits))
            x += step
        y += step
    n = len(centres)
    return {frozenset(c): best[frozenset(c)][1]
            for k in range(1, n + 1) for c in itertools.combinations(range(n), k)
            if frozenset(c) in best}


VENN = ["Trypsin", "LysC", "GluC"]


def draw_venn(c, regions, centres, radii, names, tots, fsc):
    """Area-proportional three-set Venn with region counts and circle totals."""
    colour = assign(ORDER + ["All"])
    fill_o = 0.16
    for (cx, cy), r, d in zip(centres, radii, VENN):
        c.add(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{colour[d]}" '
              f'fill-opacity="{fill_o}" stroke="{colour[d]}" stroke-width="1.8"/>')

    def composite(members):
        """The region's fill as seen, each circle's colour laid over white."""
        rgb = [255.0, 255.0, 255.0]
        for i in sorted(members):
            h = colour[VENN[i]]
            for j in range(3):
                rgb[j] = rgb[j] * (1 - fill_o) + int(h[1 + 2 * j:3 + 2 * j], 16) * fill_o
        return "#" + "".join(f"{round(v):02x}" for v in rgb)
    line = 10 * fsc * TEXT_BOOST
    for members in venn_label_points(centres, radii):
        n = regions[frozenset(VENN[i] for i in members)]
        for size in (10, 9.5, 9, 8.5, 8, 7.5, 7):    # the largest that fits
            box = (text_width(f"{n:,}", size, fsc) + 2, 0.72 * size * fsc * TEXT_BOOST)
            x, y, fits = venn_label_points(centres, radii, box, step=0.5)[members]
            if fits:
                break
        else:
            sys.exit(f"Venn: {n:,} does not fit its region; enlarge the circles")
        assert ink_on(composite(members)) == INK
        c.text(x, y + 0.35 * size * fsc * TEXT_BOOST, f"{n:,}", size, INK,
               "middle", "600")
    for i, ((cx, cy), r, name, tot) in enumerate(zip(centres, radii, names, tots)):
        if i == 2:
            x, anchor = cx + r + 12, "start"
            span = [x + t / 10 * text_width(name, 10, fsc) for t in range(11)]
            under = max((oy + (ro * ro - (sx - ox_) ** 2) ** 0.5
                         for (ox_, oy), ro in zip(centres[:2], radii[:2])
                         for sx in span if abs(sx - ox_) < ro), default=cy)
            y = max(cy, under + 0.8 * line)
        else:
            left = i == 0
            x = cx - r - 12 if left else cx + r + 12
            anchor, y = ("end" if left else "start"), cy - r * 0.35
        c.text(x, y, name, 10, INK, anchor, "600")
        c.text(x, y + line, f"{tot:,}", 10, INK_SECONDARY, anchor)


def venn_extent(centres, radii):
    """-> (x0, y0, x1, y1) of the circles."""
    return (min(cx - r for (cx, _y), r in zip(centres, radii)),
            min(cy - r for (_x, cy), r in zip(centres, radii)),
            max(cx + r for (cx, _y), r in zip(centres, radii)),
            max(cy + r for (_x, cy), r in zip(centres, radii)))


def venn_at(regions, R, x0, y0):
    """-> `venn_layout` shifted so the circles start at (x0, y0)."""
    centres, radii = venn_layout(regions, VENN, R)
    ex, ey = venn_extent(centres, radii)[:2]
    return [(x - ex + x0, y - ey + y0) for x, y in centres], radii


def panel_venn_grouped(rows, mod, out, letter=""):
    """The Venn and the grouped bars as one panel."""
    sel = [r_ for r_ in rows if r_["mod"] == mod]
    if not sel:
        sys.exit(f"no sites for {mod!r}")
    regions = venn_regions(sel, VENN)
    total = sum(regions.values())
    groups = [(d, peptide_classes(sel, d)) for d in BY_DIGEST]
    decades = len(str(int(max(max(cnt.values()) for _d, cnt in groups))))
    fsc = 1.45

    def wide(txt, size):
        return text_width(txt, size, fsc)
    W = 1010.0
    names = [display(d) for d in VENN]
    tots = [sum(v for k, v in regions.items() if d in k) for d in VENN]
    lab_w = max(wide(t, 10) for t in names + [f"{v:,}" for v in tots])
    centres, radii = venn_at(regions, 105.0, 12 + lab_w + 14, 90.0)
    x1, y1 = venn_extent(centres, radii)[2:]
    wv = x1 + 14 + lab_w + 12
    mt, ph = 84.0, 230.0
    x0 = wv
    ml = x0 + 36 + wide(f"{10 ** decades:,}", 9) + 10
    x_end = W - 20
    slot = (x_end - ml) / len(groups)
    gap = 2.0
    bw = round((0.84 * slot - (NCAP - 1) * gap) / NCAP, 1)
    H = round(max(y1 + 18, mt + ph + 46), 1)
    c = Canvas(W, H, FONT, font_scale=fsc, out_w=1215.0)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    c.text(48, 30, "Number of phosphosites detected", 11.5, INK, "start", "600")
    c.text(48, 54, f"n = {total:,} phosphosites", 9.6, INK_SECONDARY, "start")
    grouped_key(c, min(x0 + 10, x_end - grouped_key_width(fsc)), 54.0, fsc)
    draw_venn(c, regions, centres, radii, names, tots, fsc)
    draw_grouped(c, groups, decades, x0, x_end, ml, mt, ph, bw=bw, gap=gap)
    print(f"    Venn {wv:.0f} units wide; bars {bw} units, total {total:,}")
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", required=True, help="prep_phospho.py export's table")
    ap.add_argument("--mod", default="Phospho (STY)")
    ap.add_argument("--letter", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    panel_venn_grouped(load(args.sites), args.mod, args.out, args.letter)
    return 0


if __name__ == "__main__":
    sys.exit(main())
