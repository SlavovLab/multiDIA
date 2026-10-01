#!/usr/bin/env python3
"""Modification sites confirmed by peptides from more than one digest.

    # scan once, cache the sites
    python3 fig3a_phospho_sites.py scan 'data/search/*-60min-Phospho.parquet' --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta \\
        --out derived/mods/sites.tsv

    # the panel
    python3 fig3a_phospho_sites.py panel --sites derived/mods/sites.tsv --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta \\
        --letter d --out figures/fig2d.svg
"""

import argparse
import collections
import csv
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_svg import Canvas                                  # noqa: E402
from lib_palette import (AXIS, DEPTH, FONT, GRID, INK, INK_MUTED, INK_SECONDARY,
                     STRIP_FILL, UNION, assign, display,
                     ink_on, TEXT_BOOST)                                        # noqa: E402

ORDER = ["GluC", "LysC", "Trypsin"]

# Share of the example row given to the pies; MIN_PITCH (units per residue) wins.
PIE_SHARE = 0.42
MIN_PITCH = 11.0
PIE_R_ALONE = 48.0        # pie radius when the pies are drawn alone
LABEL_SEP = 0.26          # min angle (rad) between outside pie labels
NCAP = 5                  # peptide-count classes: 1..4, 5+


def pep_shade(k):
    """fill-opacity for a site placed by `k` distinct peptides (k <= NCAP)."""
    return 0.18 + 0.54 * k / NCAP


PIE_GAP = 24.0          # between the example's right edge and the first pie
PIE_PITCH_R = 2.9       # pie pitch, in radii


def text_width(txt, size, fsc):
    """Estimated drawn width of `txt` in canvas units, after TEXT_BOOST."""
    return 0.55 * size * fsc * TEXT_BOOST * len(txt)


def over_white(hex_fill, alpha):
    """`hex_fill` at `alpha` composited on white, as a hex string."""
    r, g, b = (int(hex_fill[i:i + 2], 16) for i in (1, 3, 5))
    return "#" + "".join(f"{round(v * alpha + 255 * (1 - alpha)):02x}"
                         for v in (r, g, b))


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


def scan(paths, fasta, precursor_q, out):
    import pyarrow.parquet as pq
    from lib_fasta import read_fasta
    seqs = read_fasta(fasta)
    print(f"{len(seqs):,} sequences in {os.path.basename(fasta)}")

    sites = collections.defaultdict(lambda: collections.defaultdict(set))
    unlocated = 0
    import lib_report as rp
    for r in rp.open_reports(sorted(paths), order=ORDER):
        path, dig = r.path, r.protease
        F = ["protein_groups", "modified", "peptide", "precursor_q"]
        for b in r.batches(F):
            g, m, st, q = (b[f] for f in F)
            for k in range(b["_n"]):
                if q[k] is None or q[k] > precursor_q or not (m[k] and g[k] and st[k]):
                    continue
                if "[" not in m[k]:
                    continue
                accs = [a.strip() for a in g[k].split(";") if a.strip()]
                acc = next((a for a in accs if a in seqs), None)
                if acc is None:
                    continue
                start = seqs[acc].find(st[k])
                if start < 0:
                    unlocated += 1
                    continue
                for off, mod in parse_mods(m[k]):
                    if "Carbamidomethyl" in mod:
                        continue
                    sites[(acc, start + off + 1, mod)][dig].add(st[k])
        print(f"  scanned {os.path.basename(path)[:44]:<44s} {dig}")
    if unlocated:
        print(f"  {unlocated:,} peptides not found in their representative "
              f"sequence (they belong to another member of the group)")

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "residue", "modification", "digests",
                    "n_digests", "n_distinct_peptides", "peptides_by_digest"])
        for (acc, res, mod), byd in sorted(sites.items()):
            peps = {p for s in byd.values() for p in s}
            enc = "|".join(f"{d}:" + ",".join(sorted(byd[d]))
                           for d in sorted(byd))
            w.writerow([acc, res, mod, ",".join(sorted(byd)), len(byd),
                        len(peps), enc])
    print(f"  wrote {out}  ({len(sites):,} sites)")
    return sites


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


def summarise(rows, mod):
    sel = [r for r in rows if r["mod"] == mod]
    by_n = collections.Counter(len(r["digests"]) for r in sel)
    combos = collections.Counter(tuple(sorted(r["digests"])) for r in sel)
    return sel, by_n, combos


def gained(sel):
    """-> [(digest added, sites so far, of them placed by >= 2 so far)], best order."""
    sets = [set(r["digests"]) for r in sel]
    left, have, steps = sorted({d for s in sets for d in s}), set(), []

    def score(d):
        h = have | {d}
        return (sum(bool(s & h) for s in sets), sum(len(s & h) >= 2 for s in sets))
    while left:
        d = max(left, key=lambda d_: (score(d_), d_))
        left.remove(d)
        have.add(d)
        steps.append((d,) + score(d))
    return steps


def panel_gained(rows, mod, out, font, letter="", width=1215.0, ts=1.7):
    """Sites as each digest is added, shaded by how many peptides place them."""
    sel = [r for r in rows if r["mod"] == mod]
    if not sel:
        sys.exit(f"no sites for {mod!r}")
    steps = gained(sel)

    W = 1010.0
    fsc = 1.45 * ts / 1.7

    def wide(txt, size):
        return text_width(txt, size, fsc)
    names = [display(d) if i == 0 else f"+ {display(d)}"
             for i, (d, *_rest) in enumerate(steps)]
    ml = 48 + max(wide(t, 10) for t in names) + 12
    mr = 8 + wide(f"{steps[-1][1]:,}", 10) + 6
    mt = 76.0
    bh, pitch = 30.0, 52.0
    pw = W - ml - mr
    H = mt + pitch * (len(steps) - 1) + bh + 58
    c = Canvas(W, H, font, font_scale=fsc, out_w=width)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    c.text(48, 30, "Number of phosphosites detected", 11.5, INK, "start",
           "600")
    total = steps[-1][1]

    def X(v):
        return ml + v / total * pw
    have = set()
    for i, (d, n, conf) in enumerate(steps):
        have.add(d)
        cnt = collections.Counter()
        for r in sel:
            peps = {p_ for dd in have for p_ in r["by_digest"].get(dd, [])}
            if peps:
                cnt[min(len(peps), NCAP)] += 1
        assert sum(cnt.values()) == n
        y = mt + i * pitch
        mid = y + bh / 2 + 4
        c.text(ml - 12, mid, names[i], 10, INK, "end")
        c.text(X(n) + 8, mid, f"{n:,}", 10, INK, "start", "600")
        x, outside = ml, []
        for k in range(1, NCAP + 1):
            v = cnt.get(k, 0)
            if not v:
                continue
            w = X(v) - ml
            assert ink_on(over_white(UNION, pep_shade(k))) == INK
            c.rect(x, y, w, bh, UNION, fo=pep_shade(k))
            lab = f"{k}+" if k == NCAP else f"{k}"
            size = next((sz for sz in (8.2, 6.8) if w >= wide(lab, sz) + 3), None)
            if size:
                c.text(x + w / 2, mid, lab, size, INK, "middle", "600")
            else:
                outside.append((lab, x + w / 2))
            x += w
        pos = []
        for j, (lab, cx) in enumerate(outside):
            px = cx if not j else max(cx, pos[-1] + (wide(outside[j - 1][0], 7.6)
                                                     + wide(lab, 7.6)) / 2 + 4)
            pos.append(px)
        if pos:
            over = pos[-1] + wide(outside[-1][0], 7.6) / 2 - (W - 4)
            if over > 0:
                pos = [q - over for q in pos]
        for (lab, cx), px in zip(outside, pos):
            c.line(cx, y - 0.5, px, y - 5.5, INK_MUTED, 0.7)
            c.text(px, y - 7.5, lab, 7.6, INK_SECONDARY, "middle", "600")
        print(f"    {'+ ' if i else '  '}{display(d):<8} {n:>6,} sites, "
              f"{conf:>5,} by >= 2 digests; by peptides "
              + " ".join(f"{k}:{cnt.get(k, 0)}" for k in range(1, NCAP + 1)))

    ky = mt + pitch * (len(steps) - 1) + bh + 30
    title = "distinct peptides placing the site"
    c.text(ml, ky, title, 9, INK_SECONDARY)
    kx = ml + wide(title, 9) + 4
    for k in range(1, NCAP + 1):
        lab = f"{k}+" if k == NCAP else f"{k}"
        c.rect(kx, ky - 9, 10, 10, UNION, fo=pep_shade(k))
        c.text(kx + 14, ky, lab, 9, INK_SECONDARY)
        kx += 14 + wide(lab, 9) + 16
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}")


def peptide_classes(sel, digest=None):
    """-> Counter {k: sites} by distinct peptides placing each site (all digests if None)."""
    cnt = collections.Counter()
    for r in sel:
        peps = set(r["by_digest"].get(digest, [])) if digest else set(r["peps"])
        if peps:
            cnt[min(len(peps), NCAP)] += 1
    return cnt


BY_DIGEST = ["Trypsin", "LysC", "GluC", None]       # None is All proteases


def panel_by_digest(rows, mod, out, font, letter="", width=1215.0, ts=1.7):
    """One stacked bar per digest and one for All proteases."""
    sel = [r for r in rows if r["mod"] == mod]
    if not sel:
        sys.exit(f"no sites for {mod!r}")
    colour = assign(ORDER + ["All"])
    bars = [(d, colour[d] if d else colour["All"], peptide_classes(sel, d))
            for d in BY_DIGEST]
    top = max(sum(cnt.values()) for _d, _c, cnt in bars)
    step = 1000
    ymax = (top // step + 1) * step

    W = 1010.0
    fsc = 1.45 * ts / 1.7

    def wide(txt, size):
        return text_width(txt, size, fsc)
    ml = 36 + wide(f"{ymax:,}", 9) + 10
    mt, ph = 70.0, 230.0
    x_end = 700.0
    slot = (x_end - ml) / len(bars)
    bw = 80.0
    H = mt + ph + 46
    y0 = mt + ph
    c = Canvas(W, H, font, font_scale=fsc, out_w=width)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    c.text(48, 30, "Number of phosphosites detected", 11.5, INK, "start", "600")

    def Y(v):
        return y0 - v / ymax * ph
    for t in range(0, ymax + 1, step):
        c.line(ml, Y(t), x_end, Y(t), GRID if t else AXIS, 1)
        c.text(ml - 8, Y(t) + 3.5, f"{t:,}", 9, INK_SECONDARY, "end")
    c.text(30, (mt + y0) / 2, "Number of phosphosites", 9.6, INK_SECONDARY, "middle",
           rot=-90)

    lab_h = 7.6 * fsc * TEXT_BOOST
    for i, (d, col, cnt) in enumerate(bars):
        n = sum(cnt.values())
        x = ml + slot * (i + 0.5) - bw / 2
        y, outside = y0, []
        for k in range(1, NCAP + 1):
            v = cnt.get(k, 0)
            if not v:
                continue
            h = v / ymax * ph
            fo = pep_shade(k)
            c.rect(x, y - h, bw, h, col, fo=fo)
            lab = f"{k}+" if k == NCAP else f"{k}"
            size = next((sz for sz in (8.2, 6.8)
                         if h >= 0.75 * sz * fsc * TEXT_BOOST + 3), None)
            if size:
                c.text(x + bw / 2, y - h / 2 + 0.3 * size * fsc * TEXT_BOOST, lab,
                       size, ink_on(over_white(col, fo)), "middle", "600")
            else:
                outside.append((lab, y - h / 2))
            y -= h
        pos = []
        for lab, cy in reversed(outside):
            pos.append(cy if not pos else max(cy, pos[-1] + lab_h))
        if pos and pos[-1] > y0 - 6:
            pos = [q - (pos[-1] - (y0 - 6)) for q in pos]
        xr = x + bw
        for (lab, cy), py in zip(reversed(outside), pos):
            c.line(xr, cy, xr + 7, py, INK_MUTED, 0.7)
            c.text(xr + 9, py + 3.2, lab, 7.6, INK_SECONDARY, "start", "600")
        c.text(x + bw / 2, Y(n) - 7, f"{n:,}", 10, INK, "middle", "600")
        c.text(x + bw / 2, y0 + 18, display(d) if d else "All proteases", 10, INK,
               "middle")
        print(f"    {display(d) if d else 'All':<8} {n:>6,} sites; by peptides "
              + " ".join(f"{k}:{cnt.get(k, 0)}" for k in range(1, NCAP + 1)))

    kx, ky = 740.0, mt + 14
    line = 9 * fsc * TEXT_BOOST + 1
    c.text(kx, ky, "distinct peptides", 9, INK_SECONDARY)
    c.text(kx, ky + line, "placing the site", 9, INK_SECONDARY)
    for j, k in enumerate(range(NCAP, 0, -1)):
        yy = ky + line + 10 + j * 17
        c.rect(kx, yy, 12, 12, INK, fo=pep_shade(k))
        c.text(kx + 18, yy + 10, f"{k}+" if k == NCAP else f"{k}", 9, INK_SECONDARY)
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}")


def log_height(v, decades, ph):
    """Bar height for `v` on a log axis from 1 to 10**decades drawn `ph` tall."""
    import math
    return 0.0 if v <= 1 else min(math.log10(v) / decades, 1.0) * ph


def draw_grouped(c, groups, decades, x0, x_end, ml, mt, ph, bw, gap,
                 rotate=False):
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
            if rotate:
                c.text(x + bw / 2 + 3.5, y0 - h - 4, f"{v:,}", 7.6, INK_SECONDARY,
                       "start", "600", rot=-90)
            else:
                c.text(x + bw / 2, y0 - h - 5, f"{v:,}", 7.6, INK_SECONDARY,
                       "middle", "600")
            x += bw + gap
        c.text(ml + slot * (i + 0.5), y0 + 18, display(d) if d else "All proteases",
               10, INK, "middle")
        print(f"    {display(d) if d else 'All':<8} "
              + " ".join(f"{k}:{cnt.get(k, 0)}" for k in range(1, NCAP + 1)))


def grouped_key(c, x, y, fsc):
    """The 1-5+ class key on one line, starting at `x`."""
    title = "distinct peptides placing the site"
    c.text(x, y, title, 9, INK_SECONDARY)
    kx = x + text_width(title, 9, fsc) + 10
    for k in range(1, NCAP + 1):
        lab = f"{k}+" if k == NCAP else f"{k}"
        c.rect(kx, y - 9, 10, 10, UNION, fo=pep_shade(k))
        c.text(kx + 14, y, lab, 9, INK_SECONDARY)
        kx += 14 + text_width(lab, 9, fsc) + 16
    return kx


def panel_grouped(rows, mod, out, font, letter="", width=1215.0, ts=1.7):
    """Sites per peptide-count class, grouped by digest."""
    sel = [r for r in rows if r["mod"] == mod]
    if not sel:
        sys.exit(f"no sites for {mod!r}")
    groups = [(d, peptide_classes(sel, d)) for d in BY_DIGEST]
    top = max(max(cnt.values()) for _d, cnt in groups)
    decades = len(str(int(top)))

    W = 1010.0
    fsc = 1.45 * ts / 1.7

    def wide(txt, size):
        return text_width(txt, size, fsc)
    ml = 36 + wide(f"{10 ** decades:,}", 9) + 10
    mr, mt, ph = 20.0, 84.0, 230.0
    y0 = mt + ph
    H = y0 + 46
    c = Canvas(W, H, font, font_scale=fsc, out_w=width)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    c.text(48, 30, "Number of phosphosites detected", 11.5, INK, "start", "600")
    grouped_key(c, 48, 56.0, fsc)

    draw_grouped(c, groups, decades, 0.0, W - mr, ml, mt, ph, bw=34.0, gap=3.0)
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}")


def venn_regions(sel, sets):
    """-> {frozenset of digests: sites placed by exactly those digests}."""
    out = {frozenset(c): 0 for k in range(1, len(sets) + 1)
           for c in __import__("itertools").combinations(sets, k)}
    for r in sel:
        key = frozenset(r["digests"]) & frozenset(sets)
        if key:
            out[key] += 1
    return out


def venn_label_points(centres, r, step=0.25):
    """-> {frozenset of circle indices: (x, y)} a label point inside each region."""
    import math
    ox = sum(x for x, _y in centres) / len(centres)
    oy = sum(y for _x, y in centres) / len(centres)

    def inside(x, y):
        return frozenset(i for i, (cx, cy) in enumerate(centres)
                         if (x - cx) ** 2 + (y - cy) ** 2 < r * r)
    n = len(centres)
    out = {frozenset(range(n)): (ox, oy)}
    for k in (1, 2):
        for combo in __import__("itertools").combinations(range(n), k):
            tx = sum(centres[i][0] for i in combo) / k - ox
            ty = sum(centres[i][1] for i in combo) / k - oy
            norm = math.hypot(tx, ty)
            ux, uy = tx / norm, ty / norm
            runs, cur, t = [], [], 0.0
            while t < 3 * r:
                if inside(ox + ux * t, oy + uy * t) == frozenset(combo):
                    cur.append(t)
                elif cur:
                    runs.append(cur)
                    cur = []
                t += step
            if cur:
                runs.append(cur)
            run = max(runs, key=len)
            tm = (run[0] + run[-1]) / 2
            out[frozenset(combo)] = (ox + ux * tm, oy + uy * tm)
    return out


VENN = ["Trypsin", "LysC", "GluC"]


def draw_venn(c, regions, centres, rad, names, tots, fsc):
    """Three-set Venn with region counts and per-circle totals."""
    colour = assign(ORDER + ["All"])
    fill_o = 0.16
    for (cx, cy), d in zip(centres, VENN):
        c.add(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{rad}" fill="{colour[d]}" '
              f'fill-opacity="{fill_o}" stroke="{colour[d]}" stroke-width="1.8"/>')

    def composite(members):
        """The region's fill as seen, each circle's colour laid over white."""
        rgb = [255.0, 255.0, 255.0]
        for i in sorted(members):
            h = colour[VENN[i]]
            for j in range(3):
                rgb[j] = rgb[j] * (1 - fill_o) + int(h[1 + 2 * j:3 + 2 * j], 16) * fill_o
        return "#" + "".join(f"{round(v):02x}" for v in rgb)
    points = venn_label_points(centres, rad)
    for members, (x, y) in points.items():
        n = regions[frozenset(VENN[i] for i in members)]
        assert ink_on(composite(members)) == INK
        c.text(x, y + 0.35 * 10 * fsc * TEXT_BOOST, f"{n:,}", 10, INK, "middle",
               "600")
    line = 10 * fsc * TEXT_BOOST
    for i, ((cx, cy), name, tot) in enumerate(zip(centres, names, tots)):
        if i == 2:
            x, anchor, y = cx + rad * 0.78 + 12, "start", cy + rad * 0.62
        else:
            left = i == 0
            x = cx - rad - 12 if left else cx + rad + 12
            anchor, y = ("end" if left else "start"), cy - rad * 0.35
        c.text(x, y, name, 10, INK, anchor, "600")
        c.text(x, y + line, f"{tot:,}", 10, INK_SECONDARY, anchor)


def panel_venn(rows, mod, out, font, letter="", width=1215.0, ts=1.7):
    """Sites by the digests that place them, as a three-set Venn."""
    import math
    sel = [r_ for r_ in rows if r_["mod"] == mod]
    if not sel:
        sys.exit(f"no sites for {mod!r}")
    regions = venn_regions(sel, VENN)
    total = sum(regions.values())
    fsc = 1.45 * ts / 1.7

    def wide(txt, size):
        return text_width(txt, size, fsc)
    rad, side = 92.0, 108.0
    rc = side / math.sqrt(3)
    names = [display(d) for d in VENN]
    tots = [sum(v for k, v in regions.items() if d in k) for d in VENN]
    lab_w = max(wide(t, 10) for t in names + [f"{v:,}" for v in tots])
    W = 2 * (side / 2 + rad + 14 + lab_w) + 24
    mt = 78.0
    oy = mt + rc / 2 + rad
    ox = W / 2
    centres = [(ox - side / 2, oy - rc / 2), (ox + side / 2, oy - rc / 2),
               (ox, oy + rc)]
    W, H = round(W, 1), round(oy + rc + rad + 18, 1)
    c = Canvas(W, H, font, font_scale=fsc, out_w=width * W / 1010.0)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    c.text(48, 30, "Number of phosphosites detected", 11.5, INK, "start", "600")
    c.text(48, 50, f"n = {total:,} phosphosites", 9.6, INK_SECONDARY, "start")

    draw_venn(c, regions, centres, rad, names, tots, fsc)
    for members, n in sorted(regions.items(), key=lambda kv: (len(kv[0]), sorted(kv[0]))):
        print(f"    {' + '.join(display(d) for d in VENN if d in members):<26} {n:>6,}")
    print(f"    total {total:,}; circles " +
          ", ".join(f"{nm} {t:,}" for nm, t in zip(names, tots)))
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}")


def panel_venn_grouped(rows, mod, out, font, letter="", width=1215.0, ts=1.7):
    """The Venn and the grouped bars as one panel."""
    import math
    sel = [r_ for r_ in rows if r_["mod"] == mod]
    if not sel:
        sys.exit(f"no sites for {mod!r}")
    regions = venn_regions(sel, VENN)
    total = sum(regions.values())
    groups = [(d, peptide_classes(sel, d)) for d in BY_DIGEST]
    decades = len(str(int(max(max(cnt.values()) for _d, cnt in groups))))
    fsc = 1.45 * ts / 1.7

    def wide(txt, size):
        return text_width(txt, size, fsc)
    W = 1010.0
    rad, side = 78.0, 92.0
    rc = side / math.sqrt(3)
    names = [display(d) for d in VENN]
    tots = [sum(v for k, v in regions.items() if d in k) for d in VENN]
    lab_w = max(wide(t, 10) for t in names + [f"{v:,}" for v in tots])
    wv = 2 * (side / 2 + rad + 14 + lab_w) + 24
    vmt = 84.0
    ox, oy = wv / 2, vmt + rc / 2 + rad
    centres = [(ox - side / 2, oy - rc / 2), (ox + side / 2, oy - rc / 2),
               (ox, oy + rc)]
    mt, ph = 84.0, 230.0
    x0 = wv
    ml = x0 + 36 + wide(f"{10 ** decades:,}", 9) + 10
    x_end = W - 20
    slot = (x_end - ml) / len(groups)
    gap = 2.0
    bw = round((0.84 * slot - (NCAP - 1) * gap) / NCAP, 1)
    H = round(max(oy + rc + rad + 18, mt + ph + 46), 1)
    c = Canvas(W, H, font, font_scale=fsc, out_w=width)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    c.text(48, 30, "Number of phosphosites detected", 11.5, INK, "start", "600")
    c.text(48, 54, f"n = {total:,} phosphosites", 9.6, INK_SECONDARY, "start")
    grouped_key(c, x0 + 10, 54.0, fsc)
    draw_venn(c, regions, centres, rad, names, tots, fsc)
    draw_grouped(c, groups, decades, x0, x_end, ml, mt, ph, bw=bw, gap=gap,
                 rotate=True)
    print(f"    Venn {wv:.0f} units wide; bars {bw} units, total {total:,}")
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}")


def panel(rows, seqs, mod, out, font, letter="", example=None, width=1215.0,
          ts=1.7, show=44, label=None, residue=None, flank=4, part="both"):
    """How many sites two or three digests agree on, and one of them drawn out."""
    if part not in ("both", "example", "summary"):
        sys.exit(f"--part {part!r}: expected both, example or summary")
    want_ex = part != "summary"
    want_sum = part != "example"
    sel, by_n, combos = summarise(rows, mod)
    if not sel:
        sys.exit(f"no sites for {mod!r}")
    colour = assign(ORDER + ["All"])
    total = len(sel)

    def distinct(r):
        """All three digests, and each contributed a peptide the others did not."""
        if len(r["digests"]) < 3:
            return False
        longest = {d: max(v, key=len) for d, v in r["by_digest"].items() if v}
        return len(set(longest.values())) == 3
    pool = [r for r in sel if distinct(r)]
    pick = next((r for r in sel if r["acc"] == example
                 and (residue is None or r["res"] == residue)), None) \
        if example else None
    if pick is None:
        if example:
            print(f"  {example} is not a 3-digest site with distinct peptides; "
                  f"choosing automatically")
        pick = max(pool, key=lambda r: max(map(len, r["peps"])) - min(map(len, r["peps"]))) \
            if pool else sel[0]

    ml, mr = 56.0, 34.0
    hy = 32.0 if letter else 24.0
    W = 1010.0
    acc, res = pick["acc"], pick["res"]
    seq = seqs.get(acc, "")
    where = {}
    for d in ORDER:
        for pep in pick["by_digest"].get(d, []):
            i = seq.find(pep)
            if i >= 0 and i < res <= i + len(pep):
                where.setdefault(pep, [i, []])[1].append(d)
    peps = sorted(where, key=lambda p: (where[p][0], len(p)))
    if not peps:
        sys.exit(f"{acc}: no peptide covers residue {res}")
    lo = max(0, min(where[p][0] for p in peps) - flank)
    hi = min(len(seq), max(where[p][0] + len(p) for p in peps) + flank)
    row_h = 13.0

    dist = {}
    for r in sel:
        nd = len(r["digests"])
        dist.setdefault(nd, collections.Counter())[min(len(r["peps"]), NCAP)] += 1
    span = hi - lo
    if part == "example":
        pw = W - ml - mr
        pie_w = pr = pie_pitch = 0.0
    elif part == "summary":
        pw = 0.0
        pr = PIE_R_ALONE
        pie_pitch = PIE_PITCH_R * pr
        pie_w = pie_pitch * (len(dist) - 1) + 2 * pr
        W = ml + pie_w + mr
        width = width * W / 1010.0
    else:
        avail = W - ml - mr - PIE_GAP
        pie_w = avail * PIE_SHARE
        pw = avail - pie_w
        if pw / max(span, 1) < MIN_PITCH:
            pw = min(avail, MIN_PITCH * span)
            pie_w = avail - pw
        pr = pie_w / (PIE_PITCH_R * (len(dist) - 1) + 2)
        pie_pitch = PIE_PITCH_R * pr

    bar_h, bar_gap = 5.0, 2.2
    step = bar_h + bar_gap
    half = (len(peps) + 1) // 2
    above, below = peps[:half], peps[half:]
    rule_y = hy + 12
    top = rule_y + 40
    ruler = top - 12
    sy = top + len(above) * step + 12
    bot = sy + 8
    ex_bottom = bot + len(below) * step
    if want_ex:
        py = (ruler + ex_bottom) / 2 - (30 - 14) / 2
        py = max(py, rule_y + 30 + pr)
    else:
        py = rule_y + 30 + pr
    pie_bottom = py + pr + 30
    H = round(max(ex_bottom if want_ex else 0.0,
                  pie_bottom if want_sum else 0.0) + 14)
    c = Canvas(W, H, font, font_scale=ts * 0.82, out_w=width)

    if letter:
        c.text(22, hy, letter, 13, INK, "start", "600")

    starts = sorted({where[p][0] + 1 for p in peps})
    ends = sorted({where[p][0] + len(p) for p in peps})
    c.line(ml - 12, rule_y, W - mr, rule_y, stroke=GRID, sw=1)
    if want_ex:
        c.text(48, hy, (f"{label} · {acc}" if label else acc)
               + f" · {mod.split(' (')[0].lower()} at {seq[res - 1]}{res} of "
               f"{len(seq)} aa · {len(peps)} peptides", 9.6, INK, "start", "600")
    else:
        c.text(48, hy, f"{mod} · {total:,} sites", 9.6, INK, "start", "600")

    pitch = min(pw / max(span, 1), 17.0)

    def X(pos):
        return ml + (pos - lo) * pitch

    def draw(pep, by):
        """One peptide, striped by digest where more than one made it."""
        i, digs = where[pep]
        h = bar_h / len(digs)
        for k, d in enumerate(digs):
            c.rect(X(i), by + k * h, X(i + len(pep)) - X(i), h,
                   fill=colour[d], fo=0.9, rx=0)

    if want_ex:
        c.line(ml, ruler, ml + pw, ruler, stroke=AXIS, sw=1)
        decades = [p for p in range(lo + 1, hi + 1) if p % 10 == 0]
        ticks = sorted(set(decades) | {p for p in (lo + 1, hi)
                                       if all(abs(p - q) > 2 for q in decades)})
        ticks = [p for p in ticks if abs(p - res) > 2]
        for pos in ticks:
            c.line(X(pos - 1) + pitch / 2, ruler, X(pos - 1) + pitch / 2,
                   ruler - 4, stroke=AXIS, sw=1)
            c.text(X(pos - 1) + pitch / 2, ruler - 7, str(pos), 7.4, INK_MUTED,
                   "middle")

        for k, pep in enumerate(above):
            draw(pep, top + k * step)
        for k, pep in enumerate(below):
            draw(pep, bot + k * step)

        c.rect(X(res - 1), top - 4, pitch, (bot + len(below) * step) - top + 2,
               fill=UNION, fo=0.10, stroke="none", rx=2)
        c.text(X(res - 1) + pitch / 2, top - 23, mod.split(" (")[0].lower(), 7.4,
               UNION, "middle", "600")

        for k in range(lo, hi):
            c.text(X(k) + pitch / 2, sy, seq[k], 8.6,
                   INK if k == res - 1 else INK_SECONDARY, "middle",
                   "600" if k == res - 1 else None)

    import math
    fs = ts * 0.82
    if want_sum:
        x0 = W - mr - pie_w + pr
        for i, nd in enumerate(sorted(dist)):
            cnt = dist[nd]
            tot_n = sum(cnt.values())
            cx = x0 + pie_pitch * i
            ang = -math.pi / 2
            outside = []
            for k in range(1, NCAP + 1):
                v = cnt.get(k, 0)
                if not v:
                    continue
                sweep = 2 * math.pi * v / tot_n
                fo = pep_shade(k)
                assert ink_on(over_white(UNION, fo)) == INK, \
                    f"ramp step {k} is too dark for an ink label"
                if v == tot_n:
                    c.add(f'<circle cx="{cx:.1f}" cy="{py:.1f}" r="{pr}" '
                          f'fill="{UNION}" fill-opacity="{fo:.2f}"/>')
                else:
                    x1 = cx + pr * math.cos(ang)
                    y1 = py + pr * math.sin(ang)
                    x2 = cx + pr * math.cos(ang + sweep)
                    y2 = py + pr * math.sin(ang + sweep)
                    big = 1 if sweep > math.pi else 0
                    c.add(f'<path d="M {cx:.1f} {py:.1f} L {x1:.1f} {y1:.1f} '
                          f'A {pr} {pr} 0 {big} 1 {x2:.1f} {y2:.1f} Z" '
                          f'fill="{UNION}" fill-opacity="{fo:.2f}" '
                          f'stroke="#ffffff" stroke-width="1"/>')
                lab = f"{k}+" if k == NCAP else f"{k}"
                mid = ang + sweep / 2
                want = (0.58 * 8.2 * fs * len(lab) + 4.0) / max(sweep, 1e-9)
                if want <= 0.80 * pr:
                    r_lab = max(want, 0.52 * pr)
                    c.text(cx + r_lab * math.cos(mid),
                           py + r_lab * math.sin(mid) + 3.0, lab, 8.2, INK,
                           "middle", "600")
                else:
                    outside.append((mid, lab))
                ang += sweep
            placed = []
            for mid, lab in sorted(outside):
                a = mid
                if placed and a - placed[-1][0] < LABEL_SEP:
                    a = placed[-1][0] + LABEL_SEP
                placed.append((a, mid, lab))
            for a, mid, lab in placed:
                c.line(cx + pr * math.cos(mid), py + pr * math.sin(mid),
                       cx + (pr + 9) * math.cos(a), py + (pr + 9) * math.sin(a),
                       stroke=INK_MUTED, sw=1.0)
                c.text(cx + (pr + 17) * math.cos(a),
                       py + (pr + 17) * math.sin(a) + 3.0, lab, 7.8,
                       INK_SECONDARY, "middle", "600")

            c.text(cx, py + pr + 16, f"{nd} digest" + ("" if nd == 1 else "s"),
                   8.6, INK, "middle", "600")
            c.text(cx, py + pr + 26, f"{tot_n:,} sites", 7.6, INK_MUTED, "middle")

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}")
    print(f"  {mod}: {total:,} sites; " +
          ", ".join(f"{n} digest{'s' if n > 1 else ''} {by_n.get(n, 0):,}"
                    for n in (1, 2, 3)))
    print(f"  example {acc} {seq[res - 1]}{res} of {len(seq)} aa: "
          f"{len(peps)} peptides over {lo + 1}-{hi}, "
          f"{starts} N-termini, {ends} C-termini")
    return W, H


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("scan", help="read the reports and cache the sites")
    s.add_argument("reports", nargs="+")
    s.add_argument("--fasta", required=True)
    s.add_argument("--precursor-q", type=float, default=0.01)
    s.add_argument("--out", default="derived/mods/sites.tsv")

    p = sub.add_parser("panel", help="draw the panel from a cached scan")
    p.add_argument("--sites", default="derived/mods/sites.tsv")
    p.add_argument("--fasta", required=True)
    p.add_argument("--mod", default="Acetyl (Protein N-term)")
    p.add_argument("--example", default=None, help="accession to draw")
    p.add_argument("--residue", type=int, default=None,
                   help="which site on that accession, when it has several")
    p.add_argument("--part", default="both",
                   choices=("both", "example", "summary", "gained",
                            "by-digest", "venn", "grouped", "venn-grouped"),
                   help="which block to draw")
    p.add_argument("--flank", type=int, default=4,
                   help="residues of context beyond the peptides' own extent")
    p.add_argument("--label", default=None, help="gene name for the example")
    p.add_argument("--letter", default="")
    p.add_argument("--width", type=float, default=1215.0)
    p.add_argument("--text-scale", type=float, default=1.7)
    p.add_argument("--out", default="figures/fig2d.svg")
    p.add_argument("--font", default=FONT)

    r = sub.add_parser("report", help="print the cross-digest summary")
    r.add_argument("--sites", default="derived/mods/sites.tsv")

    args = ap.parse_args(argv)
    if args.cmd == "scan":
        paths = []
        for pat in args.reports:
            paths.extend(sorted(glob.glob(pat)) or [pat])
        scan(paths, args.fasta, args.precursor_q, args.out)
        return 0

    rows = load(args.sites)
    if args.cmd == "panel" and args.part == "venn-grouped":
        panel_venn_grouped(rows, args.mod, args.out, args.font, args.letter,
                           args.width, args.text_scale)
        return 0
    if args.cmd == "panel" and args.part == "grouped":
        panel_grouped(rows, args.mod, args.out, args.font, args.letter,
                      args.width, args.text_scale)
        return 0
    if args.cmd == "panel" and args.part == "venn":
        panel_venn(rows, args.mod, args.out, args.font, args.letter,
                   args.width, args.text_scale)
        return 0
    if args.cmd == "panel" and args.part == "by-digest":
        panel_by_digest(rows, args.mod, args.out, args.font, args.letter,
                        args.width, args.text_scale)
        return 0
    if args.cmd == "panel" and args.part == "gained":
        panel_gained(rows, args.mod, args.out, args.font, args.letter,
                     args.width, args.text_scale)
        return 0
    if args.cmd == "report":
        for mod in sorted({r_["mod"] for r_ in rows}):
            sel, by_n, combos = summarise(rows, mod)
            tot = len(sel)
            multi = tot - by_n.get(1, 0)
            print(f"{mod}: {tot:,} sites, {multi:,} ({multi / tot:.1%}) by >=2 digests")
            for combo, n in combos.most_common(6):
                print(f"    {'+'.join(combo):<26s} {n:>7,}")
        return 0

    from lib_fasta import read_fasta
    seqs = read_fasta(args.fasta)
    panel(rows, seqs, args.mod, args.out, args.font, args.letter, args.example,
          args.width, args.text_scale, label=args.label,
          residue=args.residue, flank=args.flank, part=args.part)
    return 0


if __name__ == "__main__":
    sys.exit(main())
