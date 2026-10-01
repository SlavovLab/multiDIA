#!/usr/bin/env python3
"""Cohort-level panels for what the extra proteases actually buy.

    python3 extra_why_multienzyme.py 'data/search/*-60min-Phospho.parquet' --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta \
        --sample-regex 'CF_(\\d+)' --sample 2744 --outdir figures/why
"""

import argparse
import glob
import math
import os
import statistics
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_svg import Canvas                                 # noqa: E402
from lib_palette import (AXIS, DEPTH, EMPTY, FONT, GRID, INK, INK_MUTED,
                     INK_SECONDARY, STRIP_FILL, TEXT_BOOST, UNION, assign)     # noqa: E402
from fig1cd_coverage import mark, representative      # noqa: E402
from lib_fasta import read_fasta      # noqa: E402

ORDER = ["GluC", "LysC", "Trypsin"]


def load(paths, precursor_q, protein_q, sample_regex, sample):
    import re
    import pyarrow.parquet as pq
    srx = re.compile(sample_regex) if sample_regex else None
    peps = defaultdict(set)                       # digest -> stripped seqs
    grp_peps = defaultdict(lambda: defaultdict(set))   # group -> digest -> seqs
    groups = defaultdict(set)                     # digest -> protein groups
    iso_pt = defaultdict(set)                     # digest -> isoform accessions
    quant = defaultdict(lambda: defaultdict(dict))    # digest -> group -> run -> q
    import lib_report as rp
    WANT = ["run", "protein_groups", "peptide", "precursor_q",
            "protein_q_run", "proteotypic", "protein_qty"]
    for r in rp.open_reports(paths, order=ORDER):
        path, digest = r.path, r.protease
        have = [f for f in WANT if r.has(f)]
        n = 0
        for b in r.batches(have):
            d = b
            for i in range(b["_n"]):
                run = d["run"][i]
                if srx and sample:
                    m = srx.search(run)
                    if not m or (m.group(1) if m.groups() else m.group(0)) != sample:
                        continue
                grp = d["protein_groups"][i]
                if not grp:
                    continue
                eq = d["precursor_q"][i]
                if eq is not None and eq <= precursor_q:
                    pep = d["peptide"][i]
                    if pep:
                        peps[digest].add(pep)
                        grp_peps[grp][digest].add(pep)
                        n += 1
                    if d.get("proteotypic", [None] * b["_n"])[i]:
                        for a in grp.split(";"):
                            if "-" in a:
                                iso_pt[digest].add(a)
                pgq = d.get("protein_q_run", [None] * b["_n"])[i]
                if pgq is not None and pgq <= protein_q:
                    groups[digest].add(grp)
                    q = d.get("protein_qty", [None] * b["_n"])[i]
                    if q:
                        quant[digest][grp][run] = q
        print(f"  {os.path.basename(path)[:50]:<50s} {digest:<8s} {n:>9,} rows")
    return peps, grp_peps, groups, iso_pt, quant


class Axes:
    """Minimal cartesian frame: hairline grid, one baseline, no chart junk."""

    def __init__(self, c, x, y, w, h, xlim, ylim):
        self.c, self.x, self.y, self.w, self.h = c, x, y, w, h
        self.x0, self.x1 = xlim
        self.y0, self.y1 = ylim

    def X(self, v):
        return self.x + (v - self.x0) / (self.x1 - self.x0 or 1) * self.w

    def Y(self, v):
        return self.y + self.h - (v - self.y0) / (self.y1 - self.y0 or 1) * self.h

    def ygrid(self, ticks, fmt=lambda v: f"{v:,.0f}"):
        for t in ticks:
            self.c.line(self.x, self.Y(t), self.x + self.w, self.Y(t),
                        stroke=GRID, sw=1)
            self.c.text(self.x - 8, self.Y(t) + 3.5, fmt(t), 9, INK_MUTED, "end")

    def xaxis(self, ticks, fmt=lambda v: f"{v:,.0f}"):
        self.c.line(self.x, self.y + self.h, self.x + self.w, self.y + self.h,
                    stroke=AXIS, sw=1)
        for t in ticks:
            self.c.line(self.X(t), self.y + self.h, self.X(t), self.y + self.h + 4,
                        stroke=AXIS, sw=1)
            self.c.text(self.X(t), self.y + self.h + 16, fmt(t), 9,
                        INK_MUTED, "middle")

    def title(self, s, sub=None):
        # a subtitle that would overrun the canvas goes on its own line
        k = self.c.fs * TEXT_BOOST
        x_sub = self.x + len(s) * 6.6 * k + 12
        if sub and x_sub + len(sub) * 4.9 * k > self.c.w - 12:
            self.c.text(self.x, self.y - 18 - 15 * k, s, 11.5, INK, "start", "600")
            self.c.text(self.x, self.y - 16, sub, 9, INK_MUTED, "start")
            return
        self.c.text(self.x, self.y - 18, s, 11.5, INK, "start", "600")
        if sub:
            self.c.text(x_sub, self.y - 18, sub, 9, INK_MUTED, "start")


def nice(hi, target=5):
    if hi <= 0:
        return [0]
    raw = hi / target
    mag = 10 ** math.floor(math.log10(raw))
    step = next((m * mag for m in (1, 2, 2.5, 5, 10) if raw <= m * mag), 10 * mag)
    return [i * step for i in range(int(hi / step) + 2) if i * step <= hi * 1.001]


def save(c, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, name)
    with open(p, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {p}")


def panel_missed(groups, outdir, font, letter=""):
    """Proteoforms trypsin alone does not reach, and which digest supplies them."""
    g, l, t = (groups.get(d, set()) for d in ORDER)
    bars = [
        ("Any non-trypsin digest", len((g | l) - t), UNION),
        ("Glu-C and Lys-C", len((g & l) - t), UNION),
        ("Lys-C only", len(l - g - t), assign(ORDER)["LysC"]),
        ("Glu-C only", len(g - l - t), assign(ORDER)["GluC"]),
    ]
    W, H = 620, 250
    c = Canvas(W, H, font)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    ax = Axes(c, 190, 56, W - 230, H - 100, (0, max(b[1] for b in bars) * 1.12), (0, 1))
    ax.title("Proteoforms missed by trypsin",
             f"recovered from {len(t):,} trypsin proteoforms")
    bh = ax.h / len(bars) * 0.52
    for i, (name, v, col) in enumerate(bars):
        y = ax.y + ax.h * (i + 0.5) / len(bars) - bh / 2
        c.rect(ax.x, y, max(ax.X(v) - ax.x, 1), bh, fill=col, fo=0.85, rx=3)
        c.text(ax.x - 10, y + bh / 2 + 4, name, 9.5, INK, "end")
        c.text(ax.X(v) + 8, y + bh / 2 + 4, f"{v:,}", 9.5, INK_SECONDARY, "start")
    ax.xaxis(nice(max(b[1] for b in bars)))
    c.text(ax.x + ax.w / 2, H - 14, "Proteoforms", 10, INK, "middle")
    save(c, outdir, "missed_by_trypsin.svg")
    return {n: v for n, v, _ in bars}


def panel_redundancy(grp_peps, seqs, outdir, font, letter=""):
    """Residues by how many independent digests cover them."""
    tot = [0, 0, 0, 0]
    prot_multi = prot_any = 0
    for grp, byd in grp_peps.items():
        rep = representative(grp, seqs)
        if not rep:
            continue
        seq = seqs[rep]
        masks = {d: mark(seq, s)[0] for d, s in byd.items()}
        if not masks:
            continue
        prot_any += 1
        multi = False
        for i in range(len(seq)):
            k = sum(m[i] for m in masks.values())
            tot[min(k, 3)] += 1
            if k >= 2:
                multi = True
        prot_multi += multi
    cov = sum(tot[1:])
    W, H = 620, 210
    c = Canvas(W, H, font)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    x, y, w, h = 44, 74, W - 88, 44
    c.text(x, 44, "Independent digests covering each residue", 11.5, INK,
           "start", "600")
    c.text(x, y - 8, f"{cov:,} covered residues", 9, INK_MUTED, "start")
    cx = x
    for k in (1, 2, 3):
        frac = tot[k] / cov
        bw = frac * w
        c.rect(cx, y, max(bw, 1), h, fill=DEPTH[k - 1], rx=2)
        if frac > 0.05:
            c.text(cx + bw / 2, y + h / 2 + 4, f"{frac:.0%}", 10,
                   "#ffffff" if k >= 2 else INK, "middle")
        c.text(cx + bw / 2, y + h + 16,
               f"{k} digest" + ("s" if k > 1 else ""), 9, INK_MUTED, "middle")
        cx += bw
    c.text(x, y + h + 44,
           f"{(tot[2] + tot[3]) / cov:.0%} of covered residues carry peptides from "
           f"two or more independent digests", 9.5, INK_SECONDARY, "start")
    c.text(x, y + h + 60,
           f"{prot_multi:,} of {prot_any:,} proteoforms have at least one such residue",
           9.5, INK_SECONDARY, "start")
    save(c, outdir, "redundancy.svg")
    return {"residues": cov, "by_digests": tot[1:], "proteoforms": prot_any,
            "proteoforms_with_overlap": prot_multi}


def panel_length(peps, outdir, font, letter=""):
    """Peptide length distribution per digest."""
    colour = assign(ORDER)
    hist = {}
    hi = 45
    for d in ORDER:
        h = [0] * (hi + 1)
        for p in peps.get(d, ()):
            h[min(len(p), hi)] += 1
        n = sum(h) or 1
        hist[d] = [v / n for v in h]
    W, H = 620, 300
    c = Canvas(W, H, font)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    top = max(max(v) for v in hist.values())
    ax = Axes(c, 60, 60, W - 100, H - 110, (5, hi), (0, top * 1.08))
    ax.title("Peptide length by digest", "distinct sequences, one sample")
    ax.ygrid(nice(top * 1.08, 4), lambda v: f"{v:.0%}")
    ax.xaxis(list(range(10, hi + 1, 10)))
    for d in ORDER:
        pts = [(ax.X(i), ax.Y(hist[d][i])) for i in range(5, hi + 1)]
        path = "M " + " L ".join(f"{a:.1f} {b:.2f}" for a, b in pts)
        c.add(f'<path d="{path}" fill="none" stroke="{colour[d]}" '
              f'stroke-width="2" stroke-linejoin="round"/>')
        med = statistics.median([len(p) for p in peps[d]]) if peps.get(d) else 0
        c.line(ax.X(med), ax.y + ax.h, ax.X(med), ax.y + ax.h + 7,
               stroke=colour[d], sw=2)
    ly = 62
    for d in ORDER:
        med = statistics.median([len(p) for p in peps[d]]) if peps.get(d) else 0
        c.line(ax.x + ax.w - 96, ly - 4, ax.x + ax.w - 78, ly - 4,
               stroke=colour[d], sw=2)
        c.text(ax.x + ax.w - 72, ly, f"{d}  median {med:.0f}", 9, INK, "start")
        ly += 15
    c.text(ax.x + ax.w / 2, H - 14, "Peptide length (residues)", 10, INK, "middle")
    save(c, outdir, "peptide_length.svg")
    return {d: statistics.median([len(p) for p in peps[d]]) for d in ORDER if peps.get(d)}


def panel_saturation(groups, peps, outdir, font, letter=""):
    """What each additional digest adds, in the order that adds most."""
    colour = assign(ORDER + ["All"])
    steps, chosen = [], []
    have_g, have_p = set(), set()
    remaining = [d for d in ORDER if d in groups]
    while remaining:
        best = max(remaining, key=lambda d: len(have_g | groups[d]))
        chosen.append(best)
        have_g = have_g | groups[best]
        have_p = have_p | peps[best]
        steps.append((" + ".join(chosen), len(have_g), len(have_p), best))
        remaining.remove(best)
    W, H = 620, 300
    c = Canvas(W, H, font)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    ax = Axes(c, 74, 62, W - 120, H - 116, (0, 1), (0, max(s[1] for s in steps) * 1.14))
    ax.title("What each added digest contributes",
             "greedy order: the digest adding most proteoforms first")
    ax.ygrid(nice(max(s[1] for s in steps) * 1.14))
    bw = ax.w / len(steps) * 0.46
    for i, (label, ng, npep, last) in enumerate(steps):
        x = ax.x + ax.w * (i + 0.5) / len(steps) - bw / 2
        prev = steps[i - 1][1] if i else 0
        c.rect(x, ax.Y(ng), bw, ax.y + ax.h - ax.Y(ng), fill=colour[last],
               fo=0.85, rx=3)
        c.text(x + bw / 2, ax.Y(ng) - 8, f"{ng:,}", 9.5, INK, "middle")
        if i:
            c.text(x + bw / 2, ax.Y(ng) - 21, f"+{ng - prev:,}", 8.6,
                   INK_MUTED, "middle")
        c.text(x + bw / 2, ax.y + ax.h + 16, label, 9, INK, "middle")
        c.text(x + bw / 2, ax.y + ax.h + 29, f"{npep:,} peptides", 8.4,
               INK_MUTED, "middle")
    c.line(ax.x, ax.y + ax.h, ax.x + ax.w, ax.y + ax.h, stroke=AXIS, sw=1)
    c.add(f'<g transform="translate(22 {ax.y + ax.h / 2:.1f}) rotate(-90)">'
          f'<text x="0" y="0" font-size="10" text-anchor="middle" fill="{INK}">'
          f'Proteoforms</text></g>')
    save(c, outdir, "saturation.svg")
    return [(s[0], s[1], s[2]) for s in steps]


def panel_isoforms(iso_pt, outdir, font, letter=""):
    """Isoform accessions carrying at least one proteotypic peptide."""
    colour = assign(ORDER + ["All"])
    vals = [(d, len(iso_pt.get(d, set())), colour[d]) for d in ORDER]
    union = set().union(*iso_pt.values()) if iso_pt else set()
    vals.append(("All", len(union), colour["All"]))
    W, H = 520, 290
    c = Canvas(W, H, font)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    ax = Axes(c, 74, 62, W - 110, H - 116, (0, 1), (0, max(v for _, v, _ in vals) * 1.16))
    ax.title("Isoforms with proteotypic evidence",
             "accessions with an isoform suffix")
    ax.ygrid(nice(max(v for _, v, _ in vals) * 1.16))
    bw = ax.w / len(vals) * 0.46
    for i, (name, v, col) in enumerate(vals):
        x = ax.x + ax.w * (i + 0.5) / len(vals) - bw / 2
        c.rect(x, ax.Y(v), bw, ax.y + ax.h - ax.Y(v), fill=col, fo=0.85, rx=3)
        c.text(x + bw / 2, ax.Y(v) - 8, f"{v:,}", 9.5, INK, "middle")
        c.text(x + bw / 2, ax.y + ax.h + 16, name, 9.5, INK, "middle")
    c.line(ax.x, ax.y + ax.h, ax.x + ax.w, ax.y + ax.h, stroke=AXIS, sw=1)
    c.add(f'<g transform="translate(22 {ax.y + ax.h / 2:.1f}) rotate(-90)">'
          f'<text x="0" y="0" font-size="10" text-anchor="middle" fill="{INK}">'
          f'Isoform accessions</text></g>')
    save(c, outdir, "isoform_evidence.svg")
    return {n: v for n, v, _ in vals}


def panel_cv(quant, outdir, font, min_runs=4, letter=""):
    """Across-run CV of protein quantity, per digest."""
    colour = assign(ORDER)
    cvs = {}
    for d in ORDER:
        out = []
        for grp, byrun in quant.get(d, {}).items():
            v = [x for x in byrun.values() if x and x > 0]
            if len(v) >= min_runs:
                m = statistics.mean(v)
                if m > 0:
                    out.append(statistics.stdev(v) / m)
        cvs[d] = sorted(out)
    if not any(cvs.values()):
        print("  cv: not enough runs per digest in this selection — skipped")
        return {}
    hi = 1.0
    W, H = 620, 300
    c = Canvas(W, H, font)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    ax = Axes(c, 64, 62, W - 104, H - 116, (0, hi), (0, 1))
    ax.title("Quantitative precision by digest",
             f"protein groups seen in >= {min_runs} runs")
    ax.ygrid([0, .25, .5, .75, 1.0], lambda v: f"{v:.0%}")
    ax.xaxis([0, .25, .5, .75, 1.0], lambda v: f"{v:.0%}")
    for d in ORDER:
        v = cvs[d]
        if not v:
            continue
        # one point per pixel column
        step = max(1, len(v) // int(ax.w))
        pts = [(ax.X(min(v[i], hi)), ax.Y(i / len(v)))
               for i in range(0, len(v), step)]
        pts.append((ax.X(min(v[-1], hi)), ax.Y(1.0)))
        path = "M " + " L ".join(f"{a:.1f} {b:.2f}" for a, b in pts)
        c.add(f'<path d="{path}" fill="none" stroke="{colour[d]}" stroke-width="2"/>')
    ly = 62
    for d in ORDER:
        if not cvs[d]:
            continue
        med = statistics.median(cvs[d])
        c.line(ax.x + ax.w - 172, ly - 4, ax.x + ax.w - 154, ly - 4,
               stroke=colour[d], sw=2)
        c.text(ax.x + ax.w - 148, ly,
               f"{d}  median {med:.0%}  n = {len(cvs[d]):,}", 9, INK, "start")
        ly += 15
    c.text(ax.x + ax.w / 2, H - 14, "Coefficient of variation across runs", 10,
           INK, "middle")
    c.add(f'<g transform="translate(22 {ax.y + ax.h / 2:.1f}) rotate(-90)">'
          f'<text x="0" y="0" font-size="10" text-anchor="middle" fill="{INK}">'
          f'Cumulative fraction</text></g>')
    save(c, outdir, "quant_cv.svg")
    return {d: (statistics.median(v) if v else None) for d, v in cvs.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reports", nargs="+")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--outdir", default="figures/why")
    ap.add_argument("--sample-regex", default=None)
    ap.add_argument("--sample", default=None)
    ap.add_argument("--precursor-q", type=float, default=0.01)
    ap.add_argument("--protein-q", type=float, default=0.01)
    ap.add_argument("--min-runs", type=int, default=4,
                    help="runs a protein must appear in for the CV panel")
    ap.add_argument("--panels", default="all",
                    help="comma-separated: missed,redundancy,length,saturation,"
                         "isoforms,cv")
    ap.add_argument("--font", default=FONT)
    args = ap.parse_args(argv)

    paths = []
    for pat in args.reports:
        paths.extend(sorted(glob.glob(pat)) or [pat])
    want = ({"missed", "redundancy", "length", "saturation", "isoforms", "cv"}
            if args.panels == "all" else set(args.panels.split(",")))

    print(f"reading reports"
          + (f" (sample {args.sample})" if args.sample else " (all samples)") + ":")
    peps, grp_peps, groups, iso_pt, quant = load(
        paths, args.precursor_q, args.protein_q, args.sample_regex, args.sample)
    seqs = read_fasta(args.fasta)
    print(f"{len(seqs):,} sequences; "
          + "  ".join(f"{d}: {len(peps[d]):,} peptides / {len(groups[d]):,} proteoforms"
                      for d in ORDER if d in peps))

    print("\npanels:")
    if "missed" in want:
        r = panel_missed(groups, args.outdir, args.font, "a")
        print("   ", r)
    if "redundancy" in want:
        r = panel_redundancy(grp_peps, seqs, args.outdir, args.font, "b")
        print("   ", r)
    if "length" in want:
        r = panel_length(peps, args.outdir, args.font, "c")
        print("    median lengths", r)
    if "saturation" in want:
        r = panel_saturation(groups, peps, args.outdir, args.font, "d")
        for label, ng, npep in r:
            print(f"    {label:<26} {ng:>7,} proteoforms  {npep:>8,} peptides")
    if "isoforms" in want:
        r = panel_isoforms(iso_pt, args.outdir, args.font, "e")
        print("   ", r)
    if "cv" in want:
        r = panel_cv(quant, args.outdir, args.font, args.min_runs, "f")
        print("    median CV", r)
    return 0


if __name__ == "__main__":
    sys.exit(main())
