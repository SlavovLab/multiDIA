#!/usr/bin/env python3
"""Fig. 2c and Supplementary Fig. 2: peptide abundance along a protein, one row per donor, relative to the controls."""

import argparse
import collections
import csv
import glob
import math
import os
import statistics
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lib_isoform as ie                                             # noqa: E402
import lib_report as rp                                              # noqa: E402
import supp1_isoform_da as da                                        # noqa: E402
from lib_fasta import gene_map, read_fasta                           # noqa: E402
from lib_svg import Canvas, text_width                               # noqa: E402
from lib_palette import (AXIS, FIG_W, FONT, INK, INK_SECONDARY, TEXT_BOOST,  # noqa: E402
                         assign, display, pt_scale)

ORDER = ["Trypsin", "LysC", "GluC"]
CUT = {"Trypsin": "KR", "LysC": "K", "GluC": "E"}
GROUP = {"Control": "Control", "LBD": "PD"}
PCOL = assign(ORDER)
LOW_FRAC = 0.5
W = FIG_W
FS = pt_scale(W, W)
U = FS * TEXT_BOOST                     # units per printed pt
LOW, MID, HIGH = (33, 102, 172), (247, 247, 247), (178, 24, 43)
CLIP = 2.5
NOT_QUANT = "#d4d3cd"
ROW_PT = 8
LIM, COLS = 5.5, 3



def scan(paths, fasta, want, outdir):
    """Run-normalised log2 abundance per (gene, protease, peptide, donor), and quantified peptides per run."""
    os.makedirs(outdir, exist_ok=True)
    seqs, genes = read_fasta(fasta), gene_map(fasta)
    idx, med, bygene = da.read_reports(paths, seqs, genes)
    med = da.common_medians(idx, med)
    n = 0
    with open(os.path.join(outdir, "peptides.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "protease", "peptide", "donor", "log2"])
        for g in want:
            groups = sorted(bygene.get(g, ()))
            for d in ORDER:
                for p in sorted(set().union(*[idx.get((grp, d), {}).keys() for grp in groups])):
                    members = [(grp, p) for grp in groups if p in idx.get((grp, d), {})]
                    for s, v in sorted(da.sample_levels(members, d, idx, med).items()):
                        w.writerow([g, d, p, s, f"{v:.4f}"])
                        n += 1
    npep = collections.Counter((s, d) for (_g, d), peps in idx.items() for q in peps.values()
                               for s, v in q.items() if v > 0 and (s, d) in med)
    with open(os.path.join(outdir, "runs.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["donor", "protease", "peptides"])
        w.writerows([s, d, npep[(s, d)]] for s, d in sorted(med))
    print(f"  wrote {outdir}/peptides.tsv ({n:,} values) and runs.tsv ({len(med)} runs)")


def load(levels):
    """-> ({gene: {(protease, peptide): {donor: log2}}}, {(donor, protease): quantified peptides})."""
    out = collections.defaultdict(lambda: collections.defaultdict(dict))
    with open(os.path.join(levels, "peptides.tsv"), newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["gene"]][(r["protease"], r["peptide"])][r["donor"]] = float(r["log2"])
    with open(os.path.join(levels, "runs.tsv"), newline="") as fh:
        runs = {(r["donor"], r["protease"]): int(r["peptides"])
                for r in csv.DictReader(fh, delimiter="\t")}
    return out, runs


def low_runs(runs):
    """Runs with fewer than LOW_FRAC of their protease's median quantified peptides."""
    med = {d: statistics.median(n for (_s, dd), n in runs.items() if dd == d) for d in ORDER}
    return {k for k, n in runs.items() if n < LOW_FRAC * med[k[1]]}


def relative(raw, cond):
    """-> {donor: log2 minus the mean of the other controls}, for donors with at least two other controls."""
    ctrl = [s for s in raw if cond.get(s) == "Control"]
    out = {}
    for s, v in raw.items():
        others = [raw[c] for c in ctrl if c != s]
        if len(others) >= 2:
            out[s] = v - statistics.mean(others)
    return out


def place(peptide, gene, seqs, genes):
    """-> (start, end, tag): 1-based span in the gene's canonical entry, else in its first isoform holding it."""
    canon = [a for a in seqs if "-" not in a and genes.get(a) == gene]
    for a in canon + sorted(a for a in seqs if "-" in a and a.split("-")[0] in canon):
        k = seqs[a].find(peptide)
        if k >= 0:
            return k + 1, k + len(peptide), "" if a in canon else a
    return 0, 0, "?"


def missed(peptide, d):
    return sum(a in CUT[d] for a in peptide[:-1])


def columns(gene, data, seqs, genes, cond):
    """-> ((accession, length), canonical peptide columns with values relative to the controls)."""
    cols = []
    for (d, p), vals in data.get(gene, {}).items():
        a, b, tag = place(p, gene, seqs, genes)
        rel = relative(vals, cond)
        if rel and not tag:
            cols.append({"protease": d, "peptide": p, "start": a, "end": b, "values": rel})
    cols.sort(key=lambda c_: (ORDER.index(c_["protease"]), c_["start"], c_["end"]))
    acc = next(a for a in sorted(seqs) if "-" not in a and genes.get(a) == gene)
    return (acc, len(seqs[acc])), cols


def cluster_order(cols, donors):
    """Donors ordered by average-linkage clustering of their peptide profiles, optimal leaf order; missing = 0."""
    import numpy as np
    from scipy.cluster import hierarchy
    if len(donors) < 3:
        return donors
    m = np.array([[c_["values"].get(s, 0.0) for c_ in cols] for s in donors])
    z = hierarchy.linkage(m, "average", metric="euclidean")
    return [donors[i] for i in hierarchy.leaves_list(hierarchy.optimal_leaf_ordering(z, m))]


def junction_starts(acc, seqs):
    """1-based canonical residues where a proteoform's internal deletion starts or ends."""
    out = set()
    for a in seqs:
        if ie.ISO.sub("", a) == acc and a != acc:
            for op, c0, c1, i0, _i1 in ie.events(seqs[acc], seqs[a]):
                if op == "delete" and 0 < i0 < len(seqs[a]):
                    out |= {c0 + 1, c1 + 1}
    return sorted(out)



def shade(v):
    """Blue-white-red for a log2 value, clipped at ±CLIP."""
    t = max(-1.0, min(1.0, v / CLIP))
    end = HIGH if t > 0 else LOW
    t = abs(t)
    return "#" + "".join(f"{round(m + t * (e - m)):02x}" for m, e in zip(MID, end))


def tw(s, size, bold=False):
    return text_width(s, size * U, bold)


def runs_of(cells, lo, hi):
    """Contiguous [a, b) stretches of `cells` inside [lo, hi)."""
    out, i = [], lo
    while i < hi:
        if i not in cells:
            i += 1
            continue
        j = i
        while j < hi and j in cells:
            j += 1
        out.append((i, j))
        i = j
    return out


def along(gene, ref, cols, cond, out, letter=None, key=True, row_pt=ROW_PT, axis_title=True,
          numbered=False):
    """Donors by residues, each residue the mean of its covering peptides, the colour key beside the title; -> [(PD donor, row label)]."""
    acc, n = ref
    rows = [s for g in GROUP for s in cluster_order(cols, sorted(d for d in cond if cond[d] == g))]
    n_ctrl = sum(cond[s] == "Control" for s in rows)
    labels = {s: str(i - n_ctrl + 1) for i, s in enumerate(rows) if i >= n_ctrl} if numbered else {}
    top, rh, gap_g, lane = 30.0, row_pt * U, 5 * U, 4.3 * U
    x0 = 26 + 14 * U + (tw("2744", 10) + 6 * U if numbered else 0)
    x1 = W - 8 - tw(str(n), 10) / 2

    def X(r):
        return x0 + (r - 1) / n * (x1 - x0)
    y_lanes = top + (16 if key else 12) * U
    y0 = y_lanes + 3 * lane + 5 * U
    ys = [y0 + i * rh + (gap_g if i >= n_ctrl else 0) for i in range(len(rows))]
    y_end = ys[-1] + rh
    c = Canvas(W, y_end + (32 if axis_title else 19) * U, FONT, font_scale=FS, out_w=W)
    if letter:
        c.text(20, top, letter, 13, INK, "start", "600")
    title = f"{gene} ({acc}, {n} aa)"
    c.text(62, top, title, 12, INK, "start", "600")
    for k, d in enumerate(ORDER):
        for a, b in runs_of({r for c_ in cols if c_["protease"] == d
                             for r in range(c_["start"], c_["end"] + 1)}, 1, n + 1):
            c.rect(X(a), y_lanes + k * lane, X(b) - X(a), lane - 1.2 * U, fill=PCOL[d])
    if key:
        kx, ky = 62 + tw(title, 12, True) + 24 * U, top
        c.text(kx, ky, "log2", 10, INK)
        kx += tw("log2", 10) + 6 * U
        steps, kw = 40, 110 * U / 40
        for i in range(steps):
            c.rect(kx + i * kw, ky - 9 * U, kw + 0.3, 10 * U,
                   fill=shade(-CLIP + (i + 0.5) * 2 * CLIP / steps))
        for v in (-2, 0, 2):
            c.text(kx + (v + CLIP) / (2 * CLIP) * steps * kw, ky + 10 * U, f"{v:+d}" if v else "0", 10,
                   INK_SECONDARY, "middle")
        kx += steps * kw + 14 * U
        c.rect(kx, ky - 9 * U, 10 * U, 10 * U, fill=NOT_QUANT)
        c.text(kx + 13 * U, ky, "not quantified", 10, INK)
        if kx + 13 * U + tw("not quantified", 10) > W - 8:
            sys.exit(f"{gene}: the colour key does not fit beside the title")
    for s, y in zip(rows, ys):
        if s in labels:
            c.text(x0 - 6 * U, y + rh / 2 + 3.5 * U, labels[s], 10, INK, "end")
        prev, start = None, 1
        for r in range(1, n + 2):
            vals = [c_["values"][s] for c_ in cols
                    if r <= n and c_["start"] <= r <= c_["end"] and s in c_["values"]]
            fill = shade(statistics.mean(vals)) if vals else NOT_QUANT
            if fill != prev or r > n:
                if prev:
                    c.rect(X(start), y + 0.5, X(r) - X(start), rh - 1.0, fill=prev)
                prev, start = fill, r
    for grp, ya, yb in (("Control", ys[0], ys[n_ctrl - 1] + rh), ("LBD", ys[n_ctrl], y_end)):
        c.text(26 + 10 * U, (ya + yb) / 2, GROUP[grp], 10, INK, "middle", "600", rot=-90)
    step = next(s for s in (50, 100, 250, 500, 1000) if n / s <= 6)
    for t in [1] + [t for t in range(step, n, step) if n - t >= step / 2] + [n]:
        c.line(X(t + 0.5), y_end + 1, X(t + 0.5), y_end + 4 * U, stroke=INK, sw=1.0)
        c.text(X(t + 0.5), y_end + 14 * U, str(t), 10, INK_SECONDARY, "middle")
    if axis_title:
        c.text((x0 + x1) / 2, y_end + 27 * U, "residue", 10, INK_SECONDARY, "middle")
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({len(cols)} peptides along {n} residues)")
    return [(s, labels[s]) for s in rows if s in labels]


def peptides(gdata, runs, seq, cond):
    """-> canonical peptides quantified in at least 3 controls, with log2 relative to the controls per donor."""
    low = low_runs(runs)
    out = []
    for (d, p), v in gdata.items():
        k = seq.find(p)
        v = {s: x for s, x in v.items() if (s, d) not in low}
        if k < 0 or sum(cond.get(s) == "Control" for s in v) < 3:
            continue
        out.append({"d": d, "p": p, "a": k + 1, "b": k + len(p), "mc": missed(p, d),
                    "fc": relative(v, cond)})
    return out


def facets(peps, n, marks, pd_rows, out):
    """One panel per PD donor: each peptide a line over its residues at its log2 relative to the controls."""
    rows_n = math.ceil(len(pd_rows) / COLS)
    x_l, gap_x, fh, gap_y = 26 + 30 * U, 34 * U, 45 * U, 10 * U
    fw = (W - 8 - tw(str(n), 10) / 2 - x_l - gap_x * (COLS - 1)) / COLS
    top = 20 * U
    body = rows_n * fh + (rows_n - 1) * gap_y
    c = Canvas(W, top + body + 30 * U, FONT, font_scale=FS, out_w=W)
    keys = [(display(d), PCOL[d], 1.0) for d in ORDER] + [
        ("fully cleaved", INK_SECONDARY, 1.0), ("missed cleavage", INK_SECONDARY, 0.35)]
    kx = x_l
    for lab, colour, op in keys:
        c.line(kx, 10.5 * U, kx + 15 * U, 10.5 * U, stroke=colour, sw=2.4 * U, so=op)
        c.text(kx + 19 * U, 14 * U, lab, 10, INK)
        kx += 19 * U + tw(lab, 10) + 16 * U
    if kx - 16 * U > W - 8:
        sys.exit("the peptide key does not fit one row")
    c.text(26 + 8 * U, top + body / 2, "log2 PD / Control", 10, INK_SECONDARY, "middle", rot=-90)
    c.text(x_l + (W - 8 - x_l) / 2, top + body + 27 * U, "residue", 10, INK_SECONDARY, "middle")
    for i, (s, lab) in enumerate(pd_rows):
        r, k = divmod(i, COLS)
        fx0 = x_l + k * (fw + gap_x)
        fy0 = top + r * (fh + gap_y)

        def X(res):
            return fx0 + (res - 1) / (n - 1) * fw

        def Y(v):
            return fy0 + (LIM - max(-LIM, min(LIM, v))) / (2 * LIM) * fh

        c.rect(fx0, fy0, fw, fh, stroke=AXIS, sw=0.8 * U)
        c.line(fx0, Y(0), fx0 + fw, Y(0), stroke=AXIS, sw=0.8 * U)
        c.text(fx0 + 4 * U, fy0 + 12 * U, f"PD {lab}", 10, INK, "start", "600")
        for m in marks:
            c.line(X(m - 0.5), fy0, X(m - 0.5), fy0 + fh, stroke=AXIS, sw=0.6 * U)
        for v in (4, 0, -4):
            c.line(fx0 - 3 * U, Y(v), fx0, Y(v), stroke=INK, sw=0.8 * U)
            if k == 0:
                c.text(fx0 - 5 * U, Y(v) + 3.5 * U, f"{v:+d}" if v else "0", 10, INK, "end")
        if r == rows_n - 1:
            for t in (1, 50, 100, n):
                c.line(X(t), fy0 + fh, X(t), fy0 + fh + 3 * U, stroke=INK, sw=0.8 * U)
                c.text(X(t), fy0 + fh + 14 * U, str(t), 10, INK_SECONDARY, "middle")
        for q in sorted(peps, key=lambda q: (q["mc"] == 0, ORDER.index(q["d"]))):
            if s in q["fc"]:
                c.line(X(q["a"] - 0.5), Y(q["fc"][s]), X(q["b"] + 0.5), Y(q["fc"][s]),
                       stroke=PCOL[q["d"]], sw=2.2 * U, so=1.0 if q["mc"] == 0 else 0.35)
    with open(out, "w") as fh_:
        fh_.write(c.out())


def profile(gene, data, runs, seqs, genes, cond, out, letter):
    """The gene's donor heatmap, PD rows numbered, over one peptide panel per PD donor."""
    import lib_compose
    ref, cols = columns(gene, data, seqs, genes, cond)
    acc = ref[0]
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmp:
        heat, low = os.path.join(tmp, "heat.svg"), os.path.join(tmp, "facets.svg")
        pd_rows = along(gene, ref, cols, cond, heat, letter=letter, numbered=True, axis_title=False)
        facets(peptides(data[gene], runs, seqs[acc], cond), ref[1], junction_starts(acc, seqs), pd_rows, low)
        lib_compose.compose([heat, low], out, width=W)


def along_panel(genes_, data, seqs, genes, cond, out):
    """Along-the-protein panels of `genes_` stacked, the key on the first."""
    import lib_compose
    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmp:
        parts = []
        for i, g in enumerate(genes_):
            ref, cols = columns(g, data, seqs, genes, cond)
            part = os.path.join(tmp, f"along_{g}.svg")
            along(g, ref, cols, cond, part, key=i == 0, axis_title=i == len(genes_) - 1)
            parts.append(part)
        lib_compose.compose(parts, out, gap=2.0, width=W)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("mode", choices=("scan", "along", "profile"))
    ap.add_argument("items", nargs="+", help="scan: report globs; along/profile: genes")
    ap.add_argument("--genes", nargs="+", help="scan: genes to keep")
    ap.add_argument("--levels", default="derived/peptide_levels", help="scan's output folder")
    ap.add_argument("--fasta", default="data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta")
    ap.add_argument("--metadata", default="data/metadata.xlsx")
    ap.add_argument("--letter", default=None)
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    if args.mode == "scan":
        paths = sorted(p for pat in args.items for p in glob.glob(pat))
        if not paths or not args.genes:
            sys.exit("scan needs reports and --genes")
        scan(paths, args.fasta, args.genes, args.out or args.levels)
        return 0
    data, runs = load(args.levels)
    seqs, genes = read_fasta(args.fasta), gene_map(args.fasta)
    cond = rp.read_metadata(args.metadata)
    if args.mode == "along":
        along_panel(args.items, data, seqs, genes, cond, args.out)
    else:
        profile(args.items[0], data, runs, seqs, genes, cond, args.out, args.letter)
    return 0


if __name__ == "__main__":
    sys.exit(main())
