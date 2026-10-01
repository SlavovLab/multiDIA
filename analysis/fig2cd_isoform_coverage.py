#!/usr/bin/env python3
"""Sequence coverage of the regions that discriminate an isoform.

    python3 fig2cd_isoform_coverage.py scan 'data/search/*-60min-Phospho.parquet' \\
        --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta --junctions --tsv derived/da_iso/isoform_disc_coverage.tsv

    python3 fig2cd_isoform_coverage.py --junctions --letter c --out figures/fig2c_disccov.svg
"""

import argparse
import collections
import csv
import glob
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_svg import Canvas                                  # noqa: E402
from lib_palette import (AXIS, FONT, GRID, INK, INK_MUTED, INK_SECONDARY,
                     TEXT_BOOST, assign, display)                                          # noqa: E402

ORDER = ["GluC", "LysC", "Trypsin"]
UNION = "union"
FIELDS = ["isoform", "n_disc"] + [f"cov_{d}" for d in ORDER] + ["cov_union"]
# junction columns, kept apart from the residue ones
J_FIELDS = ["n_junc"] + [f"jcov_{d}" for d in ORDER] + ["jcov_union"]


def peptides_by_base(paths, precursor_q=0.01):
    """-> {canonical accession: {digest: set(peptides)}}, one gene's groups."""
    import fig2b_isoform_strip as ie
    groups = ie.collect_groups(sorted(paths), precursor_q)
    bybase = collections.defaultdict(lambda: collections.defaultdict(set))
    for g, byd in groups.items():
        accs = [a.strip() for a in g.split(";") if a.strip()]
        bases = {ie.ISO.sub("", a) for a in accs}
        if len(bases) != 1:
            continue
        base = bases.pop()
        for d, v in byd.items():
            bybase[base][d] |= v
    return bybase


def own_residues(cseq, iseq):
    """-> (canonical residues missing from the isoform, isoform-only residues)."""
    import fig2b_isoform_strip as ie
    c_own, i_own = set(), set()
    for op, c0, c1, i0, i1 in ie.events(cseq, iseq):
        if op in ("delete", "replace"):
            c_own |= set(range(c0, c1))
        if op in ("insert", "replace"):
            i_own |= set(range(i0, i1))
    return c_own, i_own


def discriminating(cseq, iseq):
    """0-based residues of the isoform that the canonical does not have."""
    return own_residues(cseq, iseq)[1]


def junctions(cseq, iseq):
    """0-based isoform positions just after each internal deletion."""
    import fig2b_isoform_strip as ie
    return sorted({i0 for op, c0, c1, i0, i1 in ie.events(cseq, iseq)
                   if op == "delete" and 0 < i0 < len(iseq)})


def junction_coverage(iseq, joins, byd):
    """-> {digest: set of junctions spanned by one peptide}."""
    cov = {d: set() for d in ORDER}
    for d, peps in byd.items():
        if d not in cov:
            continue
        for p in peps:
            k = iseq.find(p)
            while k != -1:
                cov[d] |= {j for j in joins if k <= j - 1 and j < k + len(p)}
                k = iseq.find(p, k + 1)
    return cov


def junction_peptides(iseq, joins, byd):
    """-> {digest: set of peptides spanning at least one of `joins`}."""
    out = {d: set() for d in ORDER}
    for d, peps in byd.items():
        if d not in out:
            continue
        for p in peps:
            k = iseq.find(p)
            while k != -1:
                if any(k <= j - 1 and j < k + len(p) for j in joins):
                    out[d].add(p)
                    break
                k = iseq.find(p, k + 1)
    return out


def region_coverage(iseq, disc, byd):
    """-> {digest: set of `disc` residues its peptides cover in the isoform}."""
    cov = {d: set() for d in ORDER}
    for d, peps in byd.items():
        if d not in cov:
            continue
        for p in peps:
            j = iseq.find(p)
            while j != -1:
                cov[d] |= set(range(j, j + len(p))) & disc
                j = iseq.find(p, j + 1)
    return cov


def scan(paths, fasta, out, precursor_q=0.01, with_junctions=False):
    """Per isoform: discriminating residues, and how many each protease covers."""
    import fig2b_isoform_strip as ie
    from lib_fasta import read_fasta
    seqs = read_fasta(fasta)
    bybase = peptides_by_base(paths, precursor_q)

    isoforms = collections.defaultdict(list)
    for a in seqs:
        if ie.ISO.search(a):
            isoforms[ie.ISO.sub("", a)].append(a)

    rows, deletions = [], 0
    for base, byd in bybase.items():
        if base not in seqs:
            continue
        cseq = seqs[base]
        for iso in isoforms.get(base, ()):
            iseq = seqs.get(iso)
            if not iseq or iseq == cseq:
                continue
            disc = discriminating(cseq, iseq)
            joins = junctions(cseq, iseq) if with_junctions else []
            if not disc and not joins:
                deletions += 1
                continue
            cov = region_coverage(iseq, disc, byd)
            jcov = junction_coverage(iseq, joins, byd) if joins else \
                {d: set() for d in ORDER}
            row = {"isoform": iso, "n_disc": len(disc),
                   **{f"cov_{d}": len(cov[d]) for d in ORDER},
                   "cov_union": len(set().union(*cov.values()))}
            if with_junctions:
                row.update({"n_junc": len(joins),
                            **{f"jcov_{d}": len(jcov[d]) for d in ORDER},
                            "jcov_union": len(set().union(*jcov.values()))})
            rows.append(row)

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, FIELDS + (J_FIELDS if with_junctions else []),
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    what = "a discriminating region or junction" if with_junctions else \
        "a discriminating region"
    left = "terminal-only deletions" if with_junctions else "pure deletions"
    print(f"  wrote {out}  ({len(rows):,} isoforms with {what}; "
          f"{deletions:,} {left} excluded)")
    return rows


def load(path):
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r["n_disc"] = int(r["n_disc"])
            for d in ORDER + [UNION]:
                r[f"cov_{d}"] = int(r[f"cov_{d}"])
            for k in J_FIELDS:
                if k in r:
                    r[k] = int(r[k])
            rows.append(r)
    return rows


def curves(rows, order=None, cov="cov_", size="n_disc"):
    """-> (seen, {series: [isoforms reaching >= t%]}) for t = 1..100."""
    order = order or ORDER
    seen = [r for r in rows if r.get(size, 0) > 0 and r[f"{cov}{UNION}"] > 0]
    out = {}
    for d in order + [UNION]:
        fr = [r[f"{cov}{d}"] / r[size] for r in seen]
        out[d] = [sum(1 for x in fr if x >= t / 100.0) for t in range(1, 101)]
    return seen, out


SEQ = "#c9c8c1"      # a sequence in the schematics; its changed stretch is INK


def schematic(c, x, y, w, kind, pep):
    """Canonical over isoform, with the peptide that proves the isoform."""
    a, b = x + w * 0.38, x + w * 0.62          # the changed stretch
    rows = (("canonical", y), ("isoform", y + 18))
    for name, yy in rows:
        c.text(x - 6, yy + 6, name, 7.8, INK, "end")
        c.rect(x, yy, a - x, 7, fill=SEQ, rx=1)
        c.rect(b, yy, x + w - b, 7, fill=SEQ, rx=1)
        has = (name == "isoform") == (kind == "own")
        if has:
            c.rect(a, yy, b - a, 7, fill=INK, rx=0)
        else:
            c.line(a, yy + 3.5, b, yy + 3.5, stroke=INK, sw=1.2)
    py = y + 32
    if kind == "own":
        c.rect(a + (b - a) * 0.2, py, (b - a) * 0.6, 4, fill=pep, rx=1)
    else:
        c.rect(a - w * 0.1, py, w * 0.1, 4, fill=pep, rx=1)
        c.line(a, py + 2, b, py + 2, stroke=pep, sw=1.0)
        c.rect(b, py, w * 0.1, 4, fill=pep, rx=1)


def panel(rows, out, font, letter="b", width=1215.0, ts=1.7, order=None,
          with_junctions=False):
    """Isoforms against how much of their discriminating region is covered."""
    order = order or ORDER
    seen, cur = curves(rows, order)
    jn = {}
    if with_junctions:
        # isoforms with a deletion junction spanned by one peptide
        jn = {d: sum(1 for r in rows if r.get("n_junc", 0) and r[f"jcov_{d}"])
              for d in order + [UNION]}
    colour = assign(list(order) + ["All"])
    colour[UNION] = colour["All"]

    W = 1010.0
    # right margin sized for the widest end label at the drawn type
    fs_ = 1.45 * ts / 1.7
    mr = 30.0 + 0.62 * 9.4 * fs_ * TEXT_BOOST * max(
        len(display(k)) for k in list(order) + ["All"])
    ml = 90.0
    mt = 112.0 if with_junctions else 56.0
    ph = 208.0
    bx0, bx1 = 0.0, W - 10.0
    if with_junctions:
        pw = W - ml - 300.0
        bx0 = ml + pw + 40.0
    else:
        pw = W - ml - mr
    H = mt + ph + 52
    c = Canvas(W, H, font, font_scale=1.45 * ts / 1.7, out_w=width)

    c.text(20, 30, letter, 13, INK, "start", "600")
    if not with_junctions:
        c.text(ml - 12, 30, f"n = {len(seen):,} non-canonical isoforms", 9.6,
               INK, "start")

    top = (int(max(cur[UNION] + [jn.get(UNION, 0)]) / 250) + 1) * 250

    def X(t):
        return ml + (t - 1) / 99.0 * pw

    def Y(v):
        return mt + ph - v / top * ph

    for t in range(0, top + 1, 250):
        c.line(ml, Y(t), ml + pw, Y(t), stroke=GRID, sw=1)
        if with_junctions:
            c.line(bx0, Y(t), bx1, Y(t), stroke=GRID, sw=1)
        c.text(ml - 8, Y(t) + 3, f"{t:,}", 8.4, INK_MUTED, "end")
    for t in (1, 25, 50, 75, 100):
        c.text(X(t), mt + ph + 17, f"{t}%", 8.8, INK_MUTED, "middle")
    c.text(26, mt + ph / 2, "non-canonical isoforms", 9.6, INK_SECONDARY, "middle", rot=-90)
    c.text(ml + pw / 2, mt + ph + 36,
           "of the discriminating region covered",
           9.2, INK_SECONDARY, "middle")

    for key in order + [UNION]:
        pts = " ".join(f"{X(t + 1):.1f},{Y(v):.1f}"
                       for t, v in enumerate(cur[key]))
        c.add(f'<polyline points="{pts}" fill="none" stroke="{colour[key]}" '
              f'stroke-width="{2.6 if key == UNION else 1.8}" '
              f'stroke-linejoin="round"/>')
    # Junction bars, on the curves' own y axis, in the proteases' colours.
    if with_junctions:
        x0, x1 = bx0, bx1
        slot = (x1 - x0) / (len(order) + 1)
        sw_ = 170.0
        lab_w = 0.55 * 7.8 * c.fs * TEXT_BOOST * len("canonical") + 6
        for cx_, title, kind in (
                (ml + pw / 2, "Discriminating-region peptides", "own"),
                ((x0 + x1) / 2, "Junction-spanning peptides", "junction")):
            c.text(cx_, 38, title, 9.2, INK, "middle", "600")
            schematic(c, cx_ - (sw_ - lab_w) / 2, 56, sw_, kind, colour[UNION])
        c.text(ml + pw - 6, mt + 16, f"n = {len(seen):,}", 9.2, INK, "end")
        for i, key in enumerate(order + [UNION]):
            bx, v = x0 + slot * i + slot * 0.18, jn[key]
            c.rect(bx, Y(v), slot * 0.64, Y(0) - Y(v), fill=colour[key],
                   fo=0.9, rx=1.5)
            c.text(bx + slot * 0.32, Y(v) - 5, f"{v:,}", 8.4, INK, "middle")
            c.text(bx + slot * 0.32, mt + ph + 17,
                   "All" if key == UNION else display(key), 8.4, INK, "middle")
        c.line(x0, mt + ph, x1, mt + ph, stroke=AXIS, sw=1)

    # direct labels at the right edge, nudged apart where the curves converge
    ends = [] if with_junctions else \
        sorted(((cur[k][-1], k) for k in order + [UNION]), reverse=True)
    placed = []
    for v, k in ends:
        y = Y(v) + 3.4
        # one line of the labels' own drawn type apart, boost included
        step = 9.4 * fs_ * TEXT_BOOST * 1.05
        if placed and y < placed[-1] + step:
            y = placed[-1] + step
        placed.append(y)
        lab = "All" if k == UNION else display(k)
        c.text(ml + pw + 10, y, lab, 9.4, colour[k], "start", "600")

    c.line(ml, mt + ph, ml + pw, mt + ph, stroke=AXIS, sw=1)

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({W:.0f}x{H:.0f})")
    tot = sum(r["n_disc"] for r in seen)
    for k in order + [UNION]:
        cv = sum(r[f"cov_{k}"] for r in seen)
        med = statistics.median(r[f"cov_{k}"] / r["n_disc"] for r in seen)
        print(f"    {display(k) if k in order else 'All':<10s} "
              f"{100 * cv / tot:>5.1f}% of all discriminating residues, "
              f"median per isoform {med:>5.1%}, "
              f"≥50% in {cur[k][49]:>5,} isoforms")
    for k in (order + [UNION]) if with_junctions else ():
        print(f"    {display(k) if k in order else 'All':<10s} a junction "
              f"spanned in {jn[k]:>5,} isoforms")
    return cur


STEPS = [("Trypsin",), ("Trypsin", "LysC"), ("Trypsin", "LysC", "GluC")]


def pd_steps(paths, fasta, want, precursor_q=0.01):
    """Per isoform of `want`: evidence covered by trypsin, + Lys-C, + Glu-C.

    -> [(gene, isoform, units, n_residues, n_junctions, frac per step,
         residue frac per step, junction frac per step, {digest: residue frac},
         {digest: junction frac}, {digest or "All": junction peptides},
         [(canonical start, end, {digest or "All": peptides}) per junction])]
    """
    import fig2b_isoform_strip as ie
    from lib_fasta import gene_map
    from lib_fasta import read_fasta
    seqs = read_fasta(fasta)
    gm = gene_map(fasta)
    gene = lambda a: gm.get(a) or gm.get(a.split("-")[0])
    bybase = peptides_by_base(paths, precursor_q)
    out = []
    for iso in sorted(a for a in seqs if ie.ISO.search(a) and gene(a) in want):
        base = ie.ISO.sub("", iso)
        if base not in seqs or seqs[iso] == seqs[base]:
            continue
        cseq, iseq = seqs[base], seqs[iso]
        disc = discriminating(cseq, iseq)
        joins = junctions(cseq, iseq) if not disc else \
            [j for j in junctions(cseq, iseq)]
        units = len(disc) + len(joins)
        if not units:
            continue
        byd = bybase.get(base, {})
        rc = region_coverage(iseq, disc, byd)
        jc = junction_coverage(iseq, joins, byd)
        fr, rfr, jfr = [], [], []
        for step in STEPS:
            got = set().union(*(rc[d] for d in step))
            jg = set().union(*(jc[d] for d in step))
            fr.append((len(got) + len(jg)) / units)
            rfr.append(len(got) / len(disc) if disc else 0.0)
            jfr.append(len(jg) / len(joins) if joins else 0.0)
        rd = {d: len(rc[d]) / len(disc) if disc else 0.0 for d in ORDER}
        jd = {d: len(jc[d]) / len(joins) if joins else 0.0 for d in ORDER}
        jp = junction_peptides(iseq, joins, byd) if joins else {d: set() for d in ORDER}
        jn = {d: len(jp[d]) for d in ORDER}
        jn["All"] = len(set().union(*jp.values()))
        cut = {i0: (c0 + 1, c1) for op, c0, c1, i0, i1 in ie.events(cseq, iseq)
               if op == "delete" and 0 < i0 < len(iseq)}
        per_j = []
        for j in joins:
            pj = junction_peptides(iseq, [j], byd)
            n = {d: len(pj[d]) for d in ORDER}
            n["All"] = len(set().union(*pj.values()))
            per_j.append(cut[j] + (n,))
        out.append((gene(iso), iso, units, len(disc), len(joins), fr, rfr, jfr,
                    rd, jd, jn, per_j))
    return out


def pd_path(steps, out, font, letter="", width=1215.0, title=None,
            n_total=None):
    """Cumulative coverage bars: own residues left, deletion junctions right."""
    from lib_palette import UNION as UNION_HEX
    colour = assign(list(ORDER) + ["All"])
    seg_col = [colour["Trypsin"], colour["LysC"], colour["GluC"]]
    left = [t for t in steps if t[3] > 0 and t[6][-1] > 1e-9]
    right = [t for t in steps if t[4] > 0 and t[7][-1] > 1e-9]
    left.sort(key=lambda t: (-t[6][-1], t[0], t[1]))
    right.sort(key=lambda t: (-t[7][-1], t[0], t[1]))
    blocks = [(left, 6, "Discriminating residues",
               "of the isoform's own residues covered"),
              (right, 7, "Deletion junctions",
               "of the isoform's deletion junctions spanned")]
    covered = {t[1] for t in left} | {t[1] for t in right}
    gain = len({t[1] for rows, k, _t, _x in blocks for t in rows
                if t[k][-1] > t[k][0] + 1e-9})
    n_tot = n_total if n_total is not None else len(steps)

    W, fs = 1010.0, 1.45

    def wide(txt, size):
        return 0.56 * size * fs * TEXT_BOOST * len(txt)
    mt, row, gap = 118.0, 34.0, 36.0
    n_rows = max(len(left), len(right))
    H = mt + row * n_rows + 60
    bw = (W - 26 - 16 - gap) / 2
    c = Canvas(W, H, font, font_scale=fs, out_w=width)
    c.text(20, 30, letter, 13, INK, "start", "600")
    c.text(26, 30, title or "PD-implicated isoforms: discriminating-region "
           "coverage as digests are added", 11.5, INK, "start", "600")
    c.text(26, 50, f"{len(covered)} of {n_tot} isoforms covered; {gain} gain "
           f"coverage beyond trypsin", 8.6, INK_SECONDARY, "start")
    kx, ky = 26.0, 72.0
    for kind, col, lab in [("seg", seg_col[0], "Trypsin"),
                           ("seg", seg_col[1], "+ Lys-C"),
                           ("seg", seg_col[2], "+ Glu-C"),
                           ("dot", UNION_HEX, "All three")]:
        if kind == "dot":
            c.add(f'<circle cx="{kx + 6:.1f}" cy="{ky - 4:.1f}" r="5" fill="{col}"/>')
            kx += 16
        else:
            c.rect(kx, ky - 6.5, 22, 5, fill=col, rx=0)
            kx += 28
        c.text(kx, ky, lab, 8.2, INK_SECONDARY, "start")
        kx += wide(lab, 8.2) + 22

    for bi, (rows, k, head, axis) in enumerate(blocks):
        x0 = 26 + bi * (bw + gap)
        # bold names run wider than the estimate, hence the 1.15
        acc_x = x0 + max((1.15 * wide(t[0], 9.4) for t in rows), default=0) + 10
        ml = acc_x + max((wide(t[1], 8.4) for t in rows), default=0) + 12
        pr = x0 + bw - wide("100%", 8.4) - 10     # plot's right edge

        def X(f, ml=ml, pr=pr):
            return ml + f * (pr - ml)
        c.text(x0, mt - 22, f"{head} · {len(rows)}", 9.8, INK, "start", "600")
        y_end = mt + row * len(rows) - 6
        for t in (0, 0.5, 1.0):
            c.line(X(t), mt - 10, X(t), y_end, stroke=GRID, sw=1.0)
            c.text(X(t), y_end + 18, f"{100 * t:.0f}%", 8.0, INK_SECONDARY,
                   "middle")
        c.text((ml + pr) / 2, y_end + 38, axis, 8.6, INK_SECONDARY, "middle")
        for i, t in enumerate(rows):
            g, iso, fr = t[0], t[1], t[k]
            y = mt + i * row + row / 2 - 6
            c.text(x0, y + 4, g, 9.4, INK, "start", "600")
            c.text(acc_x, y + 4, iso, 8.4, INK_SECONDARY, "start")
            c.line(X(0), y, X(1), y, stroke=GRID, sw=1.0)
            edges = [0.0] + list(fr)             # trypsin starts at 0
            for s_ in range(3):
                if edges[s_ + 1] > edges[s_] + 1e-9:
                    c.rect(X(edges[s_]), y - 3, X(edges[s_ + 1]) - X(edges[s_]),
                           6, fill=seg_col[s_], rx=0)
            c.add(f'<circle cx="{X(fr[-1]):.1f}" cy="{y:.1f}" r="5.5" '
                  f'fill="{UNION_HEX}"/>')
            c.text(pr + 10, y + 4, f"{100 * fr[-1]:.0f}%", 8.4, INK, "start")

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({len(covered)} of {n_tot} isoforms covered: "
          f"{len(left)} on residues, {len(right)} on junctions; {gain} gain "
          f"beyond trypsin)")
    for rows, k, head, _x in blocks:
        print(f"  {head}:")
        for t in rows:
            print(f"    {t[0]:8s} {t[1]:12s} {t[3]:4d} aa {t[4]} junc  "
                  + " -> ".join(f"{100 * f:.0f}%" for f in t[k]))


GRID_BINS = (0.25, 0.50, 0.75, 1.0)          # upper edges; 100% is its own bin


def grid_fill(f):
    """Fig. 3c's cell colour for a coverage fraction."""
    from fig3bc_phospho_atlas import NONE_FILL, STEPS as CELL
    if f <= 1e-9:
        return NONE_FILL
    # binned on the rounded percent, matching the quoted number
    if f >= 1.0 - 1e-9:                      # the top bin means complete
        return CELL[-1]
    pct = min(99, max(1, round(100 * f)))
    return CELL[next(i for i, b in enumerate(GRID_BINS) if pct <= round(100 * b))]


def count_fill(n):
    """Fig. 3a's peptide-count colour."""
    from fig3a_phospho_sites import NCAP, over_white, pep_shade
    from lib_palette import UNION as UNION_HEX
    from fig3bc_phospho_atlas import NONE_FILL
    return NONE_FILL if n <= 0 else over_white(UNION_HEX, pep_shade(min(n, NCAP)))


def grid_block(c, rows, k, cum, head, x, mt, fs, cw=26.0, ch=24.0, gap_all=8.0):
    """Draw one block of `pd_grid`'s cells from `x`. -> right edge."""
    def wide(txt, size):
        return 0.56 * size * fs * TEXT_BOOST * len(txt)
    rows_d = ["GluC", "LysC", "Trypsin"]
    lab_w = wide("Trypsin", 9.6) + 12
    grid_h = ch * 3 + gap_all + ch
    gx = x + lab_w
    c.text(x, mt - 16, f"{head} · {len(rows)}", 9.8, INK, "start", "600")
    for r, d in enumerate(rows_d + ["All"]):
        y = mt + r * ch + (gap_all if d == "All" else 0)
        c.text(gx - 8, y + ch / 2 + 4, display(d) if d != "All" else "All",
               9.6, INK, "end", "600" if d == "All" else None)
        for j, t in enumerate(rows):
            cx = gx + j * cw
            name = display(d) if d != "All" else "All"
            if k == 10:                      # junctions: spanning peptides
                n = t[k][d]
                tip, col = f"{n} peptide{'s' * (n != 1)}", count_fill(n)
            else:
                f = t[cum][-1] if d == "All" else t[k][d]
                tip, col = f"{100 * f:.0f}%", grid_fill(f)
            c.add(f'<g><title>{t[0]} {t[1]} · {name}: {tip}</title>')
            c.rect(cx + 1, y + 1, cw - 2, ch - 2, col)
            c.add('</g>')
    for j, t in enumerate(rows):
        lx = gx + j * cw + cw / 2
        ly = mt + grid_h + 12
        c.text(lx, ly, f"{t[0]} {t[1]}", 8.4, INK_SECONDARY, "end", rot=-55)
    return gx + len(rows) * cw


def pd_grid(steps, out, font, letter="", width=1215.0, title=None, n_total=None):
    """The PD coverage as an isoform x protease grid."""
    from fig3bc_phospho_atlas import NONE_FILL, STEPS as CELL
    left = [t for t in steps if t[3] > 0 and t[6][-1] > 1e-9]
    right = [t for t in steps if t[4] > 0 and t[7][-1] > 1e-9]
    left.sort(key=lambda t: (-t[6][-1], t[0], t[1]))
    right.sort(key=lambda t: (-t[7][-1], t[0], t[1]))
    blocks = [(left, 8, 6, "Discriminating residues"),
              (right, 10, 7, "Deletion junctions")]
    covered = {t[1] for t in left} | {t[1] for t in right}
    gain = len({t[1] for rows, _k, cum, _h in blocks for t in rows
                if t[cum][-1] > t[cum][0] + 1e-9})
    n_tot = n_total if n_total is not None else len(steps)
    rows_d = ["GluC", "LysC", "Trypsin"]

    W, fs = 1010.0, 1.45

    def wide(txt, size):
        return 0.56 * size * fs * TEXT_BOOST * len(txt)
    cw, ch, gap_all, block_gap = 26.0, 24.0, 8.0, 40.0
    lab_w = wide("Trypsin", 9.6) + 12
    mt = 110.0
    grid_h = ch * 3 + gap_all + ch
    label_h = max(wide(f"{t[0]} {t[1]}", 8.4) for t in left + right) * 0.72 + 10
    H = mt + grid_h + label_h + 92
    c = Canvas(W, H, font, font_scale=fs, out_w=width)
    c.text(20, 30, letter, 13, INK, "start", "600")
    c.text(26, 30, title or "PD-implicated isoforms: discriminating evidence "
           "covered by each protease", 11.5, INK, "start", "600")
    c.text(26, 50, f"{len(covered)} of {n_tot} isoforms covered; {gain} gain "
           f"coverage beyond trypsin", 8.6, INK_SECONDARY, "start")

    x, key_x = 26.0, []
    for rows, k, _cum, head in blocks:
        key_x.append(x)
        x = grid_block(c, rows, k, _cum, head, x, mt, fs, cw, ch, gap_all) + block_gap

    # one key per block, under it: residues on 3c's bins, junctions on counts
    ky = H - 22
    keys = [("own residues covered", ["1–25%", "26–50%", "51–75%", "76–99%",
                                      "100%"], list(CELL)),
            ("peptides spanning a junction", ["1", "2", "3", "4", "5+"],
             [count_fill(n) for n in range(1, 6)])]
    # two lines under the figure: the residue key is wider than its block
    head_w = max(wide(h, 9) for h, _l, _f in keys) + 14
    for line, (head, labels, fills) in enumerate(keys):
        ky = H - 44 + 24 * line
        c.text(26, ky + 1, head, 9, INK_SECONDARY)
        kx = 26 + head_w
        for lab, col in zip(labels, fills):
            c.rect(kx, ky - 9, 12, 12, col)
            c.text(kx + 17, ky + 1, lab, 8.6, INK_SECONDARY)
            kx += 17 + wide(lab, 8.6) + 12
        c.rect(kx, ky - 9, 12, 12, NONE_FILL)
        c.text(kx + 17, ky + 1, "none", 8.6, INK_SECONDARY)

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({len(covered)} of {n_tot} isoforms covered: "
          f"{len(left)} on residues, {len(right)} on junctions; {gain} gain "
          f"beyond trypsin)")


def pd_grid_dots(steps, out, font, letter="", width=1215.0, title=None,
                 n_total=None, lengths=False):
    """Residue bars over junction bars; `lengths` scales residues in aa."""
    from lib_palette import AXIS as AXIS_HEX, UNION as UNION_HEX
    from fig3bc_phospho_atlas import NONE_FILL, STEPS as CELL
    left = [t for t in steps if t[3] > 0 and t[6][-1] > 1e-9]
    right = [t for t in steps if t[4] > 0 and t[7][-1] > 1e-9]
    left.sort(key=lambda t: (-t[3], -t[6][-1], t[0], t[1]) if lengths
              else (-t[6][-1], t[0], t[1]))
    right.sort(key=lambda t: (-t[7][-1], t[0], t[1]))
    covered = {t[1] for t in left} | {t[1] for t in right}
    gain = len({t[1] for rows, cum in ((left, 6), (right, 7)) for t in rows
                if t[cum][-1] > t[cum][0] + 1e-9})
    n_tot = n_total if n_total is not None else len(steps)
    colour = assign(list(ORDER) + ["All"])
    cols = ["GluC", "LysC", "Trypsin", "All"]
    fs = 1.45

    def wide(txt, size):
        return 0.56 * size * fs * TEXT_BOOST * len(txt)
    cw, gap_all = 58.0, 8.0                   # a protease column; All set apart

    def col_x(i):
        return g0 + i * cw + (gap_all if cols[i] == "All" else 0)
    n_j = sum(len(t[11]) for t in right)
    lx = 26.0
    gene_w = max(1.15 * wide(t[0], 9.4) for t in left + right)
    acc_w = max(wide(t[1], 8.4) for t in left + right)
    dl_w = max(wide(f"Δ{a_}–{b_}", 8.0) for t in right for a_, b_, _n in t[11])
    # one set of protease columns for both blocks, clear of the widest labels
    g0 = lx + gene_w + 8 + acc_w + 10 + dl_w + 12
    W = round(col_x(3) + cw + 20, 1)

    rh, jh, iso_gap = 18.0, 15.0, 3.0
    mt = 96.0                                  # residue block's first row
    jt = mt + rh * len(left) + 90              # + the bars' axis under them
    k_top = jt + jh * n_j + iso_gap * (len(right) - 1) + 24   # the count axis
    H = round(k_top, 1)

    c = Canvas(W, H, font, font_scale=fs, out_w=width * W / 1010.0)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    # clear of the letter lib_compose.py stamps
    c.text(48, 30, title or "PD-implicated isoforms", 11.5, INK, "start", "600")

    def headers(y_top, head):
        c.text(lx, y_top - 34, head, 9.8, INK, "start", "600")
        for i, d in enumerate(cols):
            c.text(col_x(i) + cw / 2, y_top - 7, display(d) if d != "All" else "All",
                   9.0, INK, "middle", "600" if d == "All" else None)

    # Residues: trypsin from 0, then what Lys-C and Glu-C add.
    c.text(lx, mt - 34, f"Discriminating residues · {len(left)} isoforms", 9.8, INK,
           "start", "600")
    n_aa = max(t[3] for t in left) if lengths else 1
    end_w = max(wide(f"{round(t[6][-1] * t[3])}/{t[3]} aa", 8.4) for t in left) \
        if lengths else wide("100%", 8.4)
    bx0 = lx + gene_w + 8 + acc_w + 12
    bx1 = W - 20 - end_w - 8

    def BX(v):
        return bx0 + v / n_aa * (bx1 - bx0)
    seg = [colour["Trypsin"], colour["LysC"], colour["GluC"]]
    if lengths:
        step = next(s_ for s_ in (1, 2, 5, 10, 20, 25, 50, 100, 200, 500)
                    if n_aa / s_ <= 4)
        ticks = [(q, str(q)) for q in range(0, n_aa + 1, step)]
    else:
        ticks = [(q, f"{100 * q:.0f}%") for q in (0.0, 0.5, 1.0)]
    for q, lab in ticks:
        c.line(BX(q), mt - 4, BX(q), mt + rh * len(left), stroke=GRID, sw=1.0)
        c.text(BX(q), mt + rh * len(left) + 14, lab, 8.0, INK_SECONDARY, "middle")
    if lengths:
        c.text((bx0 + bx1) / 2, mt + rh * len(left) + 34,
               "discriminating residues (aa); grey = not covered", 8.6,
               INK_SECONDARY, "middle")
    for r, t in enumerate(left):
        y = mt + r * rh + rh / 2
        c.text(lx, y + 4, t[0], 9.4, INK, "start", "600")
        c.text(lx + gene_w + 8, y + 4, t[1], 8.4, INK_SECONDARY, "start")
        scale = t[3] if lengths else 1
        edges = [0.0] + [v * scale for v in t[6]]    # trypsin, + Lys-C, + Glu-C
        tip = " → ".join(f"{round(v * t[3])} aa" if lengths else f"{100 * v:.0f}%"
                         for v in t[6])
        c.add(f'<g><title>{t[0]} {t[1]}: trypsin, + Lys-C, + Glu-C: {tip}'
              + (f' of {t[3]} aa' if lengths else '') + '</title>')
        if lengths:
            c.rect(BX(0), y - 5, BX(t[3]) - BX(0), 10, AXIS_HEX)
        for k in range(3):
            if edges[k + 1] > edges[k] + 1e-9:
                c.rect(BX(edges[k]), y - 5, BX(edges[k + 1]) - BX(edges[k]), 10,
                       seg[k])
        c.add('</g>')
        lab = f"{round(t[6][-1] * t[3])}/{t[3]} aa" if lengths \
            else f"{100 * t[6][-1]:.0f}%"
        c.text(BX(scale) + 6 if lengths else BX(t[6][-1]) + 6, y + 4, lab, 8.4,
               INK, "start")

    # Junctions: spanning peptides from trypsin, then Lys-C, then Glu-C.
    c.text(lx, jt - 34, f"Deletion junctions · {n_j} in {len(right)} isoforms",
           9.8, INK, "start", "600")
    order = ["Trypsin", "LysC", "GluC"]
    n_max = max(sum(n[d] for d in order) for t in right for _a, _b, n in t[11])
    jb0 = g0
    jb1 = W - 20 - wide(str(n_max), 8.4) - 8

    def JX(v):
        return jb0 + v / max(n_max, 1) * (jb1 - jb0)
    y_end = jt + jh * n_j + iso_gap * (len(right) - 1)
    for q in range(n_max + 1):
        c.line(JX(q), jt - 4, JX(q), y_end, stroke=GRID, sw=1.0)
        c.text(JX(q), y_end + 14, str(q), 8.0, INK_SECONDARY, "middle")
    y = jt
    for t in right:
        y0 = y
        for ji, (a_, b_, n) in enumerate(t[11]):
            cy = y + jh / 2
            if ji == 0:
                c.text(lx, cy + 4, t[0], 9.4, INK, "start", "600")
                c.text(lx + gene_w + 8, cy + 4, t[1], 8.4, INK_SECONDARY, "start")
            c.text(g0 - 12, cy + 4, f"Δ{a_}–{b_}", 8.0, INK_SECONDARY, "end")
            total = sum(n[d] for d in order)
            c.add(f'<g><title>{t[0]} {t[1]} Δ{a_}–{b_}: ' + ", ".join(
                f"{display(d)} {n[d]}" for d in order) + '</title>')
            x = 0
            for d in order:
                if n[d]:
                    c.rect(JX(x), cy - 4.5, JX(x + n[d]) - JX(x), 9, colour[d])
                    x += n[d]
            c.add('</g>')
            c.text(JX(total) + 6, cy + 4, str(total), 8.4,
                   INK if total else INK_SECONDARY, "start")
            y += jh
        if len(t[11]) > 1:                     # a bracket down the isoform's rows
            c.line(lx - 8, y0 + 4, lx - 8, y - 4, INK_MUTED, 1.0)
        y += iso_gap

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({W:.0f} x {H:.0f} units; {len(covered)} of {n_tot} "
          f"isoforms covered: {len(left)} on residues; {n_j} junctions in "
          f"{len(right)} isoforms)")


PD_GENES = ("SNCA LRRK2 PINK1 PRKN UBB UBC UBA52 RPS27A RAB10 RAB8A RAB12 RAB29 "
            "RAB35 TH GSK3B PLK2 CSNK2A1 CSNK2A2 CSNK2B PARK7 SYNJ1 EIF4EBP1 "
            "MSN EZR RDX VPS35 MAPT GBA1 SNCB SNCG DNAJC6 VPS13C ATP13A2 FBXO7 "
            "DDC SLC6A3 ALDH1A1 CHCHD2 GCH1 UCHL1 LAMP2")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", nargs="?", default="panel",
                    choices=("panel", "scan", "pd"))
    ap.add_argument("--genes", default=PD_GENES,
                    help="pd: gene symbols, space- or comma-separated")
    ap.add_argument("--title", default=None)
    ap.add_argument("--grid-dots", action="store_true",
                    help="pd: residue heatmap over a junction dot matrix")
    ap.add_argument("--lengths", action="store_true",
                    help="pd --grid-dots: residue bars in aa over each region's length")
    ap.add_argument("--grid", action="store_true",
                    help="pd: isoform x protease grid instead of bars")
    ap.add_argument("reports", nargs="*", default=["data/search/*-60min-Phospho.parquet"])
    ap.add_argument("--fasta", default="data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta")
    ap.add_argument("--tsv", default="derived/da_iso/isoform_disc_coverage.tsv")
    ap.add_argument("--precursor-q", type=float, default=0.01)
    ap.add_argument("--letter", default="b")
    ap.add_argument("--out", default="figures/fig2c_disccov.svg")
    ap.add_argument("--font", default=FONT)
    ap.add_argument("--width", type=float, default=1215.0)
    ap.add_argument("--junctions", action="store_true",
                    help="also score internal deletion junctions")
    args = ap.parse_args(argv)

    if args.mode == "scan":
        paths = sorted(p for pat in args.reports for p in glob.glob(pat))
        if not paths:
            sys.exit("no reports matched")
        scan(paths, args.fasta, args.tsv, args.precursor_q, args.junctions)
        return 0
    if args.mode == "pd":
        paths = sorted(p for pat in args.reports for p in glob.glob(pat))
        if not paths:
            sys.exit("no reports matched")
        want = set(args.genes.replace(",", " ").split())
        steps = pd_steps(paths, args.fasta, want, args.precursor_q)
        draw = pd_grid_dots if args.grid_dots else pd_grid if args.grid \
            else pd_path
        kw = {"lengths": True} if args.lengths and args.grid_dots else {}
        draw(steps, args.out, args.font, args.letter, args.width, args.title, **kw)
        return 0
    panel(load(args.tsv), args.out, args.font, args.letter, args.width,
          with_junctions=args.junctions)
    return 0


if __name__ == "__main__":
    sys.exit(main())
