#!/usr/bin/env python3
"""Coverage of the residues and deletion junctions that discriminate each isoform.

    python3 fig2cd_isoform_coverage.py scan 'data/search/*-60min-Phospho.parquet' --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta --tsv derived/da_iso/isoform_disc_coverage.tsv
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
from lib_palette import (AXIS, FONT, GRID, INK, INK_MUTED, INK_SECONDARY,  # noqa: E402
                         TEXT_BOOST, assign, display)

ORDER = ["GluC", "LysC", "Trypsin"]
UNION = "union"
FIELDS = (["isoform", "n_disc"] + [f"cov_{d}" for d in ORDER] + ["cov_union"]
          + ["n_junc"] + [f"jcov_{d}" for d in ORDER] + ["jcov_union"])


def peptides_by_base(paths):
    """-> {canonical accession: {digest: set(peptides)}}, one gene's groups."""
    import fig2a_isoform_strip as ie
    groups = ie.collect_groups(sorted(paths), 0.01)
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
    import fig2a_isoform_strip as ie
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
    import fig2a_isoform_strip as ie
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


def scan(paths, fasta, out):
    """Per isoform: discriminating residues and junctions, and how many each protease covers."""
    import fig2a_isoform_strip as ie
    from lib_fasta import read_fasta
    seqs = read_fasta(fasta)
    bybase = peptides_by_base(paths)

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
            joins = junctions(cseq, iseq)
            if not disc and not joins:
                deletions += 1
                continue
            cov = region_coverage(iseq, disc, byd)
            jcov = junction_coverage(iseq, joins, byd) if joins else \
                {d: set() for d in ORDER}
            rows.append({"isoform": iso, "n_disc": len(disc),
                         **{f"cov_{d}": len(cov[d]) for d in ORDER},
                         "cov_union": len(set().union(*cov.values())),
                         "n_junc": len(joins),
                         **{f"jcov_{d}": len(jcov[d]) for d in ORDER},
                         "jcov_union": len(set().union(*jcov.values()))})

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {out}  ({len(rows):,} isoforms with a discriminating region or "
          f"junction; {deletions:,} terminal-only deletions excluded)")
    return rows


def load(path):
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for k in FIELDS[1:]:
                r[k] = int(r[k])
            rows.append(r)
    return rows


def curves(rows):
    """-> (seen, {series: [isoforms reaching >= t%]}) for t = 1..100."""
    seen = [r for r in rows if r["n_disc"] > 0 and r[f"cov_{UNION}"] > 0]
    out = {}
    for d in ORDER + [UNION]:
        fr = [r[f"cov_{d}"] / r["n_disc"] for r in seen]
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


def panel(rows, out, letter=""):
    """Isoforms against how much of their discriminating region is covered, beside junction bars."""
    seen, cur = curves(rows)
    jn = {d: sum(1 for r in rows if r["n_junc"] and r[f"jcov_{d}"])
          for d in ORDER + [UNION]}
    colour = assign(list(ORDER) + ["All"])
    colour[UNION] = colour["All"]

    W = 1010.0
    ml = 90.0
    mt = 112.0
    ph = 208.0
    pw = W - ml - 300.0
    bx0, bx1 = ml + pw + 40.0, W - 10.0
    H = mt + ph + 52
    c = Canvas(W, H, FONT, font_scale=1.45, out_w=1215.0)

    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")

    top = (int(max(cur[UNION] + [jn[UNION]]) / 250) + 1) * 250

    def X(t):
        return ml + (t - 1) / 99.0 * pw

    def Y(v):
        return mt + ph - v / top * ph

    for t in range(0, top + 1, 250):
        c.line(ml, Y(t), ml + pw, Y(t), stroke=GRID, sw=1)
        c.line(bx0, Y(t), bx1, Y(t), stroke=GRID, sw=1)
        c.text(ml - 8, Y(t) + 3, f"{t:,}", 8.4, INK_MUTED, "end")
    for t in (1, 25, 50, 75, 100):
        c.text(X(t), mt + ph + 17, f"{t}%", 8.8, INK_MUTED, "middle")
    c.text(26, mt + ph / 2, "non-canonical isoforms", 9.6, INK_SECONDARY, "middle", rot=-90)
    c.text(ml + pw / 2, mt + ph + 36,
           "of the discriminating region covered",
           9.2, INK_SECONDARY, "middle")

    for key in ORDER + [UNION]:
        pts = " ".join(f"{X(t + 1):.1f},{Y(v):.1f}"
                       for t, v in enumerate(cur[key]))
        c.add(f'<polyline points="{pts}" fill="none" stroke="{colour[key]}" '
              f'stroke-width="{2.6 if key == UNION else 1.8}" '
              f'stroke-linejoin="round"/>')

    # junction bars share the curves' y axis
    slot = (bx1 - bx0) / (len(ORDER) + 1)
    sw_ = 170.0
    lab_w = 0.55 * 7.8 * c.fs * TEXT_BOOST * len("canonical") + 6
    for cx_, title, kind in (
            (ml + pw / 2, "Discriminating-region peptides", "own"),
            ((bx0 + bx1) / 2, "Junction-spanning peptides", "junction")):
        c.text(cx_, 38, title, 9.2, INK, "middle", "600")
        schematic(c, cx_ - (sw_ - lab_w) / 2, 56, sw_, kind, colour[UNION])
    c.text(ml + pw - 6, mt + 16, f"n = {len(seen):,}", 9.2, INK, "end")
    for i, key in enumerate(ORDER + [UNION]):
        bx, v = bx0 + slot * i + slot * 0.18, jn[key]
        c.rect(bx, Y(v), slot * 0.64, Y(0) - Y(v), fill=colour[key],
               fo=0.9, rx=1.5)
        c.text(bx + slot * 0.32, Y(v) - 5, f"{v:,}", 8.4, INK, "middle")
        c.text(bx + slot * 0.32, mt + ph + 17,
               "All" if key == UNION else display(key), 8.4, INK, "middle")
    c.line(bx0, mt + ph, bx1, mt + ph, stroke=AXIS, sw=1)

    c.line(ml, mt + ph, ml + pw, mt + ph, stroke=AXIS, sw=1)

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({W:.0f}x{H:.0f})")
    tot = sum(r["n_disc"] for r in seen)
    for k in ORDER + [UNION]:
        cv = sum(r[f"cov_{k}"] for r in seen)
        med = statistics.median(r[f"cov_{k}"] / r["n_disc"] for r in seen)
        print(f"    {display(k) if k in ORDER else 'All':<10s} "
              f"{100 * cv / tot:>5.1f}% of all discriminating residues, "
              f"median per isoform {med:>5.1%}, "
              f"≥50% in {cur[k][49]:>5,} isoforms")
    for k in ORDER + [UNION]:
        print(f"    {display(k) if k in ORDER else 'All':<10s} a junction "
              f"spanned in {jn[k]:>5,} isoforms")


STEPS = [("Trypsin",), ("Trypsin", "LysC"), ("Trypsin", "LysC", "GluC")]


def pd_steps(paths, fasta, want):
    """-> [(gene, isoform, n_res, n_junc, res frac per STEPS, junc frac per STEPS, [(start, end, {digest: n})])]."""
    import fig2a_isoform_strip as ie
    from lib_fasta import gene_map
    from lib_fasta import read_fasta
    seqs = read_fasta(fasta)
    gm = gene_map(fasta)
    gene = lambda a: gm.get(a) or gm.get(a.split("-")[0])
    bybase = peptides_by_base(paths)
    out = []
    for iso in sorted(a for a in seqs if ie.ISO.search(a) and gene(a) in want):
        base = ie.ISO.sub("", iso)
        if base not in seqs or seqs[iso] == seqs[base]:
            continue
        cseq, iseq = seqs[base], seqs[iso]
        disc = discriminating(cseq, iseq)
        joins = junctions(cseq, iseq)
        if not disc and not joins:
            continue
        byd = bybase.get(base, {})
        rc = region_coverage(iseq, disc, byd)
        jc = junction_coverage(iseq, joins, byd)
        rfr, jfr = [], []
        for step in STEPS:
            got = set().union(*(rc[d] for d in step))
            jg = set().union(*(jc[d] for d in step))
            rfr.append(len(got) / len(disc) if disc else 0.0)
            jfr.append(len(jg) / len(joins) if joins else 0.0)
        cut = {i0: (c0 + 1, c1) for op, c0, c1, i0, i1 in ie.events(cseq, iseq)
               if op == "delete" and 0 < i0 < len(iseq)}
        per_j = []
        for j in joins:
            pj = junction_peptides(iseq, [j], byd)
            per_j.append(cut[j] + ({d: len(pj[d]) for d in ORDER},))
        out.append((gene(iso), iso, len(disc), len(joins), rfr, jfr, per_j))
    return out


def pd_grid_dots(steps, out, letter=""):
    """Residue bars in aa over each region's length, junction bars under them."""
    left = [t for t in steps if t[2] > 0 and t[4][-1] > 1e-9]
    right = [t for t in steps if t[3] > 0 and t[5][-1] > 1e-9]
    left.sort(key=lambda t: (-t[2], -t[4][-1], t[0], t[1]))
    right.sort(key=lambda t: (-t[5][-1], t[0], t[1]))
    covered = {t[1] for t in left} | {t[1] for t in right}
    colour = assign(list(ORDER) + ["All"])
    cols = ["GluC", "LysC", "Trypsin", "All"]
    fs = 1.45

    def wide(txt, size):
        return 0.56 * size * fs * TEXT_BOOST * len(txt)
    cw, gap_all = 58.0, 8.0

    def col_x(i):
        return g0 + i * cw + (gap_all if cols[i] == "All" else 0)
    n_j = sum(len(t[6]) for t in right)
    lx = 26.0
    gene_w = max(1.15 * wide(t[0], 9.4) for t in left + right)
    acc_w = max(wide(t[1], 8.4) for t in left + right)
    dl_w = max(wide(f"Δ{a_}–{b_}", 8.0) for t in right for a_, b_, _n in t[6])
    g0 = lx + gene_w + 8 + acc_w + 10 + dl_w + 12
    W = round(col_x(3) + cw + 20, 1)

    rh, jh, iso_gap = 18.0, 15.0, 3.0
    mt = 96.0                                  # residue block's first row
    jt = mt + rh * len(left) + 90              # junction block's first row
    H = round(jt + jh * n_j + iso_gap * (len(right) - 1) + 24, 1)

    c = Canvas(W, H, FONT, font_scale=fs, out_w=1215.0 * W / 1010.0)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    c.text(48, 30, "PD-implicated isoforms", 11.5, INK, "start", "600")

    # Residues: trypsin from 0, then what Lys-C and Glu-C add.
    c.text(lx, mt - 34, f"Discriminating residues · {len(left)} isoforms", 9.8, INK,
           "start", "600")
    n_aa = max(t[2] for t in left)
    end_w = max(wide(f"{round(t[4][-1] * t[2])}/{t[2]} aa", 8.4) for t in left)
    bx0 = lx + gene_w + 8 + acc_w + 12
    bx1 = W - 20 - end_w - 8

    def BX(v):
        return bx0 + v / n_aa * (bx1 - bx0)
    seg = [colour["Trypsin"], colour["LysC"], colour["GluC"]]
    step = next(s_ for s_ in (1, 2, 5, 10, 20, 25, 50, 100, 200, 500)
                if n_aa / s_ <= 4)
    for q in range(0, n_aa + 1, step):
        c.line(BX(q), mt - 4, BX(q), mt + rh * len(left), stroke=GRID, sw=1.0)
        c.text(BX(q), mt + rh * len(left) + 14, str(q), 8.0, INK_SECONDARY,
               "middle")
    c.text((bx0 + bx1) / 2, mt + rh * len(left) + 34,
           "discriminating residues (aa); grey = not covered", 8.6,
           INK_SECONDARY, "middle")
    for r, t in enumerate(left):
        y = mt + r * rh + rh / 2
        c.text(lx, y + 4, t[0], 9.4, INK, "start", "600")
        c.text(lx + gene_w + 8, y + 4, t[1], 8.4, INK_SECONDARY, "start")
        edges = [0.0] + [v * t[2] for v in t[4]]     # trypsin, + Lys-C, + Glu-C
        tip = " → ".join(f"{round(v * t[2])} aa" for v in t[4])
        c.add(f'<g><title>{t[0]} {t[1]}: trypsin, + Lys-C, + Glu-C: {tip} '
              f'of {t[2]} aa</title>')
        c.rect(BX(0), y - 5, BX(t[2]) - BX(0), 10, AXIS)
        for k in range(3):
            if edges[k + 1] > edges[k] + 1e-9:
                c.rect(BX(edges[k]), y - 5, BX(edges[k + 1]) - BX(edges[k]), 10,
                       seg[k])
        c.add('</g>')
        c.text(BX(t[2]) + 6, y + 4, f"{round(t[4][-1] * t[2])}/{t[2]} aa", 8.4,
               INK, "start")

    # Junctions: spanning peptides from trypsin, then Lys-C, then Glu-C.
    c.text(lx, jt - 34, f"Deletion junctions · {n_j} in {len(right)} isoforms",
           9.8, INK, "start", "600")
    order = ["Trypsin", "LysC", "GluC"]
    n_max = max(sum(n[d] for d in order) for t in right for _a, _b, n in t[6])
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
        for ji, (a_, b_, n) in enumerate(t[6]):
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
        if len(t[6]) > 1:                      # a bracket down the isoform's rows
            c.line(lx - 8, y0 + 4, lx - 8, y - 4, INK_MUTED, 1.0)
        y += iso_gap

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({W:.0f} x {H:.0f} units; {len(covered)} of {len(steps)} "
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
                    choices=("panel", "scan", "pd"),
                    help="panel: curves + junction bars from --tsv; scan: write --tsv; "
                         "pd: the PD-gene isoforms")
    ap.add_argument("reports", nargs="*", default=["data/search/*-60min-Phospho.parquet"])
    ap.add_argument("--fasta", default="data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta")
    ap.add_argument("--tsv", default="derived/da_iso/isoform_disc_coverage.tsv")
    ap.add_argument("--letter", default="")
    ap.add_argument("--out", default="coverage.svg")
    args = ap.parse_args(argv)

    if args.mode == "panel":
        panel(load(args.tsv), args.out, args.letter)
        return 0
    paths = sorted(p for pat in args.reports for p in glob.glob(pat))
    if not paths:
        sys.exit("no reports matched")
    if args.mode == "scan":
        scan(paths, args.fasta, args.tsv)
    else:
        pd_grid_dots(pd_steps(paths, args.fasta, set(PD_GENES.split())), args.out,
                     args.letter)
    return 0


if __name__ == "__main__":
    sys.exit(main())
