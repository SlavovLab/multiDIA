#!/usr/bin/env python3
"""Supplementary Fig. 3: tau's exon-skipping proteoforms on big tau, each with the tau-unique peptides it contains; big tau's peptides coloured by log2 PD / Control."""

import argparse
import glob
import itertools
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fig2c_supp2_donor_levels as dl                                # noqa: E402
import lib_isoform as ie                                             # noqa: E402
from lib_fasta import gene_map, read_fasta                           # noqa: E402
from lib_palette import (FIG_W, FONT, INK, INK_SECONDARY, STRIP_FILL, TEXT_BOOST,  # noqa: E402
                         assign, display, pt_scale)
from lib_report import ORDER, open_reports, read_metadata            # noqa: E402
from lib_svg import Canvas, text_width                               # noqa: E402

W = FIG_W
FS = pt_scale(W, W)
U = FS * TEXT_BOOST                     # units per printed pt
SHOW = ["Trypsin", "LysC", "GluC"]
COL = assign(SHOW)


def parse_exons(specs):
    """['2:45-73', ...] -> {'2': (45, 73)}, 1-based residues of --protein."""
    out = {}
    for s in specs:
        name, _, span = s.rpartition(":")
        a, b = (int(x) for x in span.split("-"))
        out[name] = (a, b)
    return out


def forms(seqs, protein, exons):
    """-> {entry: the --protein residues it keeps} for entries that are --protein minus some exons."""
    canon = seqs[protein]
    built = {}
    for k in range(len(exons) + 1):
        for skip in itertools.combinations(sorted(exons), k):
            gone = {p for n in skip for p in range(exons[n][0], exons[n][1] + 1)}
            keep = [p for p in range(1, len(canon) + 1) if p not in gone]
            built.setdefault("".join(canon[p - 1] for p in keep), keep)
    family = sorted(a for a in seqs if ie.ISO.sub("", a) == protein)
    return {a: built[seqs[a]] for a in family if seqs[a] in built}


def detected(paths, protein, q_cut=0.01):
    """-> {(digest, peptide)} for --protein's groups at precursor and run-wise protein q <= q_cut."""
    seen = set()
    for r in open_reports(paths, order=ORDER):
        F = ["protein_groups", "peptide", "precursor_q"]
        pg = r.has("protein_q_run")
        if pg:
            F.append("protein_q_run")
        for b in r.batches(F):
            for k in range(b["_n"]):
                g, s, q = b["protein_groups"][k], b["peptide"][k], b["precursor_q"][k]
                if not g or not s or q is None or q > q_cut:
                    continue
                if pg and (b["protein_q_run"][k] is None or b["protein_q_run"][k] > q_cut):
                    continue
                if protein in {ie.ISO.sub("", a.strip()) for a in g.split(";")}:
                    seen.add((r.protease, s))
    return seen


def tau_name(keep, exons):
    """'2N4R' style name from the exons a proteoform keeps; big tau when it keeps exon 4a."""
    has = {n: all(p in keep for p in range(a, b + 1)) for n, (a, b) in exons.items()}
    if has["4a"]:
        return "Big tau"
    return f"{int(has['2']) + int(has['3'])}N{'4' if has['10'] else '3'}R"


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


def tw(s, size, bold=False):
    return text_width(s, size * U, bold)


def pd_ratio(gdata, runs, cond):
    """-> {(protease, peptide): mean log2 PD / Control over PD donors}, each peptide against its own control mean (Fig. 2c)."""
    low = dl.low_runs(runs)
    out = {}
    for (d, p), v in gdata.items():
        v = {s: x for s, x in v.items() if (s, d) not in low}
        if sum(cond.get(s) == "Control" for s in v) < 3:
            continue
        rel = [x for s, x in dl.relative(v, cond).items() if cond.get(s) == "LBD"]
        if len(rel) >= 2:
            out[(d, p)] = sum(rel) / len(rel)
    return out


def draw(seqs, protein, exons, peps, ratio, out):
    """Proteoforms on the protein's sequence, skipped exons as thin lines, protease coverage lanes over big tau."""
    found = forms(seqs, protein, exons)
    n = len(seqs[protein])
    rows = sorted(((tau_name(set(keep), exons), acc, keep) for acc, keep in found.items()),
                  key=lambda r: (r[0] != "Big tau", r[0][2:], r[0][:2]))
    top, lab_w = 30.0, 26 + tw("Big tau", 10, True) + 6 * U + tw("P10636-8", 10) + 8 * U
    x0, x1 = lab_w, W - 10

    def X(pos):
        return x0 + (pos - 1) / n * (x1 - x0)
    bar_h, pep_h, pitch, gap, lane = 8 * U, 3.0 * U, 4.2 * U, 7 * U, 4.3 * U
    laid = []
    for name, acc, keep in rows:
        seq = seqs[acc]
        hits = []
        for d, pep in peps:
            i = seq.find(pep)
            if i >= 0:
                hits.append((keep[i], keep[i + len(pep) - 1], d, pep, [keep[j] for j in range(i, i + len(pep))]))
        hits.sort(key=lambda h: (h[0], h[1]))
        ends, placed = [], []
        for a_, b_, d, pep, pos in hits:
            for t, last in enumerate(ends):
                if a_ > last + 1:
                    ends[t] = b_
                    break
            else:
                t = len(ends)
                ends.append(b_)
            placed.append((t, d, pep, pos))
        laid.append((name, acc, keep, placed, len(ends)))
    y_ex = top + 64 * U
    y_lanes = y_ex + 8 * U
    y = y_lanes + 3 * lane + 6 * U
    H = y + sum(bar_h + 3 * U + max(nr, 1) * pitch + gap for *_r, nr in laid) + 4 * U
    c = Canvas(W, H, FONT, font_scale=FS, out_w=W)
    c.text(26, top, "Peptides detected on MAPT proteoforms", 12, INK, "start", "600")
    x = 26
    for d in SHOW:
        c.rect(x, top + 9 * U, 9 * U, 9 * U, fill=COL[d], rx=1.5)
        c.text(x + 12 * U, top + 17 * U, display(d), 10, INK)
        x += 12 * U + tw(display(d), 10) + 10 * U
    key = "big tau: log2 PD / Control, per peptide"
    c.text(26, top + 37 * U, key, 10, INK)
    x = 26 + tw(key, 10) + 6 * U
    for i in range(40):
        c.rect(x + i * 2.75 * U, top + 28 * U, 2.75 * U + 0.3, 10 * U,
               fill=dl.shade(-dl.CLIP + (i + 0.5) * 2 * dl.CLIP / 40))
    for v in (-2, 0, 2):
        c.text(x + (v + dl.CLIP) / (2 * dl.CLIP) * 110 * U, top + 48 * U, f"{v:+d}" if v else "0", 10,
               INK_SECONDARY, "middle")
    x += 110 * U + 14 * U
    c.rect(x, top + 28 * U, 10 * U, 10 * U, fill=dl.NOT_QUANT)
    c.text(x + 13 * U, top + 37 * U, "not computed", 10, INK)
    for name_, (a_, b_) in exons.items():
        c.rect(X(a_), y_ex - 9 * U, X(b_ + 1) - X(a_), 11 * U, fill=STRIP_FILL)
        c.text((X(a_) + X(b_ + 1)) / 2, y_ex, f"exon {name_}" if X(b_ + 1) - X(a_) > tw("exon 4a", 10) + 4 * U
               else name_, 10, INK, "middle")
    big_hits = next(placed for name, _a, _k, placed, _n in laid if name == "Big tau")
    for k, d in enumerate(SHOW):
        for a_, b_ in runs_of({r for _t, dd, _p, pos in big_hits if dd == d for r in pos}, 1, n + 1):
            c.rect(X(a_), y_lanes + k * lane, X(b_) - X(a_), lane - 1.2 * U, fill=COL[d])
    for name, acc, keep, placed, nr in laid:
        c.text(26, y + bar_h - 1.5 * U, name, 10, INK, "start", "600")
        c.text(26 + tw("Big tau", 10, True) + 6 * U, y + bar_h - 1.5 * U, acc, 10, INK_SECONDARY)
        c.line(X(keep[0]), y + bar_h / 2, X(keep[-1] + 1), y + bar_h / 2, stroke=INK, sw=1.0)
        for a_, b_ in runs_of(set(keep), 1, n + 1):
            c.rect(X(a_), y, X(b_) - X(a_), bar_h, fill=INK)
        for t, d, pep, pos in placed:
            yy = y + bar_h + 3 * U + t * pitch
            if name == "Big tau":
                fill = dl.shade(ratio[(d, pep)]) if (d, pep) in ratio else dl.NOT_QUANT
            else:
                fill = COL[d]
            pieces = runs_of(set(pos), 1, n + 1)
            if len(pieces) > 1:
                c.line(X(pieces[0][1]), yy + pep_h / 2, X(pieces[-1][0]), yy + pep_h / 2, stroke=fill, sw=1.0)
            for a_, b_ in pieces:
                c.rect(X(a_), yy, max(X(b_) - X(a_), 1.5), pep_h, fill=fill, rx=1)
        y += bar_h + 3 * U + max(nr, 1) * pitch + gap
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({len(laid)} proteoforms, {len(peps)} peptides, {sum(k in ratio for k in peps)} with a ratio)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("reports", nargs="+")
    ap.add_argument("--fasta", default="data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta")
    ap.add_argument("--metadata", default="data/metadata.xlsx")
    ap.add_argument("--levels", default="derived/peptide_levels", help="fig2c_supp2_donor_levels.py scan's output")
    ap.add_argument("--protein", required=True, help="entry holding every exon, e.g. P10636")
    ap.add_argument("--exon", action="append", required=True, metavar="NAME:START-END",
                    help="an alternative exon, in --protein residues")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    paths = sorted(p for pat in args.reports for p in glob.glob(pat))
    if not paths:
        sys.exit("no reports matched")
    seqs = read_fasta(args.fasta)
    other = "\n".join(sq for acc, sq in seqs.items() if ie.ISO.sub("", acc) != args.protein)
    peps = sorted((d, p) for d, p in detected(paths, args.protein) if p not in other)
    data, runs = dl.load(args.levels)
    gene = gene_map(args.fasta)[args.protein]
    ratio = pd_ratio(data[gene], runs, read_metadata(args.metadata))
    draw(seqs, args.protein, parse_exons(args.exon), peps, ratio, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
