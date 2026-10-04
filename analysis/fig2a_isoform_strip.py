#!/usr/bin/env python3
"""Canonical and isoform on one frame, with the peptides only one of them contains.

    python3 analysis/fig2a_isoform_strip.py 'data/search/*-60min-Phospho.parquet' --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta --protein P27816 --isoform P27816-3 --label MAP4 --out strip.svg
"""

import argparse
import collections
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_svg import Canvas                                  # noqa: E402
from lib_palette import (FONT, INK, INK_MUTED, INK_SECONDARY, TEXT_BOOST,  # noqa: E402
                         assign)

ORDER = ["GluC", "LysC", "Trypsin"]
ISO = re.compile(r"-\d+$")

# Residues of agreement needed to keep two differing blocks apart.
MIN_EQUAL = 3


def align(c, i):
    """-> (common-prefix length, non-equal difflib opcodes of the trimmed cores)."""
    p = 0
    while p < min(len(c), len(i)) and c[p] == i[p]:
        p += 1
    s = 0
    while s < min(len(c), len(i)) - p and c[len(c) - 1 - s] == i[len(i) - 1 - s]:
        s += 1
    cc, ii = c[p:len(c) - s], i[p:len(i) - s]
    if not cc and not ii:
        return p, []
    if not ii:
        return p, [("delete", 0, len(cc), 0, 0)]
    if not cc:
        return p, [("insert", 0, 0, 0, len(ii))]
    import difflib
    ops = [o for o in difflib.SequenceMatcher(None, cc, ii, autojunk=False)
           .get_opcodes() if o[0] != "equal"]
    return p, ops


def merge(p, ops, min_equal=MIN_EQUAL):
    """Pool differing blocks closer than `min_equal`, shifted back by the prefix `p`."""
    merged = []
    for tag, a1, a2, b1, b2 in ops:
        if merged and a1 - merged[-1][2] < min_equal:
            merged[-1][2], merged[-1][4] = a2, b2
            merged[-1][0] = "replace"
        else:
            merged.append([tag, a1, a2, b1, b2])
    out = []
    for tag, a1, a2, b1, b2 in merged:
        if a2 > a1 and b2 > b1:
            tag = "replace"
        elif a2 > a1:
            tag = "delete"
        else:
            tag = "insert"
        out.append((tag, p + a1, p + a2, p + b1, p + b2))
    return out


def events(c, i, min_equal=MIN_EQUAL):
    """Discrete differences between canonical and isoform."""
    return merge(*align(c, i), min_equal=min_equal)


def collect_groups(paths, precursor_q, bases=None, protein_q=0.01):
    """{group: {digest: set(peptides)}} in one pass, optionally filtered."""
    out = collections.defaultdict(lambda: collections.defaultdict(set))
    import lib_report as rp
    for r in rp.open_reports(paths, order=ORDER):
        path, dig = r.path, r.protease
        F = ["protein_groups", "peptide", "precursor_q"]
        has_pg = protein_q is not None and r.has("protein_q_run")
        if has_pg:
            F.append("protein_q_run")
        for b in r.batches(F):
            g = b["protein_groups"]
            s = b["peptide"]
            q = b["precursor_q"]
            pgq = b["protein_q_run"] if has_pg else None
            for k in range(b["_n"]):
                if not g[k] or not s[k] or q[k] is None or q[k] > precursor_q:
                    continue
                if pgq is not None and (pgq[k] is None or pgq[k] > protein_q):
                    continue
                if bases is not None:
                    accs = [a.strip() for a in g[k].split(";") if a.strip()]
                    if not ({ISO.sub("", a) for a in accs} & bases):
                        continue
                out[g[k]][dig].add(s[k])
        print(f"  scanned {os.path.basename(path)[:46]:<46s} {dig}")
    return out


def pack_hits(hits, mapper=lambda j: j):
    """Peptides -> ([(start, end, digest, row)], n_rows), packed by overlap."""
    out, ends = [], []
    for pos, pep, d in sorted(hits, key=lambda t: mapper(t[0])):
        a = mapper(pos)
        b = mapper(pos + len(pep) - 1) + 1
        for t, last in enumerate(ends):
            if a >= last:
                ends[t] = b
                break
        else:
            t = len(ends)
            ends.append(b)
        out.append((a, b, d, t))
    return out, len(ends)


def aligned_panel(base, iso_acc, cseq, iseq, can_peps, iso_peps, out, label=None,
                  min_equal=MIN_EQUAL, margin=146.0, bar_height=22.0):
    """Both forms in one frame, each change laid out canonical-then-isoform; absent = thin line."""
    colour = assign(ORDER)
    ev = events(cseq, iseq, min_equal)
    if not ev:
        sys.exit(f"{base} and {iso_acc} have identical sequences")
    cmap, imap, segs, u, c0, i0 = {}, {}, [], 0, 0, 0
    for _t, cs, ce, is_, ie in ev + [("", len(cseq), len(cseq), len(iseq), len(iseq))]:
        if cs > c0:
            segs.append(("shared", u, u + cs - c0))
            for r in range(cs - c0):
                cmap[c0 + r] = imap[i0 + r] = u + r
            u += cs - c0
        if ce > cs or ie > is_:
            for r in range(ce - cs):
                cmap[cs + r] = u + r
            for r in range(ie - is_):
                imap[is_ + r] = u + (ce - cs) + r
            segs.append(("change", u, u + (ce - cs) + (ie - is_), cs, ce, is_, ie))
            u += (ce - cs) + (ie - is_)
        c0, i0 = ce, ie
    W, bar_h = 1010.0, bar_height
    ml, mr = margin, 30.0
    k = (W - ml - mr) / u

    def X(v):
        return ml + v * k

    def runs(positions):
        """Union positions -> contiguous [a, b) runs."""
        out_ = []
        for v in positions:
            if out_ and v == out_[-1][1]:
                out_[-1][1] = v + 1
            else:
                out_.append([v, v + 1])
        return out_

    def hits(seq, other, peps, m):
        found = []
        for dig in ORDER:
            for pep in sorted(peps.get(dig, ())):
                if pep in other:
                    continue
                i = seq.find(pep)
                if i != -1:
                    found.append((i, pep, dig))
        found.sort(key=lambda t: (m[t[0]], len(t[1]), t[1], t[2]))
        rows, n = pack_hits(found, lambda j: m[min(j, len(seq) - 1)])
        return [(rw, runs([m[j] for j in range(i, i + len(pep))]), dig)
                for (i, pep, dig), rw in zip(sorted(found, key=lambda t: m[t[0]]),
                                             [r[3] for r in rows])], n

    can_rows, can_n = hits(cseq, iseq, can_peps, cmap)
    iso_rows, iso_n = hits(iseq, cseq, iso_peps, imap)
    pep_h, pitch = 4.0, 6.0
    yc = 62.0
    yi = yc + bar_h + 12 + max(can_n, 1) * pitch + 14
    H = yi + bar_h + 12 + max(iso_n, 1) * pitch + 14
    c = Canvas(W, H, FONT, font_scale=1.4, out_w=1215.0)
    c.text(58.0, 28, label or base, 12, INK, "start", "600")
    for name, seq, y, m, rows in ((base, cseq, yc, cmap, can_rows),
                                  (iso_acc, iseq, yi, imap, iso_rows)):
        c.text(ml - 10, y + 9, name, 9.5, INK, "end", "600")
        c.text(ml - 10, y + 9 + 9.5 * c.fs * TEXT_BOOST * 0.95, f"{len(seq)} aa", 8, INK_MUTED, "end")
        own = runs([m[r] for r in range(len(seq))])
        c.rect(X(own[0][0]), y + bar_h / 2 - 1.0, (own[-1][1] - own[0][0]) * k, 2.0,
               fill=INK, stroke="none", rx=0)
        for a_, b_ in own:
            c.rect(X(a_), y, (b_ - a_) * k, bar_h, fill=INK, stroke="none", rx=0)
        for rw, rr, dig in rows:
            yy = y + bar_h + 6 + rw * pitch
            if len(rr) > 1:
                c.line(X(rr[0][0]), yy + pep_h / 2, X(rr[-1][1]), yy + pep_h / 2,
                       stroke=colour[dig], sw=0.8, so=0.7)
            for a_, b_ in rr:
                c.rect(X(a_), yy, max((b_ - a_) * k, 1.5), pep_h, fill=colour[dig],
                       stroke="none", rx=1)
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({len(ev)} changes; {len(can_rows)} canonical-only, "
          f"{len(iso_rows)} isoform-only peptides)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reports", nargs="+")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--protein", required=True, help="canonical accession")
    ap.add_argument("--isoform", required=True)
    ap.add_argument("--label", default=None)
    ap.add_argument("--min-equal", dest="min_equal", type=int, default=MIN_EQUAL,
                    help="shared residues needed to keep two differing blocks apart")
    ap.add_argument("--margin", type=float, default=146.0,
                    help="left margin before the bars")
    ap.add_argument("--bar-height", dest="bar_height", type=float, default=22.0,
                    help="bar thickness")
    ap.add_argument("--out", default="strip.svg")
    args = ap.parse_args(argv)

    paths = []
    for pat in args.reports:
        paths.extend(sorted(glob.glob(pat)) or [pat])
    from lib_fasta import read_fasta
    seqs = read_fasta(args.fasta)
    print(f"{len(seqs):,} sequences in {os.path.basename(args.fasta)}")
    for a in (args.protein, args.isoform):
        if a not in seqs:
            sys.exit(f"{a} not in {args.fasta}")

    # every peptide of the gene, placed by sequence rather than by reported group
    pooled = collections.defaultdict(set)
    for byd in collect_groups(paths, 0.01, bases={args.protein}).values():
        for d, v in byd.items():
            pooled[d] |= v
    aligned_panel(args.protein, args.isoform, seqs[args.protein], seqs[args.isoform],
                  pooled, pooled, args.out, args.label, args.min_equal, args.margin,
                  args.bar_height)
    return 0


if __name__ == "__main__":
    sys.exit(main())
