#!/usr/bin/env python3
"""Sequence coverage per protein group from Spectronaut reports: a rank curve or a beeswarm.

    python3 fig1cd_coverage.py 'data/search/*.parquet' --fasta search_db.fasta --rank --by-sample --sample-regex 'CF_([0-9]{4})' --out figures/fig1c.svg
"""

import argparse
import csv
import glob
import math
import os
import re
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_svg import Canvas, swarm, text_width                      # noqa: E402
from lib_fasta import read_fasta                                   # noqa: E402
from lib_palette import (AXIS, FONT, GRID, INK, INK_MUTED, assign,  # noqa: E402
                         display, TEXT_BOOST)
from prep_counts import (PRECURSOR_Q, PROTEIN_Q, detect_protease,   # noqa: E402
                         parse_float, read_rows, resolve_columns)
from fig1b_depth import fmt_compact as compact                     # noqa: E402

WIDTH, HEIGHT, TEXT_SCALE = 595.0, 490.0, 2.12


def collect(paths, sample_regex=None, split_samples=False):
    """-> {protein_group: {protease: set(peptides)}}, keyed first by sample if split_samples."""
    out = (defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
           if split_samples else defaultdict(lambda: defaultdict(set)))
    srx = re.compile(sample_regex) if sample_regex else None
    for path in paths:
        probe, _, fh0 = read_rows(path)
        if fh0 is not None:
            fh0.close()
        idx0, _ = resolve_columns(probe)
        wanted = {probe[i] for i in idx0.values()}
        header, reader, fh = read_rows(path, wanted)
        try:
            idx, used = resolve_columns(header)
            for f in ("run", "protein_group", "strip_seq"):
                if f not in idx:
                    sys.exit(f"{path}: no column resolved for {f}\n"
                             f"  resolved: {used}")
            # the filename may name several enzymes; then the run name decides
            protease = detect_protease(os.path.basename(path))
            if protease and protease.startswith("AMBIGUOUS"):
                protease = None
            n = 0
            for row in reader:
                if "eg_qvalue" in idx:
                    q = parse_float(row[idx["eg_qvalue"]])
                    if q is None or q > PRECURSOR_Q:
                        continue
                if "decoy" in idx and row[idx["decoy"]].strip().lower() == "true":
                    continue
                if "pg_qvalue" in idx:
                    q = parse_float(row[idx["pg_qvalue"]])
                    if q is None or q > PROTEIN_Q:
                        continue
                pep = row[idx["strip_seq"]].strip()
                grp = row[idx["protein_group"]].strip()
                if not (pep and grp):
                    continue
                run = row[idx["run"]]
                p = protease or detect_protease(run)
                if not p or p.startswith("AMBIGUOUS"):
                    p = "?"
                if split_samples:
                    m = srx.search(run) if srx else None
                    s = (m.group(1) if m.groups() else m.group(0)) if m else None
                    out[s][grp][p].add(pep)
                else:
                    out[grp][p].add(pep)
                n += 1
            print(f"  {os.path.basename(path)[:52]:<52s} {protease or '?':<8s} "
                  f"{n:>9,} rows kept")
        finally:
            if fh is not None:
                fh.close()
    return out


def mark(seq, peptides):
    """Residue mask of everything the peptides cover."""
    cov = bytearray(len(seq))
    found = 0
    for p in peptides:
        i = seq.find(p)
        if i == -1:
            continue
        found += 1
        ones = b"\x01" * len(p)
        while i != -1:
            cov[i:i + len(p)] = ones
            i = seq.find(p, i + 1)
    return cov, found


def representative(group, seqs):
    """The first accession in a group that the FASTA holds."""
    accs = [a.strip() for a in re.split(r"[;,]", group) if a.strip()]
    accs = [a.split("|")[1] if a.count("|") >= 2 else a for a in accs]
    return next((a for a in accs if a in seqs), None)


def nice_ticks(hi, target=6):
    """Round tick steps giving roughly `target` intervals across [0, hi]."""
    if hi <= 0:
        return [0]
    raw = hi / target
    mag = 10 ** math.floor(math.log10(raw))
    step = next((m * mag for m in (1, 2, 2.5, 5, 10) if raw <= m * mag), 10 * mag)
    return [i * step for i in range(int(hi / step) + 1)]


def med(v):
    v = sorted(v)
    n = len(v)
    if not n:
        return 0
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2


def median_curve(lists):
    """Per-rank median coverage across samples, ending at the median count."""
    n = len(lists)
    if not n:
        return []
    out = []
    for i in range(max(len(v) for v in lists)):
        vals = sorted(v[i] if i < len(v) else 0.0 for v in lists)
        m = (vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2)
        out.append(m)
    # int(x + 0.5), not round(): round() is banker's rounding
    return out[:int(med([len(v) for v in lists]) + 0.5)]


def cross_check_per_sample(per, counts_csv):
    """Compare proteoform totals sample by sample against the count table."""
    want = {}
    with open(counts_csv, newline="") as fh:
        for r in csv.DictReader(fh):
            if r["unit"].strip() == "protein_isoform_groups":
                want[(r["sample"].strip(), r["protease"].strip())] = int(r["count"])
    names = sorted({k for v in per.values() for k in v})
    print(f"\n  cross-check against {os.path.basename(counts_csv)}, per sample")
    worst = 0
    for name in names:
        rows = []
        for s_ in sorted(per):
            got = len(per[s_].get(name, []))
            exp = want.get((s_, name))
            if not got and exp is None:
                continue
            d = got - (exp or 0)
            worst = max(worst, abs(d))
            if abs(d) > 2:
                rows.append(f"{s_}:{got:,} vs {exp if exp is not None else 'absent'}")
        print(f"    {name:<8s} {len(rows)} sample(s) differ by more than 2"
              + ("   " + ", ".join(rows[:5]) if rows else ""))
    return worst


def cross_check(series, counts_csv):
    """Check pooled curves against the count table's per-sample medians."""
    want = {}
    with open(counts_csv, newline="") as fh:
        for r in csv.DictReader(fh):
            if r["unit"].strip() == "protein_isoform_groups":
                want.setdefault(r["protease"].strip(), []).append(int(r["count"]))
    print(f"\n  cross-check against {os.path.basename(counts_csv)} (pooled)")
    worst = 0
    for name, vals in series.items():
        w = want.get(name)
        if not w:
            print(f"    {name:<8s} not in the counts table")
            continue
        exp = med(w)
        d = len(vals) - exp
        worst = max(worst, abs(d))
        flag = "ok" if abs(d) <= 2 else "MISMATCH"
        print(f"    {name:<8s} curve {len(vals):>7,}   counts median "
              f"{exp:>9,.0f}   diff {d:>+6.0f}   {flag}")
    return worst


def draw_beeswarm(prot, out, letter):
    """Per-protein coverage as a beeswarm plus a box, one column per series."""
    W, H, ts = WIDTH, HEIGHT, TEXT_SCALE
    names = [n for n in ("GluC", "LysC", "Trypsin", "All")
             if any(n in r for r in prot.values())]
    full = prot
    prot = {g: r for g, r in prot.items() if all(n in r for n in names)}
    if not prot:
        sys.exit("no protein is measured by every series")
    colour = assign(names)
    sz = {"letter": 13.0 * ts, "tick": 9.5 * ts, "axis": 11.0 * ts,
          "cat": 10.5 * ts, "n": 9.5 * ts}
    ml, mr = 34.0 + 46.0 * ts, 34.0
    mt, mb = 18.0 + 21.0 * ts, 26.0 + 26.0 * ts
    pw, ph = W - ml - mr, H - mt - mb
    c = Canvas(W, H, FONT)
    cw = pw / len(names)
    bin_h = 1.6 / ph

    def Y(f):
        return mt + ph - f * ph

    if letter:
        c.text(22, 18 + 12 * ts, letter, sz["letter"], INK, "start", "600")
    for f in (0, .25, .5, .75, 1.0):
        c.line(ml, Y(f), ml + pw, Y(f), stroke=GRID, sw=1)
        c.text(ml - 9, Y(f) + 3.5 * ts, f"{f:.0%}", sz["tick"], INK_MUTED, "end")

    for k, nm in enumerate(names):
        vals = sorted(r[nm] for r in prot.values() if nm in r)
        if not vals:
            continue
        cx = ml + (k + 0.5) * cw
        band = cw * 0.56
        dots = [f'<circle cx="{cx + off:.1f}" cy="{yv:.2f}" r="1.05"/>'
                for _g, off, yv in swarm(((g, r[nm]) for g, r in prot.items() if nm in r),
                                         Y, bin_h, band)]
        c.add(f'<g fill="{colour[nm]}" fill-opacity="0.55">' + "".join(dots)
              + "</g>")

        def q(p_):
            i = (len(vals) - 1) * p_
            lo = int(i)
            return vals[lo] + (i - lo) * (vals[min(lo + 1, len(vals) - 1)] - vals[lo])
        q1, q2, q3 = q(.25), q(.5), q(.75)
        iqr = q3 - q1
        inside = [v for v in vals if q1 - 1.5 * iqr <= v <= q3 + 1.5 * iqr]
        bw = cw * 0.26
        c.line(cx, Y(min(inside)), cx, Y(max(inside)), stroke=INK, sw=1.2)
        c.rect(cx - bw / 2, Y(q3), bw, Y(q1) - Y(q3), fill="#ffffff", fo=0.0,
               stroke=INK, sw=1.2, rx=0)
        c.line(cx - bw / 2, Y(q2), cx + bw / 2, Y(q2), stroke=INK, sw=2.4)
        # white backing: the label sits over the swarm
        lab, fsz = f"{med(vals):.0%}", sz["n"] * c.fs * TEXT_BOOST
        lx, lw = cx + bw / 2 + 4, text_width(lab, fsz, True) + 6
        c.rect(lx - 3, Y(q2) - 0.62 * fsz, lw, 1.24 * fsz, fill="#ffffff",
               fo=0.85, rx=2)
        c.text(lx, Y(q2) + 0.36 * fsz, lab, sz["n"], INK, "start", "600")
        # a label wider than its column wraps at its spaces
        lab = display(nm)
        lines = lab.split() if text_width(lab, sz["cat"] * TEXT_BOOST) > 0.9 * cw else [lab]
        for j, ln in enumerate(lines):
            c.text(cx, mt + ph + 16 * ts + j * 1.15 * sz["cat"] * TEXT_BOOST, ln, sz["cat"],
                   INK, "middle")
    c.line(ml, mt + ph, ml + pw, mt + ph, stroke=AXIS, sw=1)
    c.text(ml, mt - 14 * ts,
           f"n = {len(prot):,} protein groups",
           sz["n"], INK, "start")
    c.add(f'<g transform="translate({18 + 10 * ts:.1f} {mt + ph / 2:.1f}) '
          f'rotate(-90)"><text x="0" y="0" font-size="{sz["axis"]:.1f}" '
          f'text-anchor="middle" fill="{INK}">Sequence coverage</text></g>')
    write(c, out)
    for nm in names:
        vals = sorted(r[nm] for r in prot.values() if nm in r)
        alv = sorted(r[nm] for r in full.values() if nm in r)
        print(f"  {nm:<8s} common n={len(vals):>7,} median {med(vals):6.1%}"
              f"   |   all n={len(alv):>7,} median {med(alv):6.1%}")
    return W, H


def draw_rank(series, out, letter, curves=None):
    """series: {name: [coverage fractions, descending]}"""
    W, H, ts = WIDTH, HEIGHT, TEXT_SCALE
    c = Canvas(W, H, FONT)
    sz = {"letter": 13.0 * ts, "tick": 9.5 * ts, "name": 10.5 * ts,
          "sub": 9.5 * ts, "axis": 11.0 * ts}
    tb = ts * TEXT_BOOST
    row = 1.08 * sz["name"] * TEXT_BOOST
    ylab = wrap2("Fraction of sequence covered")
    pitch = 1.15 * sz["axis"] * TEXT_BOOST
    ml = 34.0 + 46.0 * ts + 0.5 * pitch * (len(ylab) - 1)
    mr = 18.0
    mt, mb = 32.0 + 4 * row, 20.0 + 22.0 * tb
    pw, ph = W - ml - mr, H - mt - mb
    order = [n for n in ("All", "Trypsin", "LysC", "GluC") if n in series] + \
            [n for n in series if n not in ("All", "Trypsin", "LysC", "GluC")]
    colour = assign(list(series))
    xmax = max(len(v) for v in series.values()) or 1

    if letter:
        c.text(22, 18 + 12 * ts, letter, sz["letter"], INK, "start", "600")

    def X(i):
        return ml + i / xmax * pw

    def Y(f):
        return mt + ph - f * ph

    for f in (0, .25, .5, .75, 1.0):
        c.line(ml, Y(f), ml + pw, Y(f), stroke=GRID, sw=1)
        c.text(ml - 9, Y(f) + 3.5 * ts, f"{f:.2f}", sz["tick"], INK_MUTED, "end")
    for t in nice_ticks(xmax):
        if t <= xmax:
            c.line(X(t), mt + ph, X(t), mt + ph + 4, stroke=AXIS, sw=1)
            c.text(X(t), mt + ph + 13 * ts, compact(t, exact=True), sz["tick"],
                   INK_MUTED, "middle")
    c.line(ml, mt + ph, ml + pw, mt + ph, stroke=AXIS, sw=1)

    for name in order:
        vals = series[name]
        if not vals:
            continue
        # one point per pixel column
        step = max(1, len(vals) // int(pw))
        pts = [(X(i), Y(vals[i])) for i in range(0, len(vals), step)]
        pts.append((X(len(vals) - 1), Y(vals[-1])))
        d = "M " + " L ".join(f"{x:.1f} {y:.2f}" for x, y in pts)
        c.add(f'<path d="{d}" fill="none" stroke="{colour[name]}" stroke-width="3"'
              f' stroke-linejoin="round" stroke-linecap="round"/>')

    # legend above the plot: one row per series, right-aligned to the plot
    subs = {}
    for name in order:
        v = series[name]
        if curves is not None:
            half = [sum(1 for x in u if x >= 0.5) for u in curves.get(name, [])]
            subs[name] = f"{compact(med(half))} at ≥ 50% · n = {len(half)}"
        else:
            subs[name] = f"{compact(sum(1 for x in v if x >= 0.5))} at ≥ 50%"
    name_w = max(text_width(display(n), sz["name"] * TEXT_BOOST, True) for n in order)
    sub_w = max(text_width(t, sz["sub"] * TEXT_BOOST) for t in subs.values())
    lx = ml + pw - (18 * ts + name_w + 30 + sub_w)
    ly = 32.0
    for name in order:
        sy = ly - 0.35 * sz["name"] * TEXT_BOOST
        c.line(lx, sy, lx + 14 * ts, sy, stroke=colour[name], sw=3)
        tx = lx + 18 * ts
        c.text(tx, ly, display(name), sz["name"], INK, "start", "600")
        c.text(tx + name_w + 30, ly, subs[name], sz["sub"], INK_MUTED, "start")
        ly += row

    c.text(ml + pw / 2, H - 9 * ts, "Protein group rank", sz["axis"], INK, "middle")
    yc = mt + ph / 2
    c.add(f'<g transform="translate({18 + 10 * ts - 5 * (len(ylab) - 1):.1f} {yc:.1f}) '
          f'rotate(-90)">' + "".join(
              f'<text x="0" y="{j * pitch:.1f}" font-size="{sz["axis"]:.1f}" '
              f'text-anchor="middle" fill="{INK}">{t}</text>' for j, t in enumerate(ylab))
          + "</g>")
    write(c, out)
    return W, H


def wrap2(label):
    """Two lines, split at the space nearest the middle."""
    cut = min((i for i, ch in enumerate(label) if ch == " "),
              key=lambda i: abs(i - len(label) / 2), default=None)
    return [label] if cut is None else [label[:cut], label[cut + 1:]]


def write(c, out):
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reports", nargs="+")
    ap.add_argument("--fasta", required=True, help="the FASTA the search used")
    ap.add_argument("--rank", action="store_true", help="rank-coverage curves")
    ap.add_argument("--beeswarm", action="store_true",
                    help="per-protein coverage as a box plus point cloud")
    ap.add_argument("--by-sample", action="store_true",
                    help="median over per-sample coverage")
    ap.add_argument("--counts", default=None, metavar="COUNTS_CSV",
                    help="cross-check proteoform totals against the count-panel table")
    ap.add_argument("--sample-regex", default=None,
                    help="regex on the run name; group 1 is the sample id")
    ap.add_argument("--letter", default="")
    ap.add_argument("--out", default="coverage.svg")
    args = ap.parse_args(argv)

    if not (args.rank or args.beeswarm):
        sys.exit("choose --rank or --beeswarm")

    paths = []
    for pat in args.reports:
        paths.extend(sorted(glob.glob(pat)) or [pat])
    from lib_report import one_search_per_run
    one_search_per_run(paths)
    seqs = read_fasta(args.fasta)
    print(f"{len(seqs):,} sequences in {os.path.basename(args.fasta)}")

    print("reading reports:")
    groups = collect(paths, args.sample_regex, split_samples=args.by_sample)
    if args.by_sample:
        print(f"{len(groups):,} samples")
    else:
        print(f"{len(groups):,} protein groups with peptides")

    def cov_one(gs):
        """One sample's per-protein coverage: {group: {series: fraction}}."""
        out = {}
        for grp, byp in gs.items():
            rep = representative(grp, seqs)
            if rep is None or not seqs[rep]:
                continue
            sq = seqs[rep]
            row, allp = {}, set()
            for pn, st in byp.items():
                m, _f = mark(sq, st)
                row[pn] = sum(m) / len(sq)
                allp |= st
            m, _ = mark(sq, allp)
            row["All"] = sum(m) / len(sq)
            out[grp] = row
        return out

    def rank_one(gs):
        """One sample's (or the pool's) sorted coverage list per series."""
        out = defaultdict(list)
        miss = unmap = 0
        for grp, byp in gs.items():
            rep = representative(grp, seqs)
            if rep is None:
                miss += 1
                continue
            sq = seqs[rep]
            if not sq:
                continue
            allp = set()
            for p, st in byp.items():
                m, f = mark(sq, st)
                unmap += len(st) - f
                out[p].append(sum(m) / len(sq))
                allp |= st
            m, _ = mark(sq, allp)
            out["All"].append(sum(m) / len(sq))
        for k in out:
            out[k].sort(reverse=True)
        return out, miss, unmap

    if args.beeswarm:
        acc = defaultdict(lambda: defaultdict(list))
        src = groups if args.by_sample else {"pooled": groups}
        for s_, gs in sorted(src.items()):
            for grp, row in cov_one(gs).items():
                for k, v in row.items():
                    acc[grp][k].append(v)
        prot = {g: {k: med(v) for k, v in r.items()} for g, r in acc.items()}
        print(f"{len(prot):,} proteins with coverage in >=1 sample")
        W, H = draw_beeswarm(prot, args.out, args.letter)
        print(f"\nwrote {args.out}  ({W}x{H}px)")
        return 0

    if args.by_sample:
        per = {s_: rank_one(gs)[0] for s_, gs in sorted(groups.items())}
        names = sorted({k for v in per.values() for k in v})
        keep = None
        if args.counts:
            keep = defaultdict(set)
            with open(args.counts, newline="") as fh:
                for r in csv.DictReader(fh):
                    if r["unit"].strip() == "protein_isoform_groups":
                        keep[r["protease"].strip()].add(r["sample"].strip())
            # only the samples the count panel counted
            dropped = {k: sorted(s_ for s_ in per if per[s_].get(k)
                                 and s_ not in keep.get(k, set())) for k in names}
            for k, ds in dropped.items():
                if ds:
                    print(f"  {k}: dropping {', '.join(ds)} -- not in the count "
                          f"table for this series")
        curves = {k: [per[s_][k] for s_ in per if per[s_].get(k)
                      and (keep is None or s_ in keep.get(k, set()))]
                  for k in names}
        series = {k: median_curve(v) for k, v in curves.items() if v}
        print("\nper-sample rank curves, median across samples:")
        for k in sorted(series, key=lambda k: -len(series[k])):
            ns = [len(v) for v in curves[k]]
            print(f"  {k:<8s} {len(curves[k]):>2} samples   "
                  f"proteoforms {min(ns):,}-{max(ns):,}   "
                  f"median curve {len(series[k]):,}")
        if args.counts:
            drawn = {s_: {k: v for k, v in per[s_].items()
                          if keep is None or s_ in keep.get(k, set())}
                     for s_ in per}
            cross_check_per_sample(drawn, args.counts)
        W, H = draw_rank(series, args.out, args.letter, curves=curves)
        print(f"\nwrote {args.out}  ({W}x{H}px)")
        for k in sorted(series, key=lambda k: -len(series[k])):
            m = med([sum(1 for x in v if x >= 0.5) for v in curves[k]])
            print(f"  {k:<8s} median {len(series[k]):>7,} proteoforms   "
                  f"median >=50% {m:>7,.0f}")
        return 0

    series, missing, unmapped = rank_one(groups)
    if missing:
        print(f"{missing:,} group(s) had no member in the FASTA -- skipped")
    if unmapped:
        print(f"{unmapped:,} peptide-group pairs did not occur in the "
              f"representative sequence")
    if args.counts:
        cross_check(series, args.counts)
    W, H = draw_rank(series, args.out, args.letter)
    print(f"\nwrote {args.out}  ({W}x{H}px)")
    for k in sorted(series, key=lambda k: -len(series[k])):
        v = series[k]
        mid = v[len(v) // 2] if v else 0
        print(f"  {k:<8s} n = {len(v):>7,}   >=50% {sum(1 for x in v if x >= 0.5):>7,}"
              f"   >=80% {sum(1 for x in v if x >= 0.8):>6,}"
              f"   (median of its own list {mid:.1%})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
