#!/usr/bin/env python3
"""Sequence coverage from Spectronaut reports: a rank curve, or one protein.

    # every proteoform, ranked by coverage, one curve per protease
    python3 fig1cd_coverage.py 'data/search/*-60min-Phospho.parquet' --fasta search_db.fasta --rank \
        --out figures/coverage_rank.svg

    # one protein, residue by residue, with an optional isoform track
    python3 fig1cd_coverage.py 'data/search/*-60min-Phospho.parquet' --fasta search_db.fasta \
        --protein P37840 --isoform P37840-3 --out figures/snca.svg
"""

import argparse
import glob
import math
import os
import re
import sys
import collections
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_svg import Canvas, esc                          # noqa: E402
from lib_fasta import read_fasta                                   # noqa: E402
from lib_palette import (AXIS, DEPTH, EMPTY, EMPTY_DARK, FONT, GRID, INK, INK_MUTED,
                     INK_SECONDARY, STRIP_FILL, SURFACE, UNION, assign,
                     display, ink_on, ramp_for, TEXT_BOOST)      # noqa: E402
from prep_counts import (detect_protease, parse_float, read_rows,  # noqa: E402
                       resolve_columns)

ISO = re.compile(r"-\d+$")


# ------------------------------------------------------------------ inputs ---


def collect(paths, precursor_q, accession=None, sample_regex=None, sample=None,
            split_samples=False, protein_q=None):
    """-> ({protein_group: {protease: set(peptides)}}, {(protease, sample): n})"""
    out = (defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
           if split_samples else defaultdict(lambda: defaultdict(set)))
    seen_samples = defaultdict(int)
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
                    if q is None or q > precursor_q:
                        continue
                if "decoy" in idx and row[idx["decoy"]].strip().lower() == "true":
                    continue
                # same protein-group cutoff as prep_counts.py
                if protein_q is not None and "pg_qvalue" in idx:
                    q = parse_float(row[idx["pg_qvalue"]])
                    if q is None or q > protein_q:
                        continue
                pep = row[idx["strip_seq"]].strip()
                grp = row[idx["protein_group"]].strip()
                if not (pep and grp):
                    continue
                if accession and accession not in grp:
                    continue
                run = row[idx["run"]]
                if srx:
                    m = srx.search(run)
                    s = (m.group(1) if m and m.groups() else
                         m.group(0)) if m else None
                else:
                    s = None
                if sample and s != sample:
                    continue
                p = protease or detect_protease(run)
                if not p or p.startswith("AMBIGUOUS"):
                    p = "?"
                if split_samples:
                    out[s][grp][p].add(pep)
                else:
                    out[grp][p].add(pep)
                seen_samples[(p, s)] += 1
                n += 1
            print(f"  {os.path.basename(path)[:52]:<52s} {protease or '?':<8s} "
                  f"{n:>9,} rows kept")
        finally:
            if fh is not None:
                fh.close()
    return out, seen_samples


# ---------------------------------------------------------------- coverage ---
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


def representative(group, seqs, mode="first"):
    """Which accession in a group the coverage is measured against."""
    accs = [a.strip() for a in re.split(r"[;,]", group) if a.strip()]
    accs = [a.split("|")[1] if a.count("|") >= 2 else a for a in accs]
    present = [a for a in accs if a in seqs]
    if not present:
        return None
    if mode == "longest":
        return max(present, key=lambda a: len(seqs[a]))
    return present[0]


def deletion_map(canon, iso):
    """Map an isoform built from `canon` by one deletion back to canonical."""
    d = len(canon) - len(iso)
    if d <= 0:
        return None
    p = 0
    while p < len(iso) and canon[p] == iso[p]:
        p += 1
    s = 0
    while s < len(iso) - p and canon[len(canon) - 1 - s] == iso[len(iso) - 1 - s]:
        s += 1
    if p + s != len(iso):
        return None
    return p, d


def iso_segments(a, b, p, d, iso_len):
    """Isoform interval [a,b) -> list of canonical intervals."""
    segs = []
    if a < p:
        segs.append((a, min(b, p)))
    if b > p:
        segs.append((max(a, p) + d, b + d))
    return [(x, y) for x, y in segs if y > x]


# ------------------------------------------------------------------- plots ---
def fmt_thousands(v):
    return f"{v:,.0f}"


def compact(v, exact=False):
    """Format a count as 377K / 14.7K / 6.4K."""
    from fig1b_depth import fmt_compact
    return fmt_compact(v, exact)


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
    # no early stop at zero coverage: a zero is a real zero
    out = []
    for i in range(max(len(v) for v in lists)):
        vals = sorted(v[i] if i < len(v) else 0.0 for v in lists)
        m = (vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2)
        out.append(m)
    # int(x + 0.5), not round(): round() is banker's rounding
    return out[:int(med([len(v) for v in lists]) + 0.5)]


def cross_check_per_sample(per, counts_csv):
    """Compare proteoform totals sample by sample against the count table."""
    import csv as _csv
    want = {}
    with open(counts_csv, newline="") as fh:
        for r in _csv.DictReader(fh):
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


def cross_check(series, counts_csv, sample):
    """Check these curves against the counts panel on the same sample."""
    import csv as _csv
    want = {}
    with open(counts_csv, newline="") as fh:
        for r in _csv.DictReader(fh):
            if r["unit"].strip() != "protein_isoform_groups":
                continue
            if sample and r["sample"].strip() != str(sample):
                continue
            want.setdefault(r["protease"].strip(), []).append(int(r["count"]))
    stat = "median" if not sample else "value"
    print("\n  cross-check against " + os.path.basename(counts_csv)
          + (f", sample {sample}" if sample else " (pooled)"))
    worst = 0
    for name, vals in series.items():
        w = want.get(name)
        if not w:
            print(f"    {name:<8s} not in the counts table")
            continue
        w = sorted(w)
        exp = (w[len(w) // 2] if len(w) % 2 else
               (w[len(w) // 2 - 1] + w[len(w) // 2]) / 2)
        d = len(vals) - exp
        worst = max(worst, abs(d))
        flag = "ok" if abs(d) <= 2 else "MISMATCH"
        print(f"    {name:<8s} curve {len(vals):>7,}   counts {stat} "
              f"{exp:>9,.0f}   diff {d:>+6.0f}   {flag}")
    if worst > 2:
        print("    ^ more than the one group whose accession is absent from the "
              "FASTA; the two panels would disagree")
    return worst


def draw_paired(prot, out, font, letter, width=1215.0, height=430.0, ts=1.7,
                against="Trypsin"):
    """Coverage against protein rank, every series read off the same protein."""
    names = [against] + [n for n in ("All",) if n != against]
    colour = assign(list({k for r in prot.values() for k in r}))
    rows = sorted((r for r in prot.values() if against in r),
                  key=lambda r: r[against])
    n = len(rows) or 1
    sz = {"letter": 13.0 * ts, "tick": 9.5 * ts, "axis": 11.0 * ts,
          "name": 10.0 * ts}
    ml, mr = 34.0 + 46.0 * ts, 34.0
    mt, mb = 30.0 + 14.0 * ts, 20.0 + 22.0 * ts
    W, H = width, height
    pw, ph = W - ml - mr, H - mt - mb
    c = Canvas(W, H, font)

    def X(i):
        return ml + i / max(n - 1, 1) * pw

    def Y(f):
        return mt + ph - f * ph

    if letter:
        c.text(22, 18 + 12 * ts, letter, sz["letter"], INK, "start", "600")
    for f in (0, .25, .5, .75, 1.0):
        c.line(ml, Y(f), ml + pw, Y(f), stroke=GRID, sw=1)
        c.text(ml - 9, Y(f) + 3.5 * ts, f"{f:.0%}", sz["tick"], INK_MUTED, "end")
    for t in nice_ticks(n):
        if t <= n:
            c.line(X(t), mt + ph, X(t), mt + ph + 4, stroke=AXIS, sw=1)
            c.text(X(t), mt + ph + 12 * ts, compact(t, exact=True), sz["tick"],
                   INK_MUTED, "middle")
    c.line(ml, mt + ph, ml + pw, mt + ph, stroke=AXIS, sw=1)

    # union first, so the sorted line draws on top
    step = max(1, n // 5200)
    dots = []
    for i in range(0, n, step):
        v = rows[i].get("All")
        if v is not None:
            dots.append(f'<circle cx="{X(i):.1f}" cy="{Y(v):.2f}" r="1.15"/>')
    c.add(f'<g fill="{colour.get("All", UNION)}" fill-opacity="0.30">'
          + "".join(dots) + "</g>")

    pts = [(X(i), Y(rows[i][against])) for i in range(0, n, max(1, n // int(pw)))]
    pts.append((X(n - 1), Y(rows[-1][against])))
    c.add('<path d="M ' + " L ".join(f"{x:.1f} {y:.2f}" for x, y in pts)
          + f'" fill="none" stroke="{colour[against]}" stroke-width="3" '
          f'stroke-linejoin="round"/>')

    c.text(ml + 14, mt + 10 + 8 * ts, "All proteases", sz["name"],
           colour.get("All", UNION), "start", "600")
    c.text(ml + pw - 14, Y(rows[max(0, int(n * 0.86))][against]) + 24 * ts,
           display(against), sz["name"], colour[against], "end", "600")
    c.text(ml + pw / 2, H - 12 * ts, f"Protein group rank, sorted by {display(against)}",
           sz["axis"], INK, "middle")
    c.add(f'<g transform="translate({18 + 10 * ts:.1f} {mt + ph / 2:.1f}) '
          f'rotate(-90)"><text x="0" y="0" font-size="{sz["axis"]:.1f}" '
          f'text-anchor="middle" fill="{INK}">Sequence coverage</text></g>')
    missing = len(prot) - n
    if missing > 0:
        c.text(ml + pw, mt - 8,
               f"{missing:,} proteins {display(against)} never saw are absent — "
               f"they have no rank on this axis",
               9.0 * ts, INK_MUTED, "end")
    write(c, out)
    gain = [r["All"] - r[against] for r in rows if "All" in r]
    print(f"  sorted on {against}; median per-protein gain "
          f"{med(gain):.1%}, mean {sum(gain)/len(gain):.1%}")
    print(f"  {missing:,} of {len(prot):,} proteins are absent from this plot "
          f"because {against} did not see them; the gain above is therefore a "
          f"lower bound on what the union adds")
    return W, H


def draw_beeswarm(prot, out, font, letter, width=1215.0, height=430.0, ts=1.7,
                  common=True):
    """Per-protein coverage as a beeswarm plus a box, one column per series."""
    import hashlib
    names = [n for n in ("GluC", "LysC", "Trypsin", "All")
             if any(n in r for r in prot.values())]
    full = prot
    if common:
        prot = {g: r for g, r in prot.items() if all(n in r for n in names)}
        if not prot:
            sys.exit("no protein is measured by every series")
    colour = assign(names)
    sz = {"letter": 13.0 * ts, "tick": 9.5 * ts, "axis": 11.0 * ts,
          "cat": 10.5 * ts, "n": 8.8 * ts}
    ml, mr = 34.0 + 46.0 * ts, 34.0
    mt, mb = 26.0 + 30.0 * ts, 26.0 + 26.0 * ts
    W, H = width, height
    pw, ph = W - ml - mr, H - mt - mb
    c = Canvas(W, H, font)
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
        # beeswarm: bin by coverage, lay points out symmetrically within each bin
        band = cw * 0.62
        rows = collections.defaultdict(list)
        for g, r in prot.items():
            if nm in r:
                rows[round(r[nm] / bin_h)].append(g)
        widest = max((len(v) for v in rows.values()), default=1)
        dots = []
        for b_, members in rows.items():
            yv = Y(b_ * bin_h)
            # width scaled to the fullest bin, capped to stay inside the column
            span = band * min(1.0, (len(members) / widest) ** 0.55)
            n_ = len(members)
            for i, g in enumerate(sorted(members)):
                off = 0.0 if n_ == 1 else (i / (n_ - 1) - 0.5) * span
                dots.append(f'<circle cx="{cx + off:.1f}" cy="{yv:.2f}" r="1.05"/>')
        c.add(f'<g fill="{colour[nm]}" fill-opacity="0.55">' + "".join(dots)
              + "</g>")

        def q(p_):
            i = (len(vals) - 1) * p_
            lo = int(i)
            return vals[lo] + (i - lo) * (vals[min(lo + 1, len(vals) - 1)] - vals[lo])
        q1, q2, q3 = q(.25), q(.5), q(.75)
        iqr = q3 - q1
        inside = [v for v in vals if q1 - 1.5 * iqr <= v <= q3 + 1.5 * iqr]
        bw = cw * 0.30
        c.line(cx, Y(min(inside)), cx, Y(max(inside)), stroke=INK, sw=1.2)
        c.rect(cx - bw / 2, Y(q3), bw, Y(q1) - Y(q3), fill="#ffffff", fo=0.0,
               stroke=INK, sw=1.2, rx=0)
        c.line(cx - bw / 2, Y(q2), cx + bw / 2, Y(q2), stroke=INK, sw=2.4)
        # white backing: the label sits over the swarm
        lab, fsz = f"{med(vals):.0%}", sz["n"] * c.fs
        lx, lw = cx + bw / 2 + 5, 0.56 * fsz * len(lab) + 6
        c.rect(lx - 3, Y(q2) - 0.62 * fsz, lw, 1.24 * fsz, fill="#ffffff",
               fo=0.85, rx=2)
        c.text(lx, Y(q2) + 0.36 * fsz, lab, sz["n"], INK, "start", "600")
        c.text(cx, mt + ph + 16 * ts, display(nm), sz["cat"], INK, "middle")
    c.line(ml, mt + ph, ml + pw, mt + ph, stroke=AXIS, sw=1)
    if common:
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
    print("  the 'all' column is not a fair comparison: each series summarises a "
          "different protein set")
    return W, H


def draw_rank(series, out, font, letter, width=1215.0, height=430.0, ts=1.7,
              sample=None, curves=None, show_curves=False):
    """series: {name: [coverage fractions, descending]}"""
    W, H = width, height
    c = Canvas(W, H, font)
    sz = {"letter": 13.0 * ts, "tick": 9.5 * ts, "name": 10.0 * ts,
          "sub": 8.5 * ts, "axis": 11.0 * ts}
    tb = ts * TEXT_BOOST
    ml = 34.0 + 46.0 * tb
    mr = 34.0
    mt, mb = 30.0 + 14.0 * tb, 20.0 + 22.0 * tb
    pw, ph = W - ml - mr, H - mt - mb
    order = [n for n in ("All", "Trypsin", "LysC", "GluC") if n in series] + \
            [n for n in series if n not in ("All", "Trypsin", "LysC", "GluC")]
    colour = assign(list(series))
    # only what is drawn sets the x range
    xmax = max([len(v) for v in series.values()]
               + ([len(v) for vs in (curves or {}).values() for v in vs]
                  if show_curves else [])) or 1

    if letter:
        c.text(22, 18 + 12 * ts, letter, sz["letter"], INK, "start", "600")
    if sample:
        c.text(ml, mt - 12, f"sample {sample}", sz["sub"], INK_MUTED, "start")

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
            c.text(X(t), mt + ph + 12 * ts, compact(t, exact=True), sz["tick"],
                   INK_MUTED, "middle")
    c.line(ml, mt + ph, ml + pw, mt + ph, stroke=AXIS, sw=1)

    def polyline(vals, col, sw, opacity=None):
        # one point per pixel column
        step = max(1, len(vals) // int(pw))
        pts = [(X(i), Y(vals[i])) for i in range(0, len(vals), step)]
        pts.append((X(len(vals) - 1), Y(vals[-1])))
        d = "M " + " L ".join(f"{x:.1f} {y:.2f}" for x, y in pts)
        op = f' stroke-opacity="{opacity}"' if opacity is not None else ""
        c.add(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="{sw}"'
              f'{op} stroke-linejoin="round" stroke-linecap="round"/>')

    if show_curves:
        for name in order:
            for vals in (curves or {}).get(name, []):
                if vals:
                    polyline(vals, colour[name], 0.9, opacity=0.28)
    for name in order:
        vals = series[name]
        if vals:
            polyline(vals, colour[name], 2.6)

    # two counts, not a median: medians of unequal-length lists don't compare
    ly = mt + 10 * ts
    lx = ml + pw - (150.0 + 26.0 * ts) * TEXT_BOOST
    for name in order:
        v = series[name]
        c.line(lx, ly - 4, lx + 20 * ts, ly - 4, stroke=colour[name], sw=2.4)
        tx = lx + 26 * ts
        c.text(tx, ly + 2, display(name), sz["name"], INK, "start", "600")
        if curves is not None:
            half = sorted(sum(1 for x in u if x >= 0.5) for u in curves.get(name, []))
            n = len(half) or 1
            hm = (half[n // 2] if n % 2 else (half[n // 2 - 1] + half[n // 2]) / 2) \
                if half else 0
            c.text(tx, ly + 10 + 8 * ts,
                   f"{compact(hm)} at \u2265 50% · n = {len(half)}",
                   sz["sub"], INK_MUTED, "start")
        else:
            c.text(tx, ly + 10 + 8 * ts,
                   f"{compact(sum(1 for x in v if x >= 0.5))} at \u2265 50%",
                   sz["sub"], INK_MUTED, "start")
        ly += 20 * ts + 16

    c.text(ml + pw / 2, H - 12 * ts, "Protein group rank", sz["axis"], INK, "middle")
    yc = mt + ph / 2
    c.add(f'<g transform="translate({18 + 10 * ts:.1f} {yc:.1f}) rotate(-90)">'
          f'<text x="0" y="0" font-size="{sz["axis"]:.1f}" text-anchor="middle" '
          f'fill="{INK}">Fraction of sequence covered</text></g>')
    write(c, out)
    return W, H


def merged(mask):
    """Contiguous covered spans of a residue mask, as (start, end)."""
    out, i, n = [], 0, len(mask)
    while i < n:
        if mask[i]:
            j = i
            while j < n and mask[j]:
                j += 1
            out.append((i, j))
            i = j
        else:
            i += 1
    return out


def locate(seq, peptides):
    """Every occurrence of every peptide, as (start, end) intervals."""
    out = []
    for p in peptides:
        i = seq.find(p)
        while i != -1:
            out.append((i, i + len(p)))
            i = seq.find(p, i + 1)
    return out


def pack(intervals, gap):
    """Greedy interval packing into tiers, like a read pileup."""
    tiers = []
    for s, e in sorted(intervals, key=lambda t: (t[0], -(t[1] - t[0]))):
        for row in tiers:
            if row[-1][1] + gap <= s:
                row.append((s, e))
                break
        else:
            tiers.append([(s, e)])
    return tiers


def depth_counts(seq, peptides):
    """How many distinct peptides cover each residue."""
    d = [0] * len(seq)
    for s, e in locate(seq, peptides):
        for i in range(s, e):
            d[i] += 1
    return d


def draw_grid(acc, seq, series, out, font, letter, label=None, sample=None,
              per_row=50, cols=1, group=10, tile=17.0, pad=0.0, gap=0.0,
              empty=EMPTY_DARK, radius=0.0):
    """Sequence laid out as tiles, each shaded by peptides-per-residue."""
    n = len(seq)
    names = list(series)
    colour = assign(names)
    rows = (n + per_row - 1) // per_row

    # row-major across `cols` columns
    gutter = 34.0
    row_w = per_row * (tile + pad) + ((per_row - 1) // group) * gap
    block_w = gutter + row_w
    col_gap = 46.0
    ml, mr = 30.0, 24.0
    W = ml + cols * block_w + (cols - 1) * col_gap + mr

    head_h, row_h, block_gap = 26.0, tile + 4.0, 26.0
    block_h = head_h + rows * row_h + block_gap
    nrows_of_blocks = (len(names) + cols - 1) // cols
    top = 74.0 if letter else 54.0
    H = top + nrows_of_blocks * block_h - block_gap + 30

    c = Canvas(W, H, font)
    if letter:
        c.text(22, 32, letter, 13, INK, "start", "600")
    title = label or acc
    tx = ml + (26 if letter else 0)
    c.text(tx, 34, title, 12.5, INK, "start", "600")
    c.text(tx + len(title) * 7.8 + 16, 34,
           f"{acc} · {n} aa" + (f" · sample {sample}" if sample else " · all samples"),
           9.5, INK_MUTED, "start")

    for bi, name in enumerate(names):
        bx = ml + (bi % cols) * (block_w + col_gap)
        by = top + (bi // cols) * block_h
        peps = series[name]
        depth = depth_counts(seq, peps)
        cov = sum(1 for d in depth if d) / n
        ramp = ramp_for(colour[name], on_dark=(empty == EMPTY_DARK))

        c.rect(bx + gutter - 13, by + 4, 9, 9, fill=colour[name], rx=2)
        c.text(bx + gutter, by + 12, name, 11, INK, "start", "600")
        c.text(bx + gutter + len(name) * 7 + 10, by + 12,
               f"{cov:.0%} · {len(peps)} pep", 9.2, INK_MUTED, "start")

        kx = bx + block_w - (len(ramp) + 1) * 21
        for k, fill in enumerate([empty] + ramp):
            c.rect(kx, by + 4, 9, 9, fill=fill,
                   stroke=AXIS if fill in (EMPTY, "#ffffff") else "none",
                   sw=0.7, rx=radius)
            c.text(kx + 11, by + 12,
                   f"{k}+" if k == len(ramp) else str(k), 8.0, INK_MUTED, "start")
            kx += 21

        for r in range(rows):
            ry = by + head_h + r * row_h
            c.text(bx + gutter - 8, ry + tile * 0.72, str(r * per_row + 1), 8.2,
                   INK_MUTED, "end")
            for col in range(per_row):
                i = r * per_row + col
                if i >= n:
                    break
                d = depth[i]
                fill = empty if d == 0 else ramp[min(d, len(ramp)) - 1]
                x = bx + gutter + col * (tile + pad) + (col // group) * gap
                # touching tiles: no stroke, or hairlines double up on shared edges
                c.rect(x, ry, tile, tile, fill=fill,
                       stroke=AXIS if (d == 0 and fill in (EMPTY, "#ffffff"))
                       else "none", sw=0.6, rx=radius)
                c.text(x + tile / 2, ry + tile * 0.72, seq[i],
                       tile * 0.53, ink_on(fill), "middle")

    write(c, out)
    return W, H


def runs(depth, cap):
    """Merge consecutive residues of equal (capped) depth -> [(start, end, d)]."""
    out, i, n = [], 0, len(depth)
    while i < n:
        d = min(depth[i], cap)
        j = i + 1
        while j < n and min(depth[j], cap) == d:
            j += 1
        out.append((i, j, d))
        i = j
    return out


def draw_line(acc, seq, series, out, font, letter, label=None, sample=None,
              strip=15.0, empty=EMPTY_DARK, width=1010.0):
    """One unbroken strip per series, the whole sequence at a common scale."""
    n = len(seq)
    names = list(series)
    colour = assign(names)
    on_dark = empty == EMPTY_DARK

    ml, mr = 104.0, 30.0
    top = 88.0 if letter else 68.0
    head, gap = 15.0, 13.0
    row_h = head + strip + gap
    W = width
    pw = W - ml - mr
    H = top + len(names) * row_h + 34.0
    c = Canvas(W, H, font)

    def x_of(i):                       # 0-based residue index -> x
        return ml + i / n * pw

    if letter:
        c.text(22, 32, letter, 13, INK, "start", "600")
    title = label or acc
    tx = 48 if letter else 30
    c.text(tx, 34, title, 12.5, INK, "start", "600")
    c.text(tx + len(title) * 7.8 + 16, 34,
           f"{acc} · {n} aa" + (f" · sample {sample}" if sample else " · all samples"),
           9.5, INK_MUTED, "start")

    ay = top - 24
    step = 25 if n <= 200 else (50 if n <= 600 else 100)
    c.line(ml, ay + 5, ml + pw, ay + 5, stroke=AXIS, sw=1)
    for r in list(range(1, n + 1, step)) + [n]:
        if r != n and n - r < step * 0.45:
            continue                   # keep the last tick from colliding
        c.line(x_of(r - 1), ay + 5, x_of(r - 1), ay + 9, stroke=AXIS, sw=1)
        c.text(x_of(r - 1), ay, str(r), 7.8, INK_MUTED, "middle")

    for bi, name in enumerate(names):
        by = top + bi * row_h
        peps = series[name]
        depth = depth_counts(seq, peps)
        cov = sum(1 for d in depth if d) / n
        ramp = ramp_for(colour[name], on_dark=on_dark)

        c.rect(ml - 13, by + 2, 9, 9, fill=colour[name], rx=2)
        c.text(ml, by + 10, name, 11, INK, "start", "600")
        c.text(ml + len(name) * 7 + 10, by + 10,
               f"{cov:.0%} · {len(peps)} pep", 9.2, INK_MUTED, "start")

        kx = ml + pw - (len(ramp) + 1) * 21
        for k, fill in enumerate([empty] + ramp):
            c.rect(kx, by + 2, 9, 9, fill=fill,
                   stroke=AXIS if fill in (EMPTY, "#ffffff") else "none",
                   sw=0.7, rx=0)
            c.text(kx + 11, by + 10,
                   f"{k}+" if k == len(ramp) else str(k), 8.0, INK_MUTED, "start")
            kx += 21

        sy = by + head
        c.rect(ml, sy, pw, strip, fill=empty, rx=0)
        for i, j, d in runs(depth, len(ramp)):
            if d == 0:
                continue
            c.rect(x_of(i), sy, max(x_of(j) - x_of(i), 0.35), strip,
                   fill=ramp[d - 1], rx=0)
        if empty in (EMPTY, "#ffffff"):
            c.rect(ml, sy, pw, strip, fill="none", stroke=AXIS, sw=0.6, rx=0)

    c.text(ml + pw / 2, H - 12, "Residue", 9, INK_SECONDARY, "middle")
    write(c, out)
    return W, H


def draw_stack(acc, seq, series, out, font, letter, label=None, sample=None,
               per_row=60, tile=14.0, track=10.0, empty=None, radius=0.0,
               union_name="All"):
    """The sequence once, with one track per digest underneath it."""
    empty = EMPTY_DARK if empty is None else empty
    on_dark = empty == EMPTY_DARK
    n = len(seq)
    names = list(series)
    colour = assign(names)
    masks = {nm: mark(seq, series[nm])[0] for nm in names}
    digests = [nm for nm in names if nm != union_name]
    ndig = [sum(masks[nm][i] for nm in digests) for i in range(n)]
    uramp = ramp_for(colour.get(union_name, UNION), on_dark=on_dark)[:max(len(digests), 1)]

    rows = (n + per_row - 1) // per_row
    gutter = 62.0
    ml, mr = 30.0, 24.0
    W = ml + gutter + per_row * tile + mr
    seq_h, tgap, block_gap = 15.0, 1.5, 20.0
    block_h = 12 + seq_h + len(names) * (track + tgap) + block_gap
    top = 74.0 if letter else 54.0
    H = top + rows * block_h - block_gap + 14

    c = Canvas(W, H, font)
    if letter:
        c.text(22, 32, letter, 13, INK, "start", "600")
    title = label or acc
    tx = ml + (26 if letter else 0)
    c.text(tx, 34, title, 12.5, INK, "start", "600")
    c.text(tx + len(title) * 7.8 + 16, 34,
           f"{acc} · {n} aa" + (f" · sample {sample}" if sample else " · all samples"),
           9.5, INK_MUTED, "start")
    # the key describes the union track only
    kx = W - mr - 4 - (len(uramp) + 1) * 21
    c.text(kx - 8, 34, f"digests covering ({union_name})", 8.4, INK_MUTED, "end")
    for k, fill in enumerate([empty] + list(uramp)):
        c.rect(kx, 26, 9, 9, fill=fill, rx=radius,
               stroke=AXIS if fill in (EMPTY, "#ffffff") else "none", sw=0.7)
        c.text(kx + 11, 34, str(k), 8.0, INK_MUTED, "start")
        kx += 21

    x0 = ml + gutter
    for r in range(rows):
        by = top + r * block_h
        lo, hi = r * per_row, min(r * per_row + per_row, n)
        c.text(x0 - 8, by + 10, str(lo + 1), 8.4, INK_MUTED, "end")
        c.text(x0 + (hi - lo) * tile + 6, by + 10, str(hi), 8.4, INK_MUTED, "start")
        for i in range(lo, hi):
            c.text(x0 + (i - lo) * tile + tile / 2, by + 10, seq[i], 8.8,
                   INK_SECONDARY, "middle")
        for k, nm in enumerate(names):
            ty = by + 14 + k * (track + tgap)
            c.text(x0 - 8, ty + track * 0.78, nm, 8.6, INK, "end")
            for i in range(lo, hi):
                if nm == union_name:
                    fill = empty if ndig[i] == 0 else \
                        uramp[min(ndig[i], len(uramp)) - 1]
                else:
                    fill = colour[nm] if masks[nm][i] else empty
                c.rect(x0 + (i - lo) * tile, ty, tile, track, fill=fill,
                       rx=radius,
                       stroke=AXIS if (fill in (EMPTY, "#ffffff")) else "none",
                       sw=0.6)
    write(c, out)
    return W, H


def draw_protein(acc, seq, peps, out, font, letter, label=None,
                 iso_acc=None, iso_seq=None, iso_peps=None, max_tiers=0,
                 sample=None):
    n = len(seq)
    names = [p for p in ("GluC", "LysC", "Trypsin") if p in peps] or sorted(peps)
    colour = assign(names + ["All"])
    masks = {p: mark(seq, peps[p])[0] for p in names}
    depth = [sum(masks[p][i] for p in names) for i in range(n)]

    TIER_H, TIER_GAP, BLOCK_GAP = 6.0, 2.5, 20.0
    W = 1010
    bx0, bx1 = 150.0, W - 34.0
    pitch = (bx1 - bx0) / n
    min_gap = max(1, round(2.0 / pitch))          # keep 2px of surface between peptides

    tiers = {name: pack(locate(seq, peps[name]), min_gap) for name in names}
    shown = {name: (t[:max_tiers] if max_tiers else t) for name, t in tiers.items()}

    seq_y = 82.0
    top = seq_y + 14
    y = top
    block_y = {}
    for name in names:
        block_y[name] = y
        y += len(shown[name]) * (TIER_H + TIER_GAP) + BLOCK_GAP
    dy = y - BLOCK_GAP + 26                       # cross-digest depth track
    iso_rows = 1 if (iso_seq and iso_peps) else 0
    bottom = (dy + 96) if iso_rows else (dy + 42)
    H = bottom + 14
    c = Canvas(W, H, font)
    if letter:
        c.text(22, 30, letter, 13, INK, "start", "600")

    def X(r):
        return bx0 + r * pitch

    title = label or acc
    c.text(bx0, 40, title, 12, INK, "start", "600")
    sub = f"{acc} · {n} aa" + (f" · sample {sample}" if sample else " · all samples")
    c.text(bx0 + len(title) * 7 + 12, 40, sub, 9.5, INK_MUTED, "start")

    ruler_y = 68.0
    tick = 10 if n <= 200 else (50 if n <= 800 else 100)
    for r in range(tick, n + 1, tick):
        c.line(X(r), ruler_y - 5, X(r), ruler_y - 1, stroke=AXIS, sw=1)
        c.text(X(r), ruler_y - 8, str(r), 8.2, INK_MUTED, "middle")

    show_seq = pitch >= 4.5
    if show_seq:
        for i, aa in enumerate(seq):
            c.text(X(i) + pitch / 2, seq_y, aa, 8.6, INK_SECONDARY, "middle")
    else:
        c.rect(bx0, seq_y - 9, bx1 - bx0, 11, fill=GRID, stroke=AXIS, sw=0.8, rx=3)

    for name in names:
        by = block_y[name]
        rows = shown[name]
        stack_h = len(rows) * (TIER_H + TIER_GAP)
        cov = sum(masks[name]) / n

        for r, row in enumerate(rows):
            ry = by + r * (TIER_H + TIER_GAP)
            for s, e in row:
                c.rect(X(s) + 0.5, ry, max((e - s) * pitch - 1.0, 1.2), TIER_H,
                       fill=colour[name], rx=2)

        mid = by + stack_h / 2 - 4
        c.rect(bx0 - 13, mid - 4, 9, 9, fill=colour[name], rx=2)
        c.text(bx0 - 19, mid + 4, name, 9.5, INK, "end")
        c.text(bx0 - 19, mid + 16, f"{len(peps[name])} peptides", 8.4,
               INK_MUTED, "end")
        c.text(bx0 - 19, mid + 27, f"{cov:.0%} covered", 8.4, INK_MUTED, "end")
        if max_tiers and len(tiers[name]) > max_tiers:
            c.text(bx1, by + stack_h + 12,
                   f"+{len(tiers[name]) - max_tiers} more tiers", 8.0,
                   INK_MUTED, "end")

    c.line(bx0, dy - 13, bx1, dy - 13, stroke=GRID, sw=1)
    c.text(bx0 - 13, dy + 11, "covered by", 9.5, INK, "end")
    i = 0
    while i < n:
        j = i
        while j < n and depth[j] == depth[i]:
            j += 1
        if depth[i]:
            c.rect(X(i), dy, X(j) - X(i), 15,
                   fill=DEPTH[min(depth[i], len(DEPTH)) - 1])
        i = j
    c.rect(bx0, dy, bx1 - bx0, 15, fill="none", stroke=AXIS, sw=0.8, rx=2)
    lx = bx0
    for k, fill in enumerate(DEPTH[:len(names)], start=1):
        c.rect(lx, dy + 25, 10, 10, fill=fill, rx=2)
        c.text(lx + 13, dy + 34, str(k), 8.6, INK, "start")
        lx += 26
    c.text(lx + 2, dy + 34, "digests", 8.6, INK_MUTED, "start")

    # optional isoform track, drawn in canonical coordinates
    if iso_rows:
        m = deletion_map(seq, iso_seq)
        iy = dy + 52
        c.text(bx0 - 13, iy + 11, iso_acc, 9.5, INK, "end")
        c.text(bx0 - 13, iy + 23, f"{len(iso_seq)} aa", 8.4, INK_MUTED, "end")
        c.rect(bx0, iy, bx1 - bx0, 14, fill=STRIP_FILL, fo=0.7, rx=3)
        if m:
            p, d = m
            c.rect(X(p), iy, d * pitch, 14, fill=INK_MUTED, fo=0.18)
            c.text(X(p + d / 2), iy + 27, f"Δ{p + 1}–{p + d}", 8.4,
                   INK_MUTED, "middle")
            for name, pset in iso_peps.items():
                for pep in sorted(pset):
                    a = iso_seq.find(pep)
                    if a == -1:
                        continue
                    segs = iso_segments(a, a + len(pep), p, d, len(iso_seq))
                    for x0, x1 in segs:
                        c.rect(X(x0) + 0.6, iy, max((x1 - x0) * pitch - 1.2, 1.0),
                               14, fill=colour.get(name, UNION), fo=0.9, rx=2.5)
                    if len(segs) > 1:
                        c.line(X(segs[0][1]), iy + 7, X(segs[1][0]), iy + 7,
                               stroke=colour.get(name, UNION), sw=1.4)
        else:
            c.text(bx0 + 6, iy + 11, "not a single-deletion isoform", 8.6,
                   INK_MUTED, "start")
    write(c, out)
    return W, H


def write(c, out):
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())


# -------------------------------------------------------------------- main ---
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reports", nargs="+")
    ap.add_argument("--fasta", required=True, help="the FASTA the search used")
    ap.add_argument("--rank", action="store_true", help="rank-coverage curves")
    ap.add_argument("--protein-q", type=float, default=0.01,
                    help="protein-group q-value cutoff for --rank, as in prep_counts.py")
    ap.add_argument("--show-samples", action="store_true",
                    help="with --by-sample, draw every sample faintly behind its median")
    ap.add_argument("--paired", action="store_true",
                    help="coverage vs protein rank, every series on the protein sorted by --against")
    ap.add_argument("--against", default="Trypsin",
                    help="series the paired plot sorts on")
    ap.add_argument("--all-proteins", action="store_true",
                    help="beeswarm: don't restrict to proteins every series measured")
    ap.add_argument("--beeswarm", action="store_true",
                    help="per-protein coverage as a box plus point cloud")
    ap.add_argument("--by-sample", action="store_true",
                    help="median curve over per-sample ranks")
    ap.add_argument("--counts", default=None, metavar="COUNTS_CSV",
                    help="cross-check proteoform totals against the count-panel table")
    ap.add_argument("--width", type=float, default=1215.0,
                    help="figure width for --rank")
    ap.add_argument("--height", type=float, default=430.0)
    ap.add_argument("--text-scale", type=float, default=1.7,
                    help="type size for --rank")
    ap.add_argument("--protein", default=None, help="accession for a single-protein map")
    ap.add_argument("--isoform", default=None, help="a second accession to draw below")
    ap.add_argument("--label", default=None, help="gene name shown on the panel")
    ap.add_argument("--representative", choices=("first", "longest"), default="first")
    ap.add_argument("--sample", default=None,
                    help="restrict to one patient")
    ap.add_argument("--sample-regex", default=None,
                    help="regex on the run name; group 1 is the sample id")
    ap.add_argument("--survey", action="store_true",
                    help="with --protein: print coverage per digest per sample and stop")
    ap.add_argument("--layout", choices=("pileup", "grid", "stack", "line"),
                    default="pileup",
                    help="single-protein map layout")
    ap.add_argument("--per-row", type=int, default=30,
                    help="residues per row in --layout grid")
    ap.add_argument("--grid-cols", type=int, default=2,
                    help="blocks per row in --layout grid")
    ap.add_argument("--tile", type=float, default=15.0,
                    help="tile size in px for --layout grid")
    ap.add_argument("--pad", type=float, default=0.0,
                    help="gap between residue tiles (0 = a continuous band)")
    ap.add_argument("--group-gap", type=float, default=0.0,
                    help="extra gap every --group residues (0 = none)")
    ap.add_argument("--strip", type=float, default=15.0,
                    help="strip height for --layout line")
    ap.add_argument("--tile-radius", type=float, default=0.0,
                    help="tile corner radius (0 = square)")
    ap.add_argument("--empty", default="black", choices=("black", "surface"),
                    help="colour for a residue no peptide reached")
    ap.add_argument("--max-tiers", type=int, default=20,
                    help="cap the pileup rows per digest (0 = no cap)")
    ap.add_argument("--precursor-q", type=float, default=0.01)
    ap.add_argument("--letter", default="")
    ap.add_argument("--font", default=FONT)
    ap.add_argument("--out", default="coverage.svg")
    args = ap.parse_args(argv)

    if not (args.rank or args.paired or args.beeswarm or args.protein):
        sys.exit("choose --rank, --paired, --beeswarm or --protein ACC")

    paths = []
    for pat in args.reports:
        paths.extend(sorted(glob.glob(pat)) or [pat])
    from lib_report import one_search_per_run
    one_search_per_run(paths)
    seqs = read_fasta(args.fasta)
    print(f"{len(seqs):,} sequences in {os.path.basename(args.fasta)}")
    if args.survey:
        if not (args.protein and args.sample_regex):
            sys.exit("--survey needs --protein and --sample-regex")
        acc = args.protein
        if acc not in seqs:
            sys.exit(f"{acc} not in {args.fasta}")
        seq = seqs[acc]
        print("reading reports:")
        per = {}
        samples = set()
        for path in paths:
            g, tally = collect([path], args.precursor_q, accession=acc,
                               sample_regex=args.sample_regex)
            # re-collect the small per-protein slice instead of re-reading the file
            for (p, s), _ in tally.items():
                samples.add(s)
        for s in sorted(x for x in samples if x):
            for path in paths:
                g, _ = collect([path], args.precursor_q, accession=acc,
                               sample_regex=args.sample_regex, sample=s)
                for grp, byp in g.items():
                    if acc not in {a.split("|")[1] if a.count("|") >= 2 else a
                                   for a in re.split(r"[;,]", grp)}:
                        continue
                    for p, pset in byp.items():
                        per.setdefault((s, p), set()).update(pset)
        prots = sorted({p for _, p in per})
        print(f"\n{acc} ({len(seq)} aa) coverage per sample:")
        head = "".join(f"{p:>22}" for p in prots)
        print(f"  {'sample':<8}{head}{'union':>22}   digests")
        for s in sorted(x for x in samples if x):
            row, allp = "", set()
            for p in prots:
                pset = per.get((s, p), set())
                allp |= pset
                m, _ = mark(seq, pset)
                row += (f"{len(pset):>4d} pep {sum(m) / len(seq):>6.1%}"
                        if pset else f"{'-':>15}") + "       "
            m, _ = mark(seq, allp)
            k = sum(1 for p in prots if per.get((s, p)))
            print(f"  {s:<8}{row}{len(allp):>4d} pep {sum(m) / len(seq):>6.1%}"
                  f"      {k}/{len(prots)}")
        return 0

    print("reading reports:")
    wants_cov = args.rank or args.paired or args.beeswarm
    by_sample = wants_cov and args.by_sample and not args.sample
    groups, _ = collect(paths, args.precursor_q,
                        accession=args.protein if args.protein else None,
                        sample_regex=args.sample_regex, sample=args.sample,
                        split_samples=by_sample,
                        protein_q=args.protein_q if wants_cov else None)
    if by_sample:
        print(f"{len(groups):,} samples")
    else:
        print(f"{len(groups):,} protein groups with peptides"
              + (f" for {args.protein}" if args.protein else "")
              + (f" in sample {args.sample}" if args.sample else ""))

    if args.protein:
        acc = args.protein
        if acc not in seqs:
            sys.exit(f"{acc} not in {args.fasta}")
        peps = defaultdict(set)
        for grp, byp in groups.items():
            members = {a.split("|")[1] if a.count("|") >= 2 else a
                       for a in re.split(r"[;,]", grp)}
            if acc in members:
                for p, s in byp.items():
                    peps[p] |= s
        if not peps:
            sys.exit(f"no peptides assigned to a group containing {acc}")
        seq = seqs[acc]
        iso_seq = iso_peps = None
        if args.isoform:
            if args.isoform not in seqs:
                sys.exit(f"{args.isoform} not in {args.fasta}")
            iso_seq = seqs[args.isoform]
            iso_peps = defaultdict(set)
            for grp, byp in groups.items():
                members = {a.strip() for a in re.split(r"[;,]", grp)}
                if args.isoform in members:
                    for p, s in byp.items():
                        iso_peps[p] |= s
        if args.layout in ("grid", "stack", "line"):
            order = [p for p in ("GluC", "LysC", "Trypsin") if p in peps] or \
                sorted(peps)
            series = {p: peps[p] for p in order}
            series["All"] = set().union(*peps.values())
            emptycol = EMPTY_DARK if args.empty == "black" else EMPTY
            if args.layout == "line":
                W, H = draw_line(acc, seq, series, args.out, args.font,
                                 args.letter, args.label, args.sample,
                                 strip=args.strip, empty=emptycol)
            elif args.layout == "stack":
                W, H = draw_stack(acc, seq, series, args.out, args.font,
                                  args.letter, args.label, args.sample,
                                  per_row=args.per_row, tile=args.tile,
                                  empty=emptycol, radius=args.tile_radius)
            else:
                W, H = draw_grid(acc, seq, series, args.out, args.font,
                                 args.letter, args.label, args.sample,
                                 args.per_row, args.grid_cols, tile=args.tile,
                                 pad=args.pad, gap=args.group_gap,
                                 empty=emptycol, radius=args.tile_radius)
        else:
            W, H = draw_protein(acc, seq, peps, args.out, args.font, args.letter,
                                args.label, args.isoform, iso_seq, iso_peps,
                                args.max_tiers, args.sample)
        print(f"\nwrote {args.out}  ({W}x{H}px)")
        allp = set().union(*peps.values())
        cov, found = mark(seq, allp)
        print(f"  {args.label or acc}  {acc}  {len(seq)} aa")
        for p in sorted(peps):
            m, f = mark(seq, peps[p])
            print(f"    {p:<8s} {len(peps[p]):>4d} peptides  "
                  f"{sum(m) / len(seq):6.1%} coverage")
        print(f"    {'union':<8s} {len(allp):>4d} peptides  "
              f"{sum(cov) / len(seq):6.1%} coverage")
        if found != len(allp):
            print(f"    {len(allp) - found} peptide(s) not found in {acc} "
                  f"(they belong to another member of the group)")
        if iso_peps:
            iso_all = set().union(*iso_peps.values())
            m, f = mark(iso_seq, iso_all)
            print(f"  {args.isoform}  {len(iso_seq)} aa  {len(iso_all)} peptides  "
                  f"{sum(m) / len(iso_seq):.1%} coverage")
            junction = [p for p in iso_all if p not in seq]
            for p in sorted(junction):
                who = sorted(k for k, v in iso_peps.items() if p in v)
                print(f"    isoform-specific: {p}  ({', '.join(who)})")
        return 0

    # rank mode
    def cov_one(gs):
        """One sample's per-protein coverage: {group: {series: fraction}}."""
        out = {}
        for grp, byp in gs.items():
            rep = representative(grp, seqs, args.representative)
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
            rep = representative(grp, seqs, args.representative)
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

    if args.paired or args.beeswarm:
        acc = defaultdict(lambda: defaultdict(list))
        src = groups if by_sample else {"pooled": groups}
        for s_, gs in sorted(src.items()):
            for grp, row in cov_one(gs).items():
                for k, v in row.items():
                    acc[grp][k].append(v)
        prot = {g: {k: med(v) for k, v in r.items()} for g, r in acc.items()}
        print(f"{len(prot):,} proteins with coverage in >=1 sample")
        if args.paired:
            W, H = draw_paired(prot, args.out, args.font, args.letter,
                               args.width, args.height, args.text_scale,
                               args.against)
        else:
            W, H = draw_beeswarm(prot, args.out, args.font, args.letter,
                                 args.width, args.height, args.text_scale,
                                 not args.all_proteins)
        print(f"\nwrote {args.out}  ({W}x{H}px)")
        return 0

    if by_sample:
        per = {s_: rank_one(gs)[0] for s_, gs in sorted(groups.items())}
        names = sorted({k for v in per.values() for k in v})
        keep = None
        if args.counts:
            import csv as _csv
            keep = defaultdict(set)
            with open(args.counts, newline="") as fh:
                for r in _csv.DictReader(fh):
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
        print(f"\nper-sample rank curves, median across samples:")
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
        W, H = draw_rank(series, args.out, args.font, args.letter,
                         args.width, args.height, args.text_scale,
                         sample=None, curves=curves,
                         show_curves=args.show_samples)
        print(f"\nwrote {args.out}  ({W}x{H}px)")
        for k in sorted(series, key=lambda k: -len(series[k])):
            half = [sum(1 for x in v if x >= 0.5) for v in curves[k]]
            half.sort()
            m = (half[len(half) // 2] if len(half) % 2 else
                 (half[len(half) // 2 - 1] + half[len(half) // 2]) / 2)
            print(f"  {k:<8s} median {len(series[k]):>7,} proteoforms   "
                  f"median >=50% {m:>7,.0f}")
        return 0

    series = defaultdict(list)
    missing = unmapped = 0
    for grp, byp in groups.items():
        rep = representative(grp, seqs, args.representative)
        if rep is None:
            missing += 1
            continue
        seq = seqs[rep]
        if not seq:
            continue
        allp = set()
        for p, s in byp.items():
            m, f = mark(seq, s)
            unmapped += len(s) - f
            series[p].append(sum(m) / len(seq))
            allp |= s
        m, _ = mark(seq, allp)
        series["All"].append(sum(m) / len(seq))
    for k in series:
        series[k].sort(reverse=True)
    if missing:
        print(f"{missing:,} group(s) had no member in the FASTA -- skipped")
    if unmapped:
        print(f"{unmapped:,} peptide-group pairs did not occur in the "
              f"representative sequence")
    if args.counts:
        cross_check(series, args.counts, args.sample)
    W, H = draw_rank(series, args.out, args.font, args.letter,
                     args.width, args.height, args.text_scale, args.sample)
    print(f"\nwrote {args.out}  ({W}x{H}px)")
    print("\n  proteoform counts are comparable across series; medians are not,\n"
          "  because each series ranks a different number of proteoforms")
    for k in sorted(series, key=lambda k: -len(series[k])):
        v = series[k]
        mid = v[len(v) // 2] if v else 0   # not `med`: that is the module function
        print(f"  {k:<8s} n = {len(v):>7,}   >=50% {sum(1 for x in v if x >= 0.5):>7,}"
              f"   >=80% {sum(1 for x in v if x >= 0.8):>6,}"
              f"   (median of its own list {mid:.1%})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
