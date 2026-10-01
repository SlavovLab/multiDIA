#!/usr/bin/env python3
"""Differential abundance between conditions, per digest and across digests.

    # once: protein x sample quantity matrices from the reports
    python3 extra_differential.py extract 'data/search/*-60min-Phospho.parquet' --outdir quant

    # then: test, as often as you like
    python3 extra_differential.py test quant --metadata data/metadata.xlsx \
        --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta --outdir da
"""

import argparse
import collections
import csv
import glob
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_svg import Canvas                                  # noqa: E402
from lib_fasta import gene_map                                     # noqa: E402
from lib_palette import (AXIS, FONT, GRID, INK, INK_MUTED, INK_SECONDARY,
                     STRIP_FILL, UNION, assign)                     # noqa: E402

ORDER = ["GluC", "LysC", "Trypsin"]
from lib_report import sample_of                                       # noqa: E402


def extract(paths, outdir, protein_q):
    import lib_report as rp
    os.makedirs(outdir, exist_ok=True)
    for r in rp.open_reports(paths, order=ORDER):
        dig = r.protease
        q = collections.defaultdict(dict)
        unresolved = collections.Counter()
        for b in r.batches(["run", "protein_groups", "protein_q_run",
                            "protein_qty"]):
            fn, pg = b["run"], b["protein_groups"]
            qv, qt = b["protein_q_run"], b["protein_qty"]
            for i in range(b["_n"]):
                if qv[i] is None or qv[i] > protein_q or not qt[i] or qt[i] <= 0:
                    continue
                s_ = sample_of(fn[i])
                if s_:
                    q[pg[i]][s_] = qt[i]
                else:
                    unresolved[fn[i]] += 1
        samples = sorted({s for d in q.values() for s in d})
        # Warn: an unparsed run name is a whole patient missing.
        if unresolved:
            print(f"  WARNING {dig}: {sum(unresolved.values()):,} rows in "
                  f"{len(unresolved)} run(s) have no patient id and were dropped")
            for run, n_ in sorted(unresolved.items()):
                print(f"      {run}  {n_:,} rows")
        out = os.path.join(outdir, f"{dig}.csv")
        with open(out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["protein_group"] + samples)
            for g in sorted(q):
                w.writerow([g] + [q[g].get(s, "") for s in samples])
        print(f"  {dig:<8} {len(q):>6,} groups x {len(samples)} samples -> {out}")


def load_matrix(path):
    rows = list(csv.DictReader(open(path)))
    samples = [c for c in rows[0] if c != "protein_group"]
    mat = {r["protein_group"]: {s: float(r[s]) for s in samples if r[s]} for r in rows}
    return mat, samples


def label_for(group, genes):
    names = []
    for a in group.split(";"):
        g = genes.get(a) or genes.get(a.split("-")[0])
        if g and g not in names:
            names.append(g)
    return ";".join(names[:2]) if names else group.split(";")[0]


def qc_runs(mat, samples, min_frac):
    """Drop runs whose quantified-protein count is far below the digest median."""
    n = {s: sum(1 for g in mat if s in mat[g]) for s in samples}
    import statistics
    med = statistics.median(n.values())
    keep = [s for s in samples if n[s] >= med * min_frac]
    dropped = {s: n[s] for s in samples if s not in keep}
    return keep, dropped, med


def outlier_runs(logs, samples, max_mad):
    """Runs whose profile sits far from the rest of the digest."""
    import numpy as np
    if max_mad <= 0 or len(samples) < 4:
        return {}
    M = np.array([[logs[g][s] for s in samples] for g in logs
                  if all(s in logs[g] for s in samples)])
    if len(M) < 50:
        return {}
    M = M - np.median(M, axis=0, keepdims=True)
    C = np.corrcoef(M.T)
    mc = np.array([np.median(np.delete(C[i], i)) for i in range(len(samples))])
    med = np.median(mc)
    mad = np.median(np.abs(mc - med)) or 1e-9
    return {s: (float(v), float((med - v) / mad))
            for s, v in zip(samples, mc) if (med - v) / mad > max_mad}


def normalise(mat, samples):
    """log2, then median-centre each run onto the grand median."""
    import statistics
    logs = {g: {s: math.log2(v) for s, v in d.items() if s in samples}
            for g, d in mat.items()}
    med = {s: statistics.median([logs[g][s] for g in logs if s in logs[g]])
           for s in samples}
    grand = statistics.median(med.values())
    for g in logs:
        for s in logs[g]:
            logs[g][s] += grand - med[s]
    return logs, med, grand


def test_digest(logs, groups, min_per_group):
    """Welch's t per protein, BH within the digest."""
    import numpy as np
    from scipy import stats
    names, lfc, pv, na, nb = [], [], [], [], []
    ga, gb = groups          # (case_samples, control_samples)
    for g, d in logs.items():
        a = [d[s] for s in ga if s in d]
        b = [d[s] for s in gb if s in d]
        if len(a) < min_per_group or len(b) < min_per_group:
            continue
        t = stats.ttest_ind(a, b, equal_var=False)
        if not np.isfinite(t.pvalue):
            continue
        names.append(g)
        lfc.append(float(np.mean(a) - np.mean(b)))
        pv.append(float(t.pvalue))
        na.append(len(a))
        nb.append(len(b))
    if not names:
        return []
    q = stats.false_discovery_control(np.array(pv), method="bh")
    return [{"group": g, "log2fc": f, "p": p, "q": float(qq), "n_case": x,
             "n_ctrl": y}
            for g, f, p, qq, x, y in zip(names, lfc, pv, q, na, nb)]


def trigamma_inverse(x):
    """Solve trigamma(y) = x for y (limma's trigammaInverse)."""
    from scipy.special import polygamma
    if x > 1e7:
        return 1.0 / math.sqrt(x)
    if x < 1e-6:
        return 1.0 / x
    y = 0.5 + 1.0 / x
    for _ in range(60):
        tri = polygamma(1, y)
        d = tri * (1 - tri / x) / polygamma(2, y)
        y += d
        if -d / y < 1e-8:
            break
    return y


def ebayes(s2, df):
    """Smyth (2004) method of moments for the prior variance."""
    import numpy as np
    from scipy.special import digamma, polygamma
    s2 = np.asarray(s2, float)
    ok = s2 > 0
    z = np.log(s2[ok])
    d = np.full(z.shape, float(df))
    e = z - digamma(d / 2) + np.log(d / 2)
    ebar = e.mean()
    G = len(e)
    var_e = (G / (G - 1.0)) * ((e - ebar) ** 2).mean() - polygamma(1, d / 2).mean()
    if var_e > 0:
        d0 = 2.0 * trigamma_inverse(var_e)
        s02 = math.exp(ebar + digamma(d0 / 2) - math.log(d0 / 2))
    else:
        d0 = float("inf")
        s02 = math.exp(ebar)
    return d0, s02


def test_digest_moderated(logs, groups, min_per_group):
    """Two-group limma-style moderated t."""
    import numpy as np
    from scipy import stats
    ga, gb = groups
    names, diff, s2, dfs, na, nb = [], [], [], [], [], []
    for g, d in logs.items():
        a = [d[s] for s in ga if s in d]
        b = [d[s] for s in gb if s in d]
        if len(a) < min_per_group or len(b) < min_per_group:
            continue
        n1, n2 = len(a), len(b)
        va = np.var(a, ddof=1) if n1 > 1 else 0.0
        vb = np.var(b, ddof=1) if n2 > 1 else 0.0
        pooled = ((n1 - 1) * va + (n2 - 1) * vb) / (n1 + n2 - 2)
        names.append(g)
        diff.append(float(np.mean(a) - np.mean(b)))
        s2.append(float(pooled))
        dfs.append(n1 + n2 - 2)
        na.append(n1)
        nb.append(n2)
    if not names:
        return []
    s2 = np.array(s2)
    diff = np.array(diff)
    dfs = np.array(dfs, float)
    # one prior per residual-df stratum, since missing values give a few strata
    d0s = np.zeros_like(s2)
    s02s = np.zeros_like(s2)
    for d in np.unique(dfs):
        m = dfs == d
        if m.sum() < 20:
            d0s[m], s02s[m] = 0.0, s2[m].mean() if m.sum() else 1.0
            continue
        d0, s02 = ebayes(s2[m], d)
        d0s[m] = 0.0 if not np.isfinite(d0) else d0
        s02s[m] = s02
    post = (d0s * s02s + dfs * s2) / (d0s + dfs)
    se = np.sqrt(post * (1.0 / np.array(na) + 1.0 / np.array(nb)))
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(se > 0, diff / se, 0.0)
    dtot = dfs + d0s
    p = 2 * stats.t.sf(np.abs(t), dtot)
    q = stats.false_discovery_control(p, method="bh")
    return [{"group": g, "log2fc": float(f), "p": float(pp), "q": float(qq),
             "n_case": int(x), "n_ctrl": int(y), "t": float(tt),
             "df": float(dd)}
            for g, f, pp, qq, x, y, tt, dd
            in zip(names, diff, p, q, na, nb, t, dtot)]


def combine_patients(all_logs, min_digests):
    """One value per protein per patient, pooling the digests."""
    import numpy as np
    per = collections.defaultdict(lambda: collections.defaultdict(list))
    seen = collections.defaultdict(set)
    for dig, logs in all_logs.items():
        for g, d in logs.items():
            if len(d) < 2:
                continue
            mu = float(np.mean(list(d.values())))
            seen[g].add(dig)
            for s, v in d.items():
                per[g][s].append(v - mu)
    return ({g: {s: float(np.mean(v)) for s, v in d.items()}
             for g, d in per.items() if len(seen[g]) >= min_digests},
            {g: sorted(v) for g, v in seen.items()})


def nice(hi, target=5):
    if hi <= 0:
        return [0]
    raw = hi / target
    mag = 10 ** math.floor(math.log10(raw))
    step = next((m * mag for m in (1, 2, 2.5, 5, 10) if raw <= m * mag), 10 * mag)
    return [i * step for i in range(int(hi / step) + 2) if i * step <= hi * 1.001]


def volcano(results, genes, outdir, font, qcut, letters="abc"):
    colour = assign(ORDER + ["All"])
    W, H = 1020, 340
    c = Canvas(W, H, font)
    ml, gap = 62.0, 46.0
    pw = (W - ml - 24 - gap * (len(ORDER) - 1)) / len(ORDER)
    ph = H - 96
    allf = [abs(r["log2fc"]) for res in results.values() for r in res]
    xlim = min(max(allf) if allf else 2, 6)
    ymax = max([-math.log10(r["p"]) for res in results.values() for r in res] or [1])
    for pi, dig in enumerate(ORDER):
        res = results.get(dig, [])
        x0 = ml + pi * (pw + gap)
        y0 = 56.0
        c.text(x0 - 40 if pi == 0 else x0 - 34, 32, letters[pi], 13, INK,
               "start", "600")
        c.text(x0, 46, dig, 11, INK, "start", "600")
        n_sig = sum(1 for r in res if r["q"] < qcut)
        c.text(x0 + len(dig) * 7 + 10, 46,
               f"{len(res):,} tested · {n_sig:,} at q<{qcut:g}", 8.8,
               INK_MUTED, "start")

        def X(v):
            return x0 + (max(-xlim, min(xlim, v)) + xlim) / (2 * xlim) * pw

        def Y(v):
            return y0 + ph - min(v, ymax) / ymax * ph

        for t in nice(ymax, 4):
            c.line(x0, Y(t), x0 + pw, Y(t), stroke=GRID, sw=1)
            if pi == 0:
                c.text(x0 - 8, Y(t) + 3.5, f"{t:.0f}", 9, INK_MUTED, "end")
        c.line(X(0), y0, X(0), y0 + ph, stroke=AXIS, sw=1)
        for v in (-xlim, 0, xlim):
            c.text(X(v), y0 + ph + 16, f"{v:.0f}", 9, INK_MUTED, "middle")
        c.line(x0, y0 + ph, x0 + pw, y0 + ph, stroke=AXIS, sw=1)

        for r in sorted(res, key=lambda r: r["q"], reverse=True):
            sig = r["q"] < qcut
            c.add(f'<circle cx="{X(r["log2fc"]):.1f}" cy="{Y(-math.log10(r["p"])):.1f}" '
                  f'r="{2.6 if sig else 1.6}" fill="{colour[dig] if sig else INK_MUTED}" '
                  f'fill-opacity="{0.9 if sig else 0.22}"/>')
        # label the strongest by p even when nothing clears FDR
        top = sorted(res, key=lambda r: r["p"])[:6]
        for k, r in enumerate(top):
            c.add(f'<circle cx="{X(r["log2fc"]):.1f}" '
                  f'cy="{Y(-math.log10(r["p"])):.1f}" r="3" fill="none" '
                  f'stroke="{colour[dig]}" stroke-width="1.6"/>')
            c.text(X(r["log2fc"]) + (7 if r["log2fc"] > 0 else -7),
                   Y(-math.log10(r["p"])) + 3,
                   label_for(r["group"], genes), 8,
                   INK_SECONDARY, "start" if r["log2fc"] > 0 else "end")
    c.text(W / 2, H - 14, "log2 fold change   (LBD − Control)", 10, INK, "middle")
    c.add(f'<g transform="translate(20 {56 + ph / 2:.1f}) rotate(-90)">'
          f'<text x="0" y="0" font-size="10" text-anchor="middle" fill="{INK}">'
          f'−log10 p</text></g>')
    p = os.path.join(outdir, "volcano.svg")
    with open(p, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {p}")


def concordance(results, genes, outdir, font, qcut, letter="d"):
    """Do the digests agree, on the proteins all of them measured?"""
    per = {d: {r["group"]: r for r in res} for d, res in results.items()}
    shared = set.intersection(*[set(v) for v in per.values()]) if per else set()
    pairs = [("Trypsin", "LysC"), ("Trypsin", "GluC"), ("LysC", "GluC")]
    colour = assign(ORDER + ["All"])
    W, H = 1020, 366
    c = Canvas(W, H, font)
    ml, gap = 66.0, 52.0
    pw = (W - ml - 30 - gap * 2) / 3
    ph = H - 132
    c.text(24, 32, letter, 13, INK, "start", "600")
    head = "Do the digests agree?"
    c.text(ml, 32, head, 11.5, INK, "start", "600")
    c.text(ml + len(head) * 7.2 + 16, 32,
           f"{len(shared):,} protein groups tested in all three", 9, INK_MUTED,
           "start")
    lim = 3.0
    for pi, (a, b) in enumerate(pairs):
        x0 = ml + pi * (pw + gap)
        y0 = 56.0
        xs = [(per[a][g]["log2fc"], per[b][g]["log2fc"]) for g in shared]
        n = len(xs)
        r = 0.0
        if n > 2:
            ma = sum(v[0] for v in xs) / n
            mb = sum(v[1] for v in xs) / n
            sa = math.sqrt(sum((v[0] - ma) ** 2 for v in xs))
            sb = math.sqrt(sum((v[1] - mb) ** 2 for v in xs))
            if sa and sb:
                r = sum((v[0] - ma) * (v[1] - mb) for v in xs) / (sa * sb)
        agree = sum(1 for u, v in xs if u * v > 0) / n if n else 0

        def X(v):
            return x0 + (max(-lim, min(lim, v)) + lim) / (2 * lim) * pw

        def Y(v):
            return y0 + ph - (max(-lim, min(lim, v)) + lim) / (2 * lim) * ph

        for t in (-2, 0, 2):
            c.line(X(t), y0, X(t), y0 + ph, stroke=GRID, sw=1)
            c.line(x0, Y(t), x0 + pw, Y(t), stroke=GRID, sw=1)
            c.text(X(t), y0 + ph + 16, f"{t:g}", 9, INK_MUTED, "middle")
            if pi == 0:
                c.text(x0 - 8, Y(t) + 3.5, f"{t:g}", 9, INK_MUTED, "end")
        c.add(f'<path d="M {X(-lim):.1f} {Y(-lim):.1f} L {X(lim):.1f} {Y(lim):.1f}" '
              f'stroke="{AXIS}" stroke-width="1" fill="none" stroke-dasharray="0"/>')
        for u, v in xs:
            sig = (per[a][g]["q"] < qcut for g in ())
            c.add(f'<circle cx="{X(u):.1f}" cy="{Y(v):.1f}" r="1.5" '
                  f'fill="{INK_MUTED}" fill-opacity="0.25"/>')
        both = [g for g in shared
                if per[a][g]["q"] < qcut and per[b][g]["q"] < qcut]
        for g in both:
            c.add(f'<circle cx="{X(per[a][g]["log2fc"]):.1f}" '
                  f'cy="{Y(per[b][g]["log2fc"]):.1f}" r="3" fill="{UNION}" '
                  f'fill-opacity="0.9"/>')
        c.text(x0 + pw / 2, 46, f"{a} vs {b}", 10, INK, "middle", "600")
        c.text(x0, y0 + ph + 32, f"r = {r:.2f}", 9, INK_SECONDARY, "start")
        c.text(x0 + pw, y0 + ph + 32, f"{agree:.0%} same direction", 9,
               INK_SECONDARY, "end")
    c.text(W / 2, H - 10, "log2 fold change, first digest → second digest", 10,
           INK, "middle")
    p = os.path.join(outdir, "concordance.svg")
    with open(p, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {p}")
    return shared


def panel_markers(all_logs, cond, case, control, genes, wanted, outdir, font,
                  letter="e"):
    """Per-patient values for named genes, in every digest that measured them."""
    import numpy as np
    colour = assign(ORDER + ["All"])
    # gene -> digest -> {sample: centred log2}
    found = collections.OrderedDict()
    for want in wanted:
        per = {}
        for dig, logs in all_logs.items():
            hit = None
            for g, d in logs.items():
                if want in label_for(g, genes).split(";") and len(d) >= 4:
                    if hit is None or len(d) > len(all_logs[dig][hit]):
                        hit = g
            if hit:
                mu = float(np.mean(list(logs[hit].values())))
                per[dig] = {s: v - mu for s, v in logs[hit].items()}
        if per:
            found[want] = per
    if not found:
        print("  markers: none of the requested genes were quantified")
        return
    ncol = min(len(found), 4)
    nrow = (len(found) + ncol - 1) // ncol
    pw, ph, gapx, gapy = 210.0, 150.0, 34.0, 62.0
    ml, mt = 58.0, 74.0
    W = ml + ncol * pw + (ncol - 1) * gapx + 20
    H = mt + nrow * (ph + gapy) + 10
    c = Canvas(W, H, font)
    c.text(22, 32, letter, 13, INK, "start", "600")
    c.text(ml, 32, f"{case} vs {control}, per patient, per digest", 11.5, INK,
           "start", "600")
    for k, (gene, per) in enumerate(found.items()):
        x0 = ml + (k % ncol) * (pw + gapx)
        y0 = mt + (k // ncol) * (ph + gapy)
        vals = [v for d in per.values() for v in d.values()]
        lo, hi = min(vals), max(vals)
        pad = (hi - lo) * 0.12 or 0.5
        lo, hi = lo - pad, hi + pad

        def Y(v):
            return y0 + ph - (v - lo) / (hi - lo) * ph

        for t in (-2, -1, 0, 1, 2):
            if lo <= t <= hi:
                c.line(x0, Y(t), x0 + pw, Y(t), stroke=GRID, sw=1)
                c.text(x0 - 6, Y(t) + 3.5, f"{t:g}", 8.4, INK_MUTED, "end")
        c.text(x0, y0 - 10, gene, 11, INK, "start", "600")
        digs = [d for d in ORDER if d in per]
        bw = pw / len(digs)
        for di, dig in enumerate(digs):
            bx = x0 + di * bw
            for gi, (grp, label) in enumerate(((control, "C"), (case, "L"))):
                xs = bx + bw * (0.28 + 0.44 * gi)
                pts = [v for s, v in per[dig].items() if cond.get(s) == grp]
                if not pts:
                    continue
                med = float(np.median(pts))
                c.line(xs - 11, Y(med), xs + 11, Y(med), stroke=colour[dig], sw=2.4)
                for j, v in enumerate(sorted(pts)):
                    jx = xs + ((j % 3) - 1) * 3.6
                    c.add(f'<circle cx="{jx:.1f}" cy="{Y(v):.1f}" r="3" '
                          f'fill="{colour[dig]}" fill-opacity="0.8" '
                          f'stroke="#ffffff" stroke-width="1.2"/>')
                c.text(xs, y0 + ph + 13, label, 8.4, INK_MUTED, "middle")
            c.text(bx + bw / 2, y0 + ph + 27, dig, 8.8, INK, "middle")
        c.line(x0, y0 + ph, x0 + pw, y0 + ph, stroke=AXIS, sw=1)
    c.add(f'<g transform="translate(20 {mt + ph / 2:.1f}) rotate(-90)">'
          f'<text x="0" y="0" font-size="10" text-anchor="middle" fill="{INK}">'
          f'log2 abundance, centred per digest</text></g>')
    path = os.path.join(outdir, "markers.svg")
    with open(path, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {path}  ({', '.join(found)})")


def panel_test(cres, genes, panel, outdir, q_cut):
    """BH within a declared panel of genes, not across the whole proteome."""
    import numpy as np
    from scipy import stats
    by = {}
    for r in cres:
        for g in label_for(r["group"], genes).split(";"):
            if g not in by or r["p"] < by[g]["p"]:
                by[g] = r
    have = [(g, by[g]) for g in panel if g in by]
    if not have:
        print("  panel: none of the requested genes were quantified")
        return []
    p = np.array([r["p"] for _, r in have])
    q = stats.false_discovery_control(p, method="bh")
    out = [(g, r, float(qq)) for (g, r), qq in zip(have, q)]
    out.sort(key=lambda t: t[1]["p"])
    path = os.path.join(outdir, "panel.tsv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "protein_group", "log2fc", "p", "q_panel"])
        for g, r, qq in out:
            w.writerow([g, r["group"], f"{r['log2fc']:.4f}", f"{r['p']:.3e}",
                        f"{qq:.4f}"])
    n_sig = sum(1 for _, _, qq in out if qq < q_cut)
    print(f"\npre-declared panel ({len(have)} of {len(panel)} quantified), "
          f"BH within the panel: {n_sig} at q < {q_cut:g}")
    print(f"  {'gene':<10}{'log2FC':>9}{'p':>12}{'q_panel':>10}")
    for g, r, qq in out:
        print(f"  {g:<10}{r['log2fc']:>+9.2f}{r['p']:>12.2e}{qq:>10.3f}"
              + ("  *" if qq < q_cut else ""))
    print(f"  wrote {path}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("extract")
    e.add_argument("reports", nargs="+")
    e.add_argument("--outdir", default="quant")
    e.add_argument("--protein-q", type=float, default=0.01)

    t = sub.add_parser("test")
    t.add_argument("quantdir")
    t.add_argument("--metadata", required=True)
    t.add_argument("--fasta", default=None)
    t.add_argument("--outdir", default="da")
    t.add_argument("--case", default="LBD")
    t.add_argument("--control", default="Control")
    t.add_argument("--min-per-group", type=int, default=3)
    t.add_argument("--min-digests", type=int, default=2,
                   help="min digests for the combined test")
    t.add_argument("--max-outlier-mad", type=float, default=3.0,
                   help="MADs below median run correlation to drop (0 disables)")
    t.add_argument("--min-run-frac", type=float, default=0.5,
                   help="min fraction of the digest's median protein count")
    t.add_argument("--q", type=float, default=0.05)
    t.add_argument("--method", choices=("moderated", "welch"), default="moderated",
                   help="limma-style moderated t, or plain Welch t")
    t.add_argument("--panel", default="",
                   help="gene symbols (comma list or file) for a panel-wide BH")
    t.add_argument("--markers", default="TH,DDC,SLC6A3,TMEM119",
                   help="gene symbols to draw per-patient, per-digest")
    t.add_argument("--font", default=FONT)
    args = ap.parse_args(argv)

    if args.cmd == "extract":
        paths = []
        for pat in args.reports:
            paths.extend(sorted(glob.glob(pat)) or [pat])
        extract(paths, args.outdir, args.protein_q)
        return 0

    from lib_report import read_metadata
    cond = read_metadata(args.metadata)
    genes = gene_map(args.fasta) if args.fasta else {}
    os.makedirs(args.outdir, exist_ok=True)
    print(f"{args.case} vs {args.control};  "
          f"{sum(1 for c in cond.values() if c == args.case)} vs "
          f"{sum(1 for c in cond.values() if c == args.control)} patients in metadata")

    results = {}
    all_logs = {}
    for dig in ORDER:
        path = os.path.join(args.quantdir, f"{dig}.csv")
        if not os.path.exists(path):
            continue
        mat, samples = load_matrix(path)
        keep, dropped, med = qc_runs(mat, samples, args.min_run_frac)
        logs, run_med, grand = normalise(mat, keep)
        bad = outlier_runs(logs, keep, args.max_outlier_mad)
        if bad:
            keep = [s for s in keep if s not in bad]
            logs, run_med, grand = normalise(mat, keep)
        ga = [s for s in keep if cond.get(s) == args.case]
        gb = [s for s in keep if cond.get(s) == args.control]
        print(f"\n{dig}: {len(mat):,} groups, {len(samples)} runs "
              f"(median {med:,.0f} quantified)")
        if dropped:
            for s, n in dropped.items():
                print(f"  QC-dropped {s} ({cond.get(s,'?')}): {n:,} quantified, "
                      f"< {args.min_run_frac:.0%} of median")
        for s, (r, z) in sorted(bad.items(), key=lambda t: -t[1][1]):
            print(f"  QC-dropped {s} ({cond.get(s,'?')}): run-to-run r = {r:.3f}, "
                  f"{z:.1f} MAD below the median — outlier profile")
        unknown = [s for s in keep if s not in cond]
        if unknown:
            print(f"  no metadata for {unknown} — excluded")
        print(f"  {args.case} n={len(ga)} {ga}")
        print(f"  {args.control} n={len(gb)} {gb}")
        all_logs[dig] = logs
        res = (test_digest_moderated(logs, (ga, gb), args.min_per_group)
               if args.method == "moderated"
               else test_digest(logs, (ga, gb), args.min_per_group))
        results[dig] = res
        sig = [r for r in res if r["q"] < args.q]
        print(f"  tested {len(res):,} groups (>= {args.min_per_group} per group); "
              f"{len(sig):,} at q < {args.q:g}; "
              f"{sum(1 for r in res if r['p'] < 0.05):,} at raw p < 0.05")
        out = os.path.join(args.outdir, f"{dig}.tsv")
        with open(out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["protein_group", "gene", "log2fc", "p", "q",
                        "n_case", "n_ctrl"])
            for r in sorted(res, key=lambda r: r["p"]):
                w.writerow([r["group"], label_for(r["group"], genes),
                            f"{r['log2fc']:.4f}", f"{r['p']:.3e}",
                            f"{r['q']:.3e}", r["n_case"], r["n_ctrl"]])
        print(f"  wrote {out}")
        for r in sorted(res, key=lambda r: r["p"])[:5]:
            print(f"    {label_for(r['group'], genes):<16} "
                  f"log2FC {r['log2fc']:+.2f}  p {r['p']:.2e}  q {r['q']:.3f}")

    # combined: pool the digests into one value per patient
    if len(all_logs) >= 2:
        comb, seen = combine_patients(all_logs, args.min_digests)
        ga = [s for s in cond if cond[s] == args.case]
        gb = [s for s in cond if cond[s] == args.control]
        cres = (test_digest_moderated(comb, (ga, gb), args.min_per_group)
                if args.method == "moderated"
                else test_digest(comb, (ga, gb), args.min_per_group))
        sig = [r for r in cres if r["q"] < args.q]
        print(f"\ncombined across digests (>= {args.min_digests} digests per protein):")
        print(f"  tested {len(cres):,} groups; {len(sig):,} at q < {args.q:g}; "
              f"{sum(1 for r in cres if r['p'] < 0.05):,} at raw p < 0.05")
        out = os.path.join(args.outdir, "combined.tsv")
        with open(out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["protein_group", "gene", "digests", "log2fc", "p", "q",
                        "n_case", "n_ctrl"])
            for r in sorted(cres, key=lambda r: r["p"]):
                w.writerow([r["group"], label_for(r["group"], genes),
                            "+".join(seen[r["group"]]), f"{r['log2fc']:.4f}",
                            f"{r['p']:.3e}", f"{r['q']:.3e}", r["n_case"],
                            r["n_ctrl"]])
        print(f"  wrote {out}")
        for r in sorted(cres, key=lambda r: r["p"])[:10]:
            star = "  *" if r["q"] < args.q else ""
            print(f"    {label_for(r['group'], genes):<16} log2FC {r['log2fc']:+.2f}  "
                  f"p {r['p']:.2e}  q {r['q']:.3f}  "
                  f"[{'+'.join(seen[r['group']])}]{star}")
        results["Combined"] = cres
        if args.panel:
            panel = (open(args.panel).read().split()
                     if os.path.exists(args.panel)
                     else [x.strip() for x in args.panel.split(",") if x.strip()])
            panel_test(cres, genes, panel, args.outdir, args.q)
        if args.markers:
            panel_markers(all_logs, cond, args.case, args.control, genes,
                          [m.strip() for m in args.markers.split(",") if m.strip()],
                          args.outdir, args.font)

    if len(results) >= 2:
        print("\nagreement across digests:")
        digest_only = {d: results[d] for d in ORDER if d in results}
        shared = concordance(digest_only, genes, args.outdir, args.font, args.q)
        volcano(digest_only, genes, args.outdir, args.font, args.q)
        per = {d: {r["group"]: r for r in res} for d, res in digest_only.items()}
        rep = collections.Counter()
        for g in shared:
            k = sum(1 for d in per if per[d][g]["q"] < args.q)
            rep[k] += 1
        for k in sorted(rep, reverse=True):
            print(f"  significant in {k} of {len(per)} digests: {rep[k]:,}")
        combined = os.path.join(args.outdir, "shared.tsv")
        with open(combined, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["protein_group", "gene"] +
                       [f"{d}_{f}" for d in per for f in ("log2fc", "q")] +
                       ["n_sig_digests", "same_direction"])
            for g in sorted(shared,
                            key=lambda g: min(per[d][g]["q"] for d in per)):
                fcs = [per[d][g]["log2fc"] for d in per]
                w.writerow([g, label_for(g, genes)] +
                           [f"{per[d][g][f]:.4g}" for d in per
                            for f in ("log2fc", "q")] +
                           [sum(1 for d in per if per[d][g]["q"] < args.q),
                            "yes" if all(x > 0 for x in fcs) or
                                     all(x < 0 for x in fcs) else "no"])
        print(f"  wrote {combined}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
