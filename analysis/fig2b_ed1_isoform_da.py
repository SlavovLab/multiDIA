#!/usr/bin/env python3
"""Isoform-specific differential abundance by difference of differences: the Extended Data Fig. 1 volcano and the Fig. 2b model plots.

    python3 analysis/fig2b_ed1_isoform_da.py 'data/search/*-60min-Phospho.parquet' --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta --metadata data/metadata.xlsx --volcano volcano.svg
"""

import argparse
import collections
import glob
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_svg import Canvas, text_width                      # noqa: E402
from lib_palette import (AXIS, DIVERGING_HIGH, FONT, GRID, INK, INK_MUTED,
                     INK_SECONDARY, TEXT_BOOST, UNION, assign, pt_scale)     # noqa: E402
from lib_report import ORDER, open_reports, read_metadata, sample_of   # noqa: E402
from lib_fasta import gene_map, read_fasta                  # noqa: E402


# width every Figure 2 panel is authored at
FIG_PANEL_W = 1010.0


# --------------------------------------------------------------------- data ---
def read_reports(paths, seqs, genes, precursor_q=0.01, min_run_peptides=1000):
    """-> (idx, med, bygene): summed quantity per group, digest, peptide and patient; log2 run medians; groups per gene."""
    idx = collections.defaultdict(lambda: collections.defaultdict(dict))
    allq = collections.defaultdict(list)
    fields = ["run", "protein_groups", "peptide", "precursor_q", "quantity"]
    keep = {}
    for r in open_reports(paths, order=ORDER):
        path, dig = r.path, r.protease
        use = ["use_for_pep"] if r.has("use_for_pep") else []
        for b in r.batches(fields + use):
            run, grp, pep, q, fg = (b[f] for f in fields)
            used = b["use_for_pep"] if use else [True] * b["_n"]
            for k in range(b["_n"]):
                if q[k] is None or q[k] > precursor_q or not used[k]:
                    continue
                if fg[k] is None or fg[k] <= 0:
                    continue
                s = sample_of(run[k] or "")
                if not s:
                    continue
                allq[(s, dig)].append(fg[k])
                g = grp[k]
                if not g:
                    continue
                ok = keep.get(g)
                if ok is None:
                    mem = [a.strip() for a in g.split(";")]
                    gn = {genes.get(a) or genes.get(a.split("-")[0])
                          for a in mem}
                    ok = keep[g] = (len(gn) == 1 and None not in gn
                                    and all(a in seqs for a in mem))
                if ok:
                    d = idx[(g, dig)][pep[k]]
                    d[s] = d.get(s, 0.0) + fg[k]
        print(f"  scanned {os.path.basename(path)[:44]:<44s} {dig}", flush=True)
    med, bad = {}, []
    for key, v in allq.items():
        if len(v) < min_run_peptides:
            bad.append((key, len(v)))
            continue
        x = sorted(math.log2(y) for y in v)
        n = len(x)
        med[key] = x[n // 2] if n % 2 else (x[n // 2 - 1] + x[n // 2]) / 2
    if bad:
        print("  gated out " + ", ".join(f"{s}/{d} ({n} precursors)"
                                         for (s, d), n in sorted(bad)))
    bygene = collections.defaultdict(set)
    for (g, _d) in idx:
        m = g.split(";")[0].strip()
        bygene[genes.get(m) or genes.get(m.split("-")[0])].add(g)
    return idx, med, dict(bygene)


def common_medians(idx, med):
    """-> run medians recomputed over the peptides every kept run of the digest quantified."""
    runs = collections.defaultdict(set)
    for s, d in med:
        runs[d].add(s)
    vals = collections.defaultdict(list)
    shared = collections.Counter()
    for (g, d), peps in idx.items():
        need = runs.get(d)
        if not need:
            continue
        for p, q in peps.items():
            if need <= {s for s, v in q.items() if v > 0}:
                shared[d] += 1
                for s in need:
                    vals[(s, d)].append(math.log2(q[s]))
    out = {k: statistics.median(vals[k]) for k in med if vals.get(k)}
    for d, n in sorted(shared.items()):
        print(f"  {d}: {n:,} peptides quantified in all {len(runs[d])} runs")
    return out


def diagnostic(gene, groups, seqs, genes, idx):
    """Split a gene's quantified peptides into canonical and per-isoform sets, by FASTA sequence."""
    canon = [a for a in seqs
             if "-" not in a and (genes.get(a) == gene
                                  or genes.get(a.split("-")[0]) == gene)]
    if len(canon) != 1:
        return None, {}
    ref = seqs[canon[0]]
    family = [a for a in seqs if "-" in a and a.split("-")[0] == canon[0]]
    diag = collections.defaultdict(set)
    base = set()
    for g in groups:
        for d in ORDER:
            for p in idx.get((g, d), {}):
                if p in ref:
                    base.add((g, p))
                    continue
                who = tuple(sorted(a for a in family if p in seqs[a]))
                if who:
                    diag[who].add((g, p))
    return (canon[0], base), dict(diag)


def sample_levels(members, digest, idx, med):
    """-> {patient: run-normalised log2 of the members' summed quantity}."""
    tot = collections.defaultdict(float)
    for g, p in members:
        d = idx.get((g, digest), {}).get(p)
        if d:
            for s, v in d.items():
                tot[s] += v
    return {s: math.log2(v) - med[(s, digest)] for s, v in tot.items()
            if (s, digest) in med and v > 0}


def column(members, idx, med, cond, case, control):
    """One proteoform's per-digest change, peptides seen in both groups, and patient-level Welch p."""
    per, peps = {}, set()
    pooled = collections.defaultdict(list)
    for d in ORDER:
        vals = sample_levels(members, d, idx, med)
        peps |= {p for g, p in members
                 if {cond.get(s) for s, v in idx.get((g, d), {}).get(p, {}).items()
                     if (s, d) in med and v > 0} >= {case, control}}
        a = [v for s, v in vals.items() if cond.get(s) == case]
        b = [v for s, v in vals.items() if cond.get(s) == control]
        if a and b:
            per[d] = statistics.mean(a) - statistics.mean(b)
            for s, v in vals.items():
                pooled[s].append(v)
    if not per:
        return None
    import numpy as np
    from scipy import stats
    avg = {s: statistics.mean(v) for s, v in pooled.items()}
    a = [v for s, v in avg.items() if cond.get(s) == case]
    b = [v for s, v in avg.items() if cond.get(s) == control]
    t = stats.ttest_ind(a, b, equal_var=False)
    return {"per_digest": per, "mean": statistics.mean(per.values()),
            "sd": statistics.stdev(list(per.values())) if len(per) > 1 else 0.0,
            "npep": len(peps),
            "p": float(t.pvalue) if np.isfinite(t.pvalue) else float("nan")}


def precursor_values(members, idx, med):
    """-> [(peptide, patient, digest, log2 value centred on the peptide's mean)]."""
    out = []
    for d in ORDER:
        for g, pep in sorted(members):
            q = idx.get((g, d), {}).get(pep)
            if not q:
                continue
            v = {s: math.log2(x) - med[(s, d)] for s, x in q.items()
                 if (s, d) in med and x > 0}
            if not v:
                continue
            mid = statistics.mean(v.values())
            out += [(pep, s, d, x - mid) for s, x in sorted(v.items())]
    return out


def dd_fit(can_members, iso_members, idx, med, cond, case="LBD", control="Control"):
    """(isoform case − control) minus (canonical case − control), Welch–Satterthwaite. -> dict, or None."""
    from scipy import stats
    rows, parts, mean = {}, [], {}
    for form, members in ((1.0, iso_members), (0.0, can_members)):
        vals = precursor_values(members, idx, med)
        a = [v for _p, s, _d, v in vals if cond.get(s) == case]
        b = [v for _p, s, _d, v in vals if cond.get(s) == control]
        if len(a) < 2 or len(b) < 2:
            return None
        mean[(form, 1.0)], mean[(form, 0.0)] = statistics.mean(a), statistics.mean(b)
        parts += [(statistics.variance(a), len(a)), (statistics.variance(b), len(b))]
        for p, s, d, v in vals:
            if cond.get(s) in (case, control):
                rows.setdefault((p, d, form), []).append((s, v, 1.0 if cond.get(s) == case else 0.0))
    se2 = sum(v / k for v, k in parts)
    if se2 <= 0:
        return None
    df = se2 ** 2 / sum((v / k) ** 2 / (k - 1) for v, k in parts)
    delta = (mean[(1.0, 1.0)] - mean[(1.0, 0.0)]) - (mean[(0.0, 1.0)] - mean[(0.0, 0.0)])
    return {"rows": rows, "levels": mean, "interaction": delta, "se": math.sqrt(se2),
            "p": float(2 * stats.t.sf(abs(delta) / math.sqrt(se2), df))}


def screen(idx, med, bygene, seqs, genes, cond, case="LBD", control="Control"):
    """Test every isoform's own peptides against its canonical. -> [row dicts] with BH q."""
    rows = []
    for gene, groups in sorted(bygene.items()):
        if gene is None:
            continue
        base, diag = diagnostic(gene, groups, seqs, genes, idx)
        if base is None or not diag:
            continue
        if column(base[1], idx, med, cond, case, control) is None:
            continue
        for who, mem in sorted(diag.items()):
            st = column(mem, idx, med, cond, case, control)
            if st is None or not st["npep"]:
                continue
            f = dd_fit(base[1], mem, idx, med, cond, case, control)
            st["mean"], st["p"] = ((f["interaction"], f["p"]) if f
                                   else (float("nan"), float("nan")))
            rows.append({"gene": gene, "isoform": ";".join(who),
                         "log2fc": st["mean"], "sd": st["sd"], "p": st["p"]})
    ok = [r for r in rows if r["p"] == r["p"]]
    if ok:
        import numpy as np
        from scipy import stats
        q = stats.false_discovery_control(
            np.array([r["p"] for r in ok]), method="bh")
        for r, qq in zip(ok, q):
            r["q"] = float(qq)
    for r in rows:
        r.setdefault("q", float("nan"))
    rows.sort(key=lambda r: -(abs(r["log2fc"]) - 1.5 * r["sd"]))
    return rows


# ---------------------------------------------------------------------- draw ---
def diamond(x, y, r):
    """-> an SVG path for a diamond centred on (x, y)."""
    return f"M {x:.1f} {y - r:.1f} L {x + r:.1f} {y:.1f} L {x:.1f} {y + r:.1f} L {x - r:.1f} {y:.1f} Z"


def _quartiles(v):
    """Tukey box: (q1, median, q3, whisker_lo, whisker_hi), linear quantiles."""
    import numpy as np
    q1, q2, q3 = (float(x) for x in np.percentile(v, [25, 50, 75]))
    iqr = q3 - q1
    lo = min(x for x in v if x >= q1 - 1.5 * iqr)
    hi = max(x for x in v if x <= q3 + 1.5 * iqr)
    return q1, q2, q3, lo, hi


def _groups(fit):
    """-> ({form: [(peptide, digest, patient, case minus control mean)]}, {form: modelled level})."""
    lv = fit["levels"]
    level = {f: lv[(f, 1.0)] - lv[(f, 0.0)] for f in (0.0, 1.0)}
    groups = {f: [] for f in (0.0, 1.0)}
    for (p, d, f), v in fit["rows"].items():
        ctl = [y for _s, y, c in v if not c]
        if ctl:
            m = statistics.mean(ctl)
            groups[f] += [(p, d, s, y - m) for s, y, c in v if c]
    return groups, level


def model_plots(plots, out, width=1215.0, case="LBD", control="Control"):
    """Beeswarm per form of each case value minus its peptide's control mean, one plot per isoform, stacked."""
    import numpy as np
    from lib_svg import beeswarm
    W = width
    c0 = Canvas(W, 10, FONT, font_scale=pt_scale(W, W), out_w=W)
    u = c0.fs * TEXT_BOOST                     # units per pt
    ml, mr, ph = 31 * u, 10.0, 190.0
    head, foot, gap = 19 * u, 26 * u, 4 * u
    pw = W - ml - mr
    half = pw * 0.14
    xs = {0.0: 0.24, 1.0: 0.62}
    H = len(plots) * (head + ph + foot) + (len(plots) - 1) * gap
    c = Canvas(W, H, FONT, font_scale=c0.fs, out_w=W)
    col = assign(ORDER)
    laid = []
    for fit, gene, iso, canon in plots:
        groups, level = _groups(fit)
        vals = sorted(y for g in groups.values() for *_r, y in g)
        lo = min(vals[int(0.01 * len(vals))], -0.5, *level.values()) - 0.2
        hi = max(vals[int(0.99 * len(vals)) - 1], 0.5, *level.values()) + 0.2
        laid.append((fit, gene, iso, canon, groups, level, lo, hi))

    def swarms(rad, groups, Y, lo, hi):
        return {f: beeswarm([Y(y) for *_r, y in g if lo <= y <= hi], rad)
                for f, g in groups.items()}

    def fits(rad):
        for _f, _g, _i, _c, groups, _l, lo, hi in laid:
            Y = lambda v, lo=lo, hi=hi: ph * (hi - v) / (hi - lo)
            if any(max(map(abs, o), default=0) > half
                   for o in swarms(rad, groups, Y, lo, hi).values()):
                return False
        return True
    # one point size across the plots: the largest every group fits at
    rad = next((r for r in (3.8, 3.4, 3.0, 2.7, 2.4, 2.1, 1.8, 1.6, 1.4, 1.2) if fits(r)), 1.2)
    print(f"  swarm radius {rad}")
    for n, (fit, gene, iso, canon, groups, level, lo, hi) in enumerate(laid):
        top = n * (head + ph + foot + gap)
        mt = top + head

        def Y(v, lo=lo, hi=hi, mt=mt):
            return mt + ph * (hi - v) / (hi - lo)
        b2 = fit["interaction"]
        accs = iso.split(";")
        iso_lab = accs[0] + "".join("/-" + a.rsplit("-", 1)[1] for a in accs[1:])
        q = fit.get("q", float("nan"))
        qtext = ((f"q = {q:.2g}" if q >= 0.001 else f"q = {q:.1e}") if q == q
                 else f"p = {fit['p']:.1e}")
        base = mt - 7 * u
        c.text(ml, base, gene, 12, INK, "start", "600")
        x = ml + text_width(gene, 12 * u, True) + 10 * u
        dtext = f"Δ = {b2:+.2f} ± {fit['se']:.2f}"
        c.text(x, base, dtext, 10, UNION, "start", "600")
        c.text(x + text_width(dtext, 10 * u, True) + 10 * u, base, qtext, 10, INK_SECONDARY, "start")
        step = 0.5 if hi - lo < 3 else 1 if hi - lo < 8 else 2
        t = math.floor(lo)
        while t <= hi:
            if t >= lo:
                c.line(ml, Y(t), ml + pw, Y(t), stroke=AXIS if t == 0 else GRID,
                       sw=1.2 if t == 0 else 0.8)
                c.text(ml - 6, Y(t) + 3.5 * u, (f"{t:+g}" if t else "0"), 10, INK_MUTED, "end")
            t += step
        c.rect(ml, mt, pw, ph, stroke=AXIS, sw=0.8)
        offs = swarms(rad, groups, Y, lo, hi)
        for f, frac in xs.items():
            cx = ml + pw * frac
            inside = [(d, y) for _p, d, _s, y in groups[f] if lo <= y <= hi]
            v = np.array([y for _d, y in inside])
            for o, (d, y) in zip(offs[f], inside):
                c.add(f'<circle cx="{cx + o:.1f}" cy="{Y(y):.1f}" r="{rad}" fill="{col[d]}" '
                      f'fill-opacity="0.8"/>')
            if len(v) >= 3:
                # an unfilled box over the points: quartiles, median, whiskers
                bw = half * 0.16
                q1, q2, q3, wl, wh = _quartiles(list(v))
                c.add(f'<rect x="{cx - bw:.1f}" y="{Y(q3):.1f}" width="{2 * bw:.1f}" '
                      f'height="{Y(q1) - Y(q3):.1f}" fill="none" stroke="{INK}" stroke-width="1.2"/>')
                for w0, w1 in ((wh, q3), (q1, wl)):
                    c.line(cx, Y(w0), cx, Y(w1), stroke=INK, sw=1.0)
                c.line(cx - bw, Y(q2), cx + bw, Y(q2), stroke=INK, sw=2.0)
            c.add(f'<path d="{diamond(cx, Y(level[f]), 5.5)}" fill="{INK}" '
                  f'stroke="white" stroke-width="1"/>')
            c.text(cx, mt + ph + 11 * u, "isoform" if f else "canonical", 10, INK, "middle")
            c.text(cx, mt + ph + 22.5 * u, iso_lab if f else canon, 10, INK, "middle")
        y0, y1 = Y(level[0.0]), Y(level[1.0])
        edge = max(max(map(abs, offs[1.0]), default=0) + rad, half * 0.16)
        bx = ml + pw * xs[1.0] + edge + 12
        tip = 4 if y1 > y0 else -4
        c.line(bx, y0, bx, y1 - tip, stroke=UNION, sw=1.8)
        c.add(f'<path d="M {bx - 4:.1f} {y1 - tip:.1f} L {bx + 4:.1f} {y1 - tip:.1f} '
              f'L {bx:.1f} {y1:.1f} Z" fill="{UNION}"/>')
        c.text(bx + 8, (y0 + y1) / 2 + 3.5 * u, f"Δ {b2:+.2f}", 11, UNION, "start", "600")
        print(f"  {gene} {iso}: Δ {b2:+.2f}, p {fit['p']:.1e}, q {q:.1e}")
    c.text(2 + 8 * u, H / 2, f"log2 {case} / {control}, per peptide", 11, INK_SECONDARY,
           "middle", rot=-90)
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({W:.0f}x{H:.0f})")


def volcano(rows, out, letter="", q_cut=0.05):
    """Every screened isoform: Δ against -log10 p, BH survivors marked and labelled."""
    ok = [r for r in rows if r["p"] == r["p"] and r["p"] > 0]
    hit = [r for r in ok if r["q"] <= q_cut]
    red = DIVERGING_HIGH[2]
    W = FIG_PANEL_W
    ml, mr, mt, ph = 84.0, 30.0, 84.0, 250.0
    H = mt + ph + 70
    c = Canvas(W, H, FONT, font_scale=pt_scale(W, 1215.0), out_w=1215.0)
    u = c.fs * TEXT_BOOST                  # canvas units per pt
    pw = W - ml - mr
    xm = max(abs(r["log2fc"]) for r in ok) * 1.08
    ym = max(-math.log10(r["p"]) for r in ok) * 1.08

    def X(v):
        return ml + (v + xm) / (2 * xm) * pw

    def Y(v):
        return mt + ph - v / ym * ph

    c.text(20, 30, letter, 13, INK, "start", "600")
    c.text(ml, 32, "Isoform-specific differential abundance", 12, INK, "start", "600")
    c.text(ml, 32 + 12 * u, f"n = {len(ok):,} isoforms · {len(hit)} with q ≤ {q_cut:g}",
           10, INK_MUTED, "start")
    step = 2.0 if ym > 8 else 1.0
    t = 0.0
    while t <= ym:
        c.line(ml, Y(t), ml + pw, Y(t), stroke=GRID, sw=1.0, so=0.9)
        c.text(ml - 8, Y(t) + 3.5 * u, f"{t:g}", 10, INK_MUTED, "end")
        t += step
    xs = 0.5 if xm < 2.5 else 1.0
    t = -math.floor(xm / xs) * xs
    while t <= xm:
        c.text(X(t), mt + ph + 4 + 8 * u, f"{t:+g}" if t else "0", 10, INK_MUTED,
               "middle")
        t += xs
    c.line(X(0), mt, X(0), mt + ph, stroke=AXIS, sw=1.0)
    c.rect(ml, mt, pw, ph, fill="none", stroke=AXIS, sw=0.8, rx=0)
    c.text(ml + pw / 2, mt + ph + 6 + 20 * u, "Δ log2 LBD / Control (isoform − canonical)", 11,
           INK_SECONDARY, "middle")
    c.text(10 + 8 * u, mt + ph / 2, "−log10 p", 11, INK_SECONDARY, "middle", rot=-90)
    for r in sorted(ok, key=lambda r: r["q"] <= q_cut):
        on = r["q"] <= q_cut
        c.add(f'<circle cx="{X(r["log2fc"]):.1f}" '
              f'cy="{Y(-math.log10(r["p"])):.1f}" r="{3.4 if on else 2.4}" '
              f'fill="{red if on else INK_MUTED}" '
              f'fill-opacity="{0.9 if on else 0.35}"/>')
    fsz = 10.5 * u                     # the size the label is drawn at
    pts = [(X(r["log2fc"]), Y(-math.log10(r["p"]))) for r in hit]
    boxes = []
    names = [r["gene"] for r in hit]
    # a gene labelled twice gets each isoform's suffix
    lab = {id(r): (f'{r["gene"]} -{r["isoform"].split(";")[0].split("-")[-1]}'
                   if names.count(r["gene"]) > 1 else r["gene"]) for r in hit}

    def free(x0, y0, x1, y1, own):
        if x0 < ml + 2 or x1 > ml + pw - 2:
            return False
        if any(x0 < bx1 and bx0 < x1 and y0 < by1 and by0 < y1
               for bx0, by0, bx1, by1 in boxes):
            return False
        return not any(x0 - 4 < px < x1 + 4 and y0 - 4 < py < y1 + 4
                       for i, (px, py) in enumerate(pts) if i != own)

    # strongest first, each label in the first free slot: outward, then inward, nearest row first
    for i, r in sorted(enumerate(hit), key=lambda t: t[1]["p"]):
        px, py = pts[i]
        tw = text_width(lab[id(r)], fsz, True)
        out_ = 1 if r["log2fc"] >= 0 else -1
        for dy in (0, -1.1 * fsz, 1.1 * fsz, -2.2 * fsz, 2.2 * fsz):
            done = False
            for side in (out_, -out_):
                x0 = px + 8 if side > 0 else px - 8 - tw
                # the box spans cap height to descender around the baseline
                y0 = py + dy - fsz * 0.45
                if free(x0, y0, x0 + tw, y0 + fsz, i):
                    boxes.append((x0, y0, x0 + tw, y0 + fsz))
                    c.text(x0, py + dy + fsz * 0.3, lab[id(r)], 10.5, INK,
                           "start", "600")
                    done = True
                    break
            if done:
                break
        else:
            print(f"  no free slot for {r['gene']}'s label; not drawn")
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({W:.0f}x{H:.0f}, {len(ok)} isoforms, "
          f"{len(hit)} at q ≤ {q_cut:g})")


# ---------------------------------------------------------------------- main ---
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reports", nargs="+")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--volcano", metavar="SVG", help="volcano of every screened isoform")
    ap.add_argument("--model-plots", dest="model_plots", metavar="GENE:ISO,...",
                    help="model plots, one per isoform, stacked into --out")
    ap.add_argument("--out", help="SVG for --model-plots")
    ap.add_argument("--width", type=float, default=1215.0, help="--model-plots width")
    args = ap.parse_args(argv)
    if not args.volcano and not args.model_plots:
        ap.error("--volcano or --model-plots")
    if args.model_plots and not args.out:
        ap.error("--model-plots needs --out")

    paths = sorted(p for pat in args.reports for p in glob.glob(pat))
    if not paths:
        sys.exit("no reports matched")
    seqs = read_fasta(args.fasta)
    genes = gene_map(args.fasta)
    cond = read_metadata(args.metadata)
    idx, med, bygene = read_reports(paths, seqs, genes)
    med = common_medians(idx, med)
    rows = screen(idx, med, bygene, seqs, genes, cond)
    if args.volcano:
        volcano(rows, args.volcano)
    if args.model_plots:
        qs = {(r["gene"], r["isoform"]): r["q"] for r in rows}
        plots = []
        for item in args.model_plots.split(","):
            gene, iso = item.split(":")
            base, diag = diagnostic(gene, bygene[gene], seqs, genes, idx)
            fit = dd_fit(base[1], diag[tuple(iso.split(";"))], idx, med, cond)
            fit["q"] = qs.get((gene, ";".join(sorted(iso.split(";")))), float("nan"))
            plots.append((fit, gene, iso, base[0]))
        model_plots(plots, args.out, args.width)
    return 0


if __name__ == "__main__":
    sys.exit(main())
