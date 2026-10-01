#!/usr/bin/env python3
"""Test whether an isoform's own (diagnostic) peptides change between conditions.

    python3 fig2a_supp1_diagnostic_peptides.py 'data/search/*-60min-Phospho.parquet' --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta \
        --metadata data/metadata.xlsx --genes MAP2,SPTB,SLC25A25,RTN1 \
        --letter c --out figures/fig2a_dots.svg --tsv derived/da_iso/diagnostic_peptides.tsv

    # every isoform with diagnostic evidence, ranked
    python3 fig2a_supp1_diagnostic_peptides.py 'data/search/*-60min-Phospho.parquet' --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta \
        --metadata data/metadata.xlsx --screen
"""

import argparse
import collections
import csv
import glob
import hashlib
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_svg import Canvas                                  # noqa: E402
from lib_palette import (AXIS, DIVERGING_HIGH, FONT, GRID, INK, INK_MUTED,
                     INK_SECONDARY, TEXT_BOOST, UNION, assign,
                     display)                                      # noqa: E402
from lib_report import ORDER, read_metadata, sample_of                 # noqa: E402
from lib_fasta import read_fasta                            # noqa: E402


# width every Figure 2 panel is authored at
FIG_PANEL_W = 1010.0


# --------------------------------------------------------------------- data ---
def read_reports(paths, seqs, genes, precursor_q=0.01, min_run_peptides=1000,
                 want=None, key="peptide", quantity="quantity"):
    """-> (idx, med, bygene). One pass, all digests."""
    import lib_report as rp
    idx = collections.defaultdict(lambda: collections.defaultdict(dict))
    allq = collections.defaultdict(list)
    prec = key == "precursor"
    FIELDS = ["run", "protein_groups", "peptide" if prec else key,
              "precursor_q", quantity]
    extra = ["modified", "charge"] if prec else []
    keep = {}
    for r in rp.open_reports(paths, order=ORDER):
        path, dig = r.path, r.protease
        # DIA-NN has no used-for-quantity flag, so every precursor counts
        use = ["use_for_pep"] if r.has("use_for_pep") else []
        fields = FIELDS + use + extra
        for b in r.batches(fields):
            run, grp, pep, q, fg = (b[f] for f in fields[:5])
            used = b["use_for_pep"] if use else [True] * b["_n"]
            if prec:
                pep = [f"{a}|{m}|{z}" for a, m, z in
                       zip(pep, b["modified"], b["charge"])]
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
                                    and all(a in seqs for a in mem)
                                    and (want is None or gn <= want))
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
    """-> run medians recomputed over the peptides every kept run of the digest quantified.

    A shallow run misses mostly weak peptides, which inflates a median over its own detections;
    a median over the shared set is independent of depth.
    """
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


def diagnostic(gene, groups, seqs, genes, idx, fasta_isoforms=False):
    """Split a gene's quantified peptides into canonical and isoform-diagnostic.

    `fasta_isoforms`: a peptide's isoforms are every isoform of the gene whose sequence
    contains it, not only those in the protein group its digest reported it under.
    """
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
        mem = family if fasta_isoforms else [a.strip() for a in g.split(";")
                                             if a.strip() in seqs]
        for d in ORDER:
            for p in idx.get((g, d), {}):
                # a precursor key leads with its stripped sequence
                sq = p.split("|")[0]
                if sq in ref:
                    base.add((g, p))
                    continue
                who = tuple(sorted(a for a in mem if sq in seqs[a]))
                if who:
                    diag[who].add((g, p))
    return (canon[0], base), dict(diag)


def sample_levels(members, digest, idx, med):
    """-> ({sample: run-centred log2 quantity}, peptides contributing)."""
    tot = collections.defaultdict(float)
    n = 0
    for g, p in members:
        d = idx.get((g, digest), {}).get(p)
        if d:
            n += 1
            for s, v in d.items():
                tot[s] += v
    return ({s: math.log2(v) - med[(s, digest)] for s, v in tot.items()
             if (s, digest) in med and v > 0}, n)


def point_patients(members, digest, idx, med, cond, case, control):
    """-> [(peptide, log2 case-minus-control, case patients, control patients)]."""
    out = []
    for g, p in sorted(members):
        d = idx.get((g, digest), {}).get(p)
        if not d:
            continue
        vals = {s: math.log2(v) - med[(s, digest)] for s, v in d.items()
                if (s, digest) in med and v > 0}
        a = {s: v for s, v in vals.items() if cond.get(s) == case}
        b = {s: v for s, v in vals.items() if cond.get(s) == control}
        if not a or not b:
            continue
        out.append((p, statistics.mean(a.values()) - statistics.mean(b.values()),
                    set(a), set(b)))
    return out


def peptide_points(members, digest, idx, med, cond, case, control):
    """-> [(peptide, log2 case-minus-control, n_case, n_ctrl)], one per peptide."""
    return [(p, v, len(a), len(b)) for p, v, a, b in
            point_patients(members, digest, idx, med, cond, case, control)]


def patient_ratios(members, idx, med, cond, control, centre=True):
    """One value per patient, built from the same per-peptide ratios as the dots."""
    per = collections.defaultdict(list)
    for d in ORDER:
        for g, pep in members:
            q = idx.get((g, d), {}).get(pep)
            if not q:
                continue
            v = {s: math.log2(x) - med[(s, d)] for s, x in q.items()
                 if (s, d) in med and x > 0}
            ref = [x for s, x in v.items() if cond.get(s) == control]
            if not ref:
                continue
            mid = statistics.mean(ref) if centre else 0.0
            for s, x in v.items():
                per[s].append(x - mid)
    return {s: statistics.mean(v) for s, v in per.items()}


def vs_canonical(can, iso):
    """Welch between the two columns' drawn dots. -> (p, effect, df, label)."""
    import numpy as np
    from scipy import stats
    a = [v for _p, _d, v in iso["points"]]
    b = [v for _p, _d, v in can["points"]]
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan"), float("nan"), f"{len(a)} v {len(b)}"
    t = stats.ttest_ind(a, b, equal_var=False)
    p = float(t.pvalue) if np.isfinite(t.pvalue) else float("nan")
    return (p, statistics.mean(a) - statistics.mean(b), float(t.df),
            f"{len(a)} v {len(b)}")


def apply_vs_canonical(can, iso):
    """Rewrite `iso`'s p and effect as the between-column contrast, in place."""
    p, eff, _df, lab = vs_canonical(can, iso)
    iso["p"], iso["mean"] = p, eff
    iso["unit"], iso["n_unit"] = "peptides", lab
    iso["delta"] = True
    # mark the canonical column too, so both bars are their own dots' mean
    can["unit"] = "peptides"
    return iso


def column(members, idx, med, cond, case, control, min_group,
           p_unit="patients"):
    """Everything the panel needs about one proteoform column."""
    per, points = {}, []
    tested_pep = set()
    box_case, box_ctrl = set(), set()
    pooled = collections.defaultdict(list)
    for d in ORDER:
        vals, npep = sample_levels(members, d, idx, med)
        pp = point_patients(members, d, idx, med, cond, case, control)
        for _p, _v, pa, pb in pp:
            box_case |= pa
            box_ctrl |= pb
        pts = [(p, v, len(pa), len(pb)) for p, v, pa, pb in pp]
        points += [(p, d, v) for p, v, _na, _nb in pts]
        # `npep` keeps the patient threshold; only the drawing is ungated
        tested_pep |= {p.split("|")[0] for p, _v, na, nb in pts
                       if na >= min_group and nb >= min_group}
        a = [v for s, v in vals.items() if cond.get(s) == case]
        b = [v for s, v in vals.items() if cond.get(s) == control]
        if len(a) >= min_group and len(b) >= min_group:
            per[d] = statistics.mean(a) - statistics.mean(b)
            for s, v in vals.items():
                pooled[s].append(v)
    if not per:
        return None
    import numpy as np
    from scipy import stats

    if p_unit == "peptides":
        vs = [v for _p, _d, v in points]
        p, eff = float("nan"), statistics.mean(per.values())
        if len(vs) > 1:
            t = stats.ttest_1samp(vs, 0.0)
            p = float(t.pvalue) if np.isfinite(t.pvalue) else float("nan")
            eff = statistics.mean(vs)
        return {"per_digest": per, "mean": eff,
                "sd": statistics.stdev(vs) if len(vs) > 1 else 0.0,
                "points": points, "npep": len(tested_pep),
                "box": (len(box_case), len(box_ctrl)),
                "n_case": len(vs), "n_ctrl": 0, "p": p,
                "unit": "peptides", "n_unit": len(vs)}

    if p_unit == "patient-ratio":
        avg = patient_ratios(members, idx, med, cond, control)
    else:
        avg = {s: statistics.mean(v) for s, v in pooled.items()}
    a = [v for s, v in avg.items() if cond.get(s) == case]
    b = [v for s, v in avg.items() if cond.get(s) == control]
    p = float("nan")
    eff = statistics.mean(per.values())
    if len(a) >= min_group and len(b) >= min_group:
        t = stats.ttest_ind(a, b, equal_var=False)
        p = float(t.pvalue) if np.isfinite(t.pvalue) else float("nan")
        if p_unit == "patient-ratio":
            # effect from the same rollup as the p
            eff = statistics.mean(a) - statistics.mean(b)
    return {"per_digest": per,
            "mean": eff,
            "sd": statistics.stdev(list(per.values())) if len(per) > 1 else 0.0,
            "points": points,
            "npep": len(tested_pep),
            "box": (len(box_case), len(box_ctrl)),
            "n_case": len(a), "n_ctrl": len(b), "p": p,
            "unit": "patients", "n_unit": f"{len(a)} v {len(b)}"}


def case_points(members, idx, med, cond, case, control, centre=True):
    """-> [log2 vs the peptide's control mean], one per case patient per peptide."""
    return [v for _p, s, _d, v in precursor_values(members, idx, med, cond, control, centre)
            if cond.get(s) == case]


def apply_case_points(can_members, iso_members, iso, idx, med, cond, case, control,
                      centre=True):
    """Rewrite `iso`'s p and effect as Welch between the two forms' case points."""
    from scipy import stats
    a = case_points(iso_members, idx, med, cond, case, control, centre)
    b = case_points(can_members, idx, med, cond, case, control, centre)
    iso["p"], iso["mean"] = float("nan"), float("nan")
    if len(a) >= 2 and len(b) >= 2:
        t = stats.ttest_ind(a, b, equal_var=False)
        iso["p"] = float(t.pvalue) if t.pvalue == t.pvalue else float("nan")
        iso["mean"] = statistics.mean(a) - statistics.mean(b)
    iso["unit"], iso["n_unit"] = "case points", f"{len(a)} v {len(b)}"
    return iso


def peptide_model_fit(can_members, iso_members, idx, med, cond, case, control):
    """OLS of log2 quantity on peptide + case + case x isoform, one row per peptide per patient.

    The peptide term is absorbed by demeaning within each peptide (per digest), so every
    peptide keeps its own baseline. -> dict with the rows, each peptide's fitted baseline,
    the case effect, the case x isoform coefficient, its standard error and p; None if unfit.
    """
    import numpy as np
    from scipy import stats
    rows = collections.defaultdict(list)
    for form, members in ((0.0, can_members), (1.0, iso_members)):
        for d in ORDER:
            for g, p in members:
                q = idx.get((g, d), {}).get(p)
                if not q:
                    continue
                for s, x in q.items():
                    if (s, d) in med and x > 0 and cond.get(s) in (case, control):
                        rows[(p, d, form)].append((s, math.log2(x) - med[(s, d)],
                                                   1.0 if cond.get(s) == case else 0.0))
    y, X = [], []
    for (_p, _d, form), v in rows.items():
        a = np.array([(yy, cc) for _s, yy, cc in v])
        a = a - a.mean(axis=0)
        y.extend(a[:, 0])
        X.extend(zip(a[:, 1], a[:, 1] * form))     # form is constant within a peptide
    n, k = len(y), len(rows)
    dof = n - k - 2
    y, X = np.array(y), np.array(X)
    if dof < 1 or np.linalg.matrix_rank(X) < 2:
        return None
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    cov = float(resid @ resid) / dof * np.linalg.inv(X.T @ X)
    se = math.sqrt(cov[1, 1])
    base = {key: statistics.mean(yy for _s, yy, _c in v)
            - beta[0] * statistics.mean(c for _s, _y, c in v)
            - beta[1] * form * statistics.mean(c for _s, _y, c in v)
            for key, v in rows.items() for form in (key[2],)}
    return {"rows": rows, "baseline": base, "case": float(beta[0]),
            "interaction": float(beta[1]), "se": se, "dof": dof,
            "p": float(2 * stats.t.sf(abs(beta[1] / se), dof)), "n": n, "k": k}


def peptide_model(can_members, iso_members, idx, med, cond, case, control):
    """-> (case x isoform coefficient, p, rows, peptides) from `peptide_model_fit`."""
    f = peptide_model_fit(can_members, iso_members, idx, med, cond, case, control)
    if f is None:
        return float("nan"), float("nan"), 0, 0
    return f["interaction"], f["p"], f["n"], f["k"]


def testable(st, can, min_pep, box_patients=0):
    """Whether an isoform column enters the screen."""
    if st is None:
        return False
    if box_patients:
        # each box's points, pooled, span >= box_patients per group
        return bool(can) and all(
            len({p.split("|")[0] for p, _d, _v in c["points"]}) >= min_pep
            and min(c["box"]) >= box_patients for c in (st, can))
    return st["npep"] >= min_pep


def screen(idx, med, bygene, seqs, genes, cond, case, control, min_group,
           min_pep, min_digests, p_unit="patients", box_patients=0, centre=True,
           fasta_isoforms=False):
    """Rank every isoform with diagnostic evidence. -> [row dicts], BH within."""
    rows = []
    for gene, groups in sorted(bygene.items()):
        if gene is None:
            continue
        base, diag = diagnostic(gene, groups, seqs, genes, idx, fasta_isoforms)
        if base is None or not diag:
            continue
        can = column(base[1], idx, med, cond, case, control, min_group,
                     p_unit)
        for who, mem in sorted(diag.items()):
            st = column(mem, idx, med, cond, case, control, min_group,
                        p_unit)
            if st and p_unit == "vs-canonical" and can:
                apply_vs_canonical(can, st)
            if st and p_unit == "peptide-model" and can:
                st["mean"], st["p"], _n, _k = peptide_model(base[1], mem, idx, med, cond,
                                                            case, control)
            if st and p_unit == "lbd-points" and can:
                apply_case_points(base[1], mem, st, idx, med, cond, case, control,
                                  centre)
            if not testable(st, can, min_pep, box_patients):
                continue
            if len(st["per_digest"]) < min_digests:
                continue
            rows.append({"gene": gene, "isoform": ";".join(who),
                         "canonical": base[0],
                         "digests": ",".join(sorted(st["per_digest"])),
                         "log2fc": st["mean"], "sd": st["sd"],
                         "n_pep": st["npep"], "p": st["p"],
                         "canonical_log2fc": can["mean"] if can else float("nan"),
                         "canonical_digests": len(can["per_digest"]) if can
                         else 0})
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
    # replication first, as in `proteoform_bands.py --screen`
    rows.sort(key=lambda r: -(abs(r["log2fc"]) - 1.5 * r["sd"]))
    return rows


# ---------------------------------------------------------------------- draw ---
def _jitter(pep, digest, width):
    """Deterministic horizontal offset."""
    h = hashlib.blake2b(f"{pep}{digest}".encode(), digest_size=4).digest()
    return (int.from_bytes(h, "big") / 2 ** 32 - 0.5) * width


def _wrap(text, chars):
    """Greedy wrap to a character budget. -> [line]."""
    out, line = [], ""
    for word in text.split():
        cand = f"{line} {word}".strip()
        if len(cand) > chars and line:
            out.append(line)
            line = word
        else:
            line = cand
    if line:
        out.append(line)
    return out


def _stat(st):
    """A column's effect and p, as the facet subtitle prints them."""
    pv = ("p < 1e-4" if st["p"] < 1e-4 else f"p = {st['p']:.1e}"
          if st["p"] < 1e-3 else f"p = {st['p']:.3f}") \
        if st["p"] == st["p"] else "p n/a"
    lead = ("Δ " if st.get("delta") else "") + f"{st['mean']:+.2f}"
    return f"{lead}  ·  {pv}"


def _label(who):
    a = who.split(";")
    return a[0] if len(a) == 1 else f"{a[0]} (+{len(a) - 1})"


def panel(picks, out, font, letter="c", width=1215.0, ts=1.7, key=True,
          ylim=None, subtitle=True, panel_w=None, panel_h=None, hue=False,
          case_label="LBD", control_label="Control", note="",
          title=None, head_text=None, label_fn=None):
    """Facet per gene; column per proteoform; a dot per peptide per digest."""
    # facets fill the panel width so composed panels share one type size
    ncol = [1 + len(g["diag"]) for g in picks]
    ph, gapx = 176.0, 26.0
    ml = 62.0
    W = panel_w or FIG_PANEL_W
    span = W - ml - 22 - gapx * (len(picks) - 1)
    fw = [span * n / sum(ncol) for n in ncol]
    fs = 1.45 * ts / 1.7
    # ~0.5 em per character, as in `test_units.py`'s overflow check
    budget = int((W - ml - 20) / (0.5 * 8.4 * fs))
    head = _wrap(head_text if head_text is not None else
                 (f"one point per peptide per protease · log2 {case_label} / "
                  f"{control_label} · each column against zero"), budget) \
        if subtitle else []
    foot = _wrap(note, budget) if note else []
    mt = 54.0 + len(head) * 13.0 + 28.0
    H = panel_h or mt + ph + (48 + len(foot) * 12.0 + 8 if foot else 30)
    c = Canvas(W, H, font, font_scale=fs, out_w=width)
    colour = assign(ORDER + ["All"])

    vals = [v for g in picks for col in [g["canon"]] + [d[1] for d in g["diag"]]
            for _p, _d, v in col["points"]]
    if ylim:
        lo, hi = ylim
    else:
        lo, hi = min(vals + [0.0]), max(vals + [0.0])
        pad = (hi - lo) * 0.10 or 0.5
        lo, hi = lo - pad, hi + pad
    # points outside a fixed scale are not drawn; they are counted on stdout
    n_off = sum(1 for v in vals if v < lo or v > hi)

    def Y(v):
        return mt + ph - (v - lo) / (hi - lo) * ph

    c.text(20, 30, letter, 13, INK, "start", "600")
    c.text(ml, 30, title or "Differential abundance of isoform-discriminating "
           "peptides", 11.5, INK, "start", "600")
    for i, line in enumerate(head):
        c.text(ml, 46 + i * 13, line, 8.4, INK_MUTED, "start")
    c.text(18, mt + ph / 2, f"log2 {case_label} / {control_label}", 9.4,
           INK_SECONDARY, "middle", rot=-90)
    if key:
        kx = ml + sum(fw) + gapx * (len(picks) - 1) - 150
        for i, d in enumerate(ORDER):
            c.add(f'<circle cx="{kx + i * 52:.1f}" cy="27" r="3.4" '
                  f'fill="{colour[d]}" fill-opacity="0.85"/>')
            c.text(kx + 7 + i * 52, 30, display(d), 7.8, INK_SECONDARY, "start")
    for i, line in enumerate(foot):
        c.text(ml, mt + ph + 60 + i * 12, line, 8.0, INK_MUTED,
               "start")

    # y grid, once, so the facets share a scale
    step = 1.0 if hi - lo > 2.5 else 0.5
    t = math.ceil(lo / step) * step
    while t <= hi:
        y = Y(t)
        c.line(ml, y, ml + sum(fw) + gapx * (len(picks) - 1), y,
               stroke=GRID if abs(t) > 1e-9 else AXIS, sw=1.0,
               so=0.9 if abs(t) > 1e-9 else None)
        c.text(ml - 6, y + 3, f"{t:g}", 8.0, INK_MUTED, "end")
        t += step

    x0 = ml
    for g, w in zip(picks, fw):
        c.rect(x0, mt, w, ph, fill="none", stroke=AXIS, sw=0.8, rx=0)
        c.text(x0 + w / 2, mt - 24, g["gene"], 10.4, INK, "middle", "600")
        best = g["diag"][0][1] if g["diag"] else None
        if g.get("sub"):
            # caller-supplied subtitle
            c.text(x0 + w / 2, mt - 7, g["sub"], 8.4, INK_SECONDARY, "middle")
        elif best is not None and not g.get("per_column"):
            c.text(x0 + w / 2, mt - 7, _stat(best), 8.4, INK_SECONDARY, "middle")
        lf = label_fn or _label
        cols = [(lf(g["canon_acc"]), g["canon"], True)]
        cols += [(lf(who), st, False) for who, st in g["diag"]]
        cw = w / len(cols)
        for j, (lab, st, is_canon) in enumerate(cols):
            cx = x0 + cw * (j + 0.5)
            for pep, dig, v in st["points"]:
                px = cx + _jitter(pep, dig, cw * 0.44)
                r_ = 2.6 if is_canon else 3.4
                if v < lo or v > hi:
                    continue
                col = colour[dig] if hue else UNION
                c.add(f'<circle cx="{px:.1f}" '
                      f'cy="{Y(v):.1f}" r="{r_:.1f}" '
                      f'fill="{col}" fill-opacity='
                      f'"{0.34 if is_canon else 0.62}" stroke="{col}" '
                      f'stroke-width="0.8" stroke-opacity="0.85"/>')
            vs = [v for _p, _d, v in st["points"]]
            mid = (statistics.mean(vs) if st.get("unit") == "peptides"
                   else statistics.median(vs)) if vs else st["mean"]
            # in `vs-canonical` mode both bars are means, so their gap is the quoted Δ
            c.line(cx - cw * 0.40, Y(mid), cx + cw * 0.40, Y(mid),
                   stroke=INK, sw=2.4)
            # `per_column`: each contrast printed over its own column
            if g.get("per_column") and not is_canon:
                # under the label, on two lines, so it clears the neighbour
                lead, pv = _stat(st).split("  ·  ")
                c.text(cx, mt + ph + 27, lead, 7.6, INK_SECONDARY, "middle")
                c.text(cx, mt + ph + 38, pv, 7.6, INK_SECONDARY, "middle")
            c.text(cx, mt + ph + 14, lab, 8.2, INK if not is_canon
                   else INK_SECONDARY, "middle", "600" if not is_canon else None)
        x0 += w + gapx

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({W:.0f}x{H:.0f}, {len(picks)} gene(s))")
    if n_off:
        print(f"  NOT DRAWN: {n_off} of {len(vals)} points fall outside the "
              f"fixed y-scale [{lo:+.1f}, {hi:+.1f}]")


def precursor_values(members, idx, med, cond, control, centre=True):
    """-> [(key, patient, digest, value)], one per precursor per patient."""
    out = []
    for d in ORDER:
        for g, pep in sorted(members):
            q = idx.get((g, d), {}).get(pep)
            if not q:
                continue
            v = {s: math.log2(x) - med[(s, d)] for s, x in q.items()
                 if (s, d) in med and x > 0}
            ref = [x for s, x in v.items() if cond.get(s) == control]
            if not ref:
                continue
            mid = statistics.mean(ref) if centre else 0.0
            out += [(pep, s, d, x - mid) for s, x in sorted(v.items())]
    return out


def run_values(members, idx, med, cond, control, centre=True):
    """-> [(digest, patient, digest, value)], one per run."""
    per = collections.defaultdict(list)
    for k, s, d, v in precursor_values(members, idx, med, cond, control,
                                       centre):
        per[(s, d)].append(v)
    return [(d, s, d, statistics.mean(v)) for (s, d), v in sorted(per.items())]


def box_columns(pick, idx, med, cond, case, control, centre=True,
                unit="patient"):
    """Points for a gene's columns, and the proteoform contrast."""
    from scipy import stats
    who, mem = pick["diag_members"][0]
    can = patient_ratios(pick["canon_members"], idx, med, cond, control,
                         centre)
    iso = patient_ratios(mem, idx, med, cond, control, centre)
    if unit == "run":
        cols = [(pick["canon_acc"], run_values(
                    pick["canon_members"], idx, med, cond, control, centre),
                 True),
                (who, run_values(mem, idx, med, cond, control, centre),
                 False)]
    elif unit != "patient":
        cols = [(pick["canon_acc"], precursor_values(
                    pick["canon_members"], idx, med, cond, control, centre),
                 True),
                (who, precursor_values(mem, idx, med, cond, control, centre),
                 False)]
    else:
        cols = [(pick["canon_acc"],
                 [(s, s, None, v) for s, v in sorted(can.items())], True),
                (who, [(s, s, None, v) for s, v in sorted(iso.items())], False)]
    diff = {s: iso[s] - can[s] for s in iso if s in can}
    a = [v for s, v in diff.items() if cond.get(s) == case]
    b = [v for s, v in diff.items() if cond.get(s) == control]
    if len(a) >= 2 and len(b) >= 2:
        t = stats.ttest_ind(a, b, equal_var=False)
        p = float(t.pvalue) if t.pvalue == t.pvalue else float("nan")
        return cols, (statistics.mean(a) - statistics.mean(b), p, len(a), len(b))
    return cols, (float("nan"), float("nan"), len(a), len(b))


def _quartiles(v):
    """Tukey box: (q1, median, q3, whisker_lo, whisker_hi), linear quantiles."""
    import numpy as np
    q1, q2, q3 = (float(x) for x in np.percentile(v, [25, 50, 75]))
    iqr = q3 - q1
    lo = min(x for x in v if x >= q1 - 1.5 * iqr)
    hi = max(x for x in v if x <= q3 + 1.5 * iqr)
    return q1, q2, q3, lo, hi


def box_panel(boxes, cond, out, font, letter="a", width=1215.0, ylim=None,
              panel_w=None, panel_h=None, case_label="LBD",
              control_label="Control", title=None, subtitle=True,
              centre=True, unit="patient", per_row=None, stat="q", box_n=None, show_stat=True):
    """Facet per gene; per proteoform column, a Control box beside an LBD box."""
    ph, gapx, ml = 176.0, 26.0, 62.0
    W = panel_w or FIG_PANEL_W
    per_row = per_row or len(boxes)
    rows = [boxes[i:i + per_row] for i in range(0, len(boxes), per_row)]
    # the gap between rows holds one row's labels and the next row's headers
    rowgap = 64.0 + (14.0 if box_n else 0.0)
    fs = 1.45
    # wrapped at the drawn size, TEXT_BOOST included
    budget = int((W - ml - 20) / (0.5 * 8.4 * fs * TEXT_BOOST))
    fc = unit == "fc"
    # patient-valued units: one point per patient, no protease key
    prec = unit not in ("patient", "protein")
    what = unit if unit in ("peptide", "precursor") else "peptide"
    how = (f"each {what} centred on its {control_label} mean" if centre
           else f"mean log2 quantity over the column's {what}s, uncentred")
    if unit == "run" and centre:
        how = f"peptides centred on {control_label}, then averaged"
    per = {"patient": "", "run": "run (patient × enzyme)"}.get(
        unit, f"{what} per patient")
    head = _wrap(f"one point per {per or 'patient'} · "
                 f"{how} · Δ = isoform minus canonical within patient, Welch "
                 f"across patients", budget) if subtitle else []
    if fc and subtitle:
        # the test is between the two columns' own points
        head = _wrap(f"one point per peptide per enzyme · log2 {case_label} / "
                     f"{control_label} · Δ = isoform minus canonical, Welch "
                     f"between the columns, BH across all tested isoforms",
                     budget)
    mt = 46.0 + len(head) * 13.0 + 28.0
    top0 = mt
    body = len(rows) * ph + (len(rows) - 1) * rowgap
    H = panel_h if panel_h and len(rows) == 1 else \
        mt + body + 44 + (14.0 if box_n else 0.0)
    c = Canvas(W, H, font, font_scale=fs, out_w=width)

    vals = [pt[3] for _g, cols, _c in boxes for _l, d, _k in cols
            for pt in d]
    if ylim:
        lo, hi = ylim
    else:
        lo, hi = min(vals + [0.0]), max(vals + [0.0])
        pad = (hi - lo) * 0.10 or 0.5
        lo, hi = lo - pad, hi + pad
    n_off = sum(1 for v in vals if v < lo or v > hi)

    def Y(v):
        return mt + ph - (v - lo) / (hi - lo) * ph

    c.text(20, 30, letter, 13, INK, "start", "600")
    if title is None:
        title = "Isoforms of PD-implicated proteins, by " + unit
    if title:        # "" (--no-title) draws none; the protease key stays put
        c.text(ml, 30, title, 11.5, INK, "start", "600")
    for i, line in enumerate(head):
        c.text(ml, 46 + i * 13, line, 8.4, INK_MUTED, "start")
    c.text(18, mt + body / 2, "log2 vs run median" if fc and not centre
           else f"log2 {case_label} / {control_label}" if fc
           else f"log2 vs {control_label} mean" if centre
           else "log2 quantity (run-median centred)", 9.4,
           INK_SECONDARY, "middle", rot=-90)

    colour = assign(ORDER + ["All"])
    keyed = prec or any(pt[2] for _g, cols, _c in boxes
                        for _l, d, _k in cols for pt in d)
    if keyed:
        # each entry advances by its own name's width at the drawn type size
        adv = [7 + 0.56 * 7.8 * c.fs * TEXT_BOOST * len(display(dg)) + 12
               for dg in ORDER]
        kx = W - 22 - sum(adv) + 12
        for dg, a in zip(ORDER, adv):
            c.add(f'<circle cx="{kx:.1f}" cy="27" r="3.4" '
                  f'fill="{colour[dg]}" fill-opacity="0.85"/>')
            c.text(kx + 7, 30, display(dg), 7.8, INK_SECONDARY, "start")
            kx += a

    for r, row in enumerate(rows):
        mt = top0 + r * (ph + rowgap)
        ncol = [len(cols) for _g, cols, _c in row]
        span = W - ml - 22 - gapx * (per_row - 1)
        # a column's width is the same in every facet
        unitw = span / (per_row * max(len(cols) for _g, cols, _c in boxes))
        fw = [unitw * n for n in ncol]
        step = 2.0 if hi - lo > 8 else 1.0 if hi - lo > 2.5 else 0.5
        t = math.ceil(lo / step) * step
        while t <= hi:
            y = Y(t)
            c.line(ml, y, ml + sum(fw) + gapx * (len(row) - 1), y,
                   stroke=GRID if abs(t) > 1e-9 else AXIS, sw=1.0,
                   so=0.9 if abs(t) > 1e-9 else None)
            c.text(ml - 6, y + 3, f"{t:g}", 8.0, INK_MUTED, "end")
            t += step
        # Control hollow, the case tinted
        style = {control_label: ("none", INK_SECONDARY),
                 case_label: (UNION, UNION)}
        x0 = ml
        for (gene, cols, (delta, p, _na, _nb)), w in zip(row, fw):
            c.rect(x0, mt, w, ph, fill="none", stroke=AXIS, sw=0.8, rx=0)
            c.text(x0 + w / 2, mt - 20, gene, 10.4, INK, "middle", "600")
            if fc and show_stat:
                # `p` holds the BH q, or the raw Welch p with --stat p
                qs = (f"{stat} = {p:.1e}" if p < 1e-3
                      else f"{stat} = {p:.3f}")
                c.text(x0 + w / 2, mt - 7, f"Δ {delta:+.2f} · {qs}", 8.0,
                       INK_SECONDARY, "middle")
            elif delta == delta and show_stat:
                # at --unit protein: LBD minus Control of the whole protein
                c.text(x0 + w / 2, mt - 7, _stat({"mean": delta, "p": p,
                                                  "delta": not unit.startswith("protein")}),
                       8.4, INK_SECONDARY, "middle")
            cw = w / len(cols)
            for j, (lab, d, is_canon) in enumerate(cols):
                cx = x0 + cw * (j + 0.5)
                bw = min(cw * 0.26, 34.0)
                if fc:
                    bw = min(cw * 0.46, 40.0)
                for k, grp in enumerate((None,) if fc
                                        else (control_label, case_label)):
                    bx = cx if fc else cx + (k - 0.5) * cw * 0.40
                    v = sorted(x for _k, s, _d, x in d
                               if fc or cond.get(s) == grp)
                    fill, stroke = (style[control_label] if is_canon
                                    else style[case_label]) if fc else style[grp]
                    pts = [(k_, s, dg, x) for k_, s, dg, x in d
                           if (fc or cond.get(s) == grp) and lo <= x <= hi]
                    # dense precursor points go under the box; the rest on top
                    dense = prec and unit not in ("run", "fc")
                    if dense:
                        for k_, s, dg, x in pts:
                            px = bx + _jitter(f"{k_}{s}", lab, bw * 0.95)
                            c.add(f'<circle cx="{px:.1f}" cy="{Y(x):.1f}" '
                                  f'r="1.5" fill="{colour[dg]}" '
                                  f'fill-opacity="0.45"/>')
                    if len(v) >= 3:
                        q1, q2, q3, wl, wh = _quartiles(v)
                        c.line(bx, Y(wh), bx, Y(q3), stroke=stroke, sw=1.0)
                        c.line(bx, Y(q1), bx, Y(wl), stroke=stroke, sw=1.0)
                        for wy in (wl, wh):
                            c.line(bx - bw * 0.25, Y(wy), bx + bw * 0.25, Y(wy),
                                   stroke=stroke, sw=1.0)
                        c.add(f'<rect x="{bx - bw / 2:.1f}" y="{Y(q3):.1f}" '
                              f'width="{bw:.1f}" height="{Y(q1) - Y(q3):.1f}" '
                              f'fill="{fill}" fill-opacity="{0.10 if dense else 0.18}" '
                              f'stroke="{stroke}" stroke-width="1.2"/>')
                        c.line(bx - bw / 2, Y(q2), bx + bw / 2, Y(q2),
                               stroke=INK, sw=2.0)
                    for k_, s, dg, x in ([] if dense else pts):
                        px = bx + _jitter(f"{k_}{s}", lab, bw * 0.7)
                        c.add(f'<circle cx="{px:.1f}" cy="{Y(x):.1f}" '
                              f'r="{2.2 if fc else 2.8}" '
                              f'fill="{colour[dg] if dg else stroke}" '
                              f'fill-opacity="{0.6 if fc else 0.75}"/>')
                    if not fc:
                        c.text(bx, mt + ph + 12, grp, 7.4, INK_MUTED, "middle")
                        if box_n:
                            n = len({s for _k, s, _d, _x in d if cond.get(s) == grp})
                            c.text(bx, mt + ph + 23, f"n = {n}", 7.0, INK_MUTED,
                                   "middle")
                # one column per facet at --unit protein, titled by its gene
                if not unit.startswith("protein"):
                    c.text(cx, mt + ph + (14 if fc else 37 if box_n else 26),
                           _label(lab), 8.2,
                           INK_SECONDARY if is_canon else INK, "middle",
                           None if is_canon else "600")
                if fc and box_n and (gene, lab) in box_n:
                    na, nb = box_n[(gene, lab)]
                    c.text(cx, mt + ph + (26 if fc else 38),
                           f"{na} {case_label} v {nb} {control_label}", 7.4,
                           INK_MUTED, "middle")
            x0 += w + gapx

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({W:.0f}x{H:.0f}, {len(boxes)} gene(s))")
    if n_off:
        print(f"  NOT DRAWN: {n_off} of {len(vals)} {what if prec else 'patient'} points fall "
              f"outside the fixed y-scale [{lo:+.1f}, {hi:+.1f}]")


def volcano(rows, out, font, q_cut=0.05, letter="", title=None,
            xlab="Δ log2 LBD / Control (isoform − canonical)", width=1215.0,
            panel_w=None, panel_h=None):
    """Every screened isoform: effect against -log10 p, BH survivors marked."""
    ok = [r for r in rows if r["p"] == r["p"] and r["p"] > 0]
    hit = [r for r in ok if r["q"] <= q_cut]
    red = DIVERGING_HIGH[2]
    W = panel_w or FIG_PANEL_W
    ml, mr, mt, ph = 62.0, 30.0, 74.0, 250.0
    H = panel_h or mt + ph + 46
    c = Canvas(W, H, font, font_scale=1.45, out_w=width)
    pw = W - ml - mr
    xm = max(abs(r["log2fc"]) for r in ok) * 1.08
    ym = max(-math.log10(r["p"]) for r in ok) * 1.08

    def X(v):
        return ml + (v + xm) / (2 * xm) * pw

    def Y(v):
        return mt + ph - v / ym * ph

    c.text(20, 30, letter, 13, INK, "start", "600")
    c.text(ml, 30, title or "Isoform-specific differential abundance",
           11.5, INK, "start", "600")
    c.text(ml, 46, f"n = {len(ok):,} isoforms · {len(hit)} with q ≤ {q_cut:g}",
           8.4, INK_MUTED, "start")
    step = 2.0 if ym > 8 else 1.0
    t = 0.0
    while t <= ym:
        c.line(ml, Y(t), ml + pw, Y(t), stroke=GRID, sw=1.0, so=0.9)
        c.text(ml - 6, Y(t) + 3, f"{t:g}", 8.0, INK_MUTED, "end")
        t += step
    xs = 0.5 if xm < 2.5 else 1.0
    t = -math.floor(xm / xs) * xs
    while t <= xm:
        c.text(X(t), mt + ph + 14, f"{t:+g}" if t else "0", 8.0, INK_MUTED,
               "middle")
        t += xs
    c.line(X(0), mt, X(0), mt + ph, stroke=AXIS, sw=1.0)
    c.rect(ml, mt, pw, ph, fill="none", stroke=AXIS, sw=0.8, rx=0)
    c.text(ml + pw / 2, mt + ph + 32, xlab, 9.4, INK_SECONDARY, "middle")
    c.text(18, mt + ph / 2, "−log10 p", 9.4, INK_SECONDARY, "middle", rot=-90)
    for r in sorted(ok, key=lambda r: r["q"] <= q_cut):
        on = r["q"] <= q_cut
        c.add(f'<circle cx="{X(r["log2fc"]):.1f}" '
              f'cy="{Y(-math.log10(r["p"])):.1f}" r="{3.4 if on else 2.4}" '
              f'fill="{red if on else INK_MUTED}" '
              f'fill-opacity="{0.9 if on else 0.35}"/>')
    # labels: first free slot of right, left, above, below
    fsz = 8.2 * c.fs * TEXT_BOOST      # the size the label is drawn at
    pts = [(X(r["log2fc"]), Y(-math.log10(r["p"]))) for r in hit]
    boxes = []
    names = [r["gene"] for r in hit]
    # a gene labelled twice gets each isoform's suffix
    lab = {id(r): (f'{r["gene"]} -{r["isoform"].split(";")[0].split("-")[-1]}'
                   if names.count(r["gene"]) > 1 else r["gene"]) for r in hit}

    def free(x0, y0, x1, y1, own):
        if any(x0 < bx1 and bx0 < x1 and y0 < by1 and by0 < y1
               for bx0, by0, bx1, by1 in boxes):
            return False
        return not any(x0 - 4 < px < x1 + 4 and y0 - 4 < py < y1 + 4
                       for i, (px, py) in enumerate(pts) if i != own)

    for i, r in sorted(enumerate(hit), key=lambda t: t[1]["p"]):
        px, py = pts[i]
        # bold caps run ~0.66 em
        tw = 0.72 * fsz * len(lab[id(r)])   # bold caps at TEXT_BOOST
        out_ = 1 if r["log2fc"] >= 0 else -1
        for dy in (0, -12, 12, -24, 24):
            done = False
            for side in (out_, -out_):
                x0 = px + 8 if side > 0 else px - 8 - tw
                # the box spans cap height to descender around the baseline
                y0 = py + dy - fsz * 0.45
                if free(x0, y0, x0 + tw, y0 + fsz, i):
                    boxes.append((x0, y0, x0 + tw, y0 + fsz))
                    c.text(x0, py + dy + fsz * 0.3, lab[id(r)], 8.2, INK,
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
    ap.add_argument("reports", nargs="*", default=["data/search/*-60min-Phospho.parquet"])
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--case", default="LBD")
    ap.add_argument("--control", default="Control")
    ap.add_argument("--genes", help="comma-separated symbols to draw")
    ap.add_argument("--screen", action="store_true",
                    help="rank every isoform with diagnostic evidence")
    ap.add_argument("--min-group", type=int, default=3)
    ap.add_argument("--box-n", dest="box_n", action="store_true",
                    help="label each box with the patients behind its points")
    ap.add_argument("--common-median", dest="common_median", action="store_true",
                    help="normalise each run by the median of peptides every run of its digest quantified")
    ap.add_argument("--fasta-isoforms", dest="fasta_isoforms", action="store_true",
                    help="assign each peptide to every isoform whose sequence contains it")
    ap.add_argument("--no-run-median", dest="no_run_median", action="store_true",
                    help="skip the run-median subtraction (runs below --min-run-peptides still dropped)")
    ap.add_argument("--per-patient", dest="per_patient", action="store_true",
                    help="--unit fc: a point per case patient per peptide, not per peptide")
    ap.add_argument("--no-stat", dest="no_stat", action="store_true",
                    help="--box: leave the per-panel statistic off")
    ap.add_argument("--box-patients", type=int, default=0,
                    help="instead of --min-group per peptide: patients per group "
                         "each box's points must span (0 = off)")
    ap.add_argument("--min-pep", type=int, default=2)
    ap.add_argument("--min-digests", type=int, default=2)
    ap.add_argument("--precursor-q", type=float, default=0.01)
    ap.add_argument("--min-run-peptides", type=int, default=1000)
    ap.add_argument("--top", type=int, default=30)
    ap.add_argument("--max-isoforms", type=int, default=1,
                    help="isoform columns per gene, best-supported first")
    ap.add_argument("--p-unit", dest="p_unit", default="patients",
                    choices=("patients", "patient-ratio", "peptides",
                             "vs-canonical", "lbd-points", "peptide-model"),
                    help="replication unit for the p-value and the quoted effect")
    ap.add_argument("--hue", action="store_true",
                    help="colour each point by the protease that produced it")
    ap.add_argument("--panel-w", type=float, default=None,
                    help="canvas width")
    ap.add_argument("--panel-h", type=float, default=None,
                    help="canvas height")
    ap.add_argument("--subtitle", default=True,
                    action=argparse.BooleanOptionalAction,
                    help="the line under the title saying what a point is")
    ap.add_argument("--ylim", default=None, metavar="LO,HI",
                    help="fix the y-scale, e.g. -2,2")
    ap.add_argument("--note", default=True, action=argparse.BooleanOptionalAction,
                    help="footnote under the panel")
    ap.add_argument("--key", default=True, action=argparse.BooleanOptionalAction,
                    help="draw the protease key")
    ap.add_argument("--box", action="store_true",
                    help="per-patient boxes, Control beside the case")
    ap.add_argument("--centre", default=True,
                    action=argparse.BooleanOptionalAction,
                    help="centre each peptide on its control mean (--no-centre: log2 vs run median only)")
    ap.add_argument("--unit", default="patient",
                    choices=("patient", "run", "peptide", "precursor", "fc",
                             "protein", "protein-peptide"),
                    help="--box only: what one point is")
    ap.add_argument("--isoforms", default=None,
                    help="comma-separated isoform columns to draw, as the screen names them")
    ap.add_argument("--stat", choices=("q", "p"), default="q",
                    help="--unit fc: print the screen's BH q or the raw Welch p")
    ap.add_argument("--max-p", dest="max_p", type=float, default=None,
                    help="--unit protein: draw only proteins with p below this")
    ap.add_argument("--per-row", dest="per_row", type=int, default=None,
                    help="--box only: facets per row; default is one row")
    ap.add_argument("--volcano", default=None, metavar="SVG",
                    help="--screen only: volcano of every screened isoform")
    ap.add_argument("--q-cut", dest="q_cut", type=float, default=0.05)
    ap.add_argument("--title", default=None,
                    help="panel title; --box only")
    ap.add_argument("--no-title", dest="title", action="store_const", const="",
                    help="draw no panel title; --box only")
    ap.add_argument("--letter", default="c")
    ap.add_argument("--out")
    ap.add_argument("--tsv")
    ap.add_argument("--font", default=FONT)
    args = ap.parse_args(argv)

    paths = sorted(p for pat in args.reports for p in glob.glob(pat))
    if not paths:
        sys.exit("no reports matched")
    from lib_fasta import gene_map
    seqs = read_fasta(args.fasta)
    genes = gene_map(args.fasta)
    cond = read_metadata(args.metadata)
    idx, med, bygene = read_reports(paths, seqs, genes, args.precursor_q,
                                    args.min_run_peptides,
                                    key="precursor" if args.unit == "precursor"
                                    else "peptide")
    if args.common_median:
        med = common_medians(idx, med)
    if args.no_run_median:
        # Spectronaut's cross-run normalisation already applied; keep only the run gate
        med = {k: 0.0 for k in med}
    rows = screen(idx, med, bygene, seqs, genes, cond, args.case, args.control,
                  args.min_group, args.min_pep, args.min_digests, args.p_unit,
                  args.box_patients, args.centre, args.fasta_isoforms)
    best_q = min((r["q"] for r in rows if r["q"] == r["q"]), default=float("nan"))
    note = (f"{len(rows)} isoforms tested at ≥ {args.min_pep} diagnostic "
            f"peptides in ≥ {args.min_digests} digests; "
            + (f"none survives Benjamini–Hochberg (best q = {best_q:.2f}), so "
               f"these are the best-supported candidates, not hits"
               if best_q == best_q and best_q > 0.05
               else f"best q = {best_q:.2f}"))

    if args.screen:
        print(f"\n  {len(rows)} isoform(s) with >= {args.min_pep} diagnostic "
              f"peptides in >= {args.min_digests} digests\n")
        print(f"  {'gene':<10}{'isoform':<22}{'log2FC':>8}{'sd':>6}{'dig':>4}"
              f"{'pep':>5}{'p':>10}{'q':>8}   {'canonical':>10}")
        for r in rows[:args.top]:
            print(f"  {r['gene']:<10}{_label(r['isoform'])[:21]:<22}"
                  f"{r['log2fc']:>+8.2f}{r['sd']:>6.2f}"
                  f"{len(r['digests'].split(',')):>4}{r['n_pep']:>5}"
                  f"{r['p']:>10.2e}{r['q']:>8.2f}   "
                  f"{r['canonical_log2fc']:>+10.2f}")
        if args.tsv:
            os.makedirs(os.path.dirname(os.path.abspath(args.tsv)) or ".",
                        exist_ok=True)
            with open(args.tsv, "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                keys = ["gene", "isoform", "canonical", "digests", "n_pep",
                        "log2fc", "sd", "p", "q", "canonical_log2fc",
                        "canonical_digests"]
                w.writerow(keys)
                for r in rows:
                    w.writerow([f"{r[k]:.4g}" if isinstance(r[k], float)
                                else r[k] for k in keys])
            print(f"\n  wrote {args.tsv}")
        if args.volcano:
            volcano(rows, args.volcano, args.font, args.q_cut, args.letter,
                    title=args.title,
                    xlab=("Δ log2 vs run median, LBD (isoform − canonical)"
                          if args.p_unit == "lbd-points" and not args.centre
                          else "Δ log2 LBD / Control (isoform − canonical)"
                          if args.p_unit in ("vs-canonical", "lbd-points", "peptide-model")
                          else "log2 LBD / Control"))
        return 0

    if not args.genes:
        sys.exit("--genes or --screen")
    if args.box and args.unit in ("protein", "protein-peptide"):
        # whole-protein abundance per patient, Welch LBD against Control
        from scipy import stats
        boxes = []
        for gene in [g.strip().upper() for g in args.genes.split(",")]:
            groups = bygene.get(gene)
            base = None
            if groups:
                base, _ = diagnostic(gene, groups, seqs, genes, idx)
            if not base or not base[1]:
                print(f"  {gene}: not quantified at protein level, skipped")
                continue
            v = patient_ratios(base[1], idx, med, cond, args.control)
            a = [x for s_, x in v.items() if cond.get(s_) == args.case]
            b = [x for s_, x in v.items() if cond.get(s_) == args.control]
            if len(a) < 3 or len(b) < 3:
                print(f"  {gene}: under 3 patients in a group, skipped")
                continue
            t = stats.ttest_ind(a, b, equal_var=False)
            d = statistics.mean(a) - statistics.mean(b)
            pep = args.unit == "protein-peptide"
            if pep:
                # every peptide observation is a point (not independent)
                obs = precursor_values(base[1], idx, med, cond, args.control)
                oa = [x for _k, s_, _d, x in obs if cond.get(s_) == args.case]
                ob = [x for _k, s_, _d, x in obs if cond.get(s_) == args.control]
                t = stats.ttest_ind(oa, ob, equal_var=False)
                d = statistics.mean(oa) - statistics.mean(ob)
            print(f"  {gene}: {base[0]} {len(base[1])} peptide units, "
                  f"LBD - Control {d:+.2f}, p {t.pvalue:.3g} "
                  + (f"({len(oa)} v {len(ob)} peptide observations)" if pep
                     else f"(n {len(a)} v {len(b)})"))
            if args.max_p is not None and not t.pvalue < args.max_p:
                print(f"    not drawn: p >= {args.max_p:g}")
                continue
            # --hue: a point per run, coloured by protease
            pts = (obs if pep else
                   run_values(base[1], idx, med, cond, args.control)
                   if args.hue else
                   [(s_, s_, None, x) for s_, x in sorted(v.items())])
            boxes.append((gene, [(base[0], pts, True)],
                          (d, float(t.pvalue), len(a), len(b))))
        ylim = tuple(float(x) for x in args.ylim.split(",")) if args.ylim \
            else None
        box_panel(boxes, cond, args.out, args.font, args.letter, ylim=ylim,
                  panel_w=args.panel_w, panel_h=args.panel_h,
                  case_label=args.case, control_label=args.control,
                  title=args.title, subtitle=False, unit=args.unit,
                  per_row=args.per_row)
        return 0
    picks = []
    for gene in [g.strip().upper() for g in args.genes.split(",")]:
        groups = bygene.get(gene)
        if not groups:
            print(f"  {gene}: no reported group passed the filters")
            continue
        base, diag = diagnostic(gene, groups, seqs, genes, idx, args.fasta_isoforms)
        if base is None:
            print(f"  {gene}: no single canonical entry, skipped")
            continue
        can = column(base[1], idx, med, cond, args.case, args.control,
                     args.min_group, args.p_unit)
        keep = []
        for who, mem in sorted(diag.items()):
            st = column(mem, idx, med, cond, args.case, args.control,
                        args.min_group, args.p_unit)
            if st and args.p_unit == "vs-canonical" and can:
                apply_vs_canonical(can, st)
            if testable(st, can, args.min_pep, args.box_patients) \
                    and len(st["per_digest"]) >= args.min_digests:
                keep.append((";".join(who), st))
        if can is None or not keep:
            print(f"  {gene}: no isoform with >= {args.min_pep} diagnostic "
                  f"peptides in >= {args.min_digests} digests, skipped")
            continue
        # best-supported first: peptide count, then agreement between digests
        keep.sort(key=lambda t: (-t[1]["npep"], t[1]["sd"]))
        if args.isoforms:
            keep = [k for k in keep if k[0] in set(args.isoforms.split(","))]
        if args.max_isoforms:
            keep = keep[:args.max_isoforms]
        members = {";".join(who): mem for who, mem in diag.items()}
        picks.append({"gene": gene, "canon_acc": base[0], "canon": can,
                      "diag": keep, "canon_members": base[1],
                      "diag_members": [(w, members[w]) for w, _st in keep]})
        print(f"  {gene}: canonical {base[0]} {can['mean']:+.2f} over "
              f"{len(can['per_digest'])} digest(s), {can['npep']} peptides")
        for who, st in keep:
            print(f"    {who:<24}{st['mean']:>+7.2f} ± {st['sd']:.2f} over "
                  f"{len(st['per_digest'])} digest(s), {st['npep']} diagnostic "
                  f"peptides, p {st['p']:.2e}")
    if not picks:
        sys.exit("nothing to draw")
    if args.tsv:
        os.makedirs(os.path.dirname(os.path.abspath(args.tsv)) or ".",
                    exist_ok=True)
        with open(args.tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["gene", "column", "kind", "digest", "peptide",
                        "log2fc_case_minus_control"])
            for g in picks:
                for pep, dig, v in g["canon"]["points"]:
                    w.writerow([g["gene"], g["canon_acc"], "canonical", dig,
                                pep, f"{v:.4f}"])
                for who, st in g["diag"]:
                    for pep, dig, v in st["points"]:
                        w.writerow([g["gene"], who, "diagnostic", dig, pep,
                                    f"{v:.4f}"])
        print(f"  wrote {args.tsv}")
    ylim = tuple(float(v) for v in args.ylim.split(",")) if args.ylim \
        else None
    if args.out and args.box and args.unit == "fc":
        if args.p_unit != "vs-canonical":
            sys.exit("--unit fc draws the vs-canonical Welch: add "
                     "--p-unit vs-canonical")
        # the screen's own q
        qof = {(r["gene"], r["isoform"]): r["q"] for r in rows}
        boxes, box_n = [], {}
        for g in picks:
            who, st = g["diag"][0]
            if args.box_n:
                box_n[(g["gene"], g["canon_acc"])] = g["canon"]["box"]
                box_n[(g["gene"], who)] = st["box"]
            pt = lambda col: [(pp, None, dg, v) for pp, dg, v in col["points"]]
            cols = [(g["canon_acc"], pt(g["canon"]), True), (who, pt(st), False)]
            if args.per_patient:
                # one point per case patient per peptide, against the peptide's control mean
                cols = []
                for lab, mem, is_can in ((g["canon_acc"], g["canon_members"], True),
                                         (who, g["diag_members"][0][1], False)):
                    vals = precursor_values(mem, idx, med, cond, args.control,
                                            args.centre)
                    cols.append((lab, [(f"{pp}{s_}", s_, dg, v) for pp, s_, dg, v in vals
                                       if cond.get(s_) == args.case], is_can))
                    box_n[(g["gene"], lab)] = (
                        len({s_ for _p, s_, _d, _v in vals if cond.get(s_) == args.case}),
                        len({s_ for _p, s_, _d, _v in vals if cond.get(s_) == args.control}))
            boxes.append((g["gene"], cols,
                          (st["mean"],
                           st["p"] if args.stat == "p"
                           else qof.get((g["gene"], who), float("nan")),
                           0, 0)))
        box_panel(boxes, cond, args.out, args.font, args.letter, ylim=ylim,
                  panel_w=args.panel_w, panel_h=args.panel_h,
                  case_label=args.case, control_label=args.control,
                  title=args.title, subtitle=args.subtitle, unit="fc",
                  per_row=args.per_row, stat=args.stat, box_n=box_n,
                  show_stat=not args.no_stat, centre=args.centre)
    elif args.out and args.box:
        boxes = []
        for g in picks:
            cols, con = box_columns(g, idx, med, cond, args.case,
                                    args.control, args.centre,
                                    args.unit)
            boxes.append((g["gene"], cols, con))
            print(f"  {g['gene']}: patient-level Δ isoform - canonical "
                  f"{con[0]:+.2f}, p {con[1]:.3g} (n {con[2]} v {con[3]})")
        box_panel(boxes, cond, args.out, args.font, args.letter, ylim=ylim,
                  panel_w=args.panel_w, panel_h=args.panel_h,
                  case_label=args.case, control_label=args.control,
                  title=args.title, subtitle=args.subtitle,
                  centre=args.centre,
                  unit=args.unit, per_row=args.per_row,
                  box_n=args.box_n, show_stat=not args.no_stat)
    elif args.out:
        panel(picks, args.out, args.font, args.letter, key=args.key, ylim=ylim,
              subtitle=args.subtitle, hue=args.hue, panel_w=args.panel_w,
              panel_h=args.panel_h, case_label=args.case,
              control_label=args.control,
              note=note if args.note else None)
    return 0


if __name__ == "__main__":
    sys.exit(main())
