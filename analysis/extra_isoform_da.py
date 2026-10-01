#!/usr/bin/env python3
"""Differential abundance of isoforms, relative to their own canonical form.

    python3 extra_isoform_da.py quant --metadata data/metadata.xlsx \
        --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta --outdir da_iso
"""

import argparse
import collections
import csv
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_svg import Canvas                                   # noqa: E402
from lib_palette import (AXIS, FONT, GRID, INK, INK_MUTED, INK_SECONDARY,
                     UNION, assign, display)                         # noqa: E402
from extra_differential import (ORDER, label_for, load_matrix, normalise,  # noqa: E402
                                outlier_runs, qc_runs, test_digest_moderated)
from lib_fasta import gene_map                                     # noqa: E402

ISO = re.compile(r"-\d+$")


def classify(group):
    """-> ('isoform'|'canonical'|None, base_accession, key)."""
    accs = [a.strip() for a in group.split(";") if a.strip()]
    if not accs:
        return None, None, None
    bases = {ISO.sub("", a) for a in accs}
    if len(bases) != 1:
        return None, None, None          # spans genes: not a clean pair
    base = bases.pop()
    if all(ISO.search(a) for a in accs):
        return "isoform", base, ";".join(sorted(accs))
    if len(accs) == 1:
        return "canonical", base, accs[0]
    return None, None, None              # canonical + isoform together: shared


def pairs_for_digest(logs, min_runs):
    """-> {(base, isoform_key): {sample: log2 isoform/canonical}}"""
    iso, can = collections.defaultdict(list), {}
    for g, d in logs.items():
        if len(d) < min_runs:
            continue
        kind, base, key = classify(g)
        if kind == "isoform":
            iso[base].append((key, d))
        elif kind == "canonical":
            can[base] = d
    out = {}
    for base, entries in iso.items():
        if base not in can:
            continue
        cd = can[base]
        for key, idict in entries:
            shared = [s for s in idict if s in cd]
            if len(shared) >= min_runs:
                out[(base, key)] = {s: idict[s] - cd[s] for s in shared}
    return out


def combine(per_digest, min_digests):
    """Average the ratio across digests per patient, centring within digest."""
    import numpy as np
    acc = collections.defaultdict(lambda: collections.defaultdict(list))
    seen = collections.defaultdict(set)
    for dig, pairs in per_digest.items():
        for k, d in pairs.items():
            if len(d) < 2:
                continue
            mu = float(np.mean(list(d.values())))
            seen[k].add(dig)
            for s, v in d.items():
                acc[k][s].append(v - mu)
    return ({k: {s: float(np.mean(v)) for s, v in d.items()}
             for k, d in acc.items() if len(seen[k]) >= min_digests},
            {k: sorted(v) for k, v in seen.items()})


def nice(hi, target=5):
    if hi <= 0:
        return [0]
    raw = hi / target
    mag = 10 ** math.floor(math.log10(raw))
    step = next((m * mag for m in (1, 2, 2.5, 5, 10) if raw <= m * mag), 10 * mag)
    return [i * step for i in range(int(hi / step) + 2) if i * step <= hi * 1.001]


def panel_access(per_digest, union_keys, outdir, font, letter="a",
                 width=1215.0, ts=1.7):
    """How many isoform/canonical pairs each digest can measure."""
    colour = assign(ORDER + ["All"])
    counts = [(d, len(per_digest.get(d, {})), colour[d]) for d in ORDER]
    tryp = set(per_digest.get("Trypsin", {}))
    only_other = len(union_keys - tryp)
    counts.append(("All", len(union_keys), colour["All"]))
    W, H = 560, 300
    c = Canvas(W, H, font, font_scale=0.78 * ts / 1.7, out_w=width)
    c.text(20, 30, letter, 13, INK, "start", "600")
    c.text(64, 30, "Isoforms measurable against their own canonical form", 11.5,
           INK, "start", "600")
    ml, mt = 64.0, 62.0
    pw, ph = W - ml - 30, H - mt - 58
    hi = max(v for _, v, _ in counts) * 1.16
    for t in nice(hi):
        y = mt + ph - t / hi * ph
        c.line(ml, y, ml + pw, y, stroke=GRID, sw=1)
        c.text(ml - 8, y + 3.5, f"{t:,.0f}", 9, INK_MUTED, "end")
    bw = pw / len(counts) * 0.46
    for i, (name, v, col) in enumerate(counts):
        x = ml + pw * (i + 0.5) / len(counts) - bw / 2
        y = mt + ph - v / hi * ph
        c.rect(x, y, bw, mt + ph - y, fill=col, fo=0.85, rx=3)
        c.text(x + bw / 2, y - 8, f"{v:,}", 9.5, INK, "middle")
        c.text(x + bw / 2, mt + ph + 16, name, 9.5, INK, "middle")
        if name == "All":
            # split the union bar: trypsin-reachable part vs second-protease-only part
            yb = mt + ph - (len(union_keys) - only_other) / hi * ph
            c.rect(x, yb, bw, mt + ph - yb, fill="#ffffff", rx=3)
            c.rect(x, yb, bw, mt + ph - yb, fill=col, fo=0.30, rx=3)
            c.rect(x, y, bw, yb - y, fill=col, fo=0.95, rx=3)
            c.line(x - 6, (y + yb) / 2, x - 2, (y + yb) / 2, stroke=INK_MUTED, sw=1)
            c.text(x - 9, (y + yb) / 2 - 2, f"{only_other:,} only with", 8.4,
                   INK_SECONDARY, "end")
            c.text(x - 9, (y + yb) / 2 + 9, "Glu-C or Lys-C", 8.4,
                   INK_SECONDARY, "end")
    c.line(ml, mt + ph, ml + pw, mt + ph, stroke=AXIS, sw=1)
    c.add(f'<g transform="translate(20 {mt + ph / 2:.1f}) rotate(-90)">'
          f'<text x="0" y="0" font-size="10" text-anchor="middle" fill="{INK}">'
          f'Isoform / canonical pairs</text></g>')
    p = os.path.join(outdir, "isoform_access.svg")
    with open(p, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {p}")
    return only_other


def panel_volcano(res, seen, genes, outdir, font, qcut, letter="b",
                  width=1215.0, ts=1.7):
    colour = assign(ORDER + ["All"])
    W, H = 560, 330
    c = Canvas(W, H, font, font_scale=0.78 * ts / 1.7, out_w=width)
    c.text(20, 30, letter, 13, INK, "start", "600")
    head = "Differential isoform usage"
    c.text(64, 30, head, 11.5, INK, "start", "600")
    c.text(64 + len(head) * 7.2 + 16, 30, f"{len(res):,} pairs tested", 9,
           INK_MUTED, "start")
    ml, mt = 64.0, 62.0
    pw, ph = W - ml - 30, H - mt - 62
    xlim = min(max((abs(r["log2fc"]) for r in res), default=2), 5)
    ymax = max((-math.log10(r["p"]) for r in res), default=1)

    def X(v):
        return ml + (max(-xlim, min(xlim, v)) + xlim) / (2 * xlim) * pw

    def Y(v):
        return mt + ph - min(v, ymax) / ymax * ph

    for t in nice(ymax, 4):
        c.line(ml, Y(t), ml + pw, Y(t), stroke=GRID, sw=1)
        c.text(ml - 8, Y(t) + 3.5, f"{t:.0f}", 9, INK_MUTED, "end")
    c.line(X(0), mt, X(0), mt + ph, stroke=AXIS, sw=1)
    c.line(ml, mt + ph, ml + pw, mt + ph, stroke=AXIS, sw=1)
    for v in (-xlim, 0, xlim):
        c.text(X(v), mt + ph + 16, f"{v:.0f}", 9, INK_MUTED, "middle")
    for r in sorted(res, key=lambda r: r["p"], reverse=True):
        multi = "Trypsin" not in seen[r["key"]]
        col = UNION if multi else INK_MUTED
        c.add(f'<circle cx="{X(r["log2fc"]):.1f}" cy="{Y(-math.log10(r["p"])):.1f}" '
              f'r="{3 if multi else 2}" fill="{col}" '
              f'fill-opacity="{0.85 if multi else 0.3}"/>')
    placed = []
    for r in sorted(res, key=lambda r: r["p"])[:6]:
        multi = "Trypsin" not in seen[r["key"]]
        ty = Y(-math.log10(r["p"])) + 3
        while any(abs(ty - q) < 10 for q in placed):
            ty += 10
        placed.append(ty)
        c.text(X(r["log2fc"]) + (7 if r["log2fc"] > 0 else -7), ty,
               label_for(r["key"][1], genes) + ("*" if multi else ""), 8,
               INK_SECONDARY, "start" if r["log2fc"] > 0 else "end")
    c.rect(ml + 6, mt + ph - 20, 9, 9, fill=UNION, rx=2)
    c.text(ml + 20, mt + ph - 12, "only measurable with Glu-C or Lys-C", 8.4,
           INK_SECONDARY, "start")
    c.text(ml + pw / 2, H - 26, "log2 change in isoform / canonical ratio", 10,
           INK, "middle")
    c.add(f'<g transform="translate(20 {mt + ph / 2:.1f}) rotate(-90)">'
          f'<text x="0" y="0" font-size="10" text-anchor="middle" fill="{INK}">'
          f'−log10 p</text></g>')
    p = os.path.join(outdir, "isoform_volcano.svg")
    with open(p, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {p}")


def panel_examples(res, seen, per_digest, cond, case, control, genes, outdir,
                   font, n=4, letter="c", width=1215.0, ts=1.7,
                   case_label=None, control_label=None):
    """Per-patient isoform ratios for the strongest multi-protease-only pairs."""
    case_label = case_label or "LBD"
    control_label = control_label or "Control"
    import numpy as np
    colour = assign(ORDER + ["All"])
    picks = [r for r in sorted(res, key=lambda r: r["p"])
             if "Trypsin" not in seen[r["key"]]][:n]
    if not picks:
        picks = sorted(res, key=lambda r: r["p"])[:n]
    pw, ph, gapx = 220.0, 150.0, 34.0
    ml, mt = 60.0, 66.0
    strip_h = 20.0
    W = ml + len(picks) * pw + (len(picks) - 1) * gapx + 20
    H = mt + ph + 26
    c = Canvas(W, H, font, font_scale=1.45 * ts / 1.7, out_w=width)
    c.text(20, 30, letter, 13, INK, "start", "600")
    c.text(ml, 30, f"Isoform / canonical ratio per patient", 11.5, INK,
           "start", "600")
    for k, r in enumerate(picks):
        key = r["key"]
        x0 = ml + k * (pw + gapx)
        digs = [d for d in ORDER if key in per_digest.get(d, {})]
        vals = [v for d in digs for v in per_digest[d][key].values()]
        lo, hi = min(vals), max(vals)
        # extra top padding: the accession and p-value sit inside the facet
        rng = (hi - lo) or 1.0
        lo, hi = lo - rng * 0.12, hi + rng * 0.40

        def Y(v):
            return mt + ph - (v - lo) / (hi - lo) * ph

        name = label_for(key[1], genes)
        c.rect(x0, mt - strip_h, pw, strip_h, fill=GRID, stroke=INK_SECONDARY,
               sw=1.0)
        c.text(x0 + pw / 2, mt - 6, name, 10.5, INK, "middle", "600")
        for t in range(int(math.floor(lo)), int(math.ceil(hi)) + 1):
            if lo <= t <= hi:
                c.line(x0, Y(t), x0 + pw, Y(t), stroke=GRID, sw=1)
                c.text(x0 - 6, Y(t) + 3.5, f"{t:g}", 8.4, INK_MUTED, "end")
        c.rect(x0, mt, pw, ph, fill="none", stroke=INK_SECONDARY, sw=1.0)
        c.text(x0 + pw / 2, mt + 14, key[1], 8.4, INK_MUTED, "middle")
        c.text(x0 + pw / 2, mt + 26, f"p = {r['p']:.1e}", 8.6, INK, "middle",
               italic=True)
        for gi, (grp, lab) in enumerate(((control, control_label),
                                         (case, case_label))):
            xs = x0 + pw * (0.30 + 0.40 * gi)
            for di, dig in enumerate(digs):
                pts = [v for s_, v in per_digest[dig][key].items()
                       if cond.get(s_) == grp]
                if not pts:
                    continue
                cx = xs + (di - (len(digs) - 1) / 2) * 13.0
                for j, v in enumerate(sorted(pts)):
                    c.add(f'<circle cx="{cx + ((j % 3) - 1) * 4.0:.1f}" '
                          f'cy="{Y(v):.1f}" r="4.0" fill="{colour[dig]}" '
                          f'fill-opacity="0.55" stroke="{darker(colour[dig])}" '
                          f'stroke-width="1.1"/>')
                c.line(cx - 13, Y(float(np.median(pts))), cx + 13,
                       Y(float(np.median(pts))), stroke=INK, sw=3.2)
            c.text(xs, mt + ph + 15, lab, 9.0, INK, "middle")
    c.add(f'<g transform="translate(20 {mt + ph / 2:.1f}) rotate(-90)">'
          f'<text x="0" y="0" font-size="10" text-anchor="middle" fill="{INK}">'
          f'log2 isoform / canonical</text></g>')
    p = os.path.join(outdir, "isoform_examples.svg")
    with open(p, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {p}  ({', '.join(label_for(r['key'][1], genes) for r in picks)})")


def darker(hex_fill, f=0.62):
    """`hex_fill` scaled toward black, for a point's outline."""
    r, g, b = (int(hex_fill[i:i + 2], 16) for i in (1, 3, 5))
    return "#" + "".join(f"{round(v * f):02x}" for v in (r, g, b))


def peptides_for(reports, bases, precursor_q):
    """-> {base: {'canonical'|'isoform:<key>'|'shared': {digest: set(peptides)}}}"""
    import pyarrow.parquet as pq
    out = collections.defaultdict(lambda: collections.defaultdict(
        lambda: collections.defaultdict(set)))
    bases = set(bases)
    import lib_report as rp
    for r in rp.open_reports(reports, order=ORDER):
        path, dig = r.path, r.protease
        for b in r.batches(["protein_groups", "peptide", "precursor_q"]):
            g = b["protein_groups"]
            s_ = b["peptide"]
            q = b["precursor_q"]
            for i in range(b["_n"]):
                if not g[i] or not s_[i] or q[i] is None or q[i] > precursor_q:
                    continue
                accs = [a.strip() for a in g[i].split(";") if a.strip()]
                hit = {ISO.sub("", a) for a in accs} & bases
                if len(hit) != 1:
                    continue
                base = hit.pop()
                if all(ISO.search(a) for a in accs):
                    slot = "isoform:" + ";".join(sorted(accs))
                elif len(accs) == 1:
                    slot = "canonical"
                else:
                    slot = "shared"
                out[base][slot][dig].add(s_[i])
        print(f"  scanned {os.path.basename(path)[:44]:<44s} {dig}")
    return out


def align_blocks(a, b):
    """Matching blocks between two sequences -> [(ia, ib, size), ...]."""
    import difflib
    sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
    return [(m.a, m.b, m.size) for m in sm.get_matching_blocks() if m.size]


def unique_spans(seq_len, blocks, which):
    """Spans of one sequence not matched in the other."""
    idx = 0 if which == "a" else 1
    covered = sorted((b[idx], b[idx] + b[2]) for b in blocks)
    out, pos = [], 0
    for s, e in covered:
        if s > pos:
            out.append((pos, s))
        pos = max(pos, e)
    if pos < seq_len:
        out.append((pos, seq_len))
    return out


def panel_strips(picks, seqs, pep_by_base, genes, outdir, font, letter="d",
                 flank=12, max_tiers=6):
    """Canonical and isoform as aligned strips, with the peptides on them."""
    colour = assign(ORDER + ["All"])
    W = 1010.0
    lab_w, call_w = 96.0, 250.0
    strip_x = lab_w + 18
    strip_w = W - strip_x - call_w - 26
    pep_h, tier_gap, strip_h = 6.0, 1.6, 14.0
    form_h = strip_h + 8 + max_tiers * (pep_h + tier_gap) + 6
    pair_gap = 48.0
    mt = 70.0
    H = mt + len(picks) * (form_h * 2 + pair_gap) + 16
    c = Canvas(W, H, font)
    c.text(20, 30, letter, 13, INK, "start", "600")
    c.text(strip_x, 30, "Which peptides tell the isoforms apart", 11.5, INK,
           "start", "600")
    lx = strip_x + 262
    for lbl, col, op in [("shared", INK_MUTED, 0.35)] + \
                        [(d, colour[d], 0.9) for d in ORDER]:
        c.rect(lx, 22, 9, 9, fill=col, fo=op, rx=2)
        c.text(lx + 13, 30, lbl, 9, INK, "start")
        lx += 58 if lbl == "shared" else 56

    def draw(seq, layers, x0, y, scale, highlight=None):
        """layers: [(peptides, colour, opacity, outline)] drawn into tiers."""
        ends, extra = [], 0
        for peps, col, op, outline in layers:
            for pep in sorted(peps, key=lambda p: -len(p)):
                i = seq.find(pep)
                if i < 0:
                    continue
                row = 0
                while row < len(ends) and ends[row] > i:
                    row += 1
                if row >= max_tiers:
                    extra += 1
                    continue
                if row == len(ends):
                    ends.append(0)
                ends[row] = i + len(pep) + 2 / max(scale, 1e-9)
                hot = highlight is not None and pep == highlight
                c.rect(x0 + i * scale, y + row * (pep_h + tier_gap),
                       max(len(pep) * scale, 1.2), pep_h, fill=col,
                       fo=1.0 if hot else op, rx=1.5,
                       stroke=INK if (hot or outline) else "none",
                       sw=1.1 if hot else 0)
        return extra

    for k, (base, iso_key, note) in enumerate(picks):
        iso_acc = iso_key.split(";")[0]
        if base not in seqs or iso_acc not in seqs:
            continue
        ca, ia = seqs[base], seqs[iso_acc]
        blocks = align_blocks(ca, ia)
        scale = strip_w / max(len(ca), len(ia))
        y0 = mt + k * (form_h * 2 + pair_gap)
        slots = pep_by_base.get(base, {})
        shared = slots.get("shared", {})
        canon_only = slots.get("canonical", {})
        iso_only = slots.get("isoform:" + iso_key, {})

        iso_unique = unique_spans(len(ia), blocks, "b")

        def informative(pep):
            """Rank: overlaps isoform-unique sequence first, then length."""
            i = ia.find(pep)
            if i < 0:
                return (-1, 0)
            ov = sum(max(0, min(i + len(pep), e) - max(i, s0))
                     for s0, e in iso_unique)
            return (ov, len(pep))

        cand = sorted({p for d in ORDER for p in iso_only.get(d, ())
                       if ia.find(p) >= 0}, key=informative, reverse=True)
        key_pep = cand[0] if cand else None

        for name, seq, which, y, spec in ((base, ca, "a", y0, canon_only),
                                          (iso_key, ia, "b", y0 + form_h,
                                           iso_only)):
            if which == "a":
                c.text(lab_w, y + 12, label_for(name, genes), 10.5, INK, "end",
                       "600")
            c.text(lab_w, y + (12 if which == "b" else 25), name.split(";")[0],
                   9, INK_SECONDARY, "end")
            n_spec = sum(len(v) for v in spec.values())
            n_sh = sum(len(v) for v in shared.values())
            c.text(lab_w, y + (24 if which == "b" else 37),
                   f"{len(seq)} aa", 8.4, INK_MUTED, "end")
            c.text(lab_w, y + (36 if which == "b" else 49),
                   f"{n_spec} specific · {n_sh} shared", 8, INK_MUTED, "end")
            c.rect(strip_x, y + 4, len(seq) * scale, strip_h, fill="#ffffff",
                   stroke=AXIS, sw=1, rx=3)
            for s0, e0 in unique_spans(len(seq), blocks, which):
                c.rect(strip_x + s0 * scale, y + 4, max((e0 - s0) * scale, 1.5),
                       strip_h, fill=INK_MUTED, fo=0.32)
            # discriminating peptides claim tiers first
            layers = []
            for d in ORDER:
                if spec.get(d):
                    layers.append((spec[d], colour[d], 0.9, False))
            if shared:
                layers.append((set().union(*shared.values()), INK_MUTED, 0.35,
                               False))
            n_extra = draw(seq, layers, strip_x, y + strip_h + 10, scale,
                           key_pep if which == "b" else None)
            if n_extra:
                c.text(lab_w, y + (46 if which == "b" else 59),
                       f"{n_extra} not drawn", 8, INK_MUTED, "end")

        yb = y0 + form_h * 2 - 2
        step = 200 if max(len(ca), len(ia)) > 600 else 50
        for r in range(0, max(len(ca), len(ia)) + 1, step):
            c.line(strip_x + r * scale, yb, strip_x + r * scale, yb + 4,
                   stroke=AXIS, sw=1)
            c.text(strip_x + r * scale, yb + 14, str(r), 8, INK_MUTED, "middle")

        cx = W - call_w - 8
        if key_pep:
            i = ia.find(key_pep)
            a0, a1 = max(0, i - flank), min(len(ia), i + len(key_pep) + flank)
            head = f"{iso_acc}  {i + 1}–{i + len(key_pep)}"
            mark = set(range(i, min(i + len(key_pep), a1)))
            dig = next(d for d in ORDER if key_pep in iso_only.get(d, ()))
        else:
            spans = unique_spans(len(ia), blocks, "b") or [(0, min(len(ia), 20))]
            s0, e0 = max(spans, key=lambda t: t[1] - t[0])
            a0, a1 = max(0, s0 - flank), min(len(ia), e0 + flank)
            head = f"{iso_acc}  {s0 + 1}–{e0}"
            mark, dig = set(range(s0, e0)), None
        # cap the window so the glyphs stay legible
        MAXW = 40
        if a1 - a0 > MAXW:
            centre = (min(mark) + max(mark)) // 2 if mark else (a0 + a1) // 2
            a0n = max(a0, centre - MAXW // 2)
            a1n = min(a1, a0n + MAXW)
            a0n = max(a0, a1n - MAXW)
            trimmed = (a0n > a0, a1n < a1)
            a0, a1 = a0n, a1n
        else:
            trimmed = (False, False)
        c.text(cx, y0 + 14, head, 9, INK, "start", "600")
        pitch = min(8.6, (call_w - 12) / max(a1 - a0, 1))
        fs = max(4.6, min(8.0, pitch * 1.18))
        for j in range(a0, a1):
            x = cx + (j - a0) * pitch
            if j in mark:
                c.rect(x, y0 + 22, pitch, 14,
                       fill=colour[dig] if dig else INK_MUTED,
                       fo=0.85 if dig else 0.3, rx=1.5)
            c.text(x + pitch / 2, y0 + 33, ia[j], fs,
                   "#ffffff" if (j in mark and dig) else INK, "middle")
        for side, on in zip(("start", "end"), trimmed):
            if on:
                c.text(cx - 4 if side == "start" else cx + (a1 - a0) * pitch + 4,
                       y0 + 33, "…", 9, INK_MUTED,
                       "end" if side == "start" else "start")
        n_iso = sum(len(v) for v in iso_only.values())
        by = " + ".join(f"{d} {len(v)}" for d, v in sorted(iso_only.items()))
        c.text(cx, y0 + 50, f"{n_iso} isoform-specific peptides ({by})", 8.4,
               INK_SECONDARY, "start")
        if note:
            c.text(cx, y0 + 63, note, 8.4, INK_MUTED, "start")

    p = os.path.join(outdir, "isoform_strips.svg")
    with open(p, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {p}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("quantdir")
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--outdir", default="derived/da_iso")
    ap.add_argument("--case", default="LBD")
    ap.add_argument("--control", default="Control")
    ap.add_argument("--min-runs", type=int, default=4,
                    help="runs a pair must be measurable in, within a digest")
    ap.add_argument("--min-per-group", type=int, default=3)
    ap.add_argument("--min-digests", type=int, default=1)
    ap.add_argument("--min-run-frac", type=float, default=0.5)
    ap.add_argument("--max-outlier-mad", type=float, default=3.0)
    ap.add_argument("--q", type=float, default=0.05)
    ap.add_argument("--reports", default=None,
                   help="report globs, for the strip panel")
    ap.add_argument("--strip-pool", type=int, default=60,
                   help="best-ranked pairs considered for the strip panel")
    ap.add_argument("--strip-min-pep", type=int, default=5,
                   help="isoform-specific peptides a pair needs to be drawn")
    ap.add_argument("--strip-trypsin-blind", action=argparse.BooleanOptionalAction,
                   default=True,
                   help="prefer pairs trypsin alone could not test")
    ap.add_argument("--strips", type=int, default=3,
                   help="how many pairs to draw as aligned strips")
    ap.add_argument("--font", default=FONT)
    args = ap.parse_args(argv)

    from lib_report import read_metadata
    cond = read_metadata(args.metadata)
    genes = gene_map(args.fasta) if args.fasta else {}
    seqs_all = {}
    if args.fasta:
        from lib_fasta import read_fasta
        seqs_all = read_fasta(args.fasta)
    os.makedirs(args.outdir, exist_ok=True)

    per_digest = {}
    for dig in ORDER:
        path = os.path.join(args.quantdir, f"{dig}.csv")
        if not os.path.exists(path):
            continue
        mat, samples = load_matrix(path)
        keep, dropped, _ = qc_runs(mat, samples, args.min_run_frac)
        logs, _, _ = normalise(mat, keep)
        bad = outlier_runs(logs, keep, args.max_outlier_mad)
        if bad:
            keep = [s for s in keep if s not in bad]
            logs, _, _ = normalise(mat, keep)
        pairs = pairs_for_digest(logs, args.min_runs)
        per_digest[dig] = pairs
        print(f"{dig}: {len(keep)} runs after QC"
              + (f" (dropped {', '.join(sorted(set(dropped) | set(bad)))})"
                 if dropped or bad else "")
              + f"; {len(pairs):,} isoform/canonical pairs measurable")

    union_keys = set().union(*[set(p) for p in per_digest.values()])
    tryp = set(per_digest.get("Trypsin", {}))
    print(f"\n{len(union_keys):,} pairs measurable in at least one digest; "
          f"{len(tryp):,} with trypsin; "
          f"{len(union_keys - tryp):,} only with Glu-C and/or Lys-C")

    comb, seen = combine(per_digest, args.min_digests)
    ga = [s for s in cond if cond[s] == args.case]
    gb = [s for s in cond if cond[s] == args.control]
    res = test_digest_moderated(comb, (ga, gb), args.min_per_group)
    for r in res:
        r["key"] = r["group"] if isinstance(r["group"], tuple) else r["group"]
    # test_digest_moderated keys on the dict key, which here is the tuple
    sig = [r for r in res if r["q"] < args.q]
    print(f"\ntested {len(res):,} isoform ratios; {len(sig):,} at q < {args.q:g}; "
          f"{sum(1 for r in res if r['p'] < 0.05):,} at raw p < 0.05")
    multi_only = [r for r in res if "Trypsin" not in seen[r["key"]]]
    print(f"  {len(multi_only):,} of the tested pairs are invisible to trypsin alone")

    out = os.path.join(args.outdir, "isoform_ratios.tsv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "isoform_group", "canonical", "digests",
                    "trypsin_can_see", "log2fc_ratio", "p", "q",
                    "n_case", "n_ctrl"])
        for r in sorted(res, key=lambda r: r["p"]):
            base, key = r["key"]
            w.writerow([label_for(key, genes), key, base,
                        "+".join(seen[r["key"]]),
                        "yes" if "Trypsin" in seen[r["key"]] else "no",
                        f"{r['log2fc']:.4f}", f"{r['p']:.3e}", f"{r['q']:.3e}",
                        r["n_case"], r["n_ctrl"]])
    print(f"  wrote {out}")
    print(f"\n  {'gene':<12}{'isoform':<18}{'log2FC':>8}{'p':>11}   digests")
    for r in sorted(res, key=lambda r: r["p"])[:12]:
        base, key = r["key"]
        mark = "  <- trypsin cannot" if "Trypsin" not in seen[r["key"]] else ""
        print(f"  {label_for(key, genes)[:11]:<12}{key[:17]:<18}"
              f"{r['log2fc']:>+8.2f}{r['p']:>11.1e}   "
              f"{'+'.join(seen[r['key']])}{mark}")

    print("\nfigures:")
    panel_access(per_digest, union_keys, args.outdir, args.font)
    panel_volcano(res, seen, genes, args.outdir, args.font, args.q)
    panel_examples(res, seen, per_digest, cond, args.case, args.control, genes,
                   args.outdir, args.font)

    if args.reports:
        import glob as _glob
        reports = []
        for pat in args.reports.split(","):
            reports.extend(sorted(_glob.glob(pat.strip())) or [pat.strip()])
        # Rank on evidence, not p alone.
        ranked = sorted(res, key=lambda r: r["p"])
        blind = [r for r in ranked if "Trypsin" not in seen[r["key"]]]
        cands = (blind or ranked) if args.strip_trypsin_blind else ranked
        cands = cands[:args.strip_pool]
        # one pass covers both the strip panel and the check below
        bases = {b for b, _ in (k for k in (union_keys - tryp))} | \
                {b for b, _ in (r["key"] for r in cands)}
        print(f"\n  one pass over the reports for {len(bases)} accessions")
        pep_by_base = peptides_for(reports, bases, 0.01)

        def iso_peps(base, key, digests=None):
            d = pep_by_base.get(base, {}).get("isoform:" + key, {})
            if digests:
                d = {k: v for k, v in d.items() if k in digests}
            return sum(len(v) for v in d.values())

        blind_keys = union_keys - tryp
        with_tryp = sum(1 for b, k in blind_keys if iso_peps(b, k, {"Trypsin"}))
        print(f"\n  of the {len(blind_keys):,} pairs trypsin cannot test, "
              f"{len(blind_keys) - with_tryp:,} have no trypsin peptide for the "
              f"isoform at all")
        print(f"  the other {with_tryp:,} do have trypsin peptides — there the "
              f"canonical side is what trypsin failed to quantify cleanly")

        if args.strips:
            good = [r for r in cands
                    if iso_peps(*r["key"]) >= args.strip_min_pep]
            if not good:
                print(f"  no candidate has >= {args.strip_min_pep} "
                      f"isoform-specific peptides; using the best supported")
                good = sorted(cands, key=lambda r: -iso_peps(*r["key"]))
            good = sorted(good, key=lambda r: r["p"])[:args.strips]
            picks = []
            for r in good:
                base, key = r["key"]
                picks.append((base, key,
                              f"{'+'.join(seen[r['key']])} · ratio "
                              f"{r['log2fc']:+.2f}, p {r['p']:.1e}"))
                print(f"    {label_for(key, genes):<12}{key:<18}"
                      f"{iso_peps(base, key):>4} isoform-specific peptides"
                      f"   p {r['p']:.1e}")
            panel_strips(picks, seqs_all, pep_by_base, genes, args.outdir,
                         args.font)
    return 0


if __name__ == "__main__":
    sys.exit(main())
