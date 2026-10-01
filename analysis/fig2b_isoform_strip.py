#!/usr/bin/env python3
"""Splice-event evidence: what distinguishes an isoform, and which peptides prove it.

    # the evidence panel for one pair
    python3 fig2b_isoform_strip.py '*.parquet' --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta \
        --protein P37840 --isoform P37840-3 --label SNCA --out figures/snca_event.svg
"""

import argparse
import collections
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_svg import Canvas, esc                             # noqa: E402
from lib_palette import (AXIS, DIVERGING_HIGH, EMPTY, FONT, GRID, INK, INK_MUTED,
                     INK_SECONDARY, STRIP_FILL, SURFACE, UNION, assign,
                     display, ink_on, TEXT_BOOST)                                 # noqa: E402

# Residue cells: black = shared sequence, red = differing; the band says which form.
SEQ_SAME = INK
SEQ_CAN = DIVERGING_HIGH[2]
SEQ_ISO = DIVERGING_HIGH[2]

ORDER = ["GluC", "LysC", "Trypsin"]
STRIP_TEXT = 1.4
ISO = re.compile(r"-\d+$")
UNMAPPED = "no Ensembl transcript"        # mirrors splice_events.UNMAPPED


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
    """The cheap half: pool differing blocks that sit too close together."""
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


def label(ev, clen, ilen):
    """A label for an already-decomposed set of changes."""
    if not ev:
        return "identical"
    if len(ev) > 1:
        return f"{min(len(ev), 6)}{'+' if len(ev) >= 6 else ''} separate changes"
    kind, cs, ce, _is, ie = ev[0]
    at_n = cs <= 3
    at_c = ce >= clen - 3 if kind != "insert" else ie >= ilen - 3
    where = "N-terminal" if at_n else ("C-terminal" if at_c else "internal")
    verb = {"delete": "segment absent", "insert": "segment added",
            "replace": "segment replaced"}[kind]
    return f"{where} {verb}"


def event_type(c, i, min_equal=MIN_EQUAL):
    """A label for the isoform's relationship to its canonical form."""
    return label(events(c, i, min_equal), len(c), len(i))


def largest_event(c, i, min_equal=MIN_EQUAL):
    """The change to draw: the one involving the most residues."""
    ev = events(c, i, min_equal)
    if not ev:
        return None
    return max(ev, key=lambda e: max(e[2] - e[1], e[4] - e[3]))


def whole_event(c, i, min_equal=MIN_EQUAL):
    """Every change as one, from the first one's start to the last one's end."""
    ev = events(c, i, min_equal)
    if not ev:
        return None
    c0, c1 = min(e[1] for e in ev), max(e[2] for e in ev)
    i0, i1 = min(e[3] for e in ev), max(e[4] for e in ev)
    op = "replace" if c1 > c0 and i1 > i0 else "delete" if c1 > c0 else "insert"
    return (op, c0, c1, i0, i1)


def collect_groups(paths, precursor_q, bases=None, protein_q=0.01):
    """{group: {digest: set(peptides)}} in one pass, optionally filtered."""
    import pyarrow.parquet as pq
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


def split_evidence(groups):
    """-> (isoform accession -> {digest: peps}, base -> {digest: peps})."""
    iso_ev = collections.defaultdict(lambda: collections.defaultdict(set))
    can_ev = collections.defaultdict(lambda: collections.defaultdict(set))
    for g, byd in groups.items():
        accs = [a.strip() for a in g.split(";") if a.strip()]
        bases = {ISO.sub("", a) for a in accs}
        if len(bases) != 1:
            continue
        base = bases.pop()
        if all(ISO.search(a) for a in accs):
            for a in accs:
                for d, v in byd.items():
                    iso_ev[a][d] |= v
        elif len(accs) == 1:
            for d, v in byd.items():
                can_ev[base][d] |= v
    return iso_ev, can_ev


def curated(cache_dir, seqs):
    """-> (types, mechanisms, spans) from the UniProt VAR_SEQ cache, or empty dicts."""
    if not cache_dir:
        return {}, {}, {}
    import prep_uniprot as uv
    path = os.path.join(cache_dir, "varseq.tsv.gz")
    if not os.path.exists(path):
        print(f"  no curation at {path}; falling back to the sequence diff")
        return {}, {}, {}
    recs = uv.parse(path)
    types, mech, spans = {}, {}, {}
    for acc, r in recs.items():
        if r["status"] == "canonical" or acc not in seqs:
            continue
        base = r["base"]
        if base not in seqs:
            continue
        c = uv.classify(r)
        types[acc] = c if isinstance(c, str) else uv.place(*c, len(seqs[base]))
        mech[acc] = r["mechanism"]
        spans[acc] = [(v, s - 1, e, k) for v, s, e, k, _n in r["vsps"]]
    print(f"  curated: {len(types):,} isoforms from {path}")
    return types, mech, spans


def pack_hits(hits, mapper=lambda j: j, order=None, group=False):
    """Peptides -> ([(start, end, digest, row)], n_rows), packed by overlap."""
    out, base_row = [], 0
    if not group:
        ends = []
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
    for dig in (order or ORDER):
        ends = []
        for pos, pep, d in hits:
            if d != dig:
                continue
            a = mapper(pos)
            b = mapper(pos + len(pep) - 1) + 1
            for t, last in enumerate(ends):
                if a >= last:
                    ends[t] = b
                    break
            else:
                t = len(ends)
                ends.append(b)
            out.append((a, b, dig, base_row + t))
        base_row += len(ends)
    return out, base_row


def evidence_panel(base, iso_acc, cseq, iseq, can_peps, iso_peps, out, font,
                   label=None, letter="", flank=9, max_rows=6, note=None,
                   letters=False, width=1215.0, flank_res=None,
                   browser=False, event=None, marks=None, marks_key=None,
                   region_label=True, dedupe=False, residues=True,
                   margin=None, bar_height=None, ruler=True,
                   gap_numbers=False):
    """Both forms as aligned bars, with the peptides that discriminate them."""
    colour = assign(ORDER + ["All"])
    ev = events(cseq, iseq)
    big = event or largest_event(cseq, iseq)
    if big is None:
        sys.exit(f"{base} and {iso_acc} have identical sequences")
    _, cs, ce, is_, ie = big
    p = cs
    kind = note or event_type(cseq, iseq)
    if len(ev) > 1 and not note and event is None:
        kind += f" (largest of {len(ev)} shown)"
    W = 1010.0
    ml, mr = (146.0 if margin is None else margin), 30.0
    barw = W - ml - mr
    scale = barw / max(len(cseq), len(iseq))
    bar_h = 14.0 if bar_height is None else bar_height
    bar_gap = 15.0
    mt = 52.0

    def picks(seq, peps, lo, hi):
        """Peptides overlapping [lo, hi) of `seq`, with position and digest."""
        found = []
        for dig in ORDER:
            for pep in sorted(peps.get(dig, ()), key=len):
                i = seq.find(pep)
                while i != -1:
                    if i < hi and i + len(pep) > lo:
                        found.append((i, pep, dig))
                        break
                    i = seq.find(pep, i + 1)
        found.sort(key=lambda t: (t[0], len(t[1])))
        return found

    can_hits = picks(cseq, can_peps, cs, max(ce, cs + 1))
    iso_hits = picks(iseq, iso_peps, max(is_ - 1, 0), max(ie, is_ + 1))
    rows = max(len(can_hits[:max_rows]), len(iso_hits[:max_rows]), 1)
    extra = 11 if (len(can_hits) > max_rows or len(iso_hits) > max_rows) else 0
    lists_y = mt + 2 * bar_h + bar_gap + 30
    if letters:
        def to_canon(j):
            """Isoform position -> the column it is drawn in."""
            if ie == is_:
                return j if j < cs else j + (ce - cs)
            return j

        ins_len, del_len = ie - is_, ce - cs
        if browser:
            def can_map(p):
                return p if p < ce else p + ins_len

            def iso_map(j):
                return j if j < is_ else j + del_len
        else:
            def can_map(p):
                return p

            iso_map = to_canon
        span = max(del_len + ins_len if browser else max(del_len, ins_len), 1)
        pad = flank_res if flank_res is not None else max(10, min(24, span))
        w0 = max(0, cs - pad)
        frame = (max(len(cseq) + ins_len, len(iseq) + del_len) if browser
                 else max(len(cseq), len(iseq)))
        w1 = min(frame, cs + span + pad)
        if dedupe:
            def unique(hits, mapper):
                spans = [(dig, max(mapper(pos), w0),
                          min(mapper(pos + len(pep) - 1) + 1, w1), (pos, pep, dig))
                         for pos, pep, dig in hits]
                spans.sort(key=lambda t: -(t[2] - t[1]))
                kept = []
                for dig, a, b, h in spans:
                    if not any(d == dig and a2 <= a and b <= b2
                               for d, a2, b2, _ in kept):
                        kept.append((dig, a, b, h))
                return sorted((h for *_, h in kept),
                              key=lambda t: (mapper(t[0]), len(t[1])))
            can_hits, iso_hits = unique(can_hits, can_map), unique(iso_hits, iso_map)
        can_rows, can_tiers = pack_hits(can_hits, can_map)
        iso_rows, iso_tiers = pack_hits(iso_hits, iso_map)
        band_c = bar_h + 14 + max(can_tiers, 1) * 6.0 + 12
        band_i = bar_h + 14 + max(iso_tiers, 1) * 6.0

        def room(sites):
            return 7.2 * max(map(len, sites.values())) + 20.0 if sites else 0.0
        h_can = room((marks or {}).get("canonical"))
        h_iso = room((marks or {}).get("isoform"))
        mt += h_can
        band_c += h_iso
        H = mt + band_c + band_i + 12 + (18 if marks_key else 0)
    else:
        H = lists_y + 14 + rows * 12 + extra + 22

    c = Canvas(W, H, font, font_scale=STRIP_TEXT, out_w=width)
    if letter:
        c.text(24, 26, letter, 13, INK, "start", "600")
    title = label or base
    tx = 146.0 - 88 + (24 if letter else 0)
    c.text(tx, 28, title, 12, INK, "start", "600")
    c.text(tx + len(title) * 7.6 * c.fs * TEXT_BOOST + 14, 28, kind, 9, INK_MUTED,
           "start")

    if letters:
        pitch = barw / max(w1 - w0, 1)

        def lollipops(y, sites, mapper, seq_):
            """Site markers above a band: one dot per digest, then the residue label."""
            placed = []
            for pos in sorted(sites):
                d = mapper(pos)
                if not (w0 <= d < w1):
                    continue
                x = ml + (d - w0 + 0.5) * pitch
                digs = sorted(sites[pos], key=ORDER.index)
                top = y - 6.0 - 7.2 * len(digs)
                c.line(x, y, x, top + 3.6, stroke=INK, sw=0.8)
                for k, dig in enumerate(digs):
                    c.add(f'<circle cx="{x:.1f}" cy="{y - 9.6 - 7.2 * k:.1f}" '
                          f'r="3.3" fill="{colour[dig]}" stroke="{SURFACE}" '
                          f'stroke-width="1"/>')
                lab = f"{seq_[pos]}{pos + 1}"
                half = 0.5 * 6.6 * 0.56 * len(lab)
                ly = top - 2.0
                while any(abs(x - px) < half + ph + 1.5 and abs(ly - py) < 7.0
                          for px, py, ph in placed):
                    ly -= 7.5
                placed.append((x, ly, half))
                c.text(x, ly, lab, 6.6, INK_SECONDARY, "middle")
        for name, y, is_iso in ((base, mt, False),
                                (iso_acc, mt + band_c, True)):
            seq_ = iseq if is_iso else cseq
            c.text(ml - 10, y + 5, name, 9.5, INK, "end", "600")
            sub_ = f"{len(seq_)} aa"
            if browser and del_len and ins_len:
                m_ = iso_map if is_iso else can_map
                rr = [r for r in range(len(seq_)) if w0 <= m_(r) < w1]
                if rr:
                    sub_ = f"{rr[0] + 1}–{rr[-1] + 1} of {len(seq_)}"
            c.text(ml - 10, y + 5 + 9.5 * c.fs * TEXT_BOOST * 0.95, sub_, 8,
                   INK_MUTED, "end")
            def cell(x, ch, tone):
                fill = SEQ_SAME if browser else \
                    {"same": SEQ_SAME, "can": SEQ_CAN, "iso": SEQ_ISO}[tone]
                ov = 0 if residues and bar_height is None else 1.0
                c.rect(x, y, pitch + ov, bar_h,
                       fill=fill, stroke="none", rx=0)
                if ch and residues:
                    fsz = min(pitch * (0.78 if bar_height is None else 1.05),
                              bar_h * 0.72)
                    by = (y + bar_h - 3.6 if bar_height is None
                          else y + bar_h / 2 + fsz * 0.36)
                    c.text(x + pitch / 2, by, ch, fsz / c.fs, ink_on(fill),
                           "middle", fit=True)

            if not is_iso:
                for pos in range(len(cseq)):
                    d = can_map(pos)
                    if not (w0 <= d < w1):
                        continue
                    cell(ml + (d - w0) * pitch, cseq[pos],
                         "can" if cs <= pos < ce else "same")
            else:
                for j in range(len(iseq)):
                    d = iso_map(j)
                    if not (w0 <= d < w1):
                        continue
                    cell(ml + (d - w0) * pitch, iseq[j],
                         "iso" if is_ <= j < ie and ie > is_ else "same")
                if region_label and (ce > cs or (browser and ie > is_)):
                    if ie == is_:
                        what = f"{ce - cs} aa absent"
                    elif ce == cs:
                        what = f"{ie - is_} aa inserted"
                    else:
                        what = f"{ce - cs} aa replaced by {ie - is_}"
                    lab = "discriminating region  ·  " + what
                    if browser:
                        e0 = max(cs if ce > cs else ce, w0)
                        e1 = min(ce if ce > cs else ce + ins_len, w1)
                        lx = ml + ((e0 + e1) / 2 - w0) * pitch
                        half = 0.5 * 7.6 * 0.56 * len(lab)
                        lx = min(max(lx, ml + half), ml + barw - half)
                    else:
                        lx = ml + barw / 2
                    c.text(lx, y - 5 - h_iso, lab, 7.6,
                           SEQ_SAME if browser else SEQ_ISO, "middle", "600")

            g0, g1 = (cs, ce) if is_iso else (ce, ce + ins_len)
            g0, g1 = max(g0, w0), min(g1, w1)
            if browser:
                if g1 > g0:
                    c.rect(ml + (g0 - w0) * pitch, y + bar_h / 2 - 1.0,
                           (g1 - g0) * pitch, 2.0, fill=SEQ_SAME,
                           stroke="none", rx=0)
                if gap_numbers and not is_iso and g1 > g0 and ce > 0:
                    c.text(ml + (g0 - w0) * pitch, y - 4, str(ce), 7.6,
                           INK_MUTED, "end")
                    if ce < len(cseq):
                        c.text(ml + (g1 - w0) * pitch, y - 4, str(ce + 1),
                               7.6, INK_MUTED, "start")

            sites = (marks or {}).get("isoform" if is_iso else "canonical")
            if sites:
                lollipops(y, sites, iso_map if is_iso else can_map, seq_)

            packed = iso_rows if is_iso else can_rows
            hits = iso_hits if is_iso else can_hits
            ty_ = y + bar_h + 3
            for a, b_, dig, t in packed:
                x0_ = ml + (max(a, w0) - w0) * pitch
                x1_ = ml + (min(b_, w1) - w0) * pitch
                if x1_ <= x0_:
                    continue
                gap_ = min(1.6, (x1_ - x0_) * 0.15)
                if browser and g1 > g0 and a < g0 < b_:
                    yb = ty_ + t * 6.0
                    for s0, s1 in ((a, g0), (g1, b_)):
                        u0 = ml + (max(s0, w0) - w0) * pitch
                        u1 = ml + (min(s1, w1) - w0) * pitch
                        if u1 > u0:
                            c.rect(u0, yb, u1 - u0, 4.4, fill=colour[dig],
                                   fo=0.9, rx=1.4)
                    k0 = ml + (max(g0, w0) - w0) * pitch
                    k1 = ml + (min(g1, w1) - w0) * pitch
                    if k1 > k0:
                        c.rect(k0, yb + 1.6, k1 - k0, 1.2, fill=colour[dig],
                               fo=0.55, stroke="none", rx=0)
                    continue
                c.rect(x0_, ty_ + t * 6.0, x1_ - x0_ - gap_, 4.4,
                       fill=colour[dig], fo=0.9, rx=1.4)

            if browser and not hits:
                c.text(ml, ty_ + 7.5,
                       "no peptide spans the junction in this form",
                       7.2, INK_MUTED, "start")

        if not ruler:
            pass
        elif browser and del_len and ins_len:
            pass
        else:
            c.text(ml, mt - h_can - 5, str(w0 + 1), 7.6, INK_MUTED, "start")
            c.text(ml + barw, mt - h_can - 5, str(w1), 7.6, INK_MUTED, "end")
            total = (len(iseq) if browser and ins_len and not del_len
                     else len(cseq))
            c.text(ml + barw / 2, mt - h_can - 5,
                   f"residues {w0 + 1}–{w1} of {total}", 7.8, INK_MUTED,
                   "middle")
    else:
        for name, seq, y in ((base, cseq, mt),
                             (iso_acc, iseq, mt + bar_h + bar_gap)):
            c.text(ml - 10, y + 7, name, 9.5, INK, "end", "600")
            c.text(ml - 10, y + 17, f"{len(seq)} aa", 8, INK_MUTED, "end")
            if name == base:
                c.rect(ml, y, len(cseq) * scale, bar_h, fill=STRIP_FILL,
                       stroke=AXIS, sw=1, rx=0)
                if ce > cs:
                    c.rect(ml + cs * scale, y, (ce - cs) * scale, bar_h,
                           fill="#8a8a8a", fo=0.9, rx=0)
                    c.text(ml + (cs + ce) / 2 * scale, y - 4, f"{cs + 1}–{ce}", 8,
                           INK_SECONDARY, "middle")
            else:
                for a, b in ((0, cs), (ce, len(cseq))):
                    if b > a:
                        c.rect(ml + a * scale, y, (b - a) * scale, bar_h,
                               fill=STRIP_FILL, stroke=AXIS, sw=1, rx=0)
                if ce > cs:
                    c.rect(ml + cs * scale, y, (ce - cs) * scale, bar_h,
                           fill=SURFACE, stroke=AXIS, sw=0.8, rx=0)
                    c.add(f'<path d="M {ml + cs * scale:.1f} {y + bar_h:.1f} '
                          f'L {ml + ce * scale:.1f} {y:.1f}" stroke="{AXIS}" '
                          f'stroke-width="0.8" fill="none"/>')
                    c.text(ml + (cs + ce) / 2 * scale, y + bar_h + 9,
                           f"{ce - cs} aa absent", 7.6, INK_MUTED, "middle")
                if ie > is_:
                    c.rect(ml + cs * scale, y, max((ie - is_) * scale, 2.0), bar_h,
                           fill=UNION, fo=0.9, rx=0)

    if letters:
        if marks_key:
            ky = H - 12
            c.line(ml + 3.3, ky + 3, ml + 3.3, ky - 3.6, stroke=INK, sw=0.8)
            c.add(f'<circle cx="{ml + 3.3:.1f}" cy="{ky - 6.9:.1f}" r="3.3" '
                  f'fill="{INK_MUTED}" stroke="{SURFACE}" stroke-width="1"/>')
            c.text(ml + 13, ky, marks_key, 7.6, INK_SECONDARY, "start")
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
        with open(out, "w") as fh:
            fh.write(c.out())
        print(f"  wrote {out}  ({kind}; {len(can_hits)} canonical-only, "
              f"{len(iso_hits)} isoform-only peptides)")
        return W, H

    ty = lists_y
    half = (W - ml + 88 - mr) / 2 - 16
    lx0 = ml - 104
    box_top = ty - 15
    box_bot = H - 26
    for col, (hits, seq, lo, hi, head, src_y) in enumerate((
            (can_hits, cseq, cs, ce, f"in {base} only", mt),
            (iso_hits, iseq, is_, ie, f"in {iso_acc} only",
             mt + bar_h + bar_gap))):
        x0 = lx0 + col * (half + 32)
        bx0, bx1 = x0 - 10, x0 + half + 4
        c.rect(bx0, box_top, bx1 - bx0, box_bot - box_top, fill=SURFACE,
               stroke=AXIS, sw=1.0, rx=6)
        if hi > lo:
            rl, rr = ml + lo * scale, ml + hi * scale
            ends = ((rl, bx0 + 12), (rr, min(bx1 - 12, bx0 + 12 + (rr - rl))))
        else:
            jx = ml + lo * scale
            ends = ((jx, bx0 + 12), (jx, bx0 + 60))
        for sx, tx_ in ends:
            c.add(f'<path d="M {sx:.1f} {src_y + bar_h:.1f} '
                  f'L {tx_:.1f} {box_top:.1f}" fill="none" '
                  f'stroke="{AXIS}" stroke-width="1" '
                  f'stroke-dasharray="3 3"/>')
        c.text(x0, ty, head, 9, INK, "start", "600")
        if not hits:
            c.text(x0, ty + 13, "no discriminating peptide", 8.2, INK_MUTED,
                   "start")
            continue
        c.text(x0 + len(head) * 6.1 + 10, ty,
               f"{len(hits)} peptide{'s' if len(hits) != 1 else ''}", 8.2,
               INK_MUTED, "start")
        shown = hits[:max_rows]
        wins = [(max(0, i - flank), min(len(seq), i + len(pep) + flank))
                for i, pep, _ in shown]
        left = max(lo - a0 for a0, _ in wins)
        right = max(a1 - lo for _, a1 in wins)
        pitch = min(7.2, (half - 26) / max(left + right, 1))
        ax = x0 + left * pitch
        for r, ((i, pep, dig), (a0, a1)) in enumerate(zip(shown, wins)):
            yy = ty + 14 + r * 12
            for j in range(a0, a1):
                x = ax + (j - lo) * pitch
                inpep = i <= j < i + len(pep)
                indiff = lo <= j < hi
                if inpep:
                    c.rect(x, yy - 7.5, pitch, 10, fill=colour[dig],
                           fo=0.9 if indiff else 0.32, rx=0)
                c.text(x + pitch / 2, yy, seq[j], 7.0,
                       "#ffffff" if (inpep and indiff) else INK, "middle")
            c.text(ax + (a1 - lo) * pitch + 7, yy, dig, 7.4, INK_MUTED, "start")
        c.line(ax, ty + 5, ax, ty + 10 + len(shown) * 12, stroke=INK, sw=0.8,
               so=0.35)
        if len(hits) > max_rows:
            c.text(x0, ty + 14 + max_rows * 12 + 1,
                   f"+{len(hits) - max_rows} more", 7.6, INK_MUTED, "start")

    c.text(lx0, H - 8,
           "Solid = inside the differing region; faded = shared flanks.", 7.8,
           INK_MUTED, "start")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {out}  ({kind}; {len(can_hits)} canonical-only, "
          f"{len(iso_hits)} isoform-only peptides)")
    return W, H


FORMAL = {
    "SE": "Skipped exon",
    "MXE": "Mutually exclusive exons",
    "A5SS": "Alternative 5' splice site",
    "A3SS": "Alternative 3' splice site",
    "RI": "Retained intron",
    "AF": "Alternative first coding exon",
    "AL": "Alternative last coding exon",
    "identical exon structure": "Identical CDS structure",
    "canonical and isoform share a transcript":
        "Cross-references conflict (shared transcript)",
    "VS_DEL": "Deletion",
    "VS_REP": "Replacement",
    "VS_MULTI": "Two or more features",
    "VS_NONE": "Sequence held in another entry",
}

CODE = {"SE", "MXE", "A5SS", "A3SS", "RI", "AF", "AL"}
ORDER_EV = ["SE", "MXE", "A5SS", "A3SS", "RI", "AF", "AL"]


def formal(k):
    """Table wording."""
    if k in CODE:
        return f"{FORMAL[k]} ({k})"
    return FORMAL.get(k, k[:1].upper() + k[1:])


def draw_pair(paths, seqs, protein, isoform, out, font, label=None, letter="",
              precursor_q=0.01, events_tsv="data/gencode/events.tsv", curation=None,
              mechanism=True, **panel):
    """One pair's evidence panel from the reports, under a line of outside claims."""
    _types, mech, spans = curation or ({}, {}, {})
    groups = collect_groups(paths, precursor_q, bases={protein})
    iso_ev, can_ev = split_evidence(groups)
    import prep_splice_events as sev
    bits = []
    codes = sev.load(events_tsv).get(isoform, ())
    typed = " + ".join(formal(k) for k in ORDER_EV if k in codes)
    if typed:
        bits.append(typed)
    vsps = spans.get(isoform, [])
    if len(vsps) == 1:
        bits.append(vsps[0][0])
    m = ", ".join(mech.get(isoform, []))
    if m and mechanism:
        bits.append(m.lower())
    note = " · ".join(bits) or None
    return evidence_panel(protein, isoform, seqs[protein], seqs[isoform],
                          can_ev.get(protein, {}), iso_ev.get(isoform, {}), out,
                          font, label, letter, note=note, **panel)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reports", nargs="+")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--protein", default=None, help="canonical accession")
    ap.add_argument("--isoform", default=None)
    ap.add_argument("--label", default=None)
    ap.add_argument("--events", default="data/gencode/events.tsv",
                    help="per-isoform splice events from prep_splice_events.py")
    ap.add_argument("--precursor-q", type=float, default=0.01)
    ap.add_argument("--letter", default="")
    ap.add_argument("--flank-res", type=int, default=None, metavar="N",
                    help="residues of context each side of the event")
    ap.add_argument("--browser", action="store_true",
                    help="genome-browser encoding (needs --letters)")
    ap.add_argument("--no-mechanism", dest="no_mechanism", action="store_true",
                    help="leave UniProt's mechanism off the header line")
    ap.add_argument("--no-ruler", dest="no_ruler", action="store_true",
                    help="drop the position ruler over the bars")
    ap.add_argument("--gap-numbers", dest="gap_numbers", action="store_true",
                    help="number the residues either side of the gap")
    ap.add_argument("--margin", type=float, default=None,
                    help="left margin before the bars (default 146)")
    ap.add_argument("--bar-height", dest="bar_height", type=float, default=None,
                    help="bar thickness (default 14)")
    ap.add_argument("--no-residues", dest="no_residues", action="store_true",
                    help="--letters layout without residue letters")
    ap.add_argument("--letters", action="store_true",
                    help="draw residue letters around the event")
    ap.add_argument("--out", default="isoform_event.svg")
    ap.add_argument("--font", default=FONT)
    ap.add_argument("--width", type=float, default=1215.0,
                    help="rendered width")
    ap.add_argument("--varseq", default="data/uniprot", metavar="DIR",
                    help="UniProt VAR_SEQ cache; '' forces the sequence diff")
    args = ap.parse_args(argv)

    paths = []
    for pat in args.reports:
        paths.extend(sorted(glob.glob(pat)) or [pat])
    from lib_fasta import read_fasta
    seqs = read_fasta(args.fasta)
    print(f"{len(seqs):,} sequences in {os.path.basename(args.fasta)}")

    types, mech, spans = curated(args.varseq, seqs)

    if not (args.protein and args.isoform):
        sys.exit("give --protein and --isoform")
    for a in (args.protein, args.isoform):
        if a not in seqs:
            sys.exit(f"{a} not in {args.fasta}")
    draw_pair(paths, seqs, args.protein, args.isoform, args.out, args.font,
              args.label, args.letter, args.precursor_q, args.events,
              (types, mech, spans), letters=args.letters, width=args.width,
              flank_res=args.flank_res, browser=args.browser,
              residues=not args.no_residues, margin=args.margin,
              bar_height=args.bar_height, ruler=not args.no_ruler,
              gap_numbers=args.gap_numbers, mechanism=not args.no_mechanism)
    return 0


if __name__ == "__main__":
    sys.exit(main())
