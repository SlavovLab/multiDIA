#!/usr/bin/env python3
"""Phosphosite x digest -log10 PEP heatmaps and per-protein modification maps.

    python3 fig3bc_phospho_atlas.py heatmap --scan derived/mods --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta \\
        --out figures/phospho/atlas_heatmap.svg

    python3 fig3bc_phospho_atlas.py heatmap --scan derived/mods --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta \\
        --genes SNCA TH@P07101-3 MAPT@P10636-8 ... \\
        --out figures/phospho/atlas_pd.svg

    python3 fig3bc_phospho_atlas.py protein 'data/search/*-60min-Phospho.parquet' \\
        --scan derived/mods --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta --isoform P10636-8 \\
        --region N1:45-73 --region N2:74-102 --region "Proline-rich:151-243" \\
        --uniprot-cache data/uniprot --rename "Tau/MAP =R" \\
        --out figures/phospho/atlas_mapt.svg
"""

import argparse
import collections
import glob
import gzip
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_svg import Canvas, esc                             # noqa: E402
from fig3a_phospho_sites import parse_mods                                   # noqa: E402
from lib_palette import (AXIS, DIVERGING_HIGH, DIVERGING_MID, FONT, INK,
                     INK_MUTED, INK_SECONDARY, PANEL_W, STRIP_FILL,
                     display, TEXT_BOOST)                                      # noqa: E402
import prep_phospho as ph                                               # noqa: E402

ORDER = ph.ORDER
NONE_FILL = "#000000"       # not detected
BREAKS = (3.0, 10.0, 20.0, 30.0)
STEPS = [DIVERGING_MID] + DIVERGING_HIGH
FIXED = "Carbamidomethyl"   # fixed modification, not an observation
TITLE_X = 36.0


def score(pep):
    """-log10 PEP, with PEP 0 put at the top bin."""
    return BREAKS[-1] if pep <= 0 else -math.log10(pep)


def fill(pep):
    if pep is None:
        return NONE_FILL
    s = score(pep)
    for i, b in enumerate(BREAKS):
        if s < b:
            return STEPS[i]
    return STEPS[-1]


def bin_labels():
    edges = [0.0] + list(BREAKS)
    out = [f"{edges[i]:g}–{edges[i + 1]:g}" for i in range(len(BREAKS))]
    return out + [f"≥ {BREAKS[-1]:g}"]


def cell(c, x, y, w, h, pep, tip):
    c.add(f'<g><title>{esc(tip)}</title>')
    c.rect(x, y, w, h, fill(pep))
    c.add('</g>')


KEY_ROW = 16.0 * TEXT_BOOST


def legend(c, x, y, size=10, max_width=None, gap=10.0):
    """Draw the PEP key; returns its right edge."""
    c.text(x, y - 5, "−log10 PEP (best precursor)", 9, INK_SECONDARY)
    entries = list(zip(STEPS, bin_labels())) + [(NONE_FILL, "not detected")]
    cx, cy, right = x, y, x
    for fill_, lab in entries:
        adv = size + 3 + 0.55 * 8.5 * TEXT_BOOST * len(lab) + gap
        if max_width and cx > x and cx + adv - gap > x + max_width:
            cx, cy = x, cy + KEY_ROW
        c.rect(cx, cy, size, size, fill_)
        c.text(cx + size + 3, cy + size - 1.5, lab, 8.5, INK_SECONDARY)
        right = max(right, cx + adv)
        cx += adv
    return right


KEY_Y = 32.0


def cell_key(c, x, y, size=14):
    """Draw the digest order of a box's cells; returns its right edge."""
    c.text(x, y - 5, "Cell order", 9, INK_SECONDARY)
    names = [display(d) for d in ORDER]
    w = 0.55 * 8.5 * TEXT_BOOST * max(map(len, names)) + 8
    for j, name in enumerate(names):
        c.rect(x + j * w, y, w - 1, size, "#ffffff", INK_MUTED, 0.8)
        c.text(x + j * w + (w - 1) / 2, y + size - 4, name, 8.5, INK_SECONDARY,
               "middle")
    return x + len(names) * w


def cohort(scan_dir):
    """Kept phospho precursor rows on the patients every digest measured."""
    patients = ph.complete_patients(ph.load_runs(scan_dir))
    keep = set(patients)
    return [r for r in ph.load(scan_dir) if r["sample"] in keep], patients


def site_cells(rows, S):
    """-> {site id: {digest: (best PEP, patients)}}."""
    out = {}
    for sid, s in S.items():
        best = {}
        for dig, sample, _, i in s["obs"]:
            p = float(rows[i]["pep"])
            b = best.setdefault(dig, [p, set()])
            b[0] = min(b[0], p)
            b[1].add(sample)
        out[sid] = {d: (p, len(n)) for d, (p, n) in best.items()}
    return out


def genes_of(fasta):
    gn = {}
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                acc = line.split("|")[1]
                m = re.search(r"\bGN=(\S+)", line)
                gn[acc] = m.group(1) if m else ""
    return gn


def base(acc):
    return ph.ISOFORM.sub("", acc)


SORT = ["Trypsin", "LysC", "GluC"]


def sort_key(cell):
    """Order by trypsin, then Lys-C, then Glu-C PEP; not detected last."""
    return tuple((0, cell[d][0]) if d in cell else (1, 0.0) for d in SORT)


def heatmap_all(rows, S, cells, out, font, key=True, title=True):
    order = sorted(S, key=lambda i: sort_key(cells[i]))
    first = [next(d for d in SORT if d in cells[i]) for i in order]
    W, H = PANEL_W, 236 if key else 192
    x0, x1 = 70.0, W - 14
    y_map, row_h = 62.0, 40.0
    c = Canvas(W, H, font)
    if title:
        ph.header(c, TITLE_X, 20, f"{len(order):,} phosphosites")
    n = len(order)
    dx = (x1 - x0) / n
    # blocks named by the first digest that detects a site
    names = {"Trypsin": "Trypsin", "LysC": "Lys-C",
             "GluC": "Glu-C"}
    k = 0
    while k < n:
        j = k
        while j + 1 < n and first[j + 1] == first[k]:
            j += 1
        xa, xb = x0 + k * dx, x0 + (j + 1) * dx
        c.line(xa + 0.5, y_map - 8, xb - 0.5, y_map - 8, INK_SECONDARY, 1.0)
        for xe in (xa + 0.5, xb - 0.5):
            c.line(xe, y_map - 8, xe, y_map - 4, INK_SECONDARY, 1.0)
        narrow = xb - xa < 120
        c.text(xb if narrow else (xa + xb) / 2, y_map - 12,
               f"{names[first[k]]} · {j - k + 1:,}", 9, INK_SECONDARY,
               "end" if narrow else "middle")
        k = j + 1
    for r, dig in enumerate(SORT):
        y = y_map + r * row_h
        c.text(x0 - 6, y + row_h / 2 + 3.5, display(dig), 10, INK, "end")
        # merge runs of one colour into one rect
        k = 0
        while k < n:
            f = fill(cells[order[k]].get(dig, (None,))[0])
            j = k
            while j + 1 < n and fill(cells[order[j + 1]]
                                     .get(dig, (None,))[0]) == f:
                j += 1
            c.rect(x0 + k * dx, y, (j - k + 1) * dx + 0.05, row_h, f)
            k = j + 1
    if key:
        legend(c, x0, y_map + 3 * row_h + 26)
    ph.save(c, os.path.dirname(out) or ".", os.path.basename(out))


def parse_items(specs):
    items = []
    for s in specs:
        s, _, acc = s.partition("@")
        label, _, genes = s.partition(":")
        items.append({"label": label,
                      "genes": [g for g in (genes or label).split(",") if g],
                      "numbering": acc or None})
    return items


def residue_in(entry_seq, rows, s, sid):
    """1-based residue of site `s` in another entry, or None."""
    for _, _, pep, i in s["obs"]:
        r = rows[i]
        for keys, _ in r["placements"]:
            if not keys & s["keys"]:
                continue
            o = next(o for o, (k, _) in zip(r["offsets"], r["placements"])
                     if k == keys)
            at = entry_seq.find(pep)
            if at >= 0:
                return at + o + 1
    return None


def heatmap_genes(rows, S, cells, items, seqs, gn, out, font,
                  key=True, width=PANEL_W, height=0.0):
    cols, empty = [], []
    for it in items:
        want = set(it["genes"])
        mine = []
        for sid, s in S.items():
            genes = {gn.get(base(k.split(":")[0]), "") for k in s["keys"]}
            if not genes & want:
                continue
            if any(genes & set(o["genes"]) for o in items[:items.index(it)]):
                continue
            hit = sorted(genes & want)
            name = (it["label"] if len(hit) != 1 or len(want) == 1
                    else hit[0])
            pos = None
            if it["numbering"]:
                pos = residue_in(seqs[it["numbering"]], rows, s, sid)
            if pos is None:
                k = sorted(k for k in s["keys"]
                           if gn.get(base(k.split(":")[0])) in want)[0]
                pos = int(k.split(":")[1])
            mine.append({"sid": sid, "name": name, "pos": pos,
                         "label": f"{name} {s['residue']}{pos}"})
        # merge site ids that land on the same residue
        merged = {}
        for m in mine:
            k = merged.setdefault(m["label"], dict(m, obs=[]))
            k["obs"] += S[m["sid"]]["obs"]
        mine = sorted(merged.values(), key=lambda m: m["pos"])
        for m in mine:
            m["cells"] = site_cells(rows, {0: m})[0]
        cols += mine
        if not mine:
            empty.append(it["label"])

    x0, row_h = 70.0, 18.0
    W = width
    cw = min(22.0, (W - x0 - 20) / max(len(cols), 1))
    y_map = 44.0
    y_lab = y_map + 3 * row_h + 6
    H = max(y_lab + 70 + (44 + KEY_ROW if key else 0), height)
    c = Canvas(W, H, font)
    ph.header(c, TITLE_X, 20, "Phosphosites on PD-implicated proteins")
    print(f"  {len(cols)} sites on {len(items) - len(empty)} of {len(items)} "
          f"items; no phosphosite on: {', '.join(empty) or 'none'}")
    for r, dig in enumerate(ORDER):
        y = y_map + r * row_h
        c.text(x0 - 6, y + row_h / 2 + 3.5, display(dig), 10, INK, "end")
        for j, m in enumerate(cols):
            pep, n = m["cells"].get(dig, (None, 0))
            tip = (f"{m['label']} · {display(dig)} · " +
                   (f"PEP {pep:.2g}, {n} patient{'s' * (n != 1)}"
                    if pep is not None else "not detected"))
            cell(c, x0 + j * cw + 1, y + 1, cw - 2, row_h - 2, pep, tip)
    prev = None
    for j, m in enumerate(cols):
        xc = x0 + j * cw + cw / 2
        c.text(xc + 3, y_lab, m["label"], 8.5, INK_SECONDARY, "end", rot=-60)
        if prev is not None and m["name"] != prev:
            c.line(x0 + j * cw, y_map - 4, x0 + j * cw, y_map + 3 * row_h + 2,
                   AXIS, 1.0)
        prev = m["name"]
    if key:
        import math
        lab_h = max(0.55 * 8.5 * TEXT_BOOST * len(m["label"]) for m in cols) \
            * math.sin(math.radians(60))
        right = legend(c, 20, y_lab + lab_h + 26, gap=9.0)
        assert right - 9.0 <= W, "PD heatmap key wider than the panel"
    ph.save(c, os.path.dirname(out) or ".", os.path.basename(out))


def heatmap(args):
    rows, patients = cohort(args.scan)
    print(f"  {len(patients)} patients measured by all three digests")
    S = ph.sites(rows)
    cells = site_cells(rows, S)
    sc = sorted(score(p) for v in cells.values() for p, _ in v.values())
    q = [sc[int(f * (len(sc) - 1))] for f in (0.1, 0.25, 0.5, 0.75, 0.9)]
    print(f"  {len(S):,} sites, {len(sc):,} site x digest cells; -log10 PEP "
          f"deciles 10/25/50/75/90: " + " / ".join(f"{v:.1f}" for v in q))
    if not args.genes:
        return heatmap_all(rows, S, cells, args.out, args.font,
                           key=not args.no_legend, title=not args.no_title)
    from lib_fasta import read_fasta
    heatmap_genes(rows, S, cells, parse_items(args.genes),
                  read_fasta(args.fasta), genes_of(args.fasta), args.out,
                  args.font, key=not args.no_legend, width=args.width,
                  height=args.height)


FEATURE = re.compile(r'(REPEAT|DOMAIN|REGION|ZN_FING|MOTIF)\s+(\d+)\.\.(\d+);'
                     r'\s*/note="([^"]*)"')


def uniprot_features(cache, acc, canon_seq, seq, rename):
    """UniProt features of `acc`, mapped onto `seq` by sequence."""
    from prep_uniprot import fetch
    path = fetch(cache, query=f"accession:{base(acc)}",
                 fields="accession,ft_repeat,ft_domain,ft_region,ft_zn_fing",
                 name=f"features_{base(acc)}.tsv.gz")
    with gzip.open(path, "rt") as fh:
        fh.readline()
        text = fh.read()
    out = []
    for kind, a, b, note in FEATURE.findall(text):
        if note == "Disordered":
            continue
        stretch = canon_seq[int(a) - 1:int(b)]
        at = seq.find(stretch)
        if at < 0:
            continue
        for old, new in rename:
            note = note.replace(old, new)
        out.append((note, at + 1, at + len(stretch), kind))
    # drop regions that only span drawn repeats/domains
    return [f for f in out if not any(
        g is not f and g[3] != "REGION" and f[1] <= g[1] and g[2] <= f[2]
        for g in out)]


def mod_label(mod, aa, pos):
    if mod.startswith("Phospho"):
        return f"p{aa}{pos}"
    if mod.startswith("Oxidation"):
        return f"{aa}{pos}ox"
    if mod.startswith("Acetyl"):
        return f"Ac-{aa}{pos}"
    return f"{mod.split(' ')[0]}-{aa}{pos}"


def protein_mods(paths, scan_dir, seq, precursor_q=0.01):
    """-> {(residue, aa, mod): {digest: [best PEP, patients]}} in `seq`."""
    import lib_report as rp
    kept = {(r["digest"], r["run"]) for r in ph.load_runs(scan_dir)
            if r["kept"] == "1"}
    keep = set(ph.complete_patients(ph.load_runs(scan_dir)))
    out = collections.defaultdict(dict)
    F = ["run", "modified", "peptide", "precursor_q", "pep_score"]
    for r in rp.open_reports(sorted(paths), order=ORDER):
        for b in r.batches(F):
            run, mod, pep, q, pe = (b[f] for f in F)
            for k in range(b["_n"]):
                if (q[k] is None or q[k] > precursor_q or not mod[k]
                        or pe[k] is None or "[" not in mod[k]
                        or (r.protease, run[k]) not in kept):
                    continue
                s = ph.sample_of(run[k] or "")
                if s not in keep:
                    continue
                at = seq.find(pep[k])
                if at < 0:
                    continue
                for o, m in parse_mods(mod[k]):
                    if m.startswith(FIXED):
                        continue
                    if "N-term" in m and at > 1:
                        continue      # protein N-term only (after Met removal)
                    key = (at + o + 1, pep[k][o], m)
                    cur = out[key].setdefault(r.protease, [1.0, set()])
                    cur[0] = min(cur[0], pe[k])
                    cur[1].add(s)
        print(f"  scanned {os.path.basename(r.path)}  ({r.protease})")
    return out


def layout(marks, px, lo, hi, box_w, gap=3.0, tiers=2):
    """Place label boxes above/below the bar -> [(mark, side, tier, box x)]."""
    ends = {(side, t): lo - gap for side in (-1, 1) for t in range(tiers)}
    placed = []
    for m in marks:
        want = min(max(px(m[0]) - box_w / 2, lo), hi - box_w)
        best = None
        for t in range(tiers):
            for side in (-1, 1):
                x = max(want, ends[(side, t)] + gap)
                cost = (abs(x - want), t)
                if best is None or cost < best[0]:
                    best = (cost, side, t, x)
        _, side, t, x = best
        ends[(side, t)] = x + box_w
        placed.append((m, side, t, x))
    # shift back anything pushed past the right edge
    for key, end in ends.items():
        over = end - hi
        if over > 0:
            row = [i for i, p in enumerate(placed) if (p[1], p[2]) == key]
            for i in reversed(row):
                m, side, t, x = placed[i]
                placed[i] = (m, side, t, x - over)
                nxt = x - over - gap
                prev = [placed[j][3] + box_w for j in row if j < i]
                over = (prev[-1] - nxt) if prev and prev[-1] > nxt else 0
                if over <= 0:
                    break
    return placed


def protein(args):
    from lib_fasta import read_fasta
    seqs = read_fasta(args.fasta)
    seq = seqs[args.isoform]
    gene = genes_of(args.fasta).get(args.isoform, args.isoform)
    paths = sorted(p for g in args.reports for p in glob.glob(g))
    mods = protein_mods(paths, args.scan, seq)
    feats = []
    if args.uniprot_cache:
        rename = [tuple(r.split("=", 1)) for r in args.rename]
        feats = uniprot_features(args.uniprot_cache, args.isoform,
                                 seqs[base(args.isoform)], seq, rename)
    for r in args.region:
        name, _, span = r.rpartition(":")
        a, b = span.split("-")
        feats.append((name, int(a), int(b), "REGION"))
    marks = sorted(((pos, aa, m) for pos, aa, m in mods),
                   key=lambda k: (k[0], k[2]))
    kinds = collections.Counter(m.split(" ")[0] for _, _, m in marks)
    print(f"  {gene} {args.isoform} ({len(seq)} aa): " +
          ", ".join(f"{n} {k}" for k, n in kinds.items()) +
          f"; {len(feats)} features")
    for pos, aa, m in marks:
        print(f"    {mod_label(m, aa, pos):>10}  " + "  ".join(
            f"{d}:{mods[(pos, aa, m)][d][0]:.1e}/{len(mods[(pos, aa, m)][d][1])}"
            if d in mods[(pos, aa, m)] else f"{d}:-" for d in ORDER))

    W = args.width
    lo, hi = 20.0, W - 20.0
    px = lambda r: lo + (r - 0.5) / len(seq) * (hi - lo)           # noqa: E731
    cs = args.cell
    box_w = 3 * cs + args.pad
    placed = layout(marks, px, lo, hi, box_w, tiers=1)
    reach = 44.0
    # alternate labels are staggered by one line
    line = args.label_size * TEXT_BOOST * 1.2
    y_bar = 48 + 14 + reach + line
    bar_h = 22.0
    H = max(y_bar + bar_h + reach + 14 + 56 + line + KEY_ROW,
            args.height)
    nth = {-1: 0, 1: 0}
    c = Canvas(W, H, font=args.font)
    ph.header(c, TITLE_X, 20, f"{gene} ({args.isoform}, {len(seq)} aa): modified "
              f"residues")
    c.rect(lo, y_bar, hi - lo, bar_h, "#ffffff", AXIS, 1.0)
    for name, a, b, _ in sorted(feats, key=lambda f: f[1]):
        c.rect(px(a - 0.5 + 0.5), y_bar, px(b + 0.5) - px(a), bar_h,
               STRIP_FILL, INK_SECONDARY, 0.8)
        mid = (px(a) + px(b)) / 2
        c.text(mid, y_bar + bar_h / 2 + 3.5, name, 9, INK, "middle")
    for r in [1] + list(range(100, len(seq), 100)) + [len(seq)]:
        c.line(px(r), y_bar + bar_h, px(r), y_bar + bar_h + 3, INK_MUTED, 0.8)
        c.text(px(r), y_bar + bar_h + 12, str(r), 7.5, INK_MUTED, "middle")
    for (pos, aa, m), side, _, x in placed:
        xs, xc = px(pos), x + box_w / 2
        stagger = line * (nth[side] % 2)
        nth[side] += 1
        if side < 0:
            ya, yb = y_bar, y_bar - reach
            yk, yl = yb + cs + 1, yb - 3 - stagger
        else:
            ya, yb = y_bar + bar_h, y_bar + bar_h + reach
            yk, yl = yb - 1, yb + cs + 9 + stagger
        bend = ya + 5 * side
        c.add(f'<path d="M {xs:.1f} {ya:.1f} L {xs:.1f} {bend:.1f} '
              f'L {xc:.1f} {yk - 5 * side:.1f} L {xc:.1f} {yk:.1f}" '
              f'fill="none" stroke="{INK_MUTED}" stroke-width="0.7"/>')
        c.text(xc, yl, mod_label(m, aa, pos), args.label_size, INK, "middle")
        cells = mods[(pos, aa, m)]
        for j, d in enumerate(ORDER):
            pep, n = (cells[d][0], len(cells[d][1])) if d in cells else (None, 0)
            tip = (f"{mod_label(m, aa, pos)} · {display(d)} · " +
                   (f"PEP {pep:.2g}, {n} patient{'s' * (n != 1)}"
                    if pep is not None else "not detected"))
            cell(c, x + args.pad / 2 + j * cs, yb, cs - 1, cs - 1, pep, tip)
    yk = H - KEY_Y - KEY_ROW
    kx = lo if args.no_legend else legend(c, lo, yk) + 10
    cell_key(c, kx, yk)
    ph.save(c, os.path.dirname(args.out) or ".", os.path.basename(args.out))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    h = sub.add_parser("heatmap", help="sites x digests, -log10 PEP cells")
    h.add_argument("--scan", required=True, help="prep_phospho.py scan's outdir")
    h.add_argument("--fasta", required=True)
    h.add_argument("--genes", nargs="*",
                   help="LABEL[:GENE,...][@ACCESSION] items, in display order")
    h.add_argument("--out", required=True)
    h.add_argument("--no-legend", action="store_true",
                   help="omit the PEP key")
    h.add_argument("--no-title", action="store_true",
                   help="omit the title")
    h.add_argument("--width", type=float, default=PANEL_W,
                   help="panel width for --genes")
    h.add_argument("--height", type=float, default=0.0,
                   help="minimum panel height, to match a row neighbour")
    h.add_argument("--font", default=FONT)
    h.set_defaults(fn=heatmap)
    p = sub.add_parser("protein", help="one protein's domains and mod cells")
    p.add_argument("reports", nargs="+", help="the phospho searches")
    p.add_argument("--scan", required=True, help="prep_phospho.py scan's outdir")
    p.add_argument("--fasta", required=True)
    p.add_argument("--isoform", required=True, help="entry to draw, e.g. P10636-8")
    p.add_argument("--region", action="append", default=[],
                   help="NAME:START-END in the drawn entry's coordinates")
    p.add_argument("--uniprot-cache", default=None,
                   help="add UniProt's repeats/domains, cached here")
    p.add_argument("--rename", action="append", default=[],
                   help="OLD=NEW substitution on UniProt feature names")
    p.add_argument("--out", required=True)
    p.add_argument("--width", type=float, default=PANEL_W, help="panel width")
    p.add_argument("--height", type=float, default=0.0,
                   help="minimum panel height, to match a row neighbour")
    p.add_argument("--cell", type=float, default=10.0, help="cell size")
    p.add_argument("--pad", type=float, default=12.0,
                   help="room around a box's three cells")
    p.add_argument("--label-size", type=float, default=8.5)
    p.add_argument("--no-legend", action="store_true",
                   help="omit the PEP key")
    p.add_argument("--font", default=FONT)
    p.set_defaults(fn=protein)
    args = ap.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
