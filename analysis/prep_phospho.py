#!/usr/bin/env python3
"""Phosphosites across the three digests.

    python3 prep_phospho.py scan 'data/search/*-Phospho.parquet' --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta \\
        --outdir derived/mods
    python3 prep_phospho.py known --cache data/uniprot
    python3 prep_phospho.py report --scan derived/mods --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta
    python3 prep_phospho.py figures --scan derived/mods --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta \\
        --outdir figures/phospho
"""

import argparse
import collections
import csv
import glob
import gzip
import math
import os
import re
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lib_svg import Canvas                                  # noqa: E402
from fig3a_phospho_sites import parse_mods                                   # noqa: E402
from lib_palette import (AXIS, FONT, INK, INK_MUTED, INK_SECONDARY, PANEL_W,
                     UNION, assign, display)                       # noqa: E402
from lib_report import sample_of                                       # noqa: E402
from lib_fasta import read_fasta                            # noqa: E402
from extra_why_multienzyme import Axes, nice, save                       # noqa: E402

ORDER = ["GluC", "LysC", "Trypsin"]
PHOSPHO = "Phospho"
STY = "STY"
ISOFORM = re.compile(r"-\d+$")
K = 7                      # shortest peptide the searches allow
STRATA = [(1, 1), (2, 3), (4, 7), (8, 10)]
MIN_STRATUM = 20           # points resting on fewer sites are not drawn


# ------------------------------------------------------------------- known ---
MOD_RES = re.compile(r'MOD_RES\s+(\d+);\s*/note="([^"]*)"'
                     r'(?:;\s*/evidence="([^"]*)")?')


def known_sites(cache_dir, force=False):
    """-> {canonical accession: {residue: evidence}} for UniProt's phospho MOD_RES."""
    from prep_uniprot import fetch
    path = fetch(cache_dir, force=force, fields="accession,ft_mod_res",
                 name="modres.tsv.gz")
    out = collections.defaultdict(dict)
    with gzip.open(path, "rt") as fh:
        fh.readline()
        for line in fh:
            acc, _, feats = line.rstrip("\n").partition("\t")
            for pos, note, ev in MOD_RES.findall(feats):
                if not note.startswith("Phospho"):
                    continue
                exp = any(c in (ev or "") for c in ("ECO:0000269", "ECO:0007744"))
                out[acc][int(pos)] = "experimental" if exp else "similarity"
    return out


# ------------------------------------------------------------------ locate ---
def locate_all(peptides, seqs):
    """-> {peptide: [(accession, 0-based start)]}, canonical entries preferred."""
    need = collections.defaultdict(list)
    for p in peptides:
        need[p[:K]].append(p)
    hits = collections.defaultdict(list)
    for acc, s in seqs.items():
        for i in range(len(s) - K + 1):
            ps = need.get(s[i:i + K])
            if ps:
                for p in ps:
                    if s.startswith(p, i):
                        hits[p].append((acc, i))
    out = {}
    for p in peptides:
        h = hits.get(p, [])
        canon = [x for x in h if not ISOFORM.search(x[0])]
        out[p] = sorted(canon or h)
    return out


def phospho_offsets(modified):
    """0-based offsets of phosphorylated residues within the stripped peptide."""
    return sorted(o for o, m in parse_mods(modified) if PHOSPHO in m)


# -------------------------------------------------------------------- scan ---
def scan(paths, fasta, outdir, precursor_q=0.01, min_run_precursors=1000):
    import lib_report as rp
    seqs = read_fasta(fasta)
    rows, depth = [], collections.defaultdict(lambda: [set(), set()])
    unresolved = collections.Counter()
    F = ["run", "modified", "peptide", "charge", "precursor_q", "pep_score",
         "quantity", "use_for_pep", "protein_groups"]
    for r in rp.open_reports(sorted(paths), order=ORDER):
        for b in r.batches(F):
            run, mod, pep, z, q, pe, qty, used, pg = (b[f] for f in F)
            for k in range(b["_n"]):
                if q[k] is None or q[k] > precursor_q or not mod[k]:
                    continue
                key = (r.protease, run[k])
                prec = hash((mod[k], z[k]))
                depth[key][0].add(prec)
                if PHOSPHO not in mod[k]:
                    continue
                depth[key][1].add(prec)
                s = sample_of(run[k] or "")
                if s is None:
                    unresolved[(r.protease, run[k])] += 1
                    continue
                rows.append([r.protease, run[k], s, mod[k], pep[k], z[k], q[k],
                             pe[k], qty[k], bool(used[k]), pg[k] or ""])
        print(f"  scanned {os.path.basename(r.path)}  ({r.protease})")
    for (dig, run), n in sorted(unresolved.items()):
        print(f"  WARNING: {dig} run {run!r} names no patient; "
              f"{n:,} phospho rows skipped")

    loc = locate_all({x[4] for x in rows}, seqs)
    os.makedirs(outdir, exist_ok=True)
    runs_path = os.path.join(outdir, "phospho_runs.tsv")
    with open(runs_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["digest", "run", "sample", "precursors",
                    "phospho_precursors", "kept"])
        for (dig, run), (allp, php) in sorted(depth.items()):
            w.writerow([dig, run, sample_of(run or "") or "", len(allp),
                        len(php), int(len(allp) >= min_run_precursors)])
    prec_path = os.path.join(outdir, "phospho_precursors.tsv")
    unplaced = 0
    with open(prec_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["digest", "run", "sample", "modified", "peptide", "charge",
                    "q", "pep", "quantity", "used_for_pep", "protein_groups",
                    "phospho_offsets", "locations"])
        for x in rows:
            where = loc.get(x[4], [])
            unplaced += not where
            w.writerow(x[:9] + [int(x[9]), x[10],
                                ",".join(map(str, phospho_offsets(x[3]))),
                                ";".join(f"{a}:{i}" for a, i in where)])
    print(f"  wrote {runs_path}  ({len(depth)} runs)")
    print(f"  wrote {prec_path}  ({len(rows):,} phospho precursors, "
          f"{len({x[4] for x in rows}):,} peptides, {unplaced:,} rows not "
          f"found in {os.path.basename(fasta)})")


# -------------------------------------------------------------------- load ---
def load_runs(scan_dir):
    with open(os.path.join(scan_dir, "phospho_runs.tsv"), newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load(scan_dir):
    """-> kept precursor rows, each with `placements`: [(site keys, residue)]."""
    kept = {(r["digest"], r["run"]) for r in load_runs(scan_dir)
            if r["kept"] == "1"}
    out = []
    path = os.path.join(scan_dir, "phospho_precursors.tsv")
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r["digest"], r["run"]) not in kept or not r["locations"]:
                continue
            where = [(a, int(i)) for a, i in
                     (x.rsplit(":", 1) for x in r["locations"].split(";"))]
            offs = [int(o) for o in r["phospho_offsets"].split(",") if o]
            r["offsets"] = offs
            r["where"] = where
            r["placements"] = [(frozenset(f"{a}:{i + o + 1}" for a, i in where),
                                r["peptide"][o]) for o in offs]
            r["q"] = float(r["q"])
            r["quantity"] = float(r["quantity"]) if r["quantity"] else None
            out.append(r)
    return out


def sites(rows):
    """Merge placements into sites by union-find over site keys."""
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for r in rows:
        for keys, _ in r["placements"]:
            ks = sorted(keys)
            for k in ks[1:]:
                a, b = find(ks[0]), find(k)
                if a != b:
                    parent[max(a, b)] = min(a, b)
    out = {}
    for i, r in enumerate(rows):
        for keys, aa in r["placements"]:
            sid = find(next(iter(keys)))
            s = out.setdefault(sid, {"keys": set(), "residue": aa, "obs": []})
            s["keys"] |= keys
            s["obs"].append((r["digest"], r["sample"], r["peptide"], i))
    return out


def complete_patients(runs):
    """Patients with a kept run in every digest, sorted, from `load_runs`."""
    have = collections.defaultdict(set)
    for r in runs:
        if r["kept"] == "1" and r["sample"]:
            have[r["sample"]].add(r["digest"])
    return sorted(s for s, d in have.items() if set(ORDER) <= d)


# ----------------------------------------------------------------- cohort ---
class Cohort:
    """Sites on the patients measured by every digest, with their evidence."""

    def __init__(self, rows, known, seqs, patients):
        self.patients = list(patients)
        keep = set(self.patients)
        self.rows = [r for r in rows if r["sample"] in keep]
        self.known = known
        self.seqs = seqs
        self.S = sites(self.rows)
        self.digests = {i: {o[0] for o in s["obs"]} for i, s in self.S.items()}
        self.seen_in = {i: {o[1] for o in s["obs"]} for i, s in self.S.items()}
        self.curated = {i: self._curated(s["keys"]) for i, s in self.S.items()}
        self.null = {i: self._null(i) for i in self.S}

    def _is_known(self, acc, pos):
        return pos in self.known.get(acc, {})

    def _curated(self, keys):
        return any(self._is_known(a, int(p))
                   for a, p in (k.rsplit(":", 1) for k in keys))

    def _peptides(self, sid):
        """Distinct (peptide, locations) carrying the site."""
        seen = {}
        for d, smp, pep, i in self.S[sid]["obs"]:
            seen.setdefault((pep, tuple(self.rows[i]["where"])), None)
        return list(seen)

    def _null(self, sid):
        vals = []
        for pep, where in self._peptides(sid):
            cand = [k for k, c in enumerate(pep) if c in STY]
            vals.append(sum(any(self._is_known(a, s + k + 1) for a, s in where)
                            for k in cand) / len(cand))
        return sum(vals) / len(vals)

    def random_residues(self, sid):
        """-> {S, T, Y: share} for a phosphate on a random S/T/Y of the site's peptides."""
        out = collections.Counter()
        peps = self._peptides(sid)
        for pep, _ in peps:
            cand = [c for c in pep if c in STY]
            for x in STY:
                out[x] += cand.count(x) / len(cand) / len(peps)
        return out

    def per_patient(self):
        """-> {patient: {digest: set(site ids)}}."""
        out = collections.defaultdict(lambda: collections.defaultdict(set))
        for sid, s in self.S.items():
            for d, smp, _, _ in s["obs"]:
                out[smp][d].add(sid)
        return out

    def greedy_order(self):
        """Digests in the order that adds most sites per patient (median)."""
        per = self.per_patient()
        chosen, left = [], list(ORDER)
        while left:
            def gain(d):
                return statistics.median(
                    len(set().union(*(per[p][x] for x in chosen + [d])))
                    for p in self.patients)
            best = max(left, key=gain)
            chosen.append(best)
            left.remove(best)
        return chosen

    def strata(self):
        """-> [(label, lo, hi, {"single"|"multi": (n, curated, null)})]."""
        out = []
        for lo, hi in STRATA:
            cell = {}
            for cls in ("single", "multi"):
                ids = [i for i in self.S
                       if lo <= len(self.seen_in[i]) <= hi
                       and (len(self.digests[i]) == 1) == (cls == "single")]
                n = len(ids)
                cell[cls] = (n, sum(self.curated[i] for i in ids) / n if n else None,
                             sum(self.null[i] for i in ids) / n if n else None)
            out.append((f"{lo}" if lo == hi else f"{lo}–{hi}", lo, hi, cell))
        return out

    def pairs(self, min_candidates=2):
        """Overlapping singly-phosphorylated peptides from two digests."""
        forms = {}
        for r in self.rows:
            if len(r["offsets"]) == 1:
                forms.setdefault((r["digest"], r["peptide"], r["offsets"][0]),
                                 r["where"])
        by_acc = collections.defaultdict(list)
        for (d, pep, off), where in forms.items():
            for acc, st in where:
                by_acc[acc].append((d, pep, off, st))
        seen, out = set(), []
        for acc, L in by_acc.items():
            seq = self.seqs[acc]
            for i in range(len(L)):
                for j in range(i + 1, len(L)):
                    a, b = L[i], L[j]
                    if a[0] == b[0]:
                        continue
                    lo = max(a[3], b[3])
                    hi = min(a[3] + len(a[1]), b[3] + len(b[1]))
                    sa, sb = a[3] + a[2], b[3] + b[2]
                    if not (lo <= sa < hi and lo <= sb < hi):
                        continue
                    key = tuple(sorted([a[:3], b[:3]]))
                    if key in seen:
                        continue
                    seen.add(key)
                    m = sum(1 for x in seq[lo:hi] if x in STY)
                    if m >= min_candidates:
                        out.append((sa == sb, m, acc))
        return out

    def region(self, acc, residue):
        """Singly-phosphorylated peptides placing the phosphate on `residue`."""
        found = collections.defaultdict(lambda: [set(), set()])
        for r in self.rows:
            if len(r["offsets"]) != 1:
                continue
            for a, st in r["where"]:
                if a == acc and st + r["offsets"][0] + 1 == residue:
                    found[(st, r["peptide"])][0].add(r["digest"])
                    found[(st, r["peptide"])][1].add(r["sample"])
        return sorted(((st, pep, [d for d in ORDER if d in ds],
                        sum(c in STY for c in pep), len(pts))
                       for (st, pep), (ds, pts) in found.items()),
                      key=lambda x: (x[0], -len(x[1])))


# ------------------------------------------------------------------ panels ---
def header(c, x, y, title, sub=None):
    c.text(x, y, title, 11.5, INK, "start", "600")
    if sub:
        c.text(x, y + 15, sub, 9, INK_MUTED, "start")


def ylabel(c, x, y, s):
    c.add(f'<g transform="translate({x:.1f} {y:.1f}) rotate(-90)">'
          f'<text x="0" y="0" font-size="10" text-anchor="middle" fill="{INK}">'
          f'{s}</text></g>')


def dot(c, x, y, r, fill, stroke="#ffffff", sw=1.0):
    c.add(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}" '
          f'stroke="{stroke}" stroke-width="{sw}"/>')


def polyline(c, pts, stroke, sw=1.5, so=None):
    d = "M " + " L ".join(f"{a:.1f} {b:.1f}" for a, b in pts)
    extra = f' stroke-opacity="{so}"' if so is not None else ""
    c.add(f'<path d="{d}" fill="none" stroke="{stroke}" stroke-width="{sw}" '
          f'stroke-linejoin="round"{extra}/>')


def dashed(c, pts, stroke, sw=1.4, on=4.0, off=3.0):
    """A dashed polyline drawn as segments; MuPDF ignores stroke-dasharray."""
    carry = 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        L = math.hypot(x1 - x0, y1 - y0)
        t = 0.0
        while t < L:
            phase = (carry + t) % (on + off)
            step = min((on if phase < on else on + off) - phase, L - t)
            if phase < on:
                a, b = t / L, (t + step) / L
                c.line(x0 + (x1 - x0) * a, y0 + (y1 - y0) * a,
                       x0 + (x1 - x0) * b, y0 + (y1 - y0) * b, stroke=stroke, sw=sw)
            t += step
        carry = (carry + L) % (on + off)


def panel_gain(co, shallow, outdir, font, letter=""):
    """Phosphosites per patient as each digest is added, greedy order."""
    order = co.greedy_order()
    colour = assign(ORDER)
    per = co.per_patient()
    steps = {p: [len(set().union(*(per[p][d] for d in order[:k + 1])))
                 for k in range(len(order))] for p in co.patients}
    med = [statistics.median(v[k] for v in steps.values())
           for k in range(len(order))]
    gain = [statistics.median(v[k] / v[0] - 1 for v in steps.values())
            for k in range(len(order))]
    W, H = PANEL_W, 310
    c = Canvas(W, H, font)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    top = max(max(v) for v in steps.values()) * 1.08
    ax = Axes(c, 86, 96, W - 86 - 190, H - 96 - 44, (-0.4, len(order) - 0.6),
              (0, top))
    header(c, ax.x, 32, "Each added protease adds phosphosites in every patient",
           f"n = {len(co.patients)} patients measured by all three digests · "
           f"median per patient above each step")
    ax.ygrid(nice(top, 4))
    c.line(ax.x, ax.y + ax.h, ax.x + ax.w, ax.y + ax.h, stroke=AXIS, sw=1)
    for p, v in steps.items():
        polyline(c, [(ax.X(k), ax.Y(n)) for k, n in enumerate(v)], INK_MUTED,
                 sw=1.0, so=0.55)
    for k, d in enumerate(order):
        for p, v in steps.items():
            dot(c, ax.X(k), ax.Y(v[k]), 3.4, colour[d])
        c.line(ax.X(k) - 22, ax.Y(med[k]), ax.X(k) - 8, ax.Y(med[k]),
               stroke=INK, sw=2)
        c.text(ax.X(k), ax.y - 16, f"{med[k]:,.0f}", 10, INK, "middle", "600")
        if k:
            c.text(ax.X(k), ax.y - 3, f"+{gain[k]:.0%}", 8.6, INK_MUTED, "middle")
        lab = display(d) if k == 0 else f"+ {display(d)}"
        c.text(ax.X(k), ax.y + ax.h + 18, lab, 9.5, INK, "middle")
    for p in shallow:
        if p in steps:
            v = steps[p]
            c.text(ax.X(len(order) - 1) + 10, ax.Y(v[-1]) + 3.5,
                   f"{p} · shallow {display(shallow[p])} run", 8.6,
                   INK_SECONDARY, "start")
    ylabel(c, 30, ax.y + ax.h / 2, "Phosphosites per patient")
    save(c, outdir, "phospho_gain.svg")
    return order, med, gain


def panel_curated(co, outdir, font, letter=""):
    """Curated share by replication, one digest against two or more."""
    st = co.strata()
    W, H = PANEL_W, 320
    c = Canvas(W, H, font)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    ax = Axes(c, 86, 78, W - 86 - 196, H - 78 - 86, (-0.3, len(st) - 0.85), (0, 1))
    header(c, ax.x, 32,
           "More patients barely move a phosphosite call from chance; "
           "a second digest does",
           "share curated in UniProt · dashed: the same peptides with the "
           "phosphate on a random S/T/Y")
    ax.ygrid([0, .25, .5, .75, 1.0], lambda v: f"{v:.0%}")
    c.line(ax.x, ax.y + ax.h, ax.x + ax.w, ax.y + ax.h, stroke=AXIS, sw=1)
    series = [("multi", UNION, "two or three digests"),
              ("single", INK_SECONDARY, "one digest")]
    ends = []
    for cls, col, name in series:
        pts = [(k, cell[cls]) for k, (_, _, _, cell) in enumerate(st)
               if cell[cls][0] >= MIN_STRATUM]
        if not pts:
            continue
        obs = [(ax.X(k), ax.Y(v[1])) for k, v in pts]
        nul = [(ax.X(k), ax.Y(v[2])) for k, v in pts]
        dashed(c, nul, col)
        polyline(c, obs, col, sw=2.2)
        for x, y in nul:
            dot(c, x, y, 3.0, "#ffffff", col, 1.3)
        for x, y in obs:
            dot(c, x, y, 3.8, col)
        ends.append((obs[-1][1], f"{name} {pts[-1][1][1]:.0%}", col, False))
        ends.append((nul[-1][1], f"random placement {pts[-1][1][2]:.0%}", col,
                     True))
    # direct labels, nudged apart where they would touch
    ends.sort()
    placed = []
    for y, name, col, is_null in ends:
        if placed and y - placed[-1] < 13:
            y = placed[-1] + 13
        placed.append(y)
        line_key(c, ax.x + ax.w + 10, y, col, dashed_key=is_null)
        c.text(ax.x + ax.w + 31, y + 3.5, name, 9,
               INK_SECONDARY if is_null else INK, "start",
               None if is_null else "600")
    for k, (lab, lo, hi, cell) in enumerate(st):
        c.text(ax.X(k), ax.y + ax.h + 18, lab, 9.5, INK, "middle")
        for row, (cls, col, _) in enumerate(series):
            n = cell[cls][0]
            c.text(ax.X(k), ax.y + ax.h + 34 + 12 * row,
                   f"{n:,}" if n >= MIN_STRATUM else f"({n})", 8.2,
                   INK_SECONDARY if n >= MIN_STRATUM else INK_MUTED, "middle")
    for row, (cls, col, _) in enumerate(series):
        line_key(c, ax.x - 26, ax.y + ax.h + 31 + 12 * row, col)
    c.text(ax.x - 30, ax.y + ax.h + 34, "sites", 8.2, INK_MUTED, "end")
    c.text(ax.x + ax.w / 2, H - 8,
           f"Patients in which the site was placed (of {len(co.patients)})",
           10, INK, "middle")
    ylabel(c, 30, ax.y + ax.h / 2, "Sites curated in UniProt")
    save(c, outdir, "phospho_curated.svg")
    return st


def panel_agreement(co, outdir, font, letter=""):
    """Agreement on the residue between two digests, against chance."""
    pairs = co.pairs()
    if not pairs:
        print("  agreement: no overlapping pairs with two or more candidates")
        return [], None, None
    bins = [(2, 2), (3, 3), (4, 5), (6, 99)]
    rows = []
    for lo, hi in bins:
        sel = [p for p in pairs if lo <= p[1] <= hi]
        if not sel:
            continue
        n = len(sel)
        rows.append((f"{lo}" if lo == hi else (f"{lo}–{hi}" if hi < 99
                                               else f"≥{lo}"),
                     n, sum(p[0] for p in sel) / n, sum(1 / p[1] for p in sel) / n))
    W, H = PANEL_W, 290
    c = Canvas(W, H, font)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    ax = Axes(c, 86, 78, W - 86 - 60, H - 78 - 70, (-0.5, len(rows) - 0.5), (0, 1))
    n_all = len(pairs)
    agree = sum(p[0] for p in pairs) / n_all
    chance = sum(1 / p[1] for p in pairs) / n_all
    header(c, ax.x, 32,
           "Two digests put the phosphate on the same residue more often than chance",
           f"{n_all:,} overlapping peptide pairs · agree {agree:.0%}, "
           f"chance {chance:.0%}")
    ax.ygrid([0, .25, .5, .75, 1.0], lambda v: f"{v:.0%}")
    c.line(ax.x, ax.y + ax.h, ax.x + ax.w, ax.y + ax.h, stroke=AXIS, sw=1)
    bw = ax.w / len(rows) * 0.42
    for k, (lab, n, a, ch) in enumerate(rows):
        x = ax.X(k) - bw / 2
        c.rect(x, ax.Y(a), bw, ax.y + ax.h - ax.Y(a), fill=UNION, fo=0.85, rx=2)
        c.text(ax.X(k), ax.Y(a) - 7, f"{a:.0%}", 9.5, INK, "middle", "600")
        c.line(x - 8, ax.Y(ch), x + bw + 8, ax.Y(ch), stroke=INK, sw=1.6)
        c.text(x + bw + 12, ax.Y(ch) + 3.5, f"chance {ch:.0%}", 8.4,
               INK_SECONDARY, "start")
        c.text(ax.X(k), ax.y + ax.h + 18, lab, 9.5, INK, "middle")
        c.text(ax.X(k), ax.y + ax.h + 32, f"{n:,} pairs", 8.2, INK_MUTED, "middle")
    c.text(ax.x + ax.w / 2, H - 8, "S/T/Y candidates in the shared stretch", 10,
           INK, "middle")
    ylabel(c, 30, ax.y + ax.h / 2, "Pairs placing the same residue")
    save(c, outdir, "phospho_agreement.svg")
    return rows, agree, chance


def residue_mix(co):
    """-> [(digests, observed Counter, random Counter, n)] for 1, 2, 3 digests."""
    out = []
    for k in (1, 2, 3):
        ids = [i for i in co.S if len(co.digests[i]) == k]
        obs = collections.Counter(co.S[i]["residue"] for i in ids)
        rnd = collections.Counter()
        for i in ids:
            rnd.update(co.random_residues(i))
        out.append((k, obs, rnd, len(ids)))
    return out


def panel_residues(co, outdir, font, letter=""):
    """pS / pT / pY share, observed against random placement, by digests."""
    mix = [m for m in residue_mix(co) if m[3]]
    shade = {"S": INK_SECONDARY, "T": INK_MUTED, "Y": AXIS}
    W = PANEL_W
    bar, gap, ggap = 17.0, 5.0, 20.0
    x0, x1 = 214.0, W - 62.0
    top = 84.0
    H = round(top + len(mix) * (2 * bar + gap + ggap) + 10)
    c = Canvas(W, H, font)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    header(c, 48, 32,
           "One-digest calls have the residue mix of random placement; "
           "agreed calls do not",
           "pS / pT / pY share of sites · random: the phosphate on a random "
           "S/T/Y of the same peptides")
    names = {1: "one digest", 2: "two digests", 3: "three digests"}
    y = top
    for k, obs, rnd, n in mix:
        c.text(48, y + bar + gap / 2 + 4, names[k], 9.5, INK, "start", "600")
        c.text(48, y + bar + gap / 2 + 17, f"{n:,} sites", 8.4, INK_MUTED, "start")
        for j, (lab, cnt) in enumerate((("observed", obs), ("random", rnd))):
            by = y + j * (bar + gap)
            tot = sum(cnt.values())
            c.text(x0 - 10, by + bar / 2 + 3.5, lab, 8.6, INK_SECONDARY, "end")
            cx = x0
            for r in STY:
                w = (x1 - x0) * cnt[r] / tot
                c.rect(cx, by, w, bar, fill=shade[r])
                if w > 34:
                    ink = "#ffffff" if r != "Y" else INK
                    c.text(cx + w / 2, by + bar / 2 + 3.5,
                           f"p{r} {cnt[r] / tot:.0%}", 8.4, ink, "middle")
                elif r == "Y":
                    c.text(x1 + 6, by + bar / 2 + 3.5,
                           f"pY {cnt[r] / tot:.0%}", 8.4, INK_SECONDARY, "start")
                cx += w
        y += 2 * bar + gap + ggap
    save(c, outdir, "phospho_residues.svg")
    return mix


def panel_example(co, acc, residue, label, outdir, font, letter="", flank=3):
    """Every peptide placing one phosphate on `residue`, candidates marked."""
    peps = co.region(acc, residue)
    if not peps:
        sys.exit(f"{acc}: no singly-phosphorylated peptide places a phosphate "
                 f"on residue {residue}")
    seq = co.seqs[acc]
    colour = assign(ORDER)
    lo = max(0, min(p[0] for p in peps) - flank)
    hi = min(len(seq), max(p[0] + len(p[1]) for p in peps) + flank)
    ml, mr = 40.0, 196.0
    W = PANEL_W
    pitch = (W - ml - mr) / (hi - lo)
    bar_h, step = 8.0, 17.0
    top = 104.0
    H = round(top + 26 + len(peps) * step + 16)
    c = Canvas(W, H, font)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    alone = sum(1 for p in peps if p[3] == 1)
    digs = {d for p in peps for d in p[2]}
    aa = seq[residue - 1]
    header(c, 48, 32,
           f"{label} {aa}{residue}: {len(peps)} peptides from {len(digs)} "
           f"digests, one placement",
           f"{aa}{residue} is the only S/T/Y in {alone} of them · each "
           f"peptide placed in {min(p[4] for p in peps)}–"
           f"{max(p[4] for p in peps)} of {len(co.patients)} patients")

    def X(i):
        return ml + (i - lo) * pitch

    ruler = top - 30
    c.line(X(lo), ruler, X(hi), ruler, stroke=AXIS, sw=1)
    for pos in range(lo + 1, hi + 1):
        if pos % 10 == 0 or pos in (lo + 1, hi):
            if pos in (lo + 1, hi) and any(abs(pos - q) <= 2 and q % 10 == 0
                                           for q in range(lo + 1, hi + 1)):
                continue
            c.line(X(pos - 0.5), ruler, X(pos - 0.5), ruler - 4, stroke=AXIS, sw=1)
            c.text(X(pos - 0.5), ruler - 7, str(pos), 7.6, INK_MUTED, "middle")
    kx = W - mr + 14
    dot(c, kx + 4, ruler - 10, 4.4, INK, "#ffffff", 1.2)
    c.text(kx + 14, ruler - 6.5, "placed phosphate", 8.4, INK_SECONDARY, "start")
    dot(c, kx + 4, ruler + 4, 3.0, "#ffffff", INK, 1.1)
    c.text(kx + 14, ruler + 7.5, "other S/T/Y", 8.4, INK_SECONDARY, "start")
    rows_top = top + 14
    c.rect(X(residue - 1), top - 12, pitch,
           rows_top + len(peps) * step - (top - 12), fill=UNION, fo=0.10, rx=2)
    for i in range(lo, hi):
        c.text(X(i + 0.5), top, seq[i], 8.6,
               INK if i == residue - 1 else INK_SECONDARY, "middle",
               "600" if i == residue - 1 else None)
    for k, (st, pep, ds, cand, npt) in enumerate(peps):
        y = rows_top + k * step
        h = bar_h / len(ds)
        for j, d in enumerate(ds):
            c.rect(X(st), y + j * h, len(pep) * pitch, h, fill=colour[d], fo=0.9)
        for off, ch in enumerate(pep):
            if ch not in STY:
                continue
            x = X(st + off + 0.5)
            if st + off + 1 == residue:
                dot(c, x, y + bar_h / 2, 4.4, INK, "#ffffff", 1.2)
            else:
                dot(c, x, y + bar_h / 2, 3.0, "#ffffff", INK, 1.1)
        name = " + ".join(display(d) for d in ds)
        c.text(W - mr + 14, y + bar_h / 2 + 3.5,
               f"{name} · " + ("only S/T/Y" if cand == 1
                                    else f"1 of {cand} S/T/Y"),
               8.6, INK_SECONDARY, "start")
    save(c, outdir, f"phospho_example_{label.lower()}_{aa.lower()}{residue}.svg")
    return peps


# ----------------------------------------------------------- differential ---
def site_quantities(rows, S):
    """-> ({digest: {site: {patient: quantity}}}, {digest: {site: peptide}})."""
    site_of = {}
    for sid, s in S.items():
        for _, _, _, i in s["obs"]:
            site_of[i] = sid
    units = collections.defaultdict(
        lambda: collections.defaultdict(lambda: collections.defaultdict(float)))
    for i, r in enumerate(rows):
        if (len(r["offsets"]) == 1 and r["used_for_pep"] == "1"
                and r["quantity"] and r["quantity"] > 0):
            units[r["digest"]][(site_of[i], r["peptide"])][r["sample"]] += r["quantity"]
    qty, ref = {}, {}
    for dig, by_unit in units.items():
        best = {}
        for (sid, pep), q in by_unit.items():
            rank = (len(q), statistics.median(q.values()))
            if sid not in best or rank > best[sid][0]:
                best[sid] = (rank, pep)
        qty[dig] = {sid: dict(by_unit[(sid, pep)]) for sid, (_, pep) in best.items()}
        ref[dig] = {sid: pep for sid, (_, pep) in best.items()}
    return qty, ref


def protein_frames(quant_dir, min_run_frac=0.5, max_mad=3.0):
    """Per digest, the protein matrix exactly as `extra_differential.py test` sees it."""
    import extra_differential as df
    out = {}
    for dig in ORDER:
        path = os.path.join(quant_dir, f"{dig}.csv")
        if not os.path.exists(path):
            continue
        mat, samples = df.load_matrix(path)
        keep, dropped, _ = df.qc_runs(mat, samples, min_run_frac)
        logs, run_med, grand = df.normalise(mat, keep)
        bad = df.outlier_runs(logs, keep, max_mad)
        if bad:
            keep = [s for s in keep if s not in bad]
            logs, run_med, grand = df.normalise(mat, keep)
        out[dig] = {"keep": keep, "dropped": sorted(set(samples) - set(keep)),
                    "logs": logs, "mat": mat,
                    "offset": {s: run_med[s] - grand for s in keep}}
    return out


def protein_of(rows, S, ref, frames):
    """-> {digest: {site: protein group}} in the unmodified search's matrix."""
    site_of = {}
    for sid, s in S.items():
        for _, _, _, i in s["obs"]:
            site_of[i] = sid
    seen = collections.defaultdict(collections.Counter)
    for i, r in enumerate(rows):
        sid = site_of.get(i)
        if sid is not None and ref.get(r["digest"], {}).get(sid) == r["peptide"]:
            seen[(r["digest"], sid)][r["protein_groups"]] += 1
    out = collections.defaultdict(dict)
    for dig, fr in frames.items():
        by_acc = collections.defaultdict(set)
        for g in fr["mat"]:
            for a in g.split(";"):
                by_acc[a].add(g)
        for sid in ref.get(dig, {}):
            groups = [g for g, _ in seen[(dig, sid)].most_common()]
            hit = next((g for g in groups if g in fr["mat"]), None)
            if hit is None:
                cand = {g for k in S[sid]["keys"]
                        for g in by_acc.get(k.rsplit(":", 1)[0], ())}
                if cand:
                    hit = max(cand, key=lambda g: (len(fr["mat"][g]), g))
            if hit is not None:
                out[dig][sid] = hit
    return out


class Differential:
    """Phosphosite LBD-vs-control tests, per digest and across digests."""

    def __init__(self, rows, S, frames, cond, case="LBD", control="Control",
                 min_per_group=3):
        self.S, self.cond = S, cond
        self.case, self.control = case, control
        self.min_per_group = min_per_group
        qty, self.ref = site_quantities(rows, S)
        self.prot = protein_of(rows, S, self.ref, frames)
        self.frames = frames
        self.level, self.occupancy = {}, {}
        for dig, fr in frames.items():
            lv, oc = {}, {}
            for sid, q in qty.get(dig, {}).items():
                v = {s: math.log2(x) - fr["offset"][s]
                     for s, x in q.items() if s in fr["offset"]}
                if not v:
                    continue
                lv[sid] = v
                g = self.prot.get(dig, {}).get(sid)
                if g is not None:
                    o = {s: x - fr["logs"][g][s] for s, x in v.items()
                         if s in fr["logs"][g]}
                    if o:
                        oc[sid] = o
            self.level[dig], self.occupancy[dig] = lv, oc

    def groups(self, labels, samples):
        return ([s for s in samples if labels.get(s) == self.case],
                [s for s in samples if labels.get(s) == self.control])

    def run(self, labels=None):
        """-> {(quantity, digest or "Combined"): results} under `labels`."""
        import extra_differential as df
        labels = labels or self.cond
        out = {}
        for name, per in (("level", self.level), ("occupancy", self.occupancy)):
            for dig, logs in per.items():
                out[(name, dig)] = df.test_digest_moderated(
                    logs, self.groups(labels, self.frames[dig]["keep"]),
                    self.min_per_group)
            comb, _ = df.combine_patients(per, 1)
            out[(name, "Combined")] = df.test_digest_moderated(
                comb, self.groups(labels, sorted(labels)), self.min_per_group)
        return out

    def permutations(self, n=200, seed=20260923):
        """Balanced label shuffles, the same shuffle applied to every digest."""
        import random
        rng = random.Random(seed)
        pats = sorted(self.cond)
        labs = [self.cond[p] for p in pats]
        seen, out = set(), []
        while len(out) < n:
            rng.shuffle(labs)
            key = tuple(labs)
            if key in seen or key == tuple(self.cond[p] for p in pats):
                continue
            seen.add(key)
            out.append(dict(zip(pats, labs)))
        return out


    def protein_tests(self, labels=None):
        """-> {digest: {group: log2FC}}, the protein test under `labels`."""
        import extra_differential as df
        labels = labels or self.cond
        return {dig: {r["group"]: r["log2fc"] for r in df.test_digest_moderated(
                    fr["logs"], self.groups(labels, fr["keep"]), self.min_per_group)}
                for dig, fr in self.frames.items()}


PAIRS = [("GluC", "LysC"), ("GluC", "Trypsin"), ("LysC", "Trypsin")]


def summarise_tests(res, q_cut=0.10):
    """-> {key: (tested, share p < 0.05, n at q < q_cut, best q)}."""
    out = {}
    for k, rs in res.items():
        if not rs:
            continue
        out[k] = (len(rs), sum(r["p"] < 0.05 for r in rs) / len(rs),
                  sum(r["q"] < q_cut for r in rs), min(r["q"] for r in rs))
    return out


def agreement(res, prot):
    """-> {(quantity, a, b): (shared, Pearson r)} between two digests' fold changes."""
    import numpy as np
    fcs = {("protein", d): v for d, v in prot.items()}
    for (name, dig), rs in res.items():
        if dig in ORDER:
            fcs[(name, dig)] = {r["group"]: r["log2fc"] for r in rs}
    out = {}
    for name in ("protein", "level", "occupancy"):
        for a, b in PAIRS:
            fa, fb = fcs.get((name, a), {}), fcs.get((name, b), {})
            shared = sorted(set(fa) & set(fb))
            if len(shared) >= 3:
                x = np.array([fa[g] for g in shared])
                y = np.array([fb[g] for g in shared])
                out[(name, a, b)] = (len(shared), float(np.corrcoef(x, y)[0, 1]))
    return out


def calibrate(D, n_perm=200):
    """Every test and every agreement, on the real labels and on `n_perm` shuffles."""
    real = D.run()
    out = {"results": real,
           "tests": (summarise_tests(real), []),
           "agreement": (agreement(real, D.protein_tests()), [])}
    for lab in D.permutations(n_perm):
        res = D.run(lab)
        out["tests"][1].append(summarise_tests(res))
        out["agreement"][1].append(agreement(res, D.protein_tests(lab)))
    return out


def perm_p(real, null):
    """One-sided permutation p for a value at least as large as `real`."""
    return (1 + sum(v >= real for v in null)) / (1 + len(null))


def line_key(c, x, y, colour, dashed_key=False):
    """A short line swatch before a direct label, so the text stays in ink."""
    if dashed_key:
        for k in range(3):
            c.line(x + k * 6, y, x + k * 6 + 3.5, y, stroke=colour, sw=1.6)
    else:
        c.line(x, y, x + 15, y, stroke=colour, sw=2.4)


def panel_da_control(D, site, label, outdir, font, letter=""):
    """One site's level against its protein's, per patient: the positive control."""
    import extra_differential as df
    sid = next((i for i, s in D.S.items() if site in s["keys"]), None)
    dig = next((d for d in ORDER if sid in D.level.get(d, {})
                and sid in D.prot.get(d, {})), None)
    if sid is None or dig is None:
        sys.exit(f"{site}: no quantified site with a protein to compare")
    fr, g = D.frames[dig], D.prot[dig][sid]
    lv = D.level[dig][sid]
    pts = [(s, fr["logs"][g][s], lv[s]) for s in sorted(lv) if s in fr["logs"][g]]
    ctrl = [p for p in pts if D.cond.get(p[0]) == D.control]
    mx = statistics.mean(p[1] for p in ctrl)
    my = statistics.mean(p[2] for p in ctrl)
    pts = [(s, x - mx, y - my, D.cond.get(s)) for s, x, y in pts]
    res = {k: next(r for r in rs if r["group"] == sid)
           for k, rs in D.run().items() if k[1] == dig
           and any(r["group"] == sid for r in rs)}
    pa = [x for _, x, _, cnd in pts if cnd == D.case]
    pb = [x for _, x, _, cnd in pts if cnd == D.control]
    fc_prot = statistics.mean(pa) - statistics.mean(pb)
    lo = math.floor(min(min(p[1], p[2]) for p in pts) - 0.5)
    hi = math.ceil(max(max(p[1], p[2]) for p in pts) + 0.5)
    W, H = PANEL_W, 360
    c = Canvas(W, H, font)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    side = H - 78 - 56
    ax = Axes(c, 96, 78, side, side, (lo, hi), (lo, hi))
    header(c, 48, 32,
           f"{label} {D.S[sid]['residue']}{site.rsplit(':', 1)[1]} falls with its "
           f"protein; its occupancy does not change",
           f"{display(dig)} · {len(pa)} LBD v {len(pb)} Control patients · log2 "
           f"relative to the Control mean")
    ticks = list(range(lo, hi + 1))
    ax.ygrid(ticks, lambda v: f"{v:+.0f}" if v else "0")
    ax.xaxis(ticks, lambda v: f"{v:+.0f}" if v else "0")
    c.line(ax.X(lo), ax.Y(lo), ax.X(hi), ax.Y(hi), stroke=INK_MUTED, sw=1.2)
    c.text(ax.X(lo + 0.35), ax.Y(lo + 0.85), "constant occupancy", 8.4, INK_MUTED,
           "start", rot=-45)
    # Control first, so cases draw on top
    for s, x, y, cnd in sorted(pts, key=lambda p: p[3] == D.case):
        case = cnd == D.case
        c.add(f'<circle cx="{ax.X(x):.1f}" cy="{ax.Y(y):.1f}" r="5" '
              f'fill="{INK if case else "#ffffff"}" stroke="{INK if not case else "#ffffff"}" '
              f'stroke-width="{1.4 if not case else 2}"><title>CF_{s} {cnd}: protein '
              f'{x:+.2f}, site {y:+.2f}</title></circle>')
    kx, ky = ax.x + ax.w + 36, ax.y + 8
    dot(c, kx + 5, ky, 5, INK, "#ffffff", 2)
    c.text(kx + 16, ky + 3.5, "LBD", 9.5, INK, "start")
    c.add(f'<circle cx="{kx + 5:.1f}" cy="{ky + 18:.1f}" r="5" fill="#ffffff" '
          f'stroke="{INK}" stroke-width="1.4"/>')
    c.text(kx + 16, ky + 21.5, "Control", 9.5, INK, "start")
    rows = [("protein", fc_prot, None),
            ("phosphosite", res[("level", dig)]["log2fc"], res[("level", dig)]["p"]),
            ("occupancy", res[("occupancy", dig)]["log2fc"],
             res[("occupancy", dig)]["p"])]
    c.text(kx, ky + 58, "LBD − Control, log2", 9, INK_MUTED, "start")
    for k, (name, fc, p) in enumerate(rows):
        y = ky + 78 + 17 * k
        c.text(kx, y, name, 9.5, INK, "start")
        c.text(kx + 100, y, f"{fc:+.2f}", 9.5, INK, "end", "600")
        if p is not None:
            c.text(kx + 108, y, f"p {p:.2g}", 8.6, INK_MUTED, "start")
    c.text(ax.x + ax.w / 2, H - 14, f"{label} protein (unmodified search)", 10,
           INK, "middle")
    ylabel(c, 34, ax.y + ax.h / 2, "Phosphosite (phospho search)")
    name = f"phospho_da_{label.lower()}.svg"
    save(c, outdir, name)
    return pts, rows


DA_ROWS = [("level", "Combined", "all digests"), ("level", "GluC", "Glu-C"),
           ("level", "LysC", "Lys-C"), ("level", "Trypsin", "Trypsin"),
           ("occupancy", "Combined", "all digests"), ("occupancy", "GluC", "Glu-C"),
           ("occupancy", "LysC", "Lys-C"), ("occupancy", "Trypsin", "Trypsin")]


def ladder(c, ax, rows, real, null, fmt, step):
    """Rows of shuffled-label values (grey) with the real value (ink) on each."""
    import random
    ys = []
    for k, (key, name, head) in enumerate(rows):
        y = ax.y + step * (k + 0.5)
        ys.append(y)
        vals = null.get(key, [])
        rng = random.Random(k)
        for v in vals:
            c.add(f'<circle cx="{ax.X(v):.1f}" cy="{y + rng.uniform(-5, 5):.1f}" '
                  f'r="2" fill="{INK_MUTED}" fill-opacity="0.35"/>')
        if key in real:
            v = real[key]
            dot(c, ax.X(v), y, 5, INK, "#ffffff", 2)
            c.text(ax.x + ax.w + 12, y + 3.5, fmt(key, v, vals), 8.8,
                   INK_SECONDARY, "start")
        c.text(ax.x - 12, y + 3.5, name, 9.5, INK, "end")
        if head:
            c.text(12, y + 3.5, head, 9.5, INK, "start", "600")
    return ys


def panel_da_calibration(cal, outdir, font, letter=""):
    """Share of sites at p < 0.05, real labels against 200 shuffles, per test."""
    real_s, null_s = cal["tests"]
    real = {k: v[1] for k, v in real_s.items()}
    null = collections.defaultdict(list)
    for n in null_s:
        for k, v in n.items():
            null[k].append(v[1])
    rows = [((q, d), f"{name} · {real_s[(q, d)][0]:,} sites",
             ("site level" if q == "level" else "occupancy") if d == "Combined"
             else "") for q, d, name in DA_ROWS if (q, d) in real_s]
    step = 26.0
    top = max(v for vs in null.values() for v in vs)
    top = max(top, max(real.values())) * 1.08
    W = PANEL_W
    H = round(92 + step * len(rows) + 52)
    c = Canvas(W, H, font)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    ax = Axes(c, 238, 92, W - 238 - 150, step * len(rows), (0, top), (0, 1))
    header(c, 48, 32, "No phosphosite test beats shuffled labels once eight are counted",
           "share of tested sites at p < 0.05 · grey: 200 balanced label shuffles "
           "· black: LBD v Control")
    ticks = [t / 100 for t in range(0, int(top * 100) + 1, 5)]
    for t in ticks:
        c.line(ax.X(t), ax.y, ax.X(t), ax.y + ax.h, stroke="#e1e0d9", sw=1)
    c.line(ax.X(0.05), ax.y - 4, ax.X(0.05), ax.y + ax.h, stroke=INK_MUTED, sw=1)
    c.text(ax.X(0.05), ax.y - 8, "5%", 8.4, INK_MUTED, "middle")
    ladder(c, ax, rows, real, null,
           lambda k, v, vals: f"{v:.1%} · perm p {perm_p(v, vals):.2f}", step)
    for k in range(1, len(rows)):
        if rows[k][2]:
            y = ax.y + step * k
            c.line(12, y, ax.x + ax.w, y, stroke="#e1e0d9", sw=1)
    ax.xaxis(ticks, lambda v: f"{v:.0%}")
    c.text(ax.x + ax.w / 2, H - 10, "Sites at p < 0.05", 10, INK, "middle")
    save(c, outdir, "phospho_da_calibration.svg")
    return rows


def panel_da_agreement(cal, outdir, font, letter=""):
    """Between-digest fold-change correlation, real labels against shuffles."""
    real_a, null_a = cal["agreement"]
    real = {k: v[1] for k, v in real_a.items()}
    null = collections.defaultdict(list)
    for n in null_a:
        for k, v in n.items():
            null[k].append(v[1])
    heads = {"protein": "proteins", "level": "site level", "occupancy": "occupancy"}
    rows = []
    for name in ("protein", "level", "occupancy"):
        for j, (a, b) in enumerate(PAIRS):
            key = (name, a, b)
            if key in real_a:
                rows.append((key, f"{display(a)} – {display(b)} · {real_a[key][0]:,}",
                             heads[name] if j == 0 else ""))
    step = 24.0
    vals = [v for vs in null.values() for v in vs] + list(real.values())
    lo = min(-0.2, math.floor(min(vals) * 10) / 10)
    hi = max(0.8, math.ceil(max(vals) * 10) / 10)
    W = PANEL_W
    H = round(92 + step * len(rows) + 52)
    c = Canvas(W, H, font)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    ax = Axes(c, 262, 92, W - 262 - 130, step * len(rows), (lo, hi), (0, 1))
    header(c, 48, 32, "Digests agree about patients, not about the disease",
           "correlation of LBD − Control fold changes between two digests · grey: "
           "200 label shuffles · black: real")
    ticks = [round(lo + 0.2 * k, 1) for k in range(int(round((hi - lo) / 0.2)) + 1)]
    for t in ticks:
        c.line(ax.X(t), ax.y, ax.X(t), ax.y + ax.h, stroke="#e1e0d9", sw=1)
    c.line(ax.X(0), ax.y, ax.X(0), ax.y + ax.h, stroke=INK_MUTED, sw=1)
    ladder(c, ax, rows, real, null,
           lambda k, v, vs: f"r {v:+.2f} · perm p {perm_p(v, vs):.2f}", step)
    for k in range(1, len(rows)):
        if rows[k][2]:
            y = ax.y + step * k
            c.line(12, y, ax.x + ax.w, y, stroke="#e1e0d9", sw=1)
    ax.xaxis(ticks, lambda v: f"{v:+.1f}" if v else "0")
    c.text(ax.x + ax.w / 2, H - 10, "Pearson r between the two digests' fold changes",
           10, INK, "middle")
    save(c, outdir, "phospho_da_agreement.svg")
    return rows


def write_da(D, cal, genes, path):
    """The table view: every tested site, level and occupancy, per test."""
    res = cal["results"]
    by = collections.defaultdict(dict)
    for k, rs in res.items():
        for r in rs:
            by[r["group"]][k] = r
    cols = [(q, d) for q, d, _ in DA_ROWS]
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["site", "gene", "residue", "reference_peptides"] +
                   [f"{q}_{d}_{f}" for q, d in cols
                    for f in ("log2fc", "p", "q", "n_case", "n_ctrl")])
        order = sorted(by, key=lambda s: min(r["p"] for r in by[s].values()))
        for sid in order:
            acc, pos = sorted(D.S[sid]["keys"])[0].rsplit(":", 1)
            peps = ";".join(f"{d}:{D.ref[d][sid]}" for d in ORDER
                            if sid in D.ref.get(d, {}))
            row = [f"{acc}:{pos}", genes.get(acc, acc),
                   f"{D.S[sid]['residue']}{pos}", peps]
            for k in cols:
                r = by[sid].get(k)
                row += ([f"{r['log2fc']:.4f}", f"{r['p']:.3e}", f"{r['q']:.3e}",
                         r["n_case"], r["n_ctrl"]] if r else [""] * 5)
            w.writerow(row)
    print(f"  wrote {path}  ({len(by):,} sites)")


# ------------------------------------------------------------------ report ---
def shallow_runs(scan_dir, frac=0.5):
    """-> {patient: digest} for kept runs below `frac` of their digest's median."""
    runs = [r for r in load_runs(scan_dir) if r["kept"] == "1"]
    out = {}
    for d in ORDER:
        n = [int(r["precursors"]) for r in runs if r["digest"] == d]
        if not n:
            continue
        m = statistics.median(n)
        for r in runs:
            if r["digest"] == d and int(r["precursors"]) < frac * m:
                out[r["sample"]] = d
    return out


def report(scan_dir, fasta, cache):
    runs = load_runs(scan_dir)
    dropped = [r for r in runs if r["kept"] != "1"]
    print(f"runs: {len(runs)} scanned, {len(dropped)} dropped by the run gate")
    for r in dropped:
        print(f"    {r['digest']:<8}{r['sample']:<6}{int(r['precursors']):>7,} "
              f"precursors  {r['run']}")
    seqs = read_fasta(fasta)
    co = Cohort(load(scan_dir), known_sites(cache), seqs, complete_patients(runs))
    S = co.S
    print(f"\n{len(S):,} sites on {len(co.patients)} patients measured by all "
          f"three digests ({', '.join(co.patients)})")
    for d in ORDER:
        print(f"  {display(d):<8}{sum(d in v for v in co.digests.values()):>7,}")

    per = co.per_patient()
    order = co.greedy_order()
    tcov = collections.defaultdict(set)
    for r in co.rows:
        if r["digest"] == "Trypsin":
            for a, st in r["where"]:
                tcov[r["sample"]].update(f"{a}:{st + k + 1}"
                                         for k in range(len(r["peptide"])))
    print(f"\nper patient, greedy order {' > '.join(order)}")
    g_all, g_out, g_cur = [], [], []
    for p in co.patients:
        t = per[p]["Trypsin"]
        u = set().union(*per[p].values())
        new = u - t
        outside = [i for i in new if not (S[i]["keys"] & tcov[p])]
        tk = sum(co.curated[i] for i in t)
        uk = sum(co.curated[i] for i in u)
        g_all.append(len(u) / len(t) - 1)
        g_out.append(len(outside) / len(t))
        g_cur.append(uk / tk - 1)
        print(f"  {p}  trypsin {len(t):>5,}  all {len(u):>5,} (+{g_all[-1]:.0%})"
              f"  outside every tryptic phosphopeptide +{g_out[-1]:.0%}"
              f"  curated {tk:,} -> {uk:,} (+{g_cur[-1]:.0%})")
    print(f"  median gain: all sites +{statistics.median(g_all):.0%}, outside "
          f"tryptic phosphopeptides +{statistics.median(g_out):.0%}, curated "
          f"sites +{statistics.median(g_cur):.0%}")

    print("\ncurated share, observed / random placement in the same peptides")
    for k in (1, 2, 3):
        ids = [i for i in S if len(co.digests[i]) == k]
        o = sum(co.curated[i] for i in ids) / len(ids)
        z = sum(co.null[i] for i in ids) / len(ids)
        print(f"  {k} digest{'s' if k > 1 else ' '}  n {len(ids):>5,}  {o:6.1%} / "
              f"{z:6.1%}   lift {(o - z) / (1 - z):5.1%}")
    combos = collections.Counter(tuple(sorted(v)) for v in co.digests.values())
    for combo, n in sorted(combos.items(), key=lambda kv: (len(kv[0]), kv[0])):
        ids = [i for i in S if tuple(sorted(co.digests[i])) == combo]
        o = sum(co.curated[i] for i in ids) / n
        print(f"    {'+'.join(display(d) for d in combo):<26}n {n:>5,}  {o:6.1%}")
    for lab, lo, hi, cell in co.strata():
        print(f"  seen in {lab:>5} patients  " + "   ".join(
            f"{cls} n {cell[cls][0]:>4,} {cell[cls][1]:6.1%} / {cell[cls][2]:6.1%}"
            if cell[cls][0] else f"{cls} n    0" for cls in ("single", "multi")))

    print("\nresidue, observed / random placement in the same peptides")
    for k, obs, rnd, n in residue_mix(co):
        print(f"  {k} digest{'s' if k > 1 else ' '}  " + "  ".join(
            f"p{x} {obs[x] / n:5.1%} / {rnd[x] / n:5.1%}" for x in STY))

    pairs = co.pairs()
    a = sum(p[0] for p in pairs) / len(pairs)
    ch = sum(1 / p[1] for p in pairs) / len(pairs)
    forced = [p for p in co.pairs(min_candidates=1) if p[1] == 1]
    top, n_top = collections.Counter(p[2] for p in pairs).most_common(1)[0]
    print(f"\n{len(pairs):,} overlapping pairs from two digests, >= 2 candidates: "
          f"agree {a:.1%}, chance {ch:.1%}; {len(forced):,} forced pairs, "
          f"{sum(p[0] for p in forced):,} agree; most pairs from one entry: "
          f"{top} {n_top / len(pairs):.0%}")

    print("\ncontrols: one digest / two or three, curated share and pY")
    qty = {i: statistics.median(co.rows[o[3]]["quantity"] or 0
                                for o in S[i]["obs"]) for i in S}
    best_q = {i: min(co.rows[o[3]]["q"] for o in S[i]["obs"]) for i in S}
    forced = {i: min(sum(c in STY for c in co.rows[o[3]]["peptide"])
                     for o in S[i]["obs"]) == 1 for i in S}
    cut = sorted(qty.values())
    t1, t2 = cut[len(cut) // 3], cut[2 * len(cut) // 3]
    tests = [("intensity, lowest third", lambda i: qty[i] < t1),
             ("intensity, middle third", lambda i: t1 <= qty[i] < t2),
             ("intensity, highest third", lambda i: qty[i] >= t2),
             ("best q < 1e-6", lambda i: best_q[i] < 1e-6),
             ("one S/T/Y in the peptide", lambda i: forced[i])]
    entry = {i: sorted(S[i]["keys"])[0].rsplit(":", 1)[0] for i in S}
    ranked = [a for a, _ in collections.Counter(entry.values()).most_common()]
    for n in (10, 30):
        drop = set(ranked[:n])
        tests.append((f"without the top {n} entries",
                      lambda i, drop=drop: entry[i] not in drop))
    for name, keep in tests:
        cells = []
        for multi in (False, True):
            ids = [i for i in S if keep(i) and (len(co.digests[i]) > 1) == multi]
            cur = sum(co.curated[i] for i in ids) / len(ids)
            py = sum(S[i]["residue"] == "Y" for i in ids) / len(ids)
            cells.append(f"n {len(ids):>5,} curated {cur:5.1%} pY {py:5.1%}")
        print(f"  {name:<28}" + "   |   ".join(cells))

    def near_miss(i, w=3):
        for k in S[i]["keys"]:
            acc, pos = k.rsplit(":", 1)
            if any(0 < abs(p - int(pos)) <= w for p in co.known.get(acc, {})):
                return True
        return False
    for multi in (False, True):
        ids = [i for i in S if (len(co.digests[i]) > 1) == multi]
        nm = sum(1 for i in ids if not co.curated[i] and near_miss(i)) / len(ids)
        print(f"  {'two or three' if multi else 'one digest'}: uncurated but "
              f"within 3 residues of a curated site {nm:.1%}")


# --------------------------------------------------------------- isoforms ---
def isoform_sites(S, rows, seqs):
    """Sites placed only on isoform entries: phosphopeptides no canonical holds."""
    from fig2cd_isoform_coverage import discriminating
    disc = {}

    def own_residues(acc):
        if acc not in disc:
            base = ISOFORM.sub("", acc)
            disc[acc] = (discriminating(seqs[base], seqs[acc])
                         if base in seqs and acc in seqs else set())
        return disc[acc]

    out = []
    for sid, s in S.items():
        keys = sorted(s["keys"])
        if not all(ISOFORM.search(k.rsplit(":", 1)[0]) for k in keys):
            continue
        own = any(int(k.rsplit(":", 1)[1]) - 1 in own_residues(k.rsplit(":", 1)[0])
                  for k in keys)
        proof = False
        for _, _, pep, i in s["obs"]:
            r = rows[i]
            for acc, st in r["where"]:
                mine = own_residues(acc)
                shared = sum(1 for k, ch in enumerate(pep)
                             if ch in STY and st + k not in mine)
                proof |= own and len(r["offsets"]) > shared
        out.append({"site": sid, "entries": keys, "residue": s["residue"],
                    "own": own, "proof": proof,
                    "digests": frozenset(o[0] for o in s["obs"]),
                    "patients": len({o[1] for o in s["obs"]})})
    return out


def iso_classes(sites_):
    """-> [(label, colour key, n, localisation-proof n)] of isoform-specific sites."""
    own = [s for s in sites_ if s["own"]]
    classes = [("trypsin only", "Trypsin", lambda d: d == {"Trypsin"}),
               ("trypsin and another digest", "All",
                lambda d: "Trypsin" in d and len(d) > 1),
               ("Lys-C only", "LysC", lambda d: d == {"LysC"}),
               ("Glu-C only", "GluC", lambda d: d == {"GluC"}),
               ("Lys-C and Glu-C", "All", lambda d: d == {"LysC", "GluC"})]
    out = []
    for name, key, test in classes:
        sel = [s for s in own if test(s["digests"])]
        if sel:
            out.append((name, key, len(sel), sum(s["proof"] for s in sel)))
    return out


def panel_iso_summary(sites_, outdir, font, letter=""):
    """How many isoform-specific phosphosites need a digest other than trypsin."""
    rows = iso_classes(sites_)
    n = sum(r[2] for r in rows)
    colour = assign(ORDER + ["All"])
    W = PANEL_W
    step = 30.0
    H = round(70 + step * len(rows) + 24)
    c = Canvas(W, H, font)
    if letter:
        c.text(20, 30, letter, 13, INK, "start", "600")
    alone = sum(r[2] for r in rows if "trypsin" not in r[0])
    header(c, 48, 32,
           f"{alone} of {n} isoform-specific phosphosites are placed only by "
           f"Lys-C or Glu-C",
           "phosphate on a residue the isoform has and its canonical does not · "
           "all patients")
    ax = Axes(c, 238, 70, W - 238 - 150, step * len(rows),
              (0, max(r[2] for r in rows) * 1.1), (0, 1))
    bh = 16.0
    for k, (name, key, v, proof) in enumerate(rows):
        y = ax.y + step * (k + 0.5) - bh / 2
        c.rect(ax.x, y, max(ax.X(v) - ax.x, 1), bh, fill=colour[key], fo=0.9, rx=3)
        c.text(ax.x - 12, y + bh / 2 + 3.5, name, 9.5, INK, "end")
        c.text(ax.X(v) + 8, y + bh / 2 + 3.5, f"{v:,}", 9.5, INK, "start", "600")
        if proof:
            c.text(ax.X(v) + 26, y + bh / 2 + 3.5,
                   f"{proof} hold whatever the localisation", 8.4, INK_MUTED,
                   "start")
    c.line(ax.x, ax.y - 4, ax.x, ax.y + ax.h + 4, stroke=AXIS, sw=1)
    save(c, outdir, "phospho_isoform_sites.svg")
    return rows


INTERNAL_EVENTS = {"SE", "MXE", "A5SS", "A3SS", "RI"}


def region_sites(S, rows, seqs, events_tsv, known):
    """Phosphosites in a stretch that tells two forms of a protein apart."""
    from fig2cd_isoform_coverage import own_residues
    events = {}
    with open(events_tsv, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            events[r["isoform"]] = set(filter(None, r["events"].split(",")))
    isos = collections.defaultdict(list)
    for a in seqs:
        if ISOFORM.search(a):
            isos[ISOFORM.sub("", a)].append(a)
    cache = {}

    def blocks(base, iso):
        if (base, iso) not in cache:
            cache[(base, iso)] = own_residues(seqs[base], seqs[iso])
        return cache[(base, iso)]

    out = []
    for sid, s in S.items():
        hit = None
        for k in sorted(s["keys"]):
            acc, pos = k.rsplit(":", 1)
            p0 = int(pos) - 1
            base = ISOFORM.sub("", acc)
            if acc != base:
                if base in seqs and acc in seqs and p0 in blocks(base, acc)[1]:
                    hit = ("isoform-own", acc, base, int(pos))
            else:
                for iso in isos.get(acc, ()):
                    if events.get(iso, set()) & INTERNAL_EVENTS and \
                            p0 in blocks(acc, iso)[0]:
                        hit = ("canonical-exon", acc, iso, int(pos))
                        break
            if hit:
                break
        if not hit:
            continue
        kind, acc, against, pos = hit
        peps = collections.defaultdict(set)
        pts = collections.defaultdict(set)
        forced = set()
        for d, smp, pep, i in s["obs"]:
            peps[d].add(pep)
            pts[d].add(smp)
            if len(rows[i]["offsets"]) == 1 and sum(c in STY for c in pep) == 1:
                forced.add(d)
        distinct = len(set().union(*peps.values()))
        out.append({"site": sid, "kind": kind, "acc": acc, "against": against,
                    "position": pos, "residue": s["residue"],
                    "events": ",".join(sorted(events.get(
                        against if kind == "canonical-exon" else acc, ()))),
                    "digests": sorted(peps, key=ORDER.index),
                    "differ": len(peps) > 1 and distinct > max(len(v) for v in
                                                               peps.values()),
                    "patients": {d: len(v) for d, v in pts.items()},
                    "forced": sorted(forced, key=ORDER.index),
                    "curated": any(int(k.rsplit(":", 1)[1]) in
                                   known.get(k.rsplit(":", 1)[0], {})
                                   for k in s["keys"])})
    return out


def pair_sites(S, seqs, bybase, known):
    """-> [dict] per (pair, site) on either form's own residues."""
    from fig2cd_isoform_coverage import own_residues
    isos = collections.defaultdict(list)
    for a in seqs:
        if ISOFORM.search(a):
            isos[ISOFORM.sub("", a)].append(a)
    cache = {}

    def pair(base, iso):
        if (base, iso) not in cache:
            cseq, iseq = seqs[base], seqs[iso]
            byd = bybase.get(base, {})
            peps = []
            for a, b in ((cseq, iseq), (iseq, cseq)):
                n = {d: sum(p in a and p not in b for p in byd.get(d, ()))
                     for d in ORDER}
                peps.append({d: v for d, v in n.items() if v})
            cache[(base, iso)] = own_residues(cseq, iseq) + tuple(peps)
        return cache[(base, iso)]

    out, seen = [], set()
    for sid, s in S.items():
        peps = collections.defaultdict(set)
        pts = collections.defaultdict(set)
        for d, smp, pep, _ in s["obs"]:
            peps[d].add(pep)
            pts[d].add(smp)
        for k in sorted(s["keys"]):
            acc, pos = k.rsplit(":", 1)
            base = ISOFORM.sub("", acc)
            for iso in ([acc] if acc != base else isos.get(base, [])):
                if base not in seqs or iso not in seqs or (base, iso, sid) in seen:
                    continue
                c_own, i_own, own_c, own_i = pair(base, iso)
                own = i_own if acc != base else c_own
                if int(pos) - 1 not in own or not own_c or not own_i:
                    continue
                seen.add((base, iso, sid))
                distinct = len(set().union(*peps.values()))
                out.append({"site": sid, "canonical": base, "isoform": iso,
                            "acc": acc, "position": int(pos),
                            "residue": s["residue"],
                            "on": "isoform" if acc != base else "canonical",
                            "own_residues": len(own),
                            "own_canonical": own_c, "own_isoform": own_i,
                            "trypsin_sees_both": "Trypsin" in own_c and
                                                 "Trypsin" in own_i,
                            "digests": sorted(peps, key=ORDER.index),
                            "differ": len(peps) > 1 and distinct > max(
                                len(v) for v in peps.values()),
                            "patients": {d: len(v) for d, v in pts.items()},
                            "curated": any(int(x.rsplit(":", 1)[1]) in
                                           known.get(x.rsplit(":", 1)[0], {})
                                           for x in s["keys"])})
    return out


def mod_sites_rows(S, mod="Phospho (STY)"):
    """The sites in `mod_sites.load`'s shape, at each site's first key."""
    out = []
    for sid, s in sorted(S.items(), key=lambda kv: sorted(kv[1]["keys"])[0]):
        acc, pos = sorted(s["keys"])[0].rsplit(":", 1)
        byd = collections.defaultdict(set)
        for d, _, pep, _ in s["obs"]:
            byd[d].add(pep)
        out.append({"site": sid, "acc": acc, "res": int(pos), "mod": mod,
                    "digests": sorted(byd),
                    "by_digest": {d: sorted(v) for d, v in sorted(byd.items())},
                    "peps": sorted(set().union(*byd.values()))})
    return out


def export_sites(S, path, mod="Phospho (STY)"):
    """Write the sites in `mod_sites.load`'s table, so its CLI draws them."""
    rows = mod_sites_rows(S, mod)
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "residue", "modification", "digests",
                    "n_digests", "n_distinct_peptides", "peptides_by_digest"])
        for r in rows:
            w.writerow([r["acc"], r["res"], mod, ",".join(r["digests"]),
                        len(r["digests"]), len(r["peps"]),
                        "|".join(f"{d}:" + ",".join(v)
                                 for d, v in r["by_digest"].items())])
    print(f"  wrote {path}  ({len(rows):,} sites)")


def isoforms(args):
    """Isoform-specific phosphosites, their regions' coverage, and the panels."""
    from lib_fasta import gene_map
    from fig2cd_isoform_coverage import discriminating, peptides_by_base, region_coverage
    seqs = read_fasta(args.fasta)
    genes = gene_map(args.fasta)
    rows = load(args.scan)
    S = sites(rows)
    found = isoform_sites(S, rows, seqs)
    print(f"{len(found):,} of {len(S):,} sites are placed only on isoform entries; "
          f"{sum(s['own'] for s in found):,} on a residue the isoform alone has")
    for name, _, n, proof in iso_classes(found):
        print(f"  {name:<28}{n:>4}   {proof} of them whatever the localisation")
    paths = [p for pat in args.reports for p in (sorted(glob.glob(pat)) or [pat])]
    bybase = peptides_by_base(paths)
    table = []
    for s in found:
        if not s["own"]:
            continue
        for k in s["entries"]:
            acc, pos = k.rsplit(":", 1)
            base = ISOFORM.sub("", acc)
            disc = discriminating(seqs[base], seqs[acc])
            if int(pos) - 1 not in disc:
                continue
            cov = region_coverage(seqs[acc], disc, bybase.get(base, {}))
            table.append((genes.get(acc, acc), acc, f"{s['residue']}{pos}",
                          "+".join(sorted(s["digests"], key=ORDER.index)),
                          s["patients"], len(disc),
                          {d: len(v) / len(disc) for d, v in cov.items()},
                          len(set().union(*cov.values())) / len(disc),
                          s["proof"]))
    os.makedirs(os.path.dirname(os.path.abspath(args.tsv)) or ".", exist_ok=True)
    with open(args.tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "isoform", "site", "digests", "patients",
                    "own_residues"] + [f"covered_{d}" for d in ORDER] +
                   ["covered_all", "whatever_localisation"])
        for g, acc, site, digs, npt, nd, cov, allc, proof in sorted(table):
            w.writerow([g, acc, site, digs, npt, nd] +
                       [f"{cov[d]:.3f}" for d in ORDER] + [f"{allc:.3f}",
                                                           int(proof)])
    print(f"  wrote {args.tsv}  ({len(table):,} isoform-site rows)")
    print("\nisoform-specific sites placed only by Lys-C or Glu-C, in >= 3 patients")
    for g, acc, site, digs, npt, nd, cov, allc, proof in sorted(table,
                                                              key=lambda t: -t[4]):
        if "Trypsin" not in digs and npt >= 3:
            print(f"  {g:<10}{acc:<12}{site:<7}{digs:<11}{npt:>2} patients  own "
                  f"residues {nd:>3}: " + " ".join(f"{display(d)} {cov[d]:.0%}"
                                                   for d in ORDER) +
                  f"  all {allc:.0%}" + ("  localisation-proof" if proof else ""))
    census = region_sites(S, rows, seqs, args.events, known_sites(args.cache))
    multi = [r for r in census if len(r["digests"]) > 1]
    print(f"\n{len(census):,} sites in a stretch that tells two forms apart: "
          f"{sum(r['kind'] == 'isoform-own' for r in census)} on an isoform's own "
          f"residues, {sum(r['kind'] == 'canonical-exon' for r in census)} in a "
          f"canonical exon an isoform skips")
    print(f"  placed by two or more digests {len(multi)}, through different "
          f"peptides {sum(r['differ'] for r in multi)}; only by Lys-C or Glu-C "
          f"{sum('Trypsin' not in r['digests'] for r in census)}")
    with open(args.regions_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "accession", "site", "kind", "against", "events",
                    "digests", "different_peptides", "patients", "forced",
                    "curated"])
        for r in sorted(census, key=lambda r: (-len(r["digests"]), -r["differ"],
                                               -max(r["patients"].values()))):
            w.writerow([genes.get(r["acc"], r["acc"]), r["acc"],
                        f"{r['residue']}{r['position']}", r["kind"], r["against"],
                        r["events"], "+".join(r["digests"]), int(r["differ"]),
                        ";".join(f"{d}:{n}" for d, n in r["patients"].items()),
                        "+".join(r["forced"]), int(r["curated"])])
    print(f"  wrote {args.regions_tsv}")
    print("  placed by three digests through different peptides, >= 5 patients:")
    for r in sorted(multi, key=lambda r: -max(r["patients"].values())):
        if len(r["digests"]) == 3 and r["differ"] and \
                max(r["patients"].values()) >= 5:
            print(f"    {genes.get(r['acc'], r['acc']):<9}{r['acc']:<11}"
                  f"{r['residue']}{r['position']:<6}{r['kind']:<15}v {r['against']:<11}"
                  f"[{r['events']}] patients "
                  + " ".join(f"{d}:{n}" for d, n in r["patients"].items())
                  + (f"  forced in {'+'.join(r['forced'])}" if r["forced"] else "")
                  + ("  curated" if r["curated"] else ""))
    drawable = pair_sites(S, seqs, bybase, known_sites(args.cache))
    blind = [r for r in drawable if not r["trypsin_sees_both"]]
    three = [r for r in blind if len(r["digests"]) == 3]

    def npairs(rs):
        return len({(r["canonical"], r["isoform"]) for r in rs})
    print(f"\n{len(drawable):,} (pair, site) combinations a pair panel can draw -- "
          f"both forms seen by peptides of their own, the site on one form's own "
          f"residues -- over {npairs(drawable)} pairs")
    print(f"  one form seen only by Lys-C or Glu-C: {len(blind)} over "
          f"{npairs(blind)} pairs; of those, placed by all three digests: "
          f"{len(three)}")
    for r in sorted(three, key=lambda r: (r["own_residues"], r["position"])):
        print(f"    {genes.get(r['canonical'], r['canonical']):<9}{r['canonical']:<9}"
              f"v {r['isoform']:<11}{r['residue']}{r['position']:<6}on the "
              f"{r['on']:<10} own residues {r['own_residues']:>5,}  patients "
              + " ".join(f"{d}:{n}" for d, n in r["patients"].items())
              + ("  curated" if r["curated"] else ""))
    with open(args.pairs_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "canonical", "isoform", "site", "on", "own_residues",
                    "own_peptides_canonical", "own_peptides_isoform",
                    "trypsin_sees_both", "digests", "different_peptides",
                    "patients", "curated"])
        for r in sorted(drawable, key=lambda r: (-len(r["digests"]),
                                                 r["trypsin_sees_both"],
                                                 r["own_residues"])):
            w.writerow([genes.get(r["canonical"], r["canonical"]), r["canonical"],
                        r["isoform"], f"{r['residue']}{r['position']}", r["on"],
                        r["own_residues"],
                        ";".join(f"{d}:{n}" for d, n in r["own_canonical"].items()),
                        ";".join(f"{d}:{n}" for d, n in r["own_isoform"].items()),
                        int(r["trypsin_sees_both"]), "+".join(r["digests"]),
                        int(r["differ"]),
                        ";".join(f"{d}:{n}" for d, n in r["patients"].items()),
                        int(r["curated"])])
    print(f"  wrote {args.pairs_tsv}")
    panel_iso_summary(found, args.outdir, args.font, "a")


def site_marks(rows, canonical, isoform, event):
    """The phosphates placed on either form's own block of `event`."""
    _, cs, ce, is_, ie_ = event
    own = {"canonical": (canonical, range(cs, ce)),
           "isoform": (isoform, range(is_, ie_))}
    marks = {form: collections.defaultdict(set) for form in own}
    pts = collections.defaultdict(set)
    for r in rows:
        for a, st in r["where"]:
            for form, (acc, span) in own.items():
                if a != acc:
                    continue
                for o in r["offsets"]:
                    if st + o in span:
                        marks[form][st + o].add(r["digest"])
                        pts[(form, st + o, r["digest"])].add(r["sample"])
    return {form: dict(v) for form, v in marks.items() if v}, pts


def strip(args):
    """Isoform strip with the phosphosites on each form's own residues."""
    import fig2b_isoform_strip as ie
    seqs = read_fasta(args.fasta)
    event = ie.whole_event(seqs[args.canonical], seqs[args.isoform])
    if event is None:
        sys.exit(f"{args.isoform}: identical to {args.canonical}")
    marks, pts = site_marks(load(args.scan), args.canonical, args.isoform, event)
    paths = [p for pat in args.reports for p in (sorted(glob.glob(pat)) or [pat])]
    ie.draw_pair(paths, seqs, args.canonical, args.isoform, args.out, args.font,
                 args.label, args.letter, events_tsv=args.events,
                 curation=ie.curated(args.varseq, seqs), letters=True,
                 browser=True, event=event, marks=marks, region_label=False,
                 dedupe=True,
                 marks_key="phosphosite on a residue only this form has · "
                           "one dot per digest that placed a phosphate there")
    for form, sites_ in marks.items():
        acc = args.canonical if form == "canonical" else args.isoform
        for pos in sorted(sites_):
            print(f"  {acc:<11}{seqs[acc][pos]}{pos + 1:<6}" + "  ".join(
                f"{display(d)} {len(pts[(form, pos, d)])} patients"
                for d in sorted(sites_[pos], key=ORDER.index)))


def one_side(canon_points, points, min_points=3):
    """True when a column's points all sit on one side of the canonical's mean."""
    mid = statistics.mean(v for _p, _d, v in canon_points)
    rel = [v - mid for _p, _d, v in points]
    return len(rel) >= min_points and (max(rel) < 0 or min(rel) > 0)


def census_picks(facets, n=4, census_only=True):
    """-> ["ACC:RESIDUE"]: census sites on one side, largest |Δ| first, one per gene."""
    out, seen = [], set()
    for f in sorted((f for f in facets if (f["census"] or not census_only)
                     and f["one_side"]),
                    key=lambda f: -abs(f["st"]["mean"])):
        if f["gene"] not in seen:
            seen.add(f["gene"])
            out.append(f"{f['acc']}:{f['label'][2:]}")
    return out[:n]


def dots(args):
    """Dot plots, one facet per phosphosite, against the gene's canonical peptides."""
    import fig2a_supp1_diagnostic_peptides as dp
    from lib_fasta import gene_map
    from lib_report import read_metadata
    seqs = read_fasta(args.fasta)
    genes = gene_map(args.fasta)
    cond = read_metadata(args.metadata)
    rows = load(args.scan)
    S = sites(rows)
    by_key = {k: sid for sid, s in S.items() for k in s["keys"]}
    if args.screen:
        todo = sorted(S)
    else:
        if not args.sites:
            sys.exit("--sites ACC:RESIDUE,... or --screen")
        todo = []
        for k in (x.strip() for x in args.sites.split(",")):
            if k not in by_key:
                sys.exit(f"{k}: no phosphosite placed there")
            todo.append(by_key[k])
    first = {sid: sorted(S[sid]["keys"])[0].rsplit(":", 1) for sid in todo}
    gene_of = {sid: genes.get(acc) or genes.get(ISOFORM.sub("", acc))
               for sid, (acc, _) in first.items()}

    def expand(pats):
        return [p for pat in pats for p in (sorted(glob.glob(pat)) or [pat])]

    def column(members, idx, med):
        return dp.column(members, idx, med, cond, args.case, args.control,
                         args.min_group, "vs-canonical")

    want = None if args.screen else set(gene_of.values())
    idx, med, bygene = dp.read_reports(expand(args.reports), seqs, genes,
                                       want=want)
    pidx, pmed, pby = dp.read_reports(expand(args.phospho), seqs, genes,
                                      want=want, key="modified")
    census = {}
    if os.path.exists(args.pairs_tsv):
        with open(args.pairs_tsv, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                acc = r["canonical"] if r["on"] == "canonical" else r["isoform"]
                census.setdefault(f"{acc}:{r['site'][1:]}", set()).add(
                    r["isoform"])
    cans, diags, facets = {}, {}, []
    for sid in todo:
        gene = gene_of[sid]
        if gene not in bygene or gene not in pby:
            continue
        if gene not in cans:
            base, diags[gene] = dp.diagnostic(gene, bygene[gene], seqs, genes,
                                              idx)
            cans[gene] = (base[0], column(base[1], idx, med)) if base else \
                (None, None)
        canon, can = cans[gene]
        mods = {rows[i]["modified"] for _d, _s, _p, i in S[sid]["obs"]}
        members = {(g, m) for g in pby[gene] for d in ORDER
                   for m in pidx.get((g, d), {}) if m in mods}
        st = column(members, pidx, pmed) if can else None
        if st is None:
            continue
        dp.apply_vs_canonical(can, st)
        acc, pos = first[sid]
        # the isoform's diagnostic peptides, from its column with the most points
        isos = set().union(*(census.get(k, set()) for k in S[sid]["keys"]))
        iso, iso_col = None, None
        if isos and not args.screen and not args.no_isoform:
            cands = [(who, column(m, idx, med)) for who, m in diags[gene].items()
                     if isos & set(who)]
            cands = [(w, c) for w, c in cands if c]
            if cands:
                who, iso_col = max(cands, key=lambda t: len(t[1]["points"]))
                dp.apply_vs_canonical(can, iso_col)
                iso = ";".join(who)
        facets.append({"site": sid, "gene": gene, "canon_acc": canon,
                       "canon": can, "acc": acc, "iso": iso, "iso_st": iso_col,
                       "label": f"p{S[sid]['residue']}{pos}", "st": st})

    if args.screen:
        ok = [f for f in facets if f["st"]["npep"] >= args.min_pep
              and len(f["st"]["per_digest"]) >= args.min_digests]
        from scipy import stats
        tested = [f for f in ok if f["st"]["p"] == f["st"]["p"]]
        for f, q in zip(tested, stats.false_discovery_control(
                [f["st"]["p"] for f in tested], method="bh")):
            f["q"] = float(q)
        # penalise shifts the digests disagree about
        ok.sort(key=lambda f: -(abs(f["st"]["mean"]) - 1.5 * f["st"]["sd"]))
        for f in ok:
            f["census"] = any(k in census for k in S[f["site"]]["keys"])
            f["one_side"] = one_side(f["canon"]["points"], f["st"]["points"])
        inc = [f for f in ok if f["census"]]
        print(f"\n{len(facets):,} sites with a column; {len(ok)} with >= "
              f"{args.min_pep} tested phosphopeptides in >= {args.min_digests} "
              f"digests, {len(inc)} of them in the pair census; best q "
              f"{min((f.get('q', 1) for f in ok), default=float('nan')):.2f}")
        for f in ok[:args.top]:
            st = f["st"]
            print(f"  {'census ' if f['census'] else '       '}"
                  f"{'one side ' if f['one_side'] else '         '}{f['gene']:<10}"
                  f"{f['acc']}:{f['label'][2:]:<7}Δ {st['mean']:+.2f}  p "
                  f"{st['p']:.3f}  q {f.get('q', float('nan')):.2f}  sd "
                  f"{st['sd']:.2f}  {len(st['points'])} points  " + " ".join(
                      f"{display(d)} {v:+.2f}" for d, v in
                      sorted(st["per_digest"].items(), key=lambda kv: ORDER.index(kv[0]))))
        print(f"  {sum(f['one_side'] for f in inc)} census sites have >= 3 points "
              f"on one side of the canonical; the four genes with the largest "
              f"|Δ|: --sites " + ",".join(census_picks(ok)))
        print(f"  {sum(f['one_side'] for f in ok)} sites of any kind do; the four "
              f"genes with the largest |Δ|: --sites "
              + ",".join(census_picks(ok, census_only=False)))
        if args.tsv:
            with open(args.tsv, "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                w.writerow(["gene", "site", "pair_census", "points_one_side",
                            "delta_vs_canonical", "p", "q", "sd_between_digests",
                            "points", "tested_phosphopeptides", "digests"])
                for f in ok:
                    st = f["st"]
                    w.writerow([f["gene"], f"{f['acc']}:{f['label'][2:]}",
                                int(f["census"]), int(f["one_side"]),
                                f"{st['mean']:.4f}",
                                f"{st['p']:.4g}", f"{f.get('q', float('nan')):.4g}",
                                f"{st['sd']:.4f}", len(st["points"]), st["npep"],
                                "+".join(sorted(st["per_digest"], key=ORDER.index))])
            print(f"  wrote {args.tsv}")
        return

    for f in facets:
        st, it = f["st"], f["iso_st"]
        print(f"  {f['gene']:<9}{f['acc']}:{f['label'][2:]:<7}Δ {st['mean']:+.2f}  "
              f"p {st['p']:.3f}  {len(st['points'])} points against "
              f"{len(f['canon']['points'])}" + (
                  f"   isoform {f['iso']} Δ {it['mean']:+.2f} p {it['p']:.3f} "
                  f"{len(it['points'])} points" if it else "   no isoform column"))
    if args.tsv:
        os.makedirs(os.path.dirname(os.path.abspath(args.tsv)) or ".",
                    exist_ok=True)
        with open(args.tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["gene", "column", "kind", "digest", "peptide",
                        "log2fc_case_minus_control"])
            for f in facets:
                cols = [(f["canon_acc"], "canonical", f["canon"])]
                if f["iso_st"]:
                    cols.append((f["iso"], "diagnostic", f["iso_st"]))
                cols.append((f"{f['acc']}:{f['label']}", "phospho", f["st"]))
                for lab, kind, st in cols:
                    for pep, dig, v in st["points"]:
                        w.writerow([f["gene"], lab, kind, dig, pep, f"{v:.4f}"])
        print(f"  wrote {args.tsv}")
    vals = [v for f in facets for st in (f["st"], f["iso_st"]) if st
            for _p, _d, v in st["points"]]
    ylim = tuple(float(v) for v in args.ylim.split(",")) if args.ylim else \
        (min(vals + [-1.75]) - 0.25, max(vals + [1.75]) + 0.25)
    dp.panel([{"gene": f["gene"], "canon_acc": f["canon_acc"], "canon": f["canon"],
               "diag": ([(f["iso"], f["iso_st"])] if f["iso_st"] else [])
               + [(f["label"], f["st"])], "per_column": True,
               "sub": (f"{f['label'][1:]} · absent from {f['iso'].split(';')[0]}"
                       if f["iso_st"] else None)} for f in facets],
             args.out, args.font, args.letter, key=False, ylim=ylim,
             subtitle=False, hue=True, panel_h=300, case_label=args.case,
             control_label=args.control,
             title=("Differential abundance of phosphosites, against their "
                    "protein" if args.no_isoform else
                    "Differential abundance of isoform-discriminating "
                    "phosphosites"))


def da(args):
    """Build, calibrate, print, write and draw the differential analysis."""
    from lib_fasta import gene_map
    from lib_report import read_metadata
    cond = read_metadata(args.metadata)
    rows = load(args.scan)
    S = sites(rows)
    frames = protein_frames(args.quant)
    D = Differential(rows, S, frames, cond)
    genes = gene_map(args.fasta)
    for dig, fr in frames.items():
        a, b = D.groups(cond, fr["keep"])
        print(f"{display(dig):<8} {len(a)} LBD v {len(b)} Control after the protein "
              f"test's run QC (dropped {', '.join(fr['dropped']) or 'none'}); "
              f"{len(D.level[dig]):,} sites with a level, "
              f"{len(D.occupancy[dig]):,} with an occupancy")
    cal = calibrate(D, args.permutations)
    real_s, null_s = cal["tests"]
    print(f"\ntests, real labels against {args.permutations} shuffles")
    for q, d, name in DA_ROWS:
        if (q, d) not in real_s:
            continue
        n, frac, nq, bq = real_s[(q, d)]
        nf = [x[(q, d)][1] for x in null_s if (q, d) in x]
        nn = [x[(q, d)][2] for x in null_s if (q, d) in x]
        print(f"  {q:<10}{name:<12}{n:>6,} sites  p<0.05 {frac:6.1%} "
              f"(null median {statistics.median(nf):.1%}, perm p {perm_p(frac, nf):.3f})"
              f"  q<0.10 {nq} (null 95th {sorted(nn)[int(0.95 * len(nn))]})"
              f"  best q {bq:.3f}")
    real_a, null_a = cal["agreement"]
    print("\nfold-change agreement between digests")
    for key, (n, r) in real_a.items():
        nr = [x[key][1] for x in null_a if key in x]
        print(f"  {key[0]:<10}{key[1]}-{key[2]:<9}{n:>6,} shared  r {r:+.3f}  "
              f"(null median {statistics.median(nr):+.3f}, perm p {perm_p(r, nr):.3f})")
    write_da(D, cal, genes, args.tsv)
    pts, ctrl_rows = panel_da_control(D, args.control_site, args.control_label,
                                      args.outdir, args.font, "a")
    print(f"\npositive control {args.control_label} ({args.control_site}): " +
          ", ".join(f"{n} {fc:+.2f}" + (f" (p {p:.2g})" if p is not None else "")
                    for n, fc, p in ctrl_rows))
    panel_da_calibration(cal, args.outdir, args.font, "b")
    panel_da_agreement(cal, args.outdir, args.font)
    for key in (("level", "Combined"), ("occupancy", "Combined")):
        print(f"\nbest {key[0]} sites, all digests")
        for r in sorted(cal["results"][key], key=lambda r: r["p"])[:6]:
            acc, pos = sorted(S[r["group"]]["keys"])[0].rsplit(":", 1)
            print(f"  {genes.get(acc, acc):<9}{S[r['group']]['residue']}{pos:<6}"
                  f"log2FC {r['log2fc']:+.2f}  p {r['p']:.1e}  q {r['q']:.2f}  "
                  f"n {r['n_case']} v {r['n_ctrl']}")


# --------------------------------------------------------------------- main ---
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("scan", help="read the reports and cache the tables")
    s.add_argument("reports", nargs="+")
    s.add_argument("--fasta", required=True)
    s.add_argument("--outdir", default="derived/mods")
    s.add_argument("--precursor-q", type=float, default=0.01)
    s.add_argument("--min-run-precursors", type=int, default=1000)
    k = sub.add_parser("known", help="download UniProt's phospho annotation")
    k.add_argument("--cache", default="data/uniprot")
    k.add_argument("--force", action="store_true")
    q = sub.add_parser("strip",
                       help="isoform strip with the phosphosites")
    q.add_argument("reports", nargs="+",
                   help="the unmodified search")
    q.add_argument("--scan", default="derived/mods")
    q.add_argument("--fasta", required=True)
    q.add_argument("--canonical", required=True)
    q.add_argument("--isoform", required=True)
    q.add_argument("--label", default=None)
    q.add_argument("--letter", default="")
    q.add_argument("--events", default="data/gencode/events.tsv",
                   help="splice-event typing from prep_splice_events.py")
    q.add_argument("--varseq", default="data/uniprot",
                   help="UniProt VAR_SEQ cache from prep_uniprot.py")
    q.add_argument("--out", required=True)
    q.add_argument("--font", default=FONT)
    o = sub.add_parser("dots", help="dot plots, a facet per phosphosite")
    o.add_argument("reports", nargs="+", help="the unmodified search")
    o.add_argument("--phospho", nargs="+", default=["data/search/*-Phospho.parquet"],
                   help="the phospho search")
    o.add_argument("--scan", default="derived/mods")
    o.add_argument("--fasta", required=True)
    o.add_argument("--metadata", required=True)
    o.add_argument("--sites", default=None, metavar="ACC:RES,...",
                   help="the facets, in order")
    o.add_argument("--screen", action="store_true",
                   help="rank every site instead of drawing")
    o.add_argument("--pairs-tsv", default="derived/da_iso/phospho_pair_sites.tsv",
                   help="pair census from `isoforms`")
    o.add_argument("--top", type=int, default=25)
    o.add_argument("--case", default="LBD")
    o.add_argument("--control", default="Control")
    o.add_argument("--min-group", type=int, default=3)
    o.add_argument("--min-pep", type=int, default=2)
    o.add_argument("--min-digests", type=int, default=2)
    o.add_argument("--ylim", default=None, metavar="LO,HI",
                   help="fix the y-scale")
    o.add_argument("--no-isoform", action="store_true",
                   help="omit the isoform column")
    o.add_argument("--letter", default="")
    o.add_argument("--out", default=None)
    o.add_argument("--tsv", default=None)
    o.add_argument("--font", default=FONT)
    a = sub.add_parser("da", help="LBD v Control phosphosite tests, calibrated")
    a.add_argument("--scan", default="derived/mods")
    a.add_argument("--quant", default="quant",
                   help="protein matrices from `extra_differential.py extract`")
    a.add_argument("--metadata", required=True)
    a.add_argument("--fasta", required=True)
    a.add_argument("--permutations", type=int, default=200)
    a.add_argument("--control-site", default="Q01959:53",
                   help="accession:residue of the positive control")
    a.add_argument("--control-label", default="DAT")
    a.add_argument("--tsv", default="da/phospho_sites.tsv")
    a.add_argument("--outdir", default="figures/phospho")
    a.add_argument("--font", default=FONT)
    x = sub.add_parser("export", help="the sites in fig3a_phospho_sites.py's table format")
    x.add_argument("--scan", default="derived/mods")
    x.add_argument("--out", default="derived/mods/sites_phospho.tsv")
    x.add_argument("--complete", action="store_true",
                   help="only patients every digest measured")
    i = sub.add_parser("isoforms",
                       help="phosphosites on residues only an isoform has")
    i.add_argument("reports", nargs="+",
                   help="the unmodified search")
    i.add_argument("--scan", default="derived/mods")
    i.add_argument("--fasta", required=True)
    i.add_argument("--tsv", default="derived/da_iso/phospho_isoform_sites.tsv")
    i.add_argument("--regions-tsv", default="derived/da_iso/phospho_region_sites.tsv")
    i.add_argument("--pairs-tsv", default="derived/da_iso/phospho_pair_sites.tsv",
                   help="sites a pair panel can draw")
    i.add_argument("--events", default="data/gencode/events.tsv",
                   help="splice-event typing from prep_splice_events.py")
    i.add_argument("--cache", default="data/uniprot")
    i.add_argument("--outdir", default="figures/phospho")
    i.add_argument("--font", default=FONT)
    for name in ("report", "figures"):
        r = sub.add_parser(name)
        r.add_argument("--scan", default="derived/mods")
        r.add_argument("--fasta", required=True)
        r.add_argument("--cache", default="data/uniprot")
        if name == "figures":
            r.add_argument("--outdir", default="figures/phospho")
            r.add_argument("--example", default="P07197",
                           help="accession of the worked example")
            r.add_argument("--residue", type=int, default=736)
            r.add_argument("--label", default="NEFM")
            r.add_argument("--font", default=FONT)
    args = ap.parse_args(argv)
    if args.cmd == "scan":
        paths = [p for pat in args.reports for p in (sorted(glob.glob(pat)) or [pat])]
        scan(paths, args.fasta, args.outdir, args.precursor_q,
             args.min_run_precursors)
    elif args.cmd == "known":
        ks = known_sites(args.cache, args.force)
        n = sum(len(v) for v in ks.values())
        exp = sum(e == "experimental" for v in ks.values() for e in v.values())
        print(f"  {n:,} phospho MOD_RES on {len(ks):,} entries, "
              f"{exp:,} with experimental evidence")
    elif args.cmd == "report":
        report(args.scan, args.fasta, args.cache)
    elif args.cmd == "da":
        da(args)
    elif args.cmd == "isoforms":
        isoforms(args)
    elif args.cmd == "strip":
        strip(args)
    elif args.cmd == "dots":
        dots(args)
    elif args.cmd == "export":
        rows = load(args.scan)
        if args.complete:
            keep = set(complete_patients(load_runs(args.scan)))
            rows = [r for r in rows if r["sample"] in keep]
        export_sites(sites(rows), args.out)
    elif args.cmd == "figures":
        co = Cohort(load(args.scan), known_sites(args.cache), read_fasta(args.fasta),
                    complete_patients(load_runs(args.scan)))
        panel_gain(co, shallow_runs(args.scan), args.outdir, args.font, "a")
        panel_curated(co, args.outdir, args.font, "b")
        panel_example(co, args.example, args.residue, args.label, args.outdir,
                      args.font, "c")
        panel_agreement(co, args.outdir, args.font)
        panel_residues(co, args.outdir, args.font)
    return 0


if __name__ == "__main__":
    sys.exit(main())
