#!/usr/bin/env python3
"""Phospho precursors placed on the FASTA, and the site table they merge into.

    python3 prep_phospho.py scan 'data/search/*-60min-Phospho.parquet' --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta \\
        --outdir derived/mods
"""

import argparse
import collections
import csv
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fig3a_phospho_sites import parse_mods                  # noqa: E402
from lib_palette import INK, INK_MUTED                      # noqa: E402
from lib_report import sample_of                            # noqa: E402
from lib_fasta import read_fasta                            # noqa: E402

ORDER = ["GluC", "LysC", "Trypsin"]
PHOSPHO = "Phospho"
ISOFORM = re.compile(r"-\d+$")
K = 7                      # shortest peptide the searches allow


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


def scan(paths, fasta, outdir, precursor_q=0.01, min_run_precursors=1000):
    """Write phospho_runs.tsv and phospho_precursors.tsv to `outdir`."""
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


def header(c, x, y, title, sub=None):
    c.text(x, y, title, 11.5, INK, "start", "600")
    if sub:
        c.text(x, y + 15, sub, 9, INK_MUTED, "start")


def save(c, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    p = os.path.join(outdir, name)
    with open(p, "w") as fh:
        fh.write(c.out())
    print(f"  wrote {p}")


def mod_sites_rows(S, mod="Phospho (STY)"):
    """The sites in `fig3a_phospho_sites.load`'s shape, at each site's first key."""
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
    """Write the sites in the table `fig3a_phospho_sites.load` reads."""
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
    x = sub.add_parser("export", help="the sites in fig3a_phospho_sites.py's table format")
    x.add_argument("--scan", default="derived/mods")
    x.add_argument("--out", default="derived/mods/sites_phospho.tsv")
    x.add_argument("--complete", action="store_true",
                   help="only patients every digest measured")
    args = ap.parse_args(argv)
    if args.cmd == "scan":
        paths = [p for pat in args.reports for p in (sorted(glob.glob(pat)) or [pat])]
        scan(paths, args.fasta, args.outdir, args.precursor_q,
             args.min_run_precursors)
    elif args.cmd == "export":
        rows = load(args.scan)
        if args.complete:
            keep = set(complete_patients(load_runs(args.scan)))
            rows = [r for r in rows if r["sample"] in keep]
        export_sites(sites(rows), args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
