#!/usr/bin/env python3
"""The reference's "differential isoform ratios", rebuilt with points coloured by protease.

    python3 extra_reported_groups.py 'data/search/*-60min-Phospho.parquet' --fasta data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta \\
        --metadata data/metadata.xlsx --genes ALDH9A1,GNAI2 \\
        --out figures/supp_ref_facets.svg
"""

import argparse
import collections
import glob
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fig2a_supp1_diagnostic_peptides as dp                                 # noqa: E402
from lib_palette import FONT, display                                # noqa: E402
from lib_report import read_metadata                                 # noqa: E402
from lib_fasta import read_fasta                          # noqa: E402

ISO_SUFFIX = "-"


def is_canonical_group(g):
    """A group naming exactly one accession, and that accession not an isoform."""
    accs = [a.strip() for a in g.split(";") if a.strip()]
    return len(accs) == 1 and ISO_SUFFIX not in accs[0]


def facet(gene, groups, idx, med, cond, case, control, min_group):
    """One gene's facet: the canonical group, then every other reported group."""
    cols = {}
    for g in groups:
        # a set: a peptide seen by two digests counts once
        mem = sorted({(g, p) for d in dp.ORDER for p in idx.get((g, d), {})})
        st = dp.column(mem, idx, med, cond, case, control, min_group)
        if st:
            cols[g] = st
    if len(cols) < 2:
        return None
    base = next((g for g in sorted(cols) if is_canonical_group(g)), None)
    if base is None:
        return None
    rest = sorted(g for g in cols if g != base)
    return {"gene": gene, "canon_acc": base, "canon": cols[base],
            "diag": [(g, cols[g]) for g in rest]}


def digest_split(st):
    """-> {digest: points drawn} for one column."""
    n = collections.Counter(d for _p, d, _v in st["points"])
    return {d: n[d] for d in dp.ORDER}


def patient_contrast(base, other, idx, med, cond, case, control):
    """-> (patient-level Welch p between columns, digests that saw both)."""
    from scipy import stats

    def pts(g):
        out = collections.defaultdict(list)
        for d in dp.ORDER:
            for p, v, _a, _b in dp.peptide_points(
                    [(g, q) for q in idx.get((g, d), {})], d, idx, med,
                    cond, case, control):
                out[d].append(v)
        return out

    pa, pb = pts(base), pts(other)
    shared = [d for d in dp.ORDER if pa[d] and pb[d]]
    if not shared:
        return float("nan"), shared

    diffs = collections.defaultdict(list)
    for d in shared:
        lv = {g: dp.sample_levels([(g, q) for q in idx.get((g, d), {})],
                                  d, idx, med)[0] for g in (base, other)}
        for s in set(lv[base]) & set(lv[other]):
            diffs[s].append(lv[other][s] - lv[base][s])
    avg = {s: statistics.mean(v) for s, v in diffs.items()}
    a = [v for s, v in avg.items() if cond.get(s) == case]
    b = [v for s, v in avg.items() if cond.get(s) == control]
    if len(a) < 2 or len(b) < 2:
        return float("nan"), shared
    return float(stats.ttest_ind(a, b, equal_var=False).pvalue), shared


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reports", nargs="*", default=["data/search/*-60min-Phospho.parquet"])
    ap.add_argument("--fasta", default="data/uniprot_sprot_2024-01-01_HUMAN_ISOFORMS.fasta")
    ap.add_argument("--metadata", default="data/metadata.xlsx")
    ap.add_argument("--genes", default="ALDH9A1,GNAI2")
    ap.add_argument("--case", default="LBD")
    ap.add_argument("--control", default="Control")
    ap.add_argument("--precursor-q", type=float, default=0.01)
    ap.add_argument("--quantity", choices=("raw", "normalised", "ms1-apex"),
                    default="raw",
                    help="precursor quantity column")
    ap.add_argument("--min-group", type=int, default=3)
    ap.add_argument("--letter", default="")
    ap.add_argument("--panel-w", type=float, default=None)
    ap.add_argument("--panel-h", type=float, default=None)
    ap.add_argument("--ylim", default=None)
    ap.add_argument("--no-key", dest="key", action="store_false",
                    help="drop the protease key")
    ap.add_argument("--no-note", dest="note", action="store_false",
                    help="drop the footnote")
    ap.add_argument("--no-subtitle", dest="subtitle", action="store_false",
                    help="drop the subtitle")
    ap.add_argument("--out", default="figures/supp_ref_facets.svg")
    ap.add_argument("--font", default=FONT)
    ap.add_argument("--width", type=float, default=1215.0)
    args = ap.parse_args(argv)

    paths = sorted(p for pat in args.reports for p in glob.glob(pat))
    if not paths:
        sys.exit("no reports matched")
    want = {g.strip() for g in args.genes.split(",") if g.strip()}
    from lib_fasta import gene_map
    seqs = read_fasta(args.fasta)
    genes = gene_map(args.fasta)
    cond = read_metadata(args.metadata)
    idx, med, bygene = dp.read_reports(
        paths, seqs, genes, args.precursor_q, want=want,
        quantity={"raw": "quantity", "normalised": "quantity_normalised",
                  "ms1-apex": "quantity_ms1_apex"}[args.quantity])

    picks = []
    for gene in sorted(want):
        f = facet(gene, bygene.get(gene, ()), idx, med, cond, args.case,
                  args.control, args.min_group)
        if not f:
            print(f"  {gene}: fewer than two reported groups, skipped")
            continue
        base, other = f["canon_acc"], f["diag"][0][0]
        dp.apply_vs_canonical(f["canon"], f["diag"][0][1])
        p_pep = f["diag"][0][1]["p"]
        eff = f["diag"][0][1]["mean"]
        p_pat, shared = patient_contrast(base, other, idx, med, cond,
                                         args.case, args.control)
        f["sub"] = ("Δ %+.2f  ·  p = %s"
                    % (eff, f"{p_pep:.0e}".replace("e-0", "e-")))
        picks.append(f)

        print(f"\n  {gene}")
        for g, st in [(base, f["canon"])] + f["diag"]:
            sp = digest_split(st)
            only = [d for d in dp.ORDER if sp[d]]
            mid = statistics.mean(v for _p, _d, v in st["points"])
            print(f"    {g:<20s} {len(st['points']):>3d} points  "
                  f"mean {mid:+.2f}   "
                  + "  ".join(f"{display(d)} {sp[d]}" for d in dp.ORDER)
                  + ("   <- ONE DIGEST ONLY" if len(only) == 1 else ""))
        print(f"    digests measuring both columns: "
              f"{', '.join(display(d) for d in shared) if shared else 'NONE'}")
        print(f"    Welch between columns, over peptides {p_pep:.1e}   "
              + ("patient-level p NOT COMPUTABLE within a digest" if not shared
                 else f"patient-level p {p_pat:.2f}"))

    if not picks:
        sys.exit("nothing to draw")
    ylim = None
    if args.ylim:
        lo, hi = args.ylim.split(",")
        ylim = (float(lo), float(hi))
    note = ("Columns are reported protein groups and points are peptides, as in "
            "the reference; each point is coloured by the protease that produced "
            "it. A column in one colour was measured by one digest, so the "
            "contrast between such columns is between separately searched "
            "experiments rather than between proteoforms.") if args.note else ""
    dp.panel(picks, args.out, args.font, args.letter, args.width, key=args.key,
             ylim=ylim, subtitle=args.subtitle, panel_w=args.panel_w,
             panel_h=args.panel_h, hue=True, case_label=args.case,
             control_label=args.control, note=note,
             label_fn=lambda who: who,
             title="Differential isoform ratios, as the reference builds them",
             head_text=f"one point per peptide per protease · log2 {args.case}"
                       f" / {args.control} · columns are reported protein "
                       f"groups, tested against each other")
    return 0


if __name__ == "__main__":
    sys.exit(main())
