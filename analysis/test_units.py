#!/usr/bin/env python3
"""Assertions for the analysis modules on hand-built fixtures.

    python3 test_units.py

Every case uses a report small enough that the right answer is countable by
hand, so a wrong count is a failure and not a judgement call.
"""

import csv
import os
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
COUNTS = os.path.join(HERE, "prep_counts.py")
FIGURE = os.path.join(HERE, "fig1b_depth.py")
WORKFLOW = os.path.join(HERE, "fig1a_workflow.py")
COVERAGE = os.path.join(HERE, "fig1cd_coverage.py")

failures = []


def check(name, got, want):
    if got == want:
        print(f"  ok   {name}")
    else:
        print(f"  FAIL {name}: got {got!r}, want {want!r}")
        failures.append(name)


def run(args, expect_rc=0):
    p = subprocess.run([sys.executable] + args, capture_output=True, text=True)
    if p.returncode != expect_rc:
        print(f"  FAIL command rc={p.returncode} (wanted {expect_rc}): {' '.join(args)}")
        print("  ---- stdout ----\n" + p.stdout[-2500:])
        print("  ---- stderr ----\n" + p.stderr[-2500:])
        failures.append(" ".join(args))
    return p


def write_report(path, header, rows, delim="\t"):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter=delim, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def load_counts(path):
    with open(path, newline="") as fh:
        return {(r["sample"], r["protease"], r["unit"]): int(r["count"])
                for r in csv.DictReader(fh)}


def cut(seq, residues, block_p=False):
    """Fully cleaved 7-30-residue peptides of `seq`, C-terminal peptide excluded."""
    sites = [0] + [i + 1 for i, a in enumerate(seq[:-1])
                   if a in residues and not (block_p and seq[i + 1] == "P")]
    return [seq[a:b] for a, b in zip(sites, sites[1:]) if 7 <= b - a <= 30]


def workflow_fasta(tmp):
    """-> (path, accession, sequence) of a pseudorandom 300-residue entry, plus a second gene
    that shares the first Lys-C peptide."""
    import random
    seq = "M" + "".join(random.Random(11).choices("ACDEFGHIKLMNPQRSTVWY", k=299))
    path = os.path.join(tmp, "workflow.fasta")
    with open(path, "w") as fh:
        fh.write(f">sp|P99999|TEST_HUMAN test GN=TEST\n{seq}\n")
        fh.write(f">sp|P99998|OTHER_HUMAN other GN=OTHER\nMW{cut(seq, 'K')[0]}WW\n")
    return path, "P99999", seq


def workflow_reports(tmp, seq):
    """-> one Spectronaut-like parquet per protease, peptides cut from `seq`; the last Glu-C
    peptide fails the precursor q-value."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    paths = []
    for k, (tag, residues, block_p) in enumerate((("G", "E", False), ("L", "K", False),
                                                   ("T", "KR", True))):
        peps = cut(seq, residues, block_p)
        q = [0.001] * len(peps)
        if tag == "G":
            q[-1] = 0.5
        p = os.path.join(tmp, f"workflow_{tag}.parquet")
        pq.write_table(pa.table({
            "R_FileName": [f"2026-01-01_CF_{1000 + k}_{tag}_test_A1_1_{9000 + k}"] * len(peps),
            "PG_ProteinGroups": ["P99999"] * len(peps),
            "PEP_StrippedSequence": peps,
            "EG_Qvalue": q,
            "FG_Quantity": [10.0] * len(peps),
            "EG_UsedForPeptideQuantity": [True] * len(peps),
            "EG_ApexRT": [1.0] * len(peps),
        }), p)
        paths.append(p)
    return paths


# ---------------------------------------------------------------------------
# 1. Filtering, union dedup and isoform collapsing on one tiny fixture.
# ---------------------------------------------------------------------------
def case_semantics(tmp):
    hdr = ["R.FileName", "R.Condition", "R.Replicate", "PG.ProteinGroups",
           "PG.Qvalue", "PEP.StrippedSequence", "EG.ModifiedSequence",
           "EG.Qvalue", "EG.IsDecoy", "FG.Charge"]
    rows = [
        # run                 cond      rep    group        pgq      pep         mod           egq     decoy  z
        ["r_S1_Trypsin", "Trypsin", "S1", "P49189",   "0.001", "PEPTIDEA", "_PEPTIDEA_", "0.001", "False", "2"],
        ["r_S1_Trypsin", "Trypsin", "S1", "P49189",   "0.001", "PEPTIDEA", "_PEPTIDEA_", "0.001", "False", "3"],
        ["r_S1_Trypsin", "Trypsin", "S1", "P49189-3", "0.001", "PEPTIDEB", "_PEPTIDEB_", "0.001", "False", "2"],
        ["r_S1_Trypsin", "Trypsin", "S1", "P04899",   "0.001", "DECOYPEP", "_DECOYPEP_", "0.001", "True",  "2"],
        ["r_S1_Trypsin", "Trypsin", "S1", "P04899",   "0.001", "BADQPEP",  "_BADQPEP_",  "0.05",  "False", "2"],
        ["r_S1_Trypsin", "Trypsin", "S1", "P04899",   "0.20",  "BADPGPEP", "_BADPGPEP_", "0.001", "False", "2"],
        ["r_S1_GluC",    "GluC",    "S1", "P49189",   "0.001", "PEPTIDEA", "_PEPTIDEA_", "0.001", "False", "2"],
        ["r_S1_GluC",    "GluC",    "S1", "P12345",   "0.001", "PEPTIDEC", "_PEPTIDEC_", "0.001", "False", "2"],
    ]
    rep = os.path.join(tmp, "sem.tsv")
    write_report(rep, hdr, rows)
    out = os.path.join(tmp, "sem")
    p = run([COUNTS, rep, "--out", out])
    c = load_counts(os.path.join(out, "counts.csv"))

    # each unit has its own q cutoff: BADPGPEP counts as a precursor, BADQPEP as a protein
    # Trypsin precursors: PEPTIDEA z2, PEPTIDEA z3, PEPTIDEB, BADPGPEP = 4
    check("trypsin precursors", c[("S1", "Trypsin", "precursors")], 4)
    check("trypsin peptides", c[("S1", "Trypsin", "peptides")], 3)
    check("gluc precursors", c[("S1", "GluC", "precursors")], 2)
    check("gluc peptides", c[("S1", "GluC", "peptides")], 2)
    # union precursors: PEPTIDEA/2 is in both digests -> 4 + 2 - 1 = 5
    check("union precursors dedups", c[("S1", "All", "precursors")], 5)
    # union peptides: {A,B,BADPGPEP} u {A,C} = 4
    check("union peptides dedups", c[("S1", "All", "peptides")], 4)
    # P49189, P49189-3 and P04899 (via BADQPEP) = 3 isoform groups, 2 canonical
    check("trypsin isoform groups", c[("S1", "Trypsin", "protein_isoform_groups")], 3)
    check("trypsin canonical", c[("S1", "Trypsin", "proteins_canonical")], 2)
    check("union canonical", c[("S1", "All", "proteins_canonical")], 3)
    check("drop counts reported", "dropped decoy" in p.stdout, True)


# ---------------------------------------------------------------------------
# 2. Header variants: other capitalisation, derived stripped sequence, precursor id.
# ---------------------------------------------------------------------------
def case_header_variants(tmp):
    hdr = ["R.Raw File Name", "R.Condition", "PG.ProteinAccessions", "PG.QValue",
           "EG.PrecursorId", "EG.QValue"]
    rows = [
        ["x_S9_LysC_y", "LysC", "sp|P11111|AAA_HUMAN", "0.002", "_AAAK[Acetyl]LM_.2", "0.002"],
        ["x_S9_LysC_y", "LysC", "sp|P11111|AAA_HUMAN", "0.002", "_AAAKLM_.3", "0.002"],
    ]
    rep = os.path.join(tmp, "hdr.csv")
    write_report(rep, hdr, rows, delim=",")
    out = os.path.join(tmp, "hdr")
    run([COUNTS, rep, "--out", out])
    c = load_counts(os.path.join(out, "counts.csv"))
    check("variant precursors", c[("x_S9_y", "LysC", "precursors")], 2)
    # both rows collapse to AAAKLM once modifications and charge are stripped
    check("variant peptides (derived strip)", c[("x_S9_y", "LysC", "peptides")], 1)
    check("variant sp| accession parsed", c[("x_S9_y", "LysC", "proteins_canonical")], 1)


# ---------------------------------------------------------------------------
# 3. An unresolvable protease must stop, not guess.
# ---------------------------------------------------------------------------
def case_unresolved(tmp):
    hdr = ["R.FileName", "PG.ProteinGroups", "EG.ModifiedSequence", "FG.Charge"]
    rows = [["run_alpha_01", "P1", "_AAA_", "2"], ["run_beta_02", "P1", "_AAB_", "2"]]
    rep = os.path.join(tmp, "unres.tsv")
    write_report(rep, hdr, rows)
    p = run([COUNTS, rep, "--out", os.path.join(tmp, "unres")], expect_rc=1)
    check("unresolved protease refuses", "UNRESOLVED" in p.stdout, True)


# ---------------------------------------------------------------------------
# 4. Missing required columns must fail loudly with the header echoed.
# ---------------------------------------------------------------------------
def case_missing_columns(tmp):
    rep = os.path.join(tmp, "bad.tsv")
    write_report(rep, ["Foo", "Bar"], [["1", "2"]])
    p = run([COUNTS, rep, "--out", os.path.join(tmp, "bad")], expect_rc=1)
    check("missing columns names header", "header seen" in p.stderr, True)


# ---------------------------------------------------------------------------
# 5. Figure geometry: valid XML, every point inside its panel, outlier drawn.
# ---------------------------------------------------------------------------
def case_figure(tmp):
    csv_path = os.path.join(tmp, "counts.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "protease", "unit", "count"])
        for i in range(1, 13):
            base = 100 * i
            w.writerow([f"S{i:02d}", "Trypsin", "peptides", base + 5000])
            w.writerow([f"S{i:02d}", "GluC", "peptides", base + 2000])
            w.writerow([f"S{i:02d}", "All", "peptides", base + 6000])
        # a failed run far below the rest -- must still be plotted
        w.writerow(["S13", "GluC", "peptides", 40])
    svg = os.path.join(tmp, "f.svg")
    p = run([FIGURE, csv_path, "--out", svg, "--units", "peptides"])
    check("figure reports the outlier", "outliers" in p.stdout, True)

    root = ET.parse(svg).getroot()          # raises if the SVG is malformed
    ns = "{http://www.w3.org/2000/svg}"
    H = float(root.get("height"))
    circles = root.findall(f"{ns}circle")
    check("all 37 points drawn", len(circles), 37)
    ys = [float(c.get("cy")) for c in circles]
    check("no point above the canvas", all(0 < y < H for y in ys), True)

    # panel band: strip bottom .. axis line. Nothing plotted may sit above it.
    rects = [r for r in root.findall(f"{ns}rect") if r.get("fill") == "#f0efec"]
    check("one facet strip", len(rects), 1)
    strip_bottom = float(rects[0].get("y")) + float(rects[0].get("height"))
    check("points below the facet strip", min(ys) > strip_bottom, True)
    grid = [float(l.get("y1")) for l in root.findall(f"{ns}line")
            if l.get("stroke") == "#e1e0d9"]
    check("gridlines below the facet strip", all(y > strip_bottom for y in grid), True)
    labels = [float(t.get("y")) for t in root.findall(f"{ns}text")
              if t.get("text-anchor") == "end"]
    check("tick labels below the facet strip", all(y > strip_bottom for y in labels), True)


# ---------------------------------------------------------------------------
# 6. The workflow schematic: valid XML, on the canvas, digest obeys specificity.
# ---------------------------------------------------------------------------
def case_workflow(tmp):
    import lib_palette as palette
    from lib_palette import INK_MUTED
    fa, acc, seq = workflow_fasta(tmp)
    reports = workflow_reports(tmp, seq)
    svg = os.path.join(tmp, "fig1a.svg")
    run([WORKFLOW] + reports + ["--out", svg, "--fasta", fa, "--protein", acc,
                                "--acquisition", "65 min · 60 windows"])

    root = ET.parse(svg).getroot()
    ns = "{http://www.w3.org/2000/svg}"
    # iter, not findall: the sequence band is drawn inside a translated <g>
    W, H = float(root.get("width")), float(root.get("height"))

    # every protease keeps the colour the count panels give it
    colours = palette.assign(["Glu-C", "Lys-C", "Trypsin", "All"])
    check("workflow colours match counts",
          [colours[k] for k in ("Glu-C", "Lys-C", "Trypsin")],
          [palette.assign(["GluC", "LysC", "Trypsin", "All"])[k]
           for k in ("GluC", "LysC", "Trypsin")])
    for name in ("Glu-C", "Lys-C", "Trypsin"):
        used = [r for r in root.iter(f"{ns}rect") if r.get("fill") == colours[name]]
        check(f"{name} coverage drawn", len(used) >= 2, True)

    texts = [t.text for t in root.iter(f"{ns}text")]
    for want in ("Glu-C", "Lys-C", "Trypsin", "dia-PASEF", "directDIA",
                 "All proteases", "All"):
        check(f"label present: {want}", want in texts, True)
    # the union row is bars where any specificity reaches, not a shaded strip
    check("the union row is bars, not a strip", "covered by" in texts, False)
    union_bars = [r for r in root.iter(f"{ns}rect")
                  if r.get("fill") == palette.UNION
                  and float(r.get("height")) == 15.0]
    check("union drawn as several bars with gaps", len(union_bars) >= 2, True)
    check("the union is named once", texts.count("All proteases"), 1)
    check("no stale union names", ("Union" in texts, "All three" in texts),
          (False, False))
    check("nothing restates what the lanes already show", "3 digests" in texts,
          False)
    check("no standalone colour key",
          [t for t in texts if t in ("protease", "used in every panel")], [])
    check("no depth key", [t for t in texts if t == "digests"], [])
    # the chips carry the colour-to-protease mapping
    for name in ("Glu-C", "Lys-C", "Trypsin"):
        check(f"{name} is still named in the workflow", texts.count(name) >= 2,
              True)

    # a long label is allowed only as a data statement, i.e. with a digit in it
    prose = [t for t in texts if t and len(t) > 24 and not any(ch.isdigit()
                                                              for ch in t)]
    check("no prose on the panel", prose, [])
    check("no label is absurdly long",
          max(len(t or "") for t in texts) <= 60, True)
    check("no full stops in labels",
          [t for t in texts if t and t.rstrip().endswith(".")], [])

    # only the panel letter may be a lone glyph
    singles = [t for t in texts if t and len(t) == 1 and t.isalpha()]
    check("no per-residue sequence letters", len(singles) <= 1, True)
    check("no in-silico footnote", [t for t in texts if t and "in silico" in t], [])

    # no filled triangle carries a protease colour; only the arrowheads are filled
    tri = [p_ for p_ in root.findall(f"{ns}path")
           if p_.get("fill") in {colours[k] for k in ("Glu-C", "Lys-C", "Trypsin")}]
    check("no cut-site markers", len(tri), 0)
    heads = [p_ for p_ in root.findall(f"{ns}path") if p_.get("fill") == INK_MUTED]
    check("workflow arrows keep their heads", len(heads) >= 4, True)

    # nothing may spill off the canvas
    xs, ys = [], []
    for el in root.iter():
        for a, bucket in (("x", xs), ("y", ys), ("cx", xs), ("cy", ys),
                          ("x1", xs), ("x2", xs), ("y1", ys), ("y2", ys)):
            if el.get(a) is not None:
                bucket.append(float(el.get(a)))
    check("no geometry left of canvas", min(xs) >= 0, True)
    check("no geometry right of canvas", max(xs) <= W, True)
    check("no geometry above canvas", min(ys) >= 0, True)
    check("no geometry below canvas", max(ys) <= H, True)

    # coverage comes from each report's own peptides, gene-specific and within the q cut
    sys.path.insert(0, HERE)
    import fig1a_workflow as wf
    seen = wf.observed(reports, fa, acc)
    check("each report's coverage lands under its protease", sorted(seen),
          ["Glu-C", "Lys-C", "Trypsin"])
    for name, residues in (("Glu-C", "E"), ("Lys-C", "K"), ("Trypsin", "KR")):
        check(f"{name} spans end at its own cleavage sites",
              all(seq[e - 1] in residues for s, e in seen[name]), True)
    gluc = cut(seq, "E")
    check("a peptide over the q cut is not drawn",
          (seq.find(gluc[-1]), seq.find(gluc[-1]) + len(gluc[-1])) in seen["Glu-C"], False)
    lysc = cut(seq, "K")
    check("a peptide another gene also contains is not drawn",
          (seq.find(lysc[0]), seq.find(lysc[0]) + len(lysc[0])) in seen["Lys-C"], False)
    check("the gene-specific Lys-C peptides are drawn", len(seen["Lys-C"]), len(lysc) - 1)
    union = wf.coverage_mask([sp for v in seen.values() for sp in v], len(seq))
    check("the pooled lane is labelled with its observed coverage",
          f"{sum(union) / len(seq):.0%}" in texts and "observed" in texts, True)


# ---------------------------------------------------------------------------
# 7. fig1cd_coverage: FASTA parsing, coverage maths and the rank curves.
# ---------------------------------------------------------------------------
def case_coverage(tmp):
    import fig1cd_coverage as pc
    from lib_fasta import read_fasta

    fa = os.path.join(tmp, "db.fasta")
    with open(fa, "w") as fh:
        fh.write(">sp|P00001|AAA_HUMAN alpha\nMKAAAWWWKDDEFGHIKLM\n")   # 19 aa
        fh.write(">P00002\nMKAAAWWWKDDEFGHIKLMNPQRST\n")                # 25 aa
        fh.write(">sp|P00001-2|AAA_HUMAN iso\nMKAAAWWWKHIKLM\n")        # del 10-14
    seqs = read_fasta(fa)
    check("fasta reads sp| accessions", sorted(seqs), ["P00001", "P00001-2", "P00002"])
    check("fasta reads bare header", len(seqs["P00002"]), 25)

    cov, found = pc.mark("MKAAAWWWKDDEFGHIKLM", {"AAAWWW", "HIKLM", "ZZZ"})
    check("coverage marks residues", sum(cov), 11)
    check("coverage counts found peptides", found, 2)

    hdr = ["R.FileName", "PG.ProteinGroups", "PEP.StrippedSequence",
           "EG.ModifiedSequence", "FG.Charge", "EG.Qvalue"]
    rows = []
    for runname, pep in (("r_Trypsin", "AAAWWWK"), ("r_Trypsin", "DDEFGHIK"),
                         ("r_GluC", "MKAAAWWWKDDE"), ("r_LysC", "HIKLM")):
        rows.append([runname, "P00001", pep, f"_{pep}_", "2", "0.001"])
    rows.append(["r_LysC", "P00001-2", "WWWKHIK", "_WWWKHIK_", "2", "0.001"])
    rows.append(["r_GluC", "P00001", "NOTINSEQ", "_NOTINSEQ_", "2", "0.001"])
    # one report per enzyme; the protease comes from each report's filename
    reps = []
    for enz in ("Trypsin", "LysC", "GluC"):
        path = os.path.join(tmp, f"PD_{enz}_Report.tsv")
        write_report(path, hdr, [r for r in rows if r[0] == f"r_{enz}"])
        reps.append(path)

    ns = "{http://www.w3.org/2000/svg}"
    svg = os.path.join(tmp, "rank.svg")
    p = run([COVERAGE] + reps + ["--fasta", fa, "--rank", "--out", svg])
    root = ET.parse(svg).getroot()
    paths = [e for e in root.findall(f"{ns}path") if e.get("fill") == "none"]
    check("one curve per protease plus union", len(paths), 4)
    check("series named per enzyme",
          all(e in p.stdout for e in ("Trypsin", "LysC", "GluC", "All")), True)
    check("rank mode reports how many reach 50%", ">=50%" in p.stdout, True)
    check("rank mode reports the peptide its group's sequence lacks",
          "1 peptide-group pairs did not occur" in p.stdout, True)


# ---------------------------------------------------------------------------
# 8. Shared-y groups: one axis and one set of tick labels per group, one letter.
# ---------------------------------------------------------------------------
def case_figure_groups(tmp):
    csv_path = os.path.join(tmp, "grouped.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "protease", "unit", "count"])
        for i in range(1, 9):
            # two units in the thousands, two in the millions: if the groups
            # leaked into each other the small pair would be flattened
            w.writerow([f"S{i}", "Trypsin", "precursors", 1_000_000 + i * 1000])
            w.writerow([f"S{i}", "Trypsin", "peptides", 900_000 + i * 1000])
            w.writerow([f"S{i}", "Trypsin", "proteins_canonical", 5_000 + i])
            w.writerow([f"S{i}", "Trypsin", "protein_isoform_groups", 6_000 + i])
    svg = os.path.join(tmp, "g.svg")
    run([FIGURE, csv_path, "--out", svg, "--letter", "b", "--units",
         "precursors,peptides|proteins_canonical,protein_isoform_groups"])
    root = ET.parse(svg).getroot()
    ns = "{http://www.w3.org/2000/svg}"

    letters = [t.text for t in root.findall(f"{ns}text")
               if t.text and len(t.text) == 1 and t.text.isalpha()]
    check("one figure-level letter", letters, ["b"])

    # tick labels are right-anchored and drawn once per y-group
    tick_x = sorted({float(t.get("x")) for t in root.findall(f"{ns}text")
                     if t.get("text-anchor") == "end"})
    check("one tick-label column per y-group", len(tick_x), 2)

    # panels in a group share gridline positions; panels across groups do not
    strips = sorted((float(r.get("x")) for r in root.findall(f"{ns}rect")
                     if r.get("fill") == "#f0efec"))
    check("four facet strips", len(strips), 4)
    grid_y = {}
    for ln in root.findall(f"{ns}line"):
        if ln.get("stroke") == "#e1e0d9":
            grid_y.setdefault(round(float(ln.get("x1"))), set()).add(
                round(float(ln.get("y1")), 1))
    cols = [grid_y[k] for k in sorted(grid_y)]
    check("panels 1 and 2 share a y scale", cols[0] == cols[1], True)
    check("panels 3 and 4 share a y scale", cols[2] == cols[3], True)
    check("the two groups differ", cols[0] != cols[2], True)

    # panel order follows --units, canonical before protein groups
    titles = [t.text for t in root.findall(f"{ns}text")
              if t.text in ("Precursors", "Peptides", "Protein (canonical)",
                            "Protein groups")]
    check("panel order follows --units", titles,
          ["Precursors", "Peptides", "Protein (canonical)", "Protein groups"])


# ---------------------------------------------------------------------------
# 10. Sequence-level event decomposition and its merge threshold.
# ---------------------------------------------------------------------------
def case_events(tmp):
    sys.path.insert(0, HERE)
    import fig2a_isoform_strip as ie

    # pseudorandom and aperiodic, so every edit has one best alignment
    import random
    A = "".join(random.Random(11).choices("ACDEFGHIKLMNPQRSTVWY", k=200))
    check("fixture has no repeated 12-mer",
          max(A.count(A[i:i + 12]) for i in range(len(A) - 12)), 1)
    check("identical sequences give no events", ie.events(A, A), [])

    # one clean internal deletion of ten residues
    delA = A[:40] + A[50:]
    check("internal deletion is one event", len(ie.events(A, delA)), 1)
    check("internal deletion is typed", ie.events(A, delA)[0][0], "delete")
    check("deletion coordinates are canonical",
          ie.events(A, delA)[0][1:3], (40, 50))

    # truncations at each end, and an insertion
    check("N-terminal truncation", ie.events(A, A[15:]), [("delete", 0, 15, 0, 0)])
    check("C-terminal truncation", ie.events(A, A[:-15]),
          [("delete", 185, 200, 185, 185)])
    check("internal insertion",
          [e[0] for e in ie.events(A, A[:40] + "WWWWWWWW" + A[40:])], ["insert"])

    # two changes far apart stay two; the same two within min_equal merge
    two = A[:20] + A[30:60] + A[70:]
    check("two distant changes stay separate", len(ie.events(A, two)), 2)
    near = A[:40] + A[43:47] + A[50:]        # 3 matching residues between them
    check("changes closer than min_equal merge",
          len(ie.events(A, near, min_equal=5)), 1)
    check("changes exactly min_equal apart stay separate",
          len(ie.events(A, near, min_equal=3)), 2)
    # both sides of the merged span are non-empty, so it is a replacement
    check("merging two deletions gives a replacement",
          ie.events(A, near, min_equal=5)[0][0], "replace")
    # flank trimming can move an endpoint by a residue, so assert containment
    cs, ce = ie.events(A, near, min_equal=5)[0][1:3]
    check("the merged span covers both changes", (cs <= 40, ce >= 47),
          (True, True))
    check("a lower threshold does not merge them",
          len(ie.events(A, near, min_equal=2)), 2)
    # align() is the threshold-independent half; merge() must reproduce events()
    for me in (0, 2, 3, 5, 12):
        check(f"align+merge == events at min_equal={me}",
              ie.merge(*ie.align(A, near), min_equal=me),
              ie.events(A, near, min_equal=me))
    for seq in (A[15:], A[:-15], A[:40] + "WWWWWWWW" + A[40:], A, delA):
        check("align+merge == events on " + seq[:6],
              ie.merge(*ie.align(A, seq)), ie.events(A, seq))

    check("five separated changes are five events",
          len(ie.events(A, A[:10] + A[12:25] + A[27:40] + A[42:55] + A[57:70]
                        + A[72:])), 5)
    check("MIN_EQUAL is 3", ie.MIN_EQUAL, 3)


# ---------------------------------------------------------------------------
# 11. The heatmap scale: a neutral midpoint, a red arm, and black for none.
# ---------------------------------------------------------------------------
def case_diverging(tmp):
    sys.path.insert(0, HERE)
    import lib_palette as palette
    from fig3bc_phospho_atlas import NONE_FILL, STEPS, fill

    check("midpoint is neutral, not a hue", max(
        int(palette.DIVERGING_MID[i:i + 2], 16) for i in (1, 3, 5)) - min(
        int(palette.DIVERGING_MID[i:i + 2], 16) for i in (1, 3, 5)) <= 8, True)
    check("no PEP means no measurement, not the lowest bin", fill(None), NONE_FILL)
    check("the no-data colour is not on the scale", NONE_FILL in STEPS, False)

    # OKLab L, which the arm was built against, not WCAG relative luminance
    def oklab_l(h):
        def lin(c):
            c = int(c, 16) / 255
            return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
        r, g, b = (lin(h[i:i + 2]) for i in (1, 3, 5))
        l = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b) ** (1 / 3)
        m = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b) ** (1 / 3)
        s_ = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b) ** (1 / 3)
        return 0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s_

    ls = [oklab_l(c) for c in palette.DIVERGING_HIGH]
    check("high arm darkens monotonically",
          all(a > b for a, b in zip(ls, ls[1:])), True)
    check("high arm steps clear the dL floor",
          all(a - b >= 0.06 for a, b in zip(ls, ls[1:])), True)


# ---------------------------------------------------------------------------
# 14. Compact number format, dropping the x-axis furniture, a pinned canvas.
# ---------------------------------------------------------------------------
def case_compact(tmp):
    sys.path.insert(0, HERE)
    from fig1b_depth import TEXT_SCALE, fmt_compact as f
    from lib_palette import AXIS, TEXT_BOOST

    check("hundreds of thousands lose the decimal", f(376550), "377K")
    check("tens of thousands keep one", f(14795), "14.8K")
    check("thousands keep one", f(6404), "6.4K")
    check("a trailing .0 is dropped on axes", f(15001, exact=True), "15K")
    check("axis steps read cleanly",
          [f(v, exact=True) for v in (0, 2500, 100000, 400000)],
          ["0", "2.5K", "100K", "400K"])
    check("millions switch suffix",
          (f(1000000, exact=True), f(1250000)), ("1M", "1.25M"))
    # axis ticks drop a trailing zero; rounded data labels keep it
    check("exact drops a trailing zero", f(5000, exact=True), "5K")
    check("rounded keeps it", f(4023), "4.0K")
    check("rounded never claims more than it has",
          (f(15001), f(9999)), ("15.0K", "10.0K"))
    check("small values are untouched", f(842), "842")
    check("negatives keep their sign", f(-14795), "-14.8K")
    # the format must never widen the label it replaces
    for v in (376550, 14795, 6404, 842, 1250000):
        check(f"compact is no wider than full for {v}",
              len(f(v)) <= len(f"{v:,.0f}"), True)

    rows = []
    for s_i, samp in enumerate(("S1", "S2", "S3")):
        for prot in ("Trypsin", "GluC"):
            for k in range(4 + s_i):
                rows.append([f"r_{samp}_{prot}", prot, samp, f"P{k:05d}", "0.001",
                             f"PEP{k}{prot}", f"_PEP{k}{prot}_", "0.001", "2"])
    rep = os.path.join(tmp, "xt.tsv")
    write_report(rep, ["R.FileName", "R.Condition", "R.Replicate",
                       "PG.ProteinGroups", "PG.Qvalue", "PEP.StrippedSequence",
                       "EG.ModifiedSequence", "EG.Qvalue", "FG.Charge"], rows)
    out = os.path.join(tmp, "xt")
    run([COUNTS, rep, "--out", out])
    counts = os.path.join(out, "counts.csv")
    ns = "{http://www.w3.org/2000/svg}"

    def texts(svg):
        return [(e.text or "").strip() for e in ET.parse(svg).getroot().iter(f"{ns}text")]

    def axis_y(svg):
        return max(float(ln.get("y1")) for ln in ET.parse(svg).getroot().iter(f"{ns}line")
                   if ln.get("stroke") == AXIS)

    on = os.path.join(tmp, "xticks_on.svg")
    run([FIGURE, counts, "--out", on, "--units", "peptides"])
    t = texts(on)
    check("default keeps the category labels", "Trypsin" in t, True)
    check("default keeps the x title", "Protease" in t, True)

    off = os.path.join(tmp, "xticks_off.svg")
    run([FIGURE, counts, "--out", off, "--units", "peptides", "--no-xticks"])
    t = texts(off)
    check("--no-xticks drops the category labels", "Trypsin" in t, False)
    check("--no-xticks drops the x title too", "Protease" in t, False)
    check("--no-xticks keeps the y title", "Count / sample" in t, True)
    # the panel grows into the bottom margin rather than leaving a white band
    check("--no-xticks gives the bottom margin to the panel",
          axis_y(off) > axis_y(on) + 30, True)

    # the canvas is the requested size, and the median labels clear the strip in it
    med = f"{9.5 * TEXT_SCALE * TEXT_BOOST:.1f}"
    for w, h in (("700", "400"), ("1215", "422")):
        fp = os.path.join(tmp, f"fix_{w}.svg")
        run([FIGURE, counts, "--out", fp, "--units", "peptides",
             "--width", w, "--height", h])
        r = ET.parse(fp).getroot()
        check(f"--width/--height pin the canvas at {w}x{h}",
              (r.get("width"), r.get("height")), (w, h))
        doc = open(fp).read()
        sbot = max(float(y) + float(hh) for _x, y, _w2, hh in
                   re.findall(r'<rect x="([\d.]+)" y="([\d.]+)" width="([\d.]+)" '
                              r'height="([\d.]+)" fill="#f0efec"', doc))
        ys = [float(y) for y in re.findall(
            rf'<text x="[\d.]+" y="([\d.-]+)"\s*font-size="{med}" '
            rf'text-anchor="middle"', doc)]
        check(f"median labels clear the strip at {w}x{h}",
              bool(ys) and all(y - float(med) * 0.78 > sbot for y in ys), True)
    p = run([FIGURE, counts, "--out", os.path.join(tmp, "tiny.svg"), "--units", "peptides",
             "--height", "120"], expect_rc=1)
    check("a canvas too small for the panel is refused", "enlarge the figure" in p.stderr, True)


# ---------------------------------------------------------------------------
# 15. Text must stay inside the canvas, measured on the rendered glyphs.
# ---------------------------------------------------------------------------
def case_text_bounds(tmp):
    sys.path.insert(0, HERE)
    try:
        import pymupdf                                    # noqa: F401
    except ImportError:
        print("  skip  pymupdf not installed")
        return

    def margins(path):
        import pymupdf
        doc = pymupdf.open(path)
        pg = doc[0]
        src = open(path).read()
        vb = [float(v) for v in
              re.search(r'viewBox="([^"]+)"', src).group(1).split()]
        W, H = vb[2], vb[3]
        sx, sy = W / pg.rect.width, H / pg.rect.height
        box = [(sp["bbox"][0] * sx, sp["bbox"][2] * sx,
                sp["bbox"][1] * sy, sp["bbox"][3] * sy)
               for blk in pg.get_text("dict")["blocks"]
               for ln in blk.get("lines", []) for sp in ln["spans"]]
        if not box:
            return None
        return (min(b[0] for b in box), W - max(b[1] for b in box),
                min(b[2] for b in box), H - max(b[3] for b in box))

    fa, acc, seq = workflow_fasta(tmp)
    wf = os.path.join(tmp, "wf.svg")
    run([WORKFLOW] + workflow_reports(tmp, seq) +
        ["--out", wf, "--letter", "a", "--fasta", fa, "--protein", acc,
         "--label", "MCCC1", "--acquisition", "65 min · 60 windows"])
    m = margins(wf)
    check("workflow text stays inside the canvas",
          m is not None and min(m) >= 0, True)
    # a renderer substituting a wider font must not push it out either
    check("workflow keeps a margin a wider font cannot eat",
          m is not None and min(m) >= 8, True)

    # a label must also stay inside the box it sits in
    def box_overflow(path):
        import pymupdf
        src = open(path).read()
        vb = [float(v) for v in
              re.search(r'viewBox="([^"]+)"', src).group(1).split()]
        W, H = vb[2], vb[3]
        # the boxes are the workflow's; the translated sequence band is skipped
        flat = re.sub(r'<g transform="translate\([^"]*\)">.*?</g>', "", src,
                      flags=re.S)
        boxes = [tuple(map(float, m)) for m in
                 re.findall(r'<rect x="([\d.]+)" y="([\d.]+)" width="([\d.]+)" '
                            r'height="([\d.]+)"[^>]*rx="[1-9]', flat)]
        pg = pymupdf.open(path)[0]
        sx, sy = W / pg.rect.width, H / pg.rect.height
        spans = [((sp["bbox"][0] + sp["bbox"][2]) / 2 * sx, sp["bbox"][0] * sx,
                  sp["bbox"][2] * sx,
                  (sp["bbox"][1] + sp["bbox"][3]) / 2 * sy, sp["text"])
                 for blk in pg.get_text("dict")["blocks"]
                 for ln in blk.get("lines", []) for sp in ln["spans"]]
        out = []
        for bx, by, bw, bh in boxes:
            for cx, x0, x1, cy, t in spans:
                # contained, not merely level with it
                if bx < cx < bx + bw and by <= cy <= by + bh:
                    if min(x0 - bx, bx + bw - x1) < 0:
                        out.append(t)
        return out

    check("no label overflows its box", box_overflow(wf), [])

    hdr = ["R.FileName", "R.Condition", "R.Replicate", "PG.ProteinGroups",
           "PG.Qvalue", "PEP.StrippedSequence", "EG.ModifiedSequence",
           "EG.Qvalue", "FG.Charge"]
    rows = [[f"r_{s_}_{p}", p, s_, f"P{k:05d}", "0.001", f"PEP{k}{p}",
             f"_PEP{k}{p}_", "0.001", "2"]
            for s_ in ("S1", "S2", "S3") for p in ("Trypsin", "GluC")
            for k in range(6)]
    rep = os.path.join(tmp, "tb.tsv")
    write_report(rep, hdr, rows)
    out = os.path.join(tmp, "tb")
    run([COUNTS, rep, "--out", out])
    counts = os.path.join(out, "counts.csv")
    for w, h in (("900", "380"), ("1215", "422")):
        f = os.path.join(tmp, f"tb_{w}.svg")
        run([FIGURE, counts, "--out", f, "--units", "peptides,precursors",
             "--letter", "b", "--no-xticks", "--width", w, "--height", h])
        m = margins(f)
        check(f"count panel text stays inside at {w}x{h}",
              m is not None and min(m) >= 8, True)


# ---------------------------------------------------------------------------
# 16. One spelling per protease across the whole figure.
# ---------------------------------------------------------------------------
def case_display_names(tmp):
    sys.path.insert(0, HERE)
    import lib_palette as palette

    check("run-name spellings become manuscript spellings",
          [palette.display(x) for x in ("GluC", "LysC", "Trypsin", "All")],
          ["Glu-C", "Lys-C", "Trypsin", "All proteases"])
    check("already-correct spellings survive",
          [palette.display(x) for x in ("Glu-C", "Lys-C", "All proteases")],
          ["Glu-C", "Lys-C", "All proteases"])
    check("display folds the same way colour does",
          {palette.display(x) for x in ("gluc", "Glu-C", "GLUC")}, {"Glu-C"})
    check("an unknown name passes through", palette.display("AspN2"), "AspN2")
    # two spellings of one protease share a colour and a label
    col = palette.assign(["GluC", "Glu-C"])
    check("one protease, one colour", len(set(col.values())), 1)


# ---------------------------------------------------------------------------
# 17. The median convention shared by the rank curves and the count panel.
# ---------------------------------------------------------------------------
def case_median(tmp):
    sys.path.insert(0, HERE)
    import fig1cd_coverage as pc

    check("odd count takes the middle", pc.med([3, 1, 2]), 2)
    check("even count averages the middle pair", pc.med([1, 2, 3, 4]), 2.5)
    check("empty is zero", pc.med([]), 0)

    # a curve must end at the median proteoform count, not one rank past it
    lists = [[1.0] * n for n in (4, 6, 8, 10)]          # median count 7
    check("median curve ends at the median count",
          len(pc.median_curve(lists)), 7)
    lists = [[1.0] * n for n in (4, 6, 9)]              # odd: median count 6
    check("odd sample count too", len(pc.median_curve(lists)), 6)

    # a sample that ran out contributes 0 to the tail
    lists = [[0.9, 0.8], [0.9, 0.8, 0.7, 0.6, 0.5]]
    c = pc.median_curve(lists)
    check("short samples pull the tail down", len(c), 4)   # median of 2 and 5
    check("the median is over every sample, absent counted as zero",
          round(c[2], 3), 0.35)
    # ties round half up rather than to even
    check("a .5 median rounds up, not to even",
          (len(pc.median_curve([[1.0] * 2, [1.0] * 3])),
           len(pc.median_curve([[1.0] * 3, [1.0] * 4]))), (3, 4))


def case_pack_hits(tmp):
    """Peptide rows on the isoform strip: packed by overlap, ends exclusive."""
    sys.path.insert(0, HERE)
    import fig2a_isoform_strip as ie

    # abutting products of one digest do not overlap and share a row
    hits = [(0, "PEPTIDEK", "Trypsin"), (8, "SAMPLER", "Trypsin")]
    rows, n = ie.pack_hits(hits)
    check("abutting peptides share a row", n, 1)
    check("and are both placed", len(rows), 2)
    check("ends are exclusive", sorted((a, b) for a, b, _d, _t in rows),
          [(0, 8), (8, 15)])

    # one residue of genuine overlap does take a second row
    rows, n = ie.pack_hits([(0, "PEPTIDEK", "Trypsin"), (7, "KSAMPLER", "Trypsin")])
    check("overlapping peptides do not", n, 2)

    # the protease is irrelevant to packing
    hits = [(0, "AAAA", "Trypsin"), (40, "CCCC", "GluC"), (80, "DDDD", "LysC")]
    rows, n = ie.pack_hits(hits)
    check("non-overlapping peptides share a row across digests", n, 1)
    check("and all three are placed", len(rows), 3)

    # overlap still forces a new row, and still does so across digests
    rows, n = ie.pack_hits([(0, "AAAAAAAA", "Trypsin"), (4, "CCCCCCCC", "GluC")])
    check("overlap across digests takes a second row", n, 2)

    # the mapper places a peptide in drawing coordinates before packing
    rows, _n = ie.pack_hits([(0, "PEPK", "Trypsin")], lambda j: j + 100)
    check("the mapper is applied to the start", rows[0][0], 100)
    check("and the end follows it", rows[0][1], 104)

    # a junction-spanning peptide's drawn end is its last residue's column
    gap = lambda j: j if j < 10 else j + 40          # noqa: E731
    rows, _n = ie.pack_hits([(6, "PEPTIDEK", "Trypsin")], gap)
    check("a junction-spanning peptide reserves the whole span",
          (rows[0][0], rows[0][1]), (6, 54))
    rows, n = ie.pack_hits([(6, "PEPTIDEK", "Trypsin"),
                            (12, "SAMPLER", "Trypsin")], gap)
    check("so a peptide inside that span is not packed under it", n, 2)

    check("no peptides, no rows", ie.pack_hits([])[1], 0)


def case_aligned_strip(tmp):
    """Both forms on one frame: each change typed, bars from the margin, title fixed."""
    import fig2a_isoform_strip as ie
    from lib_fasta import read_fasta

    can = "M" + "A" * 40 + "C" * 40
    ins = can[:41] + "WWWWWWWWWW" + can[41:]         # 10 residues inserted at 41
    dele = can[:41] + can[51:]                       # 10 residues deleted at 41
    repl = can[:41] + "WWWWW" + can[51:]             # 10 replaced by 5
    for name, iso, tag in (("an insertion", ins, "insert"),
                           ("a deletion", dele, "delete"),
                           ("a replacement", repl, "replace")):
        check(f"{name} is seen as one", [e[0] for e in ie.events(can, iso)], [tag])

    fa = os.path.join(tmp, "aligned.fasta")
    with open(fa, "w") as fh:
        for acc, seq in (("P00001", can), ("P00001-2", ins),
                         ("P00001-3", repl)):
            fh.write(f">sp|{acc}|T_HUMAN test\n{seq}\n")
    check("the fasta round-trips", read_fasta(fa)["P00001-2"], ins)

    tail = can[:41] + "WWWWW" + can[51:] + "K" * 20      # replaced, then longer
    out = os.path.join(tmp, "aligned.svg")
    ie.aligned_panel("P00001", "P00001-4", can, tail, {"Trypsin": {can[30:45]}},
                     {"Trypsin": {tail[40:60]}}, out, label="GENE", margin=80,
                     bar_height=22)
    svg = open(out).read()
    check("--margin starts the bars at the margin",
          re.search(r'<rect x="80\.0" y="[\d.]+" width="[\d.]+" height="22\.0"',
                    svg) is not None, True)
    check("... and the title does not move with them",
          re.search(r'<text x="58\.0" y="32\.0"[^>]*>GENE<', svg) is not None,
          True)
    check("each form is labelled with its length",
          (f">{len(can)} aa<" in svg, f">{len(tail)} aa<" in svg), (True, True))


def case_sample_ids(tmp):
    """Run name -> patient id, for every naming convention in the cohort, by one parser."""
    sys.path.insert(0, HERE)
    from lib_report import sample_of

    check("the 60-min naming resolves",
          sample_of("2026-03-04_CF_2522_L_new_A6_1_8022"), "2522")
    # a token between CF and the id, no separator before the protease letter
    check("a run with a token between CF and the id resolves",
          sample_of("2026-08-27_CF_PD_2522T_G-A1_1_9301"), "2522")
    check("and the 30-min naming too",
          sample_of("2026-08-27_CF_PD_1488G_30min_G-A1_1_9329"), "1488")
    # four digits, not \d+: an unanchored \d+ takes the leading year
    check("the date is not mistaken for a patient id",
          sample_of("2026-03-04_CF_1825_T_new_A6_1_8004"), "1825")
    check("an unparseable name returns None rather than guessing",
          sample_of("blank_A1_1_0001"), None)

    # no other module keeps a private copy of the rule
    copies = [n for n in sorted(os.listdir(HERE))
              if n.endswith(".py") and n not in ("lib_report.py", "test_units.py")
              and re.search(r'''(?:compile|search|match)\(\s*r?["']CF''',
                            open(os.path.join(HERE, n)).read())]
    check("no module keeps a private copy of the rule", copies, [])


def case_compose(tmp):
    """Stacking two panels: common width, summed height, nothing dropped."""
    sys.path.insert(0, HERE)
    import lib_compose as compose

    stack = 'system-ui, "Segoe UI", sans-serif'    # note the inner quotes
    a = os.path.join(tmp, "a.svg")
    b = os.path.join(tmp, "b.svg")
    with open(a, "w") as fh:                       # 100 wide in its own units
        fh.write('<svg xmlns="http://www.w3.org/2000/svg" width="200" '
                 f"height=\"100\" viewBox=\"0 0 100 50\" font-family='{stack}'>"
                 '<rect id="A"/><text>A</text></svg>')
    with open(b, "w") as fh:                       # 50 wide, half the scale
        fh.write('<svg xmlns="http://www.w3.org/2000/svg" width="100" '
                 f"height=\"60\" viewBox=\"0 0 50 30\" font-family='{stack}'>"
                 '<rect id="B"/><text>B</text></svg>')

    out = os.path.join(tmp, "stacked.svg")
    w, h = compose.compose([a, b], out, gap=10)
    src = open(out).read()

    # widest input sets the width; each panel scales to it
    check("width is the widest input", w, 200.0)
    # a: 50 * (200/100) = 100;  b: 30 * (200/50) = 120;  + one 10 gap
    check("height is the scaled sum plus gaps", h, 230.0)
    check("both bodies survive", ("A" in src) and ("B" in src), True)
    check("no nested <svg> wrappers", src.count("<svg"), 1)

    # the second panel is offset by the first's scaled height plus the gap
    check("second panel is offset by the scaled height",
          "translate(0 110.00)" in src, True)
    check("first panel sits at the origin", "translate(0 0.00)" in src, True)
    check("each panel keeps its own units", "scale(2.000000)" in src
          and "scale(4.000000)" in src, True)

    # the wrapper's font stack moves onto each panel's <g>
    check("each panel carries its font stack", src.count(f"font-family='{stack}'"), 2)
    check("and it is quoted so the inner quotes survive",
          ET.fromstring(src) is not None, True)
    g = [e for e in ET.fromstring(src)
         if e.tag.endswith("}g") or e.tag == "g"]
    check("the font is on the group, not the root",
          all(e.get("font-family") == stack for e in g), True)

    # an explicit width overrides the default and rescales both
    out2 = os.path.join(tmp, "stacked2.svg")
    w2, h2 = compose.compose([a, b], out2, gap=0, width=100.0)
    check("explicit width is honoured", w2, 100.0)
    check("and the height follows it", h2, 110.0)   # 50*1 + 30*2

    # a shared colour key takes its hexes from palette.assign
    from lib_palette import assign
    out3 = os.path.join(tmp, "keyed.svg")
    names = ["Glu-C", "Lys-C", "Trypsin"]
    compose.compose([a, b], out3, gap=0, letter="a", key=names)
    src3 = open(out3).read()
    want = assign(list(names))
    check("the key takes its colours from the palette",
          all(f'fill="{want[n]}"' in src3 for n in names), True)
    check("and names every category it colours",
          all(f">{n}</text>" in src3 for n in names), True)
    check("the composite is still well-formed XML",
          ET.fromstring(src3) is not None, True)
    # right-aligned to end inside the canvas
    xs = [float(m) for m in re.findall(r'<rect x="([\d.]+)" y="[\d.]+" '
                                       r'width="[\d.]+" height="[\d.]+" rx="2"', src3)]
    check("all three swatches are on the canvas", len(xs), 3)
    check("and none runs past the right edge", max(xs) + 9.0 <= 200.0, True)
    check("nor off the left one", min(xs) >= 0.0, True)
    check("the key does not change the height",
          compose.parse(out3)[2], 220.0)

    # a row composite has no font on its root, so a letter takes its panels' stack
    row = os.path.join(tmp, "row.svg")
    compose.compose_row([a, b], row)
    check("a row composite has no font on its root",
          ET.parse(row).getroot().get("font-family"), None)
    check("parse recovers the stack its panels share", compose.parse(row)[4], stack)
    out4 = os.path.join(tmp, "row_lettered.svg")
    compose.compose([row], out4, letter="c")
    lab = [e for e in ET.parse(out4).getroot().iter()
           if e.tag.endswith("text") and e.text == "c"]
    check("the letter stamped on a nested composite carries the font",
          [e.get("font-family") for e in lab], [stack])

    # a whole figure from lettered panels drawn at 2 px per unit, letters written two ways
    la = os.path.join(tmp, "la.svg")
    lb = os.path.join(tmp, "lb.svg")
    with open(la, "w") as fh:
        fh.write('<svg xmlns="http://www.w3.org/2000/svg" width="200" '
                 f"height=\"100\" viewBox=\"0 0 100 50\" font-family='{stack}'>"
                 '<text x="22.0" y="32.0" font-size="16.25" fill="#0b0b0b" '
                 'text-anchor="start" font-weight="600">a</text>'
                 '<text x="40" y="20" font-weight="600">Title</text></svg>')
    with open(lb, "w") as fh:
        fh.write('<svg xmlns="http://www.w3.org/2000/svg" width="100" '
                 f"height=\"60\" viewBox=\"0 0 50 30\" font-family='{stack}'>"
                 '<text x="40.8" y="36.1" font-size="22.1" font-weight="600" '
                 'fill="#0b0b0b">b</text><text x="10" y="40">x</text></svg>')
    whole = os.path.join(tmp, "whole.svg")
    w5, h5 = compose.compose([la, lb], whole, natural=True, reletter=True)
    src5 = open(whole).read()
    check("natural: the widest panel sets the width", w5, 200.0)
    check("natural: each panel keeps its own scale", h5, 160.0)   # 50*2 + 30*2
    check("natural: neither panel is stretched",
          src5.count("scale(2.000000)"), 2)
    root5 = ET.fromstring(src5)
    lets = [(e.text, e.get("x"), e.get("y"), e.get("font-size"))
            for e in root5.iter() if e.tag.endswith("text")
            and e.text in ("a", "b")]
    check("reletter: one letter per panel, at one size and corner",
          lets, [("a", "22", "32", f"{__import__('lib_palette').pt(compose.LETTER_PT):.2f}"),
                 ("b", "22", "132", f"{__import__('lib_palette').pt(compose.LETTER_PT):.2f}")])
    check("reletter: the letters are outside the panels' scaled groups",
          all(ch.text not in ("a", "b") for g in root5 if g.tag.endswith("g")
              for ch in g.iter()), True)
    check("reletter: a bold title and an unbolded letter-sized label survive",
          ">Title</text>" in src5 and ">x</text>" in src5, True)
    rown = os.path.join(tmp, "row_natural.svg")
    w6, h6 = compose.compose_row([la, lb], rown, gap=10, natural=True)
    check("row natural: widths add, the row is as tall as the tallest panel",
          (w6, h6), (310.0, 100.0))              # 200 + 10 + 100; max(100, 60)
    check("row natural: neither panel is rescaled",
          open(rown).read().count("scale(2.000000)"), 2)
    wide = os.path.join(tmp, "stretched.svg")
    compose.compose([la, lb], wide)
    check("without natural the narrow panel is stretched",
          "scale(4.000000)" in open(wide).read(), True)

    lc = os.path.join(tmp, "lc.svg")                # two candidates: refuse
    with open(lc, "w") as fh:
        fh.write('<svg xmlns="http://www.w3.org/2000/svg" width="100" '
                 'height="60" viewBox="0 0 50 30">'
                 '<text x="20" y="30" font-weight="600">c</text>'
                 '<text x="30" y="40" font-weight="600">d</text></svg>')
    try:
        compose.compose([la, lc], os.path.join(tmp, "bad.svg"), natural=True,
                        reletter=True)
        refused = False
    except SystemExit:
        refused = True
    check("reletter refuses a panel with two candidate letters", refused, True)


# ---------------------------------------------------------------------------
# 23. Isoform DA (Fig. 2b, Extended Data Fig. 1): difference of differences and its plots.
# ---------------------------------------------------------------------------
def case_isoform_da(tmp):
    sys.path.insert(0, HERE)
    import fig2b_ed1_isoform_da as dp

    # a gene that falls as a whole reads Δ 0; only an isoform-specific change reads otherwise
    cond = {"a": "LBD", "b": "LBD", "c": "LBD",
            "d": "Control", "e": "Control", "f": "Control"}
    med = {(x, "Trypsin"): 0.0 for x in cond}
    down = dict(a=1.0, b=1.1, c=0.9, d=4.0, e=4.2, f=3.8)
    flat = dict(a=4.0, b=4.1, c=3.9, d=4.0, e=4.2, f=3.8)

    def fit(iso):
        idx = {("P1", "Trypsin"): {"CANK": down},
               ("P1-2", "Trypsin"): {"ISOK": iso}}
        return dp.dd_fit([("P1", "CANK")], [("P1-2", "ISOK")], idx, med, cond)
    f = fit(down)
    check("a whole-gene fall reads Δ 0", round(f["interaction"], 6), 0.0)
    check("... with one value per patient per peptide and form",
          sorted((k, len(v)) for k, v in f["rows"].items()),
          [(("CANK", "Trypsin", 0.0), 6), (("ISOK", "Trypsin", 1.0), 6)])
    f = fit(flat)
    check("an isoform spared by the fall reads Δ > 0", f["interaction"] > 1.5, True)
    check("... and is significant", f["p"] < 0.01, True)
    check("fewer than two patients in a group is no fit",
          dp.dd_fit([("P1", "CANK")], [("P1-2", "ISOK")],
                    {("P1", "Trypsin"): {"CANK": down},
                     ("P1-2", "Trypsin"): {"ISOK": dict(a=1.0, d=4.0)}},
                    med, cond), None)

    pv = dp.precursor_values([("P1", "CANK")], {("P1", "Trypsin"): {"CANK": down}},
                             med)
    check("a precursor value per patient", len(pv), 6)
    check("... centred on the peptide's mean",
          round(sum(v for *_r, v in pv), 9), 0.0)

    ns = "{http://www.w3.org/2000/svg}"
    mp = os.path.join(tmp, "model.svg")
    dp.model_plots([(f, "G1", "P1-2", "P1")], mp)
    check("the model plot draws each case value once",
          len(list(ET.parse(mp).getroot().iter(f"{ns}circle"))), 6)

    # volcano: every row a point, only the BH survivors labelled
    vr = [{"gene": f"G{i}", "log2fc": (-1) ** i * 0.1 * i, "p": 10 ** -i,
           "q": 10 ** -i * 3} for i in range(1, 7)]
    out = os.path.join(tmp, "volcano.svg")
    dp.volcano(vr, out)
    root = ET.parse(out).getroot()
    texts = {t.text for t in root.iter(f"{ns}text")}
    check("the volcano draws a point per isoform",
          len(list(root.iter(f"{ns}circle"))), 6)
    check("... and labels only the q < 0.05 survivors",
          {f"G{i}" for i in range(1, 7)} & texts, {"G2", "G3", "G4", "G5", "G6"})
    check("a Tukey whisker stops at the last point inside 1.5 IQR",
          dp._quartiles([0.0, 1.0, 2.0, 3.0, 100.0])[4], 3.0)


def case_print_type(tmp):
    """Type size at print: every enclosing scale counts, and the figure prints PRINT_W_PT wide."""
    sys.path.insert(0, HERE)
    import audit
    from lib_palette import PRINT_W_PT, pt
    from lib_svg import text_width
    svg = os.path.join(tmp, "type.svg")
    w = 2 * PRINT_W_PT                          # 2 units per printed pt
    with open(svg, "w") as fh:
        fh.write(f'<svg xmlns="http://www.w3.org/2000/svg" width="{w:g}" height="100" '
                 f'viewBox="0 0 {w:g} 100"><text font-size="20">root</text>'
                 '<g transform="translate(5 5) scale(0.5)"><text font-size="20">half</text>'
                 '<g transform="scale(3)"><text font-size="4">nested</text></g></g>'
                 '<text font-size="30"> </text></svg>')
    h, sizes = audit.printed_type(svg)
    check("printed size follows every enclosing scale",
          sorted((t, round(s_, 2)) for s_, t in sizes),
          [("half", 5.0), ("nested", 3.0), ("root", 10.0)])
    check("printed height is in points", h, 50.0)
    check("pt() is the whole-figure size that prints at that many points",
          round(pt(10) * PRINT_W_PT / 1215.0, 6), 10.0)
    check("text_width is near a measured render (Trypsin: 2.95 em)",
          2.8 < text_width("Trypsin", 1.0) < 3.2, True)


def case_diagnostic_peptides(tmp):
    sys.path.insert(0, HERE)
    import fig2b_ed1_isoform_da as dp

    # "diagnostic" means absent from the canonical, not merely in a group naming an isoform
    canon = "MAAAKPEDR" + "GGSTVLQEK" + "YYFNMPWLR"
    ins = "HHHDDDWWWR"
    seqs = {"P1": canon,
            "P1-2": canon[:18] + ins + canon[18:],   # an insertion at 18
            "P1-3": canon[:9] + canon[18:]}          # a deletion of 10-18
    genes = {"P1": "G1", "P1-2": "G1", "P1-3": "G1"}
    shared = "GGSTVLQEK"                             # in every form
    idx = {("P1;P1-2", "Trypsin"): {shared: {"s1": 10.0}},
           ("P1-2", "Trypsin"): {ins: {"s1": 10.0}},
           ("P1-2", "GluC"): {ins: {"s1": 10.0}},
           ("P1", "Trypsin"): {"YYFNMPWLR": {"s1": 10.0}}}
    base, diag = dp.diagnostic("G1", {"P1;P1-2", "P1-2", "P1"}, seqs, genes, idx)
    check("the canonical column takes the peptides it contains",
          {p for _g, p in base[1]}, {shared, "YYFNMPWLR"})
    check("... including one reported under a group naming the isoform",
          ("P1;P1-2", shared) in base[1], True)
    check("only one isoform column, and only the inserted peptide in it",
          {who: {p for _g, p in m} for who, m in diag.items()},
          {("P1-2",): {ins}})

    # a gene with no single canonical entry is refused, not given a guessed reference
    b2, d2 = dp.diagnostic("G2", {"Q1-2"}, {"Q1-2": canon},
                           {"Q1-2": "G2"}, {})
    check("an isoform-only gene is refused", (b2, d2), (None, {}))

    # column: case minus control of each patient's summed peptides, per digest
    cond = {"a": "LBD", "b": "LBD", "c": "LBD",
            "d": "Control", "e": "Control", "f": "Control"}
    med2 = {(x, "Trypsin"): 0.0 for x in cond}
    idx2 = {("P1-2", "Trypsin"): {
        "AAAK": dict(a=4.0, b=4.1, c=3.9, d=1.0, e=1.1, f=0.9),
        "CCCK": dict(a=2.0, b=2.0, c=2.0, d=2.0, e=2.0, f=2.0)}}
    mem2 = [("P1-2", "AAAK"), ("P1-2", "CCCK")]
    st = dp.column(mem2, idx2, med2, cond, "LBD", "Control")
    check("column: one change per digest that measured both groups",
          sorted(st["per_digest"]), ["Trypsin"])
    check("... the log2 ratio of the summed peptides",
          round(st["mean"], 1), 1.0)
    check("... counting the peptides seen in both groups", st["npep"], 2)
    check("... with a patient-level p", 0 < st["p"] < 0.01, True)
    check("no digest measuring both groups is no column",
          dp.column(mem2, {("P1-2", "Trypsin"): {
              "AAAK": dict(a=4.0, b=4.0)}}, med2, cond, "LBD", "Control"), None)


def case_report(tmp):
    sys.path.insert(0, HERE)
    import lib_report as rp

    check("normalising collapses every separator",
          {rp._norm(x) for x in ("EG.Qvalue", "EG_Qvalue", "eg qvalue",
                                 "EG-QVALUE")}, {"egqvalue"})

    # -- the schema says which software, and what each column means -----------
    SN = ["R.FileName", "PG.ProteinGroups", "PEP.StrippedSequence",
          "EG.Qvalue", "FG.Quantity", "EG.UsedForPeptideQuantity",
          "EG.ApexRT", "PG.QValue (Run-Wise)"]
    DIANN = ["File.Name", "Protein.Group", "Stripped.Sequence", "Q.Value",
             "Precursor.Quantity", "RT", "Precursor.Id", "Precursor.Charge"]
    MQ = ["Raw file", "Modified sequence", "Leading razor protein",
          "Missed cleavages", "Intensity", "Reverse", "Sequence", "Charge"]
    check("Spectronaut is recognised", rp.identify(SN)[0], "Spectronaut")
    check("... and so is the underscored re-export",
          rp.identify([c.replace(".", "_").replace(" ", "_") for c in SN])[0],
          "Spectronaut")
    check("DIA-NN is recognised", rp.identify(DIANN)[0], "DIA-NN")
    check("MaxQuant is recognised", rp.identify(MQ)[0], "MaxQuant")
    check("and nothing is claimed for an unrelated table",
          rp.identify(["gene", "log2fc", "padj"])[0], None)

    # the same canonical fields come out of all three
    for name, cols in (("Spectronaut", SN), ("DIA-NN", DIANN),
                       ("MaxQuant", MQ)):
        got = rp.resolve(cols)
        check(f"{name} resolves run/peptide/quantity",
              all(f in got for f in ("run", "peptide", "quantity")), True)
    check("dots and underscores resolve identically",
          rp.resolve(SN),
          {k: v.replace("_", ".") if v.count("_") and " " not in v else v
           for k, v in rp.resolve([c.replace(".", "_") for c in SN]).items()}
          | {"protein_q_run": "PG.QValue (Run-Wise)"})
    check("an absent field is absent rather than guessed",
          "decoy" in rp.resolve(SN), False)

    # -- the content says which protease -------------------------------------
    def peps(ends, n=400):
        return {f"AAAK{i:04d}"[:8].replace("K", "A") + ends[i % len(ends)]
                for i in range(n)}
    check("all-E termini is Glu-C", rp.protease(peps("E"))[0], "GluC")
    check("all-K termini is Lys-C", rp.protease(peps("K"))[0], "LysC")
    # half K, half R is trypsin, not Lys-C
    got, ev = rp.protease(peps("KR"))
    check("half K half R is trypsin, not Lys-C", got, "Trypsin")
    check("... and it reports both residues",
          sorted(ev["per_residue"]), ["K", "R"])
    check("Lys-C is not offered as an alternative for it",
          "LysC" in ev.get("also_fits", []), False)
    check("N-terminal enzymes use the other end",
          rp.protease({f"D{i:07d}" for i in range(300)})[0], "AspN")
    bad, ev = rp.protease({f"AAAAAA{c}" for c in "ACDFGHILMNPQSTVWY" * 20})
    check("an unexplained terminus distribution is refused", bad, None)
    check("... and says what it was closest to", "closest" in ev, True)
    check("no peptides is not a protease call", rp.protease([])[0], None)

    # -- reading: a numeric column stored as text is coerced ------------------
    import pyarrow as pa
    import pyarrow.parquet as pq
    p = os.path.join(tmp, "sn_like.parquet")
    pq.write_table(pa.table({
        "R_FileName": ["2026-01-01_CF_PD_1488G_30min_A1_1_9329"] * 4,
        "PG_ProteinGroups": ["P1", "P1", "P2", "P2"],
        "PEP_StrippedSequence": ["AAAAAAE", "CCCCCCE", "DDDDDDE", "EEEEEEE"],
        # q-values as text with a capital E, as some re-exports carry them
        "EG_Qvalue": ["1.9967115436590949E-13", "0.004", "NaN", "0.5"],
        "FG_Quantity": [10.0, 20.0, 30.0, 40.0],
        "EG_UsedForPeptideQuantity": [True, True, True, True],
        "EG_ApexRT": [1.0, 2.0, 3.0, 4.0],
    }), p)
    r = rp.Report(p)
    check("a report is recognised through underscored columns",
          r.software, "Spectronaut")
    check("and its protease comes from the peptides", r.protease, "GluC")
    check("the text q-value column is noticed",
          str(r.schema.field(r.columns["precursor_q"]).type), "string")
    b = next(r.batches(["precursor_q", "quantity", "peptide"]))
    check("... and coerced to numbers on read",
          [type(x).__name__ for x in b["precursor_q"][:2]], ["float", "float"])
    check("the text NaN becomes a null, not a string",
          b["precursor_q"][2], None)
    check("so a threshold compares numerically",
          [x is not None and x <= 0.01 for x in b["precursor_q"]],
          [True, True, False, False])
    r.require("run", "peptide", "quantity")          # must not exit
    check("asking for a field the file lacks is caught",
          r.has("decoy"), False)

    # -- a declaration that contradicts the data stops the run ---------------
    try:
        rp.open_reports([p], declared={os.path.basename(p): "Trypsin"})
        check("a contradicting declaration is refused", "no exit", "SystemExit")
    except SystemExit as e:
        check("a contradicting declaration is refused",
              "declared Trypsin" in str(e), True)
    got = rp.open_reports([p], declared={os.path.basename(p): "GluC"})
    check("an agreeing declaration is accepted", got[0].protease, "GluC")

    # -- two searches of the same injections are refused ---------------------
    same = os.path.join(tmp, "sn_like_research.parquet")
    other = os.path.join(tmp, "sn_like_reinjected.parquet")
    t = pq.read_table(p)
    pq.write_table(t, same)
    pq.write_table(t.set_column(
        t.schema.get_field_index("R_FileName"), "R_FileName",
        pa.array(["2026-02-02_CF_PD_1488G_30min_A1_1_9400"] * 4)), other)
    check("a report's run names are read", rp.runs_of(p),
          {"2026-01-01_CF_PD_1488G_30min_A1_1_9329"})
    try:
        rp.open_reports([p, same])
        check("two reports sharing runs are refused", "no exit", "SystemExit")
    except SystemExit as e:
        check("two reports sharing runs are refused",
              "pooling them counts every precursor twice" in str(e), True)
    check("a re-injection with its own run names is accepted",
          len(rp.open_reports([p, other])), 2)
    tsv = os.path.join(tmp, "sn_like.tsv")
    write_report(tsv, ["R.FileName", "PEP.StrippedSequence"],
                 [["2026-01-01_CF_PD_1488G_30min_A1_1_9329", "AAAAAAE"]])
    try:
        rp.one_search_per_run([p, tsv])
        check("... across parquet and text reports alike", "no exit", "SystemExit")
    except SystemExit:
        check("... across parquet and text reports alike", True, True)


def case_isoform_coverage(tmp):
    """Discriminating residues and deletion junctions (Fig. 2c-d)."""
    sys.path.insert(0, HERE)
    import fig2cd_isoform_coverage as ic

    # AAAK|BBBK joins at 4: AKBB spans it, a peptide ending or starting there does not
    check("junction_peptides: only a peptide holding both sides of the join spans it",
          ic.junction_peptides("AAAKBBBK", [4],
                               {"LysC": {"AAAK", "BBBK"}, "Trypsin": {"AKBB"}}),
          {"GluC": set(), "LysC": set(), "Trypsin": {"AKBB"}})
    # S3 is shared; T7 S8 T9 are the isoform's own
    iseqs = {"P7": "MAASAKGSGGAKLLLLK", "P7-2": "MAASAKGTSTAKLLLLK"}
    check("the isoform's own residues are the replaced ones",
          sorted(ic.discriminating(iseqs["P7"], iseqs["P7-2"])), [7, 8, 9])
    check("each side's own residues come from one alignment",
          tuple(sorted(x) for x in ic.own_residues("MAAGGGKLL", "MAAXKLL")),
          ([3, 4, 5], [3]))
    # a pure deletion's only evidence is a peptide spanning the join in one piece
    jc, ji = "MKTAYIAKQRWWWWWWPEPTIDEK", "MKTAYIAKQRPEPTIDEK"
    check("a pure deletion has no own residue but one junction",
          (ic.discriminating(jc, ji), ic.junctions(jc, ji)), (set(), [10]))
    check("a terminal deletion is not a junction",
          ic.junctions("MKTAYIAKQRWWW", "MKTAYIAKQR"), [])
    jcov = ic.junction_coverage(ji, [10], {"Trypsin": {"QRPEPT", "AKQR"},
                                           "LysC": {"AKQR", "PEPTIDEK"}})
    check("a peptide across the join covers it; two meeting at it do not",
          (jcov["Trypsin"], jcov["LysC"]), ({10}, set()))


def case_phospho(tmp):
    sys.path.insert(0, HERE)
    import prep_phospho as ph

    # -- offsets: an N-terminal modifier precedes residue 0 and shifts nothing
    check("a phosphate on the first residue",
          ph.phospho_offsets("_S[Phospho (STY)]PEK_"), [0])
    check("an N-terminal acetyl does not shift the phospho offsets",
          ph.phospho_offsets("_[Acetyl (Protein N-term)]AS[Phospho (STY)]"
                             "PT[Phospho (STY)]K_"), [1, 3])
    check("other residue modifications are skipped",
          ph.phospho_offsets("_AAM[Oxidation (M)]S[Phospho (STY)]K_"), [3])

    # -- location: every canonical entry, isoforms only as a fallback --------
    seqs = {"P1": "MKKAESPVKEEAVAEK", "P2": "GGKAESPVKEE",
            "P1-2": "MKKAESPVKEEAVAEKWQTSPRK", "P3-2": "QQWQTSPRK"}
    loc = ph.locate_all(["KAESPVK", "WQTSPRK", "NOTTHERE"], seqs)
    check("a shared peptide is placed in every canonical entry",
          loc["KAESPVK"], [("P1", 2), ("P2", 2)])
    check("an isoform-only peptide is placed in its isoform entries",
          loc["WQTSPRK"], [("P1-2", 16), ("P3-2", 2)])
    check("an absent peptide is placed nowhere", loc["NOTTHERE"], [])

    # -- a hand-built scan: two patients, three digests ----------------------
    scan = os.path.join(tmp, "phospho_scan")
    os.makedirs(scan, exist_ok=True)
    runs = [["GluC", "CF_1111_G", "1111", 5000, 3, 1],
            ["LysC", "CF_1111_L", "1111", 5000, 3, 1],
            ["Trypsin", "CF_1111_T", "1111", 5000, 3, 1],
            ["GluC", "CF_2222_G", "2222", 5000, 1, 1],
            ["Trypsin", "CF_2222_T", "2222", 5000, 1, 1],
            ["LysC", "CF_2222_L", "2222", 12, 1, 0]]          # below the gate
    write_report(os.path.join(scan, "phospho_runs.tsv"),
                 ["digest", "run", "sample", "precursors", "phospho_precursors",
                  "kept"], runs)
    # P1 = MKKAESPVKEEAVAEK: S6 carries the phosphate
    head = ["digest", "run", "sample", "modified", "peptide", "charge", "q",
            "pep", "quantity", "used_for_pep", "protein_groups",
            "phospho_offsets", "locations"]
    prec = [
        # the same S6 through three different peptides, one per digest
        ["GluC", "CF_1111_G", "1111", "_KKAES[Phospho (STY)]PVKE_", "KKAESPVKE",
         2, 1e-4, 1e-4, 100, 1, "P1", "4", "P1:1"],
        ["LysC", "CF_1111_L", "1111", "_KAES[Phospho (STY)]PVK_", "KAESPVK",
         2, 1e-4, 1e-4, 100, 1, "P1;P2", "3", "P1:2;P2:2"],
        ["Trypsin", "CF_1111_T", "1111", "_AES[Phospho (STY)]PVKEEAVAEK_",
         "AESPVKEEAVAEK", 2, 1e-4, 1e-4, 100, 1, "P1", "2", "P1:3"],
        # the second patient: one row in a gated-out run, one Glu-C isoform placement
        ["LysC", "CF_2222_L", "2222", "_KAES[Phospho (STY)]PVK_", "KAESPVK",
         2, 1e-4, 1e-4, 100, 1, "P1", "3", "P1:2"],
        ["GluC", "CF_2222_G", "2222", "_WQTS[Phospho (STY)]PRK_", "WQTSPRK",
         2, 1e-4, 1e-4, 100, 1, "P1-2", "3", "P1-2:16"],
    ]
    write_report(os.path.join(scan, "phospho_precursors.tsv"), head, prec)
    rows = ph.load(scan)
    check("rows from a gated-out run are not loaded",
          sorted({(r["digest"], r["sample"]) for r in rows}),
          [("GluC", "1111"), ("GluC", "2222"), ("LysC", "1111"),
           ("Trypsin", "1111")])
    S = ph.sites(rows)
    s6 = [s for s in S.values() if "P1:6" in s["keys"]]
    check("three digests' placements on one residue merge into one site",
          len(s6), 1)
    check("... which carries the second entry the shared peptide sits in",
          "P2:6" in s6[0]["keys"], True)
    check("... and all three digests",
          sorted({o[0] for o in s6[0]["obs"]}), ["GluC", "LysC", "Trypsin"])
    runs_read = ph.load_runs(scan)
    check("only patients with a kept run in every digest are complete",
          ph.complete_patients(runs_read), ["1111"])
    no_phospho = [dict(r) for r in runs_read] + [
        {"digest": "LysC", "run": "CF_2222_L2", "sample": "2222",
         "precursors": "5000", "phospho_precursors": "0", "kept": "1"}]
    check("... measured is a property of the runs, not of finding phosphosites",
          ph.complete_patients(no_phospho), ["1111", "2222"])

    # -- export: fig3a_phospho_sites.py reads the same sites back ------------
    import fig3a_phospho_sites as mod_sites
    exp = os.path.join(tmp, "sites_phospho.tsv")
    ph.export_sites(S, exp)
    back = {(r["acc"], r["res"]): r for r in mod_sites.load(exp)}
    check("mod_sites reads every exported site", len(back), len(S))
    check("... with each digest's peptides for the three-digest site",
          back[("P1", 6)]["by_digest"],
          {"GluC": ["KKAESPVKE"], "LysC": ["KAESPVK"], "Trypsin": ["AESPVKEEAVAEK"]})

    # -- peptide_classes: Fig. 3a's per-digest bars -----------------------------
    # a digest counts only its own peptides, All counts every digest's, 5+ caps
    stack = [{"by_digest": {"Trypsin": ["A", "B"], "LysC": ["C"]},
              "peps": ["A", "B", "C"]},
             {"by_digest": {"Trypsin": ["D"]}, "peps": ["D"]},
             {"by_digest": {"GluC": list("EFGHIJ")}, "peps": list("EFGHIJ")}]
    check("peptide_classes: a digest counts only its own peptides",
          [dict(mod_sites.peptide_classes(stack, d)) for d in ("Trypsin", "LysC", "GluC")],
          [{2: 1, 1: 1}, {1: 1}, {5: 1}])
    check("peptide_classes: All counts every digest's, capped at 5+",
          dict(mod_sites.peptide_classes(stack)), {3: 1, 1: 1, 5: 1})
    # -- log_height: the bars' log y axis --------------------------------------
    check("log_height: 1 sits on the baseline, 10**decades at the top, 100 halfway",
          [mod_sites.log_height(v, 4, 230.0) for v in (1, 10 ** 4, 100)],
          [0.0, 230.0, 115.0])
    check("log_height: an empty class draws nothing rather than failing",
          mod_sites.log_height(0, 4, 230.0), 0.0)
    # -- venn: Fig. 3a's regions and where their counts go ---------------------
    vs = [{"digests": d} for d in (["Trypsin"], ["Trypsin", "LysC"],
                                    ["Trypsin", "LysC", "GluC"], ["LysC"],
                                    ["LysC"])]
    reg = mod_sites.venn_regions(vs, ["Trypsin", "LysC", "GluC"])
    check("venn_regions: every combination is a key, empty ones as 0",
          (len(reg), reg[frozenset(["GluC"])], reg[frozenset(["LysC"])],
           reg[frozenset(["Trypsin", "LysC", "GluC"])]), (7, 0, 2, 1))
    check("venn_regions: the regions partition the sites",
          sum(reg.values()), len(vs))
    reg3 = {frozenset(k): v for k, v in [
        (["Trypsin"], 1646), (["LysC"], 1087), (["GluC"], 263),
        (["LysC", "GluC"], 32), (["Trypsin", "GluC"], 42),
        (["Trypsin", "LysC"], 485), (["Trypsin", "LysC", "GluC"], 94)]}
    cen, rad = mod_sites.venn_layout(reg3, ["Trypsin", "LysC", "GluC"], 90.0)
    check("venn_layout: circle areas are in proportion to the totals",
          [round((r / rad[0]) ** 2, 3) for r in rad],
          [1.0, round(1698 / 2267, 3), round(431 / 2267, 3)])
    pts = mod_sites.venn_label_points(cen, rad)

    def members(x, y):
        return frozenset(i for i, ((cx, cy), r) in enumerate(zip(cen, rad))
                         if (x - cx) ** 2 + (y - cy) ** 2 < r ** 2)
    check("venn labels: each count lands inside exactly its own region",
          all(members(x, y) == k for k, (x, y, _f) in pts.items()) and len(pts) == 7,
          True)


def main():
    with tempfile.TemporaryDirectory() as tmp:
        for fn in (case_semantics, case_header_variants, case_unresolved,
                   case_missing_columns, case_figure, case_figure_groups, case_workflow,
                   case_coverage, case_events, case_diverging, case_compact,
                   case_text_bounds, case_display_names, case_median,
                   case_pack_hits, case_aligned_strip, case_sample_ids,
                   case_compose, case_print_type, case_diagnostic_peptides, case_isoform_da,
                   case_report, case_isoform_coverage, case_phospho):
            print(f"\n{fn.__name__}")
            fn(tmp)
    print()
    if failures:
        print(f"{len(failures)} FAILURE(S): " + ", ".join(failures))
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
