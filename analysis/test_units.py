#!/usr/bin/env python3
"""Assertions for prep_counts.py and fig1b_depth.py on hand-built fixtures.

    python3 test_units.py

Every case uses a report small enough that the right answer is countable by
hand, so a wrong count is a failure and not a judgement call.
"""

import csv
import glob
import math
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


# ---------------------------------------------------------------------------
# 1. Filtering, union dedup, isoform collapsing -- all on one tiny fixture.
#
# Sample S1, two proteases. PEPTIDEA is seen by both digests, so the union
# must count it once. P49189 and P49189-3 are two isoform groups but one
# canonical protein. One decoy row and one above-q row must be dropped.
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

    # Default: each unit is filtered by its own q-value. BADPGPEP clears the
    # precursor cutoff but not the protein one, so it counts as a precursor and
    # a peptide but contributes no protein group; BADQPEP is the mirror case.
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

    # --strict requires both cutoffs for every unit, dropping both odd rows
    out2 = os.path.join(tmp, "sem_strict")
    run([COUNTS, rep, "--out", out2, "--strict"])
    s = load_counts(os.path.join(out2, "counts.csv"))
    check("strict trypsin precursors", s[("S1", "Trypsin", "precursors")], 3)
    check("strict trypsin peptides", s[("S1", "Trypsin", "peptides")], 2)
    check("strict trypsin isoform groups",
          s[("S1", "Trypsin", "protein_isoform_groups")], 2)
    check("strict trypsin canonical", s[("S1", "Trypsin", "proteins_canonical")], 1)
    check("strict union precursors", s[("S1", "All", "precursors")], 4)


# ---------------------------------------------------------------------------
# 2. Header variants: different capitalisation, no stripped-sequence column
#    (must fall back to deriving it), precursor id instead of mod sequence.
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
    check("unresolved suggests regex", "--protease-regex" in p.stdout, True)

    # ... and the regex escape hatch works
    out = os.path.join(tmp, "unres2")
    run([COUNTS, rep, "--out", out,
         "--protease-regex", r"run_(alpha|beta)", "--sample-regex", r"_(\d+)$"])
    c = load_counts(os.path.join(out, "counts.csv"))
    check("regex override resolves", sorted({k[1] for k in c}), ["All", "alpha", "beta"])


# ---------------------------------------------------------------------------
# 4. Missing required columns must fail loudly with the header echoed.
# ---------------------------------------------------------------------------
def case_missing_columns(tmp):
    rep = os.path.join(tmp, "bad.tsv")
    write_report(rep, ["Foo", "Bar"], [["1", "2"]])
    p = run([COUNTS, rep, "--out", os.path.join(tmp, "bad")], expect_rc=1)
    check("missing columns names header", "header seen" in p.stderr, True)


# ---------------------------------------------------------------------------
# 5. Figure geometry: valid XML, every point inside its panel, no tick label
#    above the panel top, and the outlier is drawn rather than clipped away.
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
    p = run([FIGURE, csv_path, "--out", svg, "--letters"])
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
# 6. The workflow schematic: valid XML, nothing outside the canvas, the digest
#    actually obeys the specificities, and colours agree with the count panels.
# ---------------------------------------------------------------------------
def case_workflow(tmp):
    import lib_palette as palette
    from lib_palette import INK_MUTED
    svg = os.path.join(tmp, "fig1a.svg")
    p = run([WORKFLOW, "--out", svg, "--length", "300", "--seed", "11"])

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
        check(f"{name} peptides drawn", len(used) > 3, True)

    texts = [t.text for t in root.iter(f"{ns}text")]
    for want in ("Glu-C", "Lys-C", "Trypsin", "dia-PASEF", "directDIA",
                 "All proteases", "All"):
        check(f"label present: {want}", want in texts, True)
    # The union row is drawn like the three above it -- bars where at least one
    # specificity reaches -- not as a full-width strip shaded by digest count.
    check("the union row is bars, not a strip", "covered by" in texts, False)
    union_bars = [r for r in root.iter(f"{ns}rect")
                  if r.get("fill") == palette.UNION
                  and float(r.get("height")) == 12.0]
    check("union drawn as several bars with gaps", len(union_bars) >= 2, True)
    # one colour, one name, said exactly once. The standalone colour key was
    # removed: the labelled chips in the workflow already bind each colour to its
    # protease, so a second mapping was a restatement.
    check("the union is named once", texts.count("All proteases"), 1)
    check("no stale union names", ("Union" in texts, "All three" in texts),
          (False, False))
    check("nothing restates what the lanes already show", "3 digests" in texts,
          False)
    check("no standalone colour key",
          [t for t in texts if t in ("protease", "used in every panel")], [])
    check("no depth key", [t for t in texts if t == "digests"], [])
    # ... but the chips must still carry the mapping, or panel b is unreadable
    for name in ("Glu-C", "Lys-C", "Trypsin"):
        check(f"{name} is still named in the workflow", texts.count(name) >= 2,
              True)

    # The panel carries labels and data, not prose. A long string is allowed only
    # if it is a data statement -- it has to contain a digit. That admits
    # "94% covered - 68% of it more than once" in the terminal box and still
    # rejects the explanatory sentences this check was written to keep off.
    prose = [t for t in texts if t and len(t) > 24 and not any(ch.isdigit()
                                                              for ch in t)]
    check("no prose on the panel", prose, [])
    check("no label is absurdly long",
          max(len(t or "") for t in texts) <= 60, True)
    check("no full stops in labels",
          [t for t in texts if t and t.rstrip().endswith(".")], [])

    # The per-residue letters were removed: at ~7 px on a 95-residue window they
    # are unreadable at print size, and the panel's claim is where the cuts fall.
    # Only the panel letter should remain as a lone glyph.
    singles = [t for t in texts if t and len(t) == 1 and t.isalpha()]
    check("no per-residue sequence letters", len(singles) <= 1, True)
    check("no in-silico footnote", [t for t in texts if t and "in silico" in t], [])

    # no filled triangle should carry a protease colour any more
    tri = [p_ for p_ in root.findall(f"{ns}path")
           if p_.get("fill") in {colours[k] for k in ("Glu-C", "Lys-C", "Trypsin")}]
    # Cut-site markers were removed: on TXNDC9, 50 of 84 sat on a peptide-bar
    # edge and restated it. Only the workflow arrowheads remain.
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

    # the digest respects specificity: Lys-C peptides must end in K (or at the
    # C-terminus), and trypsin peptides in K or R
    sys.path.insert(0, HERE)
    import fig1a_workflow as wf
    seq = wf.illustrative_sequence(300, 11)
    for residues, block_p, label in ((set("K"), False, "Lys-C"),
                                     (set("KR"), True, "Trypsin"),
                                     (set("E"), False, "Glu-C")):
        peps = wf.digest(seq, residues, block_p, 0, 7, 30)
        check(f"{label} peptides obey specificity",
              all(e == len(seq) or seq[e - 1] in residues for s, e in peps), True)
        check(f"{label} peptides within length window",
              all(7 <= e - s <= 30 for s, e in peps), True)
        if block_p:
            check(f"{label} skips K/R-P bonds",
                  all(e == len(seq) or seq[e] != "P" for s, e in peps), True)
    check("stdout reports the Glu-C assumption", "buffer-dependent" in p.stdout, True)



# ---------------------------------------------------------------------------
# 7. fig1cd_coverage: FASTA parsing, coverage maths, isoform mapping, and both
#    output modes on a fixture where the right answers are hand-countable.
# ---------------------------------------------------------------------------
def case_coverage(tmp):
    import fig1cd_coverage as pc

    # --- pure functions ---------------------------------------------------
    fa = os.path.join(tmp, "db.fasta")
    with open(fa, "w") as fh:
        fh.write(">sp|P00001|AAA_HUMAN alpha\nMKAAAWWWKDDEFGHIKLM\n")   # 19 aa
        fh.write(">P00002\nMKAAAWWWKDDEFGHIKLMNPQRST\n")                # 24 aa
        fh.write(">sp|P00001-2|AAA_HUMAN iso\nMKAAAWWWKHIKLM\n")        # del 10-14
    seqs = pc.read_fasta(fa)
    check("fasta reads sp| accessions", sorted(seqs), ["P00001", "P00001-2", "P00002"])
    check("fasta reads bare header", len(seqs["P00002"]), 25)

    # pileup packing: overlapping peptides must land on different tiers, and a
    # later non-overlapping one must reuse the first
    tiers = pc.pack([(0, 10), (5, 15), (20, 30)], gap=1)
    check("overlaps go to separate tiers", tiers, [[(0, 10), (20, 30)], [(5, 15)]])
    check("nested peptide gets its own tier",
          len(pc.pack([(0, 20), (5, 10), (6, 11)], gap=1)), 3)
    check("every peptide is drawn once",
          sum(len(t) for t in pc.pack([(0, 5), (1, 6), (2, 7), (9, 12)], gap=1)), 4)
    check("merged spans", pc.merged(bytearray([1, 1, 0, 0, 1, 1, 1])),
          [(0, 2), (4, 7)])
    check("repeated peptide located twice",
          pc.locate("ABCABC", {"ABC"}), [(0, 3), (3, 6)])

    cov, found = pc.mark("MKAAAWWWKDDEFGHIKLM", {"AAAWWW", "HIKLM", "ZZZ"})
    check("coverage marks residues", sum(cov), 11)
    check("coverage counts found peptides", found, 2)

    # P00001-2 is P00001 with residues 10-14 removed
    m = pc.deletion_map(seqs["P00001"], seqs["P00001-2"])
    check("deletion located", m, (9, 5))
    p, d = m
    # a peptide spanning the junction maps to two canonical intervals
    a = seqs["P00001-2"].find("WWWKHIK")
    check("junction peptide splits",
          pc.iso_segments(a, a + 7, p, d, len(seqs["P00001-2"])),
          [(5, 9), (14, 17)])
    check("interior peptide stays whole",
          pc.iso_segments(0, 4, p, d, len(seqs["P00001-2"])), [(0, 4)])

    # --- both modes end to end -------------------------------------------
    hdr = ["R.FileName", "PG.ProteinGroups", "PEP.StrippedSequence",
           "EG.ModifiedSequence", "FG.Charge", "EG.Qvalue"]
    rows = []
    for runname, pep in (("r_Trypsin", "AAAWWWK"), ("r_Trypsin", "DDEFGHIK"),
                         ("r_GluC", "MKAAAWWWKDDE"), ("r_LysC", "HIKLM")):
        rows.append([runname, "P00001", pep, f"_{pep}_", "2", "0.001"])
    rows.append(["r_LysC", "P00001-2", "WWWKHIK", "_WWWKHIK_", "2", "0.001"])
    rows.append(["r_GluC", "P00001", "NOTINSEQ", "_NOTINSEQ_", "2", "0.001"])
    # one report per enzyme, as Spectronaut exports them: the protease comes
    # from each report's own filename
    reps = []
    for enz in ("Trypsin", "LysC", "GluC"):
        path = os.path.join(tmp, f"PD_{enz}_Report.tsv")
        write_report(path, hdr, [r for r in rows if r[0] == f"r_{enz}"])
        reps.append(path)

    svg = os.path.join(tmp, "prot.svg")
    p1 = run([COVERAGE] + reps + ["--fasta", fa, "--protein", "P00001",
              "--isoform", "P00001-2", "--out", svg])
    check("reports the unmapped peptide", "not found in P00001" in p1.stdout, True)
    check("names the isoform-specific peptide",
          "isoform-specific: WWWKHIK" in p1.stdout, True)
    root = ET.parse(svg).getroot()
    ns = "{http://www.w3.org/2000/svg}"
    H = float(root.get("height"))
    ys = [float(e.get(a)) for e in root.iter() for a in ("y", "y1", "y2")
          if e.get(a) is not None]
    check("protein panel fits its canvas", max(ys) <= H, True)
    # the pileup draws one bar per peptide occurrence, not one merged block
    import lib_palette as palette
    pal = palette.assign(["GluC", "LysC", "Trypsin", "All"])
    bars = [e for e in root.findall(f"{ns}rect")
            if e.get("fill") in {pal[k] for k in ("GluC", "LysC", "Trypsin")}]
    check("a bar per peptide plus a swatch per digest", len(bars) >= 4 + 3, True)

    svg2 = os.path.join(tmp, "rank.svg")
    p2 = run([COVERAGE] + reps + ["--fasta", fa, "--rank", "--out", svg2])
    root2 = ET.parse(svg2).getroot()
    paths = [e for e in root2.findall(f"{ns}path") if e.get("fill") == "none"]
    check("one curve per protease plus union", len(paths), 4)
    check("series named per enzyme",
          all(e in p2.stdout for e in ("Trypsin", "LysC", "GluC", "All")), True)
    check("rank mode reports comparable counts",
          ">=50%" in p2.stdout and "proteoforms" in p2.stdout, True)
    check("rank mode warns medians are not comparable",
          "medians are not" in p2.stdout, True)



# ---------------------------------------------------------------------------
# 8. Shared-y groups: panels inside a group get one axis and one set of tick
#    labels; a single --letter labels the whole figure.
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
# 9. Depth ramps must stay distinguishable from the colour that means "zero".
#    This is the check that caught the dark violet step sitting at 1.72:1
#    against black, where "covered by four peptides" and "not covered at all"
#    would have looked the same.
# ---------------------------------------------------------------------------
def case_ramp_contrast(tmp):
    import lib_palette as palette

    def lum(h):
        def ch(c):
            c /= 255.0
            return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
        r, g, b = (int(h[i:i + 2], 16) for i in (1, 3, 5))
        return 0.2126 * ch(r) + 0.7152 * ch(g) + 0.0722 * ch(b)

    def ratio(a, b):
        la, lb = lum(a) + 0.05, lum(b) + 0.05
        return max(la, lb) / min(la, lb)

    worst_dark = min(ratio(step, palette.EMPTY_DARK)
                     for ramp in palette.RAMPS_ON_DARK.values() for step in ramp)
    check("dark-surface ramps clear black by 3:1", round(worst_dark, 2) >= 3.0, True)

    # The light variant cannot clear its own empty colour on luminance alone --
    # measured 1.77-1.98:1 for the palest step against #f0efec -- so it leans on
    # a hairline outline round uncovered tiles. That is exactly why black is the
    # better default: it separates "covered once" from "not covered" by 3.2:1 or
    # more with no outline at all.
    light_pale = min(ratio(ramp[0], palette.EMPTY)
                     for ramp in palette.RAMPS.values())
    dark_deep = min(ratio(ramp[0], palette.EMPTY_DARK)
                    for ramp in palette.RAMPS_ON_DARK.values())
    check("light variant's palest step is weak against its empty",
          light_pale < 2.0, True)
    check("black separates depth 1 far better than the light empty does",
          dark_deep > light_pale * 1.5, True)

    # each ramp must be monotone in lightness, or it is not a scale
    for name, ramp in palette.RAMPS_ON_DARK.items():
        ls = [lum(x) for x in ramp]
        check(f"dark ramp {name} runs dark->light", ls == sorted(ls), True)
    for name, ramp in palette.RAMPS.items():
        ls = [lum(x) for x in ramp]
        check(f"light ramp {name} runs light->dark", ls == sorted(ls, reverse=True),
              True)

    # and the stack layout must render with square, touching tiles
    csv_path = os.path.join(tmp, "unused.csv")
    fa = os.path.join(tmp, "stack.fasta")
    with open(fa, "w") as fh:
        fh.write(">sp|P00001|A_HUMAN a\nMKAAAWWWKDDEFGHIKLMNPQRST\n")
    hdr = ["R.FileName", "PG.ProteinGroups", "PEP.StrippedSequence",
           "EG.ModifiedSequence", "FG.Charge", "EG.Qvalue"]
    reps = []
    for enz, pep in (("Trypsin", "AAAWWWK"), ("LysC", "DDEFGHIK"),
                     ("GluC", "MKAAAWWWKDDE")):
        path = os.path.join(tmp, f"PD_{enz}_Report.tsv")
        write_report(path, hdr, [[f"r_{enz}", "P00001", pep, f"_{pep}_", "2",
                                  "0.001"]])
        reps.append(path)
    svg = os.path.join(tmp, "stack.svg")
    run([COVERAGE] + reps + ["--fasta", fa, "--protein", "P00001",
                             "--layout", "stack", "--out", svg])
    root = ET.parse(svg).getroot()
    ns = "{http://www.w3.org/2000/svg}"
    tiles = [r for r in root.findall(f"{ns}rect") if r.get("rx") is None]
    check("stack tiles are square (no rx)", len(tiles) > 20, True)
    fills = {r.get("fill") for r in tiles}
    check("uncovered residues are black", palette.EMPTY_DARK in fills, True)

    # ... and with the light empty colour the outline has to be there instead
    svg2 = os.path.join(tmp, "stack_light.svg")
    run([COVERAGE] + reps + ["--fasta", fa, "--protein", "P00001",
                             "--layout", "stack", "--empty", "surface",
                             "--out", svg2])
    root2 = ET.parse(svg2).getroot()
    outlined = [r for r in root2.findall(f"{ns}rect")
                if r.get("fill") == palette.EMPTY and r.get("stroke") not in
                (None, "none")]
    check("light empty tiles carry an outline", len(outlined) > 0, True)


# ---------------------------------------------------------------------------
# 10. Sequence-level event decomposition.
#
# The classifier has two judgement calls in it -- the merge threshold that stops
# difflib splitting one change into many, and the point past which multi-event
# rows are pooled -- so both are pinned here against sequences whose answer is
# obvious by eye.
# ---------------------------------------------------------------------------
def case_events(tmp):
    sys.path.insert(0, HERE)
    import fig2b_isoform_strip as ie

    # A pseudorandom sequence, fixed seed. It must not be periodic: repeating a
    # 20-mer gives difflib several equally good alignments for the same edit, so
    # the "right" answer stops being well defined and the test would be pinning
    # a tie-break rather than the classifier.
    import random
    A = "".join(random.Random(11).choices("ACDEFGHIKLMNPQRSTVWY", k=200))
    check("fixture has no repeated 12-mer",
          max(A.count(A[i:i + 12]) for i in range(len(A) - 12)), 1)
    check("identical sequences give no events", ie.events(A, A), [])

    # one clean internal deletion of ten residues
    delA = A[:40] + A[50:]
    check("internal deletion is one event", len(ie.events(A, delA)), 1)
    check("internal deletion is typed and placed",
          ie.event_type(A, delA), "internal segment absent")
    check("deletion coordinates are canonical",
          ie.events(A, delA)[0][1:3], (40, 50))

    # truncations at each end, and an insertion
    check("N-terminal truncation", ie.event_type(A, A[15:]),
          "N-terminal segment absent")
    check("C-terminal truncation", ie.event_type(A, A[:-15]),
          "C-terminal segment absent")
    check("internal insertion", ie.event_type(A, A[:40] + "WWWWWWWW" + A[40:]),
          "internal segment added")

    # two changes far apart stay two; the same two brought within min_equal of
    # each other must merge, which is the whole point of the threshold
    two = A[:20] + A[30:60] + A[70:]
    check("two distant changes stay separate", len(ie.events(A, two)), 2)
    check("two distant changes are labelled", ie.event_type(A, two),
          "2 separate changes")
    # Thresholds are given explicitly here. These assertions are about the merge
    # rule, not about the current default -- MIN_EQUAL is a measured quantity and
    # is expected to move when it is re-scored against the curation.
    near = A[:40] + A[43:47] + A[50:]        # 3 matching residues between them
    check("changes closer than min_equal merge",
          len(ie.events(A, near, min_equal=5)), 1)
    check("changes exactly min_equal apart stay separate",
          len(ie.events(A, near, min_equal=3)), 2)
    # Two deletions three residues apart merge into a *replacement*, not a
    # deletion: the three residues between them survive, so both sides of the
    # merged span are non-empty. That is the honest label for the merged event.
    check("merging two deletions gives a replacement",
          ie.event_type(A, near, min_equal=5), "internal segment replaced")
    # The exact endpoint floats by a residue or two: flank trimming absorbs any
    # chance match at the boundary. Assert containment, not the exact span, or
    # the test pins a coincidence in the fixture.
    cs, ce = ie.events(A, near, min_equal=5)[0][1:3]
    check("the merged span covers both changes", (cs <= 40, ce >= 47),
          (True, True))
    check("a lower threshold does not merge them",
          len(ie.events(A, near, min_equal=2)), 2)
    # align() is the expensive, threshold-independent half; merge() must give
    # exactly what events() gives, or the sweep in prep_uniprot scores
    # something other than what the figures are drawn from.
    for me in (0, 2, 3, 5, 12):
        check(f"align+merge == events at min_equal={me}",
              ie.merge(*ie.align(A, near), min_equal=me),
              ie.events(A, near, min_equal=me))
    for seq in (A[15:], A[:-15], A[:40] + "WWWWWWWW" + A[40:], A, delA):
        check("align+merge == events on " + seq[:6],
              ie.merge(*ie.align(A, seq)), ie.events(A, seq))

    # the tail is pooled rather than given a row each
    many = "".join(A[i:i + 2] for i in range(0, 100, 4))
    check("seven-plus changes pool into one label",
          ie.event_type(A, many, min_equal=1), "6+ separate changes")
    check("five changes are not pooled",
          ie.event_type(A, A[:10] + A[12:25] + A[27:40] + A[42:55] + A[57:70]
                        + A[72:]), "5 separate changes")

    # largest_event picks the biggest change, not the first
    mixed = A[:10] + A[12:40] + A[60:]
    big = ie.largest_event(A, mixed)
    check("largest_event picks the 20-residue change",
          big[2] - big[1], 20)

    # Pinned to the sweep in prep_uniprot.py, not to taste. If this fails,
    # re-run `--sweep` before changing it.
    check("MIN_EQUAL is the scored optimum", ie.MIN_EQUAL, 3)


# ---------------------------------------------------------------------------
# 11. The diverging scale. Two hues plus a neutral midpoint, symmetric arms, and
#     a "no measurement" that is never a colour on the scale.
# ---------------------------------------------------------------------------
def case_diverging(tmp):
    sys.path.insert(0, HERE)
    import lib_palette as palette

    d = palette.diverging
    check("zero is the neutral midpoint", d(0.0), palette.DIVERGING_MID)
    check("small magnitudes stay neutral", d(0.4), palette.DIVERGING_MID)
    check("negative is the blue arm", d(-1.5), palette.DIVERGING_LOW[1])
    check("positive is the red arm", d(1.5), palette.DIVERGING_HIGH[1])
    check("arms saturate at the last step",
          (d(99.0), d(-99.0)),
          (palette.DIVERGING_HIGH[-1], palette.DIVERGING_LOW[-1]))
    check("None means no measurement, not zero", d(None), None)
    check("scale is symmetric in magnitude",
          [palette.DIVERGING_LOW.index(d(-v)) for v in (0.7, 1.5, 2.5, 4.0)],
          [palette.DIVERGING_HIGH.index(d(v)) for v in (0.7, 1.5, 2.5, 4.0)])
    check("arms are equal length",
          len(palette.DIVERGING_LOW), len(palette.DIVERGING_HIGH))
    check("midpoint is neutral, not a hue", max(
        int(palette.DIVERGING_MID[i:i + 2], 16) for i in (1, 3, 5)) - min(
        int(palette.DIVERGING_MID[i:i + 2], 16) for i in (1, 3, 5)) <= 8, True)
    check("no-data colour is not on the scale",
          palette.MISSING in palette.DIVERGING_LOW + palette.DIVERGING_HIGH
          + [palette.DIVERGING_MID], False)

    # Both arms must run monotonically away from the midpoint, or a bigger
    # magnitude can read as a smaller one.
    #
    # Measured in OKLab L, which is what the arms were constructed against and
    # what the dataviz validator checks -- NOT WCAG relative luminance. Those are
    # different quantities: a red and a blue of equal perceptual lightness differ
    # substantially in relative luminance, because the coefficients are 0.2126 for
    # red and 0.0722 for blue. Asserting the wrong one fails a correct palette.
    def oklab_l(h):
        def lin(c):
            c = int(c, 16) / 255
            return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
        r, g, b = (lin(h[i:i + 2]) for i in (1, 3, 5))
        l = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b) ** (1 / 3)
        m = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b) ** (1 / 3)
        s_ = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b) ** (1 / 3)
        return 0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s_

    for name, arm in (("low", palette.DIVERGING_LOW),
                      ("high", palette.DIVERGING_HIGH)):
        ls = [oklab_l(c) for c in arm]
        check(f"{name} arm darkens monotonically",
              all(a > b for a, b in zip(ls, ls[1:])), True)
        check(f"{name} arm steps clear the dL floor",
              all(a - b >= 0.06 for a, b in zip(ls, ls[1:])), True)
    # Matched steps carry equal visual weight, so neither direction shouts.
    for i, (lo, hi) in enumerate(zip(palette.DIVERGING_LOW,
                                     palette.DIVERGING_HIGH)):
        check(f"arms are lightness-matched at step {i + 1}",
              abs(oklab_l(lo) - oklab_l(hi)) < 0.02, True)


# ---------------------------------------------------------------------------
# 13. Run-length merging for the single-line coverage layout.
# ---------------------------------------------------------------------------
def case_runs(tmp):
    sys.path.insert(0, HERE)
    import fig1cd_coverage as pc

    check("empty depth gives no runs", pc.runs([], 4), [])
    check("a uniform stretch is one run", pc.runs([2, 2, 2, 2], 4),
          [(0, 4, 2)])
    check("changes split runs", pc.runs([0, 0, 1, 1, 0], 4),
          [(0, 2, 0), (2, 4, 1), (4, 5, 0)])
    # capping must happen before merging, or 4 and 9 become separate runs that
    # then paint the identical colour and hairline against each other
    check("depths above the cap merge", pc.runs([4, 9, 5], 4), [(0, 3, 4)])
    check("the cap does not merge below itself", pc.runs([3, 4], 4),
          [(0, 1, 3), (1, 2, 4)])

    # the invariants that matter: runs tile the sequence exactly once, in order,
    # and every residue keeps its own capped depth
    import random
    rnd = random.Random(3)
    depth = [rnd.choice([0, 0, 1, 1, 2, 3, 7]) for _ in range(500)]
    rs = pc.runs(depth, 4)
    check("runs start at 0 and end at n", (rs[0][0], rs[-1][1]), (0, 500))
    check("runs are contiguous",
          all(a[1] == b[0] for a, b in zip(rs, rs[1:])), True)
    check("no two adjacent runs share a depth",
          all(a[2] != b[2] for a, b in zip(rs, rs[1:])), True)
    rebuilt = [d for i, j, d in rs for _ in range(j - i)]
    check("runs reconstruct the capped depths",
          rebuilt, [min(d, 4) for d in depth])
    check("merging actually saves rects", len(rs) < len(depth), True)


# ---------------------------------------------------------------------------
# 14. Compact number format, and dropping the x-axis furniture.
# ---------------------------------------------------------------------------
def case_compact(tmp):
    sys.path.insert(0, HERE)
    from fig1b_depth import fmt_compact as f

    check("hundreds of thousands lose the decimal", f(376550), "377K")
    check("tens of thousands keep one", f(14795), "14.8K")
    check("thousands keep one", f(6404), "6.4K")
    check("a trailing .0 is dropped on axes", f(15001, exact=True), "15K")
    check("axis steps read cleanly",
          [f(v, exact=True) for v in (0, 2500, 100000, 400000)],
          ["0", "2.5K", "100K", "400K"])
    check("millions switch suffix",
          (f(1000000, exact=True), f(1250000)), ("1M", "1.25M"))
    # An axis tick at 5,000 *is* 5K, so the trailing zero is noise. A data label
    # is rounded, and dropping the decimal there turns 4,023 into "4K", which
    # reads as an exact round number rather than an approximation of one.
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

    hdr = ["R.FileName", "R.Condition", "PG.ProteinGroups", "PG.Qvalue",
           "PEP.StrippedSequence", "EG.ModifiedSequence", "EG.Qvalue", "FG.Charge"]
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

    def texts(svg):
        return [(e.text or "").strip()
                for e in ET.parse(svg).getroot().iter("{http://www.w3.org/2000/svg}text")]

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
    # dropping the labels must reclaim their space, not leave a white band
    h_on = float(ET.parse(on).getroot().get("height"))
    h_off = float(ET.parse(off).getroot().get("height"))
    check("--no-xticks reclaims the bottom margin", h_off < h_on - 30, True)

    # --text-scale must move the type AND the space around it. Scaling the
    # numbers alone would leave larger labels colliding with the panel edges.
    def sizes(svg):
        import re
        return sorted({float(v) for v in
                       re.findall(r'font-size="([\d.]+)"', open(svg).read())})

    small = os.path.join(tmp, "ts_small.svg")
    big = os.path.join(tmp, "ts_big.svg")
    run([FIGURE, counts, "--out", small, "--units", "peptides",
         "--text-scale", "1.0"])
    run([FIGURE, counts, "--out", big, "--units", "peptides",
         "--text-scale", "2.0"])
    a, b = sizes(small), sizes(big)
    check("text-scale multiplies every size", len(a), len(b))
    check("text-scale doubles them", all(abs(y - 2 * x) < 0.15    # 1-dp rounding
                                        for x, y in zip(a, b)), True)
    ra = ET.parse(small).getroot()
    rb = ET.parse(big).getroot()
    check("bigger type gets a wider left margin",
          float(rb.get("width")) > float(ra.get("width")), True)
    check("bigger type gets a taller canvas",
          float(rb.get("height")) > float(ra.get("height")), True)

    # and the labels must still clear the facet strip they sit under
    import re as _re
    doc = open(big).read()
    strip_bottom = max(float(y) + float(h) for _x, y, _w, h in
                       _re.findall(r'<rect x="([\d.]+)" y="([\d.]+)" '
                                   r'width="([\d.]+)" height="([\d.]+)" '
                                   r'fill="#f0efec"', doc))
    labels = [float(y) for _x, y in
              _re.findall(r'<text x="([\d.]+)" y="([\d.-]+)"\s*'
                          r'font-size="19\.0"', doc)]
    check("median labels clear the facet strip at 2x",
          all(y - 14 > strip_bottom for y in labels) if labels else True, True)

    # --width / --height hold the canvas while the type grows into it, which is
    # the whole point: a figure has a column to fit and the type has to be legible
    # inside it, not at the cost of it.
    fixed = []
    for scale in ("1.0", "1.8", "2.6"):
        f = os.path.join(tmp, f"fix_{scale}.svg")
        run([FIGURE, counts, "--out", f, "--units", "peptides",
             "--width", "700", "--height", "400", "--text-scale", scale])
        r = ET.parse(f).getroot()
        fixed.append((scale, float(r.get("width")), float(r.get("height")), f))
    check("--width/--height pin the canvas",
          {(w, h) for _s, w, h, _f in fixed}, {(700.0, 400.0)})
    a = sizes(fixed[0][3])
    b = sizes(fixed[-1][3])
    check("type still grows inside a pinned canvas",
          all(y > x for x, y in zip(a, b)), True)

    # the headroom for the median label is derived from the label, so it must
    # hold at every scale rather than at the one it was eyeballed on
    import re as _r
    for scale, _w, _h, f in fixed:
        doc = open(f).read()
        sbot = max(float(y) + float(hh) for _x, y, _w2, hh in
                   _r.findall(r'<rect x="([\d.]+)" y="([\d.]+)" width="([\d.]+)" '
                              r'height="([\d.]+)" fill="#f0efec"', doc))
        med = round(9.5 * float(scale), 1)
        ys = [float(y) for y in _r.findall(
            rf'<text x="[\d.]+" y="([\d.-]+)"\s*font-size="{med}" '
            rf'text-anchor="middle"', doc)]
        check(f"median labels clear the strip at scale {scale}",
              all(y - med * 0.78 > sbot for y in ys) if ys else True, True)

    # and an impossible request must say so rather than draw a broken panel
    p = run([FIGURE, counts, "--out", os.path.join(tmp, "toosmall.svg"),
             "--units", "peptides", "--width", "700", "--height", "120",
             "--text-scale", "3.0"], expect_rc=1)
    check("an impossible --height is refused with a reason",
          "shrink the type" in (p.stdout + p.stderr), True)



# ---------------------------------------------------------------------------
# 15. Text must stay inside the canvas -- measured, not asserted.
#
#     The older check looked at x/y *attributes*, which says nothing about where
#     a string actually lands: an end-anchored label at x=103 has an in-bounds
#     attribute and can still be drawn to x=-2. This rasterises and asks the
#     renderer where the glyphs went.
# ---------------------------------------------------------------------------
def case_text_bounds(tmp):
    sys.path.insert(0, HERE)
    try:
        import pymupdf                                    # noqa: F401
    except ImportError:
        print("  skip  pymupdf not installed")
        return
    import re as _r

    def margins(path):
        import pymupdf
        doc = pymupdf.open(path)
        pg = doc[0]
        src = open(path).read()
        vb = [float(v) for v in
              _r.search(r'viewBox="([^"]+)"', src).group(1).split()]
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

    # A real figure, not a fixture: the failure this guards against came from a
    # long label in a column sized for short ones, which a fixture would not have.
    wf = os.path.join(tmp, "wf.svg")
    run([WORKFLOW, "--out", wf, "--letter", "a"])
    m = margins(wf)
    check("workflow text stays inside the canvas",
          m is not None and min(m) >= 0, True)
    # a renderer substituting a wider font must not push it out either
    check("workflow keeps a margin a wider font cannot eat",
          m is not None and min(m) >= 8, True)

    # A label must also stay inside the *box* it sits in, which the canvas check
    # cannot see: "Post-mortem brain" fitted the canvas perfectly while hanging
    # out of both ends of a sample box sized for "PD brain".
    def box_overflow(path):
        import pymupdf
        src = open(path).read()
        vb = [float(v) for v in
              _r.search(r'viewBox="([^"]+)"', src).group(1).split()]
        W, H = vb[2], vb[3]
        # labelled boxes are the workflow's; the sequence band below is drawn in a
        # translated <g>, whose raw coordinates are not where it is rendered
        flat = _r.sub(r'<g transform="translate\([^"]*\)">.*?</g>', "", src,
                      flags=_r.S)
        boxes = [tuple(map(float, m)) for m in
                 _r.findall(r'<rect x="([\d.]+)" y="([\d.]+)" width="([\d.]+)" '
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
                # contained, not merely level with it -- a row label sits beside
                # its colour swatch and must not be judged against it
                if bx < cx < bx + bw and by <= cy <= by + bh:
                    if min(x0 - bx, bx + bw - x1) < 0:
                        out.append(t)
        return out

    check("no label overflows its box", box_overflow(wf), [])
    # and it must still hold when the box's text is longer than the default
    wf2 = os.path.join(tmp, "wf2.svg")
    run([WORKFLOW, "--out", wf2, "--letter", "a",
         "--tissue", "Dorsolateral prefrontal cortex",
         "--groups", "6 Lewy body dementia · 6 control"])
    check("the sample box grows with its longest line", box_overflow(wf2), [])

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
    for scale in ("1.0", "1.7", "2.4"):
        f = os.path.join(tmp, f"tb_{scale}.svg")
        run([FIGURE, counts, "--out", f, "--units", "peptides,precursors",
             "--letter", "b", "--no-xticks", "--width", "900", "--height", "380",
             "--text-scale", scale])
        m = margins(f)
        check(f"count panel text stays inside at text-scale {scale}",
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
    # display and assign must agree: two spellings of one protease share a colour
    # AND a label, or the figure contradicts itself
    col = palette.assign(["GluC", "Glu-C"])
    check("one protease, one colour", len(set(col.values())), 1)


# ---------------------------------------------------------------------------
# 17. The median convention, which is what made two panels of one figure
#     disagree: the rank curve ended at the 7th of 12 samples while the count
#     panel reported the mean of the 6th and 7th.
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

    # a sample that ran out contributes 0, so the tail is not a summary of the
    # one deepest run dressed up as a summary of all of them
    lists = [[0.9, 0.8], [0.9, 0.8, 0.7, 0.6, 0.5]]
    c = pc.median_curve(lists)
    check("short samples pull the tail down", len(c), 4)   # median of 2 and 5
    check("the median is over every sample, absent counted as zero",
          round(c[2], 3), 0.35)
    # ties round half up rather than to even, so the endpoint does not depend on
    # whether the median happens to be odd or even
    check("a .5 median rounds up, not to even",
          (len(pc.median_curve([[1.0] * 2, [1.0] * 3])),
           len(pc.median_curve([[1.0] * 3, [1.0] * 4]))), (3, 4))


def case_splice_events(tmp):
    """GENCODE CDS structures -> rMATS/SUPPA2 event codes."""
    sys.path.insert(0, HERE)
    import prep_splice_events as se

    # canonical: three coding exons; isoform: middle one skipped
    can = [(100, 200), (300, 400), (500, 600)]
    skip = [(100, 200), (500, 600)]
    check("a middle exon dropped is a skipped exon",
          se.classify(can, skip, "+"), {"SE"})
    check("and it is symmetric", se.classify(skip, can, "+"), {"SE"})

    # a shifted boundary is a splice-site change, not a skip. On the + strand
    # the genomic-high end of an exon is its donor, so moving it is A5SS.
    donor = [(100, 200), (300, 380), (500, 600)]
    check("moving the exon's 3' genomic end on + is A5SS",
          se.classify(can, donor, "+"), {"A5SS"})
    acceptor = [(100, 200), (320, 400), (500, 600)]
    check("moving its 5' genomic end on + is A3SS",
          se.classify(can, acceptor, "+"), {"A3SS"})

    # ...and both swap on the minus strand, where transcription runs the other
    # way. Getting this backwards mislabels half the genome and nothing errors.
    check("the same shift on - is A3SS",
          se.classify(can, donor, "-"), {"A3SS"})
    check("and the other is A5SS",
          se.classify(can, acceptor, "-"), {"A5SS"})

    # an intron of one transcript covered by a single exon of the other
    retained = [(100, 200), (300, 600)]
    check("a fused exon pair is a retained intron",
          "RI" in se.classify(can, retained, "+"), True)

    # first/last exon. On + the first transcribed exon is the genomic first.
    alt_first = [(10, 50), (300, 400), (500, 600)]
    check("a non-overlapping first exon on + is AF",
          "AF" in se.classify(can, alt_first, "+"), True)
    check("the same structure on - is AL, not AF",
          "AL" in se.classify(can, alt_first, "-"), True)
    check("and AF is not also claimed on -",
          "AF" in se.classify(can, alt_first, "-"), False)

    # mutually exclusive: each side takes one of a non-overlapping pair, and
    # both are internal so neither is a first/last exon change
    mxe_a = [(100, 200), (300, 340), (500, 600)]
    mxe_b = [(100, 200), (360, 400), (500, 600)]
    check("each transcript taking one of a pair is MXE",
          "MXE" in se.classify(mxe_a, mxe_b, "+"), True)

    check("identical structures give no events", se.classify(can, can, "+"), set())
    check("an empty side gives no events", se.classify(can, [], "+"), set())

    check("introns are the gaps between exons",
          se.introns([(1, 10), (21, 30), (41, 50)]), [(11, 20), (31, 40)])
    check("abutting exons leave no intron",
          se.introns([(1, 10), (11, 20)]), [])

    # the UniProt cross-reference carries the isoform tag inline
    p = os.path.join(tmp, "x.tsv.gz")
    import gzip
    with gzip.open(p, "wt") as fh:
        fh.write("Entry\tEnsembl\n")
        fh.write("Q14BN4\tENST00000449503.6 [Q14BN4-2];ENST00000659705.1 "
                 "[Q14BN4-1];\n")
        fh.write("P99999\tENST00000000001.1;\n")
    x = se.read_xrefs(p)
    check("the bracket tag routes a transcript to its isoform",
          x["Q14BN4-2"], ["ENST00000449503"])
    check("versions are stripped", x["Q14BN4-1"], ["ENST00000659705"])
    check("an untagged transcript belongs to the entry",
          x["P99999"], ["ENST00000000001"])

    # A tie between two equally parsimonious transcript pairs must go the same
    # way on every run. It used to follow set order, which changes with string
    # hashing, and 18-20 isoforms changed event type between runs on identical
    # input. T1 shifts an acceptor, T2 skips an exon: one change each, so the
    # sorted-first pair, T1's A3SS, must win under every hash seed.
    import subprocess
    probe = ("import sys; sys.path.insert(0, %r); import prep_splice_events as se; "
             "ex = {'C': ('chr1', '+', [(1, 10), (20, 30), (40, 50)]), "
             "'T1': ('chr1', '+', [(1, 10), (22, 30), (40, 50)]), "
             "'T2': ('chr1', '+', [(1, 10), (40, 50)])}; "
             "print(sorted(se.best_pair_codes({'C'}, {'T2', 'T1'}, ex)))") % HERE
    seen = {subprocess.run([sys.executable, "-c", probe], capture_output=True,
                           text=True, env={**os.environ, "PYTHONHASHSEED": str(k)}
                           ).stdout.strip() for k in range(8)}
    check("a tie between transcript pairs goes the same way under every hash seed",
          seen, {"['A3SS']"})


def case_pack_hits(tmp):
    """Peptide rows on the isoform strip: packed by overlap, ends exclusive."""
    sys.path.insert(0, HERE)
    import fig2b_isoform_strip as ie

    # abutting products of the same digest: trypsin cuts after K/R, so the next
    # peptide starts at the residue after the previous one ends. They do not
    # overlap and must share a row -- `>` instead of `>=` here is what made
    # trypsin look like it needed two bands.
    hits = [(0, "PEPTIDEK", "Trypsin"), (8, "SAMPLER", "Trypsin")]
    rows, n = ie.pack_hits(hits)
    check("abutting peptides share a row", n, 1)
    check("and are both placed", len(rows), 2)
    check("ends are exclusive", sorted((a, b) for a, b, _d, _t in rows),
          [(0, 8), (8, 15)])

    # one residue of genuine overlap does take a second row
    rows, n = ie.pack_hits([(0, "PEPTIDEK", "Trypsin"), (7, "KSAMPLER", "Trypsin")])
    check("overlapping peptides do not", n, 2)

    # By default the protease is irrelevant to packing: three peptides that do
    # not overlap share one row whichever enzyme made them. Grouping them by
    # digest put non-overlapping bars on separate rows for a reason invisible to
    # a reader, which is what a row on this panel must not mean.
    hits = [(0, "AAAA", "Trypsin"), (40, "CCCC", "GluC"), (80, "DDDD", "LysC")]
    rows, n = ie.pack_hits(hits)
    check("non-overlapping peptides share a row across digests", n, 1)
    check("and all three are placed", len(rows), 3)

    # overlap still forces a new row, and still does so across digests
    rows, n = ie.pack_hits([(0, "AAAAAAAA", "Trypsin"), (4, "CCCCCCCC", "GluC")])
    check("overlap across digests takes a second row", n, 2)

    # `group=True` keeps the old one-band-per-protease form, which the panel no
    # longer uses but which is still the right answer when the question is
    # "what did each enzyme reach"
    hits = [(0, "AAAA", "Trypsin"), (40, "CCCC", "GluC"), (80, "DDDD", "LysC")]
    rows, n = ie.pack_hits(hits, group=True)
    check("group=True gives each digest its own band", n, 3)
    band = {d: t for _a, _b, d, t in rows}
    check("and bands follow the fixed protease order",
          [band[d] for d in ("GluC", "LysC", "Trypsin")], [0, 1, 2])

    # the mapper moves a peptide into drawing coordinates before packing, so an
    # isoform peptide is tested for overlap where it is actually drawn
    rows, _n = ie.pack_hits([(0, "PEPK", "Trypsin")], lambda j: j + 100)
    check("the mapper is applied to the start", rows[0][0], 100)
    check("and the end follows it", rows[0][1], 104)

    # A peptide spanning a junction: its residues land either side of a gap, so
    # its drawn end is the last residue's column, not start + length. Getting
    # this wrong let a second peptide pack onto the same row underneath it, and
    # in browser mode it is the peptide that carries the most evidence on the
    # panel -- the one proving the form has that junction.
    gap = lambda j: j if j < 10 else j + 40          # noqa: E731
    rows, _n = ie.pack_hits([(6, "PEPTIDEK", "Trypsin")], gap)
    check("a junction-spanning peptide reserves the whole span",
          (rows[0][0], rows[0][1]), (6, 54))
    rows, n = ie.pack_hits([(6, "PEPTIDEK", "Trypsin"),
                            (12, "SAMPLER", "Trypsin")], gap)
    check("so a peptide inside that span is not packed under it", n, 2)

    check("no peptides, no rows", ie.pack_hits([])[1], 0)


def case_browser_strip(tmp):
    """The genome-browser encoding: thick where a form has sequence, thin where not.

    Thickness carries the whole message, and one union frame carries all three
    event types: shared columns, then the canonical's own residues, then the
    isoform's own. What can silently break is the frame arithmetic at the two
    ends of that range -- an insertion, where the canonical's own block is empty
    and *it* takes the gap, and a replacement, where both forms have a block and
    the drawn columns therefore belong to neither form's numbering.
    """
    import fig2b_isoform_strip as ie
    from lib_fasta import read_fasta

    can = "M" + "A" * 40 + "C" * 40
    ins = can[:41] + "WWWWWWWWWW" + can[41:]         # 10 residues inserted at 41
    dele = can[:41] + can[51:]                       # 10 residues deleted at 41
    repl = can[:41] + "WWWWW" + can[51:]             # 10 replaced by 5

    check("an insertion is seen as one", ie.largest_event(can, ins)[0], "insert")
    check("a deletion is seen as one", ie.largest_event(can, dele)[0], "delete")
    check("a replacement is seen as one", ie.largest_event(can, repl)[0],
          "replace")

    fa = os.path.join(tmp, "browser.fasta")
    with open(fa, "w") as fh:
        for acc, seq in (("P00001", can), ("P00001-2", ins),
                         ("P00001-3", repl)):
            fh.write(f">sp|{acc}|T_HUMAN test\n{seq}\n")
    check("the fasta round-trips", read_fasta(fa)["P00001-2"], ins)

    def draw(iso, extra=()):
        out = os.path.join(tmp, f"{iso}.svg")
        rc = ie.evidence_panel("P00001", iso, can, read_fasta(fa)[iso],
                               {"Trypsin": {can[30:45]}},
                               {"Trypsin": {read_fasta(fa)[iso][30:60]}},
                               out, ie.FONT, letters=True, browser=True)
        return out, rc

    # --no-residues: the same bars, no letter inside any cell; and an isoform
    # that runs past the canonical's end plus its own block is drawn whole
    # (SPTB-2 lost its last 53 residues to a frame sized on the canonical)
    tail = can[:41] + "WWWWW" + can[51:] + "K" * 20      # replaced, then longer
    bare = os.path.join(tmp, "bare.svg")
    ie.evidence_panel("P00001", "P00001-4", can, tail, {}, {}, bare, ie.FONT,
                      letters=True, browser=True, residues=False,
                      flank_res=500)
    svg = open(bare).read()
    check("--no-residues draws no residue letter",
          re.search(r">[ACMWK]</text>", svg) is None, True)
    check("... and a longer isoform is drawn to its last residue",
          f"1–{len(tail)} of {len(tail)}" in svg, True)
    # --margin / --bar-height: the bars start at the margin and thicken, while
    # the title stays clear of the panel letter, where it was
    wide = os.path.join(tmp, "wide.svg")
    ie.evidence_panel("P00001", "P00001-4", can, tail, {}, {}, wide, ie.FONT,
                      label="GENE", letters=True, browser=True, flank_res=500,
                      margin=80, bar_height=22)
    ws = open(wide).read()
    check("--margin starts the bars at the margin",
          re.search(r'<rect x="80\.0" y="[\d.]+" width="[\d.]+" height="22\.0"',
                    ws) is not None, True)
    check("... and the title does not move with them",
          re.search(r'<text x="58\.0" y="28\.0"[^>]*>GENE<', ws) is not None,
          True)

    out, _ = draw("P00001-2")
    svg = open(out).read()
    check("an insertion draws", os.path.getsize(out) > 0, True)
    check("... and says what changed", "10 aa inserted" in svg, True)
    check("... with no red left on the panel", ie.SEQ_CAN in svg, False)
    # an insertion's frame is the isoform's numbering, so the ruler's total is
    # the isoform's length (RUFY3 read "residues 1–620 of 469")
    check("... and its ruler counts the isoform's residues",
          re.search(rf"residues \d+–\d+ of {len(ins)}<", svg) is not None, True)
    # --no-ruler --gap-numbers: no ruler, and the canonical's residues either
    # side of its gap numbered (the insert follows residue 41, so 41 and 42)
    gn = os.path.join(tmp, "gapnum.svg")
    ie.evidence_panel("P00001", "P00001-2", can, ins, {}, {}, gn, ie.FONT,
                      letters=True, browser=True, residues=False,
                      flank_res=500, ruler=False, gap_numbers=True)
    gs = open(gn).read()
    check("--no-ruler drops the ruler", "residues " in gs, False)
    check("--gap-numbers numbers the canonical either side of its gap",
          (">41<" in gs, ">42<" in gs), (True, True))

    # A replacement, which an aligned frame cannot show: both forms have
    # residues at the same positions. The union frame separates them into two
    # blocks, the way a browser draws mutually exclusive exons.
    out, _ = draw("P00001-3")
    svg = open(out).read()
    check("a replacement draws too", "10 aa replaced by 5" in svg, True)
    check("... with no red on it either", ie.SEQ_CAN in svg, False)
    # and the cost of that frame is stated rather than papered over: the columns
    # are a union of two numberings, so no single residue ruler can be right.
    check("the ruler is dropped for a replacement", "residues " in svg, False)
    # and nothing replaces it. A sentence explaining the frame lived here for
    # one build and read as a title rather than a note about coordinates.
    check("... and nothing is put in its place",
          "aligned on the sequence they share" in svg, False)
    # 40 columns of canonical residues and 40 of isoform, but different ones --
    # the two ranges differ, which is the whole reason a single ruler cannot work
    check("... each band stating its own range", "27\u201366 of 81" in svg, True)
    check("... which is not the same range", "27\u201361 of 76" in svg, True)


def case_sample_ids(tmp):
    """Run name -> patient id, for every naming convention in the cohort.

    This exists because the rule was duplicated. `SAMPLE_RX` (now `lib_report.py`'s) was
    fixed to allow a token between `CF` and the number, and `extra_differential.py`
    kept a private `CF_(\\d+)_` that was never updated -- so one patient's trypsin
    run resolved everywhere except in the module that builds `quant/`, where its
    83,816 rows were dropped with no error and the patient simply had no column.
    The check is that **one** parser answers for all of them.
    """
    sys.path.insert(0, HERE)
    from lib_report import sample_of

    check("the 60-min naming resolves",
          sample_of("2026-03-04_CF_2522_L_new_A6_1_8022"), "2522")
    # the replacement injection for one patient, acquired five months later with
    # a different convention: a token between CF and the id, no separator before
    # the protease letter
    check("a run with a token between CF and the id resolves",
          sample_of("2026-08-27_CF_PD_2522T_G-A1_1_9301"), "2522")
    check("and the 30-min naming too",
          sample_of("2026-08-27_CF_PD_1488G_30min_G-A1_1_9329"), "1488")
    # four digits, not \d+: an unanchored \d+ takes the leading year
    check("the date is not mistaken for a patient id",
          sample_of("2026-03-04_CF_1825_T_new_A6_1_8004"), "1825")
    check("an unparseable name returns None rather than guessing",
          sample_of("blank_A1_1_0001"), None)

    # every module that turns run names into patients must use that one parser
    import re as _re
    src = open(os.path.join(HERE, "extra_differential.py")).read()
    check("extra_differential.py keeps no private copy of the rule",
          bool(_re.search(r'compile\(r"CF', src)), False)


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

    # the second panel is offset by the first's *scaled* height plus the gap,
    # not by its authored height -- getting this wrong overlaps the panels
    check("second panel is offset by the scaled height",
          "translate(0 110.00)" in src, True)
    check("first panel sits at the origin", "translate(0 0.00)" in src, True)
    check("each panel keeps its own units", "scale(2.000000)" in src
          and "scale(4.000000)" in src, True)

    # the font stack is an attribute of the wrapper <svg> that gets stripped, so
    # it has to be carried onto each panel's <g> or every <text> in the composite
    # falls back to the renderer's default serif -- which is what shipped once
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

    # A shared colour key, stamped on the composite rather than owned by a panel.
    # The claim worth testing is not that it draws swatches but that it cannot
    # disagree with the marks: the hexes come from `palette.assign`, so a key
    # naming a protease and a panel drawing it resolve the same colour.
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
    # right-aligned to end inside the canvas: a key that runs off the edge is
    # invisible rather than wrong, which is the harder kind of bug to notice
    xs = [float(m) for m in re.findall(r'<rect x="([\d.]+)" y="[\d.]+" '
                                       r'width="[\d.]+" height="[\d.]+" rx="2"', src3)]
    check("all three swatches are on the canvas", len(xs), 3)
    check("and none runs past the right edge", max(xs) + 9.0 <= 200.0, True)
    check("nor off the left one", min(xs) >= 0.0, True)
    check("the key does not change the height",
          compose.parse(out3)[2], 220.0)

    # A composite composed again (Figure 3c: a --row, lettered on its own). Its
    # root carries no font, since each panel's <g> does, so the letter has to
    # take the stack the groups share or it is stamped with font-family=''.
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

    # A whole figure from its lettered panels (audit.py's manuscript/figures/fig3.svg). Panel
    # la is 200 px of a 100-unit canvas; lb is 100 px of a 50-unit one, the way
    # 3a's pies are narrower than 3b. Both are drawn at 2 px per unit, so
    # --natural must keep both at x2, not stretch lb to 200 px at x4. Their
    # letters are written two different ways, as 1a's and 1b's are.
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
          lets, [("a", "22", "32", f"{18.9 * __import__('lib_palette').TEXT_BOOST:.2f}"),
                 ("b", "22", "132", f"{18.9 * __import__('lib_palette').TEXT_BOOST:.2f}")])
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
    check("without --natural the narrow panel is stretched",
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
# 23. The isoform box columns (Supp. 1): isoform minus canonical within patient.
# ---------------------------------------------------------------------------
def case_box_columns(tmp):
    sys.path.insert(0, HERE)
    import fig2a_supp1_diagnostic_peptides as dp

    # The box panel's contrast is isoform minus canonical *within patient*: a
    # gene that falls as a whole must read Δ 0, as TH does in the PD version,
    # and only an isoform-specific fall may read otherwise.
    cond = {"a": "LBD", "b": "LBD", "c": "LBD",
            "d": "Control", "e": "Control", "f": "Control"}
    med = {(x, "Trypsin"): 0.0 for x in cond}
    down = dict(a=1.0, b=1.1, c=0.9, d=4.0, e=4.2, f=3.8)
    flat = dict(a=4.0, b=4.1, c=3.9, d=4.0, e=4.2, f=3.8)

    def pick(iso):
        idx = {("P1", "Trypsin"): {"CANK": down},
               ("P1-2", "Trypsin"): {"ISOK": iso}}
        return idx, {"canon_acc": "P1", "canon_members": [("P1", "CANK")],
                     "diag_members": [("P1-2", [("P1-2", "ISOK")])]}
    idx, pk = pick(down)
    cols, (delta, _p, na, nb) = dp.box_columns(pk, idx, med, cond,
                                               "LBD", "Control")
    check("a whole-gene fall reads Δ 0 in the box panel", round(delta, 6), 0.0)
    check("... with one value per patient in each column",
          [len(d) for _l, d, _k in cols], [6, 6])
    check("... and patients, not peptides, as n", (na, nb), (3, 3))
    idx, pk = pick(flat)
    _c, (delta, p, _na, _nb) = dp.box_columns(pk, idx, med, cond,
                                              "LBD", "Control")
    check("an isoform spared by the fall reads Δ > 0", delta > 1.5, True)
    check("... and is significant at the patient level", p < 0.01, True)
    # uncentred: the same patients at their absolute level, shifted by
    # exactly the control mean the centred version subtracts
    idx, pk = pick(down)
    raw, _con = dp.box_columns(pk, idx, med, cond, "LBD", "Control",
                               centre=False)
    ctl = sum(math.log2(down[s]) for s in "def") / 3
    rawv = {s: v for _k, s, _d, v in raw[0][1]}
    check("--no-centre shifts every patient by the control mean only",
          {s: round(v - rawv[s], 9) for _k, s, _d, v in cols[0][1]},
          {s: round(-ctl, 9) for s in cond})

    # --precursors: each charge state its own point, placed by its stripped
    # sequence, and the isoform gate still counting peptides, not precursors
    two = {("P1-2", "Trypsin"): {"ISOK|ISOK|2": flat, "ISOK|ISOK|3": down}}
    pts = dp.precursor_values([("P1-2", "ISOK|ISOK|2"),
                               ("P1-2", "ISOK|ISOK|3")],
                              two, med, cond, "Control")
    check("a precursor point per charge state per patient", len(pts), 12)
    idx, pk = pick(flat)
    pc, pcon = dp.box_columns(pk, idx, med, cond, "LBD", "Control",
                              unit="peptide")
    _pc, qcon = dp.box_columns(pk, idx, med, cond, "LBD", "Control")
    check("--unit peptide: a point per peptide per patient, keyed by peptide",
          sorted({k for k, _s, _d, _v in pc[1][1]}), ["ISOK"])
    check("... and the same patient-level Δ as --unit patient",
          round(pcon[0], 9), round(qcon[0], 9))
    # --unit run: one point per patient per digest, its peptides averaged
    runs = {("P1", "Trypsin"): {"AK": flat, "CK": down},
            ("P1", "GluC"): {"DE": flat}}
    medr = {**med, **{(x, "GluC"): 0.0 for x in cond}}
    rv = dp.run_values([("P1", "AK"), ("P1", "CK"), ("P1", "DE")],
                       runs, medr, cond, "Control")
    check("--unit run: a point per patient per digest", len(rv), 12)
    pv = {(k, s, d, round(v, 9)) for k, s, d, v in dp.precursor_values(
        [("P1", "AK"), ("P1", "CK")], runs, medr, cond, "Control")}
    want = round(sum(v for _k, s, d, v in pv if s == "a") / 2, 9)
    check("... each the mean of that run's centred peptides",
          [round(v, 9) for _k, s, d, v in rv if s == "a" and d == "Trypsin"],
          [want])
    check("... each centred on its own control mean",
          round(sum(v for k, s, _d, v in pts
                    if k.endswith("|3") and cond[s] == "Control"), 9), 0.0)
    st = dp.column([("P1-2", "ISOK|ISOK|2"), ("P1-2", "ISOK|ISOK|3")],
                   two, med, cond, "LBD", "Control", 3)
    check("two charge states of one peptide gate as one peptide",
          st["npep"], 1)
    seqs = {"P1": "MAAAKCANK", "P1-2": "MAAAKISOK"}
    _b, dg = dp.diagnostic("G1", {"P1-2"}, seqs,
                           {"P1": "G1", "P1-2": "G1"}, two)
    check("a precursor key is placed by its stripped sequence",
          {w: len(m) for w, m in dg.items()}, {("P1-2",): 2})
    # --unit fc with --per-row: one box per column, wrapped rows, and the
    # canvas grows to hold them rather than cutting the second row off
    fcb = [(f"G{i}", [("P1", [("AK", None, "Trypsin", 0.1 * i)] * 3, True),
                      ("P1-2", [("BK", None, "GluC", 0.5)] * 3, False)],
            (0.4, 0.01, 0, 0)) for i in range(3)]
    heights = []
    for pr in (None, 2):
        out = os.path.join(tmp, f"fc{pr}.svg")
        dp.box_panel(fcb, {}, out, "sans-serif", unit="fc", per_row=pr,
                     panel_h=300)
        root = ET.parse(out).getroot()
        heights.append(float(root.get("viewBox").split()[3]))
    check("--per-row 2 puts three facets on two rows, on a taller canvas",
          heights[1] > heights[0] + 176, True)
    # volcano: every row a point, survivors labelled, and the BH line at the
    # weakest survivor's p, so it separates survivors from the rest exactly
    vr = [{"gene": f"G{i}", "log2fc": (-1) ** i * 0.1 * i, "p": 10 ** -i,
           "q": 10 ** -i * 3} for i in range(1, 7)]
    out = os.path.join(tmp, "volcano.svg")
    dp.volcano(vr, out, "sans-serif")
    root = ET.parse(out).getroot()
    ns = "{http://www.w3.org/2000/svg}"
    texts = {t.text for t in root.iter(f"{ns}text")}
    check("the volcano draws a point per isoform",
          len(list(root.iter(f"{ns}circle"))), 6)
    check("... and labels only the q < 0.05 survivors",
          {f"G{i}" for i in range(1, 7)} & texts, {"G2", "G3", "G4", "G5", "G6"})
    # --unit protein: one column per gene, a Control and an LBD box, and the
    # facet as wide as a two-column facet (same column width everywhere)
    pb = [(f"G{i}", [("P1", [(s, s, None, (1.0 if cond[s] == "LBD" else 0.0)
                             + 0.1 * j) for j, s in enumerate(sorted(cond))],
                      True)], (1.0, 0.001, 3, 3)) for i in range(2)]
    out = os.path.join(tmp, "protein.svg")
    dp.box_panel(pb, cond, out, "sans-serif", unit="protein", per_row=2,
                 subtitle=False)
    ps = open(out).read()
    check("--unit protein prints the LBD - Control effect without a Δ",
          ("+1.00  ·  p = 0.001" in ps, "Δ +1.00" in ps), (True, False))
    check("... and no accession under its single column", ">P1<" in ps, False)
    check("a Tukey whisker stops at the last point inside 1.5 IQR",
          dp._quartiles([0.0, 1.0, 2.0, 3.0, 100.0])[4], 3.0)


def case_diagnostic_peptides(tmp):
    sys.path.insert(0, HERE)
    import fig2a_supp1_diagnostic_peptides as dp

    # "diagnostic" must mean absent from the canonical, not merely belonging to
    # a group that names an isoform. A group like `P1;P1-2` is reported *because*
    # its peptides cannot discriminate, and the first version of the screen put
    # exactly those columns at the top of the ranking.
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

    # a gene with no single canonical entry has nothing for "absent from the
    # canonical" to mean, and must say so rather than guessing a reference
    b2, d2 = dp.diagnostic("G2", {"Q1-2"}, {"Q1-2": canon},
                           {"Q1-2": "G2"}, {})
    check("an isoform-only gene is refused", (b2, d2), (None, {}))

    # --- the p-value's replication unit ------------------------------------
    # Three units, and the panel must never quote an effect from one beside a
    # p-value from another: that mismatch is what a reader checks first.
    cond = {"a": "LBD", "b": "LBD", "c": "LBD",
            "d": "Control", "e": "Control", "f": "Control"}
    med2 = {(x, "Trypsin"): 0.0 for x in cond}
    # one peptide up 2x in cases, one flat; every patient has both
    idx2 = {("P1-2", "Trypsin"): {
        "AAAK": dict(a=4.0, b=4.0, c=4.0, d=1.0, e=1.0, f=1.0),
        "CCCK": dict(a=2.0, b=2.0, c=2.0, d=2.0, e=2.0, f=2.0)}}
    mem2 = [("P1-2", "AAAK"), ("P1-2", "CCCK")]
    got = {}
    for u in ("patients", "patient-ratio", "peptides"):
        st = dp.column(mem2, idx2, med2, cond, "LBD", "Control", 3, u)
        got[u] = st
    check("the peptide unit counts peptides",
          (got["peptides"]["unit"], got["peptides"]["n_unit"]),
          ("peptides", 2))
    check("... and the patient units count patients",
          (got["patients"]["n_unit"], got["patient-ratio"]["n_unit"]),
          ("3 v 3", "3 v 3"))
    # AAAK is +2, CCCK is 0, so the dots average +1 and that is what the
    # peptide unit must quote -- not the aggregate contrast
    check("the peptide unit's effect is the mean of the drawn dots",
          round(got["peptides"]["mean"], 6),
          round(sum(v for _p, _d, v in got["peptides"]["points"])
                / len(got["peptides"]["points"]), 6))
    check("the ratio unit's effect comes from its own rollup",
          round(got["patient-ratio"]["mean"], 6), 1.0)

    # With nothing missing all three agree, which is why the summed rollup went
    # unnoticed. Drop one peptide from one case patient -- routine in DIA -- and
    # the summed unit moves while the ratio unit does not, because a sum over
    # one peptide is not comparable with a sum over two. That is the defect
    # `patient_ratios` exists for, and on the real data it is worth 0.010 -> 0.066.
    gap = {("P1-2", "Trypsin"): {
        "AAAK": dict(a=4.0, b=4.0, c=4.0, d=1.0, e=1.0, f=1.0),
        "CCCK": dict(b=2.0, c=2.0, d=2.0, e=2.0, f=2.0)}}
    summed = dp.column(mem2, gap, med2, cond, "LBD", "Control", 3, "patients")
    ratio = dp.column(mem2, gap, med2, cond, "LBD", "Control", 3,
                      "patient-ratio")
    check("one missing value moves the summed rollup",
          round(summed["mean"], 3) == round(got["patients"]["mean"], 3), False)
    # ... and it moves the ratio rollup too, which is worth asserting rather
    # than hoping: ratios remove the *scale* artefact (a sum over one peptide
    # is not smaller than a sum over two) but not the *composition* one (a
    # patient averaged over a different subset of peptides still differs).
    # Neither rollup is a fix for missingness; only a model that estimates
    # peptide effects is, and this repo does not have one.
    check("... and the ratio rollup is NOT immune either",
          round(ratio["mean"], 6) == round(got["patient-ratio"]["mean"], 6),
          False)

    # the jitter is deterministic, or the panel is not reproducible
    check("jitter is stable across calls",
          dp._jitter("PEPTIDEK", "GluC", 10.0),
          dp._jitter("PEPTIDEK", "GluC", 10.0))
    check("... and differs between digests of one peptide",
          dp._jitter("PEPTIDEK", "GluC", 10.0)
          == dp._jitter("PEPTIDEK", "LysC", 10.0), False)
    check("... and stays inside the column",
          abs(dp._jitter("PEPTIDEK", "GluC", 10.0)) <= 5.0, True)

    # the subtitle and footnote are wrapped before the canvas is sized, so the
    # wrapper has to respect its budget -- an unwrapped line is a label drawn
    # off the viewBox, which MuPDF does not clip and so `preview.py` hides
    lines = dp._wrap("one point per peptide per protease and then some more "
                     "words to force a wrap", 24)
    check("every wrapped line fits the budget",
          max(len(x) for x in lines) <= 24, True)
    check("and no word is lost",
          " ".join(lines).split(),
          ("one point per peptide per protease and then some more words to "
           "force a wrap").split())
    check("a label longer than the budget is still emitted",
          dp._wrap("unbreakablesinglewordlongerthanbudget", 5),
          ["unbreakablesinglewordlongerthanbudget"])
    check("the column label names the extras rather than listing them",
          (dp._label("P1-2"), dp._label("P1-2;P1-5")),
          ("P1-2", "P1-2 (+1)"))

    # --- a facet with two contrasts carries each under its own column ---------
    import lib_palette as palette
    can = {"points": [("A", "GluC", 0.1), ("B", "LysC", -0.1)], "mean": 0.0,
           "p": float("nan"), "unit": "peptides"}
    d1 = {"points": [("C", "GluC", 0.5)], "mean": 0.5, "p": 0.25, "delta": True,
          "unit": "peptides"}
    d2 = {"points": [("D", "Trypsin", -0.4)], "mean": -0.4, "p": 0.5,
          "delta": True, "unit": "peptides"}

    def drawn(per_column):
        path = os.path.join(tmp, f"dots_{per_column}.svg")
        dp.panel([{"gene": "G1", "canon_acc": "P1", "canon": can,
                   "diag": [("P1-2", d1), ("P1 pS5", d2)],
                   "per_column": per_column}], path, palette.FONT, "",
                 key=False, subtitle=False, hue=True, panel_h=300)
        return [t.text for t in ET.parse(path).getroot().iter(
            "{http://www.w3.org/2000/svg}text") if t.text and
            (t.text.startswith("Δ") or t.text.startswith("p ="))]
    check("a facet quotes its first contrast once, over the facet",
          drawn(False), ["Δ +0.50  ·  p = 0.250"])
    check("... or every contrast under its own column",
          drawn(True), ["Δ +0.50", "p = 0.250", "Δ -0.40", "p = 0.500"])


def case_reported_groups(tmp):
    """The reference reconstruction: columns are group strings, not sequence.

    Two things here were wrong in the first build and neither raised anything.
    A peptide seen by two digests was listed twice in the column's members, so
    it was drawn twice -- 91 points where the data holds 74. And the console
    table quoted `column()["mean"]`, a digest-averaged protein-level contrast,
    beside a bar drawn at the points' median: two different quantities, one
    caption.

    The third check is the load-bearing one. When no digest measured both
    columns the patient-level contrast does not exist, and saying so is the
    whole point of the figure -- returning some number instead would turn "this
    comparison cannot be made" into "this comparison was not significant".
    """
    sys.path.insert(0, HERE)
    import extra_reported_groups as rg

    check("a lone canonical accession is the base column",
          rg.is_canonical_group("P04899"), True)
    check("... a group naming two is not", rg.is_canonical_group("P1;P1-2"),
          False)
    check("... nor is a lone isoform", rg.is_canonical_group("P04899-4"), False)

    # one peptide, seen by two digests: two points, one member
    st = {"points": [("PEPK", "Trypsin", 0.5), ("PEPK", "LysC", 0.4),
                     ("SAMPLER", "Trypsin", -0.2)]}
    check("the split is counted off the drawn points",
          rg.digest_split(st), {"GluC": 0, "LysC": 1, "Trypsin": 2})

    cond = {"a": "LBD", "b": "LBD", "c": "LBD",
            "d": "Control", "e": "Control", "f": "Control"}
    med = {(s, d): 0.0 for s in cond for d in ("Trypsin", "GluC")}
    # the GNAI2 shape: each column measured by one digest, and a different one
    idx = {("P1", "Trypsin"): {"AAAK": {s: 4.0 for s in cond}},
           ("P1;P1-2", "GluC"): {"CCCE": {s: 8.0 for s in cond}}}
    p_pat, shared = rg.patient_contrast("P1", "P1;P1-2", idx, med, cond,
                                        "LBD", "Control")
    check("no digest measured both columns", shared, [])
    check("... so the patient-level p is not a number", p_pat == p_pat, False)

    # and when one digest does see both, it is computable
    idx[("P1;P1-2", "Trypsin")] = {"CCCE": {s: 8.0 for s in cond}}
    _pt, shared = rg.patient_contrast("P1", "P1;P1-2", idx, med, cond,
                                      "LBD", "Control")
    check("one digest seeing both makes it computable", shared, ["Trypsin"])


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

    # the same canonical fields come out of all three, which is the point
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
    # the regression: half of trypsin's peptides end in K, so a 0.5 threshold
    # plus a smallest-residue-set tie-break called the real trypsin report LysC
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
        # exactly the shape the 30-min re-exports carry: text, capital E
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


def case_phospho(tmp):
    sys.path.insert(0, HERE)
    import lib_palette as palette
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
    # P1 = MKKAESPVKEEAVAEK: S6 has the phosphate, T is absent, so the long
    # tryptic-style peptide KAESPVKEEAVAEK (0-based start 2) has one S/T/Y.
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
        # a Glu-C-only placement on the second patient, in a dropped-run digest
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

    # -- the cohort, its null and its pairs ----------------------------------
    known = {"P1": {6: "experimental"}}
    co = ph.Cohort(rows, known, seqs, ["1111"])
    check("the cohort keeps the complete patients only", co.patients, ["1111"])
    sid = next(i for i, s in co.S.items() if "P1:6" in s["keys"])
    check("a curated residue is curated", co.curated[sid], True)
    # KKAESPVKE, KAESPVK and AESPVKEEAVAEK each hold one S/T/Y, the curated S6
    check("with one candidate per peptide the null is certainty",
          co.null[sid], 1.0)
    check("the residue mix of random placement is all S",
          dict(co.random_residues(sid)), {"S": 1.0, "T": 0.0, "Y": 0.0})
    check("a pair whose shared stretch holds one S/T/Y is forced, so dropped",
          co.pairs(), [])

    two = [dict(r) for r in rows if r["sample"] == "1111"]
    two[0] = dict(two[0], peptide="SKKAESPVKE", offsets=[5],
                  where=[("P9", 0)])
    two[2] = dict(two[2], peptide="SKKAESPVKE", offsets=[0],
                  where=[("P9", 0)])
    for r in two:
        r["placements"] = [(frozenset(f"{a}:{i + o + 1}" for a, i in r["where"]),
                            r["peptide"][o]) for o in r["offsets"]]
    co2 = ph.Cohort(two, {}, dict(seqs, P9="SKKAESPVKE"), ["1111"])
    check("two digests on one peptide, two candidates, different residues",
          co2.pairs(), [(False, 2, "P9")])
    check("a forced pair is kept only when asked for",
          sorted(p[1] for p in co.pairs(min_candidates=1)), [1, 1, 1])
    check("an uncurated protein has a zero null", set(co2.null.values()), {0.0})

    # -- panels: drawn at the shared width, dashes drawn as segments ----------
    out = os.path.join(tmp, "phospho_figs")
    ph.panel_gain(co, {}, out, palette.FONT, "a")
    ph.panel_curated(co, out, palette.FONT, "b")
    ph.panel_example(co, "P1", 6, "TEST", out, palette.FONT, "c")
    ph.panel_residues(co, out, palette.FONT)
    widths = {os.path.basename(p): float(ET.parse(p).getroot()
                                         .get("viewBox").split()[2])
              for p in glob.glob(os.path.join(out, "*.svg"))}
    check("every phospho panel is PANEL_W wide", set(widths.values()),
          {palette.PANEL_W})
    check("no panel relies on stroke-dasharray, which MuPDF drops",
          any("dasharray" in open(os.path.join(out, f)).read() for f in widths),
          False)

    # -- differential: a site that falls exactly with its protein -------------
    import collections
    import math
    pats = ["1001", "1002", "1003", "1004", "1005", "1006"]
    cond = dict(zip(pats, ["LBD"] * 3 + ["Control"] * 3))
    qdir = os.path.join(tmp, "phospho_quant")
    os.makedirs(qdir, exist_ok=True)
    # protein P1 is 4x lower in LBD; the other groups are flat filler
    groups = {"P1": [1000.0] * 3 + [4000.0] * 3}
    groups.update({f"F{i}": [500.0 + i] * 6 for i in range(8)})
    write_report(os.path.join(qdir, "Trypsin.csv"), ["protein_group"] + pats,
                 [[g] + v for g, v in groups.items()], delim=",")
    frames = ph.protein_frames(qdir, max_mad=0)
    check("protein frames keep every run that clears the QC",
          frames["Trypsin"]["keep"], pats)
    prec = []
    for p, cnd in cond.items():
        site = 100.0 if cnd == "LBD" else 400.0          # 4x down, like its protein
        prec.append(["Trypsin", f"CF_{p}_T", p, "_AAS[Phospho (STY)]PK_",
                     "AASPK", 2, 1e-4, 1e-4, site, 1, "P1", "2", "P1:0"])
    # a second peptide for the same site, seen in fewer patients
    prec.append(["Trypsin", "CF_1001_T", "1001", "_KAAS[Phospho (STY)]PK_",
                 "KAASPK", 2, 1e-4, 1e-4, 9e9, 1, "P1", "3", "P1:-1"])
    dscan = os.path.join(tmp, "phospho_da_scan")
    os.makedirs(dscan, exist_ok=True)
    write_report(os.path.join(dscan, "phospho_runs.tsv"),
                 ["digest", "run", "sample", "precursors", "phospho_precursors",
                  "kept"], [["Trypsin", f"CF_{p}_T", p, 5000, 2, 1] for p in pats])
    write_report(os.path.join(dscan, "phospho_precursors.tsv"), head, prec)
    drows = ph.load(dscan)
    dS = ph.sites(drows)
    qty, ref = ph.site_quantities(drows, dS)
    sid = next(i for i, s in dS.items() if "P1:3" in s["keys"])
    check("the reference peptide is the one placed in most patients",
          ref["Trypsin"][sid], "AASPK")
    check("... so the other peptide's quantity is not summed in",
          max(qty["Trypsin"][sid].values()), 400.0)
    D = ph.Differential(drows, dS, frames, cond)
    res = D.run()
    lv = next(r for r in res[("level", "Trypsin")] if r["group"] == sid)
    oc = next(r for r in res[("occupancy", "Trypsin")] if r["group"] == sid)
    check("a site on a 4x-lost protein falls 2 log2 units",
          round(lv["log2fc"], 6), -2.0)
    check("... and its occupancy does not change",
          round(oc["log2fc"], 6), 0.0)
    check("the protein it is scaled by is its own group", D.prot["Trypsin"][sid], "P1")
    perms = D.permutations(15)
    check("shuffles are balanced",
          {tuple(sorted(collections.Counter(p.values()).values())) for p in perms},
          {(3, 3)})
    check("... distinct, and never the real labels",
          (len({tuple(sorted(p.items())) for p in perms}),
           any(p == cond for p in perms)), (15, False))
    check("a permutation p counts the real value among its shuffles",
          ph.perm_p(0.5, [0.1, 0.5, 0.9, 0.2]), 3 / 5)
    x = {"level": {"GluC": {"a": 1.0, "b": 2.0, "c": 3.0}}}
    agree = ph.agreement({("level", "GluC"): [{"group": g, "log2fc": v}
                                              for g, v in x["level"]["GluC"].items()],
                          ("level", "LysC"): [{"group": g, "log2fc": -v}
                                              for g, v in x["level"]["GluC"].items()]},
                         {})
    check("agreement is the correlation of shared fold changes",
          round(agree[("level", "GluC", "LysC")][1], 6), -1.0)

    # -- isoform-specific sites: own residue, and whether it survives
    #    localisation. S3 is shared; T7 S8 T9 are the isoform's own.
    iseqs = {"P7": "MAASAKGSGGAKLLLLK", "P7-2": "MAASAKGTSTAKLLLLK"}
    import fig2cd_isoform_coverage as ic
    # the PD grid's bins agree with the percent a reader would quote: 16 of 63
    # residues is "25%" and binned 1-25%; only complete coverage takes 100%
    from fig3bc_phospho_atlas import NONE_FILL, STEPS as CELL
    check("grid_fill: none, then 1-25 / 26-50 / 51-75 / 76-99 / 100%, on rounded %",
          [ic.grid_fill(f) for f in (0.0, 0.001, 16 / 63, 0.255, 0.5, 0.51,
                                     0.75, 0.76, 0.996, 1.0)],
          [NONE_FILL, CELL[0], CELL[0], CELL[1], CELL[1], CELL[2], CELL[2],
           CELL[3], CELL[3], CELL[4]])
    # the PD grid's junction cells count peptides that span a join: AAAK|BBBK
    # joins at 4, so AKBB spans it, and a peptide ending or starting at the join
    # (AAAK, BBBK) does not
    check("junction_peptides: only a peptide holding both sides of the join spans it",
          ic.junction_peptides("AAAKBBBK", [4],
                               {"LysC": {"AAAK", "BBBK"}, "Trypsin": {"AKBB"}}),
          {"GluC": set(), "LysC": set(), "Trypsin": {"AKBB"}})
    check("count_fill: black for none, 3a's purple ramp for 1 to 5+",
          (ic.count_fill(0), ic.count_fill(5) == ic.count_fill(9),
           ic.count_fill(1) != ic.count_fill(2)), (NONE_FILL, True, True))
    check("the isoform's own residues are the replaced ones",
          sorted(ic.discriminating(iseqs["P7"], iseqs["P7-2"])), [7, 8, 9])
    check("each side's own residues come from one alignment",
          tuple(sorted(x) for x in ic.own_residues("MAAGGGKLL", "MAAXKLL")),
          ([3, 4, 5], [3]))
    # -- junctions (experimental 2c): a pure deletion's only evidence is a
    #    peptide spanning the join, in one piece
    jc, ji = "MKTAYIAKQRWWWWWWPEPTIDEK", "MKTAYIAKQRPEPTIDEK"
    check("a pure deletion has no own residue but one junction",
          (ic.discriminating(jc, ji), ic.junctions(jc, ji)), (set(), [10]))
    check("a terminal deletion is not a junction",
          ic.junctions("MKTAYIAKQRWWW", "MKTAYIAKQR"), [])
    jcov = ic.junction_coverage(ji, [10], {"Trypsin": {"QRPEPT", "AKQR"},
                                           "LysC": {"AKQR", "PEPTIDEK"}})
    check("a peptide across the join covers it; two meeting at it do not",
          (jcov["Trypsin"], jcov["LysC"]), ({10}, set()))
    iso_prec = [
        ["LysC", "CF_1001_L", "1001", "_AS[Phospho (STY)]AKGTSTAK_", "ASAKGTSTAK",
         2, 1e-4, 1e-4, 100, 1, "P7-2", "1", "P7-2:2"],
        ["LysC", "CF_1002_L", "1002", "_ASAKGT[Phospho (STY)]STAK_", "ASAKGTSTAK",
         2, 1e-4, 1e-4, 100, 1, "P7-2", "5", "P7-2:2"],
    ]
    iscan = os.path.join(tmp, "phospho_iso_scan")
    os.makedirs(iscan, exist_ok=True)
    write_report(os.path.join(iscan, "phospho_runs.tsv"),
                 ["digest", "run", "sample", "precursors", "phospho_precursors",
                  "kept"], [["LysC", f"CF_{p}_L", p, 5000, 1, 1]
                            for p in ("1001", "1002", "1003")])
    write_report(os.path.join(iscan, "phospho_precursors.tsv"), head, iso_prec)
    irows = ph.load(iscan)
    found = {s["entries"][0]: s for s in ph.isoform_sites(ph.sites(irows), irows,
                                                           iseqs)}
    check("a phosphate on a shared residue is not isoform-specific",
          found["P7-2:4"]["own"], False)
    check("a phosphate on the isoform's own residue is",
          found["P7-2:8"]["own"], True)
    check("... but one phosphate beside a shared S/T/Y is not localisation-proof",
          found["P7-2:8"]["proof"], False)
    iso_prec.append(["LysC", "CF_1003_L", "1003",
                     "_AS[Phospho (STY)]AKGT[Phospho (STY)]STAK_", "ASAKGTSTAK",
                     2, 1e-4, 1e-4, 100, 1, "P7-2", "1,5", "P7-2:2"])
    write_report(os.path.join(iscan, "phospho_precursors.tsv"), head, iso_prec)
    irows = ph.load(iscan)
    found = {s["entries"][0]: s for s in ph.isoform_sites(ph.sites(irows), irows,
                                                           iseqs)}
    check("two phosphates and one shared S/T/Y put one on the isoform for sure",
          found["P7-2:8"]["proof"], True)
    check("classes count isoform-specific sites by the digests placing them",
          [(n, k, v) for n, k, v, _ in ph.iso_classes(list(found.values()))],
          [("Lys-C only", "LysC", 1)])
    both = ph.pair_sites(ph.sites(irows), iseqs,
                         {"P7": {"LysC": {"ASAKGTSTAK", "ASAKGSGGAK"}}}, {})
    check("a pair panel can draw a site on one form's own residues when both "
          "forms have peptides of their own",
          [(r["isoform"], r["position"], r["on"], r["own_residues"],
            r["trypsin_sees_both"]) for r in both],
          [("P7-2", 8, "isoform", 3, False)])
    check("... and nothing when one form has no peptide of its own",
          ph.pair_sites(ph.sites(irows), iseqs, {"P7": {"LysC": {"ASAKGTSTAK"}}},
                        {}), [])
    # -- the strip: Figure 2a's isoform panel, with the phosphosites ---------
    import fig2b_isoform_strip as ie
    svg = "{http://www.w3.org/2000/svg}"
    # GFAPε's shape: one replacement at the C-terminus that a chance DGE splits
    gc = "M" + "A" * 20 + "ETSLDTKSVSE" + "DGE" + "VIKESK"
    gi = "M" + "A" * 20 + "GGKSTK" + "DGE" + "NHKVTRYLKSLTIR"
    check("two changes a chance match splits are drawn as one",
          (len(ie.events(gc, gi)), ie.whole_event(gc, gi)),
          (2, ("replace", 21, 41, 21, 44)))
    check("... which the largest alone would halve",
          ie.largest_event(gc, gi)[1:3], (35, 41))

    def strip_svg(name, **kw):
        path = os.path.join(out, name)
        ie.evidence_panel("P9", "P9-2", gc, gi, {"Trypsin": {"ETSLDTK"}},
                          {"LysC": {"AAGGK"}}, path, palette.FONT, "TEST",
                          letters=True, browser=True,
                          event=ie.whole_event(gc, gi), **kw)
        return ET.parse(path).getroot()

    plain = strip_svg("strip_plain.svg")
    letters = "".join(t.text for t in plain.iter(svg + "text")
                      if t.text and len(t.text) == 1)
    check("an isoform's own block past the canonical's end is drawn in full",
          letters.endswith("GGKSTKDGENHKVTRYLKSLTIR"), True)
    marked = strip_svg("strip_marked.svg",
                       marks={"canonical": {23: {"GluC", "Trypsin"}}},
                       marks_key="the key")
    nested = os.path.join(out, "strip_nested.svg")
    ie.evidence_panel("P9", "P9-2", gc, gi, {"Trypsin": {"ETSLDTK", "TSLDTK"}},
                      {"LysC": {"AAGGK"}}, nested, palette.FONT, "TEST",
                      letters=True, browser=True, event=ie.whole_event(gc, gi),
                      dedupe=True)
    tryp = palette.assign(ph.ORDER)["Trypsin"]
    check("dedupe leaves off a bar nested in a longer one of the same digest",
          sum(r.get("fill") == tryp for r in ET.parse(nested).getroot().iter(
              svg + "rect")), 1)
    check("the region label is there by default, and can be left off",
          ("discriminating region" in open(os.path.join(out, "strip_plain.svg")).read(),
           "discriminating region" in ET.tostring(
               strip_svg("strip_bare.svg", region_label=False)).decode()),
          (True, False))
    colour = palette.assign(ph.ORDER)
    dots = [el.get("fill") for el in marked.iter(svg + "circle")]
    texts = [t.text for t in marked.iter(svg + "text")]
    check("a site gets a dot per digest that placed it, and its residue",
          (sorted(dots[:2]) == sorted([colour["GluC"], colour["Trypsin"]]),
           "S24" in texts, "the key" in texts), (True, True, True))
    check("... and the room for them only when there are sites",
          (len(list(plain.iter(svg + "circle"))),
           float(marked.get("viewBox").split()[3]) >
           float(plain.get("viewBox").split()[3])), (0, True))
    marks, pts = ph.site_marks(irows, "P7", "P7-2",
                               ie.whole_event(iseqs["P7"], iseqs["P7-2"]))
    check("the strip marks a phosphate on a form's own residues, not a shared one",
          marks, {"isoform": {7: {"LysC"}}})
    check("... with the patients whose precursor placed it",
          sorted(pts[("isoform", 7, "LysC")]), ["1002", "1003"])

    # -- which sites 2c's dot plots draw -------------------------------------
    can_pts = [("A", "GluC", 0.0), ("B", "LysC", 0.2), ("C", "Trypsin", -0.2)]
    below = [("p1", "GluC", -1.0), ("p2", "LysC", -0.8), ("p3", "Trypsin", -0.1)]
    check("a column all below the canonical's mean is on one side",
          ph.one_side(can_pts, below), True)
    check("... but not with a point above it, whatever the column's mean",
          ph.one_side(can_pts, below[:2] + [("p3", "Trypsin", 0.9)]), False)
    check("... nor with only two points",
          ph.one_side(can_pts, below[:2]), False)

    def facet(gene, site, delta, census=True, side=True):
        return {"gene": gene, "acc": site.split(":")[0],
                "label": "pS" + site.split(":")[1], "st": {"mean": delta},
                "census": census, "one_side": side}
    check("the figure takes census sites on one side, largest |Δ|, one per gene",
          ph.census_picks([facet("A", "P1:5", -0.4), facet("B", "P2:9", 1.2),
                           facet("A", "P1:7", -0.9), facet("C", "P3:1", 2.0,
                                                           census=False),
                           facet("D", "P4:2", 1.5, side=False)], n=2),
          ["P2:9", "P1:7"])

    # -- a pure deletion: the fixture the census below uses --------------------
    cseqs = {"P8": "MKAAKGSGPEKLLLLRK", "P8-2": "MLLLLRK"}
    cscan = os.path.join(tmp, "phospho_canon_scan")
    os.makedirs(cscan, exist_ok=True)
    write_report(os.path.join(cscan, "phospho_runs.tsv"),
                 ["digest", "run", "sample", "precursors", "phospho_precursors",
                  "kept"], [["LysC", "CF_1001_L", "1001", 5000, 1, 1]])
    write_report(os.path.join(cscan, "phospho_precursors.tsv"), head,
                 [["LysC", "CF_1001_L", "1001", "_GS[Phospho (STY)]GPEK_", "GSGPEK",
                   2, 1e-4, 1e-4, 100, 1, "P8", "1", "P8:5"]])
    crows = ph.load(cscan)

    # -- the census of sites in a stretch that tells two forms apart ---------
    ev = os.path.join(tmp, "events.tsv")
    write_report(ev, ["isoform", "events", "note"], [["P8-2", "SE", ""]])
    census = ph.region_sites(ph.sites(crows), crows, cseqs, ev,
                             {"P8": {7: "experimental"}})
    check("a site in a canonical exon an isoform skips is counted",
          [(r["kind"], r["against"], r["position"], r["curated"]) for r in census],
          [("canonical-exon", "P8-2", 7, True)])
    write_report(ev, ["isoform", "events", "note"], [["P8-2", "AF", ""]])
    check("... but not when the isoform only starts elsewhere",
          ph.region_sites(ph.sites(crows), crows, cseqs, ev, {}), [])

    # -- export: fig3a_phospho_sites.py reads the same sites back ------------------------
    import fig3a_phospho_sites as mod_sites
    exp = os.path.join(tmp, "sites_phospho.tsv")
    ph.export_sites(S, exp)
    back = {(r["acc"], r["res"]): r for r in mod_sites.load(exp)}
    check("mod_sites reads every exported site", len(back), len(S))
    check("... with each digest's peptides for the three-digest site",
          back[("P1", 6)]["by_digest"],
          {"GluC": ["KKAESPVKE"], "LysC": ["KAESPVK"], "Trypsin": ["AESPVKEEAVAEK"]})

    # -- gained: Figure 3a's order and running counts, countable by hand -------
    # Trypsin places 3 sites, Lys-C 2 and Glu-C 2. Lys-C and Glu-C each add one
    # site trypsin lacks, so they tie on the total at step 2; Lys-C wins it by
    # confirming two trypsin sites where Glu-C confirms one.
    toy = [{"mod": "P", "digests": d} for d in
           (["Trypsin", "LysC"], ["Trypsin", "LysC", "GluC"], ["Trypsin"],
            ["LysC"], ["GluC"])]
    # -- peptide_classes: the experimental 3a's per-digest stacks -----------------
    # A digest's bar counts only its own peptides; All counts every digest's, so a
    # site on 2 trypsin and 1 Lys-C peptide is a 2 for trypsin, a 1 for Lys-C and
    # a 3 for All. Six peptides cap at 5+.
    stack = [{"by_digest": {"Trypsin": ["A", "B"], "LysC": ["C"]},
              "peps": ["A", "B", "C"]},
             {"by_digest": {"Trypsin": ["D"]}, "peps": ["D"]},
             {"by_digest": {"GluC": list("EFGHIJ")}, "peps": list("EFGHIJ")}]
    check("peptide_classes: a digest counts only its own peptides",
          [dict(mod_sites.peptide_classes(stack, d)) for d in ("Trypsin", "LysC", "GluC")],
          [{2: 1, 1: 1}, {1: 1}, {5: 1}])
    check("peptide_classes: All counts every digest's, capped at 5+",
          dict(mod_sites.peptide_classes(stack)), {3: 1, 1: 1, 5: 1})
    # -- log_height: the grouped experiment's y axis ---------------------------
    check("log_height: 1 sits on the baseline, 10**decades at the top, 100 halfway",
          [mod_sites.log_height(v, 4, 230.0) for v in (1, 10 ** 4, 100)],
          [0.0, 230.0, 115.0])
    check("log_height: an empty class draws nothing rather than failing",
          mod_sites.log_height(0, 4, 230.0), 0.0)
    # -- venn: the experimental 3a's regions and where their counts go ----------
    vs = [{"digests": d} for d in (["Trypsin"], ["Trypsin", "LysC"],
                                    ["Trypsin", "LysC", "GluC"], ["LysC"],
                                    ["LysC"])]
    reg = mod_sites.venn_regions(vs, ["Trypsin", "LysC", "GluC"])
    check("venn_regions: every combination is a key, empty ones as 0",
          (len(reg), reg[frozenset(["GluC"])], reg[frozenset(["LysC"])],
           reg[frozenset(["Trypsin", "LysC", "GluC"])]), (7, 0, 2, 1))
    check("venn_regions: the regions partition the sites",
          sum(reg.values()), len(vs))
    import math
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
    check("gained: most sites first, ties broken by the most confirmed",
          mod_sites.gained(toy),
          [("Trypsin", 3, 0), ("LysC", 4, 2), ("GluC", 5, 2)])


def main():
    with tempfile.TemporaryDirectory() as tmp:
        for fn in (case_semantics, case_header_variants, case_unresolved,
                   case_missing_columns, case_figure, case_figure_groups, case_workflow,
                   case_coverage, case_ramp_contrast, case_events,
                   case_diverging, case_runs, case_compact,
                   case_text_bounds, case_display_names,
                   case_median, case_splice_events, case_pack_hits,
                   case_browser_strip,
                   case_sample_ids, case_compose, case_diagnostic_peptides,
                   case_box_columns,
                   case_reported_groups,
                   case_report, case_phospho):
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
