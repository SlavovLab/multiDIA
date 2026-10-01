#!/usr/bin/env python3
"""Write the SDRF-Proteomics sample table for the PRIDE submission.

    python3 prep_sdrf.py 'data/search/*-60min-Phospho.parquet' \\
        --metadata data/metadata.xlsx --out derived/multiDIA.sdrf.tsv
"""

import argparse
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_report import CONDITION, ORDER, open_reports, read_cases, runs_of, sample_of  # noqa: E402

ORGANISM_PART = "dorsolateral striatum"
INSTRUMENT = "NT=timsTOF Ultra 2;AC=MS:1003412"
ACQUISITION = "NT=data-independent acquisition;AC=MS:1003215"
DISSOCIATION = "NT=collision-induced dissociation;AC=MS:1000133"
LABEL = "label free sample"
CLEAVAGE = {"Trypsin": "NT=Trypsin;AC=MS:1001251",
            "LysC": "NT=Lys-C;AC=MS:1001309",
            "GluC": "NT=glutamyl endopeptidase;AC=MS:1001917"}
MODIFICATIONS = ["NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed",
                 "NT=Oxidation;AC=UNIMOD:35;TA=M;MT=Variable",
                 "NT=Acetyl;AC=UNIMOD:1;PP=Protein N-term;MT=Variable",
                 "NT=Phospho;AC=UNIMOD:21;TA=S,T,Y;MT=Variable"]
DISEASE = {"LBD": "Lewy body disease", "Control": "normal"}

COLUMNS = (["source name", "characteristics[organism]", "characteristics[organism part]",
            "characteristics[disease]", "characteristics[cell type]",
            "characteristics[developmental stage]", "characteristics[individual]",
            "characteristics[age]", "characteristics[sex]",
            "characteristics[biological replicate]", "material type", "assay name",
            "technology type", "comment[proteomics data acquisition method]",
            "comment[label]", "comment[instrument]", "comment[cleavage agent details]"]
           + ["comment[modification parameters]"] * len(MODIFICATIONS)
           + ["comment[dissociation method]", "comment[fraction identifier]",
              "comment[technical replicate]", "comment[data file]", "factor value[disease]"])


def age(value):
    """Years as SDRF wants them; the de-identified '≥90' becomes the range 90Y-."""
    v = str(value).strip()
    if v.startswith(("≥", ">=")):
        return v.lstrip("≥>=").strip() + "Y-"
    return f"{int(float(v))}Y"


def rows(paths, cases):
    out = []
    for rep in open_reports(paths, order=ORDER):
        for run in sorted(runs_of(rep.path)):
            sid = sample_of(run)
            if sid not in cases:
                sys.exit(f"{run}: patient {sid} is not in the metadata")
            case = cases[sid]
            disease = DISEASE[CONDITION[str(case["NPDX1"]).strip()]]
            out.append([f"CF_{sid}", "homo sapiens", ORGANISM_PART, disease,
                        "not applicable", "adult", sid, age(case["Age at Death"]),
                        str(case["Sex"]).strip().lower(), "1", "tissue", run,
                        "proteomic profiling by mass spectrometry", ACQUISITION, LABEL,
                        INSTRUMENT, CLEAVAGE[rep.protease]] + MODIFICATIONS
                       + [DISSOCIATION, "1", "1", f"{run}.d", disease])
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reports", nargs="+")
    ap.add_argument("--metadata", default="data/metadata.xlsx")
    ap.add_argument("--out", default="derived/multiDIA.sdrf.tsv")
    args = ap.parse_args(argv)
    paths = sorted(p for pat in args.reports for p in (glob.glob(pat) or [pat]))
    table = rows(paths, read_cases(args.metadata))
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        w.writerows(table)
    print(f"  wrote {args.out}  ({len(table)} runs)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
