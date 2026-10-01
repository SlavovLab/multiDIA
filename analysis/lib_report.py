#!/usr/bin/env python3
"""Identify a search report's software, columns and protease from its content.

    python3 lib_report.py 'data/search/*.parquet'          # what am I looking at?
"""

import argparse
import collections
import glob
import os
import re
import sys


ORDER = ["GluC", "LysC", "Trypsin"]


# Primary diagnosis (NPDX1) -> condition; anything else stops the read.
CONDITION = {"Lewy body disease (LBD/DLB)": "LBD", "Control": "Control"}


def read_metadata(path):
    """-> {sample id: condition}, from `data/metadata.xlsx`."""
    from openpyxl import load_workbook
    rows = load_workbook(path, read_only=True, data_only=True).active.iter_rows(
        values_only=True)
    head = [str(h).strip() if h is not None else "" for h in next(rows)]
    try:
        i_id, i_dx = head.index("ADRC #"), head.index("NPDX1")
    except ValueError:
        sys.exit(f"{path}: needs 'ADRC #' and 'NPDX1' columns, has {head}")
    out = {}
    for r in rows:
        sid, dx = r[i_id], r[i_dx]
        if sid is None or dx is None:
            continue
        dx = str(dx).strip()
        if dx not in CONDITION:
            sys.exit(f"{path}: case {sid} has diagnosis {dx!r}, which maps to "
                     f"no condition; add it to report.CONDITION")
        out[str(sid).strip()] = CONDITION[dx]
    return out


SAMPLE_RX = r"CF[_-](?:[A-Za-z]{1,3}[_-])?(\d{4})"


def sample_of(run):
    """The patient id out of a run name. -> "2744" or None."""
    m = re.search(SAMPLE_RX, run)
    return m.group(1) if m else None


def _norm(s):
    """`EG.Qvalue`, `EG_Qvalue`, `eg qvalue` -> `egqvalue`."""
    return re.sub(r"[^a-z0-9]", "", s.strip().lstrip("\ufeff").lower())


# Canonical field -> aliases, best first.
FIELDS = {
    "run":            ["R.FileName", "File.Name", "Run", "Raw file",
                       "Spectrum File"],
    "protein_groups": ["PG.ProteinGroups", "Protein.Group", "Proteins",
                       "Protein ID", "Protein"],
    "peptide":        ["PEP.StrippedSequence", "Stripped.Sequence", "Sequence",
                       "Peptide"],
    "modified":       ["EG.ModifiedSequence", "EG.ModifiedPeptide",
                       "Modified.Sequence", "Modified sequence",
                       "Modified Peptide"],
    "precursor_q":    ["EG.Qvalue", "Q.Value", "PEP"],
    "protein_q":      ["PG.Qvalue", "PG.Q.Value", "Protein.Q.Value"],
    "protein_q_run":  ["PG.QValue (Run-Wise)", "PG.Q.Value"],
    "quantity":       ["FG.Quantity", "Precursor.Quantity",
                       "Precursor.Normalised", "Intensity"],
    "quantity_normalised": ["Precursor.Normalised"],
    "quantity_ms1_apex": ["Ms1.Apex.Area"],
    "protein_qty":    ["PG.Quantity", "PG.MaxLFQ", "Genes.MaxLFQ"],
    "total_qty":      ["EG.TotalQuantity (Settings)",
                       "EG.TargetQuantity (Settings)", "FG.Quantity"],
    "use_for_pep":    ["EG.UsedForPeptideQuantity"],
    "use_for_prot":   ["EG.UsedForProteinGroupQuantity",
                       "PEP.UsedForProteinGroupQuantity"],
    "charge":         ["FG.Charge", "Precursor.Charge", "Charge"],
    "mass":           ["FG.Mass", "Mass"],
    "mz":             ["FG.PrecMz", "Precursor.Mz", "m/z"],
    "missed":         ["PEP.NrOfMissedCleavages", "Missed cleavages"],
    "rt":             ["EG.ApexRT", "RT", "Retention time", "Retention"],
    "decoy":          ["EG.IsDecoy", "Reverse"],
    "proteotypic":    ["PEP.IsProteotypic", "Proteotypic"],
    "pep_score":      ["EG.PEP", "Posterior.Error.Probability"],
    "cscore":         ["EG.Cscore", "CScore"],
}

NUMERIC = {"precursor_q", "protein_q", "protein_q_run", "quantity",
           "quantity_normalised", "quantity_ms1_apex",
           "protein_qty", "total_qty", "charge", "missed", "rt",
           "mass", "mz", "pep_score", "cscore"}
BOOLEAN = {"use_for_pep", "use_for_prot", "decoy", "proteotypic"}

# Columns only one tool writes; identification is a majority vote.
SIGNATURES = {
    "Spectronaut": ["PEP.StrippedSequence", "EG.Qvalue", "FG.Quantity",
                    "PG.ProteinGroups", "R.FileName",
                    "EG.UsedForPeptideQuantity", "EG.ApexRT"],
    "DIA-NN":      ["Precursor.Id", "Stripped.Sequence", "Q.Value",
                    "Protein.Group", "File.Name", "Precursor.Quantity",
                    "RT"],
    "MaxQuant":    ["Raw file", "Modified sequence", "Leading razor protein",
                    "Missed cleavages", "Intensity", "Reverse"],
    "FragPipe":    ["Spectrum", "Peptide", "Protein ID", "Hyperscore",
                    "Retention", "Charge"],
}

# (residues, terminus) per enzyme.
ENZYMES = {
    "GluC":         ("E",    "C"),
    "LysC":         ("K",    "C"),
    "Trypsin":      ("KR",   "C"),
    "ArgC":         ("R",    "C"),
    "Chymotrypsin": ("FWYL", "C"),
    "AspN":         ("D",    "N"),
    "LysN":         ("K",    "N"),
}


def identify(names):
    """Which software wrote a report with these columns. -> (name, evidence)."""
    have = {_norm(n) for n in names}
    best, ev = None, {"matched": 0, "of": 0, "score": 0.0}
    for tool, sig in SIGNATURES.items():
        hit = sum(1 for c in sig if _norm(c) in have)
        score = hit / len(sig)
        if score > ev["score"]:
            best, ev = tool, {"matched": hit, "of": len(sig), "score": score}
    return (best, ev) if ev["score"] >= 0.5 else (None, ev)


def resolve(names):
    """-> {canonical field: the column that actually carries it}."""
    have = {_norm(n): n for n in names}
    out = {}
    for field, aliases in FIELDS.items():
        for a in aliases:
            if _norm(a) in have:
                out[field] = have[_norm(a)]
                break
    return out


def protease(peptides, min_frac=0.80, min_residue_frac=0.10):
    """Infer the enzyme from cleavage specificity. -> (name, evidence)."""
    term = {"C": collections.Counter(), "N": collections.Counter()}
    n = 0
    for p in peptides:
        if not p:
            continue
        n += 1
        term["C"][p[-1]] += 1
        term["N"][p[0]] += 1
    if not n:
        return None, {"n": 0}
    fits = {}
    for name, (residues, side) in ENZYMES.items():
        c = term[side]
        per = {r: c[r] / n for r in residues}
        total = sum(per.values())
        if total >= min_frac and all(v >= min_residue_frac
                                     for v in per.values()):
            fits[name] = (total, per)
    if not fits:
        best = max(ENZYMES, key=lambda e: sum(
            term[ENZYMES[e][1]][r] for r in ENZYMES[e][0]) / n)
        return None, {"n": n, "closest": best,
                      "termini": dict(term["C"].most_common(5))}
    # best explained first; a smaller residue set only breaks a tie
    name = min(fits, key=lambda e: (-round(fits[e][0], 3), len(ENZYMES[e][0])))
    total, per = fits[name]
    return name, {"n": n, "frac": total,
                  "per_residue": {k: round(v, 4) for k, v in per.items()},
                  "also_fits": sorted(set(fits) - {name})}


class Report:
    """One search report, addressed by canonical field name."""

    def __init__(self, path, declared=None, peptide_sample=200000):
        import pyarrow.parquet as pq
        self.path = path
        self._pq = pq
        self._pf = pq.ParquetFile(path)
        self.schema = self._pf.schema_arrow
        self.rows = self._pf.metadata.num_rows
        self.software, self.software_evidence = identify(self.schema.names)
        self.columns = resolve(self.schema.names)
        self.declared = declared
        self.protease, self.protease_evidence = self._protease(peptide_sample)

    def _protease(self, sample):
        if "peptide" not in self.columns:
            return self.declared, {"n": 0, "note": "no peptide column"}
        got = []
        for b in self._pf.iter_batches(batch_size=sample,
                                       columns=[self.columns["peptide"]]):
            got = b.column(0).to_pylist()
            break
        return protease(set(got))

    def has(self, *fields):
        return all(f in self.columns for f in fields)

    def require(self, *fields):
        missing = [f for f in fields if f not in self.columns]
        if missing:
            sys.exit(f"{os.path.basename(self.path)}: no column carries "
                     f"{', '.join(missing)} "
                     f"(software: {self.software or 'unrecognised'})")

    def batches(self, fields, batch_size=262144):
        """Yield {canonical field: list}, numeric fields coerced to numbers."""
        import pyarrow as pa
        import pyarrow.compute as pc
        self.require(*fields)
        cols = [self.columns[f] for f in fields]
        for b in self._pf.iter_batches(batch_size=batch_size, columns=cols):
            out = {}
            for f in fields:
                col = b.column(cols.index(self.columns[f]))
                if f in NUMERIC and pa.types.is_string(col.type):
                    col = pc.cast(col, pa.float64())
                elif f in BOOLEAN and pa.types.is_string(col.type):
                    col = pc.not_equal(pc.utf8_lower(col), "false")
                if f in NUMERIC and pa.types.is_floating(col.type):
                    # NaN -> null, so `q is None or q > cut` filters rows with no q-value
                    col = pc.if_else(pc.is_nan(col),
                                     pa.scalar(None, col.type), col)
                out[f] = col.to_pylist()
            out["_n"] = b.num_rows
            yield out


def runs_of(path):
    """The run names a report contains, parquet or delimited text."""
    if path.endswith(".parquet"):
        import pyarrow.compute as pc
        import pyarrow.parquet as pq
        col = resolve(pq.ParquetFile(path).schema_arrow.names).get("run")
        if col is None:
            return set()
        return set(pc.unique(pq.read_table(path, columns=[col]).column(0))
                   .to_pylist()) - {None}
    import csv
    import gzip
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        delim = "," if path.endswith((".csv", ".csv.gz")) else "\t"
        rows = csv.reader(fh, delimiter=delim)
        header = next(rows, [])
        col = resolve(header).get("run")
        if col is None:
            return set()
        i = header.index(col)
        return {r[i] for r in rows if len(r) > i and r[i]}


def one_search_per_run(paths):
    """Stop if two reports contain the same run."""
    owner, clash = {}, collections.defaultdict(set)
    for p in paths:
        for run in runs_of(p):
            if run in owner and owner[run] != p:
                clash[(owner[run], p)].add(run)
            owner.setdefault(run, p)
    if clash:
        lines = [f"  {os.path.basename(a)} and {os.path.basename(b)}: "
                 f"{len(runs)} run(s) in both, e.g. {sorted(runs)[0]}"
                 for (a, b), runs in sorted(clash.items())]
        sys.exit("reports share runs -- two searches of the same injections, "
                 "and pooling them counts every precursor twice:\n"
                 + "\n".join(lines)
                 + "\n  pass one search, e.g. 'data/search/*-60min-Phospho.parquet'")


def open_reports(paths, declared=None, order=None):
    """-> [Report], one per path, protease inferred and checked."""
    one_search_per_run(paths)
    reps = []
    for p in paths:
        dec = None
        if declared:
            dec = declared.get(p) or declared.get(os.path.basename(p))
        r = Report(p, declared=dec)
        if dec and r.protease and r.protease != dec:
            sys.exit(f"{os.path.basename(p)}: declared {dec} but the peptides "
                     f"say {r.protease} "
                     f"({r.protease_evidence.get('per_residue')})")
        if r.protease is None:
            r.protease = dec
        if r.protease is None:
            sys.exit(f"{os.path.basename(p)}: cannot tell which protease this "
                     f"is from its peptides ({r.protease_evidence}); pass it "
                     f"explicitly")
        reps.append(r)
    if order:
        reps.sort(key=lambda r: (order.index(r.protease)
                                 if r.protease in order else len(order)))
    dup = [k for k, v in collections.Counter(r.protease
                                             for r in reps).items() if v > 1]
    if dup:
        print(f"  note: {', '.join(dup)} appears in more than one report")
    return reps


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reports", nargs="+")
    ap.add_argument("--fields", action="store_true",
                    help="also list the canonical field mapping")
    args = ap.parse_args(argv)
    paths = sorted(p for pat in args.reports for p in glob.glob(pat))
    if not paths:
        sys.exit("no report matched")
    for p in paths:
        r = Report(p)
        ev = r.protease_evidence
        print(f"\n{os.path.basename(p)}")
        print(f"  software   {r.software or 'UNRECOGNISED'}  "
              f"({r.software_evidence['matched']}/"
              f"{r.software_evidence['of']} signature columns)")
        print(f"  rows       {r.rows:,}   columns {len(r.schema.names)}")
        print(f"  protease   {r.protease or 'UNDETERMINED'}"
              + (f"   {ev.get('frac', 0) * 100:.1f}% of {ev.get('n', 0):,} "
                 f"peptides, {ev.get('per_residue')}"
                 if r.protease else f"   {ev}"))
        missing = [f for f in FIELDS if f not in r.columns]
        print(f"  fields     {len(r.columns)}/{len(FIELDS)} resolved"
              + (f", missing: {', '.join(missing)}" if missing else ""))
        strings = [f for f in NUMERIC if f in r.columns
                   and str(r.schema.field(r.columns[f]).type) == "string"]
        if strings:
            print(f"  coerced    {', '.join(strings)} stored as text, "
                  f"cast on read")
        if args.fields:
            for f in FIELDS:
                a = r.columns.get(f)
                print(f"      {f:<16}{a or '-':<32}"
                      f"{r.schema.field(a).type if a else ''}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
