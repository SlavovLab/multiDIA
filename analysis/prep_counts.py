#!/usr/bin/env python3
"""Count identifications per run, and per-sample protease unions, from Spectronaut reports.

    python3 prep_counts.py Report.tsv --out derived/counts/
    python3 prep_counts.py 'reports/*.tsv.gz' --dry-run
"""

import argparse
import csv
import glob
import gzip
import hashlib
import io
import json
import os
import re
import statistics
import sys
from collections import defaultdict

# field -> header aliases; the first one present wins
ALIASES = {
    "run": ["R.FileName", "R.Raw File Name", "R.RawFileName", "R.Replicate Name",
            "R.ReplicateName", "R.Label", "Run", "R.Run"],
    "condition": ["R.Condition", "Condition"],
    "replicate": ["R.Replicate", "Replicate"],
    "mod_seq": ["EG.ModifiedSequence", "EG.ModifiedPeptide", "FG.LabeledSequence",
                "EG.PrecursorId", "FG.Id", "Modified.Sequence"],
    "strip_seq": ["PEP.StrippedSequence", "EG.StrippedSequence", "PEP.Sequence",
                  "Stripped.Sequence"],
    "charge": ["FG.Charge", "EG.PrecursorCharge", "Precursor.Charge"],
    "precursor_id": ["EG.PrecursorId", "FG.Id"],
    "protein_group": ["PG.ProteinGroups", "PG.ProteinAccessions", "PG.ProteinGroup",
                      "Protein.Group", "PG.Genes"],
    "eg_qvalue": ["EG.Qvalue", "EG.QValue", "EG.Q.Value", "Q.Value"],
    "pg_qvalue": ["PG.QValue (Run-Wise)", "PG.Qvalue (Run-Wise)", "PG.Qvalue",
                  "PG.QValue", "PG.Q.Value", "Protein.Q.Value"],
    "decoy": ["EG.IsDecoy", "PG.IsDecoy", "IsDecoy", "Decoy"],
    "identified": ["EG.Identified", "EG.IsIdentified"],
}

UNITS = ["precursors", "peptides", "protein_isoform_groups", "proteins_canonical"]

# protease tokens looked for in the run name / condition
PROTEASES = {
    "Trypsin": ["trypsin", "tryp", "_try_", "tryps"],
    "LysC": ["lysc", "lys-c", "lys_c", "lysargc"],
    "GluC": ["gluc", "glu-c", "glu_c", "gluce"],
    "ArgC": ["argc", "arg-c"],
    "AspN": ["aspn", "asp-n"],
    "LysN": ["lysn", "lys-n"],
    "chymotrypsin": ["chymotrypsin", "chymo"],
    "elastase": ["elastase"],
    "pepsin": ["pepsin"],
    "proalanase": ["proalanase", "pro-alanase"],
}

ISOFORM_SUFFIX = re.compile(r"-\d+$")


def norm(name):
    """Fold a header cell for case- and punctuation-insensitive matching."""
    from lib_report import _norm
    return _norm(name)


def resolve_columns(header):
    """-> ({field: index}, {field: header name actually used})."""
    lookup = {}
    for i, cell in enumerate(header):
        lookup.setdefault(norm(cell), i)
    idx, used = {}, {}
    for field, names in ALIASES.items():
        for cand in names:
            if norm(cand) in lookup:
                idx[field] = lookup[norm(cand)]
                used[field] = header[idx[field]]
                break
    return idx, used


def open_text(path):
    if path.endswith((".gz", ".gzip")):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace", newline="")


def sniff_delimiter(path):
    with open_text(path) as fh:
        first = fh.readline()
    counts = {d: first.count(d) for d in ("\t", ",", ";")}
    return max(counts, key=counts.get) if max(counts.values()) else "\t"


def read_rows(path, wanted=None):
    """-> (header, row iterator, file handle or None) for a delimited or parquet file."""
    if not path.endswith(".parquet"):
        delim = sniff_delimiter(path)
        fh = open_text(path)
        reader = csv.reader(fh, delimiter=delim)
        try:
            header = next(reader)
        except StopIteration:
            fh.close()
            sys.exit(f"{path}: empty file")
        return header, reader, fh

    try:
        import pyarrow.parquet as pq
    except ImportError:
        sys.exit(f"{path}: reading parquet needs pyarrow, which is not importable.\n"
                 "  Run with the project environment (.venv/bin/python) or "
                 "convert to TSV first.")

    pf = pq.ParquetFile(path)
    header = [f.name for f in pf.schema_arrow]
    cols = [c for c in header if c in wanted] if wanted else header
    order = [header.index(c) for c in cols]

    def rows():
        for batch in pf.iter_batches(batch_size=65_536, columns=cols):
            data = [batch.column(i).to_pylist() for i in range(batch.num_columns)]
            names = batch.schema.names
            pos = {n: i for i, n in enumerate(names)}
            picked = [pos[c] for c in cols]
            for r in range(batch.num_rows):
                row = [""] * len(header)
                for c, p in zip(order, picked):
                    v = data[p][r]
                    row[c] = "" if v is None else (
                        str(v) if not isinstance(v, bool) else str(v))
                yield row

    return header, rows(), None


def parse_float(s):
    """Parse a float, allowing locale-comma decimals; None if blank or invalid."""
    if s is None:
        return None
    s = s.strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        try:
            return float(s.replace(".", "").replace(",", "."))
        except ValueError:
            return None


def is_true(s):
    return s.strip().lower() in ("true", "1", "yes")


def detect_protease(*fields):
    """Find a protease token in any of the given strings. Ambiguity is an error."""
    blob = " ".join(f.lower() for f in fields if f)
    hits = {name for name, toks in PROTEASES.items()
            if any(t.strip("_") in blob for t in toks)}
    if len(hits) == 1:
        return hits.pop()
    if len(hits) > 1:
        return "AMBIGUOUS:" + "+".join(sorted(hits))
    return None


def strip_protease_token(run):
    """Remove the protease token so what is left identifies the sample."""
    out = run
    for toks in PROTEASES.values():
        for t in toks:
            out = re.sub(re.escape(t.strip("_")), "", out, flags=re.I)
    out = re.sub(r"[_\-.]{2,}", "_", out).strip("_-. ")
    return out or run


def derive(run, condition, replicate, prot_expr, samp_expr):
    """Return (sample, protease) for one run name."""
    if prot_expr:
        m = re.search(prot_expr, run) or (re.search(prot_expr, condition or ""))
        protease = (m.group(1) if m and m.groups() else m.group(0)) if m else None
    else:
        protease = detect_protease(condition, run)
    if samp_expr:
        m = re.search(samp_expr, run)
        sample = (m.group(1) if m and m.groups() else m.group(0)) if m else None
    else:
        sample = (replicate or "").strip() or strip_protease_token(run)
    return sample, protease


def key64(s):
    """64-bit digest of `s`, for memory-light set membership."""
    return int.from_bytes(hashlib.blake2b(s.encode("utf-8"), digest_size=8).digest(), "big")


def canonical_accessions(group_field):
    """'sp|P49189-3|..;P04899' -> {'P49189', 'P04899'} (isoform suffix dropped)."""
    out = set()
    for part in re.split(r"[;,]", group_field):
        part = part.strip()
        if not part:
            continue
        if "|" in part:  # sp|ACC|NAME
            bits = [b for b in part.split("|") if b]
            part = bits[1] if len(bits) >= 2 else bits[0]
        out.add(ISOFORM_SUFFIX.sub("", part))
    return out


def isoform_groups(group_field):
    """The protein group as reported, normalised only for separator whitespace."""
    parts = [p.strip() for p in re.split(r"[;,]", group_field) if p.strip()]
    return {";".join(sorted(parts))} if parts else set()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reports", nargs="+", help="Spectronaut report file(s); globs allowed")
    ap.add_argument("--out", default="derived/counts", help="output directory (default: counts)")
    ap.add_argument("--precursor-q", type=float, default=0.01, help="EG.Qvalue cutoff")
    ap.add_argument("--protein-q", type=float, default=0.01, help="PG.Qvalue cutoff")
    ap.add_argument("--protease-regex", default=None,
                    help="run-name regex, group 1 = protease")
    ap.add_argument("--sample-regex", default=None,
                    help="run-name regex, group 1 = sample")
    ap.add_argument("--min-run-count", type=int, default=0,
                    help="drop runs with fewer precursors than this")
    ap.add_argument("--min-union-digests", type=int, default=1,
                    help="min usable digests for a sample's 'All' union")
    ap.add_argument("--strict", action="store_true",
                    help="require both q-value cutoffs for every unit")
    ap.add_argument("--exact", action="store_true",
                    help="hold identity strings instead of 64-bit digests")
    ap.add_argument("--dry-run", action="store_true",
                    help="print resolved columns and runs, then stop")
    ap.add_argument("--max-rows", type=int, default=None, help="stop after N data rows (smoke test)")
    args = ap.parse_args(argv)

    paths = []
    for pat in args.reports:
        hits = sorted(glob.glob(pat))
        paths.extend(hits if hits else [pat])
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        sys.exit("no such report file(s): " + ", ".join(missing))
    from lib_report import one_search_per_run
    one_search_per_run(paths)

    # {(run, unit): set-of-keys}
    seen = defaultdict(set)
    rows_read = 0
    rows_kept = 0
    dropped = defaultdict(int)
    run_meta = {}          # run -> (condition, replicate, source file)
    used_cols = {}
    per_file = {}

    for path in paths:
        # probe the header first so parquet decodes only the needed columns
        probe, _, fh0 = read_rows(path)
        if fh0 is not None:
            fh0.close()
        idx0, _ = resolve_columns(probe)
        wanted = {probe[i] for i in idx0.values()}
        header, reader, fh = read_rows(path, wanted)
        try:
            idx, used = resolve_columns(header)
            if not used_cols:
                used_cols = used
            required = ["run", "protein_group"]
            have_prec = "mod_seq" in idx or ("strip_seq" in idx and "charge" in idx) \
                or "precursor_id" in idx
            lack = [r for r in required if r not in idx]
            if lack or not have_prec:
                sys.exit(f"{path}: cannot find required columns.\n"
                         f"  missing: {lack or ''}"
                         f"{'' if have_prec else '  no precursor identity column'}\n"
                         f"  header seen: {header[:25]}{' ...' if len(header) > 25 else ''}\n"
                         f"  resolved: {json.dumps(used, indent=2)}")

            file_rows = 0
            for row in reader:
                rows_read += 1
                file_rows += 1
                if args.max_rows and rows_read >= args.max_rows:
                    break
                if len(row) <= max(idx.values()):
                    dropped["short row"] += 1
                    continue
                if "decoy" in idx and is_true(row[idx["decoy"]]):
                    dropped["decoy"] += 1
                    continue
                if "identified" in idx and row[idx["identified"]].strip() and \
                        not is_true(row[idx["identified"]]):
                    dropped["not identified"] += 1
                    continue
                # each unit is filtered by its own q-value unless --strict
                pass_prec = pass_prot = True
                if "eg_qvalue" in idx:
                    q = parse_float(row[idx["eg_qvalue"]])
                    pass_prec = q is not None and q <= args.precursor_q
                if "pg_qvalue" in idx:
                    q = parse_float(row[idx["pg_qvalue"]])
                    pass_prot = q is not None and q <= args.protein_q
                if args.strict:
                    pass_prec = pass_prot = (pass_prec and pass_prot)
                if not (pass_prec or pass_prot):
                    dropped["failed both q cutoffs"] += 1
                    continue
                if not pass_prec:
                    dropped["precursor q (protein kept)"] += 1
                if not pass_prot:
                    dropped["protein q (precursor kept)"] += 1

                run = row[idx["run"]].strip()
                if run not in run_meta:
                    run_meta[run] = (row[idx["condition"]].strip() if "condition" in idx else "",
                                     row[idx["replicate"]].strip() if "replicate" in idx else "",
                                     os.path.basename(path))
                rows_kept += 1

                if args.dry_run:
                    continue

                mod = (row[idx["mod_seq"]] if "mod_seq" in idx
                       else row[idx["precursor_id"]]).strip()
                strip = (row[idx["strip_seq"]].strip() if "strip_seq" in idx
                         else re.sub(r"\[[^\]]*\]|\([^)]*\)|[^A-Z]", "", mod.upper()))
                charge = row[idx["charge"]].strip() if "charge" in idx else ""
                pg = row[idx["protein_group"]].strip()

                k = (lambda s: s) if args.exact else key64
                if pass_prec:
                    seen[(run, "precursors")].add(k(mod + "/" + charge))
                    if strip:
                        seen[(run, "peptides")].add(k(strip))
                if pass_prot:
                    for g in isoform_groups(pg):
                        seen[(run, "protein_isoform_groups")].add(k(g))
                    for a in canonical_accessions(pg):
                        seen[(run, "proteins_canonical")].add(k(a))
            per_file[path] = file_rows
        finally:
            if fh is not None:
                fh.close()
        if args.max_rows and rows_read >= args.max_rows:
            break

    if not run_meta:
        sys.exit("no rows survived filtering -- check --precursor-q/--protein-q and the "
                 "resolved columns above")

    # protease from run name/condition, else from the report filename
    table, prot_src = {}, {}
    for run, (cond, rep, src) in sorted(run_meta.items()):
        sample, protease = derive(run, cond, rep, args.protease_regex, args.sample_regex)
        prot_src[run] = "run/condition"
        if protease is None and not args.protease_regex:
            protease = detect_protease(src)
            prot_src[run] = "report filename"
        table[run] = (sample, protease)

    print(f"resolved columns ({len(used_cols)}):")
    for f in ["run", "condition", "replicate", "mod_seq", "strip_seq", "charge",
              "protein_group", "eg_qvalue", "pg_qvalue", "decoy", "identified"]:
        print(f"  {f:15s} {used_cols.get(f, '-- not present --')}")
    print(f"\nrows read {rows_read:,}  kept {rows_kept:,}")
    for reason, n in sorted(dropped.items(), key=lambda x: -x[1]):
        print(f"  dropped {reason:16s} {n:,}")

    print(f"\nderived run table ({len(table)} runs) -- CHECK THIS:")
    print(f"  {'run':<40s} {'sample':<10s} {'protease':<10s} from")
    for run, (sample, protease) in sorted(table.items(), key=lambda x: (x[1][1] or "", x[1][0] or "")):
        print(f"  {run[:40]:<40s} {str(sample)[:10]:<10s} "
              f"{protease or '?? UNRESOLVED':<10s} {prot_src[run]}")

    bad = {r: v for r, v in table.items()
           if not v[0] or not v[1] or str(v[1]).startswith("AMBIGUOUS")}
    if bad:
        print(f"\n{len(bad)} run(s) could not be resolved to one sample + one protease.")
        print("Pass --protease-regex / --sample-regex; nothing is guessed for you.")
        if not args.dry_run:
            sys.exit(1)

    proteases = sorted({p for _, p in table.values()})
    samples = sorted({s for s, _ in table.values()})
    print(f"\n{len(samples)} sample(s) x {len(proteases)} protease(s): {', '.join(proteases)}")
    if args.dry_run:
        print("\n--dry-run: stopping before counting.")
        return 0

    # run quality gate
    gate_unit = "precursors"
    failed = {run: len(seen.get((run, gate_unit), set()))
              for run in table
              if len(seen.get((run, gate_unit), set())) < args.min_run_count}
    if failed:
        print(f"\nexcluded {len(failed)} run(s) with < {args.min_run_count:,} "
              f"{gate_unit} (--min-run-count):")
        for run, n in sorted(failed.items(), key=lambda x: x[1]):
            s, p = table[run]
            print(f"  {run[:44]:<44s} {s:<8s} {p:<8s} {n:,} {gate_unit}")
    usable = {run: v for run, v in table.items() if run not in failed}

    digests_per_sample = defaultdict(set)
    for run, (sample, protease) in usable.items():
        digests_per_sample[sample].add(protease)
    thin = {s: d for s, d in digests_per_sample.items()
            if len(d) < args.min_union_digests}
    if thin:
        print(f"\nno union for {len(thin)} sample(s) with < "
              f"{args.min_union_digests} usable digests (--min-union-digests):")
        for s, d in sorted(thin.items()):
            print(f"  {s:<8s} has only {', '.join(sorted(d)) or 'none'}")

    # per-run counts and per-sample unions
    os.makedirs(args.out, exist_ok=True)
    tidy = []                                   # (sample, protease, unit, count)
    union_bound = {}
    for unit in UNITS:
        by_sample = defaultdict(list)           # sample -> [(protease, set)]
        for run, (sample, protease) in usable.items():
            s = seen.get((run, unit), set())
            tidy.append((sample, protease, unit, len(s)))
            by_sample[sample].append((protease, s))
        # union across proteases, deduplicated on the same identity key
        parts_sum = union_n = 0
        for sample, parts in by_sample.items():
            if sample in thin:
                continue
            u = set().union(*[s for _, s in parts]) if parts else set()
            tidy.append((sample, "All", unit, len(u)))
            parts_sum += sum(len(s) for _, s in parts)
            union_n += len(u)
        # union FDR upper bound: q * (sum of per-search sizes) / union size
        q_unit = args.protein_q if unit.startswith("protein") else args.precursor_q
        if union_n:
            union_bound[unit] = {"bound": round(q_unit * parts_sum / union_n, 5),
                                 "per_search_q": q_unit}

    counts_csv = os.path.join(args.out, "counts.csv")
    with open(counts_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "protease", "unit", "count"])
        for r in sorted(tidy, key=lambda r: (r[2], r[1], r[0])):
            w.writerow(r)

    print(f"\nwrote {counts_csv}  ({len(tidy)} rows)")
    print(f"\n{'unit':<24s} {'protease':<14s} {'n':>3s} {'median':>10s} {'min':>10s} {'max':>10s}")
    for unit in UNITS:
        for protease in proteases + ["All"]:
            vals = [c for s, p, u, c in tidy if u == unit and p == protease]
            if vals:
                print(f"{unit:<24s} {protease:<14s} {len(vals):>3d} "
                      f"{statistics.median(vals):>10,.0f} {min(vals):>10,} {max(vals):>10,}")

    print(f"\nunion FDR upper bound across {len(proteases)} independent searches:")
    for unit, b in union_bound.items():
        print(f"  {unit:<24s} <= {b['bound']:.2%}  "
              f"(vs {b['per_search_q']:.2%} within one search)")
    print("  Not corrected for -- reported so the 'All' bars can be qualified.")

    prov = {
        "reports": [{"path": p, "rows": n} for p, n in per_file.items()],
        "resolved_columns": used_cols,
        "precursor_q": args.precursor_q,
        "protein_q": args.protein_q,
        "identity": "exact strings" if args.exact else "blake2b-64 digests",
        "rows_read": rows_read,
        "rows_kept": rows_kept,
        "dropped": dict(dropped),
        "runs": {r: {"sample": v[0], "protease": v[1]} for r, v in table.items()},
        "samples": samples,
        "proteases": proteases,
        "union_fdr_upper_bound": union_bound,
        "units": UNITS,
    }
    prov_path = os.path.join(args.out, "provenance.json")
    with open(prov_path, "w") as fh:
        json.dump(prov, fh, indent=2, sort_keys=True)
    print(f"wrote {prov_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
