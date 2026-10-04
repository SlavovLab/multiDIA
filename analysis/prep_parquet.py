#!/usr/bin/env python3
"""Convert Spectronaut TSV reports to parquet and verify them; --delete drops the TSV.

    python3 prep_parquet.py 'data/search/*.tsv' --reference 'data/search/GluC-*.parquet' --delete
"""

import argparse
import glob
import os
import sys


NULLS = ["NaN", "", "N/A", "NA", "#N/A", "Filtered", "NULL"]


def convert(path, out, reference=None, compression="zstd"):
    """-> (rows, columns). Streams the TSV into one parquet file."""
    import pyarrow as pa
    import pyarrow.csv as csv
    import pyarrow.parquet as pq

    types = {}
    if reference:
        ref = pq.ParquetFile(reference).schema_arrow
        types = {n: ref.field(n).type for n in ref.names}

    ro = csv.ReadOptions(block_size=64 * 1024 * 1024)
    po = csv.ParseOptions(delimiter="\t")
    # pyarrow raises on a type given for an absent column.
    with open(path, newline="") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    co = csv.ConvertOptions(null_values=NULLS, strings_can_be_null=True,
                            column_types={n: t for n, t in types.items()
                                          if n in header})
    reader = csv.open_csv(path, read_options=ro, parse_options=po,
                          convert_options=co)
    writer, rows = None, 0
    try:
        for batch in reader:
            if writer is None:
                writer = pq.ParquetWriter(out, batch.schema,
                                          compression=compression)
            writer.write_table(pa.Table.from_batches([batch], batch.schema))
            rows += batch.num_rows
    finally:
        if writer is not None:
            writer.close()
        reader.close()
    if writer is None:
        sys.exit(f"{path}: no data rows")
    return rows, len(header)


def verify(path, out, rows, cols):
    """Re-open the parquet and check it against the TSV. -> [complaint]."""
    import pyarrow.parquet as pq
    bad = []
    pf = pq.ParquetFile(out)
    if pf.metadata.num_rows != rows:
        bad.append(f"parquet has {pf.metadata.num_rows:,} rows, wrote {rows:,}")
    if pf.metadata.num_columns != cols:
        bad.append(f"parquet has {pf.metadata.num_columns} columns, "
                   f"TSV header had {cols}")
    with open(path, "rb") as fh:
        lines = sum(chunk.count(b"\n") for chunk in
                    iter(lambda: fh.read(1 << 24), b""))
    if lines - 1 != rows:
        bad.append(f"TSV has {lines - 1:,} data lines, parquet {rows:,} "
                   f"(embedded newlines or a quoting difference)")
    need = {"R.FileName": "string", "PG.ProteinGroups": "string",
            "PEP.StrippedSequence": "string", "EG.Qvalue": "double",
            "FG.Quantity": "double", "EG.UsedForPeptideQuantity": "bool"}
    sch = pf.schema_arrow
    for n, want in need.items():
        if n not in sch.names:
            bad.append(f"missing column {n}")
        elif str(sch.field(n).type) != want:
            bad.append(f"{n} is {sch.field(n).type}, expected {want} "
                       f"(is `NaN` declared as a null token?)")
    return bad


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", nargs="+")
    ap.add_argument("--reference", help="parquet whose dtypes to match")
    ap.add_argument("--compression", default="zstd")
    ap.add_argument("--delete", action="store_true",
                    help="remove the TSV if verification passes")
    args = ap.parse_args(argv)

    paths = sorted(p for pat in args.tsv for p in glob.glob(pat))
    if not paths:
        sys.exit("no TSV matched")
    ref = None
    if args.reference:
        r = sorted(glob.glob(args.reference))
        if not r:
            sys.exit(f"reference not found: {args.reference}")
        ref = r[0]
        print(f"matching dtypes to {os.path.basename(ref)}")

    rc = 0
    for p in paths:
        out = os.path.splitext(p)[0] + ".parquet"
        src = os.path.getsize(p)
        print(f"\n{os.path.basename(p)}  ({src / 1e6:.0f} MB)")
        rows, cols = convert(p, out, ref, args.compression)
        dst = os.path.getsize(out)
        print(f"  -> {os.path.basename(out)}  {rows:,} rows x {cols} cols, "
              f"{dst / 1e6:.0f} MB ({src / dst:.1f}x smaller)")
        bad = verify(p, out, rows, cols)
        if ref:
            import pyarrow.parquet as pq
            a = set(pq.ParquetFile(out).schema_arrow.names)
            b = set(pq.ParquetFile(ref).schema_arrow.names)
            if a - b:
                print(f"  columns not in the reference: "
                      f"{', '.join(sorted(a - b))}")
            if b - a:
                print(f"  columns the reference has and this does not: "
                      f"{', '.join(sorted(b - a))}")
        if bad:
            rc = 1
            for m in bad:
                print(f"  FAIL  {m}")
            print("  keeping the TSV")
            continue
        print("  ok    rows, columns and the required dtypes all check out")
        if args.delete:
            os.remove(p)
            print(f"  removed {os.path.basename(p)}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
