#!/usr/bin/env python3
"""Write a synthetic Spectronaut-shaped report to exercise prep_counts.py.

    python3 make_synthetic.py --out /tmp/Report.tsv --samples 12 --scale 300
"""

import argparse
import random

HEADER = ["R.FileName", "R.Condition", "R.Replicate", "PG.ProteinGroups",
          "PG.ProteinAccessions", "PG.Qvalue", "PEP.StrippedSequence",
          "EG.ModifiedSequence", "EG.PrecursorId", "EG.Qvalue", "EG.IsDecoy",
          "FG.Charge", "FG.Quantity"]

AA = "ACDEFGHIKLMNPQRSTVWY"
PROTEASES = ["Trypsin", "LysC", "GluC"]
# Fraction of each peptide pool a digest sees.
YIELD = {"Trypsin": 0.62, "LysC": 0.52, "GluC": 0.34}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/Report.tsv")
    ap.add_argument("--samples", type=int, default=12)
    ap.add_argument("--scale", type=int, default=300,
                    help="size of the shared peptide pool")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--fail-run", default="S03/LysC",
                    help="'sample/protease' to emit at 8%% depth, or '' for none")
    args = ap.parse_args()
    rng = random.Random(args.seed)

    def peptide():
        return "".join(rng.choice(AA) for _ in range(rng.randint(7, 22)))

    shared = [peptide() for _ in range(args.scale)]
    private = {p: [peptide() for _ in range(args.scale // 2)] for p in PROTEASES}

    n_prot = max(args.scale // 6, 20)
    proteins = []
    for i in range(n_prot):
        acc = f"P{10000 + i * 7}"
        proteins.append(acc)
        if i % 5 == 0:
            proteins.append(f"{acc}-{rng.randint(2, 4)}")

    samples = [f"S{i:02d}" for i in range(1, args.samples + 1)]
    fail = tuple(args.fail_run.split("/")) if args.fail_run else ("", "")

    rows = 0
    with open(args.out, "w", newline="") as fh:
        fh.write("\t".join(HEADER) + "\n")
        for s in samples:
            for prot in PROTEASES:
                run = f"20260828_multiDIA_{s}_{prot}_dia"
                depth = YIELD[prot] * rng.uniform(0.94, 1.06)
                if (s, prot) == fail:
                    depth *= 0.08
                pool = shared + private[prot]
                for pep in pool:
                    if rng.random() > depth:
                        continue
                    pg = rng.choice(proteins)
                    if rng.random() < 0.08:
                        pg = pg + ";" + rng.choice(proteins)
                    for z in ([2, 3] if rng.random() < 0.35 else [2]):
                        decoy = rng.random() < 0.01
                        egq = round(rng.uniform(0.011, 0.4), 4) if rng.random() < 0.02 \
                            else round(rng.uniform(0, 0.0099), 5)
                        pgq = round(rng.uniform(0.011, 0.2), 4) if rng.random() < 0.01 \
                            else round(rng.uniform(0, 0.0099), 5)
                        mod = pep if rng.random() < 0.8 else \
                            pep[:2] + "[Carbamidomethyl (C)]" + pep[2:]
                        fh.write("\t".join([
                            run, prot, s, pg, pg.split(";")[0], str(pgq), pep,
                            f"_{mod}_", f"_{mod}_.{z}", str(egq),
                            "True" if decoy else "False", str(z),
                            f"{rng.uniform(1e3, 1e6):.1f}",
                        ]) + "\n")
                        rows += 1

    print(f"wrote {args.out}: {rows:,} rows, {len(samples)} samples x "
          f"{len(PROTEASES)} proteases"
          + (f", failed run {'/'.join(fail)}" if args.fail_run else ""))


if __name__ == "__main__":
    main()
