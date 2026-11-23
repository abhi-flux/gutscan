"""Command-line use:  python -m gutscan.cli sample.fasta [--csv out.csv]"""
import argparse, sys, time
from .cleaning import clean_file, CleanParams
from .refdb import RefDB
from .classify import run_classification, Params
from .report import to_csv, to_text


def main(argv=None):
    ap = argparse.ArgumentParser(prog="gutscan")
    ap.add_argument("fasta")
    ap.add_argument("--csv")
    ap.add_argument("--min-len", type=int, default=150)
    ap.add_argument("--min-reads", type=int, default=3)
    ap.add_argument("--extra-refs")
    a = ap.parse_args(argv)
    t = time.time()
    uniq, st = clean_file(a.fasta, CleanParams(min_len=a.min_len))
    db = RefDB.load(extra=a.extra_refs)
    res = run_classification(db, uniq, Params(min_reads=a.min_reads))
    print(to_text(res, st))
    if a.csv:
        to_csv(res, st, a.csv)
    print(f"[{time.time() - t:.1f}s]")


if __name__ == "__main__":
    main()
