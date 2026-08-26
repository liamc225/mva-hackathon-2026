#!/usr/bin/env python3
"""Generate a Track 1-compatible ranked CSV from a local VCF."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mva_ranker.ranker import rank_vcf, write_submission  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--max-rows", type=int, default=10)
    args = parser.parse_args()

    proband_id, candidates = rank_vcf(args.vcf, max_rows=args.max_rows)
    write_submission(proband_id, candidates, args.out)
    print(f"proband_id={proband_id}")
    for index, candidate in enumerate(candidates, 1):
        print(
            f"{index:02d} {candidate.chrom}:{candidate.pos} "
            f"{candidate.ref}>{candidate.alt} gene={candidate.gene or '?'} "
            f"score={candidate.score:.2f} vaf={candidate.vaf}"
        )
    print(f"wrote {len(candidates)} rows to {args.out}")


if __name__ == "__main__":
    main()

