#!/usr/bin/env python3
"""Rank called variants in the public MVA gene panel."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mva_ranker.panel import MVA_GENE_INTERVALS, rank_panel  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--max-rows", type=int, default=100)
    parser.add_argument("--flank", type=int, default=50_000)
    parser.add_argument("--genes", nargs="*", choices=sorted(MVA_GENE_INTERVALS))
    args = parser.parse_args()
    proband_id, candidates = rank_panel(
        args.vcf,
        max_rows=args.max_rows,
        flank=args.flank,
        genes=args.genes or None,
    )
    print(f"proband_id={proband_id}")
    for index, candidate in enumerate(candidates, 1):
        print(
            f"{index:03d} {candidate.chrom}:{candidate.pos} "
            f"{candidate.ref}>{candidate.alt} gene={candidate.gene} "
            f"score={candidate.score:.2f} dp={candidate.dp} vaf={candidate.vaf}"
        )


if __name__ == "__main__":
    main()

