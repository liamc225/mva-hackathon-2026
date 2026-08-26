#!/usr/bin/env python3
"""Print safe, non-content metadata for a local challenge VCF."""

from __future__ import annotations

import argparse
from collections import Counter

from cyvcf2 import VCF


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    args = parser.parse_args()

    vcf = VCF(args.vcf)
    print(f"samples: {list(vcf.samples)}")
    print(f"contigs: {len(vcf.seqnames)}")
    print(f"has ANN: {'ID=ANN,' in (vcf.raw_header or '')}")
    print(f"has CSQ: {'ID=CSQ,' in (vcf.raw_header or '')}")
    counts = Counter()
    total = 0
    for variant in vcf:
        total += 1
        counts[len(str(variant.REF)) == len(str(variant.ALT[0]))] += 1
    print(f"records: {total}")
    print(f"SNV-like records: {counts[True]}")
    print(f"indel-like records: {counts[False]}")


if __name__ == "__main__":
    main()

