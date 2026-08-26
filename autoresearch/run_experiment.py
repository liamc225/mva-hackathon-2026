#!/usr/bin/env python3
"""Fixed proxy checks for safe iteration on the public ranking code.

These cases are synthetic engineering fixtures, not the hackathon patient and
not evidence about the confirmed causal variant. They exist to catch obvious
regressions in consequence, rarity, depth, and mosaic handling.
"""

from __future__ import annotations

import json
from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mva_ranker.ranker import Candidate, _score  # noqa: E402


CASES = [
    Candidate("chr1", 1, "A", "G", gene="BUB1B", consequence="stop_gained", impact="HIGH", af=1e-6, dp=80, vaf=0.23, gq=60),
    Candidate("chr2", 2, "C", "T", gene="CEP57", consequence="missense_variant", impact="MODERATE", af=1e-5, dp=60, vaf=0.40, gq=55),
    Candidate("chr3", 3, "G", "A", gene="TRIP13", consequence="splice_region_variant", impact="LOW", af=1e-5, dp=45, vaf=0.14, gq=45),
    Candidate("chr4", 4, "T", "C", gene="OTHER", consequence="synonymous_variant", impact="LOW", af=0.2, dp=12, vaf=0.49, gq=20),
]


def main() -> None:
    scored = [_score(case) for case in CASES]
    ranked = sorted(scored, key=lambda item: item.score, reverse=True)
    order = [item.gene for item in ranked]
    passed = order == ["BUB1B", "CEP57", "TRIP13", "OTHER"]
    result = {
        "proxy_only": True,
        "passed": passed,
        "ranking": order,
        "scores": {item.gene: round(item.score, 4) for item in scored},
    }
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
