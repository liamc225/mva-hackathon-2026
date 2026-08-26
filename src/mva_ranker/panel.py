"""MVA-focused region selection and pair-aware ranking helpers.

The intervals are public GRCh38 gene coordinates, not patient-derived data.
They are deliberately a soft panel: calls outside these genes remain visible
when the genome-wide ranker is used.  The panel is useful for a no-annotation
VCF because it gives the downstream local GENCODE annotator a bounded set of
records to inspect without uploading the VCF to a third-party service.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

from cyvcf2 import VCF

from .ranker import Candidate, _candidate_from_variant, _score


# GRCh38, 1-based inclusive gene bounds from Ensembl public coordinates.  The
# values are intentionally kept in source so reviewers can audit the panel.
MVA_GENE_INTERVALS: dict[str, tuple[str, int, int]] = {
    "BUB1B": ("15", 40160984, 40221137),
    "CEP57": ("11", 95789965, 95837070),
    "TRIP13": ("5", 892849, 919357),
    "CENATAC": ("11", 118998051, 119015811),
    "MAD1L1": ("7", 1815787, 2233243),
    "MAD2L1BP": ("6", 43629494, 43640960),
    "CEP192": ("18", 12991283, 13125053),
    "SLF2": ("10", 100912963, 100965134),
    "SMC5": ("9", 70258270, 70354874),
    "BUB1": ("2", 110635468, 110678098),
    "BUB3": ("10", 123154365, 123313144),
}


def _has_alt_genotype(variant: Any) -> bool:
    try:
        genotype = variant.genotypes[0]
        return any(int(allele) > 0 for allele in genotype[:2] if allele is not None and allele >= 0)
    except (AttributeError, IndexError, TypeError, ValueError):
        return False


def iter_panel_candidates(
    vcf_path: str | Path,
    *,
    flank: int = 50_000,
    genes: Iterable[str] | None = None,
) -> Iterable[Candidate]:
    """Yield called variants in the public MVA panel, deduplicated by allele."""

    vcf = VCF(str(vcf_path))
    ann_fields = []
    csq_fields = []
    requested = set(genes or MVA_GENE_INTERVALS)
    seen: set[tuple[str, int, str, str]] = set()
    for gene in requested:
        if gene not in MVA_GENE_INTERVALS:
            raise ValueError(f"Unknown panel gene: {gene}")
        chrom, start, end = MVA_GENE_INTERVALS[gene]
        for variant in vcf(f"{chrom}:{max(1, start - flank)}-{end + flank}"):
            alt = str(variant.ALT[0]) if variant.ALT else ""
            if not alt or alt.startswith("<") or alt == "*" or not _has_alt_genotype(variant):
                continue
            key = (str(variant.CHROM), int(variant.POS), str(variant.REF), alt)
            if key in seen:
                continue
            seen.add(key)
            candidate = _candidate_from_variant(variant, ann_fields, csq_fields)
            candidate.gene = gene
            candidate = _score(candidate)
            candidate.evidence.append(f"within {gene} +/- {flank:,} bp public panel interval")
            yield candidate


def rank_panel(
    vcf_path: str | Path,
    *,
    max_rows: int = 100,
    flank: int = 50_000,
    genes: Iterable[str] | None = None,
) -> tuple[str, list[Candidate]]:
    """Return a quality/consequence-ranked list from the panel regions."""

    vcf = VCF(str(vcf_path))
    samples = list(vcf.samples)
    proband_id = samples[0] if samples else "proband"
    candidates = list(iter_panel_candidates(vcf_path, flank=flank, genes=genes))
    candidates.sort(key=lambda item: (item.score, item.gq or 0, item.dp or 0), reverse=True)
    return proband_id, candidates[:max_rows]


def write_pair_submission(
    proband_id: str,
    pairs: list[tuple[Candidate, Candidate | None, str, float]],
    output: str | Path,
) -> None:
    """Write the Track 1 schema, including optional two-allele hypotheses.

    ``epcr`` is a monotone prioritisation score, not a calibrated probability;
    callers should state that limitation in the methods report.
    """

    fields = [
        "proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2",
        "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes",
    ]
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for first, second, finding_type, epcr in pairs:
            row = first.as_row(proband_id)
            row.update({"epcr": f"{epcr:.3f}", "finding_type": finding_type})
            if second is not None:
                row.update({
                    "chrom_2": second.chrom,
                    "pos_2": second.pos,
                    "ref_2": second.ref,
                    "alt_2": second.alt,
                    "notes": (
                        f"candidate pair; allele 1: {first.gene} {first.consequence or 'unannotated'}; "
                        f"allele 2: {second.gene} {second.consequence or 'unannotated'}. "
                        + "; ".join(first.evidence + second.evidence)
                    ),
                })
            writer.writerow(row)

