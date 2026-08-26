"""A conservative first-pass ranker for a single rare-disease VCF.

The ranker is deliberately transparent. It is a prioritisation aid, not a
clinical diagnostic model. It works best when the VCF has consequence and
population-frequency annotations (ANN or CSQ); it still emits quality-aware
rows when those annotations are absent. ``INFO/AF`` is intentionally *not*
treated as population frequency: in many single-sample VCFs it is the sample
allele frequency and is redundant with ``FORMAT/AD``.
"""

from __future__ import annotations

import csv
import heapq
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable

if TYPE_CHECKING:
    from cyvcf2 import VCF


# Established MVA biology is used as a soft prior only. Keeping this list in
# code makes the prior explicit and reviewable.
MVA_GENE_PRIORS: dict[str, float] = {
    "BUB1B": 5.0,
    "CEP57": 4.5,
    "TRIP13": 4.0,
    # Additional MVA/related chromosomal-instability genes reported in recent
    # literature. These are intentionally lower soft priors than the original
    # three genes and never act as a hard filter.
    "CENATAC": 3.8,
    "MAD1L1": 3.8,
    "MAD2L1BP": 3.5,
    "CEP192": 3.5,
    "SLF2": 3.5,
    "SMC5": 3.5,
    "BUB1": 2.5,
    "BUB3": 2.5,
}

LOF_TERMS = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "start_lost",
    "stop_lost",
}
SPLICE_TERMS = {"splice_region_variant"}
MISSENSE_TERMS = {"missense_variant"}
SYNONYMOUS_TERMS = {"synonymous_variant"}


@dataclass
class Annotation:
    gene: str = ""
    consequence: str = ""
    impact: str = ""
    transcript: str = ""
    protein_change: str = ""
    raw: str = ""


@dataclass
class Candidate:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str = ""
    consequence: str = ""
    impact: str = ""
    transcript: str = ""
    protein_change: str = ""
    af: float | None = None
    dp: int | None = None
    vaf: float | None = None
    gq: float | None = None
    quality: float | None = None
    mq: float | None = None
    qd: float | None = None
    fs: float | None = None
    sor: float | None = None
    score: float = 0.0
    evidence: list[str] = field(default_factory=list)

    def as_row(self, proband_id: str) -> dict[str, Any]:
        note = "; ".join(self.evidence)
        return {
            "proband_id": proband_id,
            "chrom_1": self.chrom,
            "pos_1": self.pos,
            "ref_1": self.ref,
            "alt_1": self.alt,
            "chrom_2": "",
            "pos_2": "",
            "ref_2": "",
            "alt_2": "",
            "epcr": "",
            "finding_type": "primary",
            "notes": note,
        }


def _header_format(vcf: VCF, key: str) -> list[str]:
    """Extract pipe-delimited annotation field names from a VCF header."""

    raw = vcf.raw_header or ""
    for line in raw.splitlines():
        if f'ID={key},' not in line:
            continue
        match = re.search(r'(?:Format|Description)="?[^"\n]*?Format: ([^"\n]+)', line)
        if match:
            return [x.strip() for x in match.group(1).split("|")]
    return []


def _parse_annotation(value: Any, fields: list[str]) -> Annotation:
    if value is None:
        return Annotation()
    raw = str(value[0] if isinstance(value, (list, tuple)) else value)
    parts = raw.split("|")
    mapping = {fields[i]: parts[i] for i in range(min(len(fields), len(parts)))}
    consequence = mapping.get("Consequence", mapping.get("Annotation", ""))
    gene = mapping.get("SYMBOL", mapping.get("Gene_Name", mapping.get("Gene", "")))
    impact = mapping.get("IMPACT", mapping.get("Impact", ""))
    transcript = mapping.get("Feature", mapping.get("Feature_ID", ""))
    protein = mapping.get("HGVSp", mapping.get("HGVS.p", ""))
    return Annotation(
        gene=gene,
        consequence=consequence,
        impact=impact,
        transcript=transcript,
        protein_change=protein,
        raw=raw,
    )


def _as_float(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        value = value[0] if value else None
    if value in (None, "", "."):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _info_float(variant: Any, keys: Iterable[str]) -> float | None:
    for key in keys:
        try:
            value = variant.INFO.get(key)
        except Exception:
            value = None
        result = _as_float(value)
        if result is not None:
            return result
    return None


def _sample_number(variant: Any, key: str) -> float | None:
    try:
        values = variant.format(key)
    except Exception:
        return None
    if values is None or len(values) == 0:
        return None
    value = values[0]
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, (list, tuple)):
        value = value[0] if value else None
    return _as_float(value)


def _sample_vaf(variant: Any, dp: float | None) -> float | None:
    try:
        ad = variant.format("AD")
    except Exception:
        ad = None
    if ad is not None and len(ad):
        values = ad[0]
        if hasattr(values, "tolist"):
            values = values.tolist()
        if isinstance(values, (list, tuple)) and len(values) >= 2:
            ref_count = _as_float(values[0])
            alt_count = sum(x for x in (_as_float(v) for v in values[1:]) if x is not None)
            denom = (ref_count or 0.0) + alt_count
            if denom > 0:
                return alt_count / denom
    ao = _sample_number(variant, "AO")
    if ao is not None and dp and dp > 0:
        return ao / dp
    return None


def _consequence_score(consequence: str, impact: str) -> tuple[float, str]:
    terms = {term.strip().lower() for term in re.split(r"[,&]", consequence) if term.strip()}
    if terms & {x.lower() for x in LOF_TERMS} or impact.upper() == "HIGH":
        return 5.0, "high-impact/LoF consequence"
    if terms & {x.lower() for x in MISSENSE_TERMS} or impact.upper() == "MODERATE":
        return 2.5, "missense/moderate-impact consequence"
    if terms & {x.lower() for x in SPLICE_TERMS}:
        return 2.0, "splice-region consequence"
    if terms & {x.lower() for x in SYNONYMOUS_TERMS}:
        return -1.0, "synonymous consequence"
    return 0.0, "consequence unannotated or low-information"


def _score(candidate: Candidate) -> Candidate:
    score = 0.0
    evidence: list[str] = []

    gene_prior = MVA_GENE_PRIORS.get(candidate.gene.upper(), 0.0)
    if gene_prior:
        score += gene_prior
        evidence.append(f"MVA gene prior {candidate.gene.upper()} (+{gene_prior:g})")

    consequence_score, consequence_note = _consequence_score(candidate.consequence, candidate.impact)
    score += consequence_score
    if consequence_score:
        evidence.append(f"{consequence_note} ({consequence_score:+g})")

    if candidate.af is not None:
        if candidate.af <= 1e-5:
            score += 2.5
            evidence.append("very rare population frequency (+2.5)")
        elif candidate.af <= 1e-4:
            score += 2.0
            evidence.append("rare population frequency (+2)")
        elif candidate.af <= 1e-2:
            score += 0.5
            evidence.append("uncommon population frequency (+0.5)")
        else:
            score -= 3.0
            evidence.append("common population frequency (-3)")

    if candidate.vaf is not None:
        # A mosaic call can be below the usual heterozygous 0.5 expectation.
        # This is a soft bonus only; low depth remains visible in the evidence.
        if 0.03 <= candidate.vaf <= 0.45:
            score += 1.25
            evidence.append(f"mosaic-compatible VAF {candidate.vaf:.3f} (+1.25)")
        elif 0.005 <= candidate.vaf < 0.03:
            score += 0.5
            evidence.append(f"low-level alternate VAF {candidate.vaf:.3f} (+0.5)")

    if candidate.dp is not None:
        if candidate.dp >= 30:
            score += 1.0
            evidence.append(f"adequate depth DP={candidate.dp:g} (+1)")
        elif candidate.dp < 10:
            score -= 1.0
            evidence.append(f"low depth DP={candidate.dp:g} (-1)")

    if candidate.gq is not None and candidate.gq >= 30:
        score += 0.5
        evidence.append(f"good genotype quality GQ={candidate.gq:g} (+0.5)")

    candidate.score = score
    candidate.evidence = evidence
    return candidate


def _candidate_from_variant(variant: Any, ann_fields: list[str], csq_fields: list[str]) -> Candidate:
    ann = _parse_annotation(variant.INFO.get("ANN"), ann_fields)
    if not ann.gene:
        ann = _parse_annotation(variant.INFO.get("CSQ"), csq_fields)
    dp = _sample_number(variant, "DP")
    gq = _sample_number(variant, "GQ")
    return Candidate(
        chrom=str(variant.CHROM),
        pos=int(variant.POS),
        ref=str(variant.REF),
        alt=str(variant.ALT[0]),
        gene=ann.gene,
        consequence=ann.consequence,
        impact=ann.impact,
        transcript=ann.transcript,
        protein_change=ann.protein_change,
        # Never use INFO/AF as population frequency. In the challenge VCF it
        # is a single-sample allele fraction, not a gnomAD-style frequency.
        af=_info_float(variant, ("gnomAD_AF", "GNOMAD_AF")),
        dp=int(dp) if dp is not None else None,
        vaf=_sample_vaf(variant, dp),
        gq=gq,
        quality=_as_float(variant.QUAL),
        mq=_info_float(variant, ("MQ",)),
        qd=_info_float(variant, ("QD",)),
        fs=_info_float(variant, ("FS",)),
        sor=_info_float(variant, ("SOR",)),
    )


def rank_vcf(vcf_path: str | Path, max_rows: int = 10) -> tuple[str, list[Candidate]]:
    """Return the proband ID and top-ranked candidates from a VCF."""

    from cyvcf2 import VCF

    vcf = VCF(str(vcf_path))
    samples = list(vcf.samples)
    proband_id = samples[0] if samples else "proband"
    ann_fields = _header_format(vcf, "ANN")
    csq_fields = _header_format(vcf, "CSQ")

    # Keep only the best rows while streaming through the VCF. A whole-genome
    # VCF can contain millions of records; retaining every Candidate would
    # waste the Colab runtime's memory before ranking is complete.
    heap: list[tuple[tuple[float, int, int], int, Candidate]] = []
    counter = 0
    for variant in vcf:
        # Skip symbolic alleles and reference blocks for a variant-prediction
        # submission; retain SNVs/indels for the initial ranking.
        alt = str(variant.ALT[0]) if variant.ALT else ""
        if alt.startswith("<") or alt == "*":
            continue
        candidate = _score(_candidate_from_variant(variant, ann_fields, csq_fields))
        key = (candidate.score, int(bool(candidate.gene)), candidate.dp or 0)
        item = (key, counter, candidate)
        counter += 1
        if len(heap) < max_rows:
            heapq.heappush(heap, item)
        elif key > heap[0][0]:
            heapq.heapreplace(heap, item)

    ranked = [item[2] for item in heap]
    ranked.sort(key=lambda c: (c.score, c.gene != "", c.dp or 0), reverse=True)
    return proband_id, ranked[:max_rows]


def write_submission(proband_id: str, candidates: list[Candidate], output: str | Path) -> None:
    fields = [
        "proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2",
        "pos_2", "ref_2", "alt_2", "epcr", "finding_type", "notes",
    ]
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for index, candidate in enumerate(candidates):
            row = candidate.as_row(proband_id)
            # EPCR is a ranking signal, not a calibrated probability. It is
            # monotone and deliberately leaves room for near-ties.
            row["epcr"] = f"{max(0.001, 1.0 - index * 0.08):.3f}"
            writer.writerow(row)
