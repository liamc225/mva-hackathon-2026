#!/usr/bin/env python3
"""Annotate panel calls with local GENCODE CDS and public reference sequence.

This intentionally uses a public GENCODE GTF plus Ensembl's public reference
sequence endpoint. The private VCF is read locally; only public coordinates
are requested. The output is a local, ignored table for review and is not a
clinical interpretation.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import requests
from cyvcf2 import VCF

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mva_ranker.panel import MVA_GENE_INTERVALS  # noqa: E402
from mva_ranker.ranker import Candidate, _candidate_from_variant, _score  # noqa: E402


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def parse_attributes(value: str) -> dict[str, str]:
    return dict(re.findall(r'(\S+)\s+"([^"]+)"', value))


@dataclass
class Transcript:
    gene: str
    transcript: str
    chrom: str
    strand: int
    cds: list[tuple[int, int]] = field(default_factory=list)
    canonical: bool = False
    coding_sequence: str = ""

    def prepare(self, genomic_sequence: str, region_start: int) -> None:
        ordered = sorted(self.cds, reverse=self.strand == -1)
        parts = []
        for start, end in ordered:
            piece = genomic_sequence[start - region_start : end - region_start + 1]
            parts.append(piece if self.strand == 1 else reverse_complement(piece))
        self.coding_sequence = "".join(parts)

    def effect(self, pos: int, ref: str, alt: str) -> tuple[str, str] | None:
        ordered = sorted(self.cds, reverse=self.strand == -1)
        offset = 0
        for start, end in ordered:
            if start <= pos <= end:
                index = offset + (pos - start if self.strand == 1 else end - pos)
                if len(ref) != len(alt):
                    consequence = "frameshift_variant" if abs(len(ref) - len(alt)) % 3 else "inframe_indel"
                    return consequence, ""
                if len(ref) != 1 or index >= len(self.coding_sequence):
                    return "coding_sequence_variant", ""
                tx_ref = ref if self.strand == 1 else reverse_complement(ref)
                tx_alt = alt if self.strand == 1 else reverse_complement(alt)
                if self.coding_sequence[index] != tx_ref:
                    return "reference_mismatch", ""
                codon_start = index - (index % 3)
                ref_codon = self.coding_sequence[codon_start : codon_start + 3]
                if len(ref_codon) != 3:
                    return "coding_sequence_variant", ""
                alt_codon = ref_codon[: index % 3] + tx_alt + ref_codon[index % 3 + 1 :]
                aa_ref = CODON_TABLE.get(ref_codon, "X")
                aa_alt = CODON_TABLE.get(alt_codon, "X")
                aa_number = codon_start // 3 + 1
                if aa_ref == aa_alt:
                    consequence = "synonymous_variant"
                elif aa_alt == "*":
                    consequence = "stop_gained"
                elif aa_ref == "*":
                    consequence = "stop_lost"
                else:
                    consequence = "missense_variant"
                protein = f"p.{aa_ref}{aa_number}{aa_alt}"
                return consequence, protein
            offset += end - start + 1
        return None


def load_transcripts(gtf: str | Path, genes: Iterable[str]) -> dict[str, list[Transcript]]:
    wanted = set(genes)
    transcripts: dict[tuple[str, str], Transcript] = {}
    with Path(gtf).open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9 or columns[2] != "CDS":
                continue
            attrs = parse_attributes(columns[8])
            gene = attrs.get("gene_name", attrs.get("gene_id", ""))
            if gene not in wanted:
                continue
            transcript = attrs.get("transcript_id", "")
            if not transcript:
                continue
            key = (gene, transcript)
            item = transcripts.setdefault(
                key,
                Transcript(gene, transcript, columns[0].removeprefix("chr"), 1 if columns[6] == "+" else -1),
            )
            item.cds.append((int(columns[3]), int(columns[4])))
            item.canonical |= "Ensembl_canonical" in attrs.get("tag", "")
    grouped: dict[str, list[Transcript]] = {}
    for item in transcripts.values():
        item.cds = sorted(set(item.cds))
        grouped.setdefault(item.gene, []).append(item)
    return grouped


def fetch_sequence(chrom: str, start: int, end: int, session: requests.Session) -> str:
    url = f"https://rest.ensembl.org/sequence/region/human/{chrom}:{start}..{end}:1"
    response = session.get(url, headers={"Content-Type": "text/plain"}, timeout=60)
    response.raise_for_status()
    return re.sub(r"\s+", "", response.text).upper()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--max-rows", type=int, default=100)
    args = parser.parse_args()

    transcripts = load_transcripts(args.gtf, MVA_GENE_INTERVALS)
    session = requests.Session()
    for gene, models in transcripts.items():
        chrom = models[0].chrom
        start = min(start for model in models for start, _ in model.cds)
        end = max(end for model in models for _, end in model.cds)
        sequence = fetch_sequence(chrom, start, end, session)
        for model in models:
            model.prepare(sequence, start)

    vcf = VCF(args.vcf)
    ann_fields: list[str] = []
    csq_fields: list[str] = []
    rows: list[dict[str, object]] = []
    seen: set[tuple[str, int, str, str]] = set()
    for gene, (chrom, start, end) in MVA_GENE_INTERVALS.items():
        for variant in vcf(f"{chrom}:{start - 50_000}-{end + 50_000}"):
            alt = str(variant.ALT[0]) if variant.ALT else ""
            if not alt or alt.startswith("<") or alt == "*":
                continue
            key = (str(variant.CHROM), int(variant.POS), str(variant.REF), alt)
            if key in seen:
                continue
            seen.add(key)
            effects = [model.effect(int(variant.POS), str(variant.REF), alt) for model in transcripts.get(gene, [])]
            effects = [effect for effect in effects if effect is not None]
            if not effects:
                continue
            consequence, protein = max(
                effects,
                key=lambda value: {
                    "stop_gained": 5,
                    "frameshift_variant": 5,
                    "splice_acceptor_variant": 5,
                    "splice_donor_variant": 5,
                    "missense_variant": 3,
                    "inframe_indel": 2,
                    "synonymous_variant": 1,
                }.get(value[0], 0),
            )
            candidate = _score(_candidate_from_variant(variant, ann_fields, csq_fields))
            candidate.gene = gene
            candidate.consequence = consequence
            candidate.protein_change = protein
            candidate = _score(candidate)
            rows.append({
                "gene": gene,
                "chrom": str(variant.CHROM),
                "pos": int(variant.POS),
                "ref": str(variant.REF),
                "alt": alt,
                "consequence": consequence,
                "protein_change": protein,
                "score": round(candidate.score, 4),
                "dp": candidate.dp,
                "gq": candidate.gq,
                "vaf": round(candidate.vaf, 4) if candidate.vaf is not None else "",
            })

    rows.sort(key=lambda row: (float(row["score"]), row["gq"] or 0, row["dp"] or 0), reverse=True)
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["gene", "chrom", "pos", "ref", "alt", "consequence", "protein_change", "score", "dp", "gq", "vaf"]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows[: args.max_rows])
    print(f"wrote {min(len(rows), args.max_rows)} rows to {output}")


if __name__ == "__main__":
    main()
