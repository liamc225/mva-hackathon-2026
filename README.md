# MVA Hackathon 2026 — reproducible variant and mechanism analysis

This repository contains the public, data-free analysis code for the SageBio
Rare Disease, Real Kid: MVA Hackathon 2026.

The private challenge files are intentionally not stored here. They must be
downloaded directly into a local/Colab `data/` directory after the participant
has received gated access. `data/`, `runtime/`, and `outputs/` are ignored by
Git. The code emits only the files needed for a submission; it does not upload
genomic data anywhere.

## Workflow

1. Install the small analysis dependencies in Colab:

   ```bash
   pip install -r requirements.txt
   ```

2. Download only the VCF, VCF index, and clinical phenotype document from the
   gated dataset. The FASTQs are optional and should be pulled only if read-
   level validation becomes necessary.

3. Inspect the VCF and phenotype document:

   ```bash
   python scripts/inspect_dataset.py --vcf data/WGS_EX2312012_HGWCNDSX7.vcf.gz
   python scripts/extract_phenotype.py --docx data/Challenge_Clinical_Phenotype_1.docx
   ```

4. Run the deterministic first-pass ranking. The optional annotation table can
   be produced by VEP, SnpEff, or another independently documented annotator.

   ```bash
   python scripts/rank_variants.py \
     --vcf data/WGS_EX2312012_HGWCNDSX7.vcf.gz \
     --out runtime/liam_variant_predictions.csv \
     --max-rows 10
   ```

5. Review the evidence table and validate high-priority calls against the
   original VCF and, when needed, the FASTQ reads. Candidate CSVs and reports
   remain in `runtime/` until the final human review and portal submission.

## Method principles

- Treat Mosaic Variegated Aneuploidy (MVA) gene knowledge as a prior, not as a
  hard filter; keep novel or incidental findings visible.
- Reward rare, high-quality variants with biologically plausible consequences,
  but do not discard low allele-fraction calls solely because they are mosaic.
- Keep the final ranking auditable: every score is accompanied by its evidence
  components and a short rationale.
- Use the autoresearch-style loop only for reproducible engineering and proxy
  checks. Never optimize against the hidden clinical answer or repeatedly burn
  the six live submissions.
- Track 2 claims are hypotheses for follow-up, not treatment recommendations.

## Privacy and retention

The challenge rules require data deletion within 30 days after the hackathon
closes. This repository contains no patient-derived data, raw clinical notes,
or candidate output. Before the deadline, delete Colab/Drive caches and any
derived private files, then follow the organizers' deletion-confirmation
instructions.

## License

MIT for the public analysis code in this repository. Submission artifacts are
subject to the hackathon's stated CC BY 4.0 terms.

