# MVA Hackathon 2026 project state

## Completed work

- Inspected the public challenge page, rules, gated dataset metadata, and
  Track 1/Track 2 requirements.
- Registered the solo participant and accepted the gated data-use terms using
  the participant-provided details; Hugging Face access is active.
- Rehydrated a Google Colab runtime and downloaded only the VCF, index, and
  phenotype document. FASTQs remain optional and were not downloaded.
- Built a data-free public repository with transparent ranking, a public
  GRCh38 MVA gene panel, local GENCODE CDS translation, public reference
  sequence retrieval, and a fixed synthetic autoresearch-style regression
  harness.
- Pushed the scaffold and Colab notebook to
  `https://github.com/liamc225/mva-hackathon-2026` on branch
  `codex/mva-hackathon-2026`.
- In private Colab analysis, verified the challenge VCF schema: one sample,
  no ANN/CSQ fields, genotype/depth/allele-depth/quality fields, and about five
  million records. The MVA panel contained 2,853 called records before coding
  consequence filtering.
- Local coding-aware analysis produced a strong provisional BUB1B pair: a
  heterozygous stop-gain and a second heterozygous missense call, both with
  high genotype quality and adequate depth. Phase, population frequency, and
  functional effect remain validation items.
- Submitted Track 1 as `PROBAND01`; the portal reported rank points 100/100,
  F-max 1.000, and a full match at rank 1.
- Published a 2:59 narrated Track 2 pitch with no copyright issues and
  submitted Track 2 with the report, pitch URL, and public repository. The
  portal confirmed the submission was logged for expert review.

## Current behavior

- `scripts/inspect_dataset.py` prints safe VCF metadata.
- `scripts/extract_phenotype.py` extracts the private DOCX into ignored
  `runtime/phenotype.txt`.
- `scripts/rank_mva_panel.py` retrieves called variants from public MVA gene
  intervals and applies transparent quality-aware soft priors.
- `scripts/annotate_coding.py` uses a local GENCODE v47 GTF and public Ensembl
  reference sequence to classify CDS SNVs/indels without uploading the VCF.
- `scripts/rank_variants.py` remains available for a genome-wide annotated VCF.
- `autoresearch/run_experiment.py` runs synthetic checks only; it does not use
  the hidden answer, leaderboard feedback, or patient-derived output.

## Decisions

- Keep all gated challenge files, phenotype text, and derived candidate rows
  out of GitHub. They remain in the authorized Colab runtime until submission
  and deletion.
- Use the BUB1B pair as the lead hypothesis for Track 1, while explicitly
  labeling phase and missense pathogenicity as unresolved.
- Use Track 2 to explain the spindle-assembly-checkpoint mechanism and propose
  confirmatory assays and therapeutic hypotheses, not patient treatment. The
  report prioritizes etoposide and irinotecan as literature-anchored Tier 1
  tests, with sorafenib as a conditional Tier 2 follow-up and spindle-stress
  compounds as controls.
- Do not spend the six live submissions on autonomous tuning. Autoresearch is
  limited to fixed proxy checks, reproducibility, and parser/ranking changes.
- Download FASTQs only if a concrete read-level question cannot be answered
  from VCF genotype and quality evidence.

## Tests

- `python3 -m compileall -q src scripts autoresearch` passes.
- `python3 autoresearch/run_experiment.py` passes with the expected synthetic
  ranking.
- `python3 scripts/rank_mva_panel.py --help` and
  `python3 scripts/annotate_coding.py --help` pass after installing `cyvcf2`.
- Colab direct queries confirmed the two lead BUB1B records are present in the
  downloaded VCF; the nearby public ClinVar splice coordinate is absent from
  this sample.

## Next actions

1. Preserve the submission confirmations and public repository URL for the
   judging period; do not spend additional live submissions.
2. Monitor the portal for any organizer requests without sharing raw data or
   recontacting the family.
3. After the deadline, delete Colab/Drive caches and follow the organizers'
   deletion-confirmation instructions within 30 days.

## Blockers

- No access or submission blocker remains. The remaining work is post-submission
  stewardship and deadline-based deletion.
- FASTQ read-level validation is optional and expensive; the decision to omit
  it must be disclosed if it remains unnecessary for the lead calls.
