# MVA Hackathon 2026 project state

## Completed work

- Inspected the public challenge page and recorded Track 1/Track 2 submission
  requirements in the repository README.
- Created a data-free, privacy-safe Python ranking scaffold and a fixed proxy
  harness for cautious autoresearch-style iteration.
- Confirmed via public dataset metadata that the VCF is about 315 MB, the
  phenotype document is small, and the FASTQs account for roughly 85 GB.

## Current behavior

- `scripts/inspect_dataset.py` prints safe VCF metadata.
- `scripts/extract_phenotype.py` extracts a private DOCX into ignored
  `runtime/phenotype.txt`.
- `scripts/rank_variants.py` ranks up to 10 VCF records using transparent soft
  priors for established MVA genes, consequence, population frequency, depth,
  genotype quality, and mosaic-compatible VAF.
- `autoresearch/run_experiment.py` runs synthetic regression checks only.

## Decisions

- Keep all gated challenge files and derived candidate outputs out of GitHub.
- Download VCF + phenotype first; use FASTQs only for read-level validation of
  finalists.
- Treat autoresearch as an engineering/proxy loop, not a hidden-answer tuning
  loop; live submissions remain human-reviewed and capped.

## Tests

- Local proxy harness still needs to be run after the initial scaffold commit.
- Real VCF parsing and ranking are pending gated dataset access.

## Next actions

1. Complete participant registration and gated dataset access.
2. Run VCF/phenotype inspection in Colab.
3. Add annotation and HPO-matching stages based on the actual file schema.
4. Validate finalists against read evidence, then prepare the Track 1 report
   and a Track 2 mechanism/drug-repurposing report.

## Blockers

- Dataset access requires the participant's institution, city/country, intended
  use date, and explicit acceptance of the hackathon rules.
- Track 2 cannot be responsibly drafted until the actual causal variant and
  mechanism evidence are identified.

