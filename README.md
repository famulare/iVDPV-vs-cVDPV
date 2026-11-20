# iVDPV vs. cVDPV classifier

## What this project does

Polio eradication depends on telling apart two failure modes for vaccine-derived polioviruses (VDPV):

- **Circulating VDPV (cVDPV)** — OPV-derived virus (Sabin or nOPV) that has reverted toward neurovirulence and gained enough community transmission to sustain itself in an undervaccinated population. Outbreak driver.
- **Immunodeficiency-associated VDPV (iVDPV)** — OPV-derived virus that has continued evolving inside a single immunocompromised host who cannot clear the infection. Long-term shedder.

Both descend from oral polio vaccine, but they differ in *how* they evolve. Intrahost evolution looks closer to a nucleotide clock with substantial nonsynonymous drift (the host environment is unusual and there is no transmission bottleneck). Community circulation looks closer to a synonymous clock with purifying selection constraining nonsynonymous change.

This repository builds a classifier that scores a Sabin-anchored VP1 sequence on the `(d_n, d_s)` plane — codon-level nonsynonymous (including in-frame codon deletions) and synonymous Hamming differences from the matching Sabin reference (`AY184219` / `AY184220` / `AY184221`) — and returns a `cVDPVness` probability across all three serotypes (`PV1`, `PV2`, `PV3`).

The intended use is surveillance triage: given a novel post-cessation VDPV isolate of unknown origin, the score and its position relative to the trained reference contours indicate whether the sequence looks more like the cVDPV regime (outbreak source) or the iVDPV regime (long-term-shedder source). A relative-fitness layer additionally projects each sequence's `d_n` distance from the cVDPV consensus into a fitness-cost-like axis.

This repository contains only the public curated-GenBank classifier — extracted reference VP1 alignments, canonical metadata, training pipeline, baseline plots, tests, and the shared utility layer in `utils/R/common.R`. Country-specific surveillance applications that score post-cessation VP1 isolates against this classifier, and any non-public surveillance data, live in a private continuation not published here.

## Canonical data source & refresh workflow

Reference VP1 alignments and epidemiology-first metadata are extracted from [`enterovirus-genbank-curated`](https://famulare.github.io/enterovirus-genbank-curated/) (currently pinned to **Release 5.3.0**, schema `2.4.1`; full provenance and SHA-256 digests in `data/curated/source_manifest.json`):

- `data/curated/PV{1,2,3}_vp1.fasta.gz` — Sabin-anchored codon-aligned VP1 nucleotide sequences (`906` / `903` / `900` nt), with `AY184219` / `AY184220` / `AY184221` placed first in each file.
- `data/curated/poliovirus_vp1_metadata.tsv.gz` — All 26 canonical vouched poliovirus columns (`virus_type`, `poliovirus_classification`, `surveillance_stream`, `sample_origin`, `specimen_type`, `collection_date*`, `country/admin1/locality`, `engineered_or_construct`, `sequence_sha256`, etc.), linked citation fields (`pubmed_id`, `reference_title`), and VP1 alignment QC metrics.
- `data/reference/zhao_fig2_extracted_vpdv_points.csv` — Digitized inter-isolate points from Zhao et al. (2017) Fig. 2 used for the supplementary consistency check.

### Divergence & inclusion conventions

- **Codon-level Hamming counting (`utils/R/common.R::codon_divergence_stats`)**:
  - `comparable_codons`: VP1 codons where both Sabin reference and query carry three unambiguous `A/C/G/T` bases. IUPAC ambiguity codes reduce the denominator instead of creating false mismatches.
  - **In-frame deletions as difference signal**: internal deletions inside a sequence's covered VP1 span count per affected codon (`indel_codons`) toward `d_n` and `assessable_codons = comparable_codons + indel_codons`. Leading/trailing gaps are treated as missing coverage, not deletions. Sequences with frameshifts (non-multiple-of-3 internal gaps) or internal stop codons are excluded.
  - `d_n` and `d_s` are rescaled to full VP1 codon length (`vp1_codons / assessable_codons`, identically `1.0` for complete unambiguous VP1s and complete VP1s with in-frame deletions).
- **Retained canonical `poliovirus_classification` classes** (vouched, `engineered_or_construct == FALSE`, VP1 coverage ≥ 95%): `Sabin-like`, `nOPV`, `nOPV-L`, `VDPV`, `VDPV-n`, `iVDPV`, `cVDPV`, `cVDPV-n`, and `wild` (**7,539** reference isolates in Release 5.3.0).
- **Classifier training set**: `cVDPV` + `cVDPV-n` (`classifier_class = "cVDPV"`, **1,382** isolates: `PV1` 291, `PV2` 1,082 incl. 25 `cVDPV-n`, `PV3` 9) vs. `iVDPV` (`classifier_class = "iVDPV"`, **179** isolates: `PV1` 50, `PV2` 71, `PV3` 58) — **1,561** total training isolates, including 8 with in-frame 1-codon VP1 deletions. All `nOPV`-related classes retain their canonical `poliovirus_classification` label in tables and figures.

### Refreshing from a new `enterovirus-genbank-curated` release

```bash
python3 scripts/extract_curated_data.py --curated-dir /path/to/enterovirus-genbank-curated
Rscript analyses/genbank-classifier/distinguish_i_from_c.R
Rscript tests/test_classifier.R
```

## Classifier artifacts

Running `Rscript analyses/genbank-classifier/distinguish_i_from_c.R` writes to `analyses/genbank-classifier/`:

- `classifier_artifacts.rds`: `mod3` (`s(d_n) + s(d_s)`, pooled AUC = 0.945), `mod5` (`s(d_n) + s(years)`, AUC = 0.858), `mod_factor_dnds` (`s(d_n) + s(d_s) + virus_type`, overall AUC = 0.961), `mod_independent_dnds` (independent `PV1`/`PV2`/`PV3` `s(d_n) + s(d_s)` models, overall AUC = 0.965; per-serotype AUCs: `PV1` 0.987, `PV2` 0.927, `PV3` 0.882), `training_scores`, `thresholds`, `serotype_thresholds`, `dnds_contours`, `dnyears_contours`, `serotype_dnds_contours`, `color_map`, `classifier_color_map`, `relative_fitness_norm`, `model_comparison`, `serotype_model_comparison`, and `source_manifest`. *(Load `library(mgcv)` before calling `predict(artifacts$mod3, ...)`.)*
- `processed_reference_meta.rds` and `.csv`: the full processed reference metadata and codon-divergence table across all retained canonical classes.
- `classifier_training_scores.csv`, `classifier_thresholds.csv`, `classifier_serotype_thresholds.csv`, `classifier_serotype_model_comparison.csv`, `classifier_dnds_contours.csv`, `classifier_dnyears_contours.csv`, `classifier_serotype_dnds_contours.csv`, `classifier_color_map.csv`, and `relative_fitness_normalization.csv`.
- Baseline figures: `iVDPV_vs_cVDPV.png`, `relative_fitness_model.png`, `classifier_dnds.png`, `classifier_dnyears.png`, `classifier_serotype_dnds.png`, `cVDPVness.png`, and `supp_Zhao2017_consistency.png`.
