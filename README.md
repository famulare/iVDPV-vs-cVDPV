# iVDPV vs. cVDPV classifier

## What this project does

Polio eradication depends on telling apart two failure modes for vaccine-derived polioviruses (VDPV):

- **Circulating VDPV (cVDPV)** — Sabin-strain virus that has reverted toward neurovirulence and gained enough community transmission to sustain itself in an undervaccinated population. Outbreak driver.
- **Immunodeficiency-associated VDPV (iVDPV)** — Sabin-strain virus that has continued evolving inside a single immunocompromised host who cannot clear the infection. Long-term shedder.

Both descend from oral polio vaccine. They differ in *how* they evolve. Patient-host evolution looks closer to a nucleotide clock with substantial nonsynonymous drift (the host environment is unusual and there is no transmission bottleneck). Community circulation looks closer to a synonymous clock with purifying selection constraining nonsynonymous change.

This repository builds a classifier that scores a VP1 sequence on the (d_n, d_s) plane — nonsynonymous and synonymous divergence from the matching Sabin reference — and returns a `cVDPVness` probability. The reference set is GenBank-derived Sabin1/2/3 P1 alignments with published cVDPV / iVDPV / wild / Sabin / Sabin-like labels derived from carefully parsing Genbank metadata and cross-referencing to linked publications. The classifier supports all three serotypes.

The intended use is surveillance triage: given a novel post-cessation VDPV isolate of unknown origin, the score and its position relative to the trained reference contours indicate whether the sequence looks more like the cVDPV regime (outbreak source) or the iVDPV regime (long-term-shedder source). A relative-fitness layer additionally projects each sequence's d_n distance from the cVDPV consensus into a fitness-cost-like axis.

This repository contains only the GenBank-trained classifier — reference alignments, training pipeline, baseline plots, and the small utility layer the classifier uses. Country-specific surveillance applications that score post-cessation VP1 isolates against this classifier, and any non-public surveillance data, live in a private continuation that is not published here.

## Classifier artifacts

Run:

```r
Rscript analyses/genbank-classifier/distinguish_i_from_c.R
```

This writes to `analyses/genbank-classifier/`:

- `classifier_artifacts.rds`: `mod3`, `mod5`, training scores, thresholds, contours, color map, and relative-fitness normalization.
- `processed_reference_meta.rds` and `.csv`: the processed reference metadata table before classifier-only filtering.
- `classifier_training_scores.csv`, `classifier_thresholds.csv`, `classifier_dnds_contours.csv`, `classifier_dnyears_contours.csv`, `classifier_color_map.csv`, and `relative_fitness_normalization.csv`.
- Baseline figures: `iVDPV_vs_cVDPV.png`, `relative_fitness_model.png`, `classifier_dnds.png`, `classifier_dnyears.png`, `cVDPVness.png`, and `supp_Zhoa2017_consistency.png`.
