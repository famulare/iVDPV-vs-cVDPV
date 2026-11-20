# Verification suite for the curated iVDPV-vs-cVDPV classifier and utilities.
# Run from the repository root: Rscript tests/test_classifier.R
suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv)
  library(pROC)
})

source(file.path("utils", "R", "common.R"))

# --- 1. Sabin VP1 reference loader -------------------------------------------
for (st in c("PV1", "PV2", "PV3")) {
  ref <- read_sabin_vp1(st)
  stopifnot(
    identical(ref$label, unname(sabin_ref_accession[[st]])),
    identical(nchar(ref$nt), unname(vp1_nt_length[[st]])),
    identical(nchar(ref$aa), unname(vp1_nt_length[[st]]) %/% 3L),
    identical(count_ambiguous_aa(ref$aa), 0L)
  )
}
stopifnot(identical(read_sabin1_vp1(), read_sabin_vp1("PV1")))

# --- 2. Codon-level Hamming + in-frame deletion semantics --------------------
ref2 <- read_sabin_vp1("PV2")

s0 <- sequence_stats(ref2$nt, "self", ref2$nt)
stopifnot(
  s0$d_n == 0, s0$d_s == 0, s0$nt_substitutions == 0L,
  s0$indel_codons == 0L, s0$assessable_codons == 301L
)

# 1 synonymous + 1 nonsynonymous codon substitution
q_sub <- ref2$nt
c1 <- substr(ref2$nt, 1, 3)
c2 <- substr(ref2$nt, 4, 6)
syn_c1 <- names(.CODON_TABLE)[.CODON_TABLE == .CODON_TABLE[[c1]] & names(.CODON_TABLE) != c1][[1L]]
nonsyn_c2 <- names(.CODON_TABLE)[.CODON_TABLE != .CODON_TABLE[[c2]] & .CODON_TABLE != "*"][[1L]]
substr(q_sub, 1, 3) <- syn_c1
substr(q_sub, 4, 6) <- nonsyn_c2
s_sub <- sequence_stats(q_sub, "sub", ref2$nt)
stopifnot(s_sub$synonymous_codons == 1L, s_sub$nonsynonymous_sub_codons == 1L, s_sub$d_s == 1, s_sub$d_n == 1)

# 1-codon (3 nt) in-frame internal deletion counts as +1 nonsynonymous difference
q_del3 <- ref2$nt
substr(q_del3, 31, 33) <- "---"
s_del3 <- sequence_stats(q_del3, "del3", ref2$nt)
stopifnot(
  s_del3$indel_codons == 1L, s_del3$indel_events == 1L, !s_del3$has_frameshift,
  s_del3$assessable_codons == 301L, s_del3$d_n == 1, s_del3$d_s == 0,
  s_del3$n_ambiguous_aa == 0L, identical(substr(s_del3$aa_sequence, 11, 11), "-")
)

# 2-codon (6 nt) in-frame internal deletion counts as +2 nonsynonymous differences, 1 event
q_del6 <- ref2$nt
substr(q_del6, 31, 36) <- "------"
s_del6 <- sequence_stats(q_del6, "del6", ref2$nt)
stopifnot(
  s_del6$indel_codons == 2L, s_del6$indel_events == 1L, !s_del6$has_frameshift,
  s_del6$assessable_codons == 301L, s_del6$d_n == 2, s_del6$d_s == 0
)

# Non-multiple-of-3 internal gap flags frameshift
q_fs <- ref2$nt
substr(q_fs, 31, 32) <- "--"
stopifnot(sequence_stats(q_fs, "fs", ref2$nt)$has_frameshift)

# Terminal edge gaps are missing coverage, not deletions, and rescale counts
q_edge <- ref2$nt
substr(q_edge, 1, 9) <- "---------"
c10 <- substr(ref2$nt, 28, 30)
syn_c10 <- names(.CODON_TABLE)[.CODON_TABLE == .CODON_TABLE[[c10]] & names(.CODON_TABLE) != c10][[1L]]
substr(q_edge, 28, 30) <- syn_c10
s_edge <- sequence_stats(q_edge, "edge", ref2$nt)
stopifnot(
  s_edge$indel_codons == 0L, s_edge$comparable_codons == 298L, !s_edge$has_frameshift,
  abs(s_edge$d_s - 301 / 298) < 1e-9, s_edge$d_n == 0
)

# IUPAC ambiguity codes lower denominator rather than creating false mismatches
q_amb <- ref2$nt
substr(q_amb, 10, 12) <- "NNN"
substr(q_amb, 13, 13) <- "R"
s_amb <- sequence_stats(q_amb, "amb", ref2$nt)
stopifnot(
  s_amb$nt_substitutions == 0L, s_amb$d_n == 0, s_amb$d_s == 0,
  s_amb$comparable_codons == 299L
)

# --- 3. Artifact & classification invariants ---------------------------------
art <- readRDS(file.path("analyses", "genbank-classifier", "classifier_artifacts.rds"))
meta <- readRDS(file.path("analyses", "genbank-classifier", "processed_reference_meta.rds"))
tr <- art$training_scores

expected_keys <- c(
  "mod3", "mod5", "training_scores", "thresholds",
  "dnds_contours", "dnyears_contours", "color_map",
  "classifier_color_map", "relative_fitness_norm",
  "model_comparison", "source_manifest"
)
stopifnot(all(expected_keys %in% names(art)))

# Canonical Sabin-parent wild strains must remain wild (no legacy heuristic corruption)
wild_parents <- meta |>
  filter(accession %in% c("AY082689", "V01149", "AY082690", "AY082680", "AY082684", "V01132", "AY082685"))
stopifnot(nrow(wild_parents) == 7L, all(as.character(wild_parents$poliovirus_classification) == "wild"))

# nOPV-related classes are present and labeled; cVDPV-n is included in training as cVDPV
stopifnot(
  all(c("nOPV", "nOPV-L", "VDPV-n", "cVDPV-n") %in% as.character(meta$poliovirus_classification)),
  any(as.character(tr$poliovirus_classification) == "cVDPV-n"),
  all(as.character(tr$classifier_class[tr$poliovirus_classification == "cVDPV-n"]) == "cVDPV")
)

# In-frame VP1 deletions are retained and contribute to d_n
stopifnot(sum(tr$indel_codons > 0L) >= 8L, all(!tr$has_frameshift))

# Prediction works cleanly on new (d_n, d_s) inputs and polarity is P(cVDPV)
probe <- tibble(d_n = c(2, 25), d_s = c(25, 100))
p <- as.numeric(predict(art$mod3, newdata = probe, type = "response"))
stopifnot(p[[1L]] > 0.8, p[[2L]] < 0.2)

cat("All classifier and codon-divergence tests PASSED.\n")
