suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(lubridate)
  library(Biostrings)
})

# Canonical VP1 nt length and Sabin reference accession per serotype.
# All length-sensitive code derives the expected length from the Sabin
# reference (`nchar(sabin$nt)`) or from these lookups directly; there is no
# scalar global, so PV2/PV3 code paths can never silently inherit PV1 defaults.
vp1_nt_length <- c(PV1 = 906L, PV2 = 903L, PV3 = 900L)
sabin_ref_accession <- c(PV1 = "AY184219", PV2 = "AY184220", PV3 = "AY184221")

vp1_nt_length_by_serotype <- function(serotype) {
  if (length(serotype) == 1L) {
    if (!serotype %in% names(vp1_nt_length)) stop("Unknown serotype: ", serotype)
    return(unname(vp1_nt_length[[serotype]]))
  }
  unknown <- setdiff(serotype, names(vp1_nt_length))
  if (length(unknown) > 0L) stop("Unknown serotype(s): ", paste(unknown, collapse = ", "))
  unname(vp1_nt_length[serotype])
}

# Standard genetic code lookup for fast, vectorized codon-level counting.
.CODON_TABLE <- c(
  TTT = "F", TTC = "F", TTA = "L", TTG = "L",
  CTT = "L", CTC = "L", CTA = "L", CTG = "L",
  ATT = "I", ATC = "I", ATA = "I", ATG = "M",
  GTT = "V", GTC = "V", GTA = "V", GTG = "V",
  TCT = "S", TCC = "S", TCA = "S", TCG = "S",
  CCT = "P", CCC = "P", CCA = "P", CCG = "P",
  ACT = "T", ACC = "T", ACA = "T", ACG = "T",
  GCT = "A", GCC = "A", GCA = "A", GCG = "A",
  TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
  CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
  AAT = "N", AAC = "N", AAA = "K", AAG = "K",
  GAT = "D", GAC = "D", GAA = "E", GAG = "E",
  TGT = "C", TGC = "C", TGA = "*", TGG = "W",
  CGT = "R", CGC = "R", CGA = "R", CGG = "R",
  AGT = "S", AGC = "S", AGA = "R", AGG = "R",
  GGT = "G", GGC = "G", GGA = "G", GGG = "G"
)

.UNAMBIGUOUS_BASES <- c("A", "C", "G", "T")

# Rscript and interactive sessions expose script paths differently. These
# helpers keep downstream runners free of hard-coded local filesystem paths.
repo_root <- function() {
  rprojroot::find_root(rprojroot::is_git_root)
}

script_dir <- function() {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) == 1L) {
    dirname(normalizePath(sub("^--file=", "", script_arg)))
  } else {
    normalizePath(getwd())
  }
}

excel_or_text_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(janitor::excel_numeric_to_date(x))

  x_chr <- as.character(x)
  numeric_like <- suppressWarnings(as.numeric(x_chr))
  out <- as.Date(rep(NA_real_, length(x_chr)), origin = "1970-01-01")
  numeric_idx <- !is.na(numeric_like)
  out[numeric_idx] <- janitor::excel_numeric_to_date(numeric_like[numeric_idx])
  text_idx <- !numeric_idx & !is.na(x_chr)
  out[text_idx] <- as.Date(lubridate::parse_date_time(x_chr[text_idx], orders = c("ymd", "mdy", "dmy")))
  out
}

# Convert dates to continuous time for clock fits and x-axes. Leap years are
# handled by measuring each date against its own calendar year length.
decimal_year <- function(date) {
  date <- as.Date(date)
  start <- make_date(year(date), 1, 1)
  finish <- make_date(year(date) + 1, 1, 1)
  year(date) + as.numeric(date - start) / as.numeric(finish - start)
}

decimal_to_date <- function(decimal_date) {
  yr <- floor(decimal_date)
  start <- make_date(yr, 1, 1)
  finish <- make_date(yr + 1, 1, 1)
  start + round((decimal_date - yr) * as.numeric(finish - start))
}

# The classifier can return probabilities numerically equal to 0 or 1 on new
# data. Clamp before logit transformation so downstream plots/tables are finite.
safe_logit <- function(x) {
  x <- pmin(pmax(x, 1e-12), 1 - 1e-12)
  log(x / (1 - x))
}

split_chars <- function(x) {
  strsplit(as.character(x), "", fixed = TRUE)[[1L]]
}

sanitize_tip_label <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

is_env_id <- function(x) {
  startsWith(as.character(x), "ENV-")
}

# Read a plain or gzipped FASTA file into a named character vector (id -> uppercased sequence).
read_fasta_sequences <- function(path) {
  if (!file.exists(path)) stop("FASTA file not found: ", path)
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0L) stop("FASTA file is empty: ", path)

  header_idx <- which(startsWith(lines, ">"))
  if (length(header_idx) == 0L) stop("FASTA file has no headers: ", path)

  ids <- sub("\\s.*$", "", sub("^>", "", lines[header_idx]))
  starts <- header_idx + 1L
  ends <- c(header_idx[-1L] - 1L, length(lines))
  seqs <- vapply(seq_along(ids), function(i) {
    if (starts[[i]] > ends[[i]]) "" else toupper(paste0(lines[starts[[i]]:ends[[i]]], collapse = ""))
  }, character(1L))
  seqs <- chartr("U.", "T-", seqs)
  setNames(seqs, ids)
}

# Translate a Sabin-anchored VP1 sequence while preserving clean in-frame codon
# deletions (`---` -> `-`) and resolving unambiguous IUPAC codons.
translate_vp1_ambiguity_aware <- function(sequence) {
  sequence <- chartr("U.", "T-", toupper(sequence))
  if (nchar(sequence) %% 3L != 0L) stop("Sequence length is not divisible by 3")

  codon_start <- seq(1L, nchar(sequence), by = 3L)
  codons <- substring(sequence, codon_start, codon_start + 2L)
  is_del <- codons == "---"
  is_partial_gap <- !is_del & grepl("-", codons, fixed = TRUE)

  aa <- character(length(codons))
  aa[is_del] <- "-"
  aa[is_partial_gap] <- "X"

  to_translate <- which(!is_del & !is_partial_gap)
  if (length(to_translate) > 0L) {
    sub_codons <- codons[to_translate]
    exact <- unname(.CODON_TABLE[sub_codons])
    unresolved <- which(is.na(exact))
    if (length(unresolved) > 0L) {
      fuzzy_seq <- paste0(sub_codons[unresolved], collapse = "")
      fuzzy_aa <- as.character(Biostrings::translate(
        Biostrings::DNAStringSet(fuzzy_seq),
        if.fuzzy.codon = "solve",
        no.init.codon = TRUE
      ))
      exact[unresolved] <- split_chars(fuzzy_aa)
    }
    aa[to_translate] <- exact
  }

  paste0(aa, collapse = "")
}

count_mismatches <- function(reference, query) {
  sum(split_chars(reference) != split_chars(query))
}

count_ambiguous_aa <- function(aa_sequence) {
  sum(split_chars(aa_sequence) == "X")
}

# Amino-acid change strings are compact and machine-readable:
# reference residue, 1-based VP1 amino-acid position, query residue (or `-` for deleted codon).
aa_change_string <- function(reference_aa, query_aa) {
  ref <- split_chars(reference_aa)
  qry <- split_chars(query_aa)
  pos <- which(ref != qry)
  if (length(pos) == 0L) return("")
  paste0(ref[pos], pos, qry[pos], collapse = ";")
}

aa_distance <- function(a, b) {
  count_mismatches(a, b)
}

nt_distance <- function(a, b) {
  count_mismatches(a, b)
}

group_palette <- function(groups) {
  base <- c("#d95f02", "#2c7fb8", "#1b9e77", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d")
  setNames(base[seq_along(groups)], groups)
}

read_sabin1_vp1 <- function(project_dir = getwd(), data_dir = file.path("data", "curated")) {
  read_sabin_vp1("PV1", project_dir = project_dir, data_dir = data_dir)
}

# Serotype-agnostic Sabin VP1 reference loader. Anchors explicitly on the
# canonical Sabin reference accession (AY184219 / AY184220 / AY184221) inside
# data/curated/PV{1,2,3}_vp1.fasta.gz rather than assuming row ordering.
read_sabin_vp1 <- function(serotype, project_dir = getwd(),
                           data_dir = file.path("data", "curated")) {
  serotype <- match.arg(serotype, c("PV1", "PV2", "PV3"))
  canonical_length <- unname(vp1_nt_length[[serotype]])
  ref_acc <- unname(sabin_ref_accession[[serotype]])

  fasta_path <- file.path(project_dir, data_dir, sprintf("%s_vp1.fasta.gz", serotype))
  seqs <- read_fasta_sequences(fasta_path)
  if (!ref_acc %in% names(seqs)) {
    stop("Sabin reference ", ref_acc, " not found in ", fasta_path)
  }
  nt <- unname(seqs[[ref_acc]])
  if (nchar(nt) != canonical_length) {
    stop(ref_acc, " VP1 length is ", nchar(nt), ", expected ", canonical_length, " nt.")
  }
  list(label = ref_acc, nt = nt, aa = translate_vp1_ambiguity_aware(nt))
}

# Codon-level Hamming divergence against a Sabin-anchored VP1 reference.
#
# Semantics match enterovirus-genbank-curated/site/pipeline/divergence.py:
#   - Leading/trailing `-` outside [first_material, last_material] are missing
#     coverage, not deletions.
#   - Internal `-` inside [first_material, last_material] are deletions; each
#     affected codon counts as +1 nonsynonymous difference (`indel_codons`),
#     and any run whose length is not divisible by 3 sets `has_frameshift = TRUE`.
#   - `comparable_codons`: codons where both reference and query have 3
#     unambiguous A/C/G/T bases. Ambiguity codes reduce the denominator rather
#     than inventing false mismatches.
#   - `assessable_codons = comparable_codons + indel_codons` (disjoint).
#   - `d_n` and `d_s` are codon-level nonsynonymous (including in-frame indel
#     codons) and synonymous counts rescaled to full VP1 codon length so near-
#     complete or ambiguity-bearing VP1 sequences lie on the same scale as
#     complete unambiguous VP1s (where scale == 1.0 exactly).
codon_divergence_stats <- function(ref_nt, query_nt, label = "query") {
  ref_nt <- chartr("U.", "T-", toupper(ref_nt))
  query_nt <- chartr("U.", "T-", toupper(query_nt))

  vp1_nt <- nchar(ref_nt)
  if (vp1_nt %% 3L != 0L) stop("Reference VP1 length (", vp1_nt, ") is not divisible by 3.")
  if (nchar(query_nt) != vp1_nt) {
    stop(label, " has length ", nchar(query_nt), "; expected ", vp1_nt, " nt.")
  }
  bad_chars <- unique(unlist(strsplit(gsub("[ACGTRYSWKMBDHVN-]", "", query_nt), "")))
  if (length(bad_chars) > 0L) {
    stop(label, " contains unsupported nucleotide characters: ", paste(bad_chars, collapse = ", "))
  }

  vp1_codons <- vp1_nt %/% 3L
  r_chr <- split_chars(ref_nt)
  q_chr <- split_chars(query_nt)
  r_is_base <- r_chr %in% .UNAMBIGUOUS_BASES
  q_is_base <- q_chr %in% .UNAMBIGUOUS_BASES
  q_mat <- q_chr != "-"

  if (!any(q_mat)) {
    return(tibble(
      comparable_nt = 0L,
      comparable_codons = 0L,
      indel_codons = 0L,
      indel_events = 0L,
      internal_del_nt = 0L,
      has_frameshift = FALSE,
      assessable_codons = 0L,
      synonymous_codons = 0L,
      nonsynonymous_sub_codons = 0L,
      nonsynonymous_codons = 0L,
      nt_substitutions = 0L,
      nt_diff = NA_real_,
      d_n = NA_real_,
      d_s = NA_real_,
      d_n_rate = NA_real_,
      d_s_rate = NA_real_
    ))
  }

  mat_idx <- which(q_mat)
  span_start <- mat_idx[[1L]]
  span_end <- mat_idx[[length(mat_idx)]]
  pos_idx <- seq_len(vp1_nt)
  int_gap <- (pos_idx >= span_start & pos_idx <= span_end) & (q_chr == "-")

  # Contiguous internal deletion runs
  if (any(int_gap)) {
    gap_rle <- rle(int_gap)
    del_lens <- gap_rle$lengths[gap_rle$values]
    indel_events <- length(del_lens)
    internal_del_nt <- sum(del_lens)
    has_frameshift <- any(del_lens %% 3L != 0L)
  } else {
    indel_events <- 0L
    internal_del_nt <- 0L
    has_frameshift <- FALSE
  }

  # Codon-level classification
  r_cod_mat <- matrix(r_chr, nrow = 3L, ncol = vp1_codons)
  q_cod_mat <- matrix(q_chr, nrow = 3L, ncol = vp1_codons)
  comp_nt_mat <- matrix(r_is_base & q_is_base, nrow = 3L, ncol = vp1_codons)
  int_gap_mat <- matrix(int_gap, nrow = 3L, ncol = vp1_codons)

  comparable_codon <- colSums(comp_nt_mat) == 3L
  indel_codon <- colSums(int_gap_mat) > 0L

  synonymous_codons <- 0L
  nonsynonymous_sub_codons <- 0L
  if (any(comparable_codon)) {
    idx_c <- which(comparable_codon)
    r_codons <- paste0(r_cod_mat[1L, idx_c], r_cod_mat[2L, idx_c], r_cod_mat[3L, idx_c])
    q_codons <- paste0(q_cod_mat[1L, idx_c], q_cod_mat[2L, idx_c], q_cod_mat[3L, idx_c])
    differs <- q_codons != r_codons
    if (any(differs)) {
      same_aa <- unname(.CODON_TABLE[q_codons[differs]]) == unname(.CODON_TABLE[r_codons[differs]])
      synonymous_codons <- as.integer(sum(same_aa))
      nonsynonymous_sub_codons <- as.integer(sum(!same_aa))
    }
  }

  comparable_codons <- as.integer(sum(comparable_codon))
  indel_codons <- as.integer(sum(indel_codon))
  assessable_codons <- comparable_codons + indel_codons
  nonsynonymous_codons <- nonsynonymous_sub_codons + indel_codons

  comparable_nt_mask <- r_is_base & q_is_base
  comparable_nt <- as.integer(sum(comparable_nt_mask))
  nt_substitutions <- as.integer(sum(comparable_nt_mask & (q_chr != r_chr)))

  if (assessable_codons > 0L) {
    codon_scale <- vp1_codons / assessable_codons
    d_n <- nonsynonymous_codons * codon_scale
    d_s <- synonymous_codons * codon_scale
    d_n_rate <- nonsynonymous_codons / assessable_codons
    d_s_rate <- synonymous_codons / assessable_codons
  } else {
    d_n <- NA_real_
    d_s <- NA_real_
    d_n_rate <- NA_real_
    d_s_rate <- NA_real_
  }

  # Total VP1-equivalent nucleotide difference: unambiguous nt substitutions
  # scaled to comparable span plus 1 event per in-frame deleted codon.
  nt_diff <- if (comparable_nt > 0L) {
    nt_substitutions * (vp1_nt / comparable_nt) + indel_codons
  } else {
    NA_real_
  }

  tibble(
    comparable_nt = comparable_nt,
    comparable_codons = comparable_codons,
    indel_codons = indel_codons,
    indel_events = as.integer(indel_events),
    internal_del_nt = as.integer(internal_del_nt),
    has_frameshift = has_frameshift,
    assessable_codons = assessable_codons,
    synonymous_codons = synonymous_codons,
    nonsynonymous_sub_codons = nonsynonymous_sub_codons,
    nonsynonymous_codons = nonsynonymous_codons,
    nt_substitutions = nt_substitutions,
    nt_diff = as.numeric(nt_diff),
    d_n = as.numeric(d_n),
    d_s = as.numeric(d_s),
    d_n_rate = as.numeric(d_n_rate),
    d_s_rate = as.numeric(d_s_rate)
  )
}

# Batch codon-level Hamming divergence table for a named character vector of
# Sabin-anchored VP1 sequences.
codon_divergence_table <- function(ref_nt, query_seqs) {
  ids <- names(query_seqs)
  if (is.null(ids)) stop("query_seqs must be a named character vector")
  map2_dfr(unname(query_seqs), ids, function(seq, id) {
    bind_cols(tibble(accession = id), codon_divergence_stats(ref_nt, seq, label = id))
  })
}

# All classifier-style sequence features are derived here so the reference
# classifier and any downstream surveillance analyses share the exact same
# codon-level Hamming + in-frame deletion definitions of d_n, d_s, and nt_diff.
sequence_stats <- function(sequence, label, sabin_vp1_nt, sabin_vp1_aa = NULL) {
  div <- codon_divergence_stats(sabin_vp1_nt, sequence, label = label)
  aa <- translate_vp1_ambiguity_aware(sequence)
  tibble(
    aa_sequence = aa,
    vp1_aa_length = nchar(aa),
    n_ambiguous_aa = count_ambiguous_aa(aa),
    aa_changes = div$nonsynonymous_codons
  ) |>
    bind_cols(div)
}
