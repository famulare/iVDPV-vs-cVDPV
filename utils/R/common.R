suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(lubridate)
  library(Biostrings)
  library(ape)
})

# Canonical VP1 nt length differs by serotype: PV1=906, PV2=903, PV3=900. The
# per-serotype value is the length of natural VP1 in cVDPV FASTAs (cf. the
# Sabin alignments under data/genbank/, which are P1-region alignments padded
# to a common 906 nt). All length-sensitive code derives the expected length
# from the Sabin reference (`nchar(sabin$nt)`) or from this lookup directly;
# there is no scalar global, deliberately, so PV2/PV3 code paths can never
# silently inherit a PV1 default.
vp1_nt_length <- c(PV1 = 906L, PV2 = 903L, PV3 = 900L)

vp1_nt_length_by_serotype <- function(serotype) {
  if (length(serotype) == 1L) {
    if (!serotype %in% names(vp1_nt_length)) stop("Unknown serotype: ", serotype)
    return(unname(vp1_nt_length[[serotype]]))
  }
  unknown <- setdiff(serotype, names(vp1_nt_length))
  if (length(unknown) > 0L) stop("Unknown serotype(s): ", paste(unknown, collapse = ", "))
  unname(vp1_nt_length[serotype])
}

# Rscript and interactive sessions expose script paths differently. These
# helpers keep the country runners free of hard-coded local filesystem paths.
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

# The Sabin reference alignments are PHYLIP files with Sabin as the first row;
# downstream country analyses only need that first reference sequence.
read_first_phylip_sequence <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) < 2L) stop("PHYLIP file has no first sequence row: ", path)
  header_parts <- strsplit(trimws(lines[[1L]]), "\\s+")[[1L]]
  expected_length <- as.integer(header_parts[[2L]])
  first_parts <- strsplit(trimws(lines[[2L]]), "\\s+")[[1L]]
  sequence <- toupper(paste0(first_parts[-1L], collapse = ""))
  if (!is.na(expected_length) && nchar(sequence) != expected_length) {
    stop("First PHYLIP sequence length is ", nchar(sequence), ", expected ", expected_length)
  }
  list(label = first_parts[[1L]], sequence = sequence)
}

# Translate VP1 while preserving useful IUPAC information. Ambiguous codons
# become a concrete amino acid when all possible codons agree; otherwise they
# become X so the deduplication/filtering logic can make an explicit choice.
translate_vp1_ambiguity_aware <- function(sequence) {
  sequence <- toupper(sequence)
  if (nchar(sequence) %% 3L != 0L) stop("Sequence length is not divisible by 3")

  codon_start <- seq(1L, nchar(sequence), by = 3L)
  codons <- substring(sequence, codon_start, codon_start + 2L)
  codons <- if_else(str_detect(codons, "-"), "NNN", codons)
  sequence_for_translation <- str_c(codons, collapse = "")

  Biostrings::translate(
    Biostrings::DNAStringSet(sequence_for_translation),
    if.fuzzy.codon = "solve",
    no.init.codon = TRUE
  ) |>
    as.character()
}

count_mismatches <- function(reference, query) {
  sum(split_chars(reference) != split_chars(query))
}

count_ambiguous_aa <- function(aa_sequence) {
  sum(split_chars(aa_sequence) == "X")
}

# Amino-acid change strings are deliberately compact and machine-readable:
# reference residue, 1-based VP1 amino-acid position, query residue.
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

read_sabin1_vp1 <- function(project_dir, data_dir = file.path("data", "genbank")) {
  read_sabin_vp1("PV1", project_dir = project_dir, data_dir = data_dir)
}

# Serotype-agnostic Sabin VP1 reference loader. Supports PV1/PV2/PV3 by mapping
# to the matching Sabin{N}.ntP1 PHYLIP alignment in data/genbank/. Trims the
# trailing per-serotype canonical VP1 length (PV1=906, PV2=903, PV3=900) so the
# returned reference matches the natural cVDPV FASTA lengths for each serotype.
read_sabin_vp1 <- function(serotype, project_dir = getwd(),
                           data_dir = file.path("data", "genbank")) {
  serotype <- match.arg(serotype, c("PV1", "PV2", "PV3"))
  canonical_length <- unname(vp1_nt_length[[serotype]])
  sabin_n <- sub("^PV", "", serotype)
  sabin <- read_first_phylip_sequence(file.path(
    project_dir, data_dir,
    sprintf("Sabin%s.ntP1.alignment.final.phylip.reduced", sabin_n)
  ))
  nt <- substr(sabin$sequence,
               nchar(sabin$sequence) - canonical_length + 1L,
               nchar(sabin$sequence))
  list(label = sabin$label, nt = nt, aa = translate_vp1_ambiguity_aware(nt))
}

# All classifier-style sequence features are derived here so MAD, RDC, and any
# future country analysis share the same definitions of d_n, d_s, and X counts.
sequence_stats <- function(sequence, label, sabin_vp1_nt, sabin_vp1_aa) {
  expected_nt_length <- nchar(sabin_vp1_nt)
  if (nchar(sequence) != expected_nt_length) {
    stop(label, " has length ", nchar(sequence),
         "; expected ", expected_nt_length, " nt.")
  }
  bad_chars <- unique(unlist(strsplit(gsub("[ACGTRYSWKMBDHVN-]", "", sequence), "")))
  if (length(bad_chars) > 0L) {
    stop(label, " contains unsupported nucleotide characters: ", paste(bad_chars, collapse = ", "))
  }
  aa <- translate_vp1_ambiguity_aware(sequence)
  nt_substitutions <- count_mismatches(sabin_vp1_nt, sequence)
  aa_changes <- count_mismatches(sabin_vp1_aa, aa)
  tibble(
    aa_sequence = aa,
    vp1_aa_length = nchar(aa),
    n_ambiguous_aa = count_ambiguous_aa(aa),
    nt_substitutions = nt_substitutions,
    aa_changes = aa_changes,
    d_n = aa_changes,
    d_s = nt_substitutions - aa_changes
  )
}
