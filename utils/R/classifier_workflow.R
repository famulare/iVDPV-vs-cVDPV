suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(patchwork)
  library(mgcv)
  library(pROC)
  library(ggbeeswarm)
})

source(file.path("utils", "R", "common.R"))

classifier_config <- function(project_dir = getwd(),
                              data_dir = file.path("data", "curated"),
                              reference_dir = file.path("data", "reference"),
                              output_dir = file.path("analyses", "genbank-classifier")) {
  list(
    project_dir = project_dir,
    data_dir = data_dir,
    reference_dir = reference_dir,
    output_dir = output_dir,
    serotypes = tibble(
      virus_type = c("PV1", "PV2", "PV3"),
      ref_accession = unname(sabin_ref_accession[c("PV1", "PV2", "PV3")]),
      vp1_nt = unname(vp1_nt_length[c("PV1", "PV2", "PV3")])
    ),
    min_vp1_coverage = 0.95,
    allow_frameshift = FALSE,
    allow_internal_stop = FALSE,
    # Canonical poliovirus_classification values retained in processed_reference_meta
    # (Sabin anchor records serve as the zero-divergence reference and are excluded
    # from the post-Sabin reference table; nOPV-related classes are included & labeled).
    retained_classes = c(
      "Sabin-like", "nOPV", "nOPV-L",
      "VDPV", "VDPV-n",
      "iVDPV",
      "cVDPV", "cVDPV-n",
      "wild"
    ),
    cvdpv_classes = c("cVDPV", "cVDPV-n"),
    ivdpv_classes = c("iVDPV"),
    projection_classes = c("Sabin-like", "nOPV", "nOPV-L", "VDPV", "VDPV-n", "iVDPV", "cVDPV", "cVDPV-n"),
    color_map = c(
      "Sabin-like" = "#04e762",
      "nOPV"       = "#008f5a",
      "nOPV-L"     = "#57cc99",
      "VDPV"       = "#ae7766",
      "VDPV-n"     = "#d49a6a",
      "iVDPV"      = "#00a1e4",
      "cVDPV"      = "#dc0073",
      "cVDPV-n"    = "#ff6b6b",
      "wild"       = "#f5b700"
    ),
    classifier_color_map = c(
      iVDPV = "#00a1e4",
      cVDPV = "#dc0073"
    ),
    clock_rates = list(
      ivdpv_nt_per_site_yr = 0.0072,
      syn_per_site_yr = 0.01
    ),
    relative_fitness = list(
      consensus_d_n = 3,
      base = 0.95,
      genome_nt = 2200,
      vp1_nt = 900
    )
  )
}

setup_classifier_graphics <- function() {
  theme_set(theme_bw())
  if (!interactive()) {
    grDevices::pdf(NULL)
    on.exit({
      while (grDevices::dev.cur() > 1) grDevices::dev.off()
    }, add = TRUE)
  }
  invisible(NULL)
}

read_curated_metadata <- function(config) {
  meta_path <- file.path(config$project_dir, config$data_dir, "poliovirus_vp1_metadata.tsv.gz")
  if (!file.exists(meta_path)) stop("Curated metadata not found: ", meta_path)

  readr::read_tsv(
    meta_path,
    col_types = readr::cols(.default = readr::col_character()),
    na = character(),
    show_col_types = FALSE
  ) |>
    mutate(
      sequence_length_nt = as.integer(sequence_length_nt),
      vp1_ref_nt = as.integer(vp1_ref_nt),
      vp1_non_gap_nt = as.integer(vp1_non_gap_nt),
      vp1_acgt_nt = as.integer(vp1_acgt_nt),
      vp1_ambig_nt = as.integer(vp1_ambig_nt),
      vp1_coverage_frac = as.numeric(vp1_coverage_frac),
      vp1_internal_del_nt = as.integer(vp1_internal_del_nt),
      vp1_internal_del_events = as.integer(vp1_internal_del_events),
      vp1_internal_del_codons = as.integer(vp1_internal_del_codons),
      engineered_or_construct = engineered_or_construct == "TRUE",
      is_sabin_reference = is_sabin_reference == "TRUE",
      vp1_has_frameshift = vp1_has_frameshift == "TRUE",
      vp1_has_misphased_gap = vp1_has_misphased_gap == "TRUE",
      vp1_has_internal_stop = vp1_has_internal_stop == "TRUE"
    )
}

read_source_manifest <- function(config) {
  manifest_path <- file.path(config$project_dir, config$data_dir, "source_manifest.json")
  if (!file.exists(manifest_path)) return(list())
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::fromJSON(manifest_path, simplifyVector = TRUE)
  } else {
    list(manifest_path = manifest_path)
  }
}

serotype_reference_meta <- function(spec, all_meta, config) {
  vt <- spec$virus_type
  ref_acc <- spec$ref_accession
  vp1_nt <- spec$vp1_nt

  fasta_path <- file.path(config$project_dir, config$data_dir, sprintf("%s_vp1.fasta.gz", vt))
  seqs <- read_fasta_sequences(fasta_path)
  if (!ref_acc %in% names(seqs)) {
    stop("Sabin reference ", ref_acc, " missing from ", fasta_path)
  }
  ref_nt <- unname(seqs[[ref_acc]])
  if (nchar(ref_nt) != vp1_nt) {
    stop("Sabin reference ", ref_acc, " length is ", nchar(ref_nt), ", expected ", vp1_nt)
  }

  meta <- all_meta |>
    filter(
      virus_type == vt,
      curation_status == "vouched",
      !engineered_or_construct,
      !is_sabin_reference,
      poliovirus_classification %in% config$retained_classes,
      vp1_coverage_frac >= config$min_vp1_coverage,
      accession %in% names(seqs)
    )

  if (!config$allow_frameshift) {
    meta <- meta |> filter(!vp1_has_frameshift)
  }
  if (!config$allow_internal_stop) {
    meta <- meta |> filter(!vp1_has_internal_stop)
  }

  query_seqs <- seqs[meta$accession]
  distances <- codon_divergence_table(ref_nt, query_seqs)

  rates <- config$clock_rates
  meta |>
    left_join(distances, by = "accession") |>
    mutate(
      # iVDPV elapsed time follows a nucleotide clock because intrahost evolution
      # includes substantial nonsynonymous drift; community circulation / Sabin-
      # derived / wild groups follow the synonymous clock proxy.
      years = if_else(
        poliovirus_classification %in% config$ivdpv_classes,
        nt_diff / (vp1_nt * rates$ivdpv_nt_per_site_yr),
        d_s / (vp1_nt * rates$syn_per_site_yr)
      )
    )
}

prepare_reference_meta <- function(config) {
  all_meta <- read_curated_metadata(config)
  meta <- pmap_dfr(config$serotypes, \(virus_type, ref_accession, vp1_nt) {
    serotype_reference_meta(
      tibble(virus_type = virus_type, ref_accession = ref_accession, vp1_nt = vp1_nt),
      all_meta = all_meta,
      config = config
    )
  })

  classification_levels <- names(config$color_map)
  rf <- config$relative_fitness

  meta |>
    mutate(
      virus_type = factor(virus_type, levels = c("PV1", "PV2", "PV3")),
      poliovirus_classification = factor(poliovirus_classification, levels = classification_levels),
      classifier_class = case_when(
        poliovirus_classification %in% config$ivdpv_classes ~ "iVDPV",
        poliovirus_classification %in% config$cvdpv_classes ~ "cVDPV",
        TRUE ~ NA_character_
      ),
      classifier_class = factor(classifier_class, levels = c("iVDPV", "cVDPV")),
      vaccine_lineage = case_when(
        poliovirus_classification %in% c("nOPV", "nOPV-L", "VDPV-n", "cVDPV-n") ~ "nOPV",
        poliovirus_classification == "wild" ~ "wild",
        TRUE ~ "Sabin"
      ),
      d_n_from_consensus = abs(d_n - rf$consensus_d_n),
      relative_fitness = rf$base^(rf$genome_nt / rf$vp1_nt * d_n_from_consensus)
    ) |>
    mutate(
      relative_fitness = relative_fitness / mean(relative_fitness[classifier_class == "cVDPV"], na.rm = TRUE)
    )
}

classifier_training_data <- function(meta, config) {
  meta |>
    filter(!is.na(classifier_class)) |>
    ungroup() |>
    mutate(
      poliovirus_classification = droplevels(poliovirus_classification),
      classifier_class = factor(classifier_class, levels = c("iVDPV", "cVDPV"))
    )
}

plot_reference_overview <- function(meta, training, config) {
  colors <- config$color_map
  cls_colors <- config$classifier_color_map

  p_all <- ggplot(meta, aes(x = d_s, y = d_n, color = poliovirus_classification, shape = virus_type)) +
    geom_point(alpha = 0.25) +
    ylab("non-synonymous codon differences") +
    xlab("synonymous codon differences") +
    scale_color_manual(values = colors, name = "classification", drop = FALSE) +
    scale_shape_discrete(name = "virus_type", drop = FALSE) +
    guides(
      color = guide_legend(override.aes = list(alpha = 1, size = 2.2)),
      shape = guide_legend(override.aes = list(alpha = 1, size = 2.2))
    )

  p_years <- ggplot(training, aes(x = years, y = d_n)) +
    geom_point(aes(color = poliovirus_classification, shape = virus_type), alpha = 0.25) +
    geom_smooth(
      aes(group = classifier_class, color = classifier_class),
      method = "gam", fill = "black", alpha = 0.25, show.legend = FALSE
    ) +
    ylab("non-synonymous codon differences") +
    xlab("approx. years since OPV") +
    scale_x_continuous(breaks = seq(0, 35, by = 5)) +
    scale_color_manual(values = c(colors, cls_colors), guide = "none") +
    scale_shape_discrete(guide = "none")

  p_dnds <- ggplot(training, aes(x = d_s, y = d_n)) +
    geom_point(aes(color = poliovirus_classification, shape = virus_type), alpha = 0.25) +
    geom_smooth(
      aes(group = classifier_class, color = classifier_class),
      method = "gam", fill = "black", alpha = 0.25, show.legend = FALSE
    ) +
    ylab("non-synonymous codon differences") +
    xlab("synonymous codon differences") +
    scale_x_continuous(breaks = seq(0, 300, by = 20)) +
    scale_color_manual(values = c(colors, cls_colors), guide = "none") +
    scale_shape_discrete(guide = "none")

  ggsave(
    file.path(config$project_dir, config$output_dir, "iVDPV_vs_cVDPV.png"),
    p_years + p_dnds + p_all,
    units = "in", width = 14, height = 4.2
  )
}

plot_relative_fitness_model <- function(training, config) {
  sabin_df <- tibble(
    virus_type = factor(c("PV1", "PV2", "PV3"), levels = c("PV1", "PV2", "PV3")),
    relative_fitness = 2.3 / c(14, 8, 18)
  )

  p <- ggplot(training, aes(x = years, y = relative_fitness)) +
    geom_point(aes(color = poliovirus_classification, shape = virus_type), alpha = 0.25) +
    geom_smooth(
      aes(group = classifier_class, color = classifier_class),
      method = "gam", fill = "black", alpha = 0.25, show.legend = FALSE
    ) +
    geom_hline(data = sabin_df, aes(yintercept = relative_fitness, linetype = virus_type)) +
    ylab("relative fitness") +
    xlab("approx. years since OPV") +
    scale_y_continuous(breaks = seq(0, 1.2, by = 0.2)) +
    scale_color_manual(values = config$color_map, name = "classification") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2.2)))

  ggsave(
    file.path(config$project_dir, config$output_dir, "relative_fitness_model.png"),
    p, units = "in", width = 6.0, height = 3.8
  )
}

plot_zhao_consistency <- function(training, config) {
  zhao_path <- file.path(config$project_dir, config$reference_dir, "zhao_fig2_extracted_vpdv_points.csv")
  zhao <- read.csv(zhao_path)
  names(zhao) <- c("nt_identity_pct", "aa_identity_pct", "classifier_class")
  zhao <- zhao |>
    mutate(classifier_class = factor(classifier_class, levels = c("iVDPV", "cVDPV")))

  training_pct <- training |>
    mutate(
      nt_identity_pct = 100 * (1 - nt_diff / vp1_ref_nt),
      aa_identity_pct = 100 * (1 - d_n / (vp1_ref_nt / 3))
    )

  colors <- config$classifier_color_map
  p <- ggplot() +
    geom_point(
      data = training_pct,
      aes(x = nt_identity_pct, y = aa_identity_pct, color = colors[["iVDPV"]]),
      alpha = 0.4
    ) +
    geom_point(
      data = zhao,
      aes(x = nt_identity_pct, y = aa_identity_pct, color = colors[["cVDPV"]]),
      alpha = 0.6
    ) +
    facet_grid(. ~ classifier_class) +
    xlab("nucleotide identity (%)") +
    ylab("amino acid identity (%)") +
    scale_color_identity(
      name = "source",
      breaks = c(colors[["iVDPV"]], colors[["cVDPV"]]),
      labels = c("Curated distance from Sabin", "Zhao 2017\ninter-isolate distance"),
      guide = "legend"
    )

  ggsave(
    file.path(config$project_dir, config$output_dir, "supp_Zhao2017_consistency.png"),
    p, units = "in", width = 8, height = 3.5
  )
}

fit_classifier_models <- function(training) {
  mod1 <- gam(classifier_class ~ d_n + d_s, data = training, family = "binomial")
  mod2 <- gam(classifier_class ~ d_n + d_s + virus_type, data = training, family = "binomial")
  mod3 <- gam(classifier_class ~ s(d_n) + s(d_s), data = training, family = "binomial")
  mod4 <- gam(classifier_class ~ ti(d_n, d_s), data = training, family = "binomial")
  mod5 <- gam(classifier_class ~ s(d_n) + s(years), data = training, family = "binomial")

  # Serotype-stratified (d_n, d_s) models:
  #   1. factor_dnds: shared smooths s(d_n) + s(d_s) with virus_type as a factor intercept
  #   2. independent_dnds: separate s(d_n) + s(d_s) GAM fit independently per serotype
  mod_factor_dnds <- gam(classifier_class ~ s(d_n) + s(d_s) + virus_type, data = training, family = "binomial")

  serotype_levels <- levels(training$virus_type)
  mod_independent_dnds <- map(serotype_levels, function(vt) {
    sub <- training |> filter(virus_type == vt)
    k_dn <- max(3L, min(10L, n_distinct(sub$d_n) - 1L))
    k_ds <- max(3L, min(10L, n_distinct(sub$d_s) - 1L))
    gam(classifier_class ~ s(d_n, k = k_dn) + s(d_s, k = k_ds), data = sub, family = "binomial")
  }) |>
    set_names(serotype_levels)

  ind_loglik <- sum(vapply(mod_independent_dnds, function(m) as.numeric(logLik(m)), numeric(1L)))
  ind_edf <- sum(vapply(mod_independent_dnds, function(m) attr(logLik(m), "df"), numeric(1L)))
  ind_aic <- sum(vapply(mod_independent_dnds, AIC, numeric(1L)))

  serotype_comparison <- bind_rows(
    tibble(
      model = c(
        "pooled: s(d_n) + s(d_s)",
        "serotype_factor: s(d_n) + s(d_s) + virus_type",
        "independent_serotype: sum(PV1 + PV2 + PV3)"
      ),
      scope = "All",
      n = nrow(training),
      n_iVDPV = sum(training$classifier_class == "iVDPV"),
      n_cVDPV = sum(training$classifier_class == "cVDPV"),
      edf = c(attr(logLik(mod3), "df"), attr(logLik(mod_factor_dnds), "df"), ind_edf),
      logLik = c(as.numeric(logLik(mod3)), as.numeric(logLik(mod_factor_dnds)), ind_loglik),
      aic = c(AIC(mod3), AIC(mod_factor_dnds), ind_aic)
    ),
    map_dfr(serotype_levels, function(vt) {
      m <- mod_independent_dnds[[vt]]
      sub <- training |> filter(virus_type == vt)
      tibble(
        model = paste0("independent_", vt, ": s(d_n) + s(d_s)"),
        scope = vt,
        n = nrow(sub),
        n_iVDPV = sum(sub$classifier_class == "iVDPV"),
        n_cVDPV = sum(sub$classifier_class == "cVDPV"),
        edf = attr(logLik(m), "df"),
        logLik = as.numeric(logLik(m)),
        aic = AIC(m)
      )
    })
  )

  list(
    linear = mod1,
    linear_serotype = mod2,
    dnds = mod3,
    interaction = mod4,
    dnyears = mod5,
    factor_dnds = mod_factor_dnds,
    independent_dnds = mod_independent_dnds,
    model_comparison = tibble(
      model = c(
        "linear d_n+d_s",
        "linear + virus_type",
        "smooth d_n+d_s (pooled)",
        "tensor interaction",
        "smooth d_n+years",
        "smooth d_n+d_s + virus_type (factor)",
        "smooth d_n+d_s independent by virus_type"
      ),
      aic = c(AIC(mod1), AIC(mod2), AIC(mod3), AIC(mod4), AIC(mod5), AIC(mod_factor_dnds), ind_aic)
    ),
    serotype_model_comparison = serotype_comparison
  )
}

predict_independent_serotype <- function(ind_models, newdata) {
  out <- rep(NA_real_, nrow(newdata))
  for (vt in names(ind_models)) {
    idx <- which(as.character(newdata$virus_type) == vt)
    if (length(idx) > 0L) {
      out[idx] <- as.numeric(predict(ind_models[[vt]], newdata = newdata[idx, , drop = FALSE], type = "response"))
    }
  }
  out
}

score_training_data <- function(training, models) {
  training |>
    mutate(
      cVDPVness = as.numeric(predict(models$dnds, newdata = training, type = "response")),
      cVDPVness_years = as.numeric(predict(models$dnyears, newdata = training, type = "response")),
      cVDPVness_factor = as.numeric(predict(models$factor_dnds, newdata = training, type = "response")),
      cVDPVness_independent = predict_independent_serotype(models$independent_dnds, training)
    )
}

classifier_thresholds <- function(scored_training) {
  dnds_roc <- roc(
    response = scored_training$classifier_class,
    predictor = scored_training$cVDPVness,
    levels = c("iVDPV", "cVDPV"),
    direction = "<",
    quiet = TRUE
  )
  dnyears_roc <- roc(
    response = scored_training$classifier_class,
    predictor = scored_training$cVDPVness_years,
    levels = c("iVDPV", "cVDPV"),
    direction = "<",
    quiet = TRUE
  )

  # Use ret = "all" / first row if tied so output is always 1 row per operating point
  pick_one <- function(df) as_tibble(df)[1L, , drop = FALSE]
  youden <- pick_one(coords(dnds_roc, "best", best.method = "youden"))
  weighted <- pick_one(coords(dnds_roc, "best", best.method = "youden", best.weights = c(10, 0.5)))
  youden_years <- pick_one(coords(dnyears_roc, "best", best.method = "youden"))
  weighted_years <- pick_one(coords(dnyears_roc, "best", best.method = "youden", best.weights = c(10, 0.5)))

  tibble(
    specificity = c(youden$specificity, weighted$specificity, youden_years$specificity, weighted_years$specificity),
    sensitivity = c(youden$sensitivity, weighted$sensitivity, youden_years$sensitivity, weighted_years$sensitivity),
    threshold = c(youden$threshold, weighted$threshold, youden_years$threshold, weighted_years$threshold),
    auc = c(as.numeric(auc(dnds_roc)), as.numeric(auc(dnds_roc)), as.numeric(auc(dnyears_roc)), as.numeric(auc(dnyears_roc))),
    cost = c("equal", "c/i = 10x", "equal", "c/i = 10x"),
    model = c("s(dn) + s(ds)", "s(dn) + s(ds)", "s(dn) + s(years)", "s(dn) + s(years)")
  )
}

contour_table <- function(model, grid, x_col, y_col, levels, costs) {
  x_values <- sort(unique(grid[[x_col]]))
  y_values <- sort(unique(grid[[y_col]]))
  zmat <- matrix(
    as.numeric(predict(model, newdata = grid, type = "response")),
    nrow = length(x_values),
    ncol = length(y_values)
  )

  comp_id <- 0L
  map2_dfr(levels, costs, function(lvl, cst) {
    contours <- contourLines(x = x_values, y = y_values, z = zmat, levels = lvl)
    if (length(contours) == 0L) return(tibble())
    map_dfr(seq_along(contours), function(i) {
      comp_id <<- comp_id + 1L
      tibble(
        comp = comp_id,
        !!x_col := contours[[i]]$x,
        !!y_col := contours[[i]]$y,
        cost = cst
      )
    })
  })
}

build_classifier_surfaces <- function(models, thresholds, training) {
  max_ds <- max(140, ceiling(max(training$d_s, na.rm = TRUE) / 10) * 10 + 10)
  max_dn <- max(50, ceiling(max(training$d_n, na.rm = TRUE) / 5) * 5 + 5)
  max_years <- max(30, ceiling(max(training$years, na.rm = TRUE) / 5) * 5 + 5)

  dnds_grid <- expand.grid(d_s = seq(0, max_ds, by = 1), d_n = seq(0, max_dn, by = 1))
  dnyears_grid <- expand.grid(years = seq(0, max_years, by = 0.1), d_n = seq(0, max_dn, by = 1))

  list(
    dnds_contours = contour_table(
      models$dnds,
      dnds_grid,
      x_col = "d_s",
      y_col = "d_n",
      levels = thresholds$threshold[thresholds$model == "s(dn) + s(ds)"],
      costs = thresholds$cost[thresholds$model == "s(dn) + s(ds)"]
    ),
    dnyears_contours = contour_table(
      models$dnyears,
      dnyears_grid,
      x_col = "years",
      y_col = "d_n",
      levels = thresholds$threshold[thresholds$model == "s(dn) + s(years)"],
      costs = thresholds$cost[thresholds$model == "s(dn) + s(years)"]
    )
  )
}

plot_classifier_diagnostics <- function(scored_training, thresholds, surfaces, config) {
  colors <- config$color_map
  cls_colors <- config$classifier_color_map

  dnds_roc <- roc(
    response = scored_training$classifier_class,
    predictor = scored_training$cVDPVness,
    levels = c("iVDPV", "cVDPV"),
    direction = "<",
    quiet = TRUE
  )
  dnyears_roc <- roc(
    response = scored_training$classifier_class,
    predictor = scored_training$cVDPVness_years,
    levels = c("iVDPV", "cVDPV"),
    direction = "<",
    quiet = TRUE
  )

  dnds_thresh <- thresholds |> filter(model == "s(dn) + s(ds)")
  dnyears_thresh <- thresholds |> filter(model == "s(dn) + s(years)")

  p_roc_dnds <- ggroc(dnds_roc) +
    geom_point(data = dnds_thresh, aes(x = specificity, y = sensitivity, shape = cost), size = 2.6, color = "#222222") +
    geom_text(data = dnds_thresh, aes(x = specificity, y = sensitivity, label = cost), nudge_x = -0.04, nudge_y = -0.06, size = 3, hjust = 0) +
    theme_bw() +
    annotate("text", x = 0.30, y = 0.20, label = paste0("AUC = ", sprintf("%.3f", as.numeric(auc(dnds_roc))))) +
    guides(shape = "none")

  p_ecdf_dnds <- ggplot() +
    stat_ecdf(data = scored_training, aes(x = cVDPVness, color = poliovirus_classification)) +
    geom_vline(data = dnds_thresh, aes(xintercept = threshold, linetype = cost)) +
    theme_bw() +
    scale_color_manual(values = colors, guide = "none") +
    guides(linetype = "none") +
    ylab("ECDF")

  p_logit_dnds <- ggplot(
    scored_training,
    aes(x = classifier_class, y = safe_logit(cVDPVness), fill = classifier_class, group = classifier_class)
  ) +
    geom_violin(scale = "width", alpha = 0.45) +
    geom_quasirandom(aes(color = poliovirus_classification), varwidth = TRUE, width = 0.45, alpha = 0.45) +
    geom_hline(data = dnds_thresh, aes(yintercept = safe_logit(threshold), linetype = cost)) +
    theme_bw() +
    scale_fill_manual(values = cls_colors, guide = "none") +
    scale_color_manual(values = colors, name = "classification") +
    guides(
      color = guide_legend(override.aes = list(alpha = 1, size = 2.2)),
      linetype = guide_legend(title = "cost")
    ) +
    ylab("logit(cVDPVness)") +
    xlab("")

  p_surface_dnds <- ggplot(scored_training, aes(x = d_s, y = d_n)) +
    geom_jitter(aes(color = poliovirus_classification), width = 0.3, height = 0.3, alpha = 0.25) +
    geom_path(data = surfaces$dnds_contours, aes(group = interaction(comp, cost), linetype = cost)) +
    theme_bw() +
    ylab("non-synonymous codon differences") +
    xlab("synonymous codon differences") +
    scale_x_continuous(breaks = seq(0, 300, by = 20)) +
    scale_color_manual(values = colors, guide = "none") +
    guides(linetype = "none")

  p_surface_years <- ggplot(scored_training, aes(x = years, y = d_n)) +
    geom_jitter(aes(color = poliovirus_classification), width = 0.2, height = 0.3, alpha = 0.25) +
    geom_path(data = surfaces$dnyears_contours, aes(group = interaction(comp, cost), linetype = cost)) +
    theme_bw() +
    ylab("non-synonymous codon differences") +
    xlab("approx. years since OPV") +
    scale_x_continuous(breaks = seq(0, 35, by = 5)) +
    scale_color_manual(values = colors, guide = "none") +
    guides(linetype = "none")

  ggsave(
    file.path(config$project_dir, config$output_dir, "classifier_dnds.png"),
    (p_roc_dnds + p_ecdf_dnds + p_logit_dnds) / (p_surface_dnds + p_surface_years) +
      plot_layout(guides = "collect"),
    units = "in", width = 13, height = 6.5
  )

  p_roc_years <- ggroc(dnyears_roc) +
    geom_point(data = dnyears_thresh, aes(x = specificity, y = sensitivity, shape = cost), size = 2.6, color = "#222222") +
    geom_text(data = dnyears_thresh, aes(x = specificity, y = sensitivity, label = cost), nudge_x = -0.04, nudge_y = -0.06, size = 3, hjust = 0) +
    theme_bw() +
    annotate("text", x = 0.30, y = 0.20, label = paste0("AUC = ", sprintf("%.3f", as.numeric(auc(dnyears_roc))))) +
    guides(shape = "none")

  p_ecdf_years <- ggplot() +
    stat_ecdf(data = scored_training, aes(x = cVDPVness_years, color = poliovirus_classification)) +
    geom_vline(data = dnyears_thresh, aes(xintercept = threshold, linetype = cost)) +
    theme_bw() +
    scale_color_manual(values = colors, guide = "none") +
    guides(linetype = "none") +
    ylab("ECDF")

  p_logit_years <- ggplot(
    scored_training,
    aes(x = classifier_class, y = safe_logit(cVDPVness_years), fill = classifier_class, group = classifier_class)
  ) +
    geom_violin(scale = "width", alpha = 0.45) +
    geom_quasirandom(aes(color = poliovirus_classification), varwidth = TRUE, width = 0.45, alpha = 0.45) +
    geom_hline(data = dnyears_thresh, aes(yintercept = safe_logit(threshold), linetype = cost)) +
    theme_bw() +
    scale_fill_manual(values = cls_colors, guide = "none") +
    scale_color_manual(values = colors, name = "classification") +
    guides(
      color = guide_legend(override.aes = list(alpha = 1, size = 2.2)),
      linetype = guide_legend(title = "cost")
    ) +
    ylab("logit(cVDPVness)") +
    xlab("")

  ggsave(
    file.path(config$project_dir, config$output_dir, "classifier_dnyears.png"),
    (p_roc_years + p_ecdf_years + p_logit_years) + plot_layout(guides = "collect"),
    units = "in", width = 13, height = 4.0
  )
}

plot_cvdpvness_projection <- function(meta, model, dnds_contours, config) {
  meta_pred <- meta |>
    filter(poliovirus_classification %in% config$projection_classes) |>
    mutate(poliovirus_classification = droplevels(poliovirus_classification))
  meta_pred$cVDPVness <- as.numeric(predict(model, newdata = meta_pred, type = "response"))

  p_class <- ggplot(meta_pred, aes(x = d_s, y = d_n, color = poliovirus_classification)) +
    geom_jitter(width = 0.3, height = 0.3, alpha = 0.25) +
    theme_bw() +
    ylab("non-synonymous codon differences") +
    xlab("synonymous codon differences") +
    scale_x_continuous(breaks = seq(0, 300, by = 20)) +
    scale_color_manual(values = config$color_map, name = "classification") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2.2)))

  p_score <- ggplot(meta_pred, aes(x = d_s, y = d_n)) +
    geom_jitter(aes(color = cVDPVness), width = 0.3, height = 0.3, alpha = 0.3) +
    geom_path(data = dnds_contours, aes(group = interaction(comp, cost), linetype = cost)) +
    theme_bw() +
    ylab("non-synonymous codon differences") +
    xlab("synonymous codon differences") +
    scale_x_continuous(breaks = seq(0, 300, by = 40)) +
    facet_wrap(~ poliovirus_classification, nrow = 2)

  ggsave(
    file.path(config$project_dir, config$output_dir, "cVDPVness.png"),
    p_class / p_score + plot_layout(heights = c(1, 1.6)),
    units = "in", width = 13, height = 7.5
  )

  meta_pred
}

relative_fitness_normalization <- function(meta, config) {
  rf <- config$relative_fitness
  cVDPV_mean <- mean(
    rf$base^(rf$genome_nt / rf$vp1_nt * abs(meta$d_n[meta$classifier_class == "cVDPV"] - rf$consensus_d_n)),
    na.rm = TRUE
  )
  c(rf, list(cVDPV_mean = cVDPV_mean))
}

serotype_stratified_thresholds <- function(scored_training) {
  specs <- tribble(
    ~model,                  ~pred_col,
    "pooled",                "cVDPVness",
    "serotype_factor",       "cVDPVness_factor",
    "independent_serotype",  "cVDPVness_independent"
  )

  eval_one <- function(df, pred_col, model_label, scope_label) {
    r <- roc(
      response = df$classifier_class,
      predictor = df[[pred_col]],
      levels = c("iVDPV", "cVDPV"),
      direction = "<",
      quiet = TRUE
    )
    pick_one <- function(x) as_tibble(x)[1L, , drop = FALSE]
    y <- pick_one(coords(r, "best", best.method = "youden"))
    w <- pick_one(coords(r, "best", best.method = "youden", best.weights = c(10, 0.5)))

    tibble(
      model = model_label,
      scope = scope_label,
      n_iVDPV = sum(df$classifier_class == "iVDPV"),
      n_cVDPV = sum(df$classifier_class == "cVDPV"),
      auc = as.numeric(auc(r)),
      cost = c("equal", "c/i = 10x"),
      threshold = c(y$threshold, w$threshold),
      specificity = c(y$specificity, w$specificity),
      sensitivity = c(y$sensitivity, w$sensitivity)
    )
  }

  scopes <- c("All", levels(scored_training$virus_type))
  map_dfr(scopes, function(sc) {
    sub <- if (sc == "All") scored_training else scored_training |> filter(virus_type == sc)
    pmap_dfr(specs, function(model, pred_col) eval_one(sub, pred_col, model, sc))
  })
}

build_serotype_stratified_surfaces <- function(models, stratified_thresholds, training) {
  max_ds <- max(140, ceiling(max(training$d_s, na.rm = TRUE) / 10) * 10 + 10)
  max_dn <- max(50, ceiling(max(training$d_n, na.rm = TRUE) / 5) * 5 + 5)
  serotypes <- levels(training$virus_type)

  # Evaluate contours for each serotype and each model under both per-serotype
  # and global ("All") ROC thresholds.
  map_dfr(serotypes, function(vt) {
    grid <- expand.grid(d_s = seq(0, max_ds, by = 1), d_n = seq(0, max_dn, by = 1)) |>
      mutate(virus_type = factor(vt, levels = serotypes))
    x_values <- sort(unique(grid$d_s))
    y_values <- sort(unique(grid$d_n))

    pred_list <- list(
      pooled = as.numeric(predict(models$dnds, newdata = grid, type = "response")),
      serotype_factor = as.numeric(predict(models$factor_dnds, newdata = grid, type = "response")),
      independent_serotype = as.numeric(predict(models$independent_dnds[[vt]], newdata = grid, type = "response"))
    )

    map_dfr(names(pred_list), function(m_name) {
      zmat <- matrix(pred_list[[m_name]], nrow = length(x_values), ncol = length(y_values))
      map_dfr(c("per_serotype", "global"), function(t_scope) {
        sc_key <- if (t_scope == "per_serotype") vt else "All"
        t_rows <- stratified_thresholds |> filter(model == m_name, scope == sc_key)
        comp_id <- 0L
        map2_dfr(t_rows$threshold, t_rows$cost, function(lvl, cst) {
          contours <- contourLines(x = x_values, y = y_values, z = zmat, levels = lvl)
          if (length(contours) == 0L) return(tibble())
          map_dfr(seq_along(contours), function(i) {
            comp_id <<- comp_id + 1L
            tibble(
              virus_type = factor(vt, levels = serotypes),
              model = m_name,
              threshold_scope = t_scope,
              comp = comp_id,
              d_s = contours[[i]]$x,
              d_n = contours[[i]]$y,
              cost = cst
            )
          })
        })
      })
    })
  })
}

plot_serotype_stratified_diagnostics <- function(scored_training, stratified_thresholds, stratified_contours, config) {
  colors <- config$color_map
  model_colors <- c(
    pooled = "#555555",
    serotype_factor = "#e66101",
    independent_serotype = "#5e3c99"
  )
  model_labels <- c(
    pooled = "Pooled: s(dn)+s(ds)",
    serotype_factor = "Factor: s(dn)+s(ds)+type",
    independent_serotype = "Independent per type"
  )

  # 1. Decision contours per serotype (using per-serotype calibrated thresholds)
  contours_per_st <- stratified_contours |>
    filter(threshold_scope == "per_serotype") |>
    mutate(model = factor(model, levels = names(model_colors)))

  p_surfaces <- ggplot(scored_training, aes(x = d_s, y = d_n)) +
    geom_jitter(aes(color = poliovirus_classification), width = 0.3, height = 0.3, alpha = 0.3) +
    geom_path(
      data = contours_per_st,
      aes(
        group = interaction(model, comp, cost),
        color = model,
        linetype = cost
      ),
      linewidth = 0.7
    ) +
    facet_grid(. ~ virus_type) +
    theme_bw() +
    ylab("non-synonymous codon differences") +
    xlab("synonymous codon differences") +
    coord_cartesian(xlim = c(0, 120), ylim = c(0, 32)) +
    scale_x_continuous(breaks = seq(0, 120, by = 20)) +
    scale_color_manual(
      values = c(colors, model_colors),
      breaks = c(names(model_colors), "iVDPV", "cVDPV", "cVDPV-n"),
      labels = c(unname(model_labels), "iVDPV", "cVDPV", "cVDPV-n"),
      name = "model / class"
    ) +
    guides(
      color = guide_legend(override.aes = list(alpha = 1)),
      linetype = guide_legend(title = "cost")
    )

  # 2. Per-serotype ROC curves comparing the 3 models
  serotypes <- levels(scored_training$virus_type)
  roc_curve_df <- map_dfr(serotypes, function(vt) {
    sub <- scored_training |> filter(virus_type == vt)
    specs <- list(
      pooled = sub$cVDPVness,
      serotype_factor = sub$cVDPVness_factor,
      independent_serotype = sub$cVDPVness_independent
    )
    map_dfr(names(specs), function(m_name) {
      r <- roc(sub$classifier_class, specs[[m_name]], levels = c("iVDPV", "cVDPV"), direction = "<", quiet = TRUE)
      tibble(
        virus_type = factor(vt, levels = serotypes),
        model = factor(m_name, levels = names(model_colors)),
        specificity = r$specificities,
        sensitivity = r$sensitivities
      ) |>
        arrange(desc(specificity), sensitivity)
    })
  })

  auc_labels <- stratified_thresholds |>
    filter(scope != "All", cost == "equal") |>
    mutate(
      virus_type = factor(scope, levels = serotypes),
      model = factor(model, levels = names(model_colors)),
      label = sprintf("%s: AUC = %.3f", model_labels[as.character(model)], auc),
      x = 0.03,
      y = case_when(
        model == "pooled" ~ 0.28,
        model == "serotype_factor" ~ 0.18,
        TRUE ~ 0.08
      )
    )

  p_rocs <- ggplot(roc_curve_df, aes(x = specificity, y = sensitivity, color = model)) +
    geom_path(linewidth = 0.8) +
    geom_text(
      data = auc_labels,
      aes(x = x, y = y, label = label, color = model),
      hjust = 1, size = 2.9, show.legend = FALSE
    ) +
    scale_x_reverse(limits = c(1, 0)) +
    facet_grid(. ~ virus_type) +
    theme_bw() +
    scale_color_manual(values = model_colors, guide = "none") +
    xlab("specificity") +
    ylab("sensitivity")

  # 3. Per-serotype logit(cVDPVness) distributions across the 3 models
  logit_long <- scored_training |>
    select(accession, virus_type, classifier_class, poliovirus_classification,
           pooled = cVDPVness, serotype_factor = cVDPVness_factor, independent_serotype = cVDPVness_independent) |>
    pivot_longer(
      cols = c(pooled, serotype_factor, independent_serotype),
      names_to = "model",
      values_to = "score"
    ) |>
    mutate(
      model = factor(model, levels = names(model_colors), labels = unname(model_labels)),
      logit_score = safe_logit(score)
    )

  p_logits <- ggplot(logit_long, aes(x = classifier_class, y = logit_score, fill = classifier_class)) +
    geom_violin(scale = "width", alpha = 0.4) +
    geom_quasirandom(aes(color = poliovirus_classification), varwidth = TRUE, width = 0.4, alpha = 0.4, size = 1.1) +
    facet_grid(model ~ virus_type, scales = "free_y") +
    theme_bw() +
    scale_fill_manual(values = config$classifier_color_map, guide = "none") +
    scale_color_manual(values = colors, guide = "none") +
    ylab("logit(cVDPVness)") +
    xlab("")

  ggsave(
    file.path(config$project_dir, config$output_dir, "classifier_serotype_dnds.png"),
    p_surfaces / p_rocs / p_logits + plot_layout(heights = c(1.1, 1.0, 1.8), guides = "collect"),
    units = "in", width = 13.5, height = 11
  )
}

write_classifier_outputs <- function(meta, training, thresholds, surfaces, models,
                                     stratified_thresholds, stratified_contours, config) {
  out_dir <- file.path(config$project_dir, config$output_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  norm <- relative_fitness_normalization(meta, config)
  source_manifest <- read_source_manifest(config)

  artifacts <- list(
    mod3 = models$dnds,
    mod5 = models$dnyears,
    mod_factor_dnds = models$factor_dnds,
    mod_independent_dnds = models$independent_dnds,
    training_scores = training,
    thresholds = thresholds,
    serotype_thresholds = stratified_thresholds,
    dnds_contours = surfaces$dnds_contours,
    dnyears_contours = surfaces$dnyears_contours,
    serotype_dnds_contours = stratified_contours,
    color_map = config$color_map,
    classifier_color_map = config$classifier_color_map,
    relative_fitness_norm = norm,
    model_comparison = models$model_comparison,
    serotype_model_comparison = models$serotype_model_comparison,
    source_manifest = source_manifest
  )

  write_rds(artifacts, file.path(out_dir, "classifier_artifacts.rds"))
  write_rds(meta, file.path(out_dir, "processed_reference_meta.rds"))
  write.csv(meta, file.path(out_dir, "processed_reference_meta.csv"), row.names = FALSE)
  write.csv(training, file.path(out_dir, "classifier_training_scores.csv"), row.names = FALSE)
  write.csv(thresholds, file.path(out_dir, "classifier_thresholds.csv"), row.names = FALSE)
  write.csv(stratified_thresholds, file.path(out_dir, "classifier_serotype_thresholds.csv"), row.names = FALSE)
  write.csv(models$serotype_model_comparison, file.path(out_dir, "classifier_serotype_model_comparison.csv"), row.names = FALSE)
  write.csv(surfaces$dnds_contours, file.path(out_dir, "classifier_dnds_contours.csv"), row.names = FALSE)
  write.csv(surfaces$dnyears_contours, file.path(out_dir, "classifier_dnyears_contours.csv"), row.names = FALSE)
  write.csv(stratified_contours, file.path(out_dir, "classifier_serotype_dnds_contours.csv"), row.names = FALSE)
  write.csv(
    data.frame(poliovirus_classification = names(config$color_map), color = unname(config$color_map)),
    file.path(out_dir, "classifier_color_map.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(parameter = names(norm), value = unlist(norm, use.names = FALSE)),
    file.path(out_dir, "relative_fitness_normalization.csv"),
    row.names = FALSE
  )

  artifacts
}

run_classifier_workflow <- function(config = classifier_config()) {
  setup_classifier_graphics()
  dir.create(file.path(config$project_dir, config$output_dir), recursive = TRUE, showWarnings = FALSE)

  meta <- prepare_reference_meta(config)
  training <- classifier_training_data(meta, config)

  plot_reference_overview(meta, training, config)
  plot_relative_fitness_model(training, config)
  plot_zhao_consistency(training, config)

  models <- fit_classifier_models(training)
  scored_training <- score_training_data(training, models)
  thresholds <- classifier_thresholds(scored_training)
  surfaces <- build_classifier_surfaces(models, thresholds, scored_training)

  stratified_thresholds <- serotype_stratified_thresholds(scored_training)
  stratified_contours <- build_serotype_stratified_surfaces(models, stratified_thresholds, scored_training)

  plot_classifier_diagnostics(scored_training, thresholds, surfaces, config)
  plot_cvdpvness_projection(meta, models$dnds, surfaces$dnds_contours, config)
  plot_serotype_stratified_diagnostics(scored_training, stratified_thresholds, stratified_contours, config)

  artifacts <- write_classifier_outputs(
    meta, scored_training, thresholds, surfaces, models,
    stratified_thresholds, stratified_contours, config
  )
  invisible(list(
    meta = meta,
    training = scored_training,
    models = models,
    thresholds = thresholds,
    serotype_thresholds = stratified_thresholds,
    surfaces = surfaces,
    serotype_contours = stratified_contours,
    artifacts = artifacts
  ))
}
