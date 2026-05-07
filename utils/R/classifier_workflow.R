suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
  library(phangorn)
  library(lubridate)
  library(Biostrings)
  library(patchwork)
  library(mgcv)
  library(pROC)
  library(ggbeeswarm)
})

classifier_config <- function(project_dir = getwd(),
                              data_dir = file.path("data", "genbank"),
                              output_dir = file.path("analyses", "genbank-classifier")) {
  list(
    project_dir = project_dir,
    data_dir = data_dir,
    output_dir = output_dir,
    serotypes = tibble(
      serotype = 1:3,
      vp1_nt = c(906L, 903L, 900L),
      p1_end = c(2643L, 2637L, 2634L)
    ),
    retained_classes = c("cVDPV", "iVDPV", "wild", "SabinLike", "VDPV", "Sabin"),
    distance_classes = c("cVDPV", "iVDPV", "wild", "SabinLike", "VDPV"),
    classifier_classes = c("cVDPV", "iVDPV"),
    color_map = c(
      SabinLike = "#04e762",
      VDPV = "#ae7766",
      iVDPV = "#00a1e4",
      cVDPV = "#dc0073",
      wild = "#f5b700"
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

vp1_slice <- function(alignment, p1_end, vp1_nt) {
  alignment[, p1_end + seq(-vp1_nt + 1L, 0L, by = 1L)]
}

drop_indel_rows <- function(alignment) {
  chars <- as.character(alignment)
  alignment[!apply(chars == "-", 1L, any), ]
}

read_serotype_inputs <- function(spec, config) {
  k <- spec$serotype
  project_dir <- config$project_dir
  alignment_path <- file.path(project_dir, config$data_dir, paste0("Sabin", k, ".ntP1.alignment.final.phylip.reduced"))
  metadata_path <- file.path(project_dir, config$data_dir, paste0("Sabin", k, ".ntP1.alignment.final.attributes.csv"))

  alignment <- read.dna(alignment_path, format = "sequential") |>
    vp1_slice(spec$p1_end, spec$vp1_nt) |>
    drop_indel_rows()

  metadata <- read.csv(metadata_path) |>
    mutate(classification = sub("Sabin[1-3]", "Sabin", classification)) |>
    filter(accession %in% labels(alignment), classification %in% config$retained_classes)

  alignment <- alignment[labels(alignment) %in% metadata$accession, ]
  aa_alignment <- as.AAbin(translate(
    DNAStringSet(apply(as.character(alignment), 1L, paste0, collapse = "")),
    if.fuzzy.codon = "X"
  ))

  list(metadata = metadata, alignment = alignment, aa_alignment = aa_alignment)
}

class_distance_table <- function(target_class, metadata, alignment, aa_alignment, serotype, vp1_nt) {
  keep <- metadata |> filter(classification %in% c(target_class, "Sabin"))
  keep_accessions <- keep$accession

  nt_distance <- as.matrix(dist.ml(alignment[labels(alignment) %in% keep_accessions, ], model = "F81"))[1L, ] * vp1_nt
  aa_distance <- as.matrix(dist.ml(aa_alignment[labels(aa_alignment) %in% keep_accessions], model = "FLU"))[1L, ] * vp1_nt / 3
  synonymous_distance <- nt_distance - aa_distance
  synonymous_distance <- synonymous_distance - min(synonymous_distance)

  years_since_opv <- if (target_class == "iVDPV") {
    # iVDPV elapsed time is closer to a nucleotide clock because patient-level
    # evolution includes substantial nonsynonymous drift.
    nt_distance / (vp1_nt * 0.0072)
  } else {
    # For cVDPV and wild-like groups, synonymous divergence gives the cleaner
    # clock-like proxy used by the original classifier.
    synonymous_distance / (vp1_nt * 0.01)
  }

  out <- tibble(
    accession = names(years_since_opv),
    classification = c("Sabin", rep(target_class, length(years_since_opv) - 1L)),
    serotype = serotype,
    years = as.numeric(years_since_opv),
    nt_diff = as.numeric(nt_distance),
    d_n = as.numeric(aa_distance),
    d_s = as.numeric(synonymous_distance)
  )

  if (target_class == "cVDPV") out else out |> filter(classification != "Sabin")
}

serotype_reference_meta <- function(spec, config) {
  inputs <- read_serotype_inputs(spec, config)
  distances <- map_dfr(
    config$distance_classes,
    class_distance_table,
    metadata = inputs$metadata,
    alignment = inputs$alignment,
    aa_alignment = inputs$aa_alignment,
    serotype = spec$serotype,
    vp1_nt = spec$vp1_nt
  )

  inputs$metadata |>
    left_join(distances, by = c("accession", "classification")) |>
    clean_reference_classifications()
}

clean_reference_classifications <- function(metadata) {
  metadata <- metadata |>
    mutate(classification = if_else(nt_diff == 0 & !grepl("Sabin", classification), "SabinLike", classification))

  # Three historical metadata corrections from the original analysis. They are
  # intentionally narrow accession-independent rules based on dating/distance.
  metadata$classification[metadata$classification == "cVDPV" & metadata$years > 12] <- "wild"
  metadata$classification[metadata$classification == "wild" & metadata$years < 12] <- "cVDPV"
  metadata$classification[metadata$classification == "VDPV" & metadata$d_s > 200] <- "wild"
  metadata
}

prepare_reference_meta <- function(config) {
  meta <- pmap_dfr(config$serotypes, \(serotype, vp1_nt, p1_end) {
    serotype_reference_meta(tibble(serotype = serotype, vp1_nt = vp1_nt, p1_end = p1_end), config)
  })

  classification_levels <- names(config$color_map)
  rf <- config$relative_fitness

  meta |>
    filter(classification != "Sabin") |>
    mutate(
      serotype = factor(serotype, levels = c("1", "2", "3")),
      classification = factor(classification, levels = classification_levels),
      d_n_from_consensus = abs(d_n - rf$consensus_d_n),
      relative_fitness = rf$base^(rf$genome_nt / rf$vp1_nt * d_n_from_consensus)
    ) |>
    mutate(relative_fitness = relative_fitness / mean(relative_fitness[classification == "cVDPV"]))
}

classifier_training_data <- function(meta, config) {
  meta |>
    filter(classification %in% config$classifier_classes) |>
    ungroup() |>
    droplevels()
}

plot_reference_overview <- function(meta, training, config) {
  colors <- config$color_map

  p_all <- ggplot(meta, aes(x = d_s, y = d_n, color = classification, shape = serotype)) +
    geom_point(alpha = 0.2) +
    ylab("non-synonymous substitutions") +
    xlab("synonymous substitutions") +
    scale_color_manual(values = colors)

  p_years <- ggplot(training, aes(x = years, y = d_n, color = classification, fill = classification, group = classification)) +
    geom_point(alpha = 0.2) +
    geom_smooth(method = "gam", fill = "black", alpha = 0.3) +
    ylab("non-synonymous substitutions") +
    xlab("approx. years since OPV") +
    scale_x_continuous(breaks = seq(0, 30, by = 5)) +
    scale_color_manual(values = colors)

  p_dnds <- ggplot(training, aes(x = d_s, y = d_n, color = classification, fill = classification)) +
    geom_point(aes(shape = serotype), alpha = 0.2) +
    geom_smooth(method = "gam", fill = "black", alpha = 0.3) +
    ylab("non-synonymous substitutions") +
    xlab("synonymous substitutions") +
    scale_x_continuous(breaks = seq(0, 300, by = 20)) +
    scale_color_manual(values = colors)

  ggsave(
    file.path(config$project_dir, config$output_dir, "iVDPV_vs_cVDPV.png"),
    p_years + p_dnds + p_all,
    units = "in", width = 12, height = 3.5
  )
}

plot_relative_fitness_model <- function(training, config) {
  sabin_df <- tibble(serotype = c("1", "2", "3"), relative_fitness = 2.3 / c(14, 8, 18))

  p <- ggplot(training, aes(x = years, y = relative_fitness, color = classification, fill = classification)) +
    geom_point(aes(shape = serotype), alpha = 0.2) +
    geom_smooth(method = "gam", fill = "black", alpha = 0.3) +
    geom_hline(data = sabin_df, aes(yintercept = relative_fitness, linetype = serotype)) +
    ylab("relative fitness") +
    xlab("approx. years since OPV") +
    scale_y_continuous(breaks = seq(0, 1.2, by = 0.2)) +
    scale_color_manual(values = config$color_map)

  ggsave(file.path(config$project_dir, config$output_dir, "relative_fitness_model.png"), p, units = "in", width = 5, height = 3.5)
}

plot_zhao_consistency <- function(training, config) {
  zhao <- read.csv(file.path(config$project_dir, config$data_dir, "zhao_fig2_extracted_vpdv_points.csv"))
  names(zhao) <- c("nt_diff", "d_n", "classification")
  zhao <- zhao |>
    mutate(
      nt_diff = (1 - nt_diff / 100) * 903,
      d_n = (1 - d_n / 100) * 301,
      d_s = nt_diff - d_n,
      d_s = d_s - min(d_s)
    )

  colors <- config$color_map
  p <- ggplot() +
    geom_point(data = training, aes(x = 100 * (1 - nt_diff / 903), y = 1 - d_n / 301, color = colors[["iVDPV"]])) +
    geom_point(data = zhao, aes(x = 100 * (1 - nt_diff / 903), y = 1 - d_n / 301, color = colors[["cVDPV"]])) +
    facet_grid(. ~ classification) +
    xlab("nucleotide identity (%)") +
    ylab("amino acid identity (%)") +
    scale_color_identity(
      name = "source",
      breaks = c(colors[["iVDPV"]], colors[["cVDPV"]]),
      labels = c("Distance from Sabin", "Zhao 2017\ninter-isolate distance"),
      guide = "legend"
    )

  # Keep the original misspelled filename for compatibility with existing notes.
  ggsave(file.path(config$project_dir, config$output_dir, "supp_Zhoa2017_consistency.png"), p, units = "in", width = 8, height = 3.5)
}

fit_classifier_models <- function(training) {
  mod1 <- gam(classification ~ d_n + d_s, data = training, family = "binomial")
  mod2 <- gam(classification ~ d_n + d_s + serotype, data = training, family = "binomial")
  mod3 <- gam(classification ~ s(d_n) + s(d_s), data = training, family = "binomial")
  mod4 <- gam(classification ~ ti(d_n, d_s), data = training, family = "binomial")
  mod5 <- gam(classification ~ s(d_n) + s(years), data = training, family = "binomial")

  list(
    linear = mod1,
    linear_serotype = mod2,
    dnds = mod3,
    interaction = mod4,
    dnyears = mod5,
    model_comparison = tibble(
      model = c("linear d_n+d_s", "linear + serotype", "smooth d_n+d_s", "tensor interaction", "smooth d_n+years"),
      aic = c(AIC(mod1), AIC(mod2), AIC(mod3), AIC(mod4), AIC(mod5))
    )
  )
}

score_training_data <- function(training, models) {
  training |>
    mutate(
      cVDPVness = predict(models$dnds, newdata = training, type = "response"),
      cVDPVness_years = predict(models$dnyears, newdata = training, type = "response")
    )
}

classifier_thresholds <- function(scored_training) {
  dnds_roc <- roc(scored_training$classification, scored_training$cVDPVness, quiet = TRUE)
  dnyears_roc <- roc(scored_training$classification, scored_training$cVDPVness_years, quiet = TRUE)

  youden <- coords(dnds_roc, "best", best.method = "youden")
  weighted <- coords(dnds_roc, "best", best.method = "youden", best.weights = c(10, 0.5))
  youden_years <- coords(dnyears_roc, "best", best.method = "youden")
  weighted_years <- coords(dnyears_roc, "best", best.method = "youden", best.weights = c(10, 0.5))

  tibble(
    specificity = c(youden$specificity, weighted$specificity, youden_years$specificity, weighted_years$specificity),
    sensitivity = c(youden$sensitivity, weighted$sensitivity, youden_years$sensitivity, weighted_years$sensitivity),
    threshold = c(youden$threshold, weighted$threshold, youden_years$threshold, weighted_years$threshold),
    cost = c("equal", "c/i = 10x", "equal", "c/i = 10x"),
    model = c("s(dn) + s(ds)", "s(dn) + s(ds)", "s(dn) + s(years)", "s(dn) + s(years)")
  )
}

contour_table <- function(model, grid, x_col, y_col, levels, costs) {
  x_values <- sort(unique(grid[[x_col]]))
  y_values <- sort(unique(grid[[y_col]]))
  zmat <- matrix(predict(model, newdata = grid, type = "response"), nrow = length(x_values), ncol = length(y_values))
  contours <- contourLines(x = x_values, y = y_values, z = zmat, levels = levels)

  map_dfr(seq_along(contours), function(i) {
    tibble(
      comp = i,
      !!x_col := contours[[i]]$x,
      !!y_col := contours[[i]]$y,
      cost = costs[match(contours[[i]]$level, levels)]
    )
  })
}

build_classifier_surfaces <- function(models, thresholds) {
  dnds_grid <- expand.grid(d_s = seq(0, 140, by = 1), d_n = seq(0, 50, by = 1))
  dnyears_grid <- expand.grid(years = seq(0, 30, by = 0.1), d_n = seq(0, 50, by = 1))

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
  dnds_roc <- roc(scored_training$classification, scored_training$cVDPVness, quiet = TRUE)
  dnyears_roc <- roc(scored_training$classification, scored_training$cVDPVness_years, quiet = TRUE)

  p_roc_dnds <- ggroc(dnds_roc) +
    geom_point(data = thresholds |> filter(model == "s(dn) + s(ds)"), aes(x = specificity, y = sensitivity, color = cost)) +
    theme_bw() +
    annotate("text", x = 0.2, y = 0.9, label = paste0("auc=", round(auc(dnds_roc), 2)))

  p_ecdf_dnds <- ggplot() +
    stat_ecdf(data = scored_training, aes(x = cVDPVness, color = classification)) +
    geom_hline(data = thresholds |> filter(model == "s(dn) + s(ds)"), aes(yintercept = threshold, linetype = cost)) +
    theme_bw() +
    guides() +
    scale_color_manual(values = colors)

  p_logit_dnds <- ggplot(scored_training, aes(x = classification, y = log(cVDPVness / (1 - cVDPVness)), fill = classification, group = classification)) +
    geom_violin(scale = "width") +
    geom_quasirandom(varwidth = TRUE, width = 0.45, alpha = 0.3) +
    geom_hline(data = thresholds |> filter(model == "s(dn) + s(ds)"),
               aes(yintercept = log(threshold / (1 - threshold)), linetype = cost)) +
    theme_bw() +
    scale_fill_manual(values = colors) +
    guides() +
    ylab("logit(cVDPVness)") +
    xlab("")

  p_surface_dnds <- ggplot(scored_training, aes(x = d_s, y = d_n)) +
    geom_jitter(aes(color = classification, fill = classification), width = 0.3, height = 0.3, alpha = 0.2) +
    geom_line(data = surfaces$dnds_contours, aes(linetype = cost)) +
    theme_bw() +
    ylab("non-synonymous substitutions") +
    xlab("synonymous substitutions") +
    scale_x_continuous(breaks = seq(0, 300, by = 20)) +
    scale_color_manual(values = colors)

  p_surface_years <- ggplot(scored_training, aes(x = years, y = d_n)) +
    geom_jitter(aes(color = classification, fill = classification), width = 0.3, height = 0.3, alpha = 0.2) +
    geom_line(data = surfaces$dnyears_contours, aes(linetype = cost)) +
    theme_bw() +
    ylab("non-synonymous substitutions") +
    xlab("approx. years since OPV") +
    scale_x_continuous(breaks = seq(0, 30, by = 5)) +
    scale_color_manual(values = colors)

  ggsave(
    file.path(config$project_dir, config$output_dir, "classifier_dnds.png"),
    (p_roc_dnds + p_ecdf_dnds + p_logit_dnds) / (p_surface_dnds + p_surface_years),
    units = "in", width = 12, height = 6
  )

  p_roc_years <- ggroc(dnyears_roc) +
    geom_point(data = thresholds |> filter(model != "s(dn) + s(ds)"), aes(x = specificity, y = sensitivity, color = cost)) +
    theme_bw() +
    annotate("text", x = 0.2, y = 0.9, label = paste0("auc=", round(auc(dnyears_roc), 2)))

  p_ecdf_years <- ggplot() +
    stat_ecdf(data = scored_training, aes(x = cVDPVness_years, color = classification)) +
    geom_hline(data = thresholds |> filter(model != "s(dn) + s(ds)"), aes(yintercept = threshold, linetype = cost)) +
    theme_bw() +
    guides() +
    scale_color_manual(values = colors)

  p_logit_years <- ggplot(scored_training, aes(x = classification, y = log(cVDPVness_years / (1 - cVDPVness_years)), fill = classification, group = classification)) +
    geom_violin(scale = "width") +
    geom_quasirandom(varwidth = TRUE, width = 0.45, alpha = 0.3) +
    geom_hline(data = thresholds |> filter(model != "s(dn) + s(ds)"),
               aes(yintercept = log(threshold / (1 - threshold)), linetype = cost)) +
    theme_bw() +
    guides() +
    ylab("logit(cVDPVness)") +
    xlab("") +
    scale_fill_manual(values = colors)

  ggsave(
    file.path(config$project_dir, config$output_dir, "classifier_dnyears.png"),
    p_roc_years + p_ecdf_years + p_logit_years,
    units = "in", width = 12, height = 3.5
  )
}

plot_cvdpvness_projection <- function(meta, model, dnds_contours, config) {
  meta_pred <- meta |>
    filter(classification %in% c("cVDPV", "iVDPV", "VDPV", "SabinLike"))
  meta_pred$cVDPVness <- predict(model, newdata = meta_pred, type = "response")

  p_class <- ggplot(meta_pred, aes(x = d_s, y = d_n, color = classification, fill = classification)) +
    geom_jitter(width = 0.3, height = 0.3, alpha = 0.2) +
    theme_bw() +
    ylab("non-synonymous substitutions") +
    xlab("synonymous substitutions") +
    scale_x_continuous(breaks = seq(0, 300, by = 20)) +
    scale_color_manual(values = config$color_map)

  p_score <- ggplot(meta_pred, aes(x = d_s, y = d_n)) +
    geom_jitter(aes(color = cVDPVness, fill = cVDPVness), width = 0.3, height = 0.3, alpha = 0.2) +
    geom_line(data = dnds_contours, aes(linetype = cost)) +
    theme_bw() +
    ylab("non-synonymous substitutions") +
    xlab("synonymous substitutions") +
    scale_x_continuous(breaks = seq(0, 300, by = 20)) +
    guides(shape = guide_legend("serotype")) +
    facet_grid(. ~ classification)

  ggsave(
    file.path(config$project_dir, config$output_dir, "cVDPVness.png"),
    p_class / p_score + plot_layout(guides = "collect"),
    units = "in", width = 12, height = 6
  )

  meta_pred
}

relative_fitness_normalization <- function(meta, config) {
  rf <- config$relative_fitness
  cVDPV_mean <- mean(rf$base^(rf$genome_nt / rf$vp1_nt * abs(meta$d_n[meta$classification == "cVDPV"] - rf$consensus_d_n)))
  c(rf, list(cVDPV_mean = cVDPV_mean))
}

write_classifier_outputs <- function(meta, training, thresholds, surfaces, models, config) {
  out_dir <- file.path(config$project_dir, config$output_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  norm <- relative_fitness_normalization(meta, config)
  artifacts <- list(
    mod3 = models$dnds,
    mod5 = models$dnyears,
    training_scores = training,
    thresholds = thresholds,
    dnds_contours = surfaces$dnds_contours,
    dnyears_contours = surfaces$dnyears_contours,
    color_map = config$color_map,
    relative_fitness_norm = norm
  )

  write_rds(artifacts, file.path(out_dir, "classifier_artifacts.rds"))
  write_rds(meta, file.path(out_dir, "processed_reference_meta.rds"))
  write.csv(meta, file.path(out_dir, "processed_reference_meta.csv"), row.names = FALSE)
  write.csv(training, file.path(out_dir, "classifier_training_scores.csv"), row.names = FALSE)
  write.csv(thresholds, file.path(out_dir, "classifier_thresholds.csv"), row.names = FALSE)
  write.csv(surfaces$dnds_contours, file.path(out_dir, "classifier_dnds_contours.csv"), row.names = FALSE)
  write.csv(surfaces$dnyears_contours, file.path(out_dir, "classifier_dnyears_contours.csv"), row.names = FALSE)
  write.csv(
    data.frame(classification = names(config$color_map), color = unname(config$color_map)),
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
  surfaces <- build_classifier_surfaces(models, thresholds)

  plot_classifier_diagnostics(scored_training, thresholds, surfaces, config)
  plot_cvdpvness_projection(meta, models$dnds, surfaces$dnds_contours, config)

  artifacts <- write_classifier_outputs(meta, scored_training, thresholds, surfaces, models, config)
  invisible(list(meta = meta, training = scored_training, models = models, thresholds = thresholds, surfaces = surfaces, artifacts = artifacts))
}
