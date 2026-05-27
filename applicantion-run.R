RUN_APPLICATION_MAIN <- FALSE

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- sub(file_arg, "", cmd_args[grep(file_arg, cmd_args)])
  if (length(script_path) > 0 && nzchar(script_path[1])) {
    return(dirname(normalizePath(script_path[1], winslash = "/", mustWork = FALSE)))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile,
                                 winslash = "/", mustWork = FALSE)))
  }
  getwd()
}

script_dir <- get_script_dir()
revision_dir <- normalizePath(file.path(script_dir, ".."),
                              winslash = "/", mustWork = FALSE)
source(file.path(revision_dir, "code", "application-final.R"), encoding = "UTF-8")
RUN_APPLICATION_MAIN <- TRUE

######################## Complete application analysis ########################
# This caller runs the traditional strategies, JointMR, and relaxed NOME in one pass.

parse_env_list <- function(name, default) {
  x <- Sys.getenv(name, "")
  if (!nzchar(x)) {
    return(default)
  }
  out <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  out[nzchar(out)]
}

format_threshold_tag <- function(x) {
  gsub("\\.", "p", as.character(x))
}

log_main <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      " | complete-strictF-stouffer | ",
      paste0(..., collapse = ""), "\n", sep = "")
  flush.console()
}

safe_meta_summary <- function(meta_obj) {
  if (inherits(meta_obj, "try-error") || is.null(meta_obj$pval.Q)) {
    return(c(beta = NA_real_, se = NA_real_, pval = NA_real_))
  }
  if (is.finite(meta_obj$pval.Q) && meta_obj$pval.Q < 0.05) {
    c(beta = meta_obj$TE.random,
      se = meta_obj$seTE.random,
      pval = meta_obj$pval.random)
  } else {
    c(beta = meta_obj$TE.fixed,
      se = meta_obj$seTE.fixed,
      pval = meta_obj$pval.fixed)
  }
}

prefix_list_names <- function(x, prefix) {
  names(x) <- paste0(prefix, names(x))
  x
}

pack_jointmr_grid <- function(re) {
  list(
    SNPnew = re$SNPnew,
    original_beta = re$theta_original,
    original_se = re$theta_original_se,
    original_pval = re$theta_original_p_value,
    pleiotropy_beta = re$theta_hat,
    pleiotropy_se = re$theta_se,
    pleiotropy_pval = re$theta_p_value,
    relaxedNOME_beta = re$theta_relaxed_nome,
    relaxedNOME_se = re$theta_relaxed_nome_se,
    relaxedNOME_pval = re$theta_relaxed_nome_p_value,
    relaxedNOME_profile_se = re$theta_relaxed_nome_profile_se,
    relaxedNOME_se_method = re$relaxed_nome_se_method,
    relaxedNOME_bootstrap_success = re$relaxed_nome_bootstrap_success,
    relaxedNOME_bootstrap_requested = re$relaxed_nome_bootstrap_requested,
    relaxedNOME_pleiotropy_adjusted = re$relaxed_nome_pleiotropy_adjusted,
    relaxedNOME_tau2_used = re$relaxed_nome_tau2_used,
    relaxedNOME_xy_overlap_pair_count = if (!is.null(re$xy_overlap_matrix)) {
      sum(re$xy_overlap_matrix, na.rm = TRUE)
    } else {
      NA_real_
    },
    relaxedNOME_xy_overlap_pairs = format_xy_overlap_pairs(re$xy_overlap_matrix),
    relaxedNOME_xx_correlation_pairs = format_xx_correlation_pairs(re$exposure_cor_matrix),
    relaxedNOME_xy_overlap_correlation_pairs = format_xy_correlation_pairs(
      re$xy_rho_matrix,
      re$xy_overlap_matrix,
      re$xy_rho_available_matrix
    ),
    pleiotropy_global_Q = re$pleiotropy_Q,
    pleiotropy_global_df = re$pleiotropy_df,
    pleiotropy_global_pval = re$pleiotropy_p_value,
    pleiotropy_global_significant = as.integer(isTRUE(re$pleiotropy_significant)),
    heterogeneity_original_tau2_DL = re$tau2_original_DL,
    heterogeneity_original_I2_DL = re$I2_original_DL,
    heterogeneity_original_Q = re$original_DL_Q,
    heterogeneity_original_df = re$original_DL_df,
    heterogeneity_original_pval = re$original_DL_p_value,
    heterogeneity_original_significant = as.integer(is.finite(re$original_DL_p_value) &&
                                                      re$original_DL_p_value < 0.05),
    heterogeneity_pleiotropy_tau2_used = re$tau_sq,
    heterogeneity_pleiotropy_I2_DL = re$I2_DL,
    heterogeneity_pleiotropy_Q = re$DL_Q,
    heterogeneity_pleiotropy_df = re$DL_df,
    heterogeneity_pleiotropy_pval = re$DL_p_value,
    heterogeneity_pleiotropy_significant = as.integer(is.finite(re$DL_p_value) &&
                                                        re$DL_p_value < 0.05)
  )
}

name_pairwise_mr_results <- function(x, nx, ny) {
  methods <- c("dIVW", "IVW_fixed", "IVW_random", "wme")
  all_beta_names <- unlist(lapply(methods, function(method) {
    paste0("MR_", method, "_beta", seq_len(nx * ny))
  }))
  all_se_names <- unlist(lapply(methods, function(method) {
    paste0("MR_", method, "_se", seq_len(nx * ny))
  }))
  all_pval_names <- unlist(lapply(methods, function(method) {
    paste0("MR_", method, "_pval", seq_len(nx * ny))
  }))
  names(x) <- c(all_beta_names, all_se_names, all_pval_names)
  x
}

run_complete_combo <- function(exposure_list,
                               outcome_list,
                               exposure_name,
                               outcome_name,
                               base_path,
                               output_dir,
                               clump_kb,
                               clump_r2,
                               JointMR_F_threshold,
                               bootstrap_time,
                               relaxed_nome_bootstrap_time,
                               f_label = "strictF",
                               f_selection_value = "all",
                               p_combine_grid = "stouffer") {
  combo_name <- paste0(gsub("_list", "", exposure_name), "_",
                       gsub("_list", "", outcome_name))
  ex_file <- sub("\\d.*", "", exposure_name)
  nx <- length(exposure_list)
  ny <- length(outcome_list)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)

  log_main("[", combo_name, "] traditional MR-Meta started")
  results_exposure <- process_gwas_data(
    exposure_list,
    ncore = 8,
    pval_threshold = 5e-8,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    F_threshold = 10
  )
  SNP_list <- results_exposure$clumped_results
  results_traditional_mr_meta <- run_mr_analysis_no_ref(
    exposure_list = exposure_list,
    outcome_list = outcome_list,
    SNP_list = SNP_list,
    nx = nx,
    ny = ny
  )

  res_divw <- results_traditional_mr_meta$res_divw
  res_fixed <- results_traditional_mr_meta$res_fixed
  res_random <- results_traditional_mr_meta$res_random
  res_wme <- results_traditional_mr_meta$res_wme
  SNPSMR <- results_traditional_mr_meta$SNPSMR
  names(SNPSMR) <- paste0("SNPSMR_", seq_along(SNPSMR))

  pairwise_mr <- c(res_divw[, 1], res_fixed[, 1], res_random[, 1], res_wme[, 1],
                  res_divw[, 2], res_fixed[, 2], res_random[, 2], res_wme[, 2],
                  res_divw[, 3], res_fixed[, 3], res_random[, 3], res_wme[, 3])
  pairwise_mr <- name_pairwise_mr_results(pairwise_mr, nx, ny)

  meta_old_divw <- safe_meta_summary(
    try(metagen(res_divw[, 1], res_divw[, 2], sm = "MD"), silent = TRUE)
  )
  meta_old_ivw_fixed <- safe_meta_summary(
    try(metagen(res_fixed[, 1], res_fixed[, 2], sm = "MD"), silent = TRUE)
  )
  meta_old_ivw_random <- safe_meta_summary(
    try(metagen(res_random[, 1], res_random[, 2], sm = "MD"), silent = TRUE)
  )
  meta_old_wme <- safe_meta_summary(
    try(metagen(res_wme[, 1], res_wme[, 2], sm = "MD"), silent = TRUE)
  )
  mr_meta_old <- c(meta_old_divw,
                   meta_old_ivw_fixed,
                   meta_old_ivw_random,
                   meta_old_wme)
  names(mr_meta_old) <- c(
    paste0("MR_meta_old_dIVW_", c("beta", "se", "pval")),
    paste0("MR_meta_old_ivw_fixed_", c("beta", "se", "pval")),
    paste0("MR_meta_old_ivw_random_", c("beta", "se", "pval")),
    paste0("MR_meta_old_wme_", c("beta", "se", "pval"))
  )

  log_main("[", combo_name, "] GWAS-meta MR started")
  combo_dir <- file.path(base_path, combo_name)
  if (!dir.exists(combo_dir)) {
    dir.create(combo_dir, recursive = TRUE)
  }

  results_harmonise <- harmonize_and_export_efficient(
    data_list = exposure_list,
    output_dir = combo_dir,
    n_cores = 4
  )
  global_ref <- results_harmonise$global_ref

  exposure_metal_dir <- file.path(combo_dir, ex_file)
  dir.create(exposure_metal_dir, recursive = TRUE, showWarnings = FALSE)
  setwd(exposure_metal_dir)
  metal_script <- file.path(exposure_metal_dir, "metal_script.txt")
  dataset_names <- names(exposure_list)
  process_lines <- sapply(dataset_names, function(name) {
    file_path <- file.path(base_path, ex_file, paste0(name, ".txt"))
    paste0("PROCESS ", file_path)
  })
  script_lines <- trimws(c("SCHEME STDERR",
                           "MARKER SNP",
                           "ALLELE effect_allele other_allele",
                           "EFFECT beta",
                           "STDERR se",
                           "PVAL p",
                           process_lines,
                           "ANALYZE"))
  script_lines <- script_lines[nchar(script_lines) > 0]
  writeLines(script_lines, con = metal_script)
  system2(command = "/data/wusijia/generic-metal/metal",
          args = metal_script,
          stdout = TRUE,
          stderr = TRUE)

  alpha_data <- process_metal_results(
    file_path = file.path(exposure_metal_dir, "METAANALYSIS1.TBL"),
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    F_threshold = 10
  )
  a_clump_alpha_data1 <- alpha_data$clumped_data
  GWAS_Meta_exposure <- alpha_data$metal_data[
    alpha_data$metal_data$SNP %in% a_clump_alpha_data1$rsid, ]

  outcome_metal_dir <- file.path(combo_dir, "T2D")
  dir.create(outcome_metal_dir, recursive = TRUE, showWarnings = FALSE)
  harmonize_and_export_general(
    outcome_list = outcome_list,
    global_ref = global_ref,
    snp_list = a_clump_alpha_data1,
    output_dir = outcome_metal_dir,
    n_cores = 4
  )

  setwd(outcome_metal_dir)
  metal_script2 <- file.path(outcome_metal_dir, "metal_script2.txt")
  dataset_names2 <- names(outcome_list)
  process_lines2 <- sapply(seq_along(dataset_names2), function(i) {
    file_path <- file.path(outcome_metal_dir, paste0("dataset_", i, ".txt"))
    paste0("PROCESS ", file_path)
  })
  script_lines2 <- trimws(c("SCHEME STDERR",
                            "MARKER SNP",
                            "ALLELE effect_allele other_allele",
                            "EFFECT beta",
                            "STDERR se",
                            "PVAL p",
                            process_lines2,
                            "ANALYZE"))
  script_lines2 <- script_lines2[nchar(script_lines2) > 0]
  writeLines(script_lines2, con = metal_script2)
  system2(command = "/data/wusijia/generic-metal/metal",
          args = metal_script2,
          stdout = TRUE,
          stderr = TRUE)

  GWAS_Meta_outcome <- fread(file.path(outcome_metal_dir, "METAANALYSIS1.TBL"),
                             sep = "\t",
                             header = TRUE,
                             fill = TRUE)
  names(GWAS_Meta_outcome)[names(GWAS_Meta_outcome) == "P-value"] <- "pval"
  names(GWAS_Meta_outcome)[names(GWAS_Meta_outcome) == "MarkerName"] <- "SNP"
  GWAS_Meta_outcome <- subset(GWAS_Meta_outcome, SNP %in% GWAS_Meta_exposure$SNP)

  SNP_GWAS_meta <- Reduce(intersect, list(GWAS_Meta_exposure$SNP, GWAS_Meta_outcome$SNP))
  GWAS_Meta_exposure <- subset(GWAS_Meta_exposure, SNP %in% SNP_GWAS_meta)
  GWAS_Meta_outcome <- subset(GWAS_Meta_outcome, SNP %in% SNP_GWAS_meta)
  GWAS_Meta_exposure <- GWAS_Meta_exposure[order(GWAS_Meta_exposure$SNP), ]
  GWAS_Meta_outcome <- GWAS_Meta_outcome[order(GWAS_Meta_outcome$SNP), ]
  if (!identical(GWAS_Meta_exposure$SNP, GWAS_Meta_outcome$SNP)) {
    GWAS_Meta_outcome <- GWAS_Meta_outcome[
      match(GWAS_Meta_exposure$SNP, GWAS_Meta_outcome$SNP), ]
  }

  mr_object_fix <- mr_input(bx = GWAS_Meta_exposure$beta,
                            bxse = GWAS_Meta_exposure$StdErr,
                            by = GWAS_Meta_outcome$Effect,
                            byse = GWAS_Meta_outcome$StdErr)
  r_divw <- mr.divw(GWAS_Meta_exposure$beta,
                    GWAS_Meta_outcome$Effect,
                    GWAS_Meta_exposure$StdErr,
                    GWAS_Meta_outcome$StdErr)
  r_IVW_fixed <- MendelianRandomization::mr_ivw(mr_object_fix, model = "fixed")
  r_IVW_random <- MendelianRandomization::mr_ivw(mr_object_fix, model = "random")
  r_wme <- MendelianRandomization::mr_median(mr_object_fix, iterations = 20)
  gwas_meta_mr <- c(
    meta_res_dIVW_beta = r_divw$beta.hat,
    meta_res_dIVW_se = r_divw$beta.se,
    meta_res_dIVW_pval = ci(r_divw$beta.hat, r_divw$beta.se)$p,
    meta_res_IVW_fixed_beta = r_IVW_fixed@Estimate,
    meta_res_IVW_fixed_se = r_IVW_fixed@StdError,
    meta_res_IVW_fixed_pval = r_IVW_fixed@Pvalue,
    meta_res_IVW_random_beta = r_IVW_random@Estimate,
    meta_res_IVW_random_se = r_IVW_random@StdError,
    meta_res_IVW_random_pval = r_IVW_random@Pvalue,
    meta_res_wme_beta = r_wme@Estimate,
    meta_res_wme_se = r_wme@StdError,
    meta_res_wme_pval = r_wme@Pvalue
  )

  log_main("[", combo_name, "] JointMR and relaxed NOME started")
  N_list <- get_application_N_list(exposure_list, outcome_list)
  jointmr_values <- list()
  ldsc_cache <- NULL

  for (p_method in p_combine_grid) {
    processed_grid <- full_dataset_processing(
      exposure_list = exposure_list,
      outcome_list = outcome_list,
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      F_threshold = JointMR_F_threshold,
      f_selection = f_selection_value,
      p_combine_method = p_method
    )

    re_grid <- run_new_mr_analysis(
      exposure_filter_list = processed_grid$processed_exposure,
      outcome_filter_list = processed_grid$processed_outcome,
      original_exposure_list = exposure_list,
      original_outcome_list = outcome_list,
      N_list = N_list,
      ancestry = "EUR",
      bootstrap_time = bootstrap_time,
      relaxed_nome_bootstrap_time = relaxed_nome_bootstrap_time,
      ldsc_cache = ldsc_cache,
      run_relaxed_nome = TRUE,
      ldsc_scope = "all"
    )
    if (is.null(ldsc_cache)) {
      ldsc_cache <- re_grid$ldsc_cache
    }

    jointmr_values <- c(
      jointmr_values,
      prefix_list_names(pack_jointmr_grid(re_grid),
                        paste0("JointMR_", f_label, "_", p_method, "_"))
    )
  }

  exposure_dataset_values <- as.list(names(exposure_list))
  names(exposure_dataset_values) <- paste0("exposure_dataset_",
                                           seq_along(exposure_dataset_values))
  outcome_dataset_values <- as.list(names(outcome_list))
  names(outcome_dataset_values) <- paste0("outcome_dataset_",
                                          seq_along(outcome_dataset_values))

  res_values <- c(
    list(combo_name = combo_name,
         clump_kb = clump_kb,
         clump_r2 = clump_r2,
         traditional_F_threshold = 10,
         JointMR_F_threshold = JointMR_F_threshold,
         JointMR_F_selection = paste0("strictF_allX_F_gt_", JointMR_F_threshold),
         bootstrap_time_original_pleiotropy = bootstrap_time,
         bootstrap_time_relaxedNOME = relaxed_nome_bootstrap_time,
         relaxedNOME_se_method = "bootstrap_with_profile_hessian_fallback",
         nx = nx,
         ny = ny),
    as.list(SNPSMR),
    list(SNP_GWAS_meta = length(SNP_GWAS_meta)),
    exposure_dataset_values,
    outcome_dataset_values,
    as.list(gwas_meta_mr),
    as.list(pairwise_mr),
    as.list(mr_meta_old),
    jointmr_values
  )

  res_all <- as.data.frame(res_values, check.names = FALSE)
  out_file <- file.path(
    output_dir,
    paste0(combo_name,
           "_complete_strictF_stouffer_3scenario_F",
           format_threshold_tag(JointMR_F_threshold),
           "_0.001_boot",
           bootstrap_time,
           "_relaxedNOMEboot",
           relaxed_nome_bootstrap_time,
           ".csv")
  )
  write.csv(res_all, out_file, row.names = FALSE)
  log_main("[", combo_name, "] saved: ", out_file)
  invisible(res_all)
}

exposure_names <- c("HDL22_list",
                    "LDL2_list",
                    "TC2_list",
                    "TG2_list")
outcome_names <- c("T2D5_list")

exposure_names <- parse_env_list("JOINTMR_EXPOSURES", exposure_names)
exposure_names <- ifelse(grepl("_list$", exposure_names), exposure_names,
                         paste0(exposure_names, "_list"))
outcome_names <- parse_env_list("JOINTMR_OUTCOMES", outcome_names)
outcome_names <- ifelse(grepl("_list$", outcome_names), outcome_names,
                        paste0(outcome_names, "_list"))

base_path <- Sys.getenv("JOINTMR_BASE_PATH", "/data/wusijia/")
output_dir <- file.path(base_path, "0506")
clump_kb <- 1000
clump_r2 <- 0.001
JointMR_F_threshold <- as.numeric(Sys.getenv("JOINTMR_F_THRESHOLD", "10"))
if (!is.finite(JointMR_F_threshold) || JointMR_F_threshold <= 0) {
  JointMR_F_threshold <- 10
}
bootstrap_time <- as.integer(Sys.getenv("JOINTMR_BOOTSTRAP", "500"))
if (!is.finite(bootstrap_time) || bootstrap_time < 0) {
  bootstrap_time <- 500L
}
relaxed_nome_bootstrap_time <- as.integer(Sys.getenv("JOINTMR_RELAXED_NOME_BOOTSTRAP", "200"))
if (!is.finite(relaxed_nome_bootstrap_time) || relaxed_nome_bootstrap_time < 0) {
  relaxed_nome_bootstrap_time <- 200L
}

f_label <- "strictF"
f_selection_value <- "all"
p_combine_grid <- "stouffer"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

log_main("script started; exposures=", paste(gsub("_list", "", exposure_names), collapse = ","),
         "; outcomes=", paste(gsub("_list", "", outcome_names), collapse = ","),
         "; traditional_F_threshold=10",
         "; JointMR_F_threshold=", JointMR_F_threshold,
         "; bootstrap_time_jointmr=", bootstrap_time,
         "; relaxedNOME_bootstrap_time=", relaxed_nome_bootstrap_time,
         "; p_methods=", paste(p_combine_grid, collapse = ","))

for (outcome_name in outcome_names) {
  log_main("loading outcome: ", outcome_name)
  load(paste0(base_path, "T2D/", outcome_name, ".Rdata"))

  for (exposure_name in exposure_names) {
    ex_file <- sub("\\d.*", "", exposure_name)
    log_main("loading exposure: ", exposure_name, " from ", ex_file)
    load(paste0(base_path, ex_file, "/", exposure_name, ".Rdata"))

    run_complete_combo(
      exposure_list = exposure_list,
      outcome_list = outcome_list,
      exposure_name = exposure_name,
      outcome_name = outcome_name,
      base_path = base_path,
      output_dir = output_dir,
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      JointMR_F_threshold = JointMR_F_threshold,
      bootstrap_time = bootstrap_time,
      relaxed_nome_bootstrap_time = relaxed_nome_bootstrap_time,
      f_label = f_label,
      f_selection_value = f_selection_value,
      p_combine_grid = p_combine_grid
    )

    rm(exposure_list)
    gc()
  }
  rm(outcome_list)
  gc()
}
