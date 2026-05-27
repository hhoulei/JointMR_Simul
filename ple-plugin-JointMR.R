if (dir.exists("/home/houl/MR-META/20260421")) {
  setwd("/home/houl/MR-META/20260421")
}

library(MASS)
library(readr)
library(meta)
library(MendelianRandomization)
library(mr.divw)
library(foreach)
if (requireNamespace("doMC", quietly = TRUE)) {
  library(doMC)
}

ginv_solve <- function(Sigma, rhs, tol = sqrt(.Machine$double.eps)) {
  Sigma <- as.matrix(Sigma)
  Sigma <- (Sigma + t(Sigma)) / 2
  rhs <- as.matrix(rhs)
  MASS::ginv(Sigma, tol = tol) %*% rhs
}

solve_shared_re <- function(Omega, rhs, tau2 = 0, tol = sqrt(.Machine$double.eps)) {
  rhs <- as.matrix(rhs)
  if (tau2 <= 0) {
    return(ginv_solve(Omega, rhs, tol = tol))
  }
  
  n <- nrow(Omega)
  z <- matrix(1, nrow = n, ncol = 1)
  Omega_inv_rhs <- ginv_solve(Omega, rhs, tol = tol)
  Omega_inv_z <- ginv_solve(Omega, z, tol = tol)
  middle <- as.numeric(1 / tau2 + crossprod(z, Omega_inv_z))
  
  Omega_inv_rhs - Omega_inv_z %*%
    (crossprod(z, Omega_inv_rhs) / middle)
}

gls_theta_plugin <- function(base_cov_list, WR_matrix, BA_list, V_alpha,
                             tau2 = 0) {
  n_ratio <- nrow(WR_matrix)
  n_snp <- ncol(WR_matrix)
  n_total <- n_ratio * n_snp
  one <- matrix(1, nrow = n_total, ncol = 1)
  y <- matrix(as.vector(WR_matrix), nrow = n_total, ncol = 1)
  U <- do.call(rbind, BA_list)
  
  D_inv_apply <- function(X) {
    X <- as.matrix(X)
    out <- matrix(NA_real_, nrow = nrow(X), ncol = ncol(X))
    for (j in seq_len(n_snp)) {
      row_ids <- ((j - 1) * n_ratio + 1):(j * n_ratio)
      out[row_ids, ] <- solve_shared_re(
        base_cov_list[[j]],
        X[row_ids, , drop = FALSE],
        tau2 = tau2
      )
    }
    out
  }
  
  D_inv_one <- D_inv_apply(one)
  D_inv_y <- D_inv_apply(y)
  D_inv_U <- D_inv_apply(U)
  
  middle <- MASS::ginv(V_alpha) + t(U) %*% D_inv_U
  middle_inv <- MASS::ginv((middle + t(middle)) / 2)
  
  Sigma_inv_one <- D_inv_one - D_inv_U %*% middle_inv %*% t(U) %*% D_inv_one
  Sigma_inv_y <- D_inv_y - D_inv_U %*% middle_inv %*% t(U) %*% D_inv_y
  
  numerator <- as.numeric(crossprod(one, Sigma_inv_y))
  denominator <- as.numeric(crossprod(one, Sigma_inv_one))
  theta_hat <- numerator / denominator
  theta_se <- sqrt(1 / denominator)
  theta_p_value <- 2 * pnorm(-abs(theta_hat / theta_se))
  
  c(theta_hat = theta_hat, theta_se = theta_se, theta_p_value = theta_p_value)
}

make_common_exposure <- function(g, nx, bx, weakp) {
  N1 <- 5000
  N2 <- 50000
  maf <- 0.2
  var_g <- 2 * maf * (1 - maf)
  sigma_x <- 0.1
  
  weak_count <- max(1, round(g * weakp))
  strong_count <- g - weak_count
  strong_instruments_base <- runif(
    strong_count,
    min = max(0.02, bx - 0.1),
    max = bx + 0.1
  )
  weak_instruments_base <- runif(weak_count, 0.01, 0.03)
  shared_instrument_order <- sample(c(strong_instruments_base, weak_instruments_base))
  
  betaXG_initial <- lapply(seq_len(nx), function(i) {
    shared_instrument_order + rnorm(g, 0, 0.001)
  })
  
  se_x <- sigma_x / sqrt(N1 * var_g)
  se_x_p <- sigma_x / sqrt(N2 * var_g)
  betaXG <- lapply(seq_len(nx), function(i) {
    betaXG_initial[[i]] + rnorm(g, 0, se_x)
  })
  seXG <- lapply(seq_len(nx), function(i) rep(se_x, g))
  PXG <- lapply(seq_len(nx), function(i) {
    beta_x_p <- betaXG_initial[[i]] + rnorm(g, 0, se_x_p)
    2 * pnorm(-abs(beta_x_p / se_x_p))
  })
  
  list(
    gamma_true = shared_instrument_order,
    betaXG = betaXG,
    seXG = seXG,
    PXG = PXG
  )
}

make_cor_matrix <- function(nx, ny, rho_1, rho_2) {
  cor_matrix <- diag(1, nrow = nx * ny, ncol = nx * ny)
  
  for (m in seq_len(ny)) {
    for (q in seq_len(ny)) {
      row_ids <- ((m - 1) * nx + 1):(m * nx)
      col_ids <- ((q - 1) * nx + 1):(q * nx)
      if (m == q) {
        cor_matrix[row_ids, col_ids] <- matrix(1, nx, nx)
      } else {
        rho_c <- runif(1, min = rho_1, max = rho_2)
        cor_matrix[row_ids, col_ids] <- rho_c
        cor_matrix[col_ids, row_ids] <- rho_c
      }
    }
  }
  
  cor_matrix
}

DataGeneration_summary <- function(g, nx, ny, bx, beta0, tau, rho_1, rho_2,
                                   siga1, siga2, weakp) {
  exposure <- make_common_exposure(g, nx, bx, weakp)
  betaXG <- exposure$betaXG
  seXG <- exposure$seXG
  seYG <- lapply(seq_len(ny), function(i) runif(g, 0.05, 0.15))
  alpha0_y <- runif(ny, siga1, siga2)
  cor_matrix <- make_cor_matrix(nx, ny, rho_1, rho_2)
  
  rho_y <- matrix(1, nrow = ny, ncol = ny)
  for (m in seq_len(ny)) {
    for (q in seq_len(ny)) {
      rho_y[m, q] <- cor_matrix[(m - 1) * nx + 1, (q - 1) * nx + 1]
    }
  }
  
  betaYG <- lapply(seq_len(ny), function(m) rep(NA_real_, g))
  u_snp <- rnorm(g, 0, tau)
  for (j in seq_len(g)) {
    se_y_j <- sapply(seYG, function(x) x[j])
    cov_y_j <- rho_y * outer(se_y_j, se_y_j)
    e_y_j <- MASS::mvrnorm(1, rep(0, ny), cov_y_j)
    for (m in seq_len(ny)) {
      betaYG[[m]][j] <- beta0 * exposure$gamma_true[j] + alpha0_y[m] +
        u_snp[j] * exposure$gamma_true[j] + e_y_j[m]
    }
  }
  
  F_stats <- lapply(seq_len(nx), function(n) {
    list(individual_F = (betaXG[[n]] / seXG[[n]])^2)
  })
  
  list(
    betaYG = betaYG,
    betaXG = betaXG,
    seXG = seXG,
    seYG = seYG,
    PXG = exposure$PXG,
    alpha0_y = alpha0_y,
    rho_matrix = cor_matrix,
    F_statistics = F_stats
  )
}

meta_exposure_for_egger <- function(data_x_beta, data_x_se) {
  w <- 1 / data_x_se^2
  bx_meta <- rowSums(w * data_x_beta) / rowSums(w)
  bxse_meta <- sqrt(1 / rowSums(w))
  list(bx = bx_meta, bxse = bxse_meta)
}

weighted_egger_intercept <- function(bx, by, byse) {
  keep <- is.finite(bx) & is.finite(by) & is.finite(byse) & byse > 0
  bx <- bx[keep]
  by <- by[keep]
  byse <- byse[keep]
  
  if (length(by) < 3 || length(unique(bx)) < 2) {
    return(c(alpha = NA_real_, se = NA_real_))
  }
  
  w <- 1 / byse^2
  X <- cbind(1, bx)
  XtWX <- crossprod(X, X * w)
  XtWy <- crossprod(X, by * w)
  XtWX_inv <- MASS::ginv(XtWX)
  coef_hat <- as.numeric(XtWX_inv %*% XtWy)
  resid <- by - as.numeric(X %*% coef_hat)
  sigma2 <- sum(w * resid^2) / max(length(by) - 2, 1)
  se_alpha <- sqrt(max(0, sigma2 * XtWX_inv[1, 1]))
  
  c(alpha = coef_hat[1], se = se_alpha)
}

estimate_egger_intercepts <- function(data_x_beta, data_x_se, data_y_beta,
                                      data_y_se, rho_matrix, nx, ny,
                                      snp_indices = NULL) {
  if (is.null(snp_indices)) {
    snp_indices <- seq_len(nrow(data_x_beta))
  }
  
  exposure_meta <- meta_exposure_for_egger(data_x_beta, data_x_se)
  bx <- exposure_meta$bx[snp_indices]
  
  alpha_hat <- rep(NA_real_, ny)
  alpha_se <- rep(NA_real_, ny)
  for (m in seq_len(ny)) {
    by <- data_y_beta[snp_indices, m]
    byse <- data_y_se[snp_indices, m]
    fit <- weighted_egger_intercept(bx, by, byse)
    alpha_hat[m] <- fit[["alpha"]]
    alpha_se[m] <- fit[["se"]]
  }
  
  rho_y <- matrix(1, nrow = ny, ncol = ny)
  for (m in seq_len(ny)) {
    for (q in seq_len(ny)) {
      rho_y[m, q] <- rho_matrix[(m - 1) * nx + 1, (q - 1) * nx + 1]
    }
  }
  
  V_alpha <- rho_y * outer(alpha_se, alpha_se)
  list(alpha_hat = alpha_hat, alpha_se = alpha_se, V_alpha = V_alpha)
}

build_plugin_adjusted_inputs <- function(fdata, nx, ny, tau, valid_indices,
                                         egger_indices = NULL) {
  data_x_beta <- do.call(cbind, lapply(fdata$betaXG, matrix))
  data_x_se <- do.call(cbind, lapply(fdata$seXG, matrix))
  data_y_beta <- do.call(cbind, lapply(fdata$betaYG, matrix))
  data_y_se <- do.call(cbind, lapply(fdata$seYG, matrix))
  
  g <- nrow(data_x_beta)
  n_ratio <- nx * ny
  
  egger_fit <- estimate_egger_intercepts(
    data_x_beta = data_x_beta,
    data_x_se = data_x_se,
    data_y_beta = data_y_beta,
    data_y_se = data_y_se,
    rho_matrix = fdata$rho_matrix,
    nx = nx,
    ny = ny,
    snp_indices = egger_indices
  )
  
  WR_matrix_full <- matrix(NA_real_, nrow = n_ratio, ncol = g)
  WR_adj_full <- matrix(NA_real_, nrow = n_ratio, ncol = g)
  seWR_matrix_full <- matrix(NA_real_, nrow = n_ratio, ncol = g)
  
  for (m in seq_len(ny)) {
    for (n in seq_len(nx)) {
      row_id <- (m - 1) * nx + n
      WR_matrix_full[row_id, ] <- data_y_beta[, m] / data_x_beta[, n]
      WR_adj_full[row_id, ] <- (data_y_beta[, m] - egger_fit$alpha_hat[m]) /
        data_x_beta[, n]
      seWR_matrix_full[row_id, ] <- sqrt(data_y_se[, m]^2 / data_x_beta[, n]^2)
    }
  }
  
  A <- matrix(0, nrow = n_ratio, ncol = ny)
  for (m in seq_len(ny)) {
    row_ids <- ((m - 1) * nx + 1):(m * nx)
    A[row_ids, m] <- 1
  }
  
  WR_adj <- WR_adj_full[, valid_indices, drop = FALSE]
  seWR <- seWR_matrix_full[, valid_indices, drop = FALSE]
  beta_x_selected <- data_x_beta[valid_indices, , drop = FALSE]
  
  base_cov_list <- vector("list", length(valid_indices))
  BA_list <- vector("list", length(valid_indices))
  
  for (idx in seq_along(valid_indices)) {
    se_vector <- seWR[, idx]
    Omega_j <- fdata$rho_matrix * outer(se_vector, se_vector)
    
    beta_x_vec <- as.numeric(beta_x_selected[idx, ])
    B_j <- diag(rep(1 / beta_x_vec, times = ny), nrow = n_ratio)
    BA_j <- B_j %*% A
    
    Sigma_base_j <- Omega_j
    
    dimnames(Sigma_base_j) <- list(rownames(WR_adj), rownames(WR_adj))
    base_cov_list[[idx]] <- Sigma_base_j
    BA_list[[idx]] <- BA_j
  }
  
  list(
    WR_adj = WR_adj,
    base_cov_list = base_cov_list,
    BA_list = BA_list,
    tau2 = tau^2,
    alpha_hat = egger_fit$alpha_hat,
    alpha_se = egger_fit$alpha_se,
    V_alpha = egger_fit$V_alpha
  )
}

safe_triplet <- function(expr) {
  out <- try(expr, silent = TRUE)
  if ("try-error" %in% class(out) || any(!is.finite(out))) {
    return(c(beta = NA_real_, se = NA_real_, pval = NA_real_))
  }
  stats::setNames(as.numeric(out), c("beta", "se", "pval"))
}

run_mr_methods <- function(bx, bxse, by, byse) {
  if (length(bx) < 2 || length(by) < 2) {
    return(list(
      divw = c(beta = NA_real_, se = NA_real_, pval = NA_real_),
      ivw_fixed = c(beta = NA_real_, se = NA_real_, pval = NA_real_),
      ivw_random = c(beta = NA_real_, se = NA_real_, pval = NA_real_),
      wme = c(beta = NA_real_, se = NA_real_, pval = NA_real_)
    ))
  }
  
  mr_object <- try(MendelianRandomization::mr_input(
    bx = bx, bxse = bxse, by = by, byse = byse
  ), silent = TRUE)
  
  divw <- safe_triplet({
    fit <- mr.divw::mr.divw(bx, by, bxse, byse)
    c(fit$beta.hat, fit$beta.se, 2 * pnorm(-abs(fit$beta.hat / fit$beta.se)))
  })
  
  ivw_fixed <- if ("try-error" %in% class(mr_object)) {
    c(beta = NA_real_, se = NA_real_, pval = NA_real_)
  } else {
    safe_triplet({
      fit <- MendelianRandomization::mr_ivw(mr_object, model = "fixed")
      c(fit@Estimate, fit@StdError, fit@Pvalue)
    })
  }
  
  ivw_random <- if ("try-error" %in% class(mr_object)) {
    c(beta = NA_real_, se = NA_real_, pval = NA_real_)
  } else {
    safe_triplet({
      fit <- MendelianRandomization::mr_ivw(mr_object, model = "random")
      c(fit@Estimate, fit@StdError, fit@Pvalue)
    })
  }
  
  wme <- if (length(bx) < 3 || "try-error" %in% class(mr_object)) {
    c(beta = NA_real_, se = NA_real_, pval = NA_real_)
  } else {
    safe_triplet({
      fit <- MendelianRandomization::mr_median(mr_object, iterations = 20)
      c(fit@Estimate, fit@StdError, fit@Pvalue)
    })
  }
  
  list(divw = divw, ivw_fixed = ivw_fixed, ivw_random = ivw_random, wme = wme)
}

run_meta_result <- function(res_matrix) {
  res_matrix <- as.matrix(res_matrix)
  keep <- is.finite(res_matrix[, 1]) & is.finite(res_matrix[, 2]) & res_matrix[, 2] > 0
  res_clean <- res_matrix[keep, , drop = FALSE]
  if (nrow(res_clean) == 0) {
    return(c(beta = NA_real_, se = NA_real_, pval = NA_real_))
  }
  
  meta_res <- try(meta::metagen(res_clean[, 1], res_clean[, 2], sm = "MD"), silent = TRUE)
  if ("try-error" %in% class(meta_res)) {
    return(c(beta = NA_real_, se = NA_real_, pval = NA_real_))
  }
  
  if (is.finite(meta_res$pval.Q) && meta_res$pval.Q < 0.05) {
    c(beta = meta_res$TE.random, se = meta_res$seTE.random, pval = meta_res$pval.random)
  } else {
    c(beta = meta_res$TE.fixed, se = meta_res$seTE.fixed, pval = meta_res$pval.fixed)
  }
}

run_comparison_methods <- function(fdata, nx, ny, g, F_THRESHOLD, P_THRESHOLD) {
  f_stat_matrix <- do.call(cbind, lapply(fdata$F_statistics, function(x) x$individual_F))
  pval_stat_matrix <- do.call(cbind, lapply(fdata$PXG, matrix))
  data_x_beta <- do.call(cbind, lapply(fdata$betaXG, matrix))
  data_x_se <- do.call(cbind, lapply(fdata$seXG, matrix))
  data_y_beta <- do.call(cbind, lapply(fdata$betaYG, matrix))
  data_y_se <- do.call(cbind, lapply(fdata$seYG, matrix))
  
  alpha_data <- matrix(NA_real_, nrow = g, ncol = 3)
  theta_data <- matrix(NA_real_, nrow = g, ncol = 3)
  for (j in seq_len(g)) {
    retx <- try(meta::metagen(data_x_beta[j, ], data_x_se[j, ], sm = "MD"), silent = TRUE)
    rety <- try(meta::metagen(data_y_beta[j, ], data_y_se[j, ], sm = "MD"), silent = TRUE)
    if (!("try-error" %in% class(retx))) {
      if (is.finite(retx$pval.Q) && retx$pval.Q < 0.05) {
        alpha_data[j, ] <- c(retx$TE.random, retx$seTE.random, retx$pval.random)
      } else {
        alpha_data[j, ] <- c(retx$TE.fixed, retx$seTE.fixed, retx$pval.fixed)
      }
    }
    if (!("try-error" %in% class(rety))) {
      if (is.finite(rety$pval.Q) && rety$pval.Q < 0.05) {
        theta_data[j, ] <- c(rety$TE.random, rety$seTE.random, rety$pval.random)
      } else {
        theta_data[j, ] <- c(rety$TE.fixed, rety$seTE.fixed, rety$pval.fixed)
      }
    }
  }
  
  F_stat_meta <- (alpha_data[, 1] / alpha_data[, 2])^2
  valid_indices_meta <- which(
    is.finite(F_stat_meta) &
      is.finite(alpha_data[, 3]) &
      F_stat_meta > F_THRESHOLD &
      alpha_data[, 3] < P_THRESHOLD
  )
  SNP_GWAS_meta <- length(valid_indices_meta)
  
  meta_mr <- run_mr_methods(
    bx = alpha_data[valid_indices_meta, 1],
    bxse = alpha_data[valid_indices_meta, 2],
    by = theta_data[valid_indices_meta, 1],
    byse = theta_data[valid_indices_meta, 2]
  )
  gwas_meta_res <- c(meta_mr$divw, meta_mr$ivw_fixed, meta_mr$ivw_random, meta_mr$wme)
  
  pair_divw <- matrix(NA_real_, nrow = nx * ny, ncol = 3)
  pair_ivw_fixed <- matrix(NA_real_, nrow = nx * ny, ncol = 3)
  pair_ivw_random <- matrix(NA_real_, nrow = nx * ny, ncol = 3)
  pair_wme <- matrix(NA_real_, nrow = nx * ny, ncol = 3)
  SNPSMR <- rep(NA_real_, nx)
  
  row_id <- 0
  for (n in seq_len(nx)) {
    valid_indices_trad <- which(
      f_stat_matrix[, n] > F_THRESHOLD &
        pval_stat_matrix[, n] < P_THRESHOLD
    )
    SNPSMR[n] <- length(valid_indices_trad)
    
    for (m in seq_len(ny)) {
      row_id <- row_id + 1
      pair_mr <- run_mr_methods(
        bx = data_x_beta[valid_indices_trad, n],
        bxse = data_x_se[valid_indices_trad, n],
        by = data_y_beta[valid_indices_trad, m],
        byse = data_y_se[valid_indices_trad, m]
      )
      pair_divw[row_id, ] <- pair_mr$divw
      pair_ivw_fixed[row_id, ] <- pair_mr$ivw_fixed
      pair_ivw_random[row_id, ] <- pair_mr$ivw_random
      pair_wme[row_id, ] <- pair_mr$wme
    }
  }
  
  SNPMR <- round(mean(SNPSMR, na.rm = TRUE), 0)
  pairwise_res <- c(
    pair_divw[, 1], pair_ivw_fixed[, 1], pair_ivw_random[, 1], pair_wme[, 1],
    pair_divw[, 2], pair_ivw_fixed[, 2], pair_ivw_random[, 2], pair_wme[, 2],
    pair_divw[, 3], pair_ivw_fixed[, 3], pair_ivw_random[, 3], pair_wme[, 3]
  )
  
  mr_meta_res <- c(
    run_meta_result(pair_divw),
    run_meta_result(pair_ivw_fixed),
    run_meta_result(pair_ivw_random),
    run_meta_result(pair_wme)
  )
  
  list(
    SNP_GWAS_meta = SNP_GWAS_meta,
    SNPMR = SNPMR,
    gwas_meta_res = gwas_meta_res,
    pairwise_res = pairwise_res,
    mr_meta_res = mr_meta_res
  )
}

comparison_result_names <- function(nx, ny) {
  methods <- c("dIVW", "IVW_fixed", "IVW_random", "wme")
  pair_count <- nx * ny
  pair_names <- c(
    unlist(lapply(methods, function(method) {
      paste0("MR_", method, "_beta", seq_len(pair_count))
    })),
    unlist(lapply(methods, function(method) {
      paste0("MR_", method, "_se", seq_len(pair_count))
    })),
    unlist(lapply(methods, function(method) {
      paste0("MR_", method, "_pval", seq_len(pair_count))
    }))
  )
  
  c(
    paste0("meta_res_dIVW_", c("beta", "se", "pval")),
    paste0("meta_res_IVW_fixed_", c("beta", "se", "pval")),
    paste0("meta_res_IVW_random_", c("beta", "se", "pval")),
    paste0("meta_res_wme_", c("beta", "se", "pval")),
    pair_names,
    paste0("MR_meta_old_dIVW_", c("beta", "se", "pval")),
    paste0("MR_meta_old_ivw_fixed_", c("beta", "se", "pval")),
    paste0("MR_meta_old_ivw_random_", c("beta", "se", "pval")),
    paste0("MR_meta_old_wme_", c("beta", "se", "pval"))
  )
}

as_result_row <- function(values) {
  as.data.frame(as.list(values), check.names = FALSE)
}

run_once <- function(i, g, nx, ny, bx, beta0, tau, rho_1, rho_2,
                     siga1, siga2, weakp, F_THRESHOLD, P_THRESHOLD,
                     seed_base = 20260428, BOOT_B = 200) {
  set.seed(seed_base + i)
  
  fdata <- DataGeneration_summary(g, nx, ny, bx, beta0, tau, rho_1, rho_2,
                                  siga1, siga2, weakp)
  comparison <- run_comparison_methods(
    fdata = fdata,
    nx = nx,
    ny = ny,
    g = g,
    F_THRESHOLD = F_THRESHOLD,
    P_THRESHOLD = P_THRESHOLD
  )
  comparison_values <- c(
    comparison$gwas_meta_res,
    comparison$pairwise_res,
    comparison$mr_meta_res
  )
  names(comparison_values) <- comparison_result_names(nx, ny)
  
  f_stat_matrix <- do.call(cbind, lapply(fdata$F_statistics, function(x) x$individual_F))
  pval_stat_matrix <- do.call(cbind, lapply(fdata$PXG, matrix))
  
  valid_indices <- which(
    apply(f_stat_matrix, 1, function(row) all(row >= F_THRESHOLD)) &
      apply(pval_stat_matrix, 1, function(row) any(row < P_THRESHOLD))
  )
  
  base_values <- c(
    iter = i,
    nx = nx,
    ny = ny,
    g = g,
    beta0 = beta0,
    rho_mean = rho_1,
    tau = tau,
    siga1 = siga1,
    siga2 = siga2,
    SNP_gaws_meta = comparison$SNP_GWAS_meta,
    SNPMR = comparison$SNPMR,
    SNPnew = length(valid_indices)
  )
  
  if (length(valid_indices) < 2) {
    return(as_result_row(c(
      base_values,
      comparison_values,
      JointMR_plugin_beta = NA_real_,
      JointMR_plugin_se = NA_real_,
      JointMR_plugin_pval = NA_real_,
      theta_hat = NA_real_,
      theta_se = NA_real_,
      theta_p_value = NA_real_,
      alpha0_mean = NA_real_,
      alpha0_true_mean = mean(fdata$alpha0_y)
    )))
  }
  
  plugin_inputs <- build_plugin_adjusted_inputs(
    fdata = fdata,
    nx = nx,
    ny = ny,
    tau = tau,
    valid_indices = valid_indices,
    egger_indices = seq_len(g)
  )
  
  est <- gls_theta_plugin(
    base_cov_list = plugin_inputs$base_cov_list,
    WR_matrix = plugin_inputs$WR_adj,
    BA_list = plugin_inputs$BA_list,
    V_alpha = plugin_inputs$V_alpha,
    tau2 = plugin_inputs$tau2
  )
  boot_est <- bootstrap_theta(
    fdata = fdata,
    nx = nx,
    ny = ny,
    tau = tau,
    valid_indices = valid_indices,
    theta_hat = unname(est[["theta_hat"]]),
    B = BOOT_B
  )
  
  as_result_row(c(
    base_values,
    comparison_values,
    JointMR_plugin_beta = unname(est[["theta_hat"]]),
    JointMR_plugin_se = unname(boot_est[["theta_se_boot"]]),
    JointMR_plugin_pval = unname(boot_est[["theta_p_value_boot"]]),
    theta_hat = unname(est[["theta_hat"]]),
    theta_se = unname(boot_est[["theta_se_boot"]]),
    theta_p_value = unname(boot_est[["theta_p_value_boot"]]),
    alpha0_mean = mean(plugin_inputs$alpha_hat),
    alpha0_true_mean = mean(fdata$alpha0_y)
  ))
}

bootstrap_theta <- function(fdata, nx, ny, tau, valid_indices, theta_hat,
                            B = 200) {
  boot_theta <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    boot_indices <- sample(valid_indices, length(valid_indices), replace = TRUE)
    boot_inputs <- try(
      build_plugin_adjusted_inputs(
        fdata = fdata,
        nx = nx,
        ny = ny,
        tau = tau,
        valid_indices = boot_indices,
        egger_indices = boot_indices
      ),
      silent = TRUE
    )
    if ("try-error" %in% class(boot_inputs)) {
      next
    }
    
    boot_est <- try(
      gls_theta_plugin(
        base_cov_list = boot_inputs$base_cov_list,
        WR_matrix = boot_inputs$WR_adj,
        BA_list = boot_inputs$BA_list,
        V_alpha = boot_inputs$V_alpha,
        tau2 = boot_inputs$tau2
      ),
      silent = TRUE
    )
    if (!("try-error" %in% class(boot_est))) {
      boot_theta[b] <- unname(boot_est[["theta_hat"]])
    }
  }
  
  boot_se <- sd(boot_theta, na.rm = TRUE)
  boot_p_value <- if (is.finite(boot_se) && boot_se > 0) {
    2 * pnorm(-abs(theta_hat / boot_se))
  } else {
    NA_real_
  }
  
  c(theta_se_boot = boot_se, theta_p_value_boot = boot_p_value)
}

build_ccname <- function(nx, ny) {
  methods <- c("dIVW", "IVW_fixed", "IVW_random", "wme")
  all_beta_names <- unlist(lapply(methods, function(method) {
    paste0("MR_", method, "_", rep("beta", nx * ny), seq_len(nx * ny))
  }))
  
  all_se_names <- unlist(lapply(methods, function(method) {
    paste0("MR_", method, "_", rep("se", nx * ny), seq_len(nx * ny))
  }))
  
  all_pval_names <- unlist(lapply(methods, function(method) {
    paste0("MR_", method, "_", rep("pval", nx * ny), seq_len(nx * ny))
  }))
  
  c(
    "nx", "ny", "g", "beta0", "rho_mean", "tau", "siga1",
    "siga2", "SNP_gaws_meta", "SNPMR", "SNPnew",
    paste0("meta_res_dIVW_", c("beta", "se", "pval")),
    paste0("meta_res_IVW_fixed_", c("beta", "se", "pval")),
    paste0("meta_res_IVW_random_", c("beta", "se", "pval")),
    paste0("meta_res_wme_", c("beta", "se", "pval")),
    all_beta_names, all_se_names, all_pval_names,
    paste0("MR_meta_old_dIVW_", c("beta", "se", "pval")),
    paste0("MR_meta_old_ivw_fixed_", c("beta", "se", "pval")),
    paste0("MR_meta_old_ivw_random_", c("beta", "se", "pval")),
    paste0("MR_meta_old_wme_", c("beta", "se", "pval")),
    paste0("MRmeta_", c("beta", "se", "pval")), "F_tat", "weakp",
    "alpha0_hat_mean", "alpha0_se_mean", "alpha0_true_mean",
    paste0("alpha0_hat_Y", seq_len(ny)),
    paste0("alpha0_se_Y", seq_len(ny)),
    paste0("alpha0_true_Y", seq_len(ny))
  )
}

#####Comp funcrion######
Comp <- function(i, g, nx, ny, bx, beta0, tau, rho_1, rho_2, siga1, siga2,
                 weakp, F_THRESHOLD, P_THRESHOLD, BOOT_B = 200) {
  cat("i=", i, "\n")
  set.seed(20260428 + i)
  
  fdata1 <- DataGeneration_summary(g, nx, ny, bx, beta0, tau, rho_1, rho_2,
                                   siga1, siga2, weakp)
  F_tat <- mean(sapply(fdata1$F_statistics, function(x) {
    mean(x$individual_F, na.rm = TRUE)
  }), na.rm = TRUE)
  
  comparison <- run_comparison_methods(
    fdata = fdata1,
    nx = nx,
    ny = ny,
    g = g,
    F_THRESHOLD = F_THRESHOLD,
    P_THRESHOLD = P_THRESHOLD
  )
  
  f_stat_matrix <- do.call(cbind, lapply(fdata1$F_statistics, function(x) {
    x$individual_F
  }))
  pval_stat_matrix <- do.call(cbind, lapply(fdata1$PXG, matrix))
  
  valid_indices_intersect <- which(
    apply(f_stat_matrix, 1, function(row) all(row >= F_THRESHOLD)) &
      apply(pval_stat_matrix, 1, function(row) any(row < P_THRESHOLD))
  )
  SNPnew <- length(valid_indices_intersect)
  
  theta_hat <- NA_real_
  theta_se <- NA_real_
  theta_p_value <- NA_real_
  alpha0_hat <- rep(NA_real_, ny)
  alpha0_se <- rep(NA_real_, ny)
  
  if (SNPnew >= 2) {
    plugin_inputs <- build_plugin_adjusted_inputs(
      fdata = fdata1,
      nx = nx,
      ny = ny,
      tau = tau,
      valid_indices = valid_indices_intersect,
      egger_indices = seq_len(g)
    )
    
    est <- gls_theta_plugin(
      base_cov_list = plugin_inputs$base_cov_list,
      WR_matrix = plugin_inputs$WR_adj,
      BA_list = plugin_inputs$BA_list,
      V_alpha = plugin_inputs$V_alpha,
      tau2 = plugin_inputs$tau2
    )
    
    theta_hat <- unname(est[["theta_hat"]])
    boot_est <- bootstrap_theta(
      fdata = fdata1,
      nx = nx,
      ny = ny,
      tau = tau,
      valid_indices = valid_indices_intersect,
      theta_hat = theta_hat,
      B = BOOT_B
    )
    theta_se <- unname(boot_est[["theta_se_boot"]])
    theta_p_value <- unname(boot_est[["theta_p_value_boot"]])
    alpha0_hat <- plugin_inputs$alpha_hat
    alpha0_se <- plugin_inputs$alpha_se
  }
  
  res_all <- c(
    nx, ny, g, beta0, rho_1, tau, siga1, siga2,
    comparison$SNP_GWAS_meta, comparison$SNPMR, SNPnew,
    comparison$gwas_meta_res,
    comparison$pairwise_res,
    comparison$mr_meta_res,
    theta_hat, theta_se, theta_p_value, F_tat, weakp,
    mean(alpha0_hat, na.rm = TRUE),
    mean(alpha0_se, na.rm = TRUE),
    mean(fdata1$alpha0_y, na.rm = TRUE),
    alpha0_hat,
    alpha0_se,
    fdata1$alpha0_y
  )
  
  res_all <- matrix(res_all, nrow = 1)
  res_all <- as.data.frame(res_all)
  colnames(res_all) <- build_ccname(nx, ny)
  
  res_all
}

####mian function####
main <- function(NN, g, nx, ny, bx, beta0, tau, rho_1, rho_2, siga1, siga2,
                 weakp, F_THRESHOLD, P_THRESHOLD, mc, BOOT_B = 200) {
                  
  registerDoMC(mc)
  res_all <- foreach(i=1:NN,
                     .combine=rbind,
                     .errorhandling = "remove",
                     .packages = c("MASS", "meta", "MendelianRandomization",
                                   "mr.divw", "readr")) %dopar% {
    Comp(i, g, nx, ny, bx, beta0, tau, rho_1, rho_2, siga1, siga2,
         weakp, F_THRESHOLD, P_THRESHOLD, BOOT_B)
  }
  
  write_csv(res_all, file = paste0("Sum_nx_", nx, "_ny_", ny, "_g_", g,
                                   "_beta0_", beta0, "_rho1_", rho_1,
                                   "_tau_", tau, "_siga1_", siga1,
                                   "_siga2_", siga2, "_weakp_", weakp,
                                   "_F_THRESHOLD_", F_THRESHOLD,
                                   "_P_THRESHOLD_", P_THRESHOLD,
                        "_0501_bootstrap_ple_plugin_delta.csv"))
  invisible(res_all)
}

F_THRESHOLDb <- c(10)
P_THRESHOLD <- 5e-8
NN <- 1000
mc <- 50
BOOT_B <- 200
bx <- 0.2
gb <- c(200)
beta0b <- c(0, 0.05)
taub <- c(0, 0.2)
rho_1b <- c(0.5, 0.8, 0.85, 0.9, 0.95)
rho_2b <- c(0.5, 0.8, 0.85, 0.9, 0.95)
siga1b <- c(0)
siga2b <- c(0.03)
nxb <- c(3)
nyb <- c(5)
weakpb <- c(0)

for (zzh in seq_along(F_THRESHOLDb)) {
  F_THRESHOLD <- F_THRESHOLDb[zzh]
  for (yh in seq_along(taub)) {
    tau <- taub[[yh]]
    for (zh in seq_along(nxb)) {
      nx <- nxb[[zh]]
      for (qh in seq_along(nyb)) {
        ny <- nyb[[qh]]
        for (xh in seq_along(siga1b)) {
          siga1 <- siga1b[[xh]]
          siga2 <- siga2b[[xh]]
          for (tth in seq_along(weakpb)) {
            weakp <- weakpb[tth]
            for (th in seq_along(gb)) {
              g <- gb[th]
              for (uh in seq_along(beta0b)) {
                beta0 <- beta0b[[uh]]
                for (rh in seq_along(rho_1b)) {
                  rho_1 <- rho_1b[[rh]]
                  rho_2 <- rho_2b[[rh]]
                  
                  ccname <- build_ccname(nx, ny)
                  
                  main(NN, g, nx, ny, bx, beta0, tau, rho_1, rho_2,
                       siga1, siga2, weakp, F_THRESHOLD, P_THRESHOLD, mc,
                       BOOT_B)
                }
              }
            }
          }
        }
      }
    }
  }
}
