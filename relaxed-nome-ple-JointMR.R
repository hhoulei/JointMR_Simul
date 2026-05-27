setwd("/home/houl/MR-META/20260421")

library(MASS)
library(readr)
library(meta)
library(MendelianRandomization)
library(mr.divw)
library(foreach)
library(doMC)

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
  
  if (all(is.finite(V_alpha)) && max(abs(V_alpha)) <= .Machine$double.eps) {
    Sigma_inv_one <- D_inv_one
    Sigma_inv_y <- D_inv_y
  } else {
    middle <- MASS::ginv(V_alpha) + t(U) %*% D_inv_U
    middle_inv <- MASS::ginv((middle + t(middle)) / 2)
    
    Sigma_inv_one <- D_inv_one - D_inv_U %*% middle_inv %*%
      t(U) %*% D_inv_one
    Sigma_inv_y <- D_inv_y - D_inv_U %*% middle_inv %*%
      t(U) %*% D_inv_y
  }
  
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
  strong_instruments_base <- rnorm(strong_count, bx, 0.03)
  weak_instruments_base <- rnorm(weak_count, 0, 0.005)
  shared_instrument_order <- sample(c(strong_instruments_base, weak_instruments_base))
  shared_instrument_order <- shared_instrument_order * sample(c(-1, 1), g, replace = TRUE)
  
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
    betaXG_initial = betaXG_initial,
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

make_psd <- function(Sigma, eps = 1e-8) {
  Sigma <- (Sigma + t(Sigma)) / 2
  eig <- eigen(Sigma, symmetric = TRUE)
  eig$values <- pmax(eig$values, eps)
  eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
}

DataGeneration_summary <- function(g, nx, ny, bx, beta0, tau, rho_1, rho_2,
                                   siga1, siga2, weakp,
                                   exposure_rho = 0, xy_rho = 0) {
  exposure <- make_common_exposure(g, nx, bx, weakp)
  betaXG_initial <- exposure$betaXG_initial
  seXG <- exposure$seXG
  seYG <- lapply(seq_len(ny), function(i) runif(g, 0.05, 0.15))
  alpha0_y <- runif(ny, siga1, siga2)
  cor_matrix <- make_cor_matrix(nx, ny, rho_1, rho_2)
  rho_x <- make_exposure_cor(nx, exposure_rho)
  
  rho_y <- matrix(1, nrow = ny, ncol = ny)
  for (m in seq_len(ny)) {
    for (q in seq_len(ny)) {
      rho_y[m, q] <- cor_matrix[(m - 1) * nx + 1, (q - 1) * nx + 1]
    }
  }
  
  betaXG <- lapply(seq_len(nx), function(n) rep(NA_real_, g))
  betaYG <- lapply(seq_len(ny), function(m) rep(NA_real_, g))
  u_snp <- rnorm(g, 0, tau)
  for (j in seq_len(g)) {
    se_x_j <- sapply(seXG, function(x) x[j])
    se_y_j <- sapply(seYG, function(x) x[j])
    cov_x_j <- rho_x * outer(se_x_j, se_x_j)
    cov_y_j <- rho_y * outer(se_y_j, se_y_j)
    cov_xy_j <- matrix(xy_rho, nrow = nx, ncol = ny) * outer(se_x_j, se_y_j)
    cov_all_j <- rbind(
      cbind(cov_x_j, cov_xy_j),
      cbind(t(cov_xy_j), cov_y_j)
    )
    e_all_j <- MASS::mvrnorm(1, rep(0, nx + ny), make_psd(cov_all_j))
    for (n in seq_len(nx)) {
      betaXG[[n]][j] <- betaXG_initial[[n]][j] + e_all_j[n]
    }
    for (m in seq_len(ny)) {
      betaYG[[m]][j] <- beta0 * exposure$gamma_true[j] + alpha0_y[m] +
        u_snp[j] * exposure$gamma_true[j] + e_all_j[nx + m]
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
    fit <- lm(by ~ bx, weights = 1 / byse^2)
    fit_sum <- summary(fit)
    alpha_hat[m] <- coef(fit)[1]
    alpha_se[m] <- fit_sum$coefficients[1, 2]
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

extract_outcome_cor <- function(rho_matrix, nx, ny) {
  rho_y <- matrix(1, nrow = ny, ncol = ny)
  for (m in seq_len(ny)) {
    for (q in seq_len(ny)) {
      rho_y[m, q] <- rho_matrix[(m - 1) * nx + 1, (q - 1) * nx + 1]
    }
  }
  rho_y
}

make_exposure_cor <- function(nx, exposure_rho = 0) {
  rho_x <- matrix(exposure_rho, nrow = nx, ncol = nx)
  diag(rho_x) <- 1
  rho_x
}

delta_ratio_cov <- function(beta_x, se_x, beta_y, se_y, alpha_hat,
                            rho_y, rho_x, nx, ny, theta_cov = NULL,
                            xy_rho = 0) {
  n_ratio <- nx * ny
  cov_matrix <- matrix(0, nrow = n_ratio, ncol = n_ratio)
  y_adj <- beta_y - alpha_hat
  
  for (m in seq_len(ny)) {
    for (n in seq_len(nx)) {
      row_id <- (m - 1) * nx + n
      for (q in seq_len(ny)) {
        for (p in seq_len(nx)) {
          col_id <- (q - 1) * nx + p
          cov_y <- rho_y[m, q] * se_y[m] * se_y[q]
          cov_x <- rho_x[n, p] * se_x[n] * se_x[p]
          if (is.null(theta_cov)) {
            exposure_component <- y_adj[m] * y_adj[q] * cov_x /
              (beta_x[n]^2 * beta_x[p]^2)
            cross_component <- 0
          } else {
            exposure_component <- theta_cov^2 * cov_x /
              (beta_x[n] * beta_x[p])
            cov_y_xp <- xy_rho * se_y[m] * se_x[p]
            cov_xn_yq <- xy_rho * se_x[n] * se_y[q]
            cross_component <- -theta_cov * (cov_y_xp + cov_xn_yq) /
              (beta_x[n] * beta_x[p])
          }
          cov_matrix[row_id, col_id] <-
            cov_y / (beta_x[n] * beta_x[p]) +
            exposure_component +
            cross_component
        }
      }
    }
  }
  
  cov_matrix <- (cov_matrix + t(cov_matrix)) / 2
  diag(cov_matrix) <- pmax(diag(cov_matrix), .Machine$double.eps)
  cov_matrix
}

build_relaxed_nome_inputs <- function(fdata, nx, ny, tau, valid_indices,
                                      adjust_pleiotropy = TRUE,
                                      egger_indices = NULL,
                                      exposure_rho = 0,
                                      xy_rho = 0,
                                      theta_cov = NULL) {
  data_x_beta <- do.call(cbind, lapply(fdata$betaXG, matrix))
  data_x_se <- do.call(cbind, lapply(fdata$seXG, matrix))
  data_y_beta <- do.call(cbind, lapply(fdata$betaYG, matrix))
  data_y_se <- do.call(cbind, lapply(fdata$seYG, matrix))
  
  g <- nrow(data_x_beta)
  n_ratio <- nx * ny
  
  if (adjust_pleiotropy) {
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
  } else {
    egger_fit <- list(
      alpha_hat = rep(0, ny),
      alpha_se = rep(0, ny),
      V_alpha = matrix(0, nrow = ny, ncol = ny)
    )
  }
  
  WR_matrix_full <- matrix(NA_real_, nrow = n_ratio, ncol = g)
  WR_adj_full <- matrix(NA_real_, nrow = n_ratio, ncol = g)
  
  for (m in seq_len(ny)) {
    for (n in seq_len(nx)) {
      row_id <- (m - 1) * nx + n
      WR_matrix_full[row_id, ] <- data_y_beta[, m] / data_x_beta[, n]
      WR_adj_full[row_id, ] <- (data_y_beta[, m] - egger_fit$alpha_hat[m]) /
        data_x_beta[, n]
    }
  }
  
  A <- matrix(0, nrow = n_ratio, ncol = ny)
  for (m in seq_len(ny)) {
    row_ids <- ((m - 1) * nx + 1):(m * nx)
    A[row_ids, m] <- 1
  }
  
  WR_adj <- WR_adj_full[, valid_indices, drop = FALSE]
  beta_x_selected <- data_x_beta[valid_indices, , drop = FALSE]
  se_x_selected <- data_x_se[valid_indices, , drop = FALSE]
  beta_y_selected <- data_y_beta[valid_indices, , drop = FALSE]
  se_y_selected <- data_y_se[valid_indices, , drop = FALSE]
  
  base_cov_list <- vector("list", length(valid_indices))
  BA_list <- vector("list", length(valid_indices))
  rho_y <- extract_outcome_cor(fdata$rho_matrix, nx, ny)
  rho_x <- make_exposure_cor(nx, exposure_rho = exposure_rho)
  
  for (idx in seq_along(valid_indices)) {
    beta_x_vec <- as.numeric(beta_x_selected[idx, ])
    se_x_vec <- as.numeric(se_x_selected[idx, ])
    beta_y_vec <- as.numeric(beta_y_selected[idx, ])
    se_y_vec <- as.numeric(se_y_selected[idx, ])
    
    Omega_j <- delta_ratio_cov(
      beta_x = beta_x_vec,
      se_x = se_x_vec,
      beta_y = beta_y_vec,
      se_y = se_y_vec,
      alpha_hat = egger_fit$alpha_hat,
      rho_y = rho_y,
      rho_x = rho_x,
      nx = nx,
      ny = ny,
      theta_cov = theta_cov,
      xy_rho = xy_rho
    )
    
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

build_plugin_adjusted_inputs <- build_relaxed_nome_inputs

estimate_jointmr_relaxed_nome <- function(fdata, nx, ny, tau, valid_indices,
                                          adjust_pleiotropy = TRUE,
                                          egger_indices = NULL,
                                          exposure_rho = 0,
                                          xy_rho = 0,
                                          max_iter = 30,
                                          tol = 1e-7) {
  data_x_beta <- do.call(cbind, lapply(fdata$betaXG, matrix))
  data_x_se <- do.call(cbind, lapply(fdata$seXG, matrix))
  data_y_beta <- do.call(cbind, lapply(fdata$betaYG, matrix))
  data_y_se <- do.call(cbind, lapply(fdata$seYG, matrix))
  if (adjust_pleiotropy) {
    egger_fit_init <- estimate_egger_intercepts(
      data_x_beta = data_x_beta,
      data_x_se = data_x_se,
      data_y_beta = data_y_beta,
      data_y_se = data_y_se,
      rho_matrix = fdata$rho_matrix,
      nx = nx,
      ny = ny,
      snp_indices = egger_indices
    )
    alpha_init <- egger_fit_init$alpha_hat
  } else {
    alpha_init <- rep(0, ny)
  }
  theta_start_values <- c()
  for (m in seq_len(ny)) {
    for (n in seq_len(nx)) {
      bx <- data_x_beta[valid_indices, n]
      by <- data_y_beta[valid_indices, m] - alpha_init[m]
      theta_start_values <- c(theta_start_values, sum(bx * by) / sum(bx^2))
    }
  }
  theta_old <- mean(theta_start_values, na.rm = TRUE)
  if (!is.finite(theta_old)) {
    theta_old <- 0
  }
  
  init_inputs <- build_relaxed_nome_inputs(
    fdata = fdata,
    nx = nx,
    ny = ny,
    tau = tau,
    valid_indices = valid_indices,
    adjust_pleiotropy = adjust_pleiotropy,
    egger_indices = egger_indices,
    exposure_rho = exposure_rho,
    xy_rho = xy_rho,
    theta_cov = theta_old
  )
  init_est <- gls_theta_plugin(
    base_cov_list = init_inputs$base_cov_list,
    WR_matrix = init_inputs$WR_adj,
    BA_list = init_inputs$BA_list,
    V_alpha = init_inputs$V_alpha,
    tau2 = init_inputs$tau2
  )
  
  theta_old <- unname(init_est[["theta_hat"]])
  inputs <- init_inputs
  est <- init_est
  for (iter in seq_len(max_iter)) {
    inputs <- build_relaxed_nome_inputs(
      fdata = fdata,
      nx = nx,
      ny = ny,
      tau = tau,
      valid_indices = valid_indices,
      adjust_pleiotropy = adjust_pleiotropy,
      egger_indices = egger_indices,
      exposure_rho = exposure_rho,
      xy_rho = xy_rho,
      theta_cov = theta_old
    )
    est <- gls_theta_plugin(
      base_cov_list = inputs$base_cov_list,
      WR_matrix = inputs$WR_adj,
      BA_list = inputs$BA_list,
      V_alpha = inputs$V_alpha,
      tau2 = inputs$tau2
    )
    theta_new <- unname(est[["theta_hat"]])
    if (is.finite(theta_new) && abs(theta_new - theta_old) < tol) {
      break
    }
    theta_old <- theta_new
  }
  
  list(est = est, inputs = inputs)
}

stack_outcome_values <- function(beta_y_vec, nx) {
  rep(beta_y_vec, each = nx)
}

stack_exposure_values <- function(beta_x_vec, ny) {
  rep(beta_x_vec, times = ny)
}

make_alpha_map <- function(nx, ny) {
  A <- matrix(0, nrow = nx * ny, ncol = ny)
  for (m in seq_len(ny)) {
    row_ids <- ((m - 1) * nx + 1):(m * nx)
    A[row_ids, m] <- 1
  }
  A
}

association_residual_cov <- function(theta, beta_x, se_x, se_y, rho_y, rho_x,
                                     nx, ny, tau2 = 0, V_alpha = NULL,
                                     xy_rho = 0) {
  n_ratio <- nx * ny
  cov_matrix <- matrix(0, nrow = n_ratio, ncol = n_ratio)
  
  for (m in seq_len(ny)) {
    for (n in seq_len(nx)) {
      row_id <- (m - 1) * nx + n
      for (q in seq_len(ny)) {
        for (p in seq_len(nx)) {
          col_id <- (q - 1) * nx + p
          cov_y <- rho_y[m, q] * se_y[m] * se_y[q]
          cov_x <- rho_x[n, p] * se_x[n] * se_x[p]
          cov_y_xp <- xy_rho * se_y[m] * se_x[p]
          cov_xn_yq <- xy_rho * se_x[n] * se_y[q]
          cov_matrix[row_id, col_id] <-
            cov_y + theta^2 * cov_x - theta * (cov_y_xp + cov_xn_yq)
        }
      }
    }
  }
  
  if (tau2 > 0) {
    x_stack <- stack_exposure_values(beta_x, ny)
    cov_matrix <- cov_matrix + tau2 * tcrossprod(x_stack)
  }
  
  if (!is.null(V_alpha) && all(is.finite(V_alpha)) &&
      max(abs(V_alpha)) > .Machine$double.eps) {
    A <- make_alpha_map(nx, ny)
    cov_matrix <- cov_matrix + A %*% V_alpha %*% t(A)
  }
  
  cov_matrix <- (cov_matrix + t(cov_matrix)) / 2
  diag(cov_matrix) <- pmax(diag(cov_matrix), .Machine$double.eps)
  make_psd(cov_matrix)
}

prepare_jointmr_assoc_data <- function(fdata, nx, ny, valid_indices,
                                       adjust_pleiotropy = TRUE,
                                       egger_indices = NULL,
                                       exposure_rho = 0) {
  data_x_beta <- do.call(cbind, lapply(fdata$betaXG, matrix))
  data_x_se <- do.call(cbind, lapply(fdata$seXG, matrix))
  data_y_beta <- do.call(cbind, lapply(fdata$betaYG, matrix))
  data_y_se <- do.call(cbind, lapply(fdata$seYG, matrix))
  
  if (adjust_pleiotropy) {
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
  } else {
    egger_fit <- list(
      alpha_hat = rep(0, ny),
      alpha_se = rep(0, ny),
      V_alpha = matrix(0, nrow = ny, ncol = ny)
    )
  }
  
  list(
    beta_x = data_x_beta[valid_indices, , drop = FALSE],
    se_x = data_x_se[valid_indices, , drop = FALSE],
    beta_y = data_y_beta[valid_indices, , drop = FALSE],
    se_y = data_y_se[valid_indices, , drop = FALSE],
    rho_y = extract_outcome_cor(fdata$rho_matrix, nx, ny),
    rho_x = make_exposure_cor(nx, exposure_rho = exposure_rho),
    alpha_hat = egger_fit$alpha_hat,
    alpha_se = egger_fit$alpha_se,
    V_alpha = egger_fit$V_alpha
  )
}

# S4 uses an association-scale exact-Q residual:
#   betaY_jm - alpha_m - theta * betaX_jn
# rather than a Wald-ratio residual. This overrides the earlier ratio-scale
# helper while keeping the same return shape used by the simulation loop.
build_relaxed_nome_inputs <- function(fdata, nx, ny, tau, valid_indices,
                                      adjust_pleiotropy = TRUE,
                                      egger_indices = NULL,
                                      exposure_rho = 0,
                                      xy_rho = 0,
                                      theta = 0,
                                      theta_cov = NULL) {
  if (!is.null(theta_cov)) {
    theta <- theta_cov
  }
  
  prep <- prepare_jointmr_assoc_data(
    fdata = fdata,
    nx = nx,
    ny = ny,
    valid_indices = valid_indices,
    adjust_pleiotropy = adjust_pleiotropy,
    egger_indices = egger_indices,
    exposure_rho = exposure_rho
  )
  
  n_snp <- length(valid_indices)
  n_ratio <- nx * ny
  residual_matrix <- matrix(NA_real_, nrow = n_ratio, ncol = n_snp)
  base_cov_list <- vector("list", n_snp)
  
  for (idx in seq_len(n_snp)) {
    beta_x_vec <- as.numeric(prep$beta_x[idx, ])
    se_x_vec <- as.numeric(prep$se_x[idx, ])
    beta_y_vec <- as.numeric(prep$beta_y[idx, ])
    se_y_vec <- as.numeric(prep$se_y[idx, ])
    
    y_stack <- stack_outcome_values(beta_y_vec, nx)
    x_stack <- stack_exposure_values(beta_x_vec, ny)
    alpha_stack <- stack_outcome_values(prep$alpha_hat, nx)
    residual_matrix[, idx] <- y_stack - alpha_stack - theta * x_stack
    
    base_cov_list[[idx]] <- association_residual_cov(
      theta = theta,
      beta_x = beta_x_vec,
      se_x = se_x_vec,
      se_y = se_y_vec,
      rho_y = prep$rho_y,
      rho_x = prep$rho_x,
      nx = nx,
      ny = ny,
      tau2 = tau^2,
      V_alpha = prep$V_alpha,
      xy_rho = xy_rho
    )
  }
  
  list(
    residual_matrix = residual_matrix,
    WR_adj = residual_matrix,
    base_cov_list = base_cov_list,
    tau2 = tau^2,
    alpha_hat = prep$alpha_hat,
    alpha_se = prep$alpha_se,
    V_alpha = prep$V_alpha
  )
}

build_plugin_adjusted_inputs <- build_relaxed_nome_inputs

jointmr_assoc_q <- function(theta, prep, nx, ny, tau, xy_rho = 0) {
  q_value <- 0
  n_snp <- nrow(prep$beta_x)
  
  for (idx in seq_len(n_snp)) {
    beta_x_vec <- as.numeric(prep$beta_x[idx, ])
    se_x_vec <- as.numeric(prep$se_x[idx, ])
    beta_y_vec <- as.numeric(prep$beta_y[idx, ])
    se_y_vec <- as.numeric(prep$se_y[idx, ])
    
    y_stack <- stack_outcome_values(beta_y_vec, nx)
    x_stack <- stack_exposure_values(beta_x_vec, ny)
    alpha_stack <- stack_outcome_values(prep$alpha_hat, nx)
    residual <- matrix(y_stack - alpha_stack - theta * x_stack, ncol = 1)
    
    Sigma <- association_residual_cov(
      theta = theta,
      beta_x = beta_x_vec,
      se_x = se_x_vec,
      se_y = se_y_vec,
      rho_y = prep$rho_y,
      rho_x = prep$rho_x,
      nx = nx,
      ny = ny,
      tau2 = tau^2,
      V_alpha = prep$V_alpha,
      xy_rho = xy_rho
    )
    
    q_j <- as.numeric(t(residual) %*% ginv_solve(Sigma, residual))
    if (!is.finite(q_j)) {
      return(Inf)
    }
    q_value <- q_value + q_j
  }
  
  q_value
}

jointmr_profile_q <- function(theta, prep, nx, ny, tau, xy_rho = 0) {
  q_value <- 0
  n_snp <- nrow(prep$beta_x)
  
  for (idx in seq_len(n_snp)) {
    beta_x_vec <- as.numeric(prep$beta_x[idx, ])
    se_x_vec <- as.numeric(prep$se_x[idx, ])
    beta_y_vec <- as.numeric(prep$beta_y[idx, ] - prep$alpha_hat)
    se_y_vec <- as.numeric(prep$se_y[idx, ])
    
    cov_x <- prep$rho_x * outer(se_x_vec, se_x_vec)
    cov_y <- prep$rho_y * outer(se_y_vec, se_y_vec)
    if (all(is.finite(prep$V_alpha)) &&
        max(abs(prep$V_alpha)) > .Machine$double.eps) {
      cov_y <- cov_y + prep$V_alpha
    }
    cov_xy <- matrix(xy_rho, nrow = nx, ncol = ny) *
      outer(se_x_vec, se_y_vec)
    
    Sigma <- rbind(
      cbind(cov_x, cov_xy),
      cbind(t(cov_xy), cov_y)
    )
    Sigma <- make_psd(Sigma)
    
    observed <- matrix(c(beta_x_vec, beta_y_vec), ncol = 1)
    design <- matrix(c(rep(1, nx), rep(theta, ny)), ncol = 1)
    denom <- as.numeric(t(design) %*% ginv_solve(Sigma, design))
    if (!is.finite(denom) || denom <= 0) {
      return(Inf)
    }
    gamma_hat <- as.numeric(t(design) %*% ginv_solve(Sigma, observed) / denom)
    
    if (tau > 0) {
      cov_y_re <- cov_y + tau^2 * gamma_hat^2 * matrix(1, nrow = ny, ncol = ny)
      Sigma <- rbind(
        cbind(cov_x, cov_xy),
        cbind(t(cov_xy), cov_y_re)
      )
      Sigma <- make_psd(Sigma)
      denom <- as.numeric(t(design) %*% ginv_solve(Sigma, design))
      if (!is.finite(denom) || denom <= 0) {
        return(Inf)
      }
      gamma_hat <- as.numeric(
        t(design) %*% ginv_solve(Sigma, observed) / denom
      )
    }
    
    residual <- observed - design * gamma_hat
    q_j <- as.numeric(t(residual) %*% ginv_solve(Sigma, residual))
    if (!is.finite(q_j)) {
      return(Inf)
    }
    q_value <- q_value + q_j
  }
  
  q_value
}

estimate_jointmr_relaxed_nome <- function(fdata, nx, ny, tau, valid_indices,
                                          adjust_pleiotropy = TRUE,
                                          egger_indices = NULL,
                                          exposure_rho = 0,
                                          xy_rho = 0,
                                          bootstrap_B = 0) {
  prep <- prepare_jointmr_assoc_data(
    fdata = fdata,
    nx = nx,
    ny = ny,
    valid_indices = valid_indices,
    adjust_pleiotropy = adjust_pleiotropy,
    egger_indices = egger_indices,
    exposure_rho = exposure_rho
  )
  
  theta_start_values <- c()
  ratio_values <- c()
  for (m in seq_len(ny)) {
    for (n in seq_len(nx)) {
      bx <- prep$beta_x[, n]
      by <- prep$beta_y[, m] - prep$alpha_hat[m]
      keep <- is.finite(bx) & is.finite(by) & abs(bx) > .Machine$double.eps
      if (any(keep)) {
        theta_start_values <- c(
          theta_start_values,
          sum(bx[keep] * by[keep]) / sum(bx[keep]^2)
        )
        ratio_values <- c(ratio_values, by[keep] / bx[keep])
      }
    }
  }
  
  theta_start <- mean(theta_start_values, na.rm = TRUE)
  if (!is.finite(theta_start)) {
    theta_start <- 0
  }
  
  ratio_values <- ratio_values[is.finite(ratio_values)]
  if (length(ratio_values) >= 5) {
    ratio_range <- as.numeric(quantile(ratio_values, c(0.01, 0.99),
                                       na.rm = TRUE))
    lower <- min(-1, theta_start - 1, ratio_range[1] - 0.5)
    upper <- max(1, theta_start + 1, ratio_range[2] + 0.5)
  } else {
    lower <- theta_start - 1
    upper <- theta_start + 1
  }
  
  q_fun <- function(theta) {
    jointmr_profile_q(theta, prep, nx, ny, tau, xy_rho = xy_rho)
  }
  
  grid <- seq(lower, upper, length.out = 201)
  grid_q <- vapply(grid, q_fun, numeric(1))
  best <- which.min(grid_q)
  if (length(best) == 0 || !is.finite(grid_q[best])) {
    return(list(
      est = c(theta_hat = NA_real_, theta_se = NA_real_,
              theta_p_value = NA_real_),
      inputs = build_relaxed_nome_inputs(
        fdata = fdata, nx = nx, ny = ny, tau = tau,
        valid_indices = valid_indices,
        adjust_pleiotropy = adjust_pleiotropy,
        egger_indices = egger_indices,
        exposure_rho = exposure_rho,
        xy_rho = xy_rho,
        theta = theta_start
      )
    ))
  }
  
  if (best == 1 || best == length(grid)) {
    opt_interval <- c(lower, upper)
  } else {
    opt_interval <- c(grid[best - 1], grid[best + 1])
  }
  
  opt <- optimize(q_fun, interval = opt_interval)
  theta_hat <- opt$minimum
  q0 <- opt$objective
  
  h <- max(1e-4, abs(theta_hat) * 1e-4)
  q_minus <- q_fun(theta_hat - h)
  q_plus <- q_fun(theta_hat + h)
  hessian <- (q_plus - 2 * q0 + q_minus) / h^2
  if (!is.finite(hessian) || hessian <= 0) {
    h <- h * 10
    q_minus <- q_fun(theta_hat - h)
    q_plus <- q_fun(theta_hat + h)
    hessian <- (q_plus - 2 * q0 + q_minus) / h^2
  }
  
  theta_se <- if (is.finite(hessian) && hessian > 0) {
    sqrt(2 / hessian)
  } else {
    NA_real_
  }
  theta_p_value <- if (is.finite(theta_se) && theta_se > 0) {
    2 * pnorm(-abs(theta_hat / theta_se))
  } else {
    NA_real_
  }
  
  inputs <- build_relaxed_nome_inputs(
    fdata = fdata,
    nx = nx,
    ny = ny,
    tau = tau,
    valid_indices = valid_indices,
    adjust_pleiotropy = adjust_pleiotropy,
    egger_indices = egger_indices,
    exposure_rho = exposure_rho,
    xy_rho = xy_rho,
    theta = theta_hat
  )
  
  est <- c(theta_hat = theta_hat, theta_se = theta_se,
           theta_p_value = theta_p_value)
  
  if (bootstrap_B > 0) {
    boot <- bootstrap_jointmr_relaxed_nome(
      fdata = fdata,
      nx = nx,
      ny = ny,
      tau = tau,
      valid_indices = valid_indices,
      adjust_pleiotropy = adjust_pleiotropy,
      egger_indices = egger_indices,
      exposure_rho = exposure_rho,
      xy_rho = xy_rho,
      B = bootstrap_B
    )
    if (is.finite(boot[["theta_se_boot"]]) && boot[["theta_se_boot"]] > 0) {
      boot <- c(
        boot,
        theta_p_value_boot = 2 * pnorm(-abs(theta_hat / boot[["theta_se_boot"]]))
      )
    } else {
      boot <- c(boot, theta_p_value_boot = NA_real_)
    }
    est <- c(est, boot)
  }
  
  list(est = est, inputs = inputs)
}

bootstrap_jointmr_relaxed_nome <- function(fdata, nx, ny, tau, valid_indices,
                                           adjust_pleiotropy = TRUE,
                                           egger_indices = NULL,
                                           exposure_rho = 0,
                                           xy_rho = 0,
                                           B = 200) {
  boot_theta <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    boot_indices <- sample(valid_indices, length(valid_indices), replace = TRUE)
    boot_fit <- estimate_jointmr_relaxed_nome(
      fdata = fdata,
      nx = nx,
      ny = ny,
      tau = tau,
      valid_indices = boot_indices,
      adjust_pleiotropy = adjust_pleiotropy,
      egger_indices = egger_indices,
      exposure_rho = exposure_rho,
      xy_rho = xy_rho,
      bootstrap_B = 0
    )
    boot_theta[b] <- unname(boot_fit$est[["theta_hat"]])
  }
  
  boot_se <- sd(boot_theta, na.rm = TRUE)
  c(theta_se_boot = boot_se)
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
                     exposure_rho = 0, xy_rho = 0,
                     seed_base = 20260428) {
  set.seed(seed_base + i)
  
  fdata <- DataGeneration_summary(g, nx, ny, bx, beta0, tau, rho_1, rho_2,
                                  siga1, siga2, weakp,
                                  exposure_rho = exposure_rho,
                                  xy_rho = xy_rho)
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
    exposure_rho = exposure_rho,
    xy_rho = xy_rho,
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
      JointMR_RN_beta = NA_real_,
      JointMR_RN_se = NA_real_,
      JointMR_RN_pval = NA_real_,
      JointMR_RN_ple_plugin_beta = NA_real_,
      JointMR_RN_ple_plugin_se = NA_real_,
      JointMR_RN_ple_plugin_pval = NA_real_,
      alpha0_mean = NA_real_,
      alpha0_true_mean = mean(fdata$alpha0_y)
    )))
  }
  
  rn_fit <- estimate_jointmr_relaxed_nome(
    fdata = fdata,
    nx = nx,
    ny = ny,
    tau = tau,
    valid_indices = valid_indices,
    adjust_pleiotropy = FALSE,
    egger_indices = seq_len(g),
    exposure_rho = exposure_rho,
    xy_rho = xy_rho
  )
  
  rn_est <- rn_fit$est
  
  rn_ple_fit <- estimate_jointmr_relaxed_nome(
    fdata = fdata,
    nx = nx,
    ny = ny,
    tau = tau,
    valid_indices = valid_indices,
    adjust_pleiotropy = TRUE,
    egger_indices = seq_len(g),
    exposure_rho = exposure_rho,
    xy_rho = xy_rho
  )
  
  rn_ple_est <- rn_ple_fit$est
  
  as_result_row(c(
    base_values,
    comparison_values,
    JointMR_RN_beta = unname(rn_est[["theta_hat"]]),
    JointMR_RN_se = unname(rn_est[["theta_se"]]),
    JointMR_RN_pval = unname(rn_est[["theta_p_value"]]),
    JointMR_RN_ple_plugin_beta = unname(rn_ple_est[["theta_hat"]]),
    JointMR_RN_ple_plugin_se = unname(rn_ple_est[["theta_se"]]),
    JointMR_RN_ple_plugin_pval = unname(rn_ple_est[["theta_p_value"]]),
    alpha0_mean = mean(rn_ple_fit$inputs$alpha_hat),
    alpha0_true_mean = mean(fdata$alpha0_y)
  ))
}

# Optional bootstrap SE/P-value code.
# The main analysis uses the model-based SE from gls_theta_plugin().
# To use bootstrap SE for sensitivity checks, uncomment this function and call it
# inside run_once() after plugin_inputs is created.
#
# bootstrap_theta <- function(fdata, nx, ny, tau, valid_indices, theta_hat, B = 1000) {
#   boot_theta <- rep(NA_real_, B)
#   for (b in seq_len(B)) {
#     boot_indices <- sample(valid_indices, length(valid_indices), replace = TRUE)
#     boot_inputs <- build_plugin_adjusted_inputs(
#       fdata = fdata,
#       nx = nx,
#       ny = ny,
#       tau = tau,
#       valid_indices = boot_indices,
#       egger_indices = boot_indices
#     )
#     boot_est <- gls_theta_plugin(
#       base_cov_list = boot_inputs$base_cov_list,
#       WR_matrix = boot_inputs$WR_adj,
#       BA_list = boot_inputs$BA_list,
#       V_alpha = boot_inputs$V_alpha,
#       tau2 = boot_inputs$tau2
#     )
#     boot_theta[b] <- unname(boot_est[["theta_hat"]])
#   }
#   boot_se <- sd(boot_theta, na.rm = TRUE)
#   boot_p_value <- 2 * pnorm(-abs(theta_hat / boot_se))
#   c(theta_se_boot = boot_se, theta_p_value_boot = boot_p_value)
# }

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
    "nx", "ny", "g", "beta0", "rho_mean", "exposure_rho", "xy_rho", "tau", "siga1",
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
    paste0("JointMR_RN_", c("beta", "se", "pval")),
    paste0("JointMR_RN_ple_plugin_", c("beta", "se", "pval")),
    "F_tat", "weakp",
    "alpha0_hat_mean", "alpha0_se_mean", "alpha0_true_mean",
    paste0("alpha0_hat_Y", seq_len(ny)),
    paste0("alpha0_se_Y", seq_len(ny)),
    paste0("alpha0_true_Y", seq_len(ny))
  )
}

#####Comp funcrion######
Comp <- function(i, g, nx, ny, bx, beta0, tau, rho_1, rho_2, siga1, siga2,
                 weakp, F_THRESHOLD, P_THRESHOLD,
                 exposure_rho = 0, xy_rho = 0) {
  cat("i=", i, "\n")
  set.seed(20260428 + i)
  
  fdata1 <- DataGeneration_summary(g, nx, ny, bx, beta0, tau, rho_1, rho_2,
                                   siga1, siga2, weakp,
                                   exposure_rho = exposure_rho,
                                   xy_rho = xy_rho)
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
  
  rn_theta_hat <- NA_real_
  rn_theta_se <- NA_real_
  rn_theta_p_value <- NA_real_
  rn_ple_theta_hat <- NA_real_
  rn_ple_theta_se <- NA_real_
  rn_ple_theta_p_value <- NA_real_
  alpha0_hat <- rep(NA_real_, ny)
  alpha0_se <- rep(NA_real_, ny)
  
  if (SNPnew >= 2) {
    rn_fit <- estimate_jointmr_relaxed_nome(
      fdata = fdata1,
      nx = nx,
      ny = ny,
      tau = tau,
      valid_indices = valid_indices_intersect,
      adjust_pleiotropy = FALSE,
      egger_indices = seq_len(g),
      exposure_rho = exposure_rho,
      xy_rho = xy_rho
    )
    
    rn_est <- rn_fit$est
    
    rn_ple_fit <- estimate_jointmr_relaxed_nome(
      fdata = fdata1,
      nx = nx,
      ny = ny,
      tau = tau,
      valid_indices = valid_indices_intersect,
      adjust_pleiotropy = TRUE,
      egger_indices = seq_len(g),
      exposure_rho = exposure_rho,
      xy_rho = xy_rho
    )
    
    rn_ple_est <- rn_ple_fit$est
    
    rn_theta_hat <- unname(rn_est[["theta_hat"]])
    rn_theta_se <- unname(rn_est[["theta_se"]])
    rn_theta_p_value <- unname(rn_est[["theta_p_value"]])
    rn_ple_theta_hat <- unname(rn_ple_est[["theta_hat"]])
    rn_ple_theta_se <- unname(rn_ple_est[["theta_se"]])
    rn_ple_theta_p_value <- unname(rn_ple_est[["theta_p_value"]])
    alpha0_hat <- rn_ple_fit$inputs$alpha_hat
    alpha0_se <- rn_ple_fit$inputs$alpha_se
  }
  
  res_all <- c(
    nx, ny, g, beta0, rho_1, exposure_rho, xy_rho, tau, siga1, siga2,
    comparison$SNP_GWAS_meta, comparison$SNPMR, SNPnew,
    comparison$gwas_meta_res,
    comparison$pairwise_res,
    comparison$mr_meta_res,
    rn_theta_hat, rn_theta_se, rn_theta_p_value,
    rn_ple_theta_hat, rn_ple_theta_se, rn_ple_theta_p_value,
    F_tat, weakp,
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
  
  write_csv(res_all,file =paste0("Sum_nx_", nx, "_ny_", ny, "_g_", g,
   "_beta0_", beta0, "_rho1_", rho_1,
    "_XXrho_", exposure_rho, "_XYrho_", xy_rho,
    "_tau_", tau, "_siga1_", siga1, "_siga2_", siga2,
    "_weakp_", weakp, "_F_THRESHOLD_", F_THRESHOLD,
    "_P_THRESHOLD_", P_THRESHOLD,
    "_0428_relaxed_nome_ple_plugin_delta.csv"), append = TRUE)
}

####mian function####
main <- function(NN, g, nx, ny, bx, beta0, tau, rho_1, rho_2, siga1, siga2,
                 weakp, F_THRESHOLD, P_THRESHOLD, mc,
                 exposure_rho, xy_rho) {
                  
  registerDoMC(mc)
  foreach(i=1:NN,
          .combine=rbind,
          .errorhandling = "remove",
          .packages = c("MASS", "meta", "MendelianRandomization", "mr.divw", "readr")) %dopar% {
            Comp(i, g, nx, ny, bx, beta0, tau, rho_1,
                 rho_2, siga1, siga2,
                 weakp, F_THRESHOLD, P_THRESHOLD,
                 exposure_rho = exposure_rho,
                 xy_rho = xy_rho)
          }
}

F_THRESHOLDb <- c(10)
P_THRESHOLD <- 5e-8
NN <- 1000
mc <- 50
bx <- 0.2
gb <- c(200)
beta0b <- c(0, 0.05)
taub <- c(0.2, 0)
rho_1b <- c(0.9, 0.5)
rho_2b <- c(0.9, 0.5)
siga1b <- c(0, 0)
siga2b <- c(0, 0.05)
nxb <- c(3)
nyb <- c(5)
weakpb <- c(0)
exposure_rhob <- c(0.9, 0.5)
xy_rhob <- c(0, 0.05)

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
                  for (xxh in seq_along(exposure_rhob)) {
                    exposure_rho <- exposure_rhob[[xxh]]
                    for (xyh in seq_along(xy_rhob)) {
                      xy_rho <- xy_rhob[[xyh]]
                      
                      ccname <- build_ccname(nx, ny)
                      
                      main(NN, g, nx, ny, bx, beta0, tau, rho_1, rho_2,
                           siga1, siga2, weakp, F_THRESHOLD, P_THRESHOLD, mc,
                           exposure_rho = exposure_rho,
                           xy_rho = xy_rho)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
