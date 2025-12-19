setwd("/data/wusijia/simulation/weak")

library(foreach)
library(readr)
library(doParallel)
library(doMC)
library(MASS)
library(MendelianRandomization)
library(meta)
library(TwoSampleMR)
library(mr.divw)
# library(Rmisc)
library(dplyr)
library(purrr) 
library(stats4)

####mle#####
MLE_S <- function(nx,ny,cov_list,WR_matrix,SNPnew){
  n <- nx * ny
  precompute <- lapply(cov_list, function(Sigma) {
    list(
      Sigma_inv = ginv(Sigma),
      log_det = det(Sigma)
    )
  })
  
  I <- matrix(1, n, 1)
  
  neg_log_lik <- function(theta) {
    total <- 0
    for (j in 1:SNPnew) {
      residual <- WR_matrix[,j] - theta * I[,1]
      quad_form <- as.numeric(crossprod(residual, precompute[[j]]$Sigma_inv) %*% residual)
      total <- total + 0.5*(quad_form + precompute[[j]]$log_det)
    }
    total + 0.5*n*SNPnew*log(2*pi)  # 添加常数项
  }
  
  
  fit <- mle(minuslogl = neg_log_lik,
             start = list(theta = 0))
  
  return(fit@coef)
}

##########data generation##########
DataGeneration <- function(g,nx,ny,bx,beta0,tau,rho_1,rho_2,siga1,siga2,weakp){
  
  N1 <- 5000
  N2 <- 50000
  #####exposure  
  weak_count <- max(1, round(g * weakp))
  strong_count <- g - weak_count

  strong_instruments_base <- rnorm(strong_count, bx, 0.03) # 基础效应
  weak_instruments_base <- rnorm(weak_count, 0, 0.005)
  
  all_instruments_base <- c(strong_instruments_base, weak_instruments_base)
  shared_instrument_order <- sample(all_instruments_base) 
  
  betaXG_initial <- list()
  for (i in 1:nx) {
    betaXG_initial[[i]] <- shared_instrument_order + rnorm(g, 0, 0.001) 
  }
  
  ###小样本
  Gx <- list()
  for(kk in 1:nx){
    G1 <- NULL
    for(i in 1:g){
    Gg <- rbinom(N1,2,0.2)
    G1 <- cbind(G1,Gg)
    }
    Gx[[kk]] <- G1 
  }
  
  X_R <- list()
  for (i in 1:nx) {
    X_R[[i]] <- Gx[[i]] %*% betaXG_initial[[i]] + rnorm(N1,0,0.1)
  }
  
  seXG <- list()
  betaXG <- list()
  for (i in 1:nx) {
    seX1G <- NULL
    betaX1G <- NULL
    for(j in 1:g){
      FF1 <- lm(X_R[[i]]~Gx[[i]][,j])
      seX1G <- c(seX1G,summary(FF1)$coef[2,2])
      betaX1G <- c(betaX1G,summary(FF1)$coef[2,1])
    }
    seXG[[i]] <- seX1G
    betaXG[[i]] <- betaX1G
  }
  
  ###大样本
  Gx2 <- list()
  for(kk in 1:nx){
    G1 <- NULL
    for(i in 1:g){
      Gg <- rbinom(N2,2,0.2)
      G1 <- cbind(G1,Gg)
    }
    Gx2[[kk]] <- G1 
  }
  
  X_R2 <- list()
  for (i in 1:nx) {
    X_R2[[i]] <- Gx2[[i]] %*% betaXG_initial[[i]] + rnorm(N2,0,0.1)
  }
  
  PXG <- lapply(1:nx, function(i) {
    sapply(1:g, function(j) {
      model <- lm(X_R2[[i]] ~ Gx2[[i]][, j])
      summary(model)$coefficients[2, 4] # 直接提取P值
    })
  })
  
  
  # PXG <- list()
  # for (i in 1:nx) {
  #   PX1G <- NULL
  #   for(j in 1:g){
  #     FF2 <- lm(X_R2[[i]]~Gx2[[i]][,j])
  #     PX1G <- c(PX1G,summary(FF2)$coef[2,4])
  #   }
  #   PXG[[i]] <- PX1G
  # }
  
  ####outcome se
  seYG <- list()
  for (i in 1:ny) {
    seYG[[i]] <- runif(g,0.05,0.15)##0.01 0.13
  } 
  
  
  gamma <- list()
  for (j in 1:(nx*ny)) {
    gamma[[j]] <- runif(g,siga1,siga2)
  }
  
  ####wald ratio
  sewald <- list()
  for (i in 1:ny) {
    for (j in 1:nx) {
      sewald[[(i-1)*nx+j]] <- sqrt((seYG[[i]]^2)/(betaXG[[j]]^2))
    }
  }
  
  cor_matrix <- diag(1, nrow = nx * ny, ncol = nx * ny)
  
  # 遍历所有 X 和 Y 数据库，填充相关性
  for (i in 1:ny) {  
    for (j in 1:ny) { 
      # 获取当前 X 数据库的起始和结束行/列
      start_row_y <- (i - 1) * nx + 1
      end_row_y <- i * nx
      start_col_y <- (j - 1) * nx + 1
      end_col_y <- j * nx
      
      # 填充相关性矩阵
      if (i == j) {
        # 对角线块填充为全 1 矩阵
        cor_matrix[start_row_y:end_row_y, start_col_y:end_col_y] <- matrix(1, nx, nx)
      } else {
        # 非对角线块填充相关性
        rho_c <- runif(1,min = rho_1, max = rho_2)
        cor_matrix[start_row_y:end_row_y, start_col_y:end_col_y] <- rho_c
        cor_matrix[start_col_y:end_col_y, start_row_y:end_row_y] <- rho_c
      }
    }
  }
  
  waldratio <- list()
  waldratio_matrix <- NULL 
  beta_m <- rep(beta0, nx*ny)
  for(k in 1:g){
    # 初始化协方差矩阵
    cov_matrix <- matrix(0, nrow = nx*ny, ncol = nx*ny)
    
    # 填充协方差矩阵
    for (s in 1:(nx*ny)) {
      for (p in 1:(nx*ny)) {
        if (s == p) {
          # 对角线元素
          cov_matrix[s, p] <- sewald[[s]][k]^2 + tau^2
        } else {
          # 非对角线元素
          cov_matrix[s, p] <- sewald[[s]][k] * sewald[[p]][k] * cor_matrix[s,p]+ tau^2
        }
      }
    }
    
    wald_r <- mvrnorm(1,beta_m,cov_matrix)
    waldratio_matrix <- rbind(waldratio_matrix,wald_r)
    
  }
  
  # 将生成的wald存储到waldratio中
  for (n in 1:(nx*ny)) {
    waldratio[[n]] <- as.vector(waldratio_matrix[,n])
  }
  
  
  ####生成outcome beta
  betaYG <- list()
  # 将生成的betaY存储到betaYG中
  for (n in 1:ny) {
    betaY_m <- NULL
    
    for (m in 1:nx) {
      
      betaY_m1 <- (waldratio[[(n-1)*nx+m]] * betaXG[[m]]) + gamma[[(n-1)*nx+m]]
      betaY_m <- rbind(betaY_m,betaY_m1)
    }
    betaYG[[n]] <- colMeans(betaY_m)
  }
  
  #### 计算F统计量
  F_stats <- list()
  F_stats_summary <- list()
  
  for (i in 1:nx) {
    # 方法1：计算每个SNP的个体F统计量 [8](@ref)
    individual_F <- (betaXG[[i]] / seXG[[i]])^2
    
    # 方法2：计算第一阶段回归的F统计量（模拟）[6](@ref)
    # 由于我们没有原始数据，这里使用近似计算
    mean_beta <- mean(betaXG[[i]])
    mean_se <- mean(seXG[[i]])
    approximate_F <- (mean_beta / mean_se)^2
    
    # 方法3：基于R²的F统计量计算 [8](@ref)
    R_squared <- (betaXG[[i]]^2) / (betaXG[[i]]^2 + seXG[[i]]^2 * (g - 1))
    F_from_R2 <- (g - 2) * R_squared / (1 - R_squared)
    
    F_stats[[i]] <- list(
      individual_F = individual_F,
      mean_individual_F = mean(individual_F),
      approximate_F = approximate_F,
      F_from_R2 = F_from_R2,
      weak_instrument_percentage = mean(individual_F < 10) * 100  # F < 10 视为弱工具变量
    )
    
    # 生成汇总统计
    F_stats_summary[[i]] <- data.frame(
      Exposure = i,
      Mean_F_individual = mean(individual_F),
      Min_F_individual = min(individual_F),
      Max_F_individual = max(individual_F),
      Approximate_F = approximate_F,
      Weak_Instruments_Pct = mean(individual_F < 10) * 100,
      Strong_Instruments_Pct = mean(individual_F >= 10) * 100
    )
  }
  
  data_sum <- list(
    betaYG = betaYG,
    betaXG = betaXG,
    seXG = seXG,
    seYG = seYG,
    PXG = PXG,
    gamma = gamma,
    sewald = sewald,
    waldratio = waldratio,
    rho_matrix = cor_matrix,
    F_statistics = F_stats,
    F_summary = do.call(rbind, F_stats_summary)
  )
  
  return(data_sum)
}


#####Comp funcrion######
Comp <- function(i,g,nx,ny,bx,beta0,tau,rho_1,rho_2,siga1,siga2,
                 weakp, F_THRESHOLD,P_THRESHOLD){
  
  cat('i=',i,'\n')
  
  fdata1 <- DataGeneration(g,nx,ny,bx,beta0,tau,rho_1,rho_2,siga1,siga2,weakp)
  F_tat <- mean(fdata1$F_summary$Mean_F_individual)
  
  f_stat_matrix <- do.call(cbind, lapply(fdata1$F_statistics, function(x) x$individual_F))
  pval_stat_matrix <- do.call(cbind, lapply(fdata1$PXG, matrix))
  data_x_beta <- do.call(cbind, lapply(fdata1$betaXG, matrix))
  data_x_se <- do.call(cbind, lapply(fdata1$seXG, matrix))
  data_y_beta <- do.call(cbind, lapply(fdata1$betaYG, matrix))
  data_y_se <- do.call(cbind, lapply(fdata1$seYG, matrix))
  
  ####GWAS-meta-MR####
  alpha_data <- NULL
  theta_data <- NULL
  for (j in 1:g) {
    
    retx <- try(metagen(data_x_beta[j,],data_x_se[j,],sm="MD"), silent=TRUE)
    rety <- try(metagen(data_y_beta[j,],data_y_se[j,],sm="MD"), silent=TRUE)
    if ('try-error' %in% class(retx) | 'try-error' %in% class(rety)) {
      next
    }else{
      if(retx$pval.Q<0.05){
        alpha_data <- rbind(alpha_data,c(retx$TE.random,retx$seTE.random,retx$pval.random))
      }else{
        alpha_data <- rbind(alpha_data,c(retx$TE.fixed,retx$seTE.fixed,retx$pval.fixed))
      }
      
      if(rety$pval.Q<0.05){
        theta_data <- rbind(theta_data,c(rety$TE.random,rety$seTE.random,rety$pval.random))
      }else{
        theta_data <- rbind(theta_data,c(rety$TE.fixed,rety$seTE.fixed,rety$pval.fixed))
      }
    }
  }
  
  # --- 筛选策略 1: 仅使用 meta-F-stat > 10 的SNP ---
  F_stat_meta <- (alpha_data[,1] / alpha_data[,2])^2
  P_valmeta <- alpha_data[,3]
  valid_indices_meta <- which(F_stat_meta > F_THRESHOLD & P_valmeta < P_THRESHOLD)
  
  alpha_data1 <- data.frame(alpha_data[valid_indices_meta,])
  theta_data1 <- data.frame(theta_data[valid_indices_meta,])
  SNP_GWAS_meta <- length(alpha_data1[,1])
  
  ####MR
  mr_object_fix <- mr_input(bx = alpha_data1[,1],
                            bxse = alpha_data1[,2],
                            by = theta_data1[,1],
                            byse = theta_data1[,2])##fix
  
  ###mr.divw
  r_divw <- mr.divw(alpha_data1[,1],theta_data1[,1],alpha_data1[,2],theta_data1[,2])
  fix_res_divw <-  c(fix_beta_divw=r_divw$beta.hat,
                     fix_se_divw=r_divw$beta.se,
                     fix_pval_divw=ci(r_divw$beta.hat, r_divw$beta.se)$p)
  
  ##fix-ivw
  r_IVW_fixed <- MendelianRandomization::mr_ivw(mr_object_fix,model = "fixed")
  fix_res_IVW_fixed <- c(fix_beta_ivw_fixed=r_IVW_fixed@Estimate,
                         fix_se_ivw_fixed=r_IVW_fixed@StdError,
                         fix_pval_ivw_fixed=r_IVW_fixed@Pvalue)
  
  
  ##random-ivw
  r_IVW_random <- MendelianRandomization::mr_ivw(mr_object_fix,model = "random")
  fix_res_IVW_random <- c(fix_beta_ivw_random=r_IVW_random@Estimate,
                          fix_se_ivw_random=r_IVW_random@StdError,
                          fix_pval_ivw_random=r_IVW_random@Pvalue)
  
  
  #weighted median
  r_wme <- MendelianRandomization::mr_median(mr_object_fix,iterations = 20)
  fix_res_wme <- c(fix_beta_wme = r_wme@Estimate,
                   fix_se_wme = r_wme@StdError,
                   fix_pval_wme = r_wme@Pvalue)
  
  rere1 <- c(fix_res_divw,fix_res_IVW_fixed,fix_res_IVW_random,fix_res_wme)
  
  
  ####传统MR####
  mr_object <- list()
  res_divw <- NULL
  res_fixed <- NULL
  res_random <- NULL
  res_wme <- NULL
  SNPSMR <- NULL
  for (m in 1:nx){
    
    # --- 筛选策略 2: 仅使用对当前暴露 m, F-stat > 10 的SNP ---
    valid_indices_trad <- which(f_stat_matrix[, m] > F_THRESHOLD & pval_stat_matrix [,m] < P_THRESHOLD)
    SNPSMR <- c(SNPSMR, length(valid_indices_trad)) # 记录使用的SNP数量
    
    # --- 检查是否有足够的SNP进行分析 ---
    if (length(valid_indices_trad) < 2) {
      for (n in 1:ny) {
        # 如果SNP不够，填充NA
        res_divw <- rbind(res_divw, c(beta_divw=NA, se_divw=NA, pval_divw=NA))
        res_fixed <- rbind(res_fixed, c(beta_ivwf=NA, se_ivwf=NA, pval_ivwf=NA))
        res_random <- rbind(res_random, c(beta_ivwr=NA, se_ivwr=NA, pval_ivwr=NA))
        res_wme <- rbind(res_wme, c(beta_wme=NA, se_wme=NA, pval_wme=NA))
      }
      next # 跳过这个暴露
    }
    
    # (如果SNP足够，则正常运行)
    for (n in 1:ny) {
      mr_object[[(m-1)*nx+n]] <- mr_input(bx=data_x_beta[valid_indices_trad,m], 
                                          bxse=data_x_se[valid_indices_trad,m], 
                                          by=data_y_beta[valid_indices_trad,n], 
                                          byse=data_y_se[valid_indices_trad,n]) 
      
      ###mr.divw
      res1_divw <- try(mr.divw(data_x_beta[valid_indices_trad,m], data_y_beta[valid_indices_trad,n],
                               data_x_se[valid_indices_trad,m], data_y_se[valid_indices_trad,n]), silent=TRUE)
      res_divw <-  rbind(res_divw, if ('try-error' %in% class(res1_divw)) rep(NA, 3) else c(res1_divw$beta.hat, res1_divw$beta.se, ci(res1_divw$beta.hat, res1_divw$beta.se)$p))
      
      ###fix-IVW
      res1_fixed <- try(MendelianRandomization::mr_ivw(mr_object[[(m-1)*nx+n]], model = "fixed"), silent=TRUE)
      res_fixed <- rbind(res_fixed, if ('try-error' %in% class(res1_fixed)) rep(NA, 3) else c(res1_fixed@Estimate, res1_fixed@StdError, res1_fixed@Pvalue))
      
      ###random IVW
      res1_random <- try(MendelianRandomization::mr_ivw(mr_object[[(m-1)*nx+n]], model = "random"), silent=TRUE)
      res_random <- rbind(res_random, if ('try-error' %in% class(res1_random)) rep(NA, 3) else c(res1_random@Estimate, res1_random@StdError, res1_random@Pvalue))
      
      #weighted median
      if (length(valid_indices_trad) < 3) {
        res_wme <- rbind(res_wme,c(beta_wme = NA, se_wme = NA, pval_wme = NA))
      }else{
        res1_wme <- try(MendelianRandomization::mr_median(mr_object[[(m-1)*nx+n]],iterations = 20), silent=TRUE)
        res_wme <- rbind(res_wme, if ('try-error' %in% class(res1_wme)) rep(NA, 3) else c(res1_wme@Estimate, res1_wme@StdError, res1_wme@Pvalue))
      }
    } # end 'n' loop
    
    rere2 <- c(res_divw[,1],res_fixed[,1],res_random[,1],res_wme[,1],
               res_divw[,2],res_fixed[,2],res_random[,2],res_wme[,2],
               res_divw[,3],res_fixed[,3],res_random[,3],res_wme[,3])
  } 
  
  SNPMR <- round(mean(SNPSMR, na.rm=TRUE), 0) # 计算平均SNP数
  
  # --- Meta-analyze 传统MR 结果 ---
  
  run_meta <- function(res_df) {
    res_clean <- na.omit(res_df)
    if (nrow(res_clean) == 0) return(rep(NA, 3))
    
    meta_res <- try(metagen(res_clean[,1], res_clean[,2], sm="MD"), silent=TRUE)
    if ('try-error' %in% class(meta_res)) return(rep(NA, 3))
    
    if(meta_res$pval.Q < 0.05){
      return(c(meta_res$TE.random, meta_res$seTE.random, meta_res$pval.random))
    } else {
      return(c(meta_res$TE.fixed, meta_res$seTE.fixed, meta_res$pval.fixed))
    }
  }
  
  meta_old_divw1 <- run_meta(res_divw)
  meta_old_ivw_fixed1 <- run_meta(res_fixed)
  meta_old_ivw_random1 <- run_meta(res_random)
  meta_old_wme1 <- run_meta(res_wme)
  
  rere3 <- c(meta_old_divw1,meta_old_ivw_fixed1,meta_old_ivw_random1,meta_old_wme1)

  
  ####MR-meta####
  # --- 筛选策略 3: 仅使用对 *所有* 暴露, F-stat > 10 的SNP ("交集") ---
  valid_indices_intersect <- which(
    apply(f_stat_matrix, 1, function(row) all(row >= F_THRESHOLD)) &  # F统计量条件
      apply(pval_stat_matrix, 1, function(row) any(row < P_THRESHOLD))   # P值条件
  )
  SNPnew <- length(valid_indices_intersect) # 这是筛选后的SNP数量
  
  # --- 检查是否有足够的SNP进行分析 ---
  if (SNPnew < 2) {
    theta_hat <- NA
    theta_se <- NA
    theta_p_value <- NA
  } else {
    
    # --- 准备数据 (仅使用筛选后的SNP) ---
    WR_matrix_full <- do.call(rbind,fdata1$waldratio)
    seWR_matrix_full <-  do.call(rbind,fdata1$sewald)
    
    colname <- NULL
    for (n in 1:ny) {
      for (m in 1:nx) {
        colname1 <- paste0("X", m, "-Y", n)
        colname <- c(colname,colname1)
      }
    }
    rowname_full <- NULL
    for (j in 1:g) {
      rowname1 <- paste0("SNP", j)
      rowname_full <- c(rowname_full,rowname1)
    }
    
    rownames(WR_matrix_full) <- colname
    rownames(seWR_matrix_full) <- colname
    colnames(WR_matrix_full) <- rowname_full
    colnames(seWR_matrix_full) <- rowname_full
    
    # --- 关键：子集化 (Subset) 矩阵 ---
    WR_matrix <- WR_matrix_full[, valid_indices_intersect, drop = FALSE]
    seWR_matrix <- seWR_matrix_full[, valid_indices_intersect, drop = FALSE]
    
    #### 计算协方差矩阵 (仅为筛选后的SNP)
    cor_matrix <- fdata1$rho_matrix
    samples <- rownames(WR_matrix)
    
    cov_list <- lapply(colnames(seWR_matrix), function(snp) { # 迭代筛选后的SNP
      se_vector <- seWR_matrix[samples, snp]
      se_outer <- outer(se_vector, se_vector)
      cov_matrix <- cor_matrix * se_outer + tau^2
      dimnames(cov_matrix) <- list(samples, samples)
      return(cov_matrix)
    })
    
    ##step 6:计算因果效应估计值
    theta_hat <- MLE_S(nx,ny,cov_list,WR_matrix,SNPnew) # SNPnew 是筛选后的数量
    
    # --- Bootstrap (使用筛选后的数据) ---
    boot_effects <- NULL
    for (k in 1:1000){
      boot_sam_indices <- sample(1:SNPnew, size = SNPnew, replace=TRUE) # 从 1:SNPnew 中采样
      WR_matrix_sample <- WR_matrix[, boot_sam_indices, drop = FALSE]
      cov_list_sample <- cov_list[boot_sam_indices]
      
      # 传递正确的SNP数量 (SNPnew)，因为重采样的大小就是SNPnew
      bott_e <- MLE_S(nx,ny,cov_list_sample,WR_matrix_sample,SNPnew) 
      boot_effects <- as.numeric(c(boot_effects,bott_e))
    }
    
    theta_se <- sd(boot_effects)
    z_score <- theta_hat / theta_se
    theta_p_value <- 2 * pnorm(-abs(z_score))
  } 
  
  
  # --- 汇总所有结果 ---
  res_all <- c(nx,ny,g,beta0,rho_1,tau,siga1,siga2,
               SNP_GWAS_meta, SNPMR, SNPnew, # 这三个现在都是*筛选后*的SNP计数
               rere1,rere2,rere3,
               theta_hat,theta_se,theta_p_value,F_tat,weakp) 
  
  res_all <- matrix(res_all ,nrow=1)
  
  res_all <- as.data.frame(res_all)
  colnames(res_all) <- ccname
  
  write_csv(res_all,file =paste0('Sum_nx_',nx,'_ny_',ny,'_g_',g,
                                 '_beta0_',beta0,"_rho1_",rho_1,
                                 '_tau_',tau,"_siga1_",siga1,"_siga2_",siga2,
                                 '_weakp_',weakp,"_F_THRESHOLD_",F_THRESHOLD,
                                 "_P_THRESHOLD_",P_THRESHOLD,
                                 "_1128.csv"), append = TRUE)
}

####mian function####
main <- function(NN,g,nx,ny,bx,beta0,tau,rho_1,rho_2,siga1,siga2,weakp,F_THRESHOLD,P_THRESHOLD,mc){
  
  
  registerDoMC(mc)
  foreach(i=1:NN,
          .combine=rbind,
          .errorhandling = "remove",
          .packages = c("MASS")) %dopar% {
            Comp(i,g,nx,ny,bx,beta0,tau,rho_1,rho_2,siga1,siga2,weakp,F_THRESHOLD,P_THRESHOLD)
          }
}


F_THRESHOLDb <- c(10)
P_THRESHOLD <- 5e-8
NN=1000
mc=50
bx=0.2
gb <- c(200)
beta0b <- c(0,0.05,0.1)
taub <- c(0,0.2)
rho_1b <- c(0.5,0.8,0.85,0.9,0.95)
rho_2b <- c(0.5,0.8,0.85,0.9,0.95)
siga1b <- c(0,0)
siga2b <- c(0,0.01)
nxb <- c(3)
nyb <- c(3,5)
weakpb <- c(0,0.3,0.5)

for(zzh in 1:length(F_THRESHOLDb)){
  F_THRESHOLD = F_THRESHOLDb[zzh]
  for(yh in 1:length(taub)){
    tau=taub[[yh]]
    for(zh in 1:length(nxb)){
      nx=nxb[[zh]]
      for(qh in 1:length(nyb)){
        ny=nyb[[qh]]
        for(xh in 1:length(siga1b)){
          siga1=siga1b[[xh]]
          siga2=siga2b[[xh]]
          for(tth in 1:length(weakpb)){
            weakp = weakpb[tth]
            for (th in 1:length(gb)) {
              g= gb[th]
              for(uh in 1:length(beta0b)){
                beta0=beta0b[[uh]]
                for(rh in 1:length(rho_1b)){
                  rho_1=rho_1b[[rh]]
                  rho_2=rho_2b[[rh]]
    
                  ####names
                  methods <- c("dIVW", "IVW_fixed", "IVW_random", "wme")
                  all_beta_names <- unlist(lapply(methods, function(method) {
                    paste0("MR_", method, "_", rep("beta", nx*ny), seq_len(nx*ny))
                  }))
                  
                  all_se_names <- unlist(lapply(methods, function(method) {
                    paste0("MR_", method, "_", rep("se", nx*ny), seq_len(nx*ny))
                  }))
                  
                  all_pval_names <- unlist(lapply(methods, function(method) {
                    paste0("MR_", method, "_", rep("pval", nx*ny), seq_len(nx*ny))
                  }))
                  
                  ccname <- c('nx','ny','g','beta0',"rho_mean",'tau',"siga1",
                              "siga2","SNP_gaws_meta","SNPMR","SNPnew",
                              
                              paste0('meta_res_dIVW_',c('beta','se','pval')),
                              paste0('meta_res_IVW_fixed_',c('beta','se','pval')),
                              paste0('meta_res_IVW_random_',c('beta','se','pval')),
                              paste0('meta_res_wme_',c('beta','se','pval')),
                              
                              all_beta_names,all_se_names,all_pval_names,
                              
                              paste0('MR_meta_old_dIVW_',c('beta','se','pval')),
                              paste0('MR_meta_old_ivw_fixed_',c('beta','se','pval')),
                              paste0('MR_meta_old_ivw_random_',c('beta','se','pval')),
                              paste0('MR_meta_old_wme_',c('beta','se','pval')),
                              paste0('MRmeta_',c('beta','se','pval')),"F_tat","weakp")
                  
                  main(NN,g,nx,ny,bx,beta0,tau,rho_1,rho_2,siga1,siga2,weakp,F_THRESHOLD,P_THRESHOLD,mc)
                }
              }
            }
          }
        }
      }
    }
  }
}