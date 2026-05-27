library(readr)
library(MASS)
library(MendelianRandomization)
library(meta)
library(TwoSampleMR)
library(mr.divw)
library(purrr) 
library(stats4)
library(data.table)
library(dplyr)
library(stringr)
library(plinkbinr)
library(ieugwasr)
library(R.utils)
library(openxlsx)
library(foreach)
library(doParallel)
library(doMC)
library(bigstatsr)
library(metafor)
library(ldscr)

############################## 处理函数 #################################
##### 基础分析函数 #####
create_initial_exposure <- function(...) {
  # 捕获所有参数
  args <- list(...)
  
  # 获取参数名称
  arg_names <- as.character(substitute(list(...)))[-1]
  
  # 处理参数命名
  if (is.null(names(args)) || all(names(args) == "")) {
    names(args) <- arg_names
  } else {
    unnamed <- which(names(args) == "")
    names(args)[unnamed] <- arg_names[unnamed]
  }
  
  # ===== 新增：检查并处理缺失值 =====
  # 定义需要检查的列名（常见beta和se列名）
  cols_to_check <- c("beta", "se")
  
  # 1. 检查所有元素是否都是数据框
  is_df <- sapply(args, function(x) is.data.frame(x) || data.table::is.data.table(x))
  if (!all(is_df)) {
    invalid_args <- names(args)[!is_df]
    stop("以下参数不是数据框或data.table: ", paste(invalid_args, collapse = ", "))
  }
  
  # 2. 检查每个数据框的缺失值
  args <- lapply(args, function(df) {
    # 识别需要检查的列（只检查实际存在的列）
    existing_cols <- cols_to_check[cols_to_check %in% names(df)]
    
    if (length(existing_cols) > 0) {
      # 安全提取列数据（兼容data.table）
      col_data <- df[, existing_cols, with = FALSE]
      
      # 检查这些列中的缺失值
      na_rows <- apply(col_data, 1, function(row) any((row==0) | is.na(row) | is.nan(row)))
      
      # 如果有缺失值，报告并删除
      if (any(na_rows)) {
        n_na <- sum(na_rows)
        n_total <- nrow(df)
        cat("检测到", n_na, "/", n_total, "行包含缺失值，已删除\n")
        df <- df[!na_rows, ]
      }
    }
    return(df)
  })
  
  # ===== 新增：确保所有数据框的SNP列完全一致 =====
  # 1. 检查所有元素是否都是数据框且包含SNP列
  is_df_with_snp <- sapply(args, function(x) is.data.frame(x) && "SNP" %in% names(x))
  
  if (!all(is_df_with_snp)) {
    invalid_args <- names(args)[!is_df_with_snp]
    stop("以下参数不是数据框或缺少SNP列: ", paste(invalid_args, collapse = ", "))
  }
  
  # 2. 提取所有SNP集合
  snp_sets <- lapply(args, function(df) unique(df$SNP))
  
  # 3. 计算所有数据框共有的SNP（交集）
  common_snps <- Reduce(intersect, snp_sets)
  
  if (length(common_snps) == 0) {
    stop("所有数据集之间没有共同的SNP！无法创建一致的列表")
  }
  
  # 4. 使用第一个数据框的SNP顺序作为基准
  base_order <- args[[1]]$SNP
  ordered_snps <- base_order[base_order %in% common_snps]
  
  # 5. 统一所有数据框的SNP
  args <- lapply(args, function(df) {
    df <- df[df$SNP %in% common_snps, ]  # 只保留共有SNP
    df <- df[match(ordered_snps, df$SNP), ]  # 按基准顺序排序
    rownames(df) <- NULL  # 重置行名
    return(df)
  })
  # ===== 新增部分结束 =====
  
  # 从全局环境中删除原始数据对象
  objects_to_remove <- arg_names
  existing_objects <- objects_to_remove[objects_to_remove %in% ls(envir = .GlobalEnv)]
  
  if (length(existing_objects) > 0) {
    rm(list = existing_objects, envir = .GlobalEnv)
    cat("已删除以下对象:", paste(existing_objects, collapse = ", "), "\n")
  }
  
  return(args)
}

create_initial_list <- function(...) {
  # 捕获所有参数
  args <- list(...)
  
  # 获取参数名称
  arg_names <- as.character(substitute(list(...)))[-1]
  
  # 如果参数有命名，使用参数名；否则使用变量名
  if (is.null(names(args)) || all(names(args) == "")) {
    names(args) <- arg_names
  } else {
    # 处理部分命名的情况
    unnamed <- which(names(args) == "")
    names(args)[unnamed] <- arg_names[unnamed]
  }
  
  # ===== 新增：检查并处理缺失值 =====
  # 定义需要检查的列名（常见beta和se列名）
  cols_to_check <- c("beta", "se")
  
  # 1. 检查所有元素是否都是数据框
  is_df <- sapply(args, function(x) is.data.frame(x) || data.table::is.data.table(x))
  if (!all(is_df)) {
    invalid_args <- names(args)[!is_df]
    stop("以下参数不是数据框或data.table: ", paste(invalid_args, collapse = ", "))
  }
  
  # 2. 检查每个数据框的缺失值
  # 修复后的缺失值检查函数
  args <- lapply(args, function(df) {
    # 识别需要检查的列（只检查实际存在的列）
    existing_cols <- cols_to_check[cols_to_check %in% names(df)]
    
    if (length(existing_cols) > 0) {
      # 安全提取列数据（兼容data.frame和data.table）
      if (inherits(df, "data.table")) {
        col_data <- df[, ..existing_cols]  # data.table语法
      } else {
        col_data <- df[, existing_cols, drop = FALSE]  # data.frame语法
      }
      
      # 检查这些列中的缺失值
      na_rows <- apply(col_data, 1, function(row) any(row == 0 | is.na(row) | is.nan(row)))
      
      # 如果有缺失值，报告并删除
      if (any(na_rows)) {
        n_na <- sum(na_rows)
        n_total <- nrow(df)
        cat("检测到", n_na, "/", n_total, "行包含缺失值，已删除\n")
        df <- df[!na_rows, ]
      }
    }
    return(df)
  })
  
  # 从全局环境中删除原始数据对象
  objects_to_remove <- arg_names
  existing_objects <- objects_to_remove[objects_to_remove %in% ls(envir = .GlobalEnv)]
  
  if (length(existing_objects) > 0) {
    rm(list = existing_objects, envir = .GlobalEnv)
    cat("已删除以下对象:", paste(existing_objects, collapse = ", "), "\n")
  }
  
  return(args)
}

##### MR-Meta暴露X处理,求并集#####
process_gwas_data <- function(data_list, 
                              ncore,
                              pval_threshold = 5e-8, 
                              clump_kb = clump_kb, 
                              clump_r2 = clump_r2, 
                              clump_p = 5e-4,
                              bfile = "/data_200t/houlei/UKB-cor/data_maf0.01_rs_ref/data_maf0.01_rs_ref",
                              plink_bin = "/opt/R/4.3.2/lib/R/library/plinkbinr/bin/plink_Linux",
                              F_threshold = 10) {
  
  # 1. 预处理每个数据集 - 计算F统计量并过滤弱工具变量
  processed_data <- map(data_list, function(df) {
    
    # 基本数据清洗
    df <- df[!duplicated(df$SNP), ]
    df <- df[df$se > 0, ]
    
    # 重命名pval列（如果存在）
    if ("pval" %in% names(df)) {
      df <- rename(df, pval.exposure = pval)
    }
    
    # ==== 计算F统计量 ====
    # 使用标准公式: F = (beta/se)^2
    df$F_stat <- (df$beta / df$se)^2
    
    # ==== 过滤弱工具变量 ====
    # 先过滤F统计量，再过滤p值（提高效率）
    df <- df %>%
      filter(F_stat >= F_threshold) %>%  # 只保留强工具变量
      filter(pval.exposure < pval_threshold)  # 然后过滤显著SNP
    
    # 仅保留必要列（减少内存占用）
    df %>%
      select(SNP, pval.exposure, F_stat)  # 保留F_stat用于诊断
  })
  
  # 输出F统计量诊断信息
  cat("=== F统计量诊断 ===\n")
  map2(processed_data, names(data_list), function(df, name) {
    if (nrow(df) > 0) {
      cat("数据集", name, ": 保留", nrow(df), "个SNP | 平均F=", 
          round(mean(df$F_stat), 1), "| 最小F=", round(min(df$F_stat), 1), "\n")
    } else {
      cat("数据集", name, ": 无符合条件的SNP\n")
    }
  })
  
  # 2. 对每个数据集执行LD聚类
  registerDoParallel(ncore)
  clumped_data <- foreach(i = seq_along(processed_data), 
                          .packages = c("dplyr"),
                          .errorhandling = "pass") %dopar% {
                            df <- processed_data[[i]]
                            
                            if (nrow(df) > 0) {
                              # 准备数据格式
                              clump_df <- df %>%
                                rename(rsid = SNP, pval = pval.exposure)
                              
                              # 执行LD聚类
                              tryCatch({
                                ld_clump_local(
                                  dat = clump_df,
                                  clump_kb = clump_kb,
                                  clump_r2 = clump_r2,
                                  clump_p = clump_p,
                                  bfile = bfile,
                                  plink_bin = plink_bin
                                )
                              }, error = function(e) {
                                message("LD聚类失败: ", conditionMessage(e))
                                data.frame(rsid = character(), pval = numeric())
                              })
                            } else {
                              data.frame(rsid = character(), pval = numeric())
                            }
                          }
  
  # 3. 整合所有聚类结果
  all_rsid <- character(0)
  
  for (i in seq_along(clumped_data)) {
    df <- clumped_data[[i]]
    if (nrow(df) > 0 && "rsid" %in% names(df)) {
      all_rsid <- union(all_rsid, df$rsid)
    }
  }
  
  union_snps <- data.frame(rsid = all_rsid)
  
  # 4. 使用整合后的SNP筛选原始数据
  filtered_data <- map(names(data_list), function(name) {
    df <- data_list[[name]]
    df %>%
      filter(SNP %in% all_rsid)
  })
  names(filtered_data) <- names(data_list)
  
  # 5. 返回结果
  return(list(
    clumped_results = clumped_data,
    union_snps = union_snps,
    filtered_data = filtered_data
  ))
}

##### 执行传统MR-Meta分析####
align_to_exposure <- function(exp_data, out_data) {
  
  exp_data$SNP <- as.character(exp_data$SNP)
  out_data$SNP <- as.character(out_data$SNP)
  
  # 合并数据
  merged <- merge(exp_data, out_data, by = "SNP", suffixes = c("_exp", "_out"))
  
  # 识别需要翻转的SNP
  flip_condition <- (merged$effect_allele_exp == merged$other_allele_out) & 
    (merged$other_allele_exp == merged$effect_allele_out)
  
  # 应用翻转
  merged[flip_condition, "beta_out"] <- -merged[flip_condition, "beta_out"]
  merged[flip_condition, c("effect_allele_out", "other_allele_out")] <- 
    merged[flip_condition, c("other_allele_out", "effect_allele_out")]
  
  # 分离数据
  exp_aligned <- merged[, c("SNP", "beta_exp", "se_exp", "effect_allele_exp", "other_allele_exp")]
  names(exp_aligned) <- c("SNP", "beta", "se", "effect_allele", "other_allele")
  
  out_aligned <- merged[, c("SNP", "beta_out", "se_out", "effect_allele_out", "other_allele_out")]
  names(out_aligned) <- c("SNP", "beta", "se", "effect_allele", "other_allele")
  
  return(list(exposure = exp_aligned, outcome = out_aligned))
}

run_mr_analysis_no_ref <- function(exposure_list, outcome_list, SNP_list, nx, ny) {
  # 初始化结果容器
  mr_objects <- list()
  res_divw <- data.frame()
  res_fixed <- data.frame()
  res_random <- data.frame()
  res_wme <- data.frame()
  SNPSMR <- c()
  
  # 循环处理每个暴露-结局组合
  for (m in 1:nx) {
    for (n in 1:ny) {
      # 1. 获取共同SNP
      re_SNP <- Reduce(intersect, list(
        SNP_list[[m]]$rsid, 
        exposure_list[[m]]$SNP, 
        outcome_list[[n]]$SNP
      ))
      
      # 记录共同SNP数量
      SNPSMR <- c(SNPSMR, length(re_SNP))
      
      # 2. 提取暴露和结局数据
      exp_data <- exposure_list[[m]] %>% 
        filter(SNP %in% re_SNP) %>%
        distinct(SNP, .keep_all = TRUE)  # 去除重复SNP
      
      out_data <- outcome_list[[n]] %>% 
        filter(SNP %in% re_SNP) %>%
        distinct(SNP, .keep_all = TRUE)  # 去除重复SNP
      
      # 3. 确保等位基因方向一致（以暴露数据为参考）
      aligned_data <- align_to_exposure(exp_data, out_data)
      exp_aligned <- aligned_data$exposure
      out_aligned <- aligned_data$outcome
      
      # 4. 确保SNP顺序一致
      exp_aligned <- exp_aligned[order(exp_aligned$SNP), ]
      out_aligned <- out_aligned[order(out_aligned$SNP), ]
      
      if (!identical(exp_aligned$SNP, out_aligned$SNP)) {
        out_aligned <- out_aligned[match(exp_aligned$SNP, out_aligned$SNP), ]
      }
      
      # 5. 提取向量
      bx <- exp_aligned$beta
      bxse <- exp_aligned$se
      by <- out_aligned$beta
      byse <- out_aligned$se
      
      # 6. 创建MR输入对象
      mr_index <- (m - 1) * ny + n
      mr_objects[[mr_index]] <- mr_input(
        bx = bx, 
        bxse = bxse, 
        by = by, 
        byse = byse
      )
      
      # 7. 执行各种MR方法
      # MR-DIVW
      if (length(re_SNP) >= 1) {
        res_divw_current <- mr.divw(bx, by, bxse, byse)
        res_divw <- rbind(res_divw, c(
          beta_divw = res_divw_current$beta.hat,
          se_divw = res_divw_current$beta.se,
          pval_divw = ci(res_divw_current$beta.hat, res_divw_current$beta.se)$p
        ))
      } else {
        res_divw <- rbind(res_divw, c(NA, NA, NA))
      }
      
      # Fixed-effect IVW
      if (length(re_SNP) >= 1) {
        res_fixed_current <- MendelianRandomization::mr_ivw(mr_objects[[mr_index]], model = "fixed")
        res_fixed <- rbind(res_fixed, c(
          beta_ivwf = res_fixed_current@Estimate,
          se_ivwf = res_fixed_current@StdError,
          pval_ivwf = res_fixed_current@Pvalue
        ))
      } else {
        res_fixed <- rbind(res_fixed, c(NA, NA, NA))
      }
      
      # Random-effect IVW
      if (length(re_SNP) >= 1) {
        res_random_current <- MendelianRandomization::mr_ivw(mr_objects[[mr_index]], model = "random")
        res_random <- rbind(res_random, c(
          beta_ivwr = res_random_current@Estimate,
          se_ivwr = res_random_current@StdError,
          pval_ivwr = res_random_current@Pvalue
        ))
      } else {
        res_random <- rbind(res_random, c(NA, NA, NA))
      }
      
      
      # Weighted Median
      if (length(re_SNP) >= 3) {
        res_wme_current <- MendelianRandomization::mr_median(mr_objects[[mr_index]], iterations = 20)
        res_wme <- rbind(res_wme, c(
          beta_wme = res_wme_current@Estimate,
          se_wme = res_wme_current@StdError,
          pval_wme = res_wme_current@Pvalue
        ))
      } else {
        res_wme <- rbind(res_wme, c(NA, NA, NA))
      }
    }
  }
  
  # 8. 整理最终结果
  final_results <- list(
    mr_objects = mr_objects,
    res_divw = res_divw,
    res_fixed = res_fixed,
    res_random = res_random,
    res_wme = res_wme,
    SNPSMR = SNPSMR
  )
  
  return(final_results)
}

##### 协调暴露数据的等位基因####
build_global_ref_stepwise <- function(data_list) {
  
  # 1. 按SNP数量排序数据集（从最多到最少）
  dataset_sizes <- sapply(data_list, function(df) length(unique(df$SNP)))
  sorted_names <- names(sort(dataset_sizes, decreasing = TRUE))
  
  # 2. 初始化全局参考
  global_ref <- data.table(SNP = character(),
                           effect_allele_ref = character(),
                           other_allele_ref = character())
  
  # 3. 逐步处理每个数据集
  for (name in sorted_names) {
    df <- data_list[[name]]
    
    # 提取当前数据集的SNP和等位基因
    current_snps <- df[, .(SNP, effect_allele, other_allele)]
    
    # 找出当前数据集中尚未在全局参考中的SNP
    new_snps <- current_snps[!SNP %in% global_ref$SNP]
    
    # 将这些新SNP添加到全局参考
    if (nrow(new_snps) > 0) {
      new_ref <- new_snps[, .(SNP, effect_allele_ref = effect_allele, other_allele_ref = other_allele)]
      global_ref <- rbind(global_ref, new_ref)
    }
  }
  
  return(global_ref)
}

# 完整的高效协调函数
harmonize_and_export_efficient <- function(data_list, 
                                           output_dir, 
                                           n_cores) {
  # 1. 构建全局参考
  global_ref <- build_global_ref_stepwise(data_list)
  
  # 2. 并行协调所有数据集
  registerDoParallel(cores = n_cores)
  
  harmonized_list <- foreach(i = 1:length(data_list), .packages = "data.table") %dopar% {
    dataset_name <- names(data_list)[i]
    target_dt <- as.data.table(data_list[[i]])
    
    # 合并全局参考
    merged_dt <- merge(target_dt, global_ref, by = "SNP", all.x = TRUE)
    
    # 识别需要翻转的SNP（等位基因完全相反）
    flip_condition <- (merged_dt$effect_allele == merged_dt$other_allele_ref) & 
      (merged_dt$other_allele == merged_dt$effect_allele_ref)
    
    # 应用翻转：将beta值取反，并交换等位基因
    merged_dt[flip_condition, `:=`(
      beta = -beta,
      effect_allele = other_allele_ref,
      other_allele = effect_allele_ref
    )]
    
    # 移除参考列
    merged_dt[, c("effect_allele_ref", "other_allele_ref") := NULL]
    
    # 重命名pval列
    if ("pval.exposure" %in% names(merged_dt)) {
      setnames(merged_dt, "pval.exposure", "pval")
    }
    
    return(list(name = dataset_name, data = merged_dt))
  }
  
  # 3. 导出所有数据集
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  output_files <- character()
  
  for (item in harmonized_list) {
    output_file <- file.path(output_dir, paste0(item$name, ".txt"))
    fwrite(item$data, output_file, sep = "\t", quote = FALSE)
    output_files <- c(output_files, output_file)
  }
  
  # 4. 返回协调后的数据列表
  harmonized_data <- lapply(harmonized_list, function(x) x$data)
  names(harmonized_data) <- sapply(harmonized_list, function(x) x$name)
  
  return(list(
    global_ref = global_ref,
    harmonized_data = harmonized_data,
    output_files = output_files
  ))
}

##### 读入meta-x结果####
process_metal_results <- function(file_path, 
                                  pval_threshold = 5e-8, 
                                  clump_kb = clump_kb, 
                                  clump_r2 = clump_r2, 
                                  clump_p = 5e-4,
                                  bfile = "/data_200t/houlei/UKB-cor/data_maf0.01_rs_ref/data_maf0.01_rs_ref",
                                  plink_bin = "/opt/R/4.3.2/lib/R/library/plinkbinr/bin/plink_Linux",
                                  F_threshold = 10) {  # 新增F统计量阈值参数
  
  # 1. 读取METAL结果文件
  metal_data <- fread(file_path, 
                      sep = "\t", 
                      header = TRUE,
                      fill = TRUE,
                      na.strings = c("NA", "N/A", ""))
  
  # 2. 数据清洗和预处理
  # 处理列名（兼容不同METAL版本）
  names(metal_data)[names(metal_data) == "P-value"] <- "pval.exposure"
  names(metal_data)[names(metal_data) == "MarkerName"] <- "SNP"
  
  # 转换P值为数值型
  metal_data$pval.exposure <- as.numeric(metal_data$pval.exposure)
  
  # ==== 新增: 计算F统计量 ====
  # 确保有必要的列
  if (!"beta" %in% names(metal_data)) {
    if ("Effect" %in% names(metal_data)) {
      names(metal_data)[names(metal_data) == "Effect"] <- "beta"
    } else {
      stop("无法找到效应值列（beta或Effect）")
    }
  }
  
  if (!"StdErr" %in% names(metal_data)) {
    stop("无法找到标准误列（StdErr）")
  }
  
  # 计算F统计量
  metal_data$F_stat <- (metal_data$beta / metal_data$StdErr)^2
  
  # ==== 新增: 过滤弱工具变量 ====
  # 筛选显著SNP并确保标准误大于0
  metal_data <- metal_data[
    pval.exposure < pval_threshold & 
      StdErr > 0 &
      F_stat >= F_threshold,  # 只保留强工具变量
  ]
  
  # 输出F统计量诊断信息
  cat("=== F统计量诊断 ===\n")
  if (nrow(metal_data) > 0) {
    cat("保留", nrow(metal_data), "个SNP | 平均F=", 
        round(mean(metal_data$F_stat), 1), "| 最小F=", 
        round(min(metal_data$F_stat), 1), "\n")
  } else {
    cat("无符合条件的SNP\n")
  }
  
  # 3. 准备LD聚类数据
  if (nrow(metal_data) > 0) {
    ld_data <- metal_data[, .(rsid = SNP, pval = pval.exposure)]
    
    # 4. 执行LD聚类
    clumped_data <- tryCatch({
      ld_clump_local(
        dat = ld_data,
        clump_kb = clump_kb,
        clump_r2 = clump_r2,
        clump_p = clump_p,
        bfile = bfile,
        plink_bin = plink_bin
      )
    }, error = function(e) {
      message("LD聚类失败: ", conditionMessage(e))
      data.frame(rsid = character(), pval = numeric())
    })
  } else {
    clumped_data <- data.frame(rsid = character(), pval = numeric())
  }
  
  # 5. 返回结果 - 添加F统计量信息
  return(list(
    metal_data = metal_data,
    clumped_data = clumped_data,
    num_clumped = nrow(clumped_data),
    F_threshold = F_threshold,  # 返回使用的阈值
    mean_F_stat = if (nrow(metal_data) > 0) mean(metal_data$F_stat) else NA
  ))
}
##### 筛选结局，生成结局txt####
harmonize_and_export_general <- function(outcome_list, 
                                         global_ref, 
                                         snp_list, 
                                         output_dir, 
                                         n_cores = 4) {
  # 1. 验证输入列表
  names(outcome_list) <- paste0("dataset_", seq_along(outcome_list))
  
  # 2. 筛选包含指定SNP的行
  filtered_datasets <- lapply(outcome_list, function(df) {
    if (!"SNP" %in% names(df)) {
      stop("数据集缺少SNP列")
    }
    df[df$SNP %in% snp_list$rsid, ]
  })
  
  # 3. 并行协调等位基因
  registerDoParallel(cores = n_cores)
  
  harmonized_tables <- foreach(i = 1:length(filtered_datasets), .packages = "data.table") %dopar% {
    dataset_name <- names(filtered_datasets)[i]
    target_data <- as.data.table(filtered_datasets[[i]])
    
    # 合并全局参考数据
    merged_data <- merge(target_data, global_ref, by = "SNP", suffixes = c("", "_ref"))
    
    # 识别需要翻转的SNP（等位基因完全相反）
    flip_condition <- (merged_data$effect_allele == merged_data$other_allele_ref) & 
      (merged_data$other_allele == merged_data$effect_allele_ref)
    
    # 应用翻转：将beta值取反，并交换等位基因
    merged_data[flip_condition, `:=`(
      beta = -beta,
      effect_allele = other_allele_ref,
      other_allele = effect_allele_ref
    )]
    
    # 移除参考列
    merged_data[, c("effect_allele_ref", "other_allele_ref") := NULL]
    
    return(list(name = dataset_name, data = merged_data))
  }
  
  # 4. 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 5. 导出协调后的数据
  output_files <- character()
  
  for (item in harmonized_tables) {
    output_file <- file.path(output_dir, paste0(item$name, ".txt"))
    fwrite(item$data, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    output_files <- c(output_files, output_file)
  }
  
  # 6. 返回协调后的数据
  harmonized_data <- lapply(harmonized_tables, function(x) x$data)
  names(harmonized_data) <- sapply(harmonized_tables, function(x) x$name)
  
  return(list(
    harmonized_data = harmonized_data,
    output_files = output_files
  ))
}

##### 处理Meta MR数据####
process_datasets_Meta_MR <- function(exposure_data_list, 
                                     outcome_data_list, 
                                     exposure_union, 
                                     outcome_union = NULL,
                                     output_dir = NULL) {
  
  # 1. 准备数据集
  # 确保所有数据集都有SNP列
  check_snp_column <- function(df_list) {
    lapply(df_list, function(df) {
      if (!"SNP" %in% names(df)) {
        stop("数据集缺少SNP列")
      }
      return(df)
    })
  }
  
  exposure_data_list <- check_snp_column(exposure_data_list)
  outcome_data_list <- check_snp_column(outcome_data_list)
  
  # 2. 筛选数据（使用union_snps）
  filter_by_union <- function(df_list, union_df) {
    lapply(df_list, function(df) {
      df %>% 
        filter(SNP %in% union_df$rsid) %>%
        filter(se > 0)  # 确保标准误大于0
    })
  }
  
  filtered_exposure <- filter_by_union(exposure_data_list, exposure_union)
  
  # 如果没有提供outcome_union，使用exposure_union
  if (is.null(outcome_union)) {
    outcome_union <- exposure_union
  }
  filtered_outcome <- filter_by_union(outcome_data_list, outcome_union)
  
  # 3. 识别重叠SNP
  all_snps <- c(
    lapply(filtered_exposure, function(df) df$SNP),
    lapply(filtered_outcome, function(df) df$SNP)
  )
  
  overlap_snps <- reduce(all_snps, intersect)
  
  # 4. 处理每个数据集（去重、排序）
  process_data <- function(df, overlap_snps) {
    df %>%
      filter(SNP %in% overlap_snps) %>%
      distinct(SNP, .keep_all = TRUE) %>%
      arrange(SNP)
  }
  
  processed_exposure <- map(filtered_exposure, process_data, overlap_snps = overlap_snps)
  processed_outcome <- map(filtered_outcome, process_data, overlap_snps = overlap_snps)
  
  # 5. 等位基因协调
  # 以第一个暴露数据集为参考
  ref_data <- processed_exposure[[1]][, c("SNP", "effect_allele", "other_allele")]
  
  # 协调暴露数据集
  processed_exposure <- map(processed_exposure, function(df) {
    harmonize_to_reference(df, ref_data)
  })
  
  # 协调结局数据集
  processed_outcome <- map(processed_outcome, function(df) {
    harmonize_to_reference(df, ref_data)
  })
  
  # 6. 返回结果
  return(list(
    processed_exposure = processed_exposure,
    processed_outcome = processed_outcome,
    overlap_snps = overlap_snps,
    num_overlap_snps = length(overlap_snps)
  ))
}

# 等位基因协调函数
harmonize_to_reference <- function(target_data, ref_data) {
  # 合并参考数据
  merged <- merge(target_data, ref_data, by = "SNP", suffixes = c("", "_ref"))
  
  # 识别需要翻转的SNP（等位基因完全相反）
  flip_condition <- (merged$effect_allele == merged$other_allele_ref) & 
    (merged$other_allele == merged$effect_allele_ref)
  
  # 应用翻转：将beta值取反，并交换等位基因
  merged[flip_condition, "beta"] <- -merged[flip_condition, "beta"]
  merged[flip_condition, c("effect_allele", "other_allele")] <- 
    merged[flip_condition, c("other_allele_ref", "effect_allele_ref")]
  
  # 移除参考列
  merged <- merged[, !grepl("_ref$", names(merged))]
  
  return(merged)
}

###combine P
Correlated_Stouffer <- function(exposure_list_ordered, exposure_cor = NULL) {
  z_matrix <- do.call(cbind, lapply(exposure_list_ordered, function(df) {
    as.numeric(df$beta) / as.numeric(df$se)
  }))
  
  if (is.null(exposure_cor)) {
    exposure_cor <- stats::cor(z_matrix, use = "pairwise.complete.obs")
  }
  
  exposure_cor <- as.matrix(exposure_cor)
  exposure_cor[!is.finite(exposure_cor)] <- 0
  diag(exposure_cor) <- 1
  exposure_cor[exposure_cor > 0.99] <- 0.99
  exposure_cor[exposure_cor < -0.99] <- -0.99
  diag(exposure_cor) <- 1
  
  k <- ncol(z_matrix)
  cor_sum <- sum(exposure_cor[upper.tri(exposure_cor)])
  denominator <- sqrt(k + 2 * cor_sum)
  if (!is.finite(denominator) || denominator <= 0) {
    denominator <- sqrt(k)
  }
  
  z_meta <- rowSums(z_matrix, na.rm = FALSE) / denominator
  p_meta <- 2 * pnorm(-abs(z_meta))
  p_meta[!is.finite(p_meta)] <- NA_real_
  
  list(
    p = p_meta,
    z = z_meta,
    exposure_cor = exposure_cor
  )
}

full_dataset_processing <- function(exposure_list, 
                                    outcome_list, 
                                    output_dir = NULL,
                                    assign_to_global = TRUE,
                                    clump_kb = clump_kb,
                                    clump_r2 = clump_r2,
                                    clump_p = 5e-4,
                                    F_threshold = 10,
                                    f_selection = c("all", "any"),
                                    p_combine_method = "stouffer",
                                    exposure_cor = NULL,
                                    bfile = "/data_200t/houlei/UKB-cor/data_maf0.01_rs_ref/data_maf0.01_rs_ref",
                                    plink_bin = "/opt/R/4.3.2/lib/R/library/plinkbinr/bin/plink_Linux") {
  
  f_selection <- match.arg(f_selection)
  p_combine_method <- match.arg(p_combine_method, choices = "stouffer")
  log_filter <- function(...) {
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | SNP-filter | ",
        paste0(..., collapse = ""), "\n", sep = "")
    flush.console()
  }
  log_filter("started; f_selection=", f_selection,
             ", p_combine_method=", p_combine_method,
             ", clump_kb=", clump_kb,
             ", clump_r2=", clump_r2,
             ", F_threshold=", F_threshold)
  required_cols <- c("SNP", "beta", "se", "pval", "effect_allele", "other_allele")
  
  normalize_dataset_list <- function(x, list_name) {
    if (is.data.frame(x)) {
      x <- list(x)
      names(x) <- list_name
    }
    if (!is.list(x) || is.null(x)) {
      stop(list_name, "必须是data.frame列表，当前类型为: ", paste(class(x), collapse = "/"))
    }
    if (is.null(names(x))) {
      names(x) <- paste0(list_name, "_", seq_along(x))
    } else {
      empty_names <- is.na(names(x)) | names(x) == ""
      names(x)[empty_names] <- paste0(list_name, "_", which(empty_names))
    }
    
    for (i in seq_along(x)) {
      if (!is.data.frame(x[[i]])) {
        stop(list_name, "[[", i, "]] (", names(x)[i], ") 必须是data.frame，当前类型为: ",
             paste(class(x[[i]]), collapse = "/"))
      }
      missing_cols <- setdiff(required_cols, names(x[[i]]))
      if (length(missing_cols) > 0 && nrow(x[[i]]) > 0) {
        stop(list_name, "[[", i, "]] (", names(x)[i], ") 缺少必要列: ",
             paste(missing_cols, collapse = ", "))
      }
      x[[i]] <- as.data.frame(x[[i]])
    }
    x
  }
  
  exposure_list <- normalize_dataset_list(exposure_list, "exposure_list")
  outcome_list <- normalize_dataset_list(outcome_list, "outcome_list")
  log_filter("datasets normalized; exposures=", length(exposure_list),
             ", outcomes=", length(outcome_list))
  
  # --- 1.SNP筛选以满足所有三个条件 ---
  F_THRESHOLD <- F_threshold
  P_THRESHOLD <- 5e-8
  
  # --- 条件 1: 必须存在于所有 N (暴露) + M (结局) 个数据集中 ---
  
  all_datasets <- c(exposure_list, outcome_list)
  
  # 过滤掉空的 data.frame
  valid_datasets <- all_datasets[sapply(all_datasets, function(df) is.data.frame(df) && nrow(df) > 0)]
  if (length(valid_datasets) < length(all_datasets)) {
    warning("一个或多个数据集为空或不是data.frame。")
  }
  
  # 获取所有SNP列表
  all_snp_lists <- lapply(valid_datasets, function(df) df$SNP)
  if (length(all_snp_lists) == 0) {
    stop("所有数据集都为空，无法找到共同SNP。")
  }
  
  # 取“N+M 全集”
  snps_in_all_datasets <- Reduce(intersect, all_snp_lists)
  
  if (length(snps_in_all_datasets) == 0) {
    warning("在所有 N (暴露) + M (结局) 数据集中没有共同的SNP。函数将停止。")
    # (根据您函数的结构，您可能需要返回一个空列表)
    # return(list(processed_exposure = list(), processed_outcome = list(), ...))
  }
  cat("--- 条件 1: 在所有N+M个数据集中共同的SNP: ", length(snps_in_all_datasets), "个\n")
  flush.console()
  
  
  # --- 预处理暴露列表 (仅保留共同的SNP) ---
  exposure_list_common <- lapply(exposure_list, function(df) {
    if (is.data.frame(df) && nrow(df) > 0) {
      df[df$SNP %in% snps_in_all_datasets, , drop = FALSE]
    } else {
      data.frame() # 返回一个空框
    }
  })
  
  
  # --- 条件 2: 按指定规则筛选 F > F_threshold ---
  
  # 2a. 计算 F 统计量
  exposure_list_common_F <- lapply(exposure_list_common, function(df) {
    if (is.data.frame(df) && nrow(df) > 0) {
      df$F_stat <- (df$beta / df$se)^2
    } else {
      df$F_stat <- numeric(0)
    }
    return(df)
  })
  
  # 2b. 获取每个暴露中 F > F_threshold 的SNP列表
  list_of_strong_snps <- lapply(exposure_list_common_F, function(df) {
    if (is.data.frame(df) && nrow(df) > 0) {
      df$SNP[df$F_stat > F_THRESHOLD & !is.na(df$F_stat)]
    } else {
      character(0) # 返回空字符向量
    }
  })
  
  # 2c. 严格规则取“F > 10 交集”；放松规则取“F > 10 并集”
  if (length(list_of_strong_snps) < length(exposure_list)) {
    warning("一个或多个暴露列表为空。'F > F_threshold'筛选结果将为空。")
    snps_f_selected <- character(0)
  } else {
    # 检查是否有空列表 (这在R中是必要的)
    if (f_selection == "all" && any(sapply(list_of_strong_snps, length) == 0)) {
      warning("一个或多个暴露数据库没有 F > F_threshold 的SNP。'F > F_threshold 交集'将为空。")
      snps_f_selected <- character(0)
    } else {
      snps_f_selected <- if (f_selection == "all") {
        Reduce(intersect, list_of_strong_snps)
      } else {
        unique(unlist(list_of_strong_snps))
      }
    }
  }
  cat("--- 条件 2: F > ", F_THRESHOLD, " '", f_selection, "' SNP: ", length(snps_f_selected), "个\n", sep = "")
  flush.console()
  
  
  # --- 条件 3: 在 *至少一个* 暴露中 P < 5e-8 ---
  # (这里仍然使用每个 X 自己的 P 值进行筛选，之后 clumping 才使用相关性校正联合P)
  
  # 3a. 获取每个暴露中 P < 5e-8 的SNP列表
  list_of_sig_snps <- lapply(exposure_list_common, function(df) {
    if (is.data.frame(df) && nrow(df) > 0) {
      df$SNP[df$pval < P_THRESHOLD & !is.na(df$pval)]
    } else {
      character(0)
    }
  })
  
  # 3b. 取“P < 5e-8 并集”
  snps_p_any_sig <- unique(unlist(list_of_sig_snps))
  cat("--- 条件 3: P < 5e-8 '并集' SNP: ", length(snps_p_any_sig), "个\n")
  flush.console()
  
  # clumping 用的联合P：对并集SNP按相关性校正 Stouffer Z 重新计算
  prepare_exposure_for_p_combine <- function(exposure_common_list) {
    if (!is.list(exposure_common_list) || is.data.frame(exposure_common_list)) {
      stop("prepare_exposure_for_p_combine需要data.frame列表，当前类型为: ",
           paste(class(exposure_common_list), collapse = "/"))
    }
    invalid_index <- which(!sapply(exposure_common_list, is.data.frame))
    if (length(invalid_index) > 0) {
      invalid_names <- names(exposure_common_list)[invalid_index]
      invalid_names[is.na(invalid_names) | invalid_names == ""] <- invalid_index[is.na(invalid_names) | invalid_names == ""]
      stop("prepare_exposure_for_p_combine发现非data.frame元素: ",
           paste(invalid_names, collapse = ", "))
    }
    
    ref_df <- exposure_common_list[[1]]
    ref_df <- ref_df[!duplicated(ref_df$SNP), , drop = FALSE]
    ref_df <- ref_df[order(ref_df$SNP), , drop = FALSE]
    ref_data <- ref_df[, c("SNP", "effect_allele", "other_allele")]
    ref_snps <- ref_df$SNP
    
    lapply(exposure_common_list, function(df) {
      df <- df[!duplicated(df$SNP), , drop = FALSE]
      df <- harmonize_to_reference(df, ref_data)
      df <- df[match(ref_snps, df$SNP), , drop = FALSE]
      df
    })
  }
  
  
  # --- 最终筛选: 合并条件 2 和 3 ---
  # (条件1 已经通过使用 exposure_list_common 隐式满足了)
  
  final_candidate_snps <- intersect(snps_f_selected, snps_p_any_sig)
  
  cat("--- 最终结果: 满足共同存在、P并集和F规则(", f_selection, ")的候选SNP: ",
      length(final_candidate_snps), "个\n", sep = "")
  flush.console()
  
  # --- 4. 创建用于 Clumping 的 new_exposure_list ---
  
  new_exposure_list <- lapply(exposure_list, function(df) {
    df[df$SNP %in% final_candidate_snps, , drop = FALSE]
  })
  
  # --- 5. Clumping 准备 ---
  
  # ① 对每个暴露数据集进行独立的LD clumping
  # 获取参考数据集的SNP顺序
  ref_df <- new_exposure_list[[1]]
  ref_snps <- ref_df$SNP
  
  #② 确保所有数据集按照参考数据集的SNP顺序排列
  ordered_new_exposure_list <- list()
  for (i in seq_along(new_exposure_list)) {
    df <- new_exposure_list[[i]]
    df_name <- names(new_exposure_list)[i]
    
    # 按照参考数据集的SNP顺序重新排列
    df_ordered <- df[match(ref_snps, df$SNP), , drop = FALSE]
    
    # 检查是否有SNP缺失
    missing_snps <- ref_snps[!ref_snps %in% df$SNP]
    if (length(missing_snps) > 0) {
      warning("数据集", df_name, "缺少", length(missing_snps), "个SNP: ", 
              paste(head(missing_snps, 5), collapse = ", "), 
              ifelse(length(missing_snps) > 5, "...", ""))
    }
    
    ordered_new_exposure_list[[df_name]] <- df_ordered
  }
  
  #④ 生成 clump data：使用相关性校正 Stouffer 合并P值
  log_filter("combining P values for global clump; method=", p_combine_method,
             ", candidate SNPs=", length(ref_snps))
  exposure_list_for_clump_p <- prepare_exposure_for_p_combine(new_exposure_list)
  combined_p_raw <- Correlated_Stouffer(exposure_list_for_clump_p,
                                        exposure_cor = exposure_cor)$p
  combined_p_by_snp <- stats::setNames(combined_p_raw,
                                       exposure_list_for_clump_p[[1]]$SNP)
  combined_p_values <- combined_p_by_snp[ref_snps]
  clump_combine_new <- data.frame(rsid = ref_snps,
                                  pval = combined_p_values)
  clump_combine_new <- na.omit(clump_combine_new)
  log_filter("global clump input ready; non-missing SNPs=", nrow(clump_combine_new))
  
  a_clump <- if (nrow(clump_combine_new) > 0) {
    log_filter("calling ld_clump_local")
    tryCatch(
      ld_clump_local(
        dat = clump_combine_new, 
        clump_kb = clump_kb,
        clump_r2 = clump_r2,
        clump_p = clump_p,
        bfile = bfile,
        plink_bin = plink_bin
      ),
      error = function(e) {
        warning("数据集clump失败: ", conditionMessage(e))
        data.frame(rsid = character(), pval = numeric())
      }
    )
  } else {
    warning("没有可用于clump的候选SNP")
    data.frame(rsid = character(), pval = numeric())
  }
  
  if (is.data.frame(a_clump) && nrow(a_clump) > 0) {
    for (k in 1:length(new_exposure_list)) {
      new_exposure_list[[k]] <- new_exposure_list[[k]][new_exposure_list[[k]]$SNP %in% a_clump$rsid, , drop = FALSE]
    }
    cat("数据集clump后保留", nrow(a_clump), "个独立SNP\n")
    flush.console()
  } else {
    warning("数据集clump后无SNP保留")
    for (k in 1:length(new_exposure_list)) {
      new_exposure_list[[k]] <- data.frame()
    }
  }
  
  
  # ⑤ 筛选数据（a_clump）
  filter_by_union <- function(df, a_clump) {
    # 确保输入是数据框
    if (!is.data.frame(df) || nrow(df) == 0) {
      return(data.frame())
    }
    
    # 确保返回数据框
    result <- df[df$SNP %in% a_clump$rsid & df$se > 0, , drop = FALSE]
    
    # 检查是否为空
    if (nrow(result) == 0) {
      warning("筛选后没有数据，返回空数据框")
      return(data.frame())
    }
    
    return(result)
  }
  
  filtered_exposure <- lapply(new_exposure_list, filter_by_union, a_clump = a_clump)
  filtered_outcome <- lapply(outcome_list, filter_by_union, a_clump = a_clump)
  
  # ⑥. 识别重叠SNP
  all_snps <- c(
    lapply(filtered_exposure, function(df) {
      if (is.data.frame(df) && nrow(df) > 0) df$SNP else character(0)
    }),
    lapply(filtered_outcome, function(df) {
      if (is.data.frame(df) && nrow(df) > 0) df$SNP else character(0)
    })
  )
  
  overlap_snps <- Reduce(intersect, all_snps)
  log_filter("post-clump overlap SNPs across exposure/outcome datasets: ",
             length(overlap_snps))
  
  # 6. 处理每个数据集（去重、排序）
  process_data <- function(df) {
    # 检查是否为空
    if (!is.data.frame(df) || nrow(df) == 0) {
      return(data.frame())
    }
    
    # 确保返回数据框
    result <- df[df$SNP %in% overlap_snps, , drop = FALSE]
    if (nrow(result) == 0) {
      return(data.frame())
    }
    
    result <- result[!duplicated(result$SNP), , drop = FALSE]
    result <- result[order(result$SNP), , drop = FALSE]
    
    return(result)
  }
  
  processed_exposure <- lapply(filtered_exposure, process_data)
  processed_outcome <- lapply(filtered_outcome, process_data)
  
  # 7. 等位基因协调
  # 检查是否有数据可处理
  valid_exposure <- which(sapply(processed_exposure, function(df) is.data.frame(df) && nrow(df) > 0))
  
  if (length(valid_exposure) == 0) {
    warning("没有有效的暴露数据可用于等位基因协调")
    results <- list(
      processed_exposure = processed_exposure,
      processed_outcome = processed_outcome,
      overlap_snps = overlap_snps,
      num_overlap_snps = length(overlap_snps)
    )
    return(results)
  }
  
  # 以第一个有效的暴露数据集为参考
  ref_index <- valid_exposure[1]
  ref_data <- processed_exposure[[ref_index]][, c("SNP", "effect_allele", "other_allele")]
  
  # 协调所有数据集
  harmonize_to_reference <- function(target_data, ref_data) {
    
    # 合并参考数据
    merged <- merge(target_data, ref_data, by = "SNP", suffixes = c("", "_ref"))
    
    # 识别需要翻转的SNP
    flip_condition <- (merged$effect_allele == merged$other_allele_ref) & 
      (merged$other_allele == merged$effect_allele_ref)
    
    # 应用翻转
    merged[flip_condition, "beta"] <- -merged[flip_condition, "beta"]
    merged <- merged[, !names(merged) %in% c("effect_allele", "other_allele"), drop = FALSE]
    names(merged)[names(merged) == "effect_allele_ref"] <-"effect_allele"
    names(merged)[names(merged) == "other_allele_ref"] <-"other_allele"
    
    return(merged)
  }
  
  # 协调暴露数据集
  processed_exposure <- lapply(processed_exposure, function(df) {
    harmonize_to_reference(df, ref_data)
  })
  
  # 协调结局数据集
  processed_outcome <- lapply(processed_outcome, function(df) {
    harmonize_to_reference(df, ref_data)
  })
  
  # 8. 返回结果
  results <- list(
    processed_exposure = processed_exposure,
    processed_outcome = processed_outcome,
    overlap_snps = overlap_snps,
    num_overlap_snps = length(overlap_snps),
    f_selection = f_selection,
    p_combine_method = p_combine_method,
    num_final_candidate_snps = length(final_candidate_snps),
    num_clumped_snps = if (is.data.frame(a_clump)) nrow(a_clump) else 0
  )
  
  return(results)
}

##### NEW method #####
run_new_mr_analysis <- function(exposure_filter_list, outcome_filter_list, original_outcome_list,
                                N_list, ancestry = "EUR", bootstrap_time,
                                relaxed_nome_bootstrap_time = NULL,
                                original_exposure_list = NULL,
                                ldsc_cache = NULL,
                                run_relaxed_nome = TRUE,
                                ldsc_scope = c("all", "yy")) {
  ldsc_scope <- match.arg(ldsc_scope)
  if (isTRUE(run_relaxed_nome)) {
    ldsc_scope <- "all"
  }
  log_step <- function(...) {
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | JointMR | ",
        paste0(..., collapse = ""), "\n", sep = "")
    flush.console()
  }
  should_report_progress <- function(k, total) {
    if (!is.finite(total) || total <= 0) {
      return(FALSE)
    }
    k == 1 || k == total || k %% max(1, floor(total / 5)) == 0
  }
  log_step("run_new_mr_analysis started; run_relaxed_nome=", run_relaxed_nome,
           ", ldsc_scope=", ldsc_scope,
           ", bootstrap_time=", bootstrap_time,
           ", relaxed_nome_bootstrap_time=",
           ifelse(is.null(relaxed_nome_bootstrap_time), bootstrap_time,
                  relaxed_nome_bootstrap_time))
  
  ginv_solve <- function(Sigma, rhs, tol = sqrt(.Machine$double.eps)) {
    Sigma <- as.matrix(Sigma)
    Sigma <- (Sigma + t(Sigma)) / 2
    rhs <- as.matrix(rhs)
    MASS::ginv(Sigma, tol = tol) %*% rhs
  }
  
  make_psd <- function(Sigma, eps = sqrt(.Machine$double.eps)) {
    Sigma <- as.matrix(Sigma)
    Sigma[!is.finite(Sigma)] <- 0
    Sigma <- (Sigma + t(Sigma)) / 2
    eig <- eigen(Sigma, symmetric = TRUE)
    values <- pmax(eig$values, eps)
    eig$vectors %*% diag(values, length(values)) %*% t(eig$vectors)
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
                                        data_y_se, cor_matrix, nx, ny,
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
    
    alpha_hat[!is.finite(alpha_hat)] <- 0
    alpha_se[!is.finite(alpha_se)] <- 0
    
    rho_y <- matrix(1, nrow = ny, ncol = ny)
    for (m in seq_len(ny)) {
      for (q in seq_len(ny)) {
        rho_y[m, q] <- cor_matrix[(m - 1) * nx + 1, (q - 1) * nx + 1]
      }
    }
    
    V_alpha <- rho_y * outer(alpha_se, alpha_se)
    V_alpha[!is.finite(V_alpha)] <- 0
    
    list(alpha_hat = alpha_hat, alpha_se = alpha_se, V_alpha = V_alpha)
  }
  
  build_plugin_adjusted_inputs <- function(data_x_beta, data_x_se, data_y_beta,
                                           data_y_se, cor_matrix, nx, ny,
                                           valid_indices,
                                           egger_indices = NULL) {
    n_ratio <- nx * ny
    
    egger_fit <- estimate_egger_intercepts(
      data_x_beta = data_x_beta,
      data_x_se = data_x_se,
      data_y_beta = data_y_beta,
      data_y_se = data_y_se,
      cor_matrix = cor_matrix,
      nx = nx,
      ny = ny,
      snp_indices = egger_indices
    )
    
    g <- nrow(data_x_beta)
    WR_adj_full <- matrix(NA_real_, nrow = n_ratio, ncol = g)
    seWR_matrix_full <- matrix(NA_real_, nrow = n_ratio, ncol = g)
    
    for (m in seq_len(ny)) {
      for (n in seq_len(nx)) {
        row_id <- (m - 1) * nx + n
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
      Omega_j <- cor_matrix * outer(se_vector, se_vector)
      
      beta_x_vec <- as.numeric(beta_x_selected[idx, ])
      B_j <- diag(rep(1 / beta_x_vec, times = ny), nrow = n_ratio)
      BA_j <- B_j %*% A
      
      dimnames(Omega_j) <- list(rownames(WR_adj), rownames(WR_adj))
      base_cov_list[[idx]] <- Omega_j
      BA_list[[idx]] <- BA_j
    }
    
    list(
      WR_adj = WR_adj,
      seWR = seWR,
      base_cov_list = base_cov_list,
      BA_list = BA_list,
      alpha_hat = egger_fit$alpha_hat,
      alpha_se = egger_fit$alpha_se,
      V_alpha = egger_fit$V_alpha
    )
  }
  
  MLE_S <- function(nx, ny, cov_list, WR_matrix, SNPnew) {
    n <- nx * ny
    I <- matrix(1, n, 1)
    precision_sum <- 0
    weighted_sum <- 0
    
    for (j in seq_len(SNPnew)) {
      Sigma <- cov_list[[j]]
      Sigma <- as.matrix(Sigma)
      Sigma <- (Sigma + t(Sigma)) / 2
      eig <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
      min_eig <- min(eig, na.rm = TRUE)
      if (!is.finite(min_eig) || min_eig <= sqrt(.Machine$double.eps)) {
        jitter <- if (is.finite(min_eig)) abs(min_eig) + sqrt(.Machine$double.eps) else sqrt(.Machine$double.eps)
        Sigma <- Sigma + diag(jitter, nrow(Sigma))
      }
      Sigma_inv <- MASS::ginv(Sigma)
      precision_sum <- precision_sum + as.numeric(crossprod(I, Sigma_inv %*% I))
      weighted_sum <- weighted_sum + as.numeric(crossprod(I, Sigma_inv %*% WR_matrix[, j, drop = FALSE]))
    }
    
    theta <- weighted_sum / precision_sum
    attr(theta, "se") <- sqrt(1 / precision_sum)
    theta
  }
  
  build_wald_ratio_inputs <- function(data_x_beta, data_y_beta, data_y_se,
                                      cor_matrix, nx, ny, valid_indices,
                                      alpha_hat = rep(0, ny),
                                      tau2 = 0) {
    n_ratio <- nx * ny
    g <- nrow(data_x_beta)
    WR_full <- matrix(NA_real_, nrow = n_ratio, ncol = g)
    seWR_full <- matrix(NA_real_, nrow = n_ratio, ncol = g)
    
    for (m in seq_len(ny)) {
      for (n in seq_len(nx)) {
        row_id <- (m - 1) * nx + n
        WR_full[row_id, ] <- (data_y_beta[, m] - alpha_hat[m]) / data_x_beta[, n]
        seWR_full[row_id, ] <- sqrt(data_y_se[, m]^2 / data_x_beta[, n]^2)
      }
    }
    
    WR <- WR_full[, valid_indices, drop = FALSE]
    seWR <- seWR_full[, valid_indices, drop = FALSE]
    beta_x_selected <- data_x_beta[valid_indices, , drop = FALSE]
    base_cov_list <- vector("list", length(valid_indices))
    
    for (idx in seq_along(valid_indices)) {
      se_vector <- seWR[, idx]
      Omega_j <- cor_matrix * outer(se_vector, se_vector) + tau2
      dimnames(Omega_j) <- list(rownames(WR), rownames(WR))
      base_cov_list[[idx]] <- Omega_j
    }
    
    list(
      WR = WR,
      seWR = seWR,
      beta_x_selected = beta_x_selected,
      base_cov_list = base_cov_list
    )
  }
  
  ####tau2: overlap-aware GLS + DL across SNPs
  collapse_gls_snp <- function(theta_vec, Sigma) {
    Sigma <- make_psd(Sigma)
    one <- rep(1, length(theta_vec))
    Sigma_inv <- MASS::ginv(Sigma)
    precision <- as.numeric(crossprod(one, Sigma_inv %*% one))
    if (!is.finite(precision) || precision <= 0) {
      stop("Invalid GLS precision in collapse_gls_snp: ", precision)
    }
    
    theta_hat <- as.numeric(crossprod(one, Sigma_inv %*% theta_vec) / precision)
    se_hat <- sqrt(1 / precision)
    if (!is.finite(theta_hat) || !is.finite(se_hat) || se_hat <= 0) {
      stop("Invalid collapsed SNP theta/se in collapse_gls_snp: theta=",
           theta_hat, ", se=", se_hat)
    }
    
    c(theta_hat = theta_hat, se_hat = se_hat)
  }
  
  calculate_tau2_DL <- function(theta_j, se_j) {
    keep <- is.finite(theta_j) & is.finite(se_j) & se_j > 0
    if (sum(keep) < length(theta_j)) {
      warning("Dropping ", length(theta_j) - sum(keep),
              " SNPs with non-finite/non-positive collapsed theta or se for tau2 DL.")
    }
    theta_j <- theta_j[keep]
    se_j <- se_j[keep]
    if (length(theta_j) < 2) {
      stop("Fewer than 2 valid SNP-level estimates remain for tau2 DL.")
    }
    w_j <- 1 / (se_j^2)
    theta_fixed <- sum(w_j * theta_j) / sum(w_j)
    Q <- sum(w_j * (theta_j - theta_fixed)^2)
    k <- length(theta_j)
    df <- k - 1
    C <- sum(w_j) - sum(w_j^2) / sum(w_j)
    if (!is.finite(Q) || !is.finite(C) || C <= 0 || df <= 0) {
      stop("Invalid tau2 DL components: Q=", Q,
           ", C=", C, ", df=", df)
    }
    tau2 <- max(0, (Q - df) / C)
    I2 <- if (Q > 0) max(0, (Q - df) / Q) * 100 else 0
    
    list(
      theta_fixed = theta_fixed,
      tau2 = tau2,
      Q = Q,
      df = df,
      p_value = 1 - pchisq(Q, df),
      I2 = I2
    )
  }
  
  get_sample_size <- function(df, name) {
    if (!is.null(N_list[[name]]) && is.finite(as.numeric(N_list[[name]][1]))) {
      return(as.numeric(N_list[[name]][1]))
    }
    candidate_cols <- c("N", "n", "samplesize", "sample_size", "N_total",
                        "n_total", "samplesize.exposure", "samplesize.outcome")
    hit <- candidate_cols[candidate_cols %in% names(df)]
    if (length(hit) > 0) {
      vals <- suppressWarnings(as.numeric(df[[hit[1]]]))
      vals <- vals[is.finite(vals) & vals > 0]
      if (length(vals) > 0) {
        return(stats::median(vals))
      }
    }
    NA_real_
  }
  
  prepare_ldsc_sumstats <- function(data_list, label) {
    if (is.null(data_list)) {
      return(list())
    }
    out <- list()
    for (name in names(data_list)) {
      df <- as.data.frame(data_list[[name]])
      N_val <- get_sample_size(df, name)
      if (!is.finite(N_val) || N_val <= 0) {
        warning("No valid sample size found for ", label, " dataset ", name,
                "; LDSC pairwise correlation involving this dataset will be set to 0.")
        next
      }
      munged_df <- data.frame(
        SNP = df$SNP,
        A1 = df$effect_allele,
        A2 = df$other_allele,
        Z = df$beta / df$se,
        N = N_val
      )
      munged_df <- munged_df[is.finite(munged_df$Z), , drop = FALSE]
      if (nrow(munged_df) > 0) {
        out[[name]] <- munged_df
      }
    }
    out
  }
  
  extract_ldsc_value <- function(ldsc_res, value_type = c("rg", "intercept")) {
    value_type <- match.arg(value_type)
    if (is.null(ldsc_res)) {
      return(NA_real_)
    }
    if (value_type == "rg") {
      if (is.null(ldsc_res$rg)) {
        return(NA_real_)
      }
      rg_tab <- as.data.frame(ldsc_res$rg)
      if ("rg" %in% names(rg_tab)) {
        return(as.numeric(rg_tab$rg[1]))
      }
      return(NA_real_)
    }
    
    # For XY sample-overlap-related covariance, use the off-diagonal element
    # of the LDSC intercept matrix rather than the genetic correlation rg.
    if (!is.null(ldsc_res$raw) && !is.null(ldsc_res$raw$I)) {
      I_mat <- as.matrix(ldsc_res$raw$I)
      if (nrow(I_mat) >= 2 && ncol(I_mat) >= 2 &&
          is.finite(as.numeric(I_mat[1, 2]))) {
        return(as.numeric(I_mat[1, 2]))
      }
    }
    
    # Fallback for LDSC wrapper versions that expose the cross-trait intercept
    # directly in the rg table.
    if (!is.null(ldsc_res$rg)) {
      rg_tab <- as.data.frame(ldsc_res$rg)
      intercept_cols <- grep("intercept|gcov.*int|int.*gcov|cov.*int|int.*cov",
                             names(rg_tab), ignore.case = TRUE, value = TRUE)
      if (length(intercept_cols) > 0) {
        return(as.numeric(rg_tab[[intercept_cols[1]]][1]))
      }
    }
    
    NA_real_
  }
  
  run_ldsc_pair <- function(sumstats_i, sumstats_j, name_i, name_j) {
    munge_pair <- list(sumstats_i, sumstats_j)
    names(munge_pair) <- c(name_i, name_j)
    tryCatch(
      ldscr::ldsc_rg(
        munged_sumstats = munge_pair,
        ancestry = ancestry
      ),
      error = function(e) {
        warning(paste("ldsc_rg failed for pair", name_i, "vs", name_j,
                      ". Error:", e$message))
        NULL
      }
    )
  }
  
  infer_dataset_source <- function(dataset_name) {
    dataset_name <- toupper(gsub("[^A-Z0-9]+", "_", dataset_name))
    if (grepl("UKB|UK_BIOBANK|LOH", dataset_name)) {
      return("UKB")
    }
    if (grepl("MVP|MILLION", dataset_name)) {
      return("MVP")
    }
    if (grepl("FINN|FINNGEN", dataset_name)) {
      return("FINNGEN")
    }
    if (grepl("EPIC", dataset_name)) {
      return("EPIC")
    }
    if (grepl("SWED", dataset_name)) {
      return("SWEDISH")
    }
    if (grepl("GLGC|WILLER", dataset_name)) {
      return("GLGC")
    }
    NA_character_
  }
  
  has_xy_sample_overlap <- function(exposure_name, outcome_name) {
    exposure_source <- infer_dataset_source(exposure_name)
    outcome_source <- infer_dataset_source(outcome_name)
    !is.na(exposure_source) &&
      !is.na(outcome_source) &&
      exposure_source == outcome_source
  }
  
  compute_ldsc_matrices <- function(original_exposure_list, original_outcome_list,
                                    scope = c("all", "yy")) {
    scope <- match.arg(scope)
    if (is.null(original_exposure_list)) {
      original_exposure_list <- filtered_exposure_list
    }
    exposure_names <- names(original_exposure_list)
    outcome_names <- names(original_outcome_list)
    nx0 <- length(exposure_names)
    ny0 <- length(outcome_names)
    
    exposure_sumstats <- if (scope == "all") {
      prepare_ldsc_sumstats(original_exposure_list, "exposure")
    } else {
      list()
    }
    outcome_sumstats <- prepare_ldsc_sumstats(original_outcome_list, "outcome")
    
    rg_x <- matrix(0, nrow = nx0, ncol = nx0,
                   dimnames = list(exposure_names, exposure_names))
    diag(rg_x) <- 1
    if (scope == "all" && nx0 > 1) {
      for (i in seq_len(nx0 - 1)) {
        for (j in (i + 1):nx0) {
          name_i <- exposure_names[i]
          name_j <- exposure_names[j]
          if (!is.null(exposure_sumstats[[name_i]]) &&
              !is.null(exposure_sumstats[[name_j]])) {
            cat("  Calculating XX rg between:", name_i, "and", name_j, "\n")
            ldsc_pair <- run_ldsc_pair(exposure_sumstats[[name_i]],
                                       exposure_sumstats[[name_j]],
                                       name_i, name_j)
            rg_val <- extract_ldsc_value(ldsc_pair, "rg")
            if (is.finite(rg_val)) {
              rg_x[name_i, name_j] <- rg_val
              rg_x[name_j, name_i] <- rg_val
            }
          }
        }
      }
    }
    
    rg_y <- matrix(0, nrow = ny0, ncol = ny0,
                   dimnames = list(outcome_names, outcome_names))
    diag(rg_y) <- 1
    if (ny0 > 1) {
      for (i in seq_len(ny0 - 1)) {
        for (j in (i + 1):ny0) {
          name_i <- outcome_names[i]
          name_j <- outcome_names[j]
          if (!is.null(outcome_sumstats[[name_i]]) &&
              !is.null(outcome_sumstats[[name_j]])) {
            cat("  Calculating YY rg between:", name_i, "and", name_j, "\n")
            ldsc_pair <- run_ldsc_pair(outcome_sumstats[[name_i]],
                                       outcome_sumstats[[name_j]],
                                       name_i, name_j)
            rg_val <- extract_ldsc_value(ldsc_pair, "rg")
            if (is.finite(rg_val)) {
              rg_y[name_i, name_j] <- rg_val
              rg_y[name_j, name_i] <- rg_val
            }
          }
        }
      }
    }
    
    xy_rho <- matrix(0, nrow = nx0, ncol = ny0,
                     dimnames = list(exposure_names, outcome_names))
    xy_overlap <- matrix(FALSE, nrow = nx0, ncol = ny0,
                         dimnames = list(exposure_names, outcome_names))
    xy_rho_available <- matrix(FALSE, nrow = nx0, ncol = ny0,
                               dimnames = list(exposure_names, outcome_names))
    if (scope == "all") {
      for (i in seq_len(nx0)) {
        for (j in seq_len(ny0)) {
          name_x <- exposure_names[i]
          name_y <- outcome_names[j]
          if (!has_xy_sample_overlap(name_x, name_y)) {
            cat("  No XY sample overlap by dataset name:", name_x, "vs", name_y, "\n")
            next
          }
          xy_overlap[name_x, name_y] <- TRUE
          if (!is.null(exposure_sumstats[[name_x]]) &&
              !is.null(outcome_sumstats[[name_y]])) {
            cat("  Calculating XY LDSC intercept between:", name_x, "and", name_y, "\n")
            ldsc_pair <- run_ldsc_pair(exposure_sumstats[[name_x]],
                                       outcome_sumstats[[name_y]],
                                       name_x, name_y)
            intercept_val <- extract_ldsc_value(ldsc_pair, "intercept")
            if (is.finite(intercept_val)) {
              xy_rho[name_x, name_y] <- intercept_val
              xy_rho_available[name_x, name_y] <- TRUE
            } else {
              warning("XY LDSC overlap intercept was not available for ",
                      name_x, " vs ", name_y,
                      "; xy_rho was set to 0 for this pair.")
            }
          }
        }
      }
    }
    
    rg_x[!is.finite(rg_x)] <- 0
    rg_y[!is.finite(rg_y)] <- 0
    xy_rho[!is.finite(xy_rho)] <- 0
    rg_x[rg_x > 0.99] <- 0.99
    rg_x[rg_x < -0.99] <- -0.99
    rg_y[rg_y > 0.99] <- 0.99
    rg_y[rg_y < -0.99] <- -0.99
    xy_rho[xy_rho > 0.99] <- 0.99
    xy_rho[xy_rho < -0.99] <- -0.99
    diag(rg_x) <- 1
    diag(rg_y) <- 1
    
    list(rg_x = rg_x, rg_y = rg_y, xy_rho = xy_rho,
         xy_overlap = xy_overlap, xy_rho_available = xy_rho_available)
  }
  
  build_ratio_cor_matrix <- function(rg_y, nx, ny, outcome_names) {
    n_total <- nx * ny
    out <- matrix(0, nrow = n_total, ncol = n_total)
    for (i in seq_len(ny)) {
      for (j in seq_len(ny)) {
        start_row_i <- (i - 1) * nx + 1
        end_row_i <- i * nx
        start_col_j <- (j - 1) * nx + 1
        end_col_j <- j * nx
        if (i == j) {
          out[start_row_i:end_row_i, start_col_j:end_col_j] <- 1
        } else {
          out[start_row_i:end_row_i, start_col_j:end_col_j] <-
            rg_y[outcome_names[i], outcome_names[j]]
        }
      }
    }
    out
  }
  
  stack_outcome_values <- function(beta_y_vec, nx) {
    rep(beta_y_vec, each = nx)
  }
  
  stack_exposure_values <- function(beta_x_vec, ny) {
    rep(beta_x_vec, times = ny)
  }
  
  jointmr_profile_q <- function(theta, beta_x, se_x, beta_y, se_y,
                                alpha_hat, V_alpha, rho_y, rho_x,
                                xy_rho, tau2 = 0) {
    q_value <- 0
    n_snp <- nrow(beta_x)
    nx_local <- ncol(beta_x)
    ny_local <- ncol(beta_y)
    
    for (idx in seq_len(n_snp)) {
      beta_x_vec <- as.numeric(beta_x[idx, ])
      se_x_vec <- as.numeric(se_x[idx, ])
      beta_y_vec <- as.numeric(beta_y[idx, ] - alpha_hat)
      se_y_vec <- as.numeric(se_y[idx, ])
      
      cov_x <- rho_x * outer(se_x_vec, se_x_vec)
      cov_y <- rho_y * outer(se_y_vec, se_y_vec)
      if (!is.null(V_alpha) && all(is.finite(V_alpha)) &&
          max(abs(V_alpha)) > .Machine$double.eps) {
        cov_y <- cov_y + V_alpha
      }
      cov_xy <- xy_rho * outer(se_x_vec, se_y_vec)
      
      Sigma <- rbind(
        cbind(cov_x, cov_xy),
        cbind(t(cov_xy), cov_y)
      )
      Sigma <- make_psd(Sigma)
      
      observed <- matrix(c(beta_x_vec, beta_y_vec), ncol = 1)
      design <- matrix(c(rep(1, nx_local), rep(theta, ny_local)), ncol = 1)
      denom <- as.numeric(t(design) %*% ginv_solve(Sigma, design))
      if (!is.finite(denom) || denom <= 0) {
        return(Inf)
      }
      gamma_hat <- as.numeric(t(design) %*% ginv_solve(Sigma, observed) / denom)
      
      if (tau2 > 0) {
        cov_y_re <- cov_y + tau2 * gamma_hat^2 *
          matrix(1, nrow = ny_local, ncol = ny_local)
        Sigma <- rbind(
          cbind(cov_x, cov_xy),
          cbind(t(cov_xy), cov_y_re)
        )
        Sigma <- make_psd(Sigma)
        denom <- as.numeric(t(design) %*% ginv_solve(Sigma, design))
        if (!is.finite(denom) || denom <= 0) {
          return(Inf)
        }
        gamma_hat <- as.numeric(t(design) %*% ginv_solve(Sigma, observed) / denom)
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
  
  estimate_relaxed_nome <- function(beta_x, se_x, beta_y, se_y,
                                    alpha_hat, V_alpha, rho_y, rho_x,
                                    xy_rho, tau2 = 0) {
    theta_start_values <- c()
    ratio_values <- c()
    nx_local <- ncol(beta_x)
    ny_local <- ncol(beta_y)
    
    for (m in seq_len(ny_local)) {
      for (n in seq_len(nx_local)) {
        bx <- beta_x[, n]
        by <- beta_y[, m] - alpha_hat[m]
        keep <- is.finite(bx) & is.finite(by) & abs(bx) > .Machine$double.eps
        if (any(keep)) {
          theta_start_values <- c(theta_start_values,
                                  sum(bx[keep] * by[keep]) / sum(bx[keep]^2))
          ratio_values <- c(ratio_values, by[keep] / bx[keep])
        }
      }
    }
    
    theta_start <- mean(theta_start_values, na.rm = TRUE)
    if (!is.finite(theta_start)) theta_start <- 0
    
    ratio_values <- ratio_values[is.finite(ratio_values)]
    if (length(ratio_values) >= 5) {
      ratio_range <- as.numeric(stats::quantile(ratio_values, c(0.01, 0.99),
                                                na.rm = TRUE))
      lower <- min(-1, theta_start - 1, ratio_range[1] - 0.5)
      upper <- max(1, theta_start + 1, ratio_range[2] + 0.5)
    } else {
      lower <- theta_start - 1
      upper <- theta_start + 1
    }
    
    q_fun <- function(theta) {
      jointmr_profile_q(theta, beta_x, se_x, beta_y, se_y,
                        alpha_hat, V_alpha, rho_y, rho_x, xy_rho,
                        tau2 = tau2)
    }
    
    grid <- seq(lower, upper, length.out = 201)
    grid_q <- vapply(grid, q_fun, numeric(1))
    best <- which.min(grid_q)
    if (length(best) == 0 || !is.finite(grid_q[best])) {
      return(c(theta_hat = NA_real_, theta_se = NA_real_,
               theta_p_value = NA_real_))
    }
    
    opt_interval <- if (best == 1 || best == length(grid)) {
      c(lower, upper)
    } else {
      c(grid[best - 1], grid[best + 1])
    }
    
    opt <- stats::optimize(q_fun, interval = opt_interval)
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
    
    theta_se <- if (is.finite(hessian) && hessian > 0) sqrt(2 / hessian) else NA_real_
    theta_p_value <- if (is.finite(theta_se) && theta_se > 0) {
      2 * pnorm(-abs(theta_hat / theta_se))
    } else {
      NA_real_
    }
    
    c(theta_hat = theta_hat, theta_se = theta_se,
      theta_p_value = theta_p_value)
  }
  
  
  # 主分析流程
  filtered_exposure_list <- exposure_filter_list
  filtered_outcome_list <- outcome_filter_list
  
  nx <- length(exposure_filter_list)
  ny <- length(outcome_filter_list)
  
  data_x_beta <- do.call(cbind, lapply(filtered_exposure_list, function(df) df$beta))
  data_x_se <- do.call(cbind, lapply(filtered_exposure_list, function(df) df$se))
  data_y_beta <- do.call(cbind, lapply(filtered_outcome_list, function(df) df$beta))
  data_y_se <- do.call(cbind, lapply(filtered_outcome_list, function(df) df$se))
  log_step("data matrices built; SNP rows=", nrow(data_x_beta),
           ", nx=", ncol(data_x_beta), ", ny=", ncol(data_y_beta))
  
  if (!all(dim(data_x_beta)[1] == c(nrow(data_x_se), nrow(data_y_beta), nrow(data_y_se)))) {
    stop("Exposure and outcome matrices have inconsistent SNP counts after harmonisation.")
  }
  
  snp_ids <- if (nx > 0 && !is.null(filtered_exposure_list[[1]]$SNP)) {
    filtered_exposure_list[[1]]$SNP
  } else {
    paste0("SNP_", seq_len(nrow(data_x_beta)))
  }
  
  # 创建行名
  colname <- NULL
  for (n in 1:ny) {
    for (m in 1:nx) {
      colname1 <- paste0("X", m, "-Y", n)
      colname <- c(colname, colname1)
    }
  }
  
  nx <- length(exposure_filter_list)
  ny <- length(outcome_filter_list)
  outcome_names <- names(original_outcome_list)
  exposure_names_for_ldsc <- names(if (!is.null(original_exposure_list)) {
    original_exposure_list
  } else {
    exposure_filter_list
  })
  
  if (is.null(ldsc_cache)) {
    if (ldsc_scope == "all") {
      log_step("Preparing LDSC matrices for XX/YY/XY correlations")
    } else {
      log_step("Preparing LDSC matrices for YY correlations only")
    }
    ldsc_cache <- compute_ldsc_matrices(original_exposure_list,
                                        original_outcome_list,
                                        scope = ldsc_scope)
    log_step("LDSC matrices finished")
  } else {
    log_step("Using cached LDSC matrices")
  }
  
  rg_matrix <- ldsc_cache$rg_y[outcome_names, outcome_names, drop = FALSE]
  rg_x_matrix <- ldsc_cache$rg_x[exposure_names_for_ldsc, exposure_names_for_ldsc, drop = FALSE]
  xy_rho_matrix <- ldsc_cache$xy_rho[exposure_names_for_ldsc, outcome_names, drop = FALSE]
  xy_overlap_matrix <- if (!is.null(ldsc_cache$xy_overlap)) {
    ldsc_cache$xy_overlap[exposure_names_for_ldsc, outcome_names, drop = FALSE]
  } else {
    matrix(FALSE, nrow = length(exposure_names_for_ldsc), ncol = length(outcome_names),
           dimnames = list(exposure_names_for_ldsc, outcome_names))
  }
  xy_rho_available_matrix <- if (!is.null(ldsc_cache$xy_rho_available)) {
    ldsc_cache$xy_rho_available[exposure_names_for_ldsc, outcome_names, drop = FALSE]
  } else {
    is.finite(xy_rho_matrix) & xy_overlap_matrix
  }
  
  cor_matrix <- build_ratio_cor_matrix(rg_matrix, nx, ny, outcome_names)
  
  data_x_ok <- apply(data_x_beta, 1, function(x) all(is.finite(x) & x != 0))
  data_x_se_ok <- apply(data_x_se, 1, function(x) all(is.finite(x) & x > 0))
  data_y_ok <- apply(data_y_beta, 1, function(x) all(is.finite(x)))
  data_y_se_ok <- apply(data_y_se, 1, function(x) all(is.finite(x) & x > 0))
  valid_indices <- which(data_x_ok & data_x_se_ok & data_y_ok & data_y_se_ok)
  log_step("valid SNPs after finite/se checks: ", length(valid_indices))
  
  if (length(valid_indices) < 3) {
    stop("Fewer than 3 valid SNPs remain for the pleiotropy-adjusted new method.")
  }
  
  original_inputs_tau0 <- build_wald_ratio_inputs(
    data_x_beta = data_x_beta,
    data_y_beta = data_y_beta,
    data_y_se = data_y_se,
    cor_matrix = cor_matrix,
    nx = nx,
    ny = ny,
    valid_indices = valid_indices
  )
  WR_original <- original_inputs_tau0$WR
  rownames(WR_original) <- colname
  colnames(WR_original) <- snp_ids[valid_indices]
  
  original_gls_summary <- t(sapply(seq_len(ncol(WR_original)), function(j) {
    collapse_gls_snp(as.numeric(WR_original[, j]), original_inputs_tau0$base_cov_list[[j]])
  }))
  
  original_dl_res <- tryCatch(
    calculate_tau2_DL(original_gls_summary[, "theta_hat"], original_gls_summary[, "se_hat"]),
    error = function(e) {
      stop("original tau2 DL failed: ", conditionMessage(e))
    }
  )
  original_tau_sq <- original_dl_res$tau2
  log_step("original tau2_DL finished; tau2=", signif(original_tau_sq, 4),
           ", Q=", signif(original_dl_res$Q, 4),
           ", p=", signif(original_dl_res$p_value, 4))
  
  original_inputs <- build_wald_ratio_inputs(
    data_x_beta = data_x_beta,
    data_y_beta = data_y_beta,
    data_y_se = data_y_se,
    cor_matrix = cor_matrix,
    nx = nx,
    ny = ny,
    valid_indices = valid_indices,
    tau2 = original_tau_sq
  )
  WR_original <- original_inputs$WR
  rownames(WR_original) <- colname
  colnames(WR_original) <- snp_ids[valid_indices]
  
  n_bootstrap <- suppressWarnings(as.integer(bootstrap_time))
  if (!is.finite(n_bootstrap) || n_bootstrap < 0) {
    n_bootstrap <- 0L
  }
  n_relaxed_nome_bootstrap <- if (is.null(relaxed_nome_bootstrap_time)) {
    n_bootstrap
  } else {
    suppressWarnings(as.integer(relaxed_nome_bootstrap_time))
  }
  if (!is.finite(n_relaxed_nome_bootstrap) || n_relaxed_nome_bootstrap < 0) {
    n_relaxed_nome_bootstrap <- 0L
  }
  
  theta_original_fit <- tryCatch(
    MLE_S(nx, ny, original_inputs$base_cov_list, WR_original, ncol(WR_original)),
    error = function(e) {
      warning("weak-1128原始MLE算法失败: ", conditionMessage(e))
      NA_real_
    }
  )
  theta_original <- as.numeric(theta_original_fit)
  theta_original_model_se <- suppressWarnings(as.numeric(attr(theta_original_fit, "se")))
  log_step("original theta estimated; theta=", signif(theta_original, 4))
  
  boot_effects_original <- numeric(0)
  if (n_bootstrap > 0) {
    log_step("original bootstrap started; B=", n_bootstrap)
    boot_effects_original <- numeric(n_bootstrap)
    for (k in seq_len(n_bootstrap)) {
      boot_cols <- sample(seq_len(ncol(WR_original)), ncol(WR_original), replace = TRUE)
      boot_effects_original[k] <- tryCatch({
        as.numeric(MLE_S(
          nx = nx,
          ny = ny,
          cov_list = original_inputs$base_cov_list[boot_cols],
          WR_matrix = WR_original[, boot_cols, drop = FALSE],
          SNPnew = ncol(WR_original)
        ))
      }, error = function(e) NA_real_)
      if (should_report_progress(k, n_bootstrap)) {
        log_step("original bootstrap progress: ", k, "/", n_bootstrap)
      }
    }
    log_step("original bootstrap finished")
  } else {
    log_step("original bootstrap skipped; using model-based SE")
  }
  boot_effects_original <- boot_effects_original[is.finite(boot_effects_original)]
  theta_original_se <- if (length(boot_effects_original) > 1) {
    sd(boot_effects_original)
  } else {
    theta_original_model_se
  }
  theta_original_p_value <- if (is.finite(theta_original) &&
                                is.finite(theta_original_se) &&
                                theta_original_se > 0) {
    2 * pnorm(-abs(theta_original / theta_original_se))
  } else {
    NA_real_
  }
  log_step("original SE/p-value finished; se=", signif(theta_original_se, 4),
           ", p=", signif(theta_original_p_value, 4))
  
  log_step("building pleiotropy-adjusted inputs")
  plugin_inputs <- build_plugin_adjusted_inputs(
    data_x_beta = data_x_beta,
    data_x_se = data_x_se,
    data_y_beta = data_y_beta,
    data_y_se = data_y_se,
    cor_matrix = cor_matrix,
    nx = nx,
    ny = ny,
    valid_indices = valid_indices,
    egger_indices = valid_indices
  )
  
  alpha_hat <- plugin_inputs$alpha_hat
  alpha_se <- plugin_inputs$alpha_se
  alpha_z <- rep(NA_real_, ny)
  alpha_p_value <- rep(NA_real_, ny)
  valid_alpha <- is.finite(alpha_hat) & is.finite(alpha_se) & alpha_se > 0
  alpha_z[valid_alpha] <- alpha_hat[valid_alpha] / alpha_se[valid_alpha]
  alpha_p_value[valid_alpha] <- 2 * pnorm(-abs(alpha_z[valid_alpha]))
  names(alpha_hat) <- outcome_names
  names(alpha_se) <- outcome_names
  names(alpha_z) <- outcome_names
  names(alpha_p_value) <- outcome_names
  
  pleiotropy_Q <- NA_real_
  pleiotropy_df <- sum(valid_alpha)
  pleiotropy_p_value <- NA_real_
  pleiotropy_significant <- NA
  if (pleiotropy_df > 0) {
    alpha_vec <- matrix(alpha_hat[valid_alpha], ncol = 1)
    V_alpha_valid <- plugin_inputs$V_alpha[valid_alpha, valid_alpha, drop = FALSE]
    pleiotropy_Q <- tryCatch(
      as.numeric(t(alpha_vec) %*% MASS::ginv(V_alpha_valid) %*% alpha_vec),
      error = function(e) NA_real_
    )
    if (is.finite(pleiotropy_Q)) {
      pleiotropy_p_value <- 1 - pchisq(pleiotropy_Q, df = pleiotropy_df)
      pleiotropy_significant <- pleiotropy_p_value < 0.05
    }
  }
  log_step("pleiotropy test finished; Q=", signif(pleiotropy_Q, 4),
           ", df=", pleiotropy_df,
           ", p=", signif(pleiotropy_p_value, 4),
           ", significant=", isTRUE(pleiotropy_significant))
  
  WR_matrix <- plugin_inputs$WR_adj
  seWR_matrix <- plugin_inputs$seWR
  rownames(WR_matrix) <- colname
  rownames(seWR_matrix) <- colname
  colnames(WR_matrix) <- snp_ids[valid_indices]
  colnames(seWR_matrix) <- snp_ids[valid_indices]
  
  SNPnew <- ncol(WR_matrix)
  
  ####estimate tau2
  gls_summary <- t(sapply(seq_len(SNPnew), function(j) {
    collapse_gls_snp(as.numeric(WR_matrix[, j]), plugin_inputs$base_cov_list[[j]])
  }))
  
  dl_res <- tryCatch(
    calculate_tau2_DL(gls_summary[, "theta_hat"], gls_summary[, "se_hat"]),
    error = function(e) {
      stop("pleiotropy-adjusted tau2 DL failed: ", conditionMessage(e))
    }
  )
  log_step("pleiotropy-adjusted tau2_DL finished; tau2=", signif(dl_res$tau2, 4),
           ", Q=", signif(dl_res$Q, 4),
           ", p=", signif(dl_res$p_value, 4))
  
  tau_sq <- dl_res$tau2
  
  est <- gls_theta_plugin(
    base_cov_list = plugin_inputs$base_cov_list,
    WR_matrix = WR_matrix,
    BA_list = plugin_inputs$BA_list,
    V_alpha = plugin_inputs$V_alpha,
    tau2 = tau_sq
  )
  
  theta_hat <- unname(est[["theta_hat"]])
  log_step("pleiotropy-adjusted theta estimated; theta=", signif(theta_hat, 4))
  
  boot_effects <- numeric(0)
  # Bootstrap计算标准误：每次重抽SNP并重新估计多效性截距
  if (n_bootstrap > 0) {
    log_step("pleiotropy-adjusted bootstrap started; B=", n_bootstrap)
    boot_effects <- numeric(n_bootstrap)
    for (k in seq_len(n_bootstrap)) {
      boot_sam <- sample(valid_indices, length(valid_indices), replace = TRUE)
      
      boot_effects[k] <- tryCatch({
        boot_inputs <- build_plugin_adjusted_inputs(
          data_x_beta = data_x_beta,
          data_x_se = data_x_se,
          data_y_beta = data_y_beta,
          data_y_se = data_y_se,
          cor_matrix = cor_matrix,
          nx = nx,
          ny = ny,
          valid_indices = boot_sam,
          egger_indices = boot_sam
        )
        
        gls_theta_plugin(
          base_cov_list = boot_inputs$base_cov_list,
          WR_matrix = boot_inputs$WR_adj,
          BA_list = boot_inputs$BA_list,
          V_alpha = boot_inputs$V_alpha,
          tau2 = tau_sq
        )[["theta_hat"]]
      }, error = function(e) NA_real_)
      if (should_report_progress(k, n_bootstrap)) {
        log_step("pleiotropy-adjusted bootstrap progress: ", k, "/", n_bootstrap)
      }
    }
    log_step("pleiotropy-adjusted bootstrap finished")
  } else {
    log_step("pleiotropy-adjusted bootstrap skipped; using GLS SE")
  }
  
  # 移除NA值
  boot_effects <- boot_effects[is.finite(boot_effects)]
  
  # 计算标准误和p值
  theta_se <- if (length(boot_effects) > 1) sd(boot_effects) else unname(est[["theta_se"]])
  z_score <- theta_hat / theta_se
  theta_p_value <- 2 * pnorm(-abs(z_score))
  log_step("pleiotropy-adjusted SE/p-value finished; se=", signif(theta_se, 4),
           ", p=", signif(theta_p_value, 4))
  
  relaxed_nome_adjusted <- NA_integer_
  relaxed_nome_tau2 <- NA_real_
  relaxed_nome_profile_se <- NA_real_
  relaxed_nome_bootstrap_success <- 0L
  relaxed_nome_se_method <- NA_character_
  relaxed_nome_est <- c(theta_hat = NA_real_, theta_se = NA_real_,
                        theta_p_value = NA_real_)
  if (isTRUE(run_relaxed_nome)) {
    relaxed_nome_adjusted <- isTRUE(pleiotropy_significant)
    relaxed_nome_alpha <- if (relaxed_nome_adjusted) alpha_hat else rep(0, ny)
    relaxed_nome_V_alpha <- if (relaxed_nome_adjusted) plugin_inputs$V_alpha else matrix(0, nrow = ny, ncol = ny)
    relaxed_nome_tau2 <- if (relaxed_nome_adjusted) tau_sq else original_tau_sq
    log_step("relaxed NOME profile started; pleiotropy_adjusted=",
             relaxed_nome_adjusted,
             ", tau2_used=", signif(relaxed_nome_tau2, 4))
    relaxed_nome_est <- tryCatch(
      estimate_relaxed_nome(
        beta_x = data_x_beta[valid_indices, , drop = FALSE],
        se_x = data_x_se[valid_indices, , drop = FALSE],
        beta_y = data_y_beta[valid_indices, , drop = FALSE],
        se_y = data_y_se[valid_indices, , drop = FALSE],
        alpha_hat = relaxed_nome_alpha,
        V_alpha = relaxed_nome_V_alpha,
        rho_y = as.matrix(rg_matrix),
        rho_x = as.matrix(rg_x_matrix),
        xy_rho = as.matrix(xy_rho_matrix),
        tau2 = relaxed_nome_tau2
      ),
      error = function(e) {
        warning("relaxed NOME profile estimator failed: ", conditionMessage(e))
        c(theta_hat = NA_real_, theta_se = NA_real_, theta_p_value = NA_real_)
      }
    )
    log_step("relaxed NOME profile finished; theta=",
             signif(unname(relaxed_nome_est[["theta_hat"]]), 4),
             ", se=", signif(unname(relaxed_nome_est[["theta_se"]]), 4),
             ", p=", signif(unname(relaxed_nome_est[["theta_p_value"]]), 4))
    relaxed_nome_profile_se <- unname(relaxed_nome_est[["theta_se"]])
    relaxed_nome_se_method <- if (n_relaxed_nome_bootstrap > 0) {
      "profile_hessian_fallback"
    } else {
      "profile_hessian_no_bootstrap"
    }
    
    relaxed_nome_boot_effects <- numeric(0)
    if (n_relaxed_nome_bootstrap > 0 && is.finite(unname(relaxed_nome_est[["theta_hat"]]))) {
      log_step("relaxed NOME bootstrap started; B=", n_relaxed_nome_bootstrap)
      relaxed_nome_boot_effects <- numeric(n_relaxed_nome_bootstrap)
      for (k in seq_len(n_relaxed_nome_bootstrap)) {
        boot_sam <- sample(valid_indices, length(valid_indices), replace = TRUE)
        
        relaxed_nome_boot_effects[k] <- tryCatch({
          boot_alpha <- relaxed_nome_alpha
          boot_V_alpha <- relaxed_nome_V_alpha
          if (isTRUE(relaxed_nome_adjusted)) {
            boot_inputs <- build_plugin_adjusted_inputs(
              data_x_beta = data_x_beta,
              data_x_se = data_x_se,
              data_y_beta = data_y_beta,
              data_y_se = data_y_se,
              cor_matrix = cor_matrix,
              nx = nx,
              ny = ny,
              valid_indices = boot_sam,
              egger_indices = boot_sam
            )
            boot_alpha <- boot_inputs$alpha_hat
            boot_V_alpha <- boot_inputs$V_alpha
          }
          
          unname(estimate_relaxed_nome(
            beta_x = data_x_beta[boot_sam, , drop = FALSE],
            se_x = data_x_se[boot_sam, , drop = FALSE],
            beta_y = data_y_beta[boot_sam, , drop = FALSE],
            se_y = data_y_se[boot_sam, , drop = FALSE],
            alpha_hat = boot_alpha,
            V_alpha = boot_V_alpha,
            rho_y = as.matrix(rg_matrix),
            rho_x = as.matrix(rg_x_matrix),
            xy_rho = as.matrix(xy_rho_matrix),
            tau2 = relaxed_nome_tau2
          )[["theta_hat"]])
        }, error = function(e) NA_real_)
        
        if (should_report_progress(k, n_relaxed_nome_bootstrap)) {
          log_step("relaxed NOME bootstrap progress: ", k, "/", n_relaxed_nome_bootstrap)
        }
      }
      relaxed_nome_boot_effects <- relaxed_nome_boot_effects[is.finite(relaxed_nome_boot_effects)]
      relaxed_nome_bootstrap_success <- length(relaxed_nome_boot_effects)
      log_step("relaxed NOME bootstrap finished; valid=", relaxed_nome_bootstrap_success,
               "/", n_relaxed_nome_bootstrap)
      
      if (relaxed_nome_bootstrap_success > 1) {
        relaxed_nome_boot_se <- sd(relaxed_nome_boot_effects)
        relaxed_nome_est[["theta_se"]] <- relaxed_nome_boot_se
        relaxed_nome_est[["theta_p_value"]] <- if (is.finite(relaxed_nome_boot_se) &&
                                                   relaxed_nome_boot_se > 0) {
          2 * pnorm(-abs(unname(relaxed_nome_est[["theta_hat"]]) / relaxed_nome_boot_se))
        } else {
          NA_real_
        }
        relaxed_nome_se_method <- "bootstrap"
        log_step("relaxed NOME bootstrap SE/p-value finished; se=",
                 signif(unname(relaxed_nome_est[["theta_se"]]), 4),
                 ", p=", signif(unname(relaxed_nome_est[["theta_p_value"]]), 4))
      } else {
        log_step("relaxed NOME bootstrap had too few valid replicates; using profile Hessian SE")
      }
    } else if (n_relaxed_nome_bootstrap > 0) {
      log_step("relaxed NOME bootstrap skipped because profile estimate is not finite")
    } else {
      log_step("relaxed NOME bootstrap skipped; using profile Hessian SE")
    }
  } else {
    log_step("relaxed NOME profile skipped")
  }
  
  # 返回结果
  return(list(
    theta_original = theta_original,
    theta_original_se = theta_original_se,
    theta_original_p_value = theta_original_p_value,
    theta_hat = theta_hat,
    theta_se = theta_se,
    theta_p_value = theta_p_value,
    theta_relaxed_nome = unname(relaxed_nome_est[["theta_hat"]]),
    theta_relaxed_nome_se = unname(relaxed_nome_est[["theta_se"]]),
    theta_relaxed_nome_p_value = unname(relaxed_nome_est[["theta_p_value"]]),
    theta_relaxed_nome_profile_se = relaxed_nome_profile_se,
    relaxed_nome_se_method = relaxed_nome_se_method,
    relaxed_nome_bootstrap_success = relaxed_nome_bootstrap_success,
    relaxed_nome_bootstrap_requested = n_relaxed_nome_bootstrap,
    relaxed_nome_pleiotropy_adjusted = as.integer(relaxed_nome_adjusted),
    relaxed_nome_tau2_used = relaxed_nome_tau2,
    SNPnew = SNPnew,
    cor_matrix =  rg_matrix,
    exposure_cor_matrix = rg_x_matrix,
    xy_rho_matrix = xy_rho_matrix,
    xy_overlap_matrix = xy_overlap_matrix,
    xy_rho_available_matrix = xy_rho_available_matrix,
    ldsc_cache = ldsc_cache,
    alpha_hat = alpha_hat,
    alpha_se = alpha_se,
    alpha_z = alpha_z,
    alpha_p_value = alpha_p_value,
    pleiotropy_Q = pleiotropy_Q,
    pleiotropy_df = pleiotropy_df,
    pleiotropy_p_value = pleiotropy_p_value,
    pleiotropy_significant = pleiotropy_significant,
    theta_original_DL_FE = original_dl_res$theta_fixed,
    tau2_original_DL = original_dl_res$tau2,
    I2_original_DL = original_dl_res$I2,
    original_DL_p_value = original_dl_res$p_value,
    original_DL_Q = original_dl_res$Q,
    original_DL_df = original_dl_res$df,
    tau_sq = tau_sq,
    theta0_DL_FE = dl_res$theta_fixed,
    I2_DL = dl_res$I2,
    DL_p_value = dl_res$p_value,
    DL_Q = dl_res$Q,
    DL_df = dl_res$df
  ))
}

############################## 实例分析 #################################
application_sample_size_map <- list(
  # Exposure sample sizes from the application data sheet.
  UKB_HDL = 400754,
  MVP_HDL = 404121,
  UKB_LDL = 431167,
  GLGC_LDL = 93982,
  UKB_TC = 437878,
  GLGC_TC = 93982,
  UKB_TG = 437532,
  GLGC_TG = 94595,
  # Outcome sample sizes.
  UKB_T2D = 468298,
  finn_T2D = 486367,
  MVP_T2D = 432648,
  EPIC_T2D = 22326,
  Swedish_T2D = 12230
)

get_application_N_list <- function(exposure_list = NULL, outcome_list = NULL) {
  dataset_names <- unique(c(names(exposure_list), names(outcome_list)))
  if (length(dataset_names) == 0) {
    return(application_sample_size_map)
  }
  out <- list()
  for (dataset_name in dataset_names) {
    if (!is.null(application_sample_size_map[[dataset_name]])) {
      out[[dataset_name]] <- application_sample_size_map[[dataset_name]]
    } else {
      warning("No application sample size was mapped for dataset: ", dataset_name)
    }
  }
  out
}

format_xy_overlap_pairs <- function(xy_overlap_matrix) {
  if (is.null(xy_overlap_matrix) || length(xy_overlap_matrix) == 0) {
    return("")
  }
  idx <- which(xy_overlap_matrix, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return("")
  }
  paste0(rownames(xy_overlap_matrix)[idx[, 1]], "-",
         colnames(xy_overlap_matrix)[idx[, 2]], collapse = ";")
}

format_xx_correlation_pairs <- function(exposure_cor_matrix, digits = 6) {
  if (is.null(exposure_cor_matrix) || length(exposure_cor_matrix) == 0) {
    return("")
  }
  exposure_cor_matrix <- as.matrix(exposure_cor_matrix)
  if (nrow(exposure_cor_matrix) < 2 || ncol(exposure_cor_matrix) < 2) {
    return("")
  }
  idx <- which(upper.tri(exposure_cor_matrix), arr.ind = TRUE)
  vals <- exposure_cor_matrix[idx]
  keep <- is.finite(vals)
  if (!any(keep)) {
    return("")
  }
  idx <- idx[keep, , drop = FALSE]
  vals <- vals[keep]
  paste0(rownames(exposure_cor_matrix)[idx[, 1]], "-",
         colnames(exposure_cor_matrix)[idx[, 2]], ":",
         formatC(vals, digits = digits, format = "fg"), collapse = ";")
}

format_xy_correlation_pairs <- function(xy_rho_matrix, xy_overlap_matrix = NULL,
                                        xy_rho_available_matrix = NULL,
                                        digits = 6) {
  if (is.null(xy_rho_matrix) || length(xy_rho_matrix) == 0) {
    return("")
  }
  xy_rho_matrix <- as.matrix(xy_rho_matrix)
  if (is.null(xy_overlap_matrix)) {
    xy_overlap_matrix <- matrix(TRUE, nrow = nrow(xy_rho_matrix),
                                ncol = ncol(xy_rho_matrix),
                                dimnames = dimnames(xy_rho_matrix))
  } else {
    xy_overlap_matrix <- as.matrix(xy_overlap_matrix)
  }
  if (is.null(xy_rho_available_matrix)) {
    xy_rho_available_matrix <- matrix(TRUE, nrow = nrow(xy_rho_matrix),
                                      ncol = ncol(xy_rho_matrix),
                                      dimnames = dimnames(xy_rho_matrix))
  } else {
    xy_rho_available_matrix <- as.matrix(xy_rho_available_matrix)
  }
  idx <- which(xy_overlap_matrix, arr.ind = TRUE)
  if (nrow(idx) == 0) {
    return("")
  }
  vals <- vapply(seq_len(nrow(idx)), function(k) {
    value <- xy_rho_matrix[idx[k, 1], idx[k, 2]]
    available <- isTRUE(xy_rho_available_matrix[idx[k, 1], idx[k, 2]])
    if (!available || !is.finite(value)) {
      return("NA")
    }
    formatC(value, digits = digits, format = "fg")
  }, character(1))
  paste0(rownames(xy_rho_matrix)[idx[, 1]], "-",
         colnames(xy_rho_matrix)[idx[, 2]], ":",
         vals, collapse = ";")
}

RUN_APPLICATION_MAIN <- if (exists("RUN_APPLICATION_MAIN", inherits = TRUE)) {
  get("RUN_APPLICATION_MAIN", inherits = TRUE)
} else {
  TRUE
}

if (isTRUE(RUN_APPLICATION_MAIN)) {
  exposure_names <- c("HDL22_list", 
                      "LDL2_list",  
                      "TC2_list",
                      "TG2_list" ##UKB
  )
  
  outcome_names <- c("T2D5_list")
  
  base_path <- "/data/wusijia/"
  clump_kb = 1000
  clump_r2 = 0.001
  JointMR_F_threshold <- 10
  
  for (outcome_name in outcome_names) {
    
    load(paste0(base_path,"T2D/",outcome_name,".Rdata"))
    
    for (exposure_name in exposure_names) {
      
      ##### 加载数据集####
      ex_file <- sub("\\d.*", "", exposure_name)
      load(paste0(base_path,ex_file,"/",exposure_name,".Rdata"))
      
      # 生成组合名称
      combo_name <- paste0(gsub("_list", "", exposure_name), "_", 
                           gsub("_list", "", outcome_name))
      cat("开始分析:", combo_name, "\n")
      
      nx = length(exposure_list)
      ny = length(outcome_list)
      
      gc()
      
      ##### 分析部分####
      ##### MR-Meta暴露X处理,求并集#####
      results_exposure <- process_gwas_data(exposure_list,
                                            ncore = 8,  # 使用8个核心
                                            pval_threshold = 5e-8,  # GWAS显著阈值
                                            clump_kb = clump_kb,
                                            clump_r2 = clump_r2,
                                            F_threshold = 10)  # 传统方法弱工具变量过滤阈值)
      union_snps <- results_exposure$union_snps
      cat("整合后的SNP数量:", nrow(union_snps), "\n")
      
      ##### 执行传统MR-Meta分析####
      SNP_list <- results_exposure$clumped_results
      
      results_traditional_mr_meta <- run_mr_analysis_no_ref(
        exposure_list = exposure_list,
        outcome_list = outcome_list,
        SNP_list = SNP_list,
        nx = length(exposure_list),
        ny = length(outcome_list)
      )
      
      res_divw <- results_traditional_mr_meta$res_divw
      res_fixed <- results_traditional_mr_meta$res_fixed
      res_random <- results_traditional_mr_meta$res_random
      res_wme <- results_traditional_mr_meta$res_wme
      SNPSMR <- results_traditional_mr_meta$SNPSMR
      
      rere2 <- c(res_divw[,1],res_fixed[,1],res_random[,1],res_wme[,1],
                 res_divw[,2],res_fixed[,2],res_random[,2],res_wme[,2],
                 res_divw[,3],res_fixed[,3],res_random[,3],res_wme[,3])
      
      meta_old_divw <- try(metagen(res_divw[,1],res_divw[,2],sm="MD"), silent=TRUE)
      meta_old_ivw_fixed <- try(metagen(res_fixed[,1],res_fixed[,2],sm="MD"), silent=TRUE)
      meta_old_ivw_random <- try(metagen(res_random[,1],res_random[,2],sm="MD"), silent=TRUE)
      meta_old_wme <- try(metagen(res_wme[,1],res_wme[,2],sm="MD"), silent=TRUE)
      
      if(meta_old_divw$pval.Q<0.05){
        meta_old_divw1 <- c(meta_old_divw$TE.random,meta_old_divw$seTE.random,meta_old_divw$pval.random)
      }else{
        meta_old_divw1 <- c(meta_old_divw$TE.fixed,meta_old_divw$seTE.fixed,meta_old_divw$pval.fixed)
      }
      
      if(meta_old_ivw_fixed$pval.Q<0.05){
        meta_old_ivw_fixed1 <- c(meta_old_ivw_fixed$TE.random,meta_old_ivw_fixed$seTE.random,meta_old_ivw_fixed$pval.random)
      }else{
        meta_old_ivw_fixed1 <- c(meta_old_ivw_fixed$TE.fixed,meta_old_ivw_fixed$seTE.fixed,meta_old_ivw_fixed$pval.fixed)
      }
      
      if(meta_old_ivw_random$pval.Q<0.05){
        meta_old_ivw_random1 <- c(meta_old_ivw_random$TE.random,meta_old_ivw_random$seTE.random,meta_old_ivw_random$pval.random)
      }else{
        meta_old_ivw_random1 <- c(meta_old_ivw_random$TE.fixed,meta_old_ivw_random$seTE.fixed,meta_old_ivw_random$pval.fixed)
      }
      
      
      if(meta_old_wme$pval.Q<0.05){
        meta_old_wme1 <- c(meta_old_wme$TE.random,meta_old_wme$seTE.random,meta_old_wme$pval.random)
      }else{
        meta_old_wme1 <- c(meta_old_wme$TE.fixed,meta_old_wme$seTE.fixed,meta_old_wme$pval.fixed)
      }
      
      rere3 <- c(meta_old_divw1,meta_old_ivw_fixed1,meta_old_ivw_random1,meta_old_wme1)
      
      ##### 协调暴露数据的等位基因####
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
      
      ##### 进行暴露的METAL分析####
      dir.create(paste0(combo_dir,"/",ex_file), recursive = TRUE)
      setwd(paste0(combo_dir,"/",ex_file))
      metal_script <- paste0(combo_dir,"/",ex_file,"/metal_script.txt")
      
      clean_script <- function(lines) {
        # 移除首尾空格
        lines <- trimws(lines)
        # 移除空行
        lines <- lines[nchar(lines) > 0]
        return(lines)
      }
      
      # 生成PROCESS行
      dataset_names <- names(exposure_list)
      process_lines <- sapply(dataset_names, function(name) {
        file_path <- file.path(base_path, ex_file, paste0(name, ".txt"))
        paste0("PROCESS ", file_path)
      })
      
      script_lines <- clean_script(c(
        "SCHEME STDERR",
        "MARKER SNP",
        "ALLELE effect_allele other_allele",
        "EFFECT beta",
        "STDERR se",
        "PVAL p",
        process_lines,
        "ANALYZE"
      ))
      
      writeLines(script_lines, con = metal_script)
      
      # 2. 运行METAL并捕获输出
      output_log <- system2(
        command = "/data/wusijia/generic-metal/metal",
        args = metal_script,
        stdout = TRUE,
        stderr = TRUE
      )
      
      
      ##### 读入meta-x结果####
      alpha_data <- process_metal_results(file_path = paste0(combo_dir,"/",ex_file,"/METAANALYSIS1.TBL"),
                                          clump_kb = clump_kb,
                                          clump_r2 = clump_r2,
                                          F_threshold = 10)
      
      a_clump_alpha_data1 <- alpha_data$clumped_data
      GWAS_Meta_exposure <- alpha_data$metal_data[alpha_data$metal_data$SNP%in%a_clump_alpha_data1$rsid,]
      
      ##### 筛选结局，生成结局txt####
      dir.create(paste0(combo_dir,"/T2D"), recursive = TRUE)
      results <- harmonize_and_export_general(
        outcome_list = outcome_list,
        global_ref = global_ref,
        snp_list = a_clump_alpha_data1,
        output_dir = paste0(combo_dir,"/T2D"),
        n_cores = 4
      )
      
      ##### 进行outcome METAL分析####
      setwd(paste0(combo_dir,"/T2D"))
      metal_script2 <- paste0(combo_dir,"/T2D/metal_script2.txt")
      
      clean_script <- function(lines) {
        # 移除首尾空格
        lines <- trimws(lines)
        # 移除空行
        lines <- lines[nchar(lines) > 0]
        return(lines)
      }
      
      dataset_names2 <- names(outcome_list)
      process_lines2 <- sapply(seq_along(dataset_names2), function(i) {
        file_path <- file.path(combo_dir, "T2D", paste0("dataset_", i, ".txt"))
        paste0("PROCESS ", file_path)
      })
      
      script_lines <- clean_script(c(
        "SCHEME STDERR",
        "MARKER SNP",
        "ALLELE effect_allele other_allele",
        "EFFECT beta",
        "STDERR se",
        "PVAL p",
        process_lines2,
        "ANALYZE"
      ))
      
      writeLines(script_lines, con = metal_script2)
      
      # 2. 运行METAL并捕获输出
      output_log <- system2(
        command = "/data/wusijia/generic-metal/metal",
        args = metal_script2,
        stdout = TRUE,
        stderr = TRUE
      )
      
      ##### 执行GWAS-meta-MR分析####
      GWAS_Meta_outcome <- fread(paste0(combo_dir,"/T2D/METAANALYSIS1.TBL"),
                                 sep = "\t",
                                 header = TRUE,
                                 fill = TRUE)
      names(GWAS_Meta_outcome)[names(GWAS_Meta_outcome) == "P-value"] <-"pval"
      names(GWAS_Meta_outcome)[names(GWAS_Meta_outcome) == "MarkerName"] <-"SNP"
      GWAS_Meta_outcome <- subset(GWAS_Meta_outcome, SNP %in% GWAS_Meta_exposure$SNP)
      
      SNP_GWAS_meta <- Reduce(intersect, list(GWAS_Meta_exposure$SNP, GWAS_Meta_outcome$SNP))
      GWAS_Meta_exposure <- subset(GWAS_Meta_exposure, SNP %in% SNP_GWAS_meta)
      GWAS_Meta_outcome <- subset(GWAS_Meta_outcome, SNP %in% SNP_GWAS_meta)
      
      GWAS_Meta_exposure <- GWAS_Meta_exposure[order(GWAS_Meta_exposure$SNP), ]
      GWAS_Meta_outcome <- GWAS_Meta_outcome[order(GWAS_Meta_outcome$SNP), ]
      
      # 确保SNP顺序一致
      if (!identical(GWAS_Meta_exposure$SNP, GWAS_Meta_outcome$SNP)) {
        # 如果排序后仍然不一致，使用match函数强制对齐
        GWAS_Meta_outcome <- GWAS_Meta_outcome[match(GWAS_Meta_exposure$SNP, GWAS_Meta_outcome$SNP), ]
      }
      
      ####MR
      mr_object_fix <- mr_input(bx = GWAS_Meta_exposure$beta,
                                bxse = GWAS_Meta_exposure$StdErr,
                                by = GWAS_Meta_outcome$Effect,
                                byse = GWAS_Meta_outcome$StdErr)##fix
      
      ###mr.divw
      r_divw <- mr.divw(GWAS_Meta_exposure$beta,GWAS_Meta_outcome$Effect,GWAS_Meta_exposure$StdErr,GWAS_Meta_outcome$StdErr)
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
      
      ##### New method #####
      # XX/YY/XY LDSC需要每个暴露和结局数据库的样本量。
      # 暴露样本量来自本应用数据表，并按 exposure_list/outcome_list 的数据集名称匹配。
      N_list <- get_application_N_list(exposure_list, outcome_list)
      
      pack_jointmr_grid <- function(re) {
        c(SNPnew = re$SNPnew,
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
                                                              re$DL_p_value < 0.05))
      }
      
      f_selection_grid <- c(strictF = "all", relaxedF = "any")
      p_combine_grid <- "stouffer"
      jointmr_grid_values <- c()
      default_re <- NULL
      default_relaxed_re <- NULL
      ldsc_cache <- NULL
      
      for (f_label in names(f_selection_grid)) {
        for (p_method in p_combine_grid) {
          cat("Running JointMR grid:", f_label, "/", p_method, "\n")
          results_Meta_MR_grid <- full_dataset_processing(
            exposure_list = exposure_list,
            outcome_list = outcome_list,
            clump_kb = clump_kb,
            clump_r2 = clump_r2,
            F_threshold = JointMR_F_threshold,
            f_selection = f_selection_grid[[f_label]],
            p_combine_method = p_method
          )
          
          re_grid <- run_new_mr_analysis(
            exposure_filter_list = results_Meta_MR_grid$processed_exposure,
            outcome_filter_list = results_Meta_MR_grid$processed_outcome,
            original_exposure_list = exposure_list,
            original_outcome_list = outcome_list,
            N_list = N_list,
            ancestry = "EUR",
            bootstrap_time = 500,
            relaxed_nome_bootstrap_time = 200,
            ldsc_cache = ldsc_cache
          )
          if (is.null(ldsc_cache)) {
            ldsc_cache <- re_grid$ldsc_cache
          }
          
          grid_values <- pack_jointmr_grid(re_grid)
          names(grid_values) <- paste0("JointMR_", f_label, "_", p_method, "_", names(grid_values))
          jointmr_grid_values <- c(jointmr_grid_values, grid_values)
          
          if (f_label == "strictF" && p_method == "stouffer") {
            new_re <- re_grid
          }
          if (f_label == "relaxedF" && p_method == "stouffer") {
            new_re_relaxed <- re_grid
          }
        }
      }
      
      theta_original <- new_re$theta_original
      theta_original_se <- new_re$theta_original_se
      theta_original_p_value <- new_re$theta_original_p_value
      theta_hat <- new_re$theta_hat
      theta_se <- new_re$theta_se
      theta_p_value <- new_re$theta_p_value
      SNPnew <- new_re$SNPnew
      pleiotropy_global <- c(new_re$pleiotropy_Q,
                             new_re$pleiotropy_df,
                             new_re$pleiotropy_p_value,
                             as.integer(isTRUE(new_re$pleiotropy_significant)))
      pleiotropy_by_outcome <- c(new_re$alpha_hat,
                                 new_re$alpha_se,
                                 new_re$alpha_p_value)
      heterogeneity_original <- c(new_re$tau2_original_DL,
                                  new_re$I2_original_DL,
                                  new_re$original_DL_Q,
                                  new_re$original_DL_df,
                                  new_re$original_DL_p_value,
                                  as.integer(is.finite(new_re$original_DL_p_value) &&
                                               new_re$original_DL_p_value < 0.05))
      heterogeneity_pleiotropy <- c(new_re$tau_sq,
                                    new_re$I2_DL,
                                    new_re$DL_Q,
                                    new_re$DL_df,
                                    new_re$DL_p_value,
                                    as.integer(is.finite(new_re$DL_p_value) &&
                                                 new_re$DL_p_value < 0.05))
      theta_original_relaxed <- new_re_relaxed$theta_original
      theta_original_se_relaxed <- new_re_relaxed$theta_original_se
      theta_original_p_value_relaxed <- new_re_relaxed$theta_original_p_value
      theta_hat_relaxed <- new_re_relaxed$theta_hat
      theta_se_relaxed <- new_re_relaxed$theta_se
      theta_p_value_relaxed <- new_re_relaxed$theta_p_value
      SNPnew_relaxed <- new_re_relaxed$SNPnew
      pleiotropy_global_relaxed <- c(new_re_relaxed$pleiotropy_Q,
                                     new_re_relaxed$pleiotropy_df,
                                     new_re_relaxed$pleiotropy_p_value,
                                     as.integer(isTRUE(new_re_relaxed$pleiotropy_significant)))
      pleiotropy_by_outcome_relaxed <- c(new_re_relaxed$alpha_hat,
                                         new_re_relaxed$alpha_se,
                                         new_re_relaxed$alpha_p_value)
      heterogeneity_original_relaxed <- c(new_re_relaxed$tau2_original_DL,
                                          new_re_relaxed$I2_original_DL,
                                          new_re_relaxed$original_DL_Q,
                                          new_re_relaxed$original_DL_df,
                                          new_re_relaxed$original_DL_p_value,
                                          as.integer(is.finite(new_re_relaxed$original_DL_p_value) &&
                                                       new_re_relaxed$original_DL_p_value < 0.05))
      heterogeneity_pleiotropy_relaxed <- c(new_re_relaxed$tau_sq,
                                            new_re_relaxed$I2_DL,
                                            new_re_relaxed$DL_Q,
                                            new_re_relaxed$DL_df,
                                            new_re_relaxed$DL_p_value,
                                            as.integer(is.finite(new_re_relaxed$DL_p_value) &&
                                                         new_re_relaxed$DL_p_value < 0.05))
      # theta_hat
      # theta_se
      # theta_p_value
      # new_re$cor_matrix
      
      ####存储结果####
      res_all <- c(clump_kb,clump_r2,nx,ny,SNPSMR,SNPnew,SNPnew_relaxed,
                   length(SNP_GWAS_meta),names(exposure_list),
                   names(outcome_list),rere1,rere2,rere3,
                   theta_original,theta_original_se,theta_original_p_value,
                   theta_hat,theta_se,theta_p_value,
                   theta_original_relaxed,theta_original_se_relaxed,theta_original_p_value_relaxed,
                   theta_hat_relaxed,theta_se_relaxed,theta_p_value_relaxed,
                   heterogeneity_original,heterogeneity_pleiotropy,
                   heterogeneity_original_relaxed,heterogeneity_pleiotropy_relaxed,
                   jointmr_grid_values,
                   pleiotropy_global,pleiotropy_by_outcome,
                   pleiotropy_global_relaxed,pleiotropy_by_outcome_relaxed) 
      res_all <- matrix(res_all ,nrow=1)
      res_all <- data.frame(res_all)
      
      ####colnames####
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
      
      SNPSMR <- NULL
      for (i in 1:(nx * ny)) {
        var_name <- paste0("SNPSMR_", i)
        SNPSMR <- c(SNPSMR,var_name)
      }
      
      ccname <- c('clump_kb','clump_r2','nx','ny',SNPSMR,'SNPnew','SNPnew_relaxed','SNP_GWAS_meta',
                  names(exposure_list),
                  names(outcome_list),
                  paste0('meta_res_dIVW_',c('beta','se','pval')),
                  paste0('meta_res_IVW_fixed_',c('beta','se','pval')),
                  paste0('meta_res_IVW_random_',c('beta','se','pval')),
                  paste0('meta_res_wme_',c('beta','se','pval')),
                  
                  all_beta_names,all_se_names,all_pval_names,
                  
                  paste0('MR_meta_old_dIVW_',c('beta','se','pval')),
                  paste0('MR_meta_old_ivw_fixed_',c('beta','se','pval')),
                  paste0('MR_meta_old_ivw_random_',c('beta','se','pval')),
                  paste0('MR_meta_old_wme_',c('beta','se','pval')),
                  paste0('MRmeta_original_',c('beta','se','pval')),
                  paste0('MRmeta_',c('beta','se','pval')),
                  paste0('MRmeta_original_relaxed_',c('beta','se','pval')),
                  paste0('MRmeta_relaxed_',c('beta','se','pval')),
                  paste0('heterogeneity_original_',c('tau2_DL','I2_DL','Q','df','pval','significant')),
                  paste0('heterogeneity_pleiotropy_',c('tau2_used','I2_DL','Q','df','pval','significant')),
                  paste0('heterogeneity_original_relaxed_',c('tau2_DL','I2_DL','Q','df','pval','significant')),
                  paste0('heterogeneity_pleiotropy_relaxed_',c('tau2_used','I2_DL','Q','df','pval','significant')),
                  names(jointmr_grid_values),
                  paste0('pleiotropy_global_',c('Q','df','pval','significant')),
                  paste0('pleiotropy_alpha_',names(outcome_list)),
                  paste0('pleiotropy_alpha_se_',names(outcome_list)),
                  paste0('pleiotropy_alpha_pval_',names(outcome_list)),
                  paste0('pleiotropy_global_relaxed_',c('Q','df','pval','significant')),
                  paste0('pleiotropy_alpha_relaxed_',names(outcome_list)),
                  paste0('pleiotropy_alpha_se_relaxed_',names(outcome_list)),
                  paste0('pleiotropy_alpha_pval_relaxed_',names(outcome_list)))
      colnames(res_all) <- ccname
      write.csv(res_all,paste0(base_path,"0506/",combo_name,"_re_0508_new_0.001_500.csv"))
      rm(exposure_list)
      gc()
    }
  }
}
