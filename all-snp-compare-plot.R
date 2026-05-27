source('F://wsj-doctor//2021-跨种族MR//ago//20230927--no pleo//QQ.plot_20230926.R')
source('F://wsj-doctor//2021-跨种族MR//ago//20230927--no pleo//histplot_20230926.R')

library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
library(patchwork)
library(outliers)
library(dplyr)

######### all valid SNPs, nx=3, ny=5, fixed-effect scenario #########
setwd("F://wsj-doctor//2024-MR-META//genemo一修//simulation addition//all_snps")

gb <- c(100)
beta0b <- c(0, 0.05, 0.1)
taub <- c(0)
rho_1b <- c(0.5, 0.8, 0.85, 0.9, 0.95)
rho_2b <- c(0.5, 0.8, 0.85, 0.9, 0.95)
siga1b <- c(0)
siga2b <- c(0)
nxb <- c(3)
nyb <- c(5)
weakp <- 0

tiquname <- c('nx', 'ny', 'g', 'beta0', "rho_mean", 'tau', "siga1",
              "siga2", "SNP_gaws_meta", "SNPMR", "SNPnew",
              paste0('meta_res_dIVW_', c('beta', 'se', 'pval')),
              paste0('meta_res_IVW_fixed_', c('beta', 'se', 'pval')),
              paste0('meta_res_IVW_random_', c('beta', 'se', 'pval')),
              paste0('meta_res_wme_', c('beta', 'se', 'pval')),
              paste0('MR_meta_old_dIVW_', c('beta', 'se', 'pval')),
              paste0('MR_meta_old_ivw_fixed_', c('beta', 'se', 'pval')),
              paste0('MR_meta_old_ivw_random_', c('beta', 'se', 'pval')),
              paste0('MR_meta_old_wme_', c('beta', 'se', 'pval')),
              paste0('MRmeta_', c('beta', 'se', 'pval')),
              "datanum", "F_tat", "weakp")

make_colnames <- function(nx, ny) {
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

  c('nx', 'ny', 'g', 'beta0', "rho_mean", 'tau', "siga1",
    "siga2", "SNP_gaws_meta", "SNPMR", "SNPnew",
    paste0('meta_res_dIVW_', c('beta', 'se', 'pval')),
    paste0('meta_res_IVW_fixed_', c('beta', 'se', 'pval')),
    paste0('meta_res_IVW_random_', c('beta', 'se', 'pval')),
    paste0('meta_res_wme_', c('beta', 'se', 'pval')),
    all_beta_names, all_se_names, all_pval_names,
    paste0('MR_meta_old_dIVW_', c('beta', 'se', 'pval')),
    paste0('MR_meta_old_ivw_fixed_', c('beta', 'se', 'pval')),
    paste0('MR_meta_old_ivw_random_', c('beta', 'se', 'pval')),
    paste0('MR_meta_old_wme_', c('beta', 'se', 'pval')),
    paste0('MRmeta_', c('beta', 'se', 'pval')), "F_tat", "weakp")
}

read_one_result <- function(nx, ny, g, beta0, rho_1, tau, siga1, siga2) {
  ccname <- make_colnames(nx, ny)
  filename <- paste0("Sum_nx_", nx, "_ny_", ny, "_g_", g,
                     "_beta0_", beta0, "_rho1_", rho_1,
                     "_tau_", tau, "_siga1_", siga1,
                     "_siga2_", siga2, "_all_valid_snps_1128.csv")

  if (!file.exists(filename)) {
    stop("Missing input file: ", filename)
  }

  once <- read.csv(filename, header = FALSE, stringsAsFactors = FALSE)
  once <- once[, seq_along(ccname), drop = FALSE]
  colnames(once) <- ccname

  # write_csv(..., append=TRUE) may leave header-like rows depending on how files were created.
  once$nx_numeric <- suppressWarnings(as.numeric(once$nx))
  once <- once[!is.na(once$nx_numeric), ]
  once$nx_numeric <- NULL

  once <- na.omit(once)
  once$datanum <- paste0(once$nx, "-", once$ny)
  once[tiquname]
}

result_all <- NULL
for (yh in seq_along(taub)) {
  tau <- taub[[yh]]
  for (th in seq_along(gb)) {
    g <- gb[[th]]
    for (xh in seq_along(siga1b)) {
      siga1 <- siga1b[[xh]]
      siga2 <- siga2b[[xh]]
      for (uh in seq_along(beta0b)) {
        beta0 <- beta0b[[uh]]
        for (rh in seq_along(rho_1b)) {
          rho_1 <- rho_1b[[rh]]
          rho_2 <- rho_2b[[rh]]
          for (zh in seq_along(nxb)) {
            nx <- nxb[[zh]]
            for (qh in seq_along(nyb)) {
              ny <- nyb[[qh]]
              once <- read_one_result(nx, ny, g, beta0, rho_1,
                                      tau, siga1, siga2)
              result_all <- rbind(result_all, once)
            }
          }
        }
      }
    }
  }
}

ress <- data.frame(beta = c(result_all$meta_res_dIVW_beta,
                            result_all$meta_res_IVW_fixed_beta,
                            result_all$meta_res_IVW_random_beta,
                            result_all$meta_res_wme_beta,
                            result_all$MR_meta_old_dIVW_beta,
                            result_all$MR_meta_old_ivw_fixed_beta,
                            result_all$MR_meta_old_ivw_random_beta,
                            result_all$MR_meta_old_wme_beta,
                            result_all$MRmeta_beta),
                   pval = c(result_all$meta_res_dIVW_pval,
                            result_all$meta_res_IVW_fixed_pval,
                            result_all$meta_res_IVW_random_pval,
                            result_all$meta_res_wme_pval,
                            result_all$MR_meta_old_dIVW_pval,
                            result_all$MR_meta_old_ivw_fixed_pval,
                            result_all$MR_meta_old_ivw_random_pval,
                            result_all$MR_meta_old_wme_pval,
                            result_all$MRmeta_pval),
                   beta1 = rep(result_all$beta0, 9),
                   rho = rep(result_all$rho_mean, 9),
                   tau = rep(result_all$tau, 9),
                   datanum = rep(result_all$datanum, 9),
                   weakp = rep(result_all$weakp, 9),
                   methods = c(rep(c('meta_dIVW', 'meta_FIVW', 'meta_RIVW',
                                     'meta_WME', "dIVW_meta", 'FIVW_meta',
                                     'RIVW_meta', 'WME_meta', 'JointMR'),
                                   each = nrow(result_all))))

ress$beta <- as.numeric(ress$beta)
ress$pval <- as.numeric(ress$pval)
ress$beta1 <- as.numeric(ress$beta1)
ress$rho <- as.numeric(ress$rho)
ress$tau <- as.numeric(ress$tau)

ress$beta0 <- factor(ress$beta1,
                     levels = beta0b,
                     ordered = TRUE,
                     labels = c("causal effect = 0.00",
                                "causal effect = 0.05",
                                "causal effect = 0.10"))
ress$rho_f <- factor(ress$rho,
                     levels = rho_1b,
                     ordered = TRUE,
                     labels = c("0.5", "0.8", "0.85", "0.9", "0.95"))
ress$methods <- factor(ress$methods,
                       levels = c('meta_dIVW', 'meta_FIVW', 'meta_RIVW',
                                  'meta_WME', "dIVW_meta", 'FIVW_meta',
                                  'RIVW_meta', 'WME_meta', 'JointMR'))

coloo <- c(brewer.pal(9, "Set3"))[c(1:3, 5:9, 4)]

####################### boxplot #################################
box_y_breaks <- sort(unique(c(beta0b - 0.1, beta0b, beta0b + 0.1)))
effect_line_data <- data.frame(
  beta0 = factor(levels(ress$beta0),
                 levels = levels(ress$beta0),
                 ordered = TRUE),
  beta1 = beta0b
)

p_box <- ggplot(data = ress, aes(x = rho_f, y = beta)) +
  geom_boxplot(aes(fill = methods), outlier.alpha = 0.5) +
  scale_fill_manual(values = coloo) +
  geom_hline(data = effect_line_data,
             aes(yintercept = beta1),
             linetype = "dashed", color = "darkred") +
  theme_bw() +
  scale_y_continuous(limits = c(min(beta0b - 0.1), max(beta0b + 0.1)),
                     breaks = box_y_breaks,
                     labels = sprintf("%.2f", box_y_breaks)) +
  xlab("dataset correlation") +
  ylab("Estimation") +
  labs(tag = "(A)") +
  facet_grid(beta0 ~ ., scales = "fixed") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        strip.text.y = element_text(size = 30),
        plot.title = element_text(size = 25, hjust = 0.5, vjust = 0.5),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 30),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        plot.tag = element_text(size = 40),
        legend.spacing = unit(2, "cm"),
        legend.key.size = unit(3, "lines")) +
  guides(fill = guide_legend(title = "Method"))

ggsave("all_snp_fixed_3-5_Estimation_boxplot_all_effects.pdf",
       plot = p_box, width = 30, height = 18)

####################### Q-Q plot #################################
method_names <- unique(ress$methods)
rho_meanb <- rho_1b
rho_meanname <- c("data correlation is 0.5",
                  "data correlation is 0.8",
                  "data correlation is 0.85",
                  "data correlation is 0.9",
                  "data correlation is 0.95")

for (m in 1:1) {
  plot_list <- list()
  for (n in seq_along(rho_meanb)) {
    once <- ress[(ress$beta1 == beta0b[m] &
                    ress$rho == rho_meanb[n]), ]

    p_value <- list()
    for (i in seq_along(method_names)) {
      p_value[[i]] <- as.numeric(once$pval[once$methods == method_names[i]])
    }

    plot <- qq.plot(p_value, legend_label = method_names) +
      labs(title = paste0(rho_meanname[n])) +
      theme(plot.title = element_text(size = 25),
            plot.tag = element_text(size = 23),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20))

    plot_list[[length(plot_list) + 1]] <- plot
  }

  final_plot <- wrap_plots(plot_list, ncol = length(rho_meanb)) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right",
          legend.spacing = unit(1.5, "cm"),
          legend.key.size = unit(3, "lines"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20))
  final_plot <- ggdraw(final_plot) +
    draw_label("(B)", x = 0.005, y = 0.995,
               hjust = 0, vjust = 1, size = 40)

  pdf(paste0("all_snp_fixed_3-5_qq.pdf"), width = 30, height = 7)
  print(final_plot)
  dev.off()
}

################# power plot ########################
thres_df <- expand.grid(methods = method_names, rho = rho_meanb)
thres_df$thres <- NA

for (r in rho_meanb) {
  once_H0 <- ress[ress$beta1 == beta0b[1] &
                    ress$rho == r, ]
  for (mth in method_names) {
    pvals <- once_H0$pval[once_H0$methods == mth]
    if (length(pvals) > 0) {
      pvals_sorted <- sort(pvals)
      index <- floor(0.05 * length(pvals_sorted)) + 1
      thres_value <- if (index > length(pvals_sorted)) {
        pvals_sorted[length(pvals_sorted)]
      } else {
        pvals_sorted[index]
      }
      thres_df$thres[thres_df$methods == mth &
                       thres_df$rho == r] <- thres_value
    } else {
      stop(sprintf("No H0 data for method %s, rho=%s", mth, r))
    }
  }
}

power_effect_rows <- list()
power_effect_names <- c("causal effect = 0.05",
                        "causal effect = 0.10")

for (m in 2:length(beta0b)) {
  plot_list <- list()

  for (n in seq_along(rho_meanb)) {
    once <- ress[(ress$beta1 == beta0b[m] &
                    ress$rho == rho_meanb[n]), ]

    powt <- NULL
    for (i in seq_along(method_names)) {
      once1 <- once[once$methods == method_names[i], ]
      thres_value <- thres_df$thres[thres_df$methods == method_names[i] &
                                      thres_df$rho == rho_meanb[n]]
      pow <- round(mean(once1$pval < thres_value), digits = 2)
      powt <- c(powt, pow)
    }

    plot <- hist.plot(power = powt,
                      legend_label = method_names) +
      labs(title = paste0(rho_meanname[n])) +
      theme(plot.title = element_text(size = 25),
            plot.tag = element_text(size = 23),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            legend.position = "none")

    plot_list[[length(plot_list) + 1]] <- plot
  }

  big_plot <- wrap_plots(plot_list, ncol = length(rho_meanb)) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right",
          legend.spacing = unit(1.5, "cm"),
          legend.key.size = unit(3, "lines"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20))

  title_plot <- ggplot() +
    theme_void() +
    annotate("text", x = 0, y = 0, label = power_effect_names[m - 1],
             size = 15) +
    theme(plot.margin = margin(b = -15, unit = "pt")) +
    labs(tag = paste0("(", LETTERS[m + 1], ")"))

  power_effect_rows[[length(power_effect_rows) + 1]] <- title_plot + big_plot +
    plot_layout(ncol = 1, heights = c(0.3, 1)) &
    theme(plot.tag = element_text(size = 40))
}

final_power_plot <- wrap_plots(power_effect_rows, ncol = 1)

pdf("all_snp_fixed_3-5_power_all_effects.pdf", width = 30, height = 16)
print(final_power_plot)
dev.off()
