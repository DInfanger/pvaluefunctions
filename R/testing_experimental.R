# #=========================================================================================
# # P-value curves for the paper
# # Author: Denis Infanger
# # Date (dd.mm.yyyy): 22.09.2018
# #=========================================================================================
# #-----------------------------------------------------------------------------------------
# # Load packages
# #-----------------------------------------------------------------------------------------
#
# library(ggpubr)
# library(metafor)
#
# #-----------------------------------------------------------------------------------------
# # Set paths
# #-----------------------------------------------------------------------------------------
#
# setwd("F:/DSBG/Papers/pvalue_function/pvaluefunctions")
#
# # gpath <- "./output/graphics"
# # tpath <- "./output/text"
# scripts_path <- "./R/"
#
# #-----------------------------------------------------------------------------------------
# # Load function
# #-----------------------------------------------------------------------------------------
#
# source(paste(scripts_path, "confidence_distributions_experimental.R", sep = "/"))
#
# #=========================================================================================
# # Graphics
# #=========================================================================================
# #-----------------------------------------------------------------------------------------
# # Stamatakis et al. (2017): br J Sports Med. Dichotomous illustration (Poole 1987)
# #-----------------------------------------------------------------------------------------
#
# plot_dat <- data.frame(
#   HR = seq(0.7, 1.9, length.out = 100)
# )
#
# plot_dat$pvalue <- NA
# plot_dat$pvalue[plot_dat$HR < 0.92 | plot_dat$HR > 1.55] <- 0
# plot_dat$pvalue[plot_dat$HR >= 0.92 & plot_dat$HR <= 1.55] <- 1
#
#
# theme_set(theme_bw())
# p_stama_dichotom <- ggplot(data = plot_dat, aes(x = HR, y = pvalue)) +
#   # geom_line(size = 1.5) +
#   geom_segment(aes(x = 0.7, xend = 0.92, y = 0, yend = 0), size = 1.5) +
#   geom_segment(aes(x = 1.55, xend = 1.9, y = 0, yend = 0), size = 1.5) +
#   geom_segment(aes(x = 0.92, xend = 1.55, y = 1, yend = 1), size = 1.5) +
#   geom_segment(aes(x = 0.92, xend = 0.92, y = 0, yend = 1), size = 1.5) +
#   geom_segment(aes(x = 1.55, xend = 1.55, y = 0, yend = 1), size = 1.5) +
#   # geom_segment(aes(x = 1.19, xend = 1.19, y = 0, yend = 1), size = 1.1, linetype = 2) +
#   geom_point(aes(x = 1.19, y = 1), size = 4.2, pch = 21, colour = "black", fill = "black", stroke = 1.7) +
#   geom_segment(aes(x = 1, xend = 1, y = 0, yend = 1), size = 1, linetype = 2) +
#   xlab("HR") +
#   ylab("Relation\nto\n95% confidence\ninterval") +
#   scale_y_continuous(breaks = c(0, 1), labels = c("Incompatible/\nsignificant", "Compatible/\nnot significant")) +
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 10), trans = "log") +
#   theme(
#     axis.title.y=element_text(colour = "black", size = 17, hjust = 0.5, vjust = 0.5, margin=margin(t = 0, r = -60, b = 0, l = 0), angle = 0),
#     axis.title.x=element_text(colour = "black", size = 17),
#     # axis.title.y=element_text(size=15,hjust=0.5, vjust=1),
#     axis.text.x=element_text(colour = "black", size=15),
#     axis.text.y=element_text(colour = "black", size=15),
#     # plot.margin=unit(c(2,2,2,2,2),"line"),
#     legend.position=c(0.9, 0.9),
#     legend.title = element_blank(),
#     legend.text=element_text(size=14),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     # panel.grid.major = element_line(colour=grey(0.8), size=0.5),
#     legend.key=element_blank(),
#     plot.title = element_text(size = 15, face = "bold"),
#     strip.text.x=element_text(size=15)
#   )
#
# p_stama_dichotom
#
#
# #-----------------------------------------------------------------------------------------
# # Stamatakis et al. (2017): br J Sports Med
# #-----------------------------------------------------------------------------------------
#
# stama <- conf_dist(
#   estimate = c(
#     0.1730998 # Stamatakis total sitting time
#     # , log(1.10) # Petersen total sitting time (both sexes)
#   )
#   , stderr = c(
#     0.1336387
#   )
#   , type = "coxreg"
#   , plot_type = "cdf"
#   , n_values = 1e4L
#   , conf_level = c(0.95)
#   , null_values = log(c(1))
#   , trans = "identity"
#   , alternative = "two_sided"
#   , log_yaxis = FALSE
#   , cut_logyaxis = 0.05
#   , xlab = "HR"
#   , xlim = log(c(0.5, 2.5))
#   , together = FALSE
#   , plot_p_limit = 1-0.999
#   , plot_counternull = TRUE
# )
#
# estimate = c(
#   0.1730998 # Stamatakis total sitting time
#   # , log(1.10) # Petersen total sitting time (both sexes)
# )
# stderr = c(
#   0.1336387
# )
# type = "coxreg"
# plot_type = "p_val"
# n_values = 1e4L
# conf_level = c(0.95)
# null_values = log(c(1))
# trans = "exp"
# alternative = "two_sided"
# log_yaxis = FALSE
# cut_logyaxis = 0.05
# xlab = "HR"
# xlim = log(c(0.7, 1.9))
# together = FALSE
# plot_p_limit = 1-0.999
# plot_counternull = TRUE
#
# #-----------------------------------------------------------------------------------------
# # Stamatakis et al. (2017) + Petersen et al. (2016)
# #-----------------------------------------------------------------------------------------
#
# stama_petersen <- conf_dist(
#   estimate = c(
#     0.1730998 # Stamatakis total sitting time
#     , log(1.10) # Petersen total sitting time (both sexes)
#     , 0.1148205 # Meta-analytic estimate
#   )
#   , stderr = c(
#     0.1336387 # Stamatakis total sitting time
#     , 0.0773228 # Petersen total sitting time (both sexes)
#     , 0.06692738 # Meta-analytic standard error
#   )
#   , est_names = c(
#     "Stamatakis et al. (2017)"
#     , "Petersen et al. (2016)"
#     , "Meta-analysis"
#   )
#   , type = "coxreg"
#   , plot_type = "cdf"
#   , n_values = 1e4L
#   , conf_level = c(0.95)
#   , null_values = log(c(1))
#   , trans = "exp"
#   , alternative = "one_sided"
#   , log_yaxis = TRUE
#   , cut_logyaxis = 0.05
#   , xlab = "HR"
#   , xlim = log(c(1.1, 1.4))
#   , together = FALSE
#   , plot_p_limit = 0.0001
#   , plot_counternull = FALSE
# )
#
# stama_petersen$plot$layers[[1]] <- NULL
#
# stama_petersen$plot <- stama_petersen$plot +
#   geom_line(aes(x = values, y = p_two, linetype = variable), size = 1.25, inherit.aes = FALSE) +
#   # scale_colour_manual(values = c("black", "black"), name = "") +
#   scale_linetype_manual(values = c(1, 2, 3), name = "") +
#   theme(
#     legend.title = element_blank()
#     , legend.spacing.x = unit(0.5, "cm")
#     , legend.key.size = unit(1.6, "line")
#     , legend.position= c(0.815, 0.87)
#
#   ) +
#   guides(
#     linetype = guide_legend(nrow = 3, keywidth = 2.5)
#   ) +
#   annotate(geom = "text", label = "B", x = 0.7, y = 1, size = 10, fontface = "bold")
#
#
# # ggsave(paste(gpath, "stamatakis_petersen.png", sep = "/"), stama_petersen$plot, width = 19*0.6, height = 11*0.6, dpi = 300, type = "cairo-png")
# # ggsave(paste(gpath, "stamatakis_petersen.tiff", sep = "/"), stama_petersen$plot, width = 19*0.6, height = 11*0.6, dpi = 1200, compression = "lzw")
#
#
# stama_petersen$conf_frame
# stama_petersen$counternull_frame
#
# # Combine those
#
# stama_peter_comb <- ggarrange(
#   stama$plot
#   , stama_petersen$plot
#   , ncol = 1
#   , nrow = 2
#   , heights = c(1, 1.1)
# )
#
# # ggsave(paste(gpath, "stamatakis_petersen_comb.png", sep = "/"), stama_peter_comb, width = 17.5*0.6, height = 19*0.6, dpi = 300, type = "cairo-png")
# # ggsave(paste(gpath, "Fig1.tiff", sep = "/"), stama_peter_comb, width = 17.5*0.6, height = 19*0.6, dpi = 400, compression = "lzw")
# # ggsave(paste(gpath, "Figure1.eps", sep = "/"), stama_peter_comb, width = 17.5*0.6, height = 19*0.6, dpi = 1200)
# # ggsave(paste(gpath, "Figure1.pdf", sep = "/"), stama_peter_comb, width = 17.5*0.6, height = 19*0.6, dpi = 1200, device=cairo_pdf)
#
# # Meta-Analysis
#
# meta_mod <- rma(
#   yi = c(
#     0.1730998 # Stamatakis total sitting time
#     , log(1.10) # Petersen total sitting time (both sexes)
#   )
#   , vi = c(
#     0.1336387 # Stamatakis total sitting time
#     , 0.0773228 # Petersen total sitting time (both sexes)
#   )^2
#   , method = "FE"
# )
#
# summary(meta_mod)
# meta_mod$se
# meta_mod$b
# forest(meta_mod, transf = exp)
# predict(meta_mod, transf = exp)
#
# #-----------------------------------------------------------------------------------------
# # Critical-value plot
# #-----------------------------------------------------------------------------------------
#
# stama_petersen$res_frame$quant_vals <- qnorm(stama_petersen$res_frame$p_one)
#
# theme_set(theme_bw())
# p <- ggplot(data = stama_petersen$res_frame, aes(x = values, y = quant_vals)) +
#   geom_line(aes(linetype = variable), size = 1.5) +
#   geom_hline(yintercept = qnorm(0.05/2), linetype = 2) +
#   xlab("HR") +
#   ylab("Critical value") +
#   scale_x_continuous(trans = "log", limits = c(0.7, 1.9), breaks = scales::pretty_breaks(n = 10)) +
#   scale_y_continuous(breaks = seq(-4, 0, 0.5)) +
#   scale_linetype_manual(name = "", values = c(1, 2, 3)) +
#   theme(
#     legend.title = element_blank()
#     , legend.spacing.x = unit(0.5, "cm")
#     , legend.key.size = unit(1.6, "line")
#     , legend.position= c(0.82, 0.87)
#     , axis.title.y.left=element_text(colour = "black", size = 17, hjust = 0.5, margin = margin(0, 10, 0, 0))
#     , axis.title.y.right=element_text(colour = "black", size = 17, hjust = 0.5, margin = margin(0, 0, 0, 10))
#     , axis.text.x=element_text(colour = "black", size=15)
#     , axis.text.y=element_text(colour = "black", size=15)
#     , axis.title.x=element_text(colour = "black", size = 17)
#     , panel.grid.minor.y = element_blank()
#     , plot.title = element_text(face = "bold")
#     , legend.text=element_text(size=15)
#   ) +
#   guides(
#     linetype = guide_legend(nrow = 3, keywidth = 2.5)
#   ) +
#   annotate(geom = "text", label = "C", x = 0.7, y = 0, size = 10, fontface = "bold") +
#   geom_label(
#     data = data.frame(theor_values = 0, p_value = qnorm(0.05/2), label = "0.05")
#     , mapping = aes(x = theor_values, y = p_value, label = label)
#     , inherit.aes = FALSE
#     , label.size = NA
#     , parse = TRUE
#     , size = 5.5
#     , hjust = "inward"
#   )
#
# # p
#
# # Combine those
#
# stama_peter_comb <- ggarrange(
#   stama$plot
#   , stama_petersen$plot
#   , p
#   , ncol = 1
#   , nrow = 3
#   # , heights = c(1, 1.1)
#   , align = c("hv")
# )
#
#
# # ggsave(paste(gpath, "Fig2.tiff", sep = "/"), stama_peter_comb, width = 17.5*0.6, height = 27*0.6, dpi = 400, compression = "lzw")
#
#
# #-----------------------------------------------------------------------------------------
# # Canning et al. (2014)
# #-----------------------------------------------------------------------------------------
#
# canning <- conf_dist(
#   estimate = c(
#     log(0.73)
#   )
#   , stderr = c(
#     0.23849
#   )
#   , type = "general_z"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , conf_level = c(0.95)
#   , null_values = c(0)
#   , trans = "exp"
#   , alternative = "two_sided"
#   , log_yaxis = FALSE
#   , cut_logyaxis = 0.05
#   , xlab = "IRR"
#   , xlim = log(c(0.35, 1.71))
#   , together = FALSE
#   , plot_p_limit = 1-0.999
#   , plot_counternull = TRUE
# )
#
# canning$plot <- canning$plot +
#   geom_vline(xintercept = 0.6, linetype = 2, size = 0.5) +
#   scale_x_continuous(trans = "log", breaks = c(seq(0.3, 1, 0.1), seq(1, 1.8, 0.2)))
# # geom_segment(aes(x = 0.5329, y = -0.5, xend = 0.5329, yend = 0.1851), linetype = 3) +
# # geom_point(aes(x = 0.5329, y = 0.1851), size = 4, pch = 21, colour = "black", fill = "white", stroke = 1.7)
#
#
# # ggsave(paste(gpath, "canning.png", sep = "/"), canning$plot, width = 19*0.6, height = 11*0.6, dpi = 300, type = "cairo-png")
# # ggsave(paste(gpath, "Fig3.tiff", sep = "/"), canning$plot, width = 19*0.6, height = 11*0.6, dpi = 400, compression = "lzw")
# # ggsave(paste(gpath, "Figure2.eps", sep = "/"), canning$plot, width = 19*0.6, height = 11*0.6, dpi = 1200)
# # ggsave(paste(gpath, "Figure2.pdf", sep = "/"), canning$plot, width = 17.5*0.6, height = 19*0.6, dpi = 1200, device=cairo_pdf)
#
#
# canning$conf_frame
# canning$counternull_frame
#
# # Assume clinically relevant reduction of risk of 40%, IRR = 0.6
# # Proportional distance
#
# log(0.73) - log(0.60)
# log(1) - log(0.73)
#
# # Proportional distance to confidence intervals
#
# log(1.17) - log(1) # Distance upper limits to 1
# log(0.60) - log(0.45) # Distance lower limit to 0.6
#
# #-----------------------------------------------------------------------------
# # Höchsmann et al. (2017)
# #-----------------------------------------------------------------------------
#
# # (1.36 - -3.5)/abs(c(1.36/qt(0.327/2, df = 29))) > qt(0.05, 29, lower.tail = FALSE) # Chow et al. (2008), p. 52
# # 2*pt(-abs(1.36 - -3.5)/abs(c(1.36/qt(0.327/2, df = 29))), df = 29)
#
# hochsmann <- conf_dist(
#   estimate = c(1.36)
#   # , n = c(30)
#   , df = c(29)
#   , stderr = abs(c(1.36/qt(0.327/2, df = 29)))
#   # , tstat = c(1.5751)
#   , type = "linreg"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , est_names = c("VO2peak")
#   , conf_level = c(0.95)
#   , null_values = c(-3.5)
#   , log_yaxis = FALSE
#   , cut_logyaxis = 0.05
#   , trans = "identity"
#   , alternative = "one_sided"
#   , xlab = expression("Adjusted difference in VO"[2*"peak"]*" [ml/(kg" %*% " min)]")
#   , xlim = c(-4, 6.5)
#   , together = FALSE
#   , plot_p_limit = 1 - 0.9999
# )
#
# hochsmann$plot$layers[[1]] <- NULL
#
# hochsmann$plot <- hochsmann$plot +
#   geom_line(aes(x = values, y = p_two, linetype = hypothesis), size = 1.25, inherit.aes = FALSE) +
#   # scale_colour_manual(values = c("black", "black"), name = "") +
#   scale_linetype_manual(values = c(1, 2), name = "", labels = c(
#     expression("H"[1]*": "*beta*" < "*theta)
#     , expression("H"[1]*": "*beta*" > "*theta))
#   ) +
#   # scale_linetype_manual(values = c(1, 2), name = "", labels = c(
#   #   expression("H"[1]*": "*beta*" < "*delta)
#   #   , expression("H"[1]*": "*beta*" > -"*delta))
#   # ) +
#   scale_x_continuous(limits = c(-4, 6.5), breaks = seq(-100, 100, 1)) +
#   theme(
#     legend.title = element_blank()
#     , legend.text.align = 0
#     , legend.text=element_text(size=17)
#     , legend.spacing.x = unit(0.5, "cm")
#     , legend.key.size = unit(2.1, "line")
#     , legend.position= c(0.9, 0.9)
#     , axis.title.x=element_text(colour = "black", size = 17, margin=margin(11, 0, 0, 0))
#   )
#
#
# # ggsave(paste(gpath, "hochsmann.png", sep = "/"), hochsmann$plot, width = 19*0.6, height = 11*0.6, dpi = 300, type = "cairo-png")
# # ggsave(paste(gpath, "Fig4.tiff", sep = "/"), hochsmann$plot, width = 19*0.6, height = 11*0.6, dpi = 400, compression = "lzw")
# # ggsave(paste(gpath, "Figure3.eps", sep = "/"), hochsmann$plot, width = 19*0.6, height = 11*0.6, dpi = 1200)
# # ggsave(paste(gpath, "Figure3.pdf", sep = "/"), hochsmann$plot, width = 19*0.6, height = 11*0.6, dpi = 1200, device=cairo_pdf)
#
# hochsmann$conf_frame
# hochsmann$counternull_frame
#
# # s_grid <- expand.grid(
# #   beta = seq(1.35, 1.45, length.out = 1e3)
# #   , df = seq(25, 17+15 - 3, length.out = 1e3)
# # )
# #
# # ci_up <- s_grid$beta + qt(0.975, df = s_grid$df)*abs((s_grid$beta/qt(0.327/2, df = s_grid$df)))
# # ci_low <- s_grid$beta - qt(0.975, df = s_grid$df)*abs((s_grid$beta/qt(0.327/2, df = s_grid$df)))
# #
# # up_ind <- which(round(ci_up, 1) == 4.1)
# # low_ind <- which(round(ci_low, 1) == -1.4)
# #
# # intersect(up_ind, low_ind)
# #
# # tail(s_grid[intersect(up_ind, low_ind), ])
#
#
# pool_sd <- function(m = NULL, s = NULL, n = NULL) {
#
#   pooled_mean <- weighted.mean(m, n)
#
#   pooled_sd <- sqrt(sum((n - 1)*s^2 + n*(m - pooled_mean)^2)/(sum(n) - 1))
#
#   list(mean = pooled_mean, sd = pooled_sd)
#
# }
#
# pool_sd(m = c(35.7, 36.4), s = c(5.8, 7.3), n = c(16, 13))$sd/2
#
#
# #-----------------------------------------------------------------------------
# # Examples
# #-----------------------------------------------------------------------------
#
# # 1 Proportion
#
# res <- conf_dist(
#   estimate = c(0.5)
#   , n = c(50)
#   # , df = 29
#   # , stderr = se_d
#   # , tstat = 31.504
#   , type = "prop"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , est_names = c("p1")
#   , log_yaxis = TRUE
#   , cut_logyaxis = 0.05
#   , conf_level = c(0.99, 0.90, 0.80, 0.5, 0.3)
#   , null_values = c(0.5)
#   , trans = "identity"
#   , alternative = "two_sided"
#   # , xlab = "Proportion"
#   # , xlim = c(0.2, 0.7)
#   , together = FALSE
#   , plot_p_limit = 1 - 0.999
#   , plot_counternull = TRUE
# )
#
# res$conf_frame
# res$counternull_frame
# res$point_est
#
# # binom::binom.confint(0.25*30, n =30, 0.90, "wilson")
# #
# # p <- seq(0, 1, 0.001)
# # p <- p[!p %in% c(0, 1)]
# #
# # pval <- sapply(p, FUN = function(x, n, p){binom.test(x = x, n = n, p = p)$p.value}, x = 25-3, n = 50)
# #
# # plot(pval~p, type = "l")
# # abline(h = c(0.05))
# # abline(v = 0.5)
# #
# # which.min(abs(pval[p<0.44] - 0.4799))
# # p[p<0.44][391]
# #
# # binom.test(25-3, 50, p = 0.5)
# # binom.test(25-3, 50, p = 0.391)
#
# # General z (Bender et al. 2005)
#
# n11 <- 30
# n12 <- 63
# n22 <- 63
# ntot <- 156
#
# p11 <- n11/ntot
# p12 <- n12/ntot
# p1 <- (n11 + n12)/ntot
#
# d <- p1 - p11/(p1)
#
# se_d <- sqrt((1/ntot)*(((p11*p12)/p1^3) + (p1*(1 - p1))))
#
# res <- conf_dist(
#   estimate = d
#   # , n = 30
#   # , df = 29
#   , stderr = se_d
#   # , tstat = 31.504
#   , type = "general_z"
#   , plot_type = "p_val"
#   , cut_logyaxis = 0.05
#   , log_yaxis = TRUE
#   , n_values = 1e4L
#   , est_names = "Diff in proportions"
#   , conf_level = c(0.95, 0.90, 0.80)
#   , null_values = c(0)
#   , trans = "identity"
#   , alternative = "two_sided"
#   # , xlab
#   , xlim = c(-1,  1)
#   , together = FALSE
#   , plot_p_limit = 1-0.999
# )
#
# res$conf_frame
#
# # Coxreg (Shakespeare et al. 2001)
#
# res <- conf_dist(
#   estimate = log(0.72)
#   # , n = 30
#   # , df = 29
#   , stderr = 0.187615
#   # , tstat = 31.504
#   , type = "coxreg"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , est_names = "HR"
#   , conf_level = c(0.95, 0.90, 0.80)
#   , null_values = c(0)
#   , trans = "exp"
#   , alternative = "one_sided"
#   # , xlab
#   # , xlim = c(90, 110)
#   , together = FALSE
# )
#
# res$conf_frame
#
# # T-tests: one-sample
#
# set.seed(142857)
# n <- 30
# x <- rnorm(n, 100, 15)
#
# t.test(x, alternative = "two.sided")
# t.test(x, alternative = "less")
# t.test(x, alternative = "greater")
#
# conf_dist(
#   estimate = 100.0357
#   , n = 30
#   , df = 29
#   # , stderr
#   , tstat = 31.504
#   , type = "ttest"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   # , est_names = "mean"
#   # , conf_level = c(0.95, 0.90, 0.80)
#   # , null_values = c(95, 100)
#   , trans = "identity"
#   , alternative = "one_sided"
#   # , xlab
#   # , xlim = c(90, 110)
#   , together = FALSE
# )
#
# set.seed(142857)
# n <- 30
# x1 <- rnorm(n, 100, 15)
# x2 <- rnorm(n/2, 90, 10)
#
# t.test(x1, alternative = "two.sided")
# t.test(x2, alternative = "two.sided")
#
# res <- conf_dist(
#   estimate = c(100.0357, 88.5764, 110)
#   # , n = c(30, 15, 27)
#   , df = c(29, 14, 26)
#   # , stderr
#   , tstat = c(31.504, 45.351, 30)
#   , type = "ttest"
#   , plot_type = "s_val"
#   , log_yaxis = TRUE
#   , cut_logyaxis = 0.05
#   , n_values = 1e4L
#   , est_names = c("mean1", "mean2", "mean3")
#   , conf_level = c(0.95, 0.90, 0.80)
#   , null_values = c(95, 100)
#   , trans = "identity"
#   , alternative = "two_sided"
#   , xlab = "Means"
#   # , xlim = c(90, 110)
#   , together = TRUE
#   , plot_p_limit = 1 - 0.999
#   , plot_counternull = TRUE
# )
#
# res$point_est
#
# # T-tests: two-samples unpaired
#
# set.seed(142857)
# n <- 30
# x1 <- rnorm(n, 100, 15)
# x2 <- rnorm(round(n/3), 95, 10)
#
# t.test(x1, x2, paired = FALSE, alternative = "two.sided", conf.level = 0.95)
#
# res <- conf_dist(
#   estimate = c(8.45763)
#   # , n = c(3000)
#   , df = c(33.5)
#   # , stderr
#   , tstat = c(2.0757)
#   , type = "ttest"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , est_names = c("mean")
#   , conf_level = c(0.95, 0.90, 0.80)
#   , null_values = c(0)
#   , trans = "identity"
#   , alternative = "two_sided"
#   , xlab = "Difference in means"
#   # , xlim = c(90, 110)
#   , together = FALSE
# )
#
# res$conf_frame
# res$counternull_frame
#
# # T-tests: two-samples paired
#
# set.seed(142857)
# n <- 30
# x1 <- rnorm(n, 100, 15)
# x2 <- rnorm(n, 95, 16)
#
# t.test(x1, x2, paired = TRUE, alternative = "two.sided")
# t.test(x1, x2, paired = TRUE, alternative = "less")
#
# res <- conf_dist(
#   estimate = c(6.224147)
#   , n = c(30)
#   , df = c(29)
#   # , stderr
#   , tstat = c(1.5751)
#   , type = "ttest"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , est_names = c("mean")
#   , log_yaxis = TRUE
#   , cut_logyaxis = 0.05
#   , conf_level = c(0.95, 0.90, 0.80)
#   , null_values = c(0)
#   , trans = "identity"
#   , alternative = "one_sided"
#   , xlab = "Mean difference"
#   # , xlim = c(90, 110)
#   , together = FALSE
# )
#
# which(sort(res$res_frame$p_two[res$res_frame$values < 5], decreasing = TRUE) <= 0.001)
#
# sort(res$res_frame$value[res$res_frame$values < 5], decreasing = TRUE)[2477]
#
# which.max(res$res_frame$p_two[res$res_frame$values > 5] <= 0.001)
# res$res_frame$values[res$res_frame$values > 5][2935]
#
# # Linear regression
#
# n <- 30
# beta0 <- 1
# beta1 <- -2.5
# beta2 <- 1.75
# beta3 <- 0.75
# sigma <- 30
#
# set.seed(142857)
# x1 <- runif(n, 0, 100)
# # x2 <- runif(n, 0, 100)
# x2 <- rbinom(n, 1, p = 0.5)
#
# y <- beta0 + beta1*x1 + beta2*x2 + rnorm(n, 0, sigma)
# # y <- beta0 + beta1*x1 + beta2*x2 + beta3*x2^2 + rnorm(n, 0, sigma)
#
# dat <- data.frame(
#   y = y
#   , x1 = x1
#   # , x2 = x2
#   , x2 = factor(x2)
# )
#
# mod <- lm(y~x1 + x2, data = dat)
#
# conf_dist(
#   estimate = c(coef(mod)[2:3])
#   # , n = c(30)
#   , df = c(mod$df.residual, mod$df.residual)
#   , stderr = c(0.1872, 11.2759)
#   # , tstat = c(1.5751)
#   , type = "linreg"
#   , plot_type = "pdf"
#   , n_values = 1e4L
#   , est_names = c("x1", "x2")
#   , conf_level = c(0.95, 0.90, 0.80)
#   , null_values = c(0)
#   , trans = "identity"
#   , alternative = "two_sided"
#   , xlab = "Coefficients"
#   # , xlim = c(90, 110)
#   , together = FALSE
# )
#
#
# res <- conf_dist(
#   estimate = c(1.4)
#   # , n = c(30)
#   , df = c(16 + 13 - 3)
#   , stderr = abs(c(1.4/qt(0.327/2, df = 16 + 13 - 3)))
#   # , tstat = c(1.5751)
#   , type = "linreg"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , est_names = c("BMI")
#   , conf_level = c(0.95, 0.90, 0.80)
#   , null_values = c(0)
#   , trans = "identity"
#   , alternative = "one_sided"
#   , xlab = "Coefficients"
#   # , xlim = c(90, 110)
#   , together = FALSE
# )
#
# res$plot + geom_vline(xintercept = 3, linetype = 2, colour = "red")
#
# res$conf_frame
# res$counternull_frame
#
# # Linear Regression: Gioia
#
# res <- conf_dist(
#   estimate = c(-0.523)
#   # , n = c(30)
#   , df = c(193.53007)
#   , stderr = (0.27814)
#   # , tstat = c(1.5751)
#   , type = "linreg"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , est_names = c("VO2peak")
#   , conf_level = c(0.95, 0.90, 0.80)
#   , null_values = c(0)
#   , trans = "identity"
#   , alternative = "two_sided"
#   , xlab = "Coefficients"
#   # , xlim = c(90, 110)
#   , together = FALSE
# )
#
# res$conf_frame
# res$counternull_frame
#
# # Stamatakis et al. (2017): Whitehall II
# # Biswas et al. (2015)
# # Ford et al. (2009)
# # Petersen et al. (2016)
#
# res <- conf_dist(
#   estimate = c(
#     0.1730998 # Stamatakis total sitting time
#     # , log(1.10) # Petersen total sitting time (both sexes)
#     # , log(1.91)
#     # log(1.14) # Ford
#     # , log(1.31) # Stamatakis highest tv time
#   )
#   # , n = c()
#   # , df = c()
#   , stderr = c(
#     0.1336387
#     # , 0.0773228
#     # , 0.0777587
#     # 0.176129
#     # , 0.150659
#   )
#   # , tstat = c(1.5751)
#   , type = "coxreg"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   # , est_names = c(
#   #   "Stamatakis et al."
#   #   , "Petersen et al."
#   #   # , "Biswas"
#   #    # "Ford"
#   #   # , "Stamatakis TV time"
#   # )
#   , conf_level = c(0.95, 0.9)
#   , null_values = c(0)
#   , trans = "exp"
#   , alternative = "two_sided"
#   , xlab = "HR"
#   , xlim = log(c(0.7, 1.8))
#   , together = TRUE
# )
#
# 2*pnorm(-0.1730998/0.1336387)
#
# res$conf_frame
# res$counternull_frame
#
# # sgrid <- expand.grid(
# #   lhr = seq(log(1.185), log(1.1945), length.out = 1e3)
# #   , se = seq(0.13, 0.14, length.out = 5e3)
# # )
# #
# # cifun_up <- function(var1, var2){
# #   exp(var1 + 1.959964*var2)
# # }
# #
# # cifun_lwr <- function(var1, var2){
# #   exp(var1 - 1.959964*var2)
# # }
# #
# # pfun <- function(var1, var2){
# #   2*pnorm(-var1/var2)
# # }
# #
# # hrs <- exp(sgrid$lhr)
# # ci_ups <- mapply(cifun_up, var1 = sgrid$lhr, var2=sgrid$se)
# # ci_lwrs <- mapply(cifun_lwr, var1 = sgrid$lhr, var2=sgrid$se)
# # pees <- mapply(pfun, var1 = sgrid$lhr, var2=sgrid$se)
# #
# # # hrs_ind <- which(round(hrs, 2) == 1.19)
# # ci_upr_ind <- which(round(ci_ups, 2) == 1.55)
# # ci_lwr_ind <- which(round(ci_lwrs, 2) == 0.92)
# # pees_ind <- which(round(pees, 2) == 0.22)
# #
# # intersect(intersect(ci_lwr_ind, ci_upr_ind), pees_ind)
# #
# # summary(sgrid[intersect(ci_lwr_ind, ci_upr_ind),])
# #
# # which.max(pees[intersect(ci_lwr_ind, ci_upr_ind)])
# #
# # sgrid[intersect(ci_lwr_ind, ci_upr_ind), ][342459, ]
#
# # Canning et al. (2014): IRR
#
# res <- conf_dist(
#   estimate = c(log(0.73))
#   # , n = c(20)
#   # , df = c(229)
#   , stderr = c(0.23849)
#   # , tstat = c(-4.687)
#   , type = "general_z"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , est_names = c("IRR")
#   , conf_level = c(0.95)
#   , null_values = c(0)
#   , trans = "exp"
#   , alternative = "two_sided"
#   , xlab = "IRR"
#   , xlim = c(0.1, 1.5)
#   # , together = FALSE
#   , plot_p_limit = 1e-16
# )
#
# res$conf_frame
# res$counternull_frame
#
# # Logistic regression
#
# n <- 100
#
# set.seed(142857)
#
# x <- runif(n, 18, 60)
# sex <- rbinom(n, 1, 0.5)
# # z <- -2 + 0.035*x
# # z <- -20 + (1/2)*x
# z <- -(65/7) + (5/21)*x + 0.1*sex
# pr <- 1/(1 + exp(-z))
# y <- rbinom(n, 1, pr)
#
# dat <- data.frame(
#   y = y
#   , x = x
#   , sex = factor(sex, levels = c(0, 1), labels = c("male", "female"))
# )
#
# glm_mod <- glm(y~x + sex, family = binomial, data = dat)
#
# res <- conf_dist(
#   estimate = c(coef(glm_mod)[2:3])
#   # , n = c(30)
#   # , df = c(mod$df.residual, mod$df.residual)
#   , stderr = c(0.04728, 0.64783)
#   # , tstat = c(1.5751)
#   , type = "logreg"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , est_names = c("x", "sex")
#   , conf_level = c(0.95, 0.90, 0.80)
#   , null_values = c(log(0.75))
#   , trans = "exp"
#   , xlab = "Odds ratio"
#   # , xlim = c(90, 110)
#   , together = TRUE
#   , plot_p_limit = 0.001
# )
#
# res$conf_frame
# res$counternull_frame
#
# res$plot + geom_vline(xintercept = 0.233114)
#
# # Correlations
#
# cor.test(swiss$Examination, swiss$Catholic, method = "pearson", conf.level = 0.95)
# set.seed(142857)
# cor.test(rnorm(20, 100, 15), rnorm(20, 100, 15), method = "pearson", conf.level = 0.95)
#
#
# res <- conf_dist(
#   estimate = c(0.9, 0.65)
#   , n = c(123, 9)
#   # , df = c(45)
#   # , stderr = c(0.04728, 0.64783)
#   # , tstat = c(-4.687, -5)
#   , type = "spearman"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , log_yaxis = TRUE
#   , cut_logyaxis = 0.05
#   , est_names = c("corr1", "corr2")
#   , conf_level = c(0.95, 0.90, 0.80, 0.5)
#   , null_values = c(0.5, 0.85)
#   , trans = "identity"
#   , alternative = "one_sided"
#   , xlab = "Correlation"
#   , xlim = c(0, 1)
#   , together = FALSE
#   , plot_p_limit = 1 - 0.999
#   , plot_counternull = TRUE
# )
#
# res$conf_frame
# res$counternull_frame
# res$plot + geom_vline(xintercept = 0.6842105)
# res$point_est
#
# # Greenland (2012)
#
# res <- conf_dist(
#   estimate = c(log(1.5))
#   # , n = c(47)
#   # , df = c(45)
#   , stderr = c(sqrt(601/30)/20)
#   # , tstat = c(-4.687)
#   , type = "general_z"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , est_names = c("RR")
#   , conf_level = c(0.95, 0.90, 0.80)
#   , null_values = c(0, log(2))
#   , trans = "exp"
#   , alternative = "two_sided"
#   , xlab = "RR"
#   , xlim = log(c(0.6, 3.2))
#   , together = FALSE
#   , plot_p_limit = 1-0.999
#   , plot_counternull = TRUE
# )
#
# res$conf_frame
# res$counternull_frame
#
# # Standard deviation
#
# set.seed(142857)
# est_1 <- var(x <- rnorm(200, 100, sqrt(250)))
# sd(x)
#
# res <- conf_dist(
#   estimate = c(est_1, 300)
#   , n = c(length(x), 250)
#   # , df = c(45)
#   # , stderr = c(sqrt(601/30)/20)
#   # , tstat = c(-4.687)
#   , type = "var"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , est_names = c("Var1", "Var2")
#   , log_yaxis = TRUE
#   , cut_logyaxis = 0.05
#   , conf_level = c(0.95)
#   , null_values = c(15^2, 18^2)
#   , trans = "identity"
#   , alternative = "two_sided"
#   , xlab = "Var"
#   # , xlim = c(200, 300)
#   , together = TRUE
#   , plot_p_limit = 1 - 0.9999
#   , plot_counternull = TRUE
# )
#
# res$conf_frame
# res$counternull_frame
# res$point_est
#
# res$plot + geom_vline(xintercept = c(250), linetype = 2)
#
# res$plot + geom_vline(xintercept = as.numeric(res$point_est[1, 1:3]), linetype = 2)
#
# (291.5835*(19))/(2*Rgamma.inv(19/2, 1/2))
#
# # Concurve meta-analysis replication
#
# res <- conf_dist(
#   estimate = c(-0.13)
#   # , n = c(20, 100)
#   # , df = c(45)
#   , stderr = c(0.224494)
#   # , tstat = c(-4.687)
#   , type = "general_z"
#   , plot_type = "s_val"
#   , n_values = 1e4L
#   , est_names = c("")
#   , log_yaxis = FALSE
#   , cut_logyaxis = 0.05
#   , conf_level = c(0.95)
#   # , null_values = c(15^2, 18^2)
#   , trans = "identity"
#   , alternative = "two_sided"
#   , xlab = "Var"
#   , xlim = c(-1, 1)
#   , together = TRUE
#   , plot_p_limit = 1 - 0.9999
# )
#
# # Meta-analysis by Köchli et al. 2018, CRAE
#
# res <- conf_dist(
#   estimate = c(-0.37)
#   , stderr = c(0.0663277)
#   # , tstat = c(-4.687)
#   , type = "general_z"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , est_names = c("Coefficient CRAE vs. BMI")
#   , log_yaxis = TRUE
#   , cut_logyaxis = 0.05
#   , conf_level = c(0.95)
#   , null_values = c(0)
#   , trans = "identity"
#   , alternative = "two_sided"
#   , xlab = "Var"
#   , xlim = c(-0.75, -0.05)
#   , together = TRUE
#   , plot_p_limit = 1 - 0.99999
# )
#
# res$conf_frame
# res$counternull_frame
#
# # Difference in proporions
#
# res <- conf_dist(
#   estimate = c(49/100, 40/90)
#   , n = c(100, 90)
#   , type = "propdiff"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , est_names = c("Test")
#   , log_yaxis = TRUE
#   , cut_logyaxis = 0.05
#   , conf_level = c(0.95, 0.99)
#   , null_values = c(0, 0.3)
#   , trans = "identity"
#   , alternative = "two_sided"
#   # , xlab = "Difference of proportions"
#   # , xlim = c(-0.75, 0.5)
#   , together = FALSE
#   , plot_p_limit = 1 - 0.9999
#   , plot_counternull = TRUE
# )
#
# res$counternull_frame
# res$conf_frame
# BinomDiffCI(49, 100, 40, 90, conf.level = c(0.95, 0.99), method = "scorecc")
# BinomDiffCI(49, 100, 40, 90, conf.level = c(0.95, 0.99), method = "scorecc", sides = "left")
# BinomDiffCI(49, 100, 40, 90, conf.level = c(0.95, 0.99), method = "scorecc", sides = "right")
#
# # Difference in proporions (Agresti-Caffo correction)
#
# # First proportion
#
# x1 <- 8
# n1 <- 40
#
# # Second proportion
#
# x2 <- 11
# n2 <- 30
#
# # Apply the correction
#
# p1hat <- (x1 + 1)/(n1 + 2)
# p2hat <- (x2 + 1)/(n2 + 2)
#
# # The estimator
#
# est0 <- (x1/n1) - (x2/n2)
#
# # The estimator and its standard error using the correction
#
# est <- p1hat - p2hat
# se <- sqrt(((p1hat*(1 - p1hat))/(n1 + 2)) + ((p2hat*(1 - p2hat))/(n2 + 2)))
#
# res <- conf_dist(
#   estimate = c(est)
#   , stderr = c(se)
#   , type = "general_z"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   # , est_names = c("Estimate")
#   , log_yaxis = TRUE
#   , cut_logyaxis = 0.05
#   , conf_level = c(0.95, 0.99)
#   , null_values = c(0, 0.3)
#   , trans = "identity"
#   , alternative = "two_sided"
#   , xlab = "Difference of proportions"
#   # , xlim = c(-0.75, 0.5)
#   , together = FALSE
#   , plot_p_limit = 1 - 0.9999
# )
#
# res$conf_frame
#
#
# estimate = c(0.804037549)
# stderr = c(0.331819298)
# type = "logreg"
# plot_type = "p_val"
# n_values = 1e4L
# est_names = c("GPA")
# conf_level = c(0.95, 0.90, 0.80)
# null_values = c(log(1)) # null value on the log-odds scale
# trans = "exp"
# alternative = "two_sided"
# log_yaxis = TRUE
# cut_logyaxis = 0.05
# xlab = "Odds Ratio (GPA)"
# xlim = log(c(0.4, 5.2)) # axis limits on the log-odds scale
# together = FALSE
# plot_p_limit = 1 - 0.999
# plot_counternull = TRUE
#
#
# res <- conf_dist(
#   estimate = c(0.804037549)
#   , stderr = c(0.331819298)
#   , type = "logreg"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , est_names = c("GPA")
#   , conf_level = c(0.95, 0.90, 0.80)
#   # , null_values = c(log(1), log(6)) # null value on the log-odds scale
#   , trans = "exp"
#   , alternative = "two_sided"
#   , log_yaxis = TRUE
#   , cut_logyaxis = 0.05
#   , xlab = "Odds Ratio (GPA)"
#   , xlim = log(c(1.5, 5.2)) # axis limits on the log-odds scale
#   , together = FALSE
#   , plot_p_limit = 1 - 0.999
#   , plot_counternull = TRUE
# )
#
# res <- conf_dist(
#   estimate = c(0.804037549)
#   , stderr = c(0.331819298)
#   , type = "logreg"
#   , plot_type = "p_val"
#   , n_values = 1e4L
#   , est_names = c("GPA")
#   # , conf_level = c(0.95, 0.90, 0.80)
#   # , null_values = c(0, 6) # null value on the log-odds scale
#   # , trans = "exp"
#   , alternative = "two_sided"
#   , log_yaxis = TRUE
#   , cut_logyaxis = 0.05
#   , xlab = "Log Odds Ratio (GPA)"
#   , xlim = (c(-5, 5.2)) # axis limits on the log-odds scale
#   , together = FALSE
#   , plot_p_limit = 1 - 0.999
#   , plot_counternull = TRUE
# )
