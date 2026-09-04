## =====================================================================
## 08_figures.R —— 论文全部图（v2 数据，output_v2）
##   Figure 1  PRISMA 式流程图
##   Figure 2  年度报告数 + 途径构成
##   Figure 4  核心 PT 森林图（总体 + 三途径分层）
##   Figure 5  TTO 累积分布（总体/分途径）+ 参数
##   Figure 5B FAERS vs Canada ROR 对比
##   Figure 7  头对头配对森林图
## 运行：LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 R_LIBS_USER=~/Rlibs Rscript 08_figures.R
## =====================================================================
suppressMessages({library(data.table); library(ggplot2)})
OUT <- "output_v2"
FIG <- file.path(OUT, "fig"); dir.create(FIG, showWarnings = FALSE)
th <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill="grey92"))

## ---------- Figure 1 流程图 ----------
dedup <- fread(file.path(OUT, "00_全局去重日志.csv"))
getv <- function(k) dedup[metric==k]$value
N_RAW <- as.numeric(getv("raw_demo_rows")); N_DEDUP <- as.numeric(getv("post_dedup_reports"))
N_REMOVED <- as.numeric(getv("removed_rows"))
logq <- fread(file.path(OUT, "00_季度处理日志.csv"))
N_SCOP <- NA  # 从 faers_target.rds 取
tgt <- readRDS(file.path(OUT, "faers_target.rds"))
N_SCOP <- uniqueN(tgt[[3]][drug_group=="scopolamine"]$primaryid)

box <- function(x, y, w, h, label, col="#eef3f9", cex=0.85) {
  rect(x - w/2, y - h/2, x + w/2, y + h/2, col = col, border = "#2b5d8a", lwd = 1.4)
  text(x, y, label, cex = cex, xpd = NA)
}
arr <- function(x0, y0, y1) arrows(x0, y0, x0, y1, length = 0.09, lwd = 1.4)
png(file.path(FIG, "fig1_flow.png"), width = 2000, height = 1250, res = 220)
par(mar = c(1,1,2,1)); plot(NA, xlim=c(0,10), ylim=c(0,10), axes=FALSE, xlab="", ylab="", main="Figure 1. Study flow (FAERS 2012Q4–2026Q2)", cex.main=1.05)
box(3, 9.0, 5.6, 1.15, sprintf("FAERS quarterly ASCII files, 55 quarters (2012Q4–2026Q2)\n%d report records", N_RAW))
arr(3, 8.4, 7.6)
box(3, 7.0, 5.6, 1.15, sprintf("Cross-quarter deduplication (one report per caseid:\nhighest caseversion, latest FDA receipt date)\n%d unique reports (%d duplicates removed)", N_DEDUP, N_REMOVED))
arr(3, 6.4, 5.6)
box(3, 5.0, 5.6, 1.0, sprintf("Reports listing scopolamine (any drug role)\nn = %s", format(N_SCOP, big.mark=",")))
arr(3, 4.5, 3.6)
box(3, 3.0, 5.6, 1.0, "Main analysis set (exclusion of non-event PTs:\nOFF LABEL USE, DRUG INEFFECTIVE, etc.)")
## 右侧分支
for (yy in c(7.0, 5.0, 3.0)) segments(5.8, yy, 6.8, yy, lwd=1.2, col="#666666")
box(8.15, 7.0, 3.4, 1.05, "Disproportionality analysis\n(ROR primary; PRR/IC parallel)", col="#f3fbf1", cex=0.8)
box(8.15, 5.0, 3.4, 1.05, "TTO Weibull + survreg re-fit;\nroute-stratified analyses", col="#f3fbf1", cex=0.8)
box(8.15, 3.0, 3.4, 1.05, "External validation (Canada Vigilance);\nlabel comparison; head-to-head", col="#f3fbf1", cex=0.8)
dev.off(); cat("[fig] fig1_flow.png\n")

## ---------- Figure 2 年度 + 途径 ----------
logq[, year := as.integer(substr(quarter, 1, 4))]
ann <- logq[, .(n = sum(n_reports)), by=year]
ann_all <- ann; ann_scop <- logq[, .(n = sum(n_target_reports)), by=year]
p2a <- ggplot(ann, aes(year, n)) +
  geom_col(fill="#4C82B6", alpha=.9) +
  labs(x="Calendar year", y="Library reports (deduplicated)",
       title="Figure 2A. FAERS library size by year") + th
p2b <- ggplot(ann_scop, aes(year, n)) +
  geom_col(fill="#C0504D", alpha=.9) +
  labs(x="Calendar year", y="Scopolamine reports",
       title="Figure 2B. Scopolamine reports by year") + th
rt <- fread(file.path(OUT, "01b_给药途径分布.csv"))[drug_group=="scopolamine"]
rt[, route_std := factor(route_std, levels=c("Transdermal","Oral","Parenteral","Other","Unknown/Missing"))]
p2c <- ggplot(rt, aes(x="", y=n_reports, fill=route_std)) +
  geom_col(width=1, color="white") + coord_polar("y") +
  scale_fill_brewer(palette="Set2", name="Route") +
  labs(x=NULL, y=NULL, title="Figure 2C. Scopolamine route composition") + th +
  theme(axis.text=element_blank(), axis.ticks=element_blank())
ggsave(file.path(FIG,"fig2a_annual_library.png"), p2a, width=6.5, height=4.2, dpi=320)
ggsave(file.path(FIG,"fig2b_annual_scop.png"),    p2b, width=6.5, height=4.2, dpi=320)
ggsave(file.path(FIG,"fig2c_route_pie.png"),      p2c, width=5.2, height=4.6, dpi=320)
cat("[fig] fig2a/2b/2c\n")

## ---------- Figure 4 核心 PT 森林图（总体+途径） ----------
rs <- fread(file.path(OUT, "25_途径分层_核心PT.csv"))
ov <- fread(file.path(OUT, "22b_核心PT_总体.csv"))
d4 <- rbind(ov[, .(pt, grp="Overall", a, ROR, lcl, ucl)],
            rs[, .(pt, grp=route, a, ROR, lcl, ucl)])
top_pt <- ov[order(-a)][1:20]$pt
d4 <- d4[pt %in% top_pt]
d4[, grp := factor(grp, levels=c("Overall","Transdermal","Oral","Parenteral"))]
d4[, pt := factor(pt, levels = rev(top_pt))]
p4 <- ggplot(d4[!is.na(ROR)], aes(ROR, pt, color=grp)) +
  geom_point(position=position_dodge(width=.6), size=1.7) +
  geom_errorbarh(aes(xmin=lcl, xmax=ucl), height=0, position=position_dodge(width=.6), linewidth=.45) +
  geom_vline(xintercept=1, linetype="dashed", color="grey40") +
  scale_x_log10() + scale_color_brewer(palette="Set1", name="") +
  labs(x="ROR (95% CI, log scale)", y=NULL,
       title="Figure 4. Disproportionality signals for core PTs, overall and by route",
       caption="ROR = reporting odds ratio vs all other reports in the deduplicated library; a≥3 required. Parenteral n small — interpret with caution.") + th
ggsave(file.path(FIG,"fig4_forest.png"), p4, width=9.5, height=8.5, dpi=320)
cat("[fig] fig4_forest.png\n")

## ---------- Figure 5 TTO ----------
tto <- readRDS(file.path(OUT, "23b_tto_报告级.rds"))
ci  <- fread(file.path(OUT, "23_TTO_重拟_CI.csv"))
## 主方法：与 survreg 差异>30% 时取 survreg，否则 fitdistr（与 Methods 一致）
wk <- ci[label=="overall"]
pick <- if (abs(wk[method=="survreg"]$shape - wk[method=="fitdistr"]$shape) /
            wk[method=="survreg"]$shape > 0.3) "survreg" else "fitdistr"
sh <- wk[method==pick]
ecdf_all <- tto[, .(tto)]
ecdf_all[, grp := "Overall"]
ecdf_rt <- tto[route_std %in% c("Transdermal","Oral","Parenteral"), .(tto, grp=route_std)]
e5 <- rbind(ecdf_all, ecdf_rt)
p5 <- ggplot(e5, aes(tto, colour=grp)) +
  stat_ecdf(linewidth=.8, pad_point = FALSE) +
  coord_cartesian(xlim=c(0, 365)) +
  scale_colour_brewer(palette="Set1", name="") +
  annotate("text", x=200, y=.35, hjust=0, size=3.4,
           label=sprintf("Weibull shape = %.3f (95%% CI %.3f–%.3f, %s)\nmedian TTO = %.0f days\nEarly-failure pattern (shape CI < 1)",
                         sh$shape, sh$shape_lcl, sh$shape_ucl, pick, sh$median_tto)) +
  labs(x="Time-to-onset (days, capped display at 365)", y="Cumulative proportion",
       title="Figure 5. Time-to-onset distributions") + th
ggsave(file.path(FIG,"fig5_tto.png"), p5, width=8.5, height=5.5, dpi=320)
cat("[fig] fig5_tto.png\n")

## ---------- Figure 5B 外部验证 ----------
ev <- fread(file.path(OUT, "21_两库对比_剔除非事件.csv"))
evs <- ev[faers_signal==TRUE & !is.na(canada_ror)]
agree <- evs[verdict %in% c("replicated","directional"), .N]/evs[, .N]
p5b <- ggplot(evs, aes(faers_ror, canada_ror, color=verdict)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
  geom_point(size=2, alpha=.85) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values=c(replicated="#1A7F37", directional="#E6A23C",
                              discordant="#B42C3F", `not evaluable`="grey60")) +
  annotate("text", x=min(evs$faers_ror), y=max(evs$canada_ror), hjust=0, size=3.4,
           label=sprintf("Directional agreement: %d/%d (%.0f%%)", evs[verdict %in% c("replicated","directional"), .N], evs[, .N], 100*agree)) +
  labs(x="FAERS ROR (log)", y="Canada Vigilance ROR (log)",
       title="Figure 5B. External validation in Canada Vigilance") + th
ggsave(file.path(FIG,"fig5b_extval.png"), p5b, width=7.2, height=6.2, dpi=320)
cat("[fig] fig5b_extval.png\n")

## ---------- Figure 7 头对头 ----------
hh <- fread(file.path(OUT, "26_头对头_tab_headtohead.csv"))
h_l <- melt(hh, id.vars="pt",
            measure.vars=c("ROR_scopolamine","ROR_promethazine","ROR_dimenhydrinate"),
            variable.name="drug", value.name="ROR")
h_l[, drug := factor(drug, levels=c("ROR_scopolamine","ROR_promethazine","ROR_dimenhydrinate"),
                     labels=c("Scopolamine","Promethazine","Dimenhydrinate"))]
keep_pt <- hh[!is.na(ROR_promethazine) & !is.na(ROR_dimenhydrinate)][order(-ROR_scopolamine)][1:15]$pt
p7 <- ggplot(h_l[pt %in% keep_pt], aes(ROR, drug, color=drug)) +
  facet_wrap(~pt, scales="free_x", ncol=5) +
  geom_vline(xintercept=1, linetype="dashed", color="grey55") +
  geom_point(size=1.6, show.legend=FALSE) +
  scale_x_log10() + scale_color_brewer(palette="Set1") +
  labs(x="ROR (log)", y=NULL,
       title="Figure 7. Head-to-head comparison of ROR estimates",
       caption="Same pipeline; comparator drug identification by name regular expressions.") + th
ggsave(file.path(FIG,"fig7_headtohead.png"), p7, width=12, height=6.5, dpi=320)
cat("[fig] fig7_headtohead.png\n")
cat("[完成] 08_figures.R\n")
