## =====================================================================
## 08d_full_figures.R —— 出版级全家桶 v2（9 图 26 面板）
## 设计系统：Okabe-Ito；图内无标题；面板字母；600dpi；双栏 180mm≈7.09in
## =====================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(survival)})
coalesce2 <- function(a,b) ifelse(is.na(a), b, a)
OUT <- "output_v2"; FIG <- file.path(OUT, "fig_pub"); dir.create(FIG, showWarnings = FALSE)
PAL <- c(Overall="#000000", Transdermal="#0072B2", Oral="#E69F00", Parenteral="#009E73",
         replicated="#0072B2", `direction-consistent`="#E69F00", discordant="#D55E00")
theme_pub <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = .22, colour = "grey90"),
        axis.title = element_text(size = 9), axis.text = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 7.5), legend.title = element_text(size = 8),
        plot.tag = element_text(face = "bold", size = 11),
        plot.margin = margin(3, 5, 3, 3))
gsave6 <- function(nm, p, w = 7.09, h = 5.5) ggsave(file.path(FIG, nm), p, width = w, height = h, dpi = 600)

## ================= F1 研究设计 + 数据流 =================
dedup <- fread(file.path(OUT, "00_全局去重日志.csv")); gv <- function(k) as.numeric(dedup[metric==k]$value)
N_RAW <- gv("raw_demo_rows"); N_DED <- gv("post_dedup_reports"); N_REM <- gv("removed_rows")
tgt <- readRDS(file.path(OUT, "faers_target.rds"))
N_SCOP <- uniqueN(tgt[[3]][drug_group=="scopolamine"]$primaryid)

## ---- A. 设计总览（四列网格） ----
trk <- data.frame(
  x = rep(c(12.5, 37.5, 62.5, 87.5), each = 3),
  y = rep(c(80, 55, 30), 4),
  w = 23, h = 19,
  fill = c(rep("#DDEBF7",3), rep("#E2EFDA",3), rep("#FFF2CC",3), rep("#FCE4D6",3)),
  lab = c(
    "FAERS\n55 quarters\n(2012Q4-2026Q2)\n20,536,224 records",
    "Canada Vigilance\n1965-2024\n1,153,422 reports",
    "DailyMed label\nTransderm Scop\n(SPL b877a694)",
    "Cross-quarter dedup\nFDA algorithm\n17,684,540\nunique reports",
    "Scopolamine reports\n7,883 (any role)\nCanada: 431\n(114 suspect)",
    "Head-to-head set\npromethazine +\ndimenhydrinate",
    "Disproportionality\nROR primary\nPRR/IC parallel\n(a >= 3)",
    "Route-stratified\nWeibull TTO\n(fitdistr +\nsurvreg)",
    "Subgroup; label\ncomparison;\nhead-to-head",
    "External validation\n94/106 (88.7%)\ndirectional\nagreement",
    "Label-external\nincremental\nsignals;\ncode on GitHub",
    "Signal spectrum\n856 positive signals\nacross 6 clinical\nclusters"))
trk_head <- data.frame(x=c(12.5,37.5,62.5,87.5), y=95,
  lab=c("Data sources","Cohort build","Analysis modules","Validation & outputs"))
p1A <- ggplot() +
  geom_rect(data=trk, aes(xmin=x-w/2, xmax=x+w/2, ymin=y-h/2, ymax=y+h/2, fill=fill),
            colour="#5B7C99", linewidth=.35) +
  geom_text(data=trk, aes(x=x,y=y,label=lab), size=2.05, lineheight=.95) +
  geom_text(data=trk_head, aes(x=x,y=y,label=lab), fontface="bold", size=2.5) +
  scale_fill_identity() + xlim(0,100) + ylim(15,100) + theme_void()

## ---- B. PRISMA 式数据流 ----
main <- data.frame(
  x = 32, w = 52,
  y = c(90, 66, 42, 16),
  h = c(12, 10, 10, 11),
  lab = c(
    "FAERS quarterly data extracts\n55 quarters (2012Q4 - 2026Q2)\n20,536,224 report records",
    "Reports screened after\ncross-quarter deduplication\nn = 17,684,540",
    "Reports listing scopolamine\n(any drug role)\nn = 7,883",
    "Primary analysis set\n856 clinical-event signals\n(non-event PTs excluded)"))
side <- data.frame(
  x = 82, w = 34,
  y = c(77.5, 66, 29),
  h = c(12, 10, 16),
  lab = c(
    "Duplicate report records removed\n(cross-quarter deduplication:\nhighest case version per case ID)\nn = 2,851,684",
    "Reports not listing scopolamine\nexcluded\nn = 17,676,657",
    "Non-event preferred terms excluded\nfrom the primary signal set\n(OFF LABEL USE, DRUG INEFFECTIVE,\nWRONG PATIENT RECEIVED PRODUCT, etc.)"))
sa <- data.frame(x=32, y=c(84,61,37), xend=32, yend=c(71,47,21.5))
ha <- data.frame(x=c(32,58,32), y=c(77.5,66,29), xend=c(65,65,65), yend=c(77.5,66,29))
p1B <- ggplot() +
  geom_rect(data=main, aes(xmin=x-w/2,xmax=x+w/2,ymin=y-h/2,ymax=y+h/2),
            fill="#FFFFFF", colour="#2B5D8A", linewidth=.45) +
  geom_text(data=main, aes(x=x,y=y,label=lab), size=2.5, lineheight=1.0) +
  geom_rect(data=side, aes(xmin=x-w/2,xmax=x+w/2,ymin=y-h/2,ymax=y+h/2),
            fill="#F2F6FA", colour="#2B5D8A", linewidth=.4) +
  geom_text(data=side, aes(x=x,y=y,label=lab), size=2.25, lineheight=1.0) +
  geom_segment(data=sa, aes(x=x,y=y,xend=x,yend=yend),
               arrow=arrow(length=unit(2,"mm"),type="closed"), linewidth=.4, colour="grey20") +
  geom_segment(data=ha, aes(x=x,y=y,xend=xend,yend=yend),
               arrow=arrow(length=unit(2,"mm"),type="closed"), linewidth=.4, colour="grey20") +
  xlim(0,100) + ylim(5,100) + theme_void()

p1 <- p1A / p1B + plot_annotation(tag_levels="A") + plot_layout(heights=c(0.9,1.25))
gsave6("fig1_design_flow.png", p1, h=7.6)

## ================= F2 报告全景（4 面板） =================
logq <- fread(file.path(OUT, "00_季度处理日志.csv")); logq[, year := as.integer(substr(quarter,1,4))]
ann  <- logq[, .(n=sum(n_reports)), by=year]
tgt2 <- readRDS(file.path(OUT, "faers_target.rds"))
demo2 <- tgt2[[1]]; sids2 <- unique(tgt2[[3]][drug_group=="scopolamine"]$primaryid)
sq <- unique(demo2[primaryid %in% sids2, .(primaryid, quarter)])
anns <- sq[, .(n=.N), by=.(year=as.integer(substr(quarter,1,4)))]
rt   <- fread(file.path(OUT, "01b_给药途径分布.csv"))[drug_group=="scopolamine"][order(n_reports)]
rt[, route_std := factor(route_std, c("Parenteral","Other","Oral","Transdermal","Unknown/Missing"))]
rt[, unk := route_std=="Unknown/Missing"]
t1 <- fread(file.path(OUT, "22_Table1_基线特征.csv"))
age <- t1[指标 %like% "^Age group", .(g=gsub("Age group ","",指标), n=例数)]
age[, g := factor(g, c("<18","18-64",">=65"))]
sx  <- t1[指标 %like% "^Sex," & !取值 %like% "UNKNOWN", .(g=gsub("Sex, ","",指标), n=例数)]
cn  <- t1[指标 %like% "^Country,", .(g=gsub("Country, ","",指标), n=例数)][1:8][order(n)]

pA <- ggplot(ann, aes(year, n/1e6)) + geom_col(fill="#0072B2", width=.72) +
  scale_y_continuous(expand=expansion(mult=c(0,.05))) +
  labs(x=NULL, y="Library reports (millions)") + theme_pub + theme(axis.text.x=element_blank())
pB <- ggplot(anns, aes(year, n)) + geom_col(fill="#D55E00", width=.72) +
  scale_y_continuous(expand=expansion(mult=c(0,.05))) +
  labs(x="Calendar year", y="Scopolamine reports") + theme_pub
pC <- ggplot(rt, aes(n_reports, route_std, fill=unk)) + geom_col(width=.68) +
  geom_text(aes(label=sprintf("%s (%.1f%%)", format(n_reports,big.mark=","), 100*n_reports/sum(n_reports))),
            hjust=-.06, size=2.3) +
  scale_fill_manual(values=c(`FALSE`="#009E73", `TRUE`="grey72"), guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0,.40))) +
  labs(x="Scopolamine reports, by route", y=NULL) + theme_pub
pD <- ggplot(cn, aes(n, reorder(g,n))) + geom_col(fill="#8250DF", width=.68) +
  scale_x_continuous(expand=expansion(mult=c(0,.22))) +
  geom_text(aes(label=sprintf("%s (%.1f%%)", format(n,big.mark=","), 100*n/N_SCOP)), hjust=-.06, size=2.3) +
  labs(x="Reports", y=NULL) + theme_pub + theme(axis.text.y=element_text(size=7))
p2 <- (pA + plot_layout(heights=1)) / pB / (pC | pD) + plot_annotation(tag_levels="A")
gsave6("fig2_landscape.png", p2, h=8.6)

## ================= 簇分类（共用） =================
CENTRAL <- c("CONFUSIONAL STATE","HALLUCINATION","DELIRIUM","AGITATION","AMNESIA","DISORIENTATION","SOMNOLENCE","DIZZINESS","LETHARGY","SEDATION")
PERIPH  <- c("DRY MOUTH","VISION BLURRED","MYDRIASIS","URINARY RETENTION","CONSTIPATION","TACHYCARDIA","DRY EYE","ACCOMMODATION DISORDER")
HYPER   <- c("ANGIOEDEMA","SWOLLEN TONGUE","SWELLING FACE","DYSPHONIA","HYPERSENSITIVITY","ANAPHYLACTIC REACTION","URTICARIA","RASH","PRURITUS","ERYTHEMA","SWOLLEN LIP","SWELLING","STRIDOR")
HEAT    <- c("HYPERTHERMIA","PYREXIA","HEAT STROKE","HEAT EXHAUSTION","HEAT CRAMPS","HEAT ILLNESS")
cluster_of <- function(pt) {
  ifelse(pt %in% HEAT, "Heat-related",
  ifelse(grepl("^APPLICATION SITE|^PRODUCT (ADHESION|QUALITY)|DERMATITIS CONTACT|WRONG TECHNIQUE", pt), "Local/product",
  ifelse(pt %in% HYPER, "Hypersensitivity",
  ifelse(pt %in% CENTRAL, "Central anticholinergic",
  ifelse(pt %in% PERIPH, "Peripheral anticholinergic", "General/other")))))
}
CLUST_COL <- c("Central anticholinergic"="#0072B2","Peripheral anticholinergic"="#56B4E9",
               "Hypersensitivity"="#D55E00","Local/product"="#009E73","Heat-related"="#E69F00",
               "General/other"="#999999")

## ================= F3 信号谱（3 面板） =================
ov <- fread(file.path(OUT, "22b_核心PT_总体.csv"))
ov[, cluster := cluster_of(pt)]
top15 <- ov[order(-ROR)][1:15]
top15[, pt := factor(pt, levels=rev(top15$pt))]
pA <- ggplot(top15, aes(ROR, pt)) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey45", linewidth=.35) +
  geom_pointrange(aes(xmin=lcl, xmax=ucl), colour="#0072B2", linewidth=.4, size=.35) +
  scale_x_log10(breaks=c(.5,1,2,5,10,20,50), expand=expansion(mult=c(.02,.08))) +
  labs(x="ROR (95% CI, log scale)", y=NULL) + theme_pub
pB <- ggplot(ov, aes(a, ROR, colour=cluster)) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey45", linewidth=.35) +
  geom_point(size=1.7, alpha=.85) +
  ggrepel::geom_text_repel(data=ov[ (rank(-ROR)<=5 | rank(-a)<=5) & pt!="ACCOMMODATION DISORDER"],
                           aes(label=pt), size=2.1, seed=7, max.iter=20000, min.segment.length=0) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values=CLUST_COL) +
  labs(x="Reports, a", y="ROR (log scale)", colour=NULL) +
  theme_pub + theme(legend.position="right", legend.key.height=unit(.24,"cm"), legend.text=element_text(size=6.5))
cc <- ov[, .(n=.N), by=cluster][order(n)]
pC <- ggplot(cc, aes(n, reorder(cluster,n), fill=cluster)) + geom_col(width=.68) +
  geom_text(aes(label=n), hjust=-.3, size=2.4) +
  scale_fill_manual(values=CLUST_COL, guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0,.14))) +
  labs(x="Core preferred terms with positive signals", y=NULL) + theme_pub
p3 <- (pA | pB) / pC + plot_annotation(tag_levels="A") + plot_layout(heights=c(1.15,0.55))
gsave6("fig3_spectrum.png", p3, h=7.4)

## ================= F4 途径分层（热图 + 森林） =================
xw <- fread(file.path(OUT, "27_途径xPT类别交叉表.csv"))
setnames(xw, "cat", "category")
m <- melt(xw, id.vars="category", variable.name="route", value.name="n")
rn <- c(Transdermal=1751, Oral=510, Parenteral=141)
m[, pct := 100*n/rn[as.character(route)]]
m[, category := fcase(category=="Local/application-site", "Local/product",
                      category=="Hypersensitivity/systemic", "Hypersensitivity",
                      category=="Other/general", "General/other",
                      default=category)]
m[, category := factor(category, c("Local/product","Heat-related","Hypersensitivity","Central anticholinergic","Peripheral anticholinergic","General/other"))]
m[, route := factor(route, c("Transdermal","Oral","Parenteral"))]
pA <- ggplot(m, aes(route, category, fill=pct)) +
  geom_tile(colour="white", linewidth=.6) +
  geom_text(aes(label=sprintf("%d\n(%.1f%%)", n, pct)), size=2.4, lineheight=.9) +
  scale_fill_gradient(low="#F7FBFF", high="#08519C", name="% of route reports") +
  labs(x=NULL, y=NULL) + theme_pub +
  theme(axis.text.x=element_text(size=8, angle=20, hjust=1), axis.text.y=element_text(size=8),
        legend.title=element_text(size=7.5), legend.text=element_text(size=7),
        panel.grid=element_blank())
rs4 <- fread(file.path(OUT, "25_途径分层_核心PT.csv"))
keep4 <- ov[order(-ROR)][1:15]$pt
d4 <- rbind(ov[pt %in% keep4, .(pt, grp="Overall", ROR, lcl, ucl)],
            rs4[pt %in% keep4, .(pt, grp=route, ROR, lcl, ucl)])
d4[, grp := factor(grp, c("Overall","Transdermal","Oral","Parenteral"))]
d4[, pt := factor(pt, levels=rev(keep4))]
pB <- ggplot(d4[!is.na(ROR)], aes(ROR, pt, colour=grp)) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey45", linewidth=.35) +
  geom_point(position=position_dodge(width=.62), size=1.3) +
  geom_linerange(aes(xmin=lcl, xmax=ucl), position=position_dodge(width=.62), linewidth=.38) +
  scale_x_log10(breaks=c(.5,2,10,50), labels=c("0.5","2","10","50"), expand=expansion(mult=c(.02,.08))) +
  scale_colour_manual(values=PAL) +
  labs(x="ROR (95% CI, log scale)", y=NULL, colour=NULL) +
  theme_pub + theme(legend.position="bottom", legend.key.height=unit(.26,"cm"),
                    legend.text=element_text(size=7), panel.grid.major.y=element_blank())
p4 <- pA + pB + plot_annotation(tag_levels="A") + plot_layout(widths=c(1,1.25))
gsave6("fig4_route.png", p4, h=4.4)

## ================= 簇分类（共用） =================
CENTRAL <- c("CONFUSIONAL STATE","HALLUCINATION","DELIRIUM","AGITATION","AMNESIA","DISORIENTATION","SOMNOLENCE","DIZZINESS","LETHARGY","SEDATION")
PERIPH  <- c("DRY MOUTH","VISION BLURRED","MYDRIASIS","URINARY RETENTION","CONSTIPATION","TACHYCARDIA","DRY EYE","ACCOMMODATION DISORDER")
HYPER   <- c("ANGIOEDEMA","SWOLLEN TONGUE","SWELLING FACE","DYSPHONIA","HYPERSENSITIVITY","ANAPHYLACTIC REACTION","URTICARIA","RASH","PRURITUS","ERYTHEMA","SWOLLEN LIP","SWELLING","STRIDOR")
HEAT    <- c("HYPERTHERMIA","PYREXIA","HEAT STROKE","HEAT EXHAUSTION","HEAT CRAMPS","HEAT ILLNESS")
cluster_of <- function(pt) {
  ifelse(pt %in% HEAT, "Heat-related",
  ifelse(grepl("^APPLICATION SITE|^PRODUCT (ADHESION|QUALITY)|DERMATITIS CONTACT|WRONG TECHNIQUE", pt), "Local/product",
  ifelse(pt %in% HYPER, "Hypersensitivity",
  ifelse(pt %in% CENTRAL, "Central anticholinergic",
  ifelse(pt %in% PERIPH, "Peripheral anticholinergic", "General/other")))))
}
CLUST_COL <- c("Central anticholinergic"="#0072B2","Peripheral anticholinergic"="#56B4E9",
               "Hypersensitivity"="#D55E00","Local/product"="#009E73","Heat-related"="#E69F00",
               "General/other"="#999999")

## ================= F3 信号谱（3 面板） =================
ov <- fread(file.path(OUT, "22b_核心PT_总体.csv"))
ov[, cluster := cluster_of(pt)]
top15 <- ov[order(-ROR)][1:15]
top15[, pt := factor(pt, levels=rev(top15$pt))]
pA <- ggplot(top15, aes(ROR, pt)) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey45", linewidth=.35) +
  geom_pointrange(aes(xmin=lcl, xmax=ucl), colour="#0072B2", linewidth=.4, size=.35) +
  scale_x_log10(breaks=c(.5,1,2,5,10,20,50), expand=expansion(mult=c(.02,.08))) +
  labs(x="ROR (95% CI, log scale)", y=NULL) + theme_pub
pB <- ggplot(ov, aes(a, ROR, colour=cluster)) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey45", linewidth=.35) +
  geom_point(size=1.7, alpha=.85) +
  ggrepel::geom_text_repel(data=ov[ (rank(-ROR)<=5 | rank(-a)<=5) & pt!="ACCOMMODATION DISORDER"],
                           aes(label=pt), size=2.1, seed=7, max.iter=20000, min.segment.length=0) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values=CLUST_COL) +
  labs(x="Reports, a", y="ROR (log scale)", colour=NULL) +
  theme_pub + theme(legend.position="right", legend.key.height=unit(.24,"cm"), legend.text=element_text(size=6.5))
cc <- ov[, .(n=.N), by=cluster][order(n)]
pC <- ggplot(cc, aes(n, reorder(cluster,n), fill=cluster)) + geom_col(width=.68) +
  geom_text(aes(label=n), hjust=-.3, size=2.4) +
  scale_fill_manual(values=CLUST_COL, guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0,.14))) +
  labs(x="Core preferred terms with positive signals", y=NULL) + theme_pub
p3 <- (pA | pB) / pC + plot_annotation(tag_levels="A") + plot_layout(heights=c(1.15,0.55))
gsave6("fig3_spectrum.png", p3, h=7.4)

## ================= F4 途径分层（热图 + 森林） =================
xw <- fread(file.path(OUT, "27_途径xPT类别交叉表.csv"))
setnames(xw, "cat", "category")
m <- melt(xw, id.vars="category", variable.name="route", value.name="n")
rn <- c(Transdermal=1751, Oral=510, Parenteral=141)
m[, pct := 100*n/rn[as.character(route)]]
m[, category := fcase(category=="Local/application-site", "Local/product",
                      category=="Hypersensitivity/systemic", "Hypersensitivity",
                      category=="Other/general", "General/other",
                      default=category)]
m[, category := factor(category, c("Local/product","Heat-related","Hypersensitivity","Central anticholinergic","Peripheral anticholinergic","General/other"))]
m[, route := factor(route, c("Transdermal","Oral","Parenteral"))]
pA <- ggplot(m, aes(route, category, fill=pct)) +
  geom_tile(colour="white", linewidth=.6) +
  geom_text(aes(label=sprintf("%d\n(%.1f%%)", n, pct)), size=2.4, lineheight=.9) +
  scale_fill_gradient(low="#F7FBFF", high="#08519C", name="% of route reports") +
  labs(x=NULL, y=NULL) + theme_pub +
  theme(axis.text.x=element_text(size=8, angle=20, hjust=1), axis.text.y=element_text(size=8),
        legend.title=element_text(size=7.5), legend.text=element_text(size=7),
        panel.grid=element_blank())
rs4 <- fread(file.path(OUT, "25_途径分层_核心PT.csv"))
keep4 <- ov[order(-ROR)][1:15]$pt
d4 <- rbind(ov[pt %in% keep4, .(pt, grp="Overall", ROR, lcl, ucl)],
            rs4[pt %in% keep4, .(pt, grp=route, ROR, lcl, ucl)])
d4[, grp := factor(grp, c("Overall","Transdermal","Oral","Parenteral"))]
d4[, pt := factor(pt, levels=rev(keep4))]
pB <- ggplot(d4[!is.na(ROR)], aes(ROR, pt, colour=grp)) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey45", linewidth=.35) +
  geom_point(position=position_dodge(width=.62), size=1.3) +
  geom_linerange(aes(xmin=lcl, xmax=ucl), position=position_dodge(width=.62), linewidth=.38) +
  scale_x_log10(breaks=c(.5,2,10,50), labels=c("0.5","2","10","50"), expand=expansion(mult=c(.02,.08))) +
  scale_colour_manual(values=PAL) +
  labs(x="ROR (95% CI, log scale)", y=NULL, colour=NULL) +
  theme_pub + theme(legend.position="bottom", legend.key.height=unit(.26,"cm"),
                    legend.text=element_text(size=7), panel.grid.major.y=element_blank())
p4 <- pA + pB + plot_annotation(tag_levels="A") + plot_layout(widths=c(1,1.25))
gsave6("fig4_route.png", p4, h=4.4)

## ================= F5 TTO（ECDF+CI 带 + 箱线） =================
tto <- readRDS(file.path(OUT, "23b_tto_报告级.rds"))
tto[, route_std := factor(fifelse(is.na(route_std),"Overall",route_std), c("Overall","Transdermal","Oral","Parenteral"))]
grps <- c("Overall","Transdermal","Oral","Parenteral")
tto_all <- rbind(tto[, .(tto, route_std = factor("Overall", levels = grps))], tto, fill = TRUE)
fits <- lapply(grps, function(g){
  v <- tto_all[route_std==g]$tto
  f <- survival::survreg(survival::Surv(v)~1, dist="weibull")
  data.table(grp=g, mu=unname(coef(f)), sigma=f$scale,
             V11=f$var[1,1], V12=f$var[1,2], V22=f$var[2,2], med=median(v))
})
fits <- rbindlist(fits)
tt <- data.table(expand.grid(t=1:365, grp=grps))
cdf <- fits[tt, on=.(grp), allow.cartesian=TRUE]
## pointwise 95% CI band via delta method on (mu, log sigma)
cdf[, lt := log(t)]
cdf[, y  := (lt - mu)/sigma]
cdf[, g1 := -1/sigma]
cdf[, g2 := -(lt - mu)/sigma]
cdf[, se := sqrt(g1^2*V11 + g2^2*V22 + 2*g1*g2*V12)]
cdf[, est := 1-exp(-exp(y))]
cdf[, lo  := 1-exp(-exp(y-1.96*se))]
cdf[, hi  := 1-exp(-exp(y+1.96*se))]
pA <- ggplot(cdf, aes(t, est, colour=grp)) +
  geom_ribbon(aes(ymin=lo, ymax=hi, fill=grp), alpha=.12, colour=NA) +
  geom_line(linewidth=.7) +
  coord_cartesian(xlim=c(0,365), expand=FALSE) +
  scale_x_continuous(breaks=seq(0,365,90)) +
  scale_colour_manual(values=PAL, breaks=grps) + scale_fill_manual(values=PAL, breaks=grps) +
  labs(x="Time-to-onset (days)", y="Cumulative proportion (Weibull)", colour=NULL, fill=NULL) +
  theme_pub + theme(legend.position=c(.6,.28), legend.background=element_rect(fill="white", colour="grey85", linewidth=.2),
                    legend.key.height=unit(.26,"cm"))
pB <- ggplot(tto_all[!is.na(route_std)], aes(route_std, tto, fill=route_std)) +
  geom_boxplot(outlier.size=.4, outlier.alpha=.35, linewidth=.35, width=.55) +
  scale_y_log10(breaks=c(1,7,30,90,180,365,1000,3000)) +
  scale_fill_manual(values=PAL, breaks=grps, guide="none") +
  labs(x=NULL, y="Time-to-onset (days, log scale)") + theme_pub
p5 <- pA / pB + plot_annotation(tag_levels="A")
gsave6("fig5_tto.png", p5, h=6.4)

## ================= F6 晕动症亚组（散点 + 棒棒糖） =================
ms  <- fread(file.path(OUT, "04_敏感性分析_晕动症适应证.csv"))[drug %like% "^scopolamine"]
m6 <- merge(ms[, .(pt, ms_a=a, ms_ROR=ROR, ms_sig=signal)], ov[, .(pt, ov_a=a, ov_ROR=ROR)], by="pt")
top6 <- m6[order(-ms_a)][1:6]$pt
pA <- ggplot(m6, aes(ov_ROR, ms_ROR)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50", linewidth=.35) +
  geom_point(aes(size=ms_a, colour=ms_sig), alpha=.85) +
  ggrepel::geom_text_repel(data=m6[pt %in% top6], aes(label=pt), size=2.3, seed=2, max.overlaps=20) +
  scale_x_log10(limits=c(.8,400)) + scale_y_log10(limits=c(.8,400)) +
  scale_colour_manual(values=c(`TRUE`="#0072B2", `FALSE`="grey62"),
                      labels=c(`TRUE`="Positive signal", `FALSE`="Not significant")) +
  scale_size_continuous(range=c(1.2,4.5), guide="none") +
  annotate("text", x=1.0, y=280, hjust=0, size=2.6,
           label=sprintf("Shared PTs: %d\nSubgroup reports: n = 554 (7.0%%)", nrow(m6))) +
  labs(x="ROR, overall (log scale)", y="ROR, motion sickness subgroup (log scale)", colour=NULL) +
  theme_pub + theme(legend.position=c(.04,.14),
                    legend.background=element_rect(fill="white", colour="grey85", linewidth=.2),
                    legend.key.height=unit(.26,"cm"))
top10 <- ms[signal==TRUE][order(-a)][1:10]
top10[, pt := factor(pt, levels=rev(top10$pt))]
pB <- ggplot(top10, aes(ROR, pt)) +
  geom_segment(aes(x=0, xend=ROR, yend=pt), colour="grey75", linewidth=.5) +
  geom_point(aes(size=a), colour="#0072B2") +
  geom_text(data=top10[ROR < 150], aes(label=sprintf("ROR %.1f", ROR)), hjust=-.15, size=2.3, colour="grey25") +
  geom_text(data=top10[ROR >= 150], aes(label=sprintf("ROR %.1f", ROR)), hjust=1.08, size=2.3, colour="grey25") +
  scale_x_log10(limits=c(1,320), breaks=c(1,10,50,200)) +
  scale_size_continuous(range=c(1.5,4), guide="none") +
  labs(x="ROR in motion sickness subgroup (log scale)", y=NULL) +
  theme_pub + theme(legend.position="none")
p6 <- (pA | pB) + plot_annotation(tag_levels="A")
gsave6("fig6_mssubgroup.png", p6, h=4.6)

## ================= F7 外验（散点 + 哑铃） =================
ev <- fread(file.path(OUT, "21_两库对比_剔除非事件.csv"))
evs <- ev[faers_signal==TRUE & !is.na(canada_ror)]
agree <- evs[verdict %in% c("replicated","direction-consistent"), .N]/evs[, .N]
rep10 <- evs[verdict=="replicated"][order(-faers_a)][1:10]
dm <- melt(rep10[, .(pt, faers_ror, canada_ror)], id.vars="pt", variable.name="db", value.name="ROR")
dm[, db := factor(db, c("canada_ror","faers_ror"), c("Canada Vigilance","FAERS"))]
dm[, pt := factor(pt, levels=rep10$pt)]
pA <- ggplot(evs, aes(faers_ror, canada_ror, colour=verdict)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50", linewidth=.35) +
  geom_point(aes(size=faers_a), alpha=.85) +
  ggrepel::geom_text_repel(data=evs[pt %in% c("MYDRIASIS","DRY MOUTH","DYSPHONIA","ANGIOEDEMA","NAUSEA")],
                           aes(label=pt), size=2.3, seed=1, max.overlaps=20, show.legend=FALSE) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values=PAL, breaks=c("replicated","direction-consistent","discordant"),
                      labels=c("Replicated","Direction-consistent","Discordant")) +
  scale_size_continuous(range=c(1.2,4.2), guide="none") +
  annotate("text", x=min(evs$faers_ror)*1.1, y=max(evs$canada_ror)*.75, hjust=0, size=2.6,
           label=sprintf("Directional agreement:\n%d/%d (%.0f%%)", round(agree*evs[,.N]), evs[,.N], 100*agree)) +
  labs(x="FAERS ROR (log scale)", y="Canada Vigilance ROR (log scale)", colour=NULL) +
  theme_pub + theme(legend.position="bottom", legend.key.height=unit(.26,"cm"), legend.text=element_text(size=7))
pB <- ggplot(dm, aes(ROR, pt, colour=db)) +
  geom_line(aes(group=pt), colour="grey78", linewidth=.8) +
  geom_point(aes(shape=db), size=2.1, alpha=.92) +
  scale_x_log10(breaks=c(1,5,10,20,50)) +
  scale_colour_manual(values=c(`Canada Vigilance`="#D55E00", FAERS="#0072B2")) +
  scale_shape_manual(values=c(`Canada Vigilance`=17, FAERS=16)) +
  labs(x="ROR (log scale)", y=NULL, colour=NULL, shape=NULL) +
  theme_pub + theme(legend.position="top", legend.key.height=unit(.26,"cm"),
                    panel.grid.major.y=element_blank())
p7 <- (pA | pB) + plot_annotation(tag_levels="A")
gsave6("fig7_extval.png", p7, h=4.6)

## ================= F8 头对头 =================
hh <- fread(file.path(OUT, "26_头对头_tab_headtohead.csv"))
h12 <- hh[!is.na(ROR_promethazine) & !is.na(ROR_dimenhydrinate)][order(-ROR_scopolamine)][1:12]
h_l <- rbind(h12[, .(pt, ROR=ROR_scopolamine,   drug="Scopolamine")],
             h12[, .(pt, ROR=ROR_promethazine, drug="Promethazine")],
             h12[, .(pt, ROR=ROR_dimenhydrinate, drug="Dimenhydrinate")])
h_l[, pt := factor(pt, levels=rev(h12$pt))]
h_l[, drug := factor(drug, c("Promethazine","Dimenhydrinate","Scopolamine"))]
p8 <- ggplot(h_l, aes(ROR, pt, colour=drug, shape=drug)) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey45", linewidth=.35) +
  geom_point(size=1.9, alpha=.92) +
  scale_x_log10(breaks=c(1,5,10,20,50)) +
  scale_colour_manual(values=c(Scopolamine="#000000", Promethazine="#E69F00", Dimenhydrinate="#009E73")) +
  scale_shape_manual(values=c(Scopolamine=16, Promethazine=1, Dimenhydrinate=2)) +
  labs(x="ROR (log scale)", y=NULL, colour=NULL, shape=NULL) +
  theme_pub + theme(legend.position="top", legend.key.height=unit(.26,"cm"),
                    panel.grid.major.y=element_blank())
gsave6("fig8_headtohead.png", p8, w=6.4, h=5.0)

## ================= F9 label（点图 + 簇×状态堆叠） =================
lab <- fread("材料包_EN/tab_label_XuCai_260904.csv", encoding="UTF-8")
setnames(lab, c("PT","FAERS ROR (95% CI)","Label mention (Y/N/U)"), c("pt","ror_txt","lab"))
lab[, ror := as.numeric(sub("^([0-9.]+).*", "\\1", ror_txt))]
lab <- lab[!is.na(ror)]
lab[, cluster := cluster_of(pt)]
lab[, lab := factor(lab, c("Yes","Unclear","No"))]
lab[, pt := factor(pt, levels = pt[order(lab, -ror)])]
nlab <- lab[, .N, by=lab]; llab <- sprintf("%s (n=%d)", nlab$lab, nlab$N); names(llab) <- as.character(nlab$lab)
pA <- ggplot(lab, aes(ror, pt, colour=lab)) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey45", linewidth=.35) +
  geom_point(size=1.8, alpha=.9) +
  scale_x_log10(breaks=c(1,2,5,10,20,50)) +
  scale_colour_manual(values=c(Yes="#0072B2", Unclear="#999999", No="#D55E00"), labels=llab, name="Label status") +
  labs(x="ROR (log scale)", y=NULL) +
  theme_pub + theme(legend.position="top", axis.text.y=element_text(size=6.2),
                    legend.key.height=unit(.26,"cm"))
cs <- unique(lab[, .(pt, cluster, lab)])
cs_sum <- cs[, .(n=.N), by=.(cluster, lab)]
cs_sum[, cluster := factor(cluster, names(CLUST_COL))]
pB <- ggplot(cs_sum, aes(cluster, n, fill=lab)) +
  geom_col(position="stack", width=.7, colour="white", linewidth=.3) +
  scale_fill_manual(values=c(Yes="#0072B2", Unclear="#999999", No="#D55E00"), name="Label status") +
  scale_x_discrete(labels=function(x) gsub(" anticholinergic","",x)) +
  labs(x=NULL, y="Preferred terms") +
  theme_pub + theme(legend.position="top", axis.text.x=element_text(size=7, angle=20, hjust=1),
                    legend.key.height=unit(.26,"cm"))
p9 <- pA / pB + plot_annotation(tag_levels="A") + plot_layout(heights=c(2.5,1))
gsave6("fig9_label.png", p9, h=8.0)

cat("[完成] 9 张出版级全家桶 →", FIG, "\n")
