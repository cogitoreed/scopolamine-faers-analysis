## =====================================================================
## 07_w1prime_tables.R —— W1' 补表（v3.1）
##   1) 重读 DEMO（按各季度自身表头，修复 faers_target.rds 的列错位）
##   2) Table 1 基线特征
##   3) TTO 重建 + survreg 交叉验证（v3.1 I 组）+ heat 类 PT 的 TTO（v3.1 C 组）
##   4) 途径分层核心 PT 信号（v3.1 E 组）+ 途径×PT 类别交叉表
##   5) 头对头（v3.1 F 组，基于 02 号全药物文件）
##   6) 摘要数字核对 facts
## 运行：LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 R_LIBS_USER=~/Rlibs Rscript 07_w1prime_tables.R
## =====================================================================
suppressMessages({library(data.table); library(MASS); library(survival)})
coalesce2 <- function(a, b) ifelse(is.na(a), b, a)
OUT <- "output_v2"
QDIR <- "data/extracted"

x <- readRDS(file.path(OUT, "faers_target.rds"))
demo_rds <- x[[1]]; pt_rds <- x[[2]]; drug_rds <- x[[3]]; outc_rds <- x[[4]]
ther_rds <- x[[5]]; indi_rds <- x[[6]]; libpt_rds <- x[[7]]

pt_total <- fread(file.path(OUT, "00_全库PT计数.csv"))
names(pt_total) <- c("pt", "n_total_pt")
## 修复表头被当数据的行（首行 pt=="PT", n=55 为文件写入时的表头残行）
pt_total <- pt_total[pt != "PT" | n_total_pt != 55]
pt_total[, pt := toupper(pt)]
setkey(pt_total, pt)

N_LIB <- sum(fread(file.path(OUT, "00_季度处理日志.csv"))$n_reports)  ## 去重后全库（v2 为全局去重口径）
stopifnot(N_LIB > 0)

## ---------------------------------------------------------------------
## 1. 重读 DEMO：按各季度文件自身表头
## ---------------------------------------------------------------------
demo_files <- list.files(QDIR, pattern="^demo.*\\.txt$", ignore.case=TRUE,
                         recursive=TRUE, full.names=TRUE)
stopifnot(length(demo_files) == 55)
target_ids <- demo_rds$primaryid  ## 112,773 目标药报告（5 药）

read_demo_one <- function(f) {
  hdr <- strsplit(toupper(sub("\\s+", "", readLines(f, n=1, warn=FALSE))), "$", fixed=TRUE)[[1]]
  old <- "GNDR_COD" %in% hdr
  want <- if (old) c("primaryid","event_dt","age","age_cod","gndr_cod","occp_cod","reporter_country","occr_country")
          else     c("primaryid","event_dt","age","age_cod","sex","occp_cod","reporter_country","occr_country")
  sel <- intersect(want, hdr)
  dt <- fread(f, sep="$", quote="", header=TRUE, select=sel,
              colClasses="character", na.strings=c(""), showProgress=FALSE)
  if (old) setnames(dt, "gndr_cod", "sex")
  dt[]
}
demo_fix <- rbindlist(lapply(demo_files, read_demo_one), fill=TRUE)
demo_fix <- demo_fix[primaryid %in% target_ids]
cat(sprintf("[DEMO 重读] 目标报告 %d / %d（应相等）\n", nrow(demo_fix), length(target_ids)))
stopifnot(nrow(demo_fix) == length(target_ids))
demo_fix[, event_dt := as.Date(event_dt, format="%Y%m%d")]
demo_fix[, age_n := suppressWarnings(as.numeric(age))]
demo_fix[, sex_u := ifelse(toupper(ifelse(is.na(sex),"",sex)) %in% c("M","F","U","NS"), toupper(sex), NA_character_)]
## 年龄统一换算为岁（handbook A2：YR×1、MON/12、WK/52、DY/365、DEC×10、HR/8760）
demo_fix[, age_yr := fcase(
  toupper(ifelse(is.na(age_cod),"",age_cod)) == "YR", age_n,
  toupper(ifelse(is.na(age_cod),"",age_cod)) == "MON", age_n/12,
  toupper(ifelse(is.na(age_cod),"",age_cod)) == "WK",  age_n/52,
  toupper(ifelse(is.na(age_cod),"",age_cod)) == "DY",  age_n/365,
  toupper(ifelse(is.na(age_cod),"",age_cod)) == "DEC", age_n*10,
  toupper(ifelse(is.na(age_cod),"",age_cod)) == "HR",  age_n/8760,
  is.na(age_n), NA_real_,
  default = NA_real_)]
demo_fix[!is.na(age_yr) & (age_yr < 0 | age_yr > 110), age_yr := NA_real_]
cat(sprintf("[DEMO 修复率] age 可用 %.1f%%  sex 可用 %.1f%%\n",
            100*mean(!is.na(demo_fix$age_yr)), 100*mean(!is.na(demo_fix$sex_u))))

## ---------------------------------------------------------------------
## 2. 东莨菪碱主分析集（any-role，n=9,548）+ 每报告途径 + 晕动症适应证
## ---------------------------------------------------------------------
scop_ids <- unique(drug_rds[drug_group=="scopolamine"]$primaryid)
N_SCOP <- length(scop_ids)
cat(sprintf("[主分析集] 东莨菪碱 any-role 报告 n=%d（v1 未全局去重时为 9,548）\n", N_SCOP))
role_dist <- unique(drug_rds[drug_group=="scopolamine", .(primaryid, role_cod)])[,
              .(n_reports = uniqueN(primaryid)), by=role_cod][order(-n_reports)]

## 每报告途径：与主管线 std_route/RT 逐字一致（unique 取首条记录）
std_route <- function(r) {
  r <- toupper(trimws(as.character(r)))
  out <- rep("Unknown/Missing", length(r))
  out[grepl("TRANSDERM|PATCH", r)] <- "Transdermal"
  out[grepl("TOPICAL|CUTANEOUS", r) & out == "Unknown/Missing"] <- "Transdermal"
  out[grepl("ORAL|MOUTH|SUBLINGUAL|^PO$", r) & out == "Unknown/Missing"] <- "Oral"
  out[grepl("INTRAVEN|INTRAMUSC|SUBCUT|INJECT|INFUSION|^IV$|^IM$|^SC$", r) &
      out == "Unknown/Missing"] <- "Parenteral"
  out[r != "" & r != "UNKNOWN" & out == "Unknown/Missing"] <- "Other"
  out
}
RT <- unique(drug_rds[drug_group=="scopolamine", .(primaryid, route)], by="primaryid")
RT[, route_std := std_route(route)]
RT_u <- RT
rt_map <- RT_u[, .(primaryid, route_std)]
route_chk <- RT_u[, .(n_reports = .N), by=route_std][order(-n_reports)]
cat("[途径核对]\n"); print(route_chk)
ref <- fread(file.path(OUT, "01b_给药途径分布.csv"))[drug_group=="scopolamine"]
rc <- copy(route_chk); setkey(rc, route_std); ref2 <- copy(ref); setkey(ref2, route_std)
cmp <- ref2[rc[, .(route_std, n_reports)]]
if (any(cmp$n_reports != cmp$i.n_reports)) { cat("[警告] 途径与 01b 不一致：\n"); print(cmp) } else cat("[途径与 01b 完全一致]\n")

MS_PAT <- "MOTION SICKNESS|MOTION SICKNESS PROPHYLAXIS|TRAVEL SICKNESS|SEASICKNESS|KINETOSIS"
ms_ids <- unique(indi_rds[grepl(MS_PAT, indi_pt)]$primaryid)
N_MS <- sum(scop_ids %in% ms_ids)
cat(sprintf("[晕动症适应证] 东莨菪碱报告 n=%d（%.1f%%）\n", N_MS, 100*N_MS/N_SCOP))

## 严重结局（断言 vs 07 文件）
oc <- unique(outc_rds[primaryid %in% scop_ids], by=c("primaryid","outc_cod"))
ser <- oc[, .(death = any(outc_cod=="DE"), hosp = any(outc_cod=="HO"),
              disab = any(outc_cod=="DS"), life = any(outc_cod=="CA")), by=primaryid]
cat(sprintf("[严重结局] death %d/%d (%.2f%%) hosp %d (%.2f%%) life %d disab %d（与 07_严重结局构成.csv 核对）\n",
            sum(ser$death), N_SCOP, 100*sum(ser$death)/N_SCOP, sum(ser$hosp), 100*sum(ser$hosp)/N_SCOP,
            sum(ser$life), sum(ser$disab)))

demo_scop <- demo_fix[primaryid %in% scop_ids]

## ---------------------------------------------------------------------
## 3. Table 1
## ---------------------------------------------------------------------
med_iqr <- function(v) { v <- v[!is.na(v)]; if (!length(v)) return(c(NA,NA,NA,0))
  c(round(median(v),1), round(quantile(v,.25,na.rm=TRUE),1), round(quantile(v,.75,na.rm=TRUE),1), length(v)) }
a_m <- med_iqr(demo_scop$age_yr)

t1 <- list()
t1[[1]] <- data.table(指标="Reports (scopolamine, any role)", 取值=sprintf("%d", N_SCOP), 例数=N_SCOP, 占比=100)
t1[[2]] <- data.table(指标="Age, years, median (Q1, Q3)", 取值=sprintf("%s (%s, %s)", a_m[1], a_m[2], a_m[3]), 例数=a_m[4], 占比=round(100*a_m[4]/N_SCOP,1))
t1[[3]] <- data.table(指标="Age missing", 取值=sprintf("%d", N_SCOP-a_m[4]), 例数=N_SCOP-a_m[4], 占比=round(100*(N_SCOP-a_m[4])/N_SCOP,1))
ag <- demo_scop[!is.na(age_yr), .(n=.N), by=.(grp=cut(age_yr, c(-Inf,18,65,Inf), labels=c("<18","18-64",">=65"), right=FALSE))]
for (i in seq_len(nrow(ag))) t1[[length(t1)+1]] <- data.table(指标=paste0("Age group ", ag$grp[i]), 取值=as.character(ag$n[i]), 例数=ag$n[i], 占比=round(100*ag$n[i]/N_SCOP,1))
sx <- demo_scop[, .(n=.N), by=.(s=toupper(ifelse(is.na(sex_u),"Unknown/missing",sex_u)))][order(-n)]
for (i in seq_len(nrow(sx))) t1[[length(t1)+1]] <- data.table(指标=paste0("Sex, ", sx$s[i]), 取值=as.character(sx$n[i]), 例数=sx$n[i], 占比=round(100*sx$n[i]/N_SCOP,1))
for (i in seq_len(nrow(route_chk))) t1[[length(t1)+1]] <- data.table(指标=paste0("Route, ", route_chk$route_std[i]), 取值=as.character(route_chk$n_reports[i]), 例数=route_chk$n_reports[i], 占比=round(100*route_chk$n_reports[i]/N_SCOP,1))
demo_scop[, cc1 := ifelse(is.na(occr_country), reporter_country, occr_country)]
demo_scop[, cc1 := toupper(ifelse(is.na(cc1),"UNKNOWN",cc1))]
cn <- demo_scop[, .(n=.N), by=.(c=cc1)][order(-n)][1:10]
for (i in seq_len(nrow(cn))) t1[[length(t1)+1]] <- data.table(指标=paste0("Country, ", cn$c[i]), 取值=as.character(cn$n[i]), 例数=cn$n[i], 占比=round(100*cn$n[i]/N_SCOP,1))
ocp <- demo_scop[, .(n=.N), by=.(o=toupper(ifelse(is.na(occp_cod),"NA",occp_cod)))][order(-n)][1:8]
for (i in seq_len(nrow(ocp))) t1[[length(t1)+1]] <- data.table(指标=paste0("Reporter, ", ocp$o[i]), 取值=as.character(ocp$n[i]), 例数=ocp$n[i], 占比=round(100*ocp$n[i]/N_SCOP,1))
yr <- demo_scop[!is.na(event_dt), .(n=.N), by=.(y=year(event_dt))][order(y)]
for (i in seq_len(nrow(yr))) t1[[length(t1)+1]] <- data.table(指标=paste0("Year of onset, ", yr$y[i]), 取值=as.character(yr$n[i]), 例数=yr$n[i], 占比=round(100*yr$n[i]/N_SCOP,1))
t1[[length(t1)+1]] <- data.table(指标="Serious outcome: death", 取值=as.character(sum(ser$death)), 例数=sum(ser$death), 占比=round(100*sum(ser$death)/N_SCOP,1))
t1[[length(t1)+1]] <- data.table(指标="Serious outcome: hospitalization", 取值=as.character(sum(ser$hosp)), 例数=sum(ser$hosp), 占比=round(100*sum(ser$hosp)/N_SCOP,1))
t1[[length(t1)+1]] <- data.table(指标="Serious outcome: life-threatening", 取值=as.character(sum(ser$life)), 例数=sum(ser$life), 占比=round(100*sum(ser$life)/N_SCOP,1))
t1[[length(t1)+1]] <- data.table(指标="Serious outcome: disability", 取值=as.character(sum(ser$disab)), 例数=sum(ser$disab), 占比=round(100*sum(ser$disab)/N_SCOP,1))
t1[[length(t1)+1]] <- data.table(指标="Motion sickness indication (INDI)", 取值=as.character(N_MS), 例数=N_MS, 占比=round(100*N_MS/N_SCOP,1))
t1[[length(t1)+1]] <- data.table(指标="Scopolamine drug role: PS/SS/C/I records", 取值=paste(role_dist$role_cod, role_dist$n_reports, sep="=", collapse="; "), 例数=sum(role_dist$n_reports), 占比=NA)
tab1 <- rbindlist(t1, fill=TRUE)
fwrite(tab1, file.path(OUT, "22_Table1_基线特征.csv"))
cat("[写出] 22_Table1_基线特征.csv\n")

## ---------------------------------------------------------------------
## 4. TTO 重建 + fitdistr vs survreg（v3.1 I 组）+ heat PT TTO（C 组）
## ---------------------------------------------------------------------
## TTO（v2 正确口径）：锁定东莨菪碱自身记录的 drug_seq ↔ THER.dsg_drug_seq，
## 取该药最早 start_dt（支持 8/6/4 位部分日期，同主管线 parse_fda_date），
## tto = event_dt - start_dt，0 < tto ≤ 3650 天（与主管线上限一致）
## 日期处理（v2.1 对齐 Ozawa 2022 惯例）：YYYYMMDD 原样；YYYYMM → 当月15日（月中插补）；
## 仅 YYYY 的记录剔除（避免长达 ±364 天的插补误差）
parse_date_mid <- function(x) {
  x <- trimws(as.character(x)); res <- as.Date(rep(NA_character_, length(x)))
  ok <- !is.na(x) & nzchar(x)
  if (any(ok)) {
    v <- x[ok]; nc <- nchar(v)
    v[nc == 8] <- v[nc == 8]
    v[nc == 6] <- paste0(v[nc == 6], "15")
    d <- as.Date(v[nc %in% c(8,6)], format = "%Y%m%d")
    d[!is.na(d) & (d < as.Date("1900-01-01") | d > Sys.Date())] <- NA
    res[ok][nc %in% c(8,6)] <- d
  }
  res
}
scop_rec <- unique(drug_rds[drug_group=="scopolamine", .(primaryid, drug_seq)])
th <- ther_rds[primaryid %in% scop_ids]
th <- merge(scop_rec, th, by.x=c("primaryid","drug_seq"), by.y=c("primaryid","dsg_drug_seq"))
th[, sdt := parse_date_mid(start_dt)]
th_u <- th[!is.na(sdt), .(sdt = min(sdt)), by=primaryid]   ## 目标药最早 start_dt
tto_dt <- merge(th_u, demo_scop[, .(primaryid, event_dt)], by="primaryid")
tto_dt[, tto := as.numeric(event_dt - sdt)]
tto_dt <- tto_dt[!is.na(tto) & tto > 0 & tto <= 3650]
tto_dt <- merge(tto_dt, rt_map, by="primaryid", all.x=TRUE)   ## 每报告途径
tto_v <- tto_dt$tto
saveRDS(tto_dt, file.path(OUT, "23b_tto_报告级.rds"))  ## 供画 TTO 曲线
cat(sprintf("[TTO] 锁定 drug_seq 口径可用 n=%d，median=%.0f 天\n", length(tto_v), median(tto_v)))

fit_w <- function(v) {
  r <- tryCatch({
    f <- MASS::fitdistr(v, "weibull")
    b <- unname(f$estimate["shape"]); se_b <- unname(f$sd["shape"])
    list(method="fitdistr", n=length(v), shape=b, shape_lcl=b-1.96*se_b, shape_ucl=b+1.96*se_b,
         scale=unname(f$estimate["scale"]))
  }, error=function(e) list(method="fitdistr", n=length(v), shape=NA_real_, shape_lcl=NA_real_,
                            shape_ucl=NA_real_, scale=NA_real_))
  r
}
fit_sr <- function(v) {
  d <- data.table(tto=v)
  f <- survreg(Surv(tto) ~ 1, data=d, dist="weibull")
  sc <- f$scale; b <- 1/sc; se_ls <- sqrt(f$var[1,1])
  list(method="survreg", n=length(v), shape=b, shape_lcl=b*exp(-1.96*se_ls), shape_ucl=b*exp(1.96*se_ls),
       scale=sc)
}
cmp_fit <- function(v, label) {
  a <- fit_w(v); b <- fit_sr(v)
  rel <- abs(a$shape_lcl - a$shape_ucl - (b$shape_lcl - b$shape_ucl)) / max(a$shape_lcl - a$shape_ucl, b$shape_lcl - b$shape_ucl)
  r <- rbind(as.data.table(a), as.data.table(b))
  r[, label := label]; r[, median_tto := round(median(v),1)]
  r[, ci_width := shape_ucl - shape_lcl]
  r[]
}
res_all <- cmp_fit(tto_v, "overall")
res_route <- rbindlist(lapply(c("Transdermal","Oral","Parenteral"), function(g) {
  v <- tto_dt[route_std==g & tto > 0]$tto
  if (length(v) >= 20) cmp_fit(v, g) else NULL
}))
tto_ci <- rbind(res_all, res_route, fill=TRUE)
fwrite(tto_ci, file.path(OUT, "23_TTO_重拟_CI.csv"))
cat("[写出] 23_TTO_重拟_CI.csv\n"); print(tto_ci)

## heat 类 PT（C 组）：计数、ROR、TTO
HEAT_PTS <- c("HYPERTHERMIA","PYREXIA","HEAT STROKE","HEAT EXHAUSTION","HEAT CRAMPS","HEAT ILLNESS")
scop_pt <- unique(pt_rds[primaryid %in% scop_ids], by=c("primaryid","pt"))
heat_tab <- rbindlist(lapply(HEAT_PTS, function(p) {
  a <- sum(scop_pt$pt == p)
  c <- coalesce2(pt_total[.(p)]$n_total_pt, 0)
  if (a < 3 || c <= a) return(data.table(pt=p, a=a, ROR=NA_real_, lcl=NA_real_, ucl=NA_real_))
  b <- N_SCOP - a; d <- N_LIB - N_SCOP - (c - a); cc <- c - a
  ror <- (a/b)/(cc/d); se <- sqrt(1/a+1/b+1/cc+1/d)
  data.table(pt=p, a=a, ROR=ror, lcl=exp(log(ror)-1.96*se), ucl=exp(log(ror)+1.96*se))
}))
heat_ids <- scop_pt[pt %in% HEAT_PTS]$primaryid
v_heat <- tto_dt[primaryid %in% heat_ids & tto > 0]$tto
heat_extra <- data.table(
  n_reports_with_heat_pt = uniqueN(heat_ids),
  n_tto_usable = length(v_heat),
  pct_le_3days = ifelse(length(v_heat), round(100*mean(v_heat <= 3),1), NA),
  median_tto = ifelse(length(v_heat), round(median(v_heat),1), NA))
fwrite(heat_tab, file.path(OUT, "24_heat_PT_ROR.csv"))
fwrite(heat_extra, file.path(OUT, "24b_heat_PT_TTO汇总.csv"))
cat("[写出] 24_heat_PT_ROR.csv / 24b\n"); print(heat_tab); print(heat_extra)
if (length(v_heat) >= 20) {
  h <- rbind(cmp_fit(v_heat, "heat-PT set")); fwrite(h, file.path(OUT, "24c_heat_TTO_拟合.csv"))
  print(h)
} else cat(sprintf("[heat TTO] n=%d < 20，仅报告描述统计（≤72h 占比 %.1f%%）\n", length(v_heat), heat_extra$pct_le_3days))

## ---------------------------------------------------------------------
## 5. 核心 PT 清单 + 途径分层（v3.1 E 组）+ 途径×类别交叉表
## ---------------------------------------------------------------------
NON_EVENT <- c("OFF LABEL USE","DRUG INEFFECTIVE","PRODUCT USE IN UNAPPROVED INDICATION")
sig03 <- fread(file.path(OUT, "03_东莨菪碱_阳性信号.csv"))
t20 <- fread(file.path(OUT, "20_Table2_主分析信号_Top25.csv"))
cat("非事件 PT 实际集合:", paste(setdiff(sig03$pt, fread(file.path(OUT,"21_两库对比_剔除非事件.csv"))$pt), collapse=" | "), "\n")

top25 <- t20$`Preferred term (PT)`
CENTRAL <- c("CONFUSIONAL STATE","HALLUCINATION","DELIRIUM","AGITATION","AMNESIA","DISORIENTATION","SOMNOLENCE","DIZZINESS","LETHARGY","SEDATION")
PERIPH  <- c("DRY MOUTH","VISION BLURRED","MYDRIASIS","URINARY RETENTION","CONSTIPATION","TACHYCARDIA","DRY EYE","ACCOMMODATION DISORDER")
HYPER   <- c("ANGIOEDEMA","SWOLLEN TONGUE","SWELLING FACE","DYSPHONIA","HYPERSENSITIVITY","ANAPHYLACTIC REACTION","URTICARIA","RASH","PRURITUS","ERYTHEMA","SWOLLEN LIP","SWELLING","SWOLLEN ARM","STRIDOR")
core <- unique(c(top25, CENTRAL, PERIPH, HYPER, HEAT_PTS))
core <- core[core %in% unique(scop_pt$pt) | scop_pt[, any(pt %in% core)]]  ## 保留存在的
core <- core[core %in% unique(scop_pt$pt)]
cat(sprintf("[核心 PT] %d 个\n", length(core)))

## 途径分层：子集 vs 全库其余（同主分析四格表逻辑）
rt_map <- unique(RT_u[, .(primaryid, route_std)])
calc_g <- function(g, pts) {
  ids_g <- rt_map[route_std==g]$primaryid
  n_g <- length(ids_g)
  sp <- scop_pt[primaryid %in% ids_g & pt %in% pts][, .(a=.N), by=pt]
  rbindlist(lapply(seq_len(nrow(sp)), function(i) {
    p <- sp$pt[i]; a <- sp$a[i]
    c_all <- coalesce2(pt_total[.(p)]$n_total_pt, 0); cc <- c_all - a
    if (a < 3 || cc <= 0) return(NULL)
    b <- n_g - a; d <- N_LIB - n_g - cc
    ror <- (a/b)/(cc/d); se <- sqrt(1/a+1/b+1/cc+1/d)
    data.table(pt=p, route=g, n_route=n_g, a=a, ROR=ror, lcl=exp(log(ror)-1.96*se), ucl=exp(log(ror)+1.96*se))
  }))
}
route_sig <- rbindlist(lapply(c("Transdermal","Oral","Parenteral"), calc_g, pts=core), fill=TRUE)
route_sig[, signal := lcl > 1]
fwrite(route_sig[order(pt, route)], file.path(OUT, "25_途径分层_核心PT.csv"))
cat(sprintf("[写出] 25_途径分层_核心PT.csv（%d 行）\n", nrow(route_sig)))

## 总体口径的核心 PT 表（主分析一致性：any-role vs 全库）
overall_sig <- rbindlist(lapply(core, function(p) {
  a <- sum(scop_pt$pt == p); if (a < 3) return(NULL)
  c_all <- coalesce2(pt_total[.(p)]$n_total_pt, 0); cc <- c_all - a; if (cc<=0) return(NULL)
  b <- N_SCOP - a; d <- N_LIB - N_SCOP - cc
  ror <- (a/b)/(cc/d); se <- sqrt(1/a+1/b+1/cc+1/d)
  data.table(pt=p, a=a, ROR=ror, lcl=exp(log(ror)-1.96*se), ucl=exp(log(ror)+1.96*se))
}))
overall_sig[, signal := lcl > 1]
fwrite(overall_sig[order(-a)], file.path(OUT, "22b_核心PT_总体.csv"))
ps_ids <- unique(drug_rds[drug_group=="scopolamine" & role_cod %in% c("PS","SS")]$primaryid)
scop_pt_ps <- unique(pt_rds[primaryid %in% ps_ids], by=c("primaryid","pt"))
ps_sig <- rbindlist(lapply(core, function(p) {
  a <- sum(scop_pt_ps$pt == p); if (a < 3) return(NULL)
  c_all <- coalesce2(pt_total[.(p)]$n_total_pt, 0); cc <- c_all - a; if (cc<=0) return(NULL)
  b <- length(ps_ids) - a; d <- N_LIB - length(ps_ids) - cc
  ror <- (a/b)/(cc/d); se <- sqrt(1/a+1/b+1/cc+1/d)
  data.table(pt=p, a=a, n_ps_ss=length(ps_ids), ROR=ror, lcl=exp(log(ror)-1.96*se), ucl=exp(log(ror)+1.96*se))
}))
ps_sig[overall_sig, on="pt", ROR_any := i.ROR]
ps_sig[, direction_consistent := sign(log(ROR)-1) == sign(log(ROR_any)-1)]
fwrite(ps_sig, file.path(OUT, "30_敏感性_PSound核心PT.csv"))
cat(sprintf("[写出] 30_敏感性_PSound核心PT.csv（PS/SS n=%d；方向一致 %d/%d）\n",
            length(ps_ids), sum(ps_sig$direction_consistent, na.rm=TRUE), sum(!is.na(ps_sig$direction_consistent))))

## 途径×PT 类别交叉表（v3.1 E3）
cat_map <- function(p) {
  if (p %in% HEAT_PTS) return("Heat-related")
  if (grepl("^APPLICATION SITE|^PRODUCT (ADHESION|QUALITY)|DERMATITIS CONTACT", p)) return("Local/application-site")
  if (p %in% HYPER) return("Hypersensitivity/systemic")
  if (p %in% CENTRAL) return("Central anticholinergic")
  if (p %in% PERIPH) return("Peripheral anticholinergic")
  "Other/general"
}
cross <- scop_pt[pt %in% core][rt_map, on="primaryid", nomatch=0]
cross[, cat := sapply(pt, cat_map)]
xt <- cross[route_std %in% c("Transdermal","Oral","Parenteral"),
            .(n_reports = uniqueN(primaryid)), by=.(route_std, cat)]
xw <- dcast(xt, cat ~ route_std, value.var="n_reports", fill=0)
fwrite(xw, file.path(OUT, "27_途径xPT类别交叉表.csv"))
cat("[写出] 27_途径xPT类别交叉表.csv\n"); print(xw)

## ---------------------------------------------------------------------
## 6. 头对头（v3.1 F 组）：scopolamine vs promethazine vs dimenhydrinate
## ---------------------------------------------------------------------
all_sig <- fread(file.path(OUT, "02_信号检测_全部药物.csv"))
hh <- all_sig[pt %in% core & drug %in% c("scopolamine","promethazine","dimenhydrinate"),
              .(pt, drug, a, ROR)]
hw <- dcast(hh, pt ~ drug, value.var="ROR")
setnames(hw, c("scopolamine","promethazine","dimenhydrinate"), c("ROR_scopolamine","ROR_promethazine","ROR_dimenhydrinate"))
ah <- dcast(hh, pt ~ drug, value.var="a")
setnames(ah, c("scopolamine","promethazine","dimenhydrinate"), c("a_scopolamine","a_promethazine","a_dimenhydrinate"))
hh_t <- merge(hw, ah, by="pt")
hh_t[, dlogROR_prom := log(ROR_scopolamine) - log(ROR_promethazine)]
hh_t[, dlogROR_dmh  := log(ROR_scopolamine) - log(ROR_dimenhydrinate)]
lab <- function(x) ifelse(is.na(x),"NA",ifelse(x > log(1.25),"scopolamine higher",ifelse(x < -log(1.25),"scopolamine lower","comparable")))
hh_t[, verdict_vs_promethazine := lab(dlogROR_prom)]
hh_t[, verdict_vs_dimenhydrinate := lab(dlogROR_dmh)]
fwrite(hh_t[order(-ROR_scopolamine)], file.path(OUT, "26_头对头_tab_headtohead.csv"))
cat("[写出] 26_头对头_tab_headtohead.csv\n")

## ---------------------------------------------------------------------
## 7. 摘要数字核对 facts
## ---------------------------------------------------------------------
ev <- fread(file.path(OUT, "21_两库对比_剔除非事件.csv"))
ms <- fread(file.path(OUT, "04_敏感性分析_晕动症适应证.csv"))
fact <- c(
  sprintf("全库报告（55 季去重后）: %d", N_LIB),
  sprintf("东莨菪碱报告（any role）: %d；PS/SS 口径: %d；角色分布: %s", N_SCOP, length(ps_ids), paste(role_dist$role_cod, role_dist$n_reports, sep="=", collapse="; ")),
  sprintf("阳性信号: 全部 %d → 剔除非事件 PT 后主分析 %d", nrow(sig03), nrow(sig03)-sum(sig03$pt %in% NON_EVENT)),
  sprintf("途径: Transdermal %d / Oral %d / Parenteral %d / Other %d / Unknown %d（可判定 %.1f%%）",
          route_chk[route_std=="Transdermal"]$n_reports, route_chk[route_std=="Oral"]$n_reports,
          route_chk[route_std=="Parenteral"]$n_reports, route_chk[route_std=="Other"]$n_reports,
          route_chk[route_std=="Unknown/Missing"]$n_reports,
          100*(1 - route_chk[route_std=="Unknown/Missing"]$n_reports/N_SCOP)),
  sprintf("晕动症适应证: n=%d (%.1f%%)；亚组阳性信号 %d", N_MS, 100*N_MS/N_SCOP, sum(ms$signal==TRUE)),
  sprintf("严重结局: 死亡 %d (%.2f%%)、住院 %d (%.2f%%)、危及生命 %d、致残 %d",
          sum(ser$death), 100*sum(ser$death)/N_SCOP, sum(ser$hosp), 100*sum(ser$hosp)/N_SCOP, sum(ser$life), sum(ser$disab)),
  sprintf("TTO: 可用 n=%d，median %.0f 天", length(tto_v), median(tto_v)),
  sprintf("外部验证（剔非事件后）: 可评 %d；verdict 计数: %s", sum(ev$faers_signal==TRUE),
          paste(names(table(ev$verdict)), table(ev$verdict), sep="=", collapse="; ")),
  sprintf("heat 类 PT: %s", paste(heat_tab[, paste0(pt, "=", a)], collapse=", ")),
  sprintf("TTO 重拟: fitdistr shape=%.4f (%.4f–%.4f) vs survreg shape=%.4f (%.4f–%.4f)",
          res_all[method=="fitdistr"]$shape, res_all[method=="fitdistr"]$shape_lcl, res_all[method=="fitdistr"]$shape_ucl,
          res_all[method=="survreg"]$shape, res_all[method=="survreg"]$shape_lcl, res_all[method=="survreg"]$shape_ucl)
)
writeLines(fact, file.path(OUT, "29_摘要数字核对_facts.txt"))
cat(paste(fact, collapse="\n"))
cat("\n[完成] 07_w1prime_tables.R\n")
