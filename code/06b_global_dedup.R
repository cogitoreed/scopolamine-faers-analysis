## =====================================================================
## 06b_global_dedup.R —— 全库跨季度去重（FDA 官方算法），产出存活 primaryid 集
##   背景：faers_scopolamine.R 原实现在季度内去重，跨季度 caseid 更新版未去重
##   （实测：东莨菪碱 9,548 行仅 7,919 唯一 caseid）。本脚本读全部 55 季 DEMO，
##   按 caseid 取 caseversion 最大者（平手取 fda_dt 最新），输出存活 primaryid。
##   同时修复 DEMO 列名（新旧两代表头按文件自身 header 读，header=TRUE）。
## 运行：LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 R_LIBS_USER=~/Rlibs Rscript 06b_global_dedup.R
## =====================================================================
suppressMessages(library(data.table))
QDIR <- "data/extracted"
OUT  <- "output_v2"
dir.create(OUT, showWarnings = FALSE)

demo_files <- sort(list.files(QDIR, pattern="^demo.*\\.txt$", ignore.case=TRUE,
                              recursive=TRUE, full.names=TRUE))
stopifnot(length(demo_files) == 55)

OLD_COLS <- c("primaryid","caseid","caseversion","i_f_code","event_dt","mfr_dt",
              "init_fda_dt","fda_dt","rept_cod","mfr_num","mfr_sndr","age","age_cod",
              "gndr_cod","e_sub","wt","wt_cod","rept_dt","to_mfr","occp_cod",
              "reporter_country","occr_country")
NEW_COLS <- c("primaryid","caseid","caseversion","i_f_code","event_dt","mfr_dt",
              "init_fda_dt","fda_dt","rept_cod","auth_num","mfr_num","mfr_sndr",
              "lit_ref","age","age_cod","age_grp","sex","e_sub","wt","wt_cod",
              "rept_dt","to_mfr","occp_cod","reporter_country","occr_country")

acc <- list()
raw_rows <- 0
t0 <- Sys.time()
for (i in seq_along(demo_files)) {
  f <- demo_files[i]
  hdr <- strsplit(toupper(sub("\\s+", "", readLines(f, n=1, warn=FALSE))), "$", fixed=TRUE)[[1]]
  is_old <- "GNDR_COD" %in% hdr
  want <- c("primaryid","caseid","caseversion","event_dt","fda_dt")
  dt <- fread(f, sep="$", quote="", header=TRUE, select=intersect(want, hdr),
              colClasses="character", na.strings=c(""), showProgress=FALSE)
  ## 表头行若被 fread 重复读入则按 primaryid=="primaryid" 过滤（header=TRUE 时不会，保险起见）
  dt <- dt[primaryid != "primaryid"]
  dt[, caseversion_n := suppressWarnings(as.numeric(caseversion))]
  dt[, fda_n := suppressWarnings(as.numeric(fda_dt))]
  dt[, ev_n := suppressWarnings(as.numeric(event_dt))]
  acc[[i]] <- dt[, .(primaryid, caseid, caseversion_n, ev_n, fda_n)]
  raw_rows <- raw_rows + nrow(dt)
  if (i %% 10 == 0 || i == length(demo_files))
    message(sprintf("[%2d/55] 累计原始行 %d，耗时 %.1f 分钟", i, raw_rows, as.numeric(difftime(Sys.time(), t0, units="mins"))))
}
all_demo <- rbindlist(acc)
rm(acc); invisible(gc())

n_raw <- nrow(all_demo)
n_caseid <- uniqueN(all_demo$caseid)
cat(sprintf("[合并] 原始 DEMO 行 %d；唯一 caseid %d（跨季度重复 %.1f%%）\n",
            n_raw, n_caseid, 100*(1 - n_caseid/n_raw)))

## FDA 算法：caseid 取 caseversion 最大；仍重复取 fda_dt 最新（再平手取 event_dt 最新）
setorder(all_demo, caseid, -caseversion_n, -fda_n, -ev_n)
surv <- unique(all_demo, by="caseid")
n_surv <- nrow(surv)
cat(sprintf("[去重后] 存活报告 %d（剔除 %d，%.1f%%）\n", n_surv, n_raw - n_surv, 100*(n_raw-n_surv)/n_raw))

saveRDS(surv$primaryid, file.path(OUT, "global_surviving_primaryids.rds"))
fwrite(data.table(
  metric = c("raw_demo_rows","unique_caseids","post_dedup_reports","removed_rows","pct_removed",
             "quarters","generated_at"),
  value  = c(n_raw, n_caseid, n_surv, n_raw-n_surv, round(100*(n_raw-n_surv)/n_raw,2),
             55, format(Sys.time())),
  stringsAsFactors = FALSE), file.path(OUT, "00_全局去重日志.csv"))
cat("[完成] 06b_global_dedup.R → output_v2/global_surviving_primaryids.rds\n")
