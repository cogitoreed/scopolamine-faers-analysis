###############################################################################
## FAERS / AEMS 抗晕动病药物不良事件信号挖掘
##
## 拟定题目：
##   Adverse event profiles of scopolamine and other motion sickness
##   prophylactics: a disproportionality analysis of the FDA Adverse Event
##   Reporting System
##
## 编制日期：2026-09-02
## 环境：R >= 4.2，需联网下载数据，建议内存 >= 8 GB
##
## 使用方法：见同目录《操作手册_FAERS从零到投稿.md》
##   1. 先设置下方 USER CONFIG 区的路径与季度范围
##   2. 全脚本运行一次到底，结果写入 output/
##
## 设计要点：
##   * 内存优化：不把 2000 万条全库明细读进内存。只累积"全库各 PT 的报告数"
##     作为信号检测的对照分母，目标药物的报告明细才落地保存。
##   * 兼容新旧格式：FAERS 在 2015Q1 前后 DRUG 表的 drugname / role_cod 列顺序
##     相反，DEMO 表 2014 年前缺少 occp_cod，脚本按位置自动识别。
##   * 适应证分层：用 INDI 表筛出"用于晕动症"的报告做敏感性分析，解决
##     东莨菪碱适应证混杂（术后恶心呕吐、临终镇静、胃肠解痉）的问题。
###############################################################################


## ============================================================================
## 0. USER CONFIG —— 只需要改这一段
## ============================================================================

WORK_DIR   <- "/Users/xiaocaixu/WorkBuddy/学位自救计划/FAERS分析"
DATA_DIR   <- file.path(WORK_DIR, "data")            # 存放下载的季度 zip
EXTRACT_DIR<- file.path(WORK_DIR, "data", "extracted") # 解压目录（可删）
OUT_DIR    <- file.path(WORK_DIR, "output_v2")

## 季度范围。
## 注意：FDA 官方目前只开放 2012Q4–2026Q2 共 55 个季度（更早的返回 HTTP 500），
##      故研究窗口为 2013–2026 年，约 3.5 GB。这已远超 FAERS 类论文的常规年限要求。
## 首次跑建议先用 TEST_MODE <- TRUE 只跑最近 2 年验证流程，再跑全量。
TEST_MODE    <- FALSE
YEARS_FULL   <- 2013:2026   # 2012 仅取 Q4，由 make_quarters 自动补入
YEARS_TEST   <- 2025:2026

## 目标药物与对照药物（大小写不敏感，正则）
## 注意：hyoscine butylbromide / buscopan 为丁溴东莨菪碱，主要作胃肠解痉、
##       中枢抗胆碱作用弱，单独归入 EXCLUDE 列表，不混入东莨菪碱主分析。
DRUG_LIST <- list(
  scopolamine   = "SCOPOLAMINE|HYOSCINE|TRANSDERM[ -]?SCOP|TRANSDERM[ -]?V|MALDEMAR|ISOPto hyoscine",
  dimenhydrinate= "DIMENHYDRINATE|DRAMAMINE|GRAVOL|GRAVOL|NAUSICALM",
  meclizine     = "MECLIZINE|MECLOZINE|ANTIVERT|BONINE|POSTAFENE",
  cinnarizine   = "CINNARIZINE|CINNARAZINE|STUGERON|STUNARONE",
  promethazine  = "PROMETHAZINE|PHENERGAN"
)

## 排除丁溴东莨菪碱等季铵类（不透过 BBB，ADE 谱完全不同）
## 实测 2026Q1：METHSCOPOLAMINE BROMIDE 存在，必须一并排除（季铵，同 Buscopan 道理）
EXCLUDE_PATTERN <- "BUTYLBROMIDE|BUTYLSCOPOLAMINE|BUSCOPAN|N-?BUTYL|METHSCOPOLAMINE"

## 晕动症相关适应证（用于敏感性分析，MedDRA PT 术语）
MS_INDI_PATTERN <- "MOTION SICKNESS|MOTION SICKNESS PROPHYLAXIS|TRAVEL SICKNESS|SEASICKNESS|KINETOSIS"

## 重点关注的系统器官分类（SOC 关键词，用于结果聚焦展示）
FOCUS_PT_PATTERN <- paste0(
  "SOMNOLENCE|DROWSINESS|SEDATION|LETHARGY|CONFUSIONAL STATE|DISORIENTATION|",
  "HALLUCINATION|DELIRIUM|AGITATION|AMNESIA|DIZZINESS|VERTIGO|HEADACHE|",
  "VISION BLURRED|MYDRIASIS|GLAUCOMA|ACCOMMODATION DISORDER|DRY EYE|DRY MOUTH|",
  "URINARY RETENTION|DYSURIA|TACHYCARDIA|BRADYCARDIA|HYPOTENSION|",
  "DERMATITIS CONTACT|ERYTHEMA|PRURITUS|RASH|APPLICATION SITE",
  collapse = ""
)

## 信号判定标准
MIN_CASES   <- 3      # 最少报告数
PRR_THRESH  <- 2      # PRR 阈值
CHISQ_THRESH<- 4      # PRR 卡方阈值


## ============================================================================
## 1. 环境准备
## ============================================================================

dir.create(DATA_DIR,    recursive = TRUE, showWarnings = FALSE)
dir.create(EXTRACT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR,     recursive = TRUE, showWarnings = FALSE)

pkgs <- c("data.table", "stringr", "ggplot2", "survival", "writexl")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) {
  message("正在安装缺失的包：", paste(to_install, collapse = ", "))
  install.packages(to_install, repos = "https://cloud.r-project.org")
}
suppressMessages({
  library(data.table); library(stringr); library(ggplot2); library(survival)
})
if (requireNamespace("writexl", quietly = TRUE)) library(writexl)

setDTthreads(0)  # 使用全部 CPU 核心

YEARS <- if (TEST_MODE) YEARS_TEST else YEARS_FULL
message("[CONFIG] 测试模式 = ", TEST_MODE, " | 年份范围 = ",
        min(YEARS), "-", max(YEARS))


## ============================================================================
## 2. 工具函数
## ============================================================================

## 生成季度标签（四位年份，如 "2013Q1"，与下载脚本的文件名一致），
## 自动补入 2012Q4（FDA 开放的最早季度）、截到 2026Q2（目前最新）
make_quarters <- function(years) {
  qs <- c()
  for (y in years) {
    for (q in 1:4) {
      if (y == 2026 && q > 2) next
      qs <- c(qs, paste0(y, "Q", q))
    }
  }
  if (min(years) <= 2013) qs <- c("2012Q4", qs)
  sort(unique(qs))
}

## 读取 FAERS ASCII：$ 分隔、无表头、Latin-1 编码、可能有引号
read_faers_txt <- function(path) {
  if (is.null(path) || length(path) == 0 || !file.exists(path)) return(NULL)
  dt <- tryCatch(
    fread(path, sep = "$", header = FALSE, quote = "\"",
          encoding = "Latin-1", fill = TRUE, na.strings = c("", "NA"),
          showProgress = FALSE, colClasses = "character"),
    error = function(e) { warning("读取失败：", basename(path), " -> ", e$message); NULL }
  )
  dt
}

## 补齐列名（按位置，取前 n 个已知列）
name_cols <- function(dt, cols) {
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  n <- min(length(cols), ncol(dt))
  setnames(dt, 1:n, cols[1:n])
  if (ncol(dt) > n) dt[, (n + 1):ncol(dt) := NULL]
  dt
}

## 统一日期：8位 yyyymmdd / 6位 yyyymm / 4位 yyyy 都转成 Date
## 注意：DEMO 的 event_dt 含 NA/空串，必须先过滤再取 nchar，
##       否则 NA 逻辑下标触发 "NAs are not allowed in subscripted assignments"（已实测踩坑）
parse_fda_date <- function(x) {
  x   <- str_trim(as.character(x))
  res <- as.Date(rep(NA_character_, length(x)))
  ok  <- !is.na(x) & nzchar(x)
  if (any(ok)) {
    v  <- x[ok]
    nc <- nchar(v)
    v[nc == 6] <- paste0(v[nc == 6], "01")
    v[nc == 4] <- paste0(v[nc == 4], "0101")
    d <- as.Date(v, format = "%Y%m%d")
    d[!is.na(d) & (d < as.Date("1900-01-01") | d > Sys.Date())] <- NA
    res[ok] <- d
  }
  res
}

## 年龄统一折算为"岁"（向量化版本：switch 不接受向量，
## 而 data.table 按整列调用，必须用逻辑索引，已实测踩坑）
age_to_year <- function(age, cod) {
  age  <- suppressWarnings(as.numeric(age))
  f    <- rep(1, length(age))
  c_u  <- toupper(cod)
  f[!is.na(c_u) & c_u == "MON"] <- 1/12
  f[!is.na(c_u) & c_u == "WK"]  <- 1/52.25
  f[!is.na(c_u) & c_u == "DY"]  <- 1/365.25
  f[!is.na(c_u) & c_u == "HR"]  <- 1/8766
  age * f
}


## ============================================================================
## 3. 主循环：逐季度读取、去重、汇总
## ============================================================================
##
## 每个季度做三件事：
##   a) 汇总 REAC 表 -> 累加到全库 PT 计数（作为信号检测的对照分母 c+d 侧）
##   b) 记录该季度去重后的总报告数
##   c) 从 DRUG 表找出目标药物报告，落地保存其 DEMO / REAC / OUTC / THER / INDI
##
## ============================================================================

quarters <- make_quarters(YEARS)
message("[STEP 3] 共 ", length(quarters), " 个季度待处理")

## 累积容器
total_reports   <- 0L                      # 全库去重后总报告数
pt_total_count  <- list()                  # list[pt] = 全库该 PT 的报告数
drug_demo <- list(); drug_reac <- list()
drug_outc <- list(); drug_ther <- list(); drug_indi <- list()
drug_drug <- list()
quarter_log <- list()

DEMO_COLS <- c("primaryid","caseid","caseversion","i_f_code","event_dt","mfr_dt",
               "init_fda_dt","fda_dt","rept_cod","auth_num","age","age_cod",
               "age_grp","sex","wt","wt_cod","reporter_country","occr_country",
               "occp_cod")
REAC_COLS <- c("primaryid","caseid","pt","drug_rec_act")
OUTC_COLS <- c("primaryid","caseid","outc_cod")
THER_COLS <- c("primaryid","caseid","dsg_drug_seq","start_dt","end_dt","dur","dur_cod")
## INDI 实测 2026Q1：仅 4 列（primaryid$caseid$indi_drug_seq$indi_pt），
## 新版格式已无 drug_seq 列；按 5 列命名会错位导致 indi_pt 不存在（已实测踩坑）
INDI_COLS <- c("primaryid","caseid","indi_drug_seq","indi_pt")
RPSR_COLS <- c("primaryid","caseid","rpsr_cod")

## ---- v2 修复（2026-09-04）------------------------------------------------
## ① DEMO 新旧两代表头不同（2014Q3 起增 auth_num/lit_ref/age_grp，gndr_cod→sex）。
##    原 header=FALSE + 固定列名按位置命名 → 错位（age 列读到药厂名）+ 表头行被当数据。
##    改为按文件自身 header=TRUE 读、按列名取列。
read_demo_era <- function(path) {
  if (is.null(path) || length(path) == 0 || !file.exists(path)) return(NULL)
  hdr <- strsplit(sub("\\s+", "", readLines(path, n = 1, warn = FALSE)),
                  "$", fixed = TRUE)[[1]]
  hdr_u <- toupper(hdr)
  want <- c("primaryid","caseid","caseversion","event_dt","age","age_cod",
            if ("GNDR_COD" %in% hdr_u) "gndr_cod" else "sex",
            "occp_cod","reporter_country","occr_country")
  sel <- intersect(want, tolower(hdr))
  dt <- fread(path, sep = "$", quote = "", header = TRUE, select = sel,
              colClasses = "character", na.strings = c(""), showProgress = FALSE)
  if ("gndr_cod" %in% names(dt)) setnames(dt, "gndr_cod", "sex")
  dt[]
}
## ② 跨季度全局去重（FDA 官方算法）：06b_global_dedup.R 预先产出存活 primaryid 集。
##    原实现只在季度内去重，跨季度 caseid 更新版未去重（实测东莨菪碱 9,548 行仅 7,919 唯一 caseid）。
GLOBAL_IDS_FILE <- file.path(WORK_DIR, "output_v2", "global_surviving_primaryids.rds")
GLOBAL_IDS <- if (file.exists(GLOBAL_IDS_FILE)) readRDS(GLOBAL_IDS_FILE) else NULL
if (is.null(GLOBAL_IDS)) stop("缺少 output_v2/global_surviving_primaryids.rds —— 请先运行 06b_global_dedup.R")

for (i in seq_along(quarters)) {
  q   <- quarters[i]
  ## 下载脚本统一存为小写文件名，这里大小写都找一遍
  zip <- file.path(DATA_DIR, paste0("faers_ascii_", tolower(q), ".zip"))
  if (!file.exists(zip)) zip <- file.path(DATA_DIR, paste0("faers_ascii_", q, ".zip"))
  if (!file.exists(zip)) {
    warning("缺少数据文件，已跳过：faers_ascii_", q, ".zip"); next
  }

  ## --- 解压到独立目录 -------------------------------------------------------
  exdir <- file.path(EXTRACT_DIR, q)
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
  if (length(list.files(exdir)) == 0) {
    utils::unzip(zip, exdir = exdir, junkpaths = TRUE)
  }
  f <- function(pfx) {
    hits <- list.files(exdir, pattern = paste0("^", pfx, ".*\\.TXT$"),
                       ignore.case = TRUE, full.names = TRUE, recursive = TRUE)
    if (length(hits) == 0) NULL else hits[1]
  }

  ## --- 读 DEMO（v2：按各季表头读 + 全局去重过滤），再季度内去重（此时已近无重复）--
  demo <- read_demo_era(f("DEMO"))
  if (is.null(demo)) { warning(q, " DEMO 读取失败，跳过"); next }
  demo <- demo[primaryid %in% GLOBAL_IDS]
  demo[, caseversion := suppressWarnings(as.numeric(caseversion))]

  ## 去重规则 1：同一 caseid 保留版本最高的记录
  setorder(demo, caseid, -caseversion)
  demo <- unique(demo, by = "caseid")

  ## 去重规则 2：caseid + 事件日 + 年龄 + 性别 + 报告国 相同的判为重复报告
  demo[, event_dt_c := as.character(event_dt)]
  demo <- demo[!duplicated(demo, by = c("caseid","event_dt_c","age","sex","reporter_country"))]
  demo[, event_dt_c := NULL]

  n_q <- nrow(demo)
  total_reports <- total_reports + n_q
  demo[, event_dt := parse_fda_date(event_dt)]
  demo[, age_yr   := age_to_year(age, age_cod)]
  demo[, quarter  := q]

  ids <- demo$primaryid

  ## --- 读 REAC：汇总 PT 计数 + 落地目标药部分 -------------------------------
  reac <- name_cols(read_faers_txt(f("REAC")), REAC_COLS)
  if (!is.null(reac)) {
    reac_q <- reac[primaryid %in% ids, .(primaryid, pt = toupper(str_trim(pt)))]
    reac_q <- reac_q[!is.na(pt) & pt != ""]
    ## 累积全库 PT 计数（按 primaryid 去重，一个报告对同一 PT 只计一次）
    cnt <- unique(reac_q, by = c("primaryid","pt"))[, .N, by = pt]
    for (k in seq_len(nrow(cnt))) {
      p <- cnt$pt[k]; v <- cnt$N[k]
      pt_total_count[[p]] <- (if (is.null(pt_total_count[[p]])) 0 else pt_total_count[[p]]) + v
    }
  } else reac_q <- data.table(primaryid = character(), pt = character())

  ## --- 读 DRUG：识别目标药物 ------------------------------------------------
  drug <- name_cols(read_faers_txt(f("DRUG")), c("primaryid","caseid","drug_seq",
                                                 "c4","c5","prod_ai","val_vbm","route"))
  target_ids <- character(0)
  drug_hit   <- data.table()
  if (!is.null(drug) && nrow(drug) > 0) {
    ## 兼容新旧格式：判断第4列是否为 role_cod
    sample_vals <- toupper(na.omit(unique(head(drug$c4, 5000))))
    is_new <- any(sample_vals %in% c("PS","SS","C","I"))
    if (is_new) { setnames(drug, "c4", "role_cod"); setnames(drug, "c5", "drugname") }
    else        { setnames(drug, "c4", "drugname"); setnames(drug, "c5", "role_cod") }

    drug[, drugname_u := toupper(str_trim(drugname))]
    drug[, is_target := FALSE]

    for (nm in names(DRUG_LIST)) {
      pat  <- DRUG_LIST[[nm]]
      hit  <- str_detect(drug$drugname_u, regex(pat, ignore_case = TRUE))
      ## 排除丁溴东莨菪碱等
      hit  <- hit & !str_detect(drug$drugname_u, regex(EXCLUDE_PATTERN, ignore_case = TRUE))
      drug[hit, `:=`(is_target = TRUE, drug_group = nm)]
    }

    drug_hit <- drug[is_target == TRUE & primaryid %in% ids,
                     .(primaryid, drug_seq, drugname_u, role_cod = toupper(role_cod),
                       route = toupper(str_trim(route)), drug_group)]
    if (nrow(drug_hit)) target_ids <- unique(drug_hit$primaryid)
  }

  ## --- 落地目标药报告的其余表 ----------------------------------------------
  if (length(target_ids)) {
    drug_demo[[q]] <- demo[primaryid %in% target_ids]
    drug_reac[[q]] <- reac_q[primaryid %in% target_ids]
    drug_drug[[q]] <- unique(drug_hit, by = c("primaryid","drug_seq"))

    outc <- name_cols(read_faers_txt(f("OUTC")), OUTC_COLS)
    if (!is.null(outc))
      drug_outc[[q]] <- outc[primaryid %in% target_ids,
                             .(primaryid, outc_cod = toupper(str_trim(outc_cod)))]

    ther <- name_cols(read_faers_txt(f("THER")), THER_COLS)
    if (!is.null(ther))
      drug_ther[[q]] <- ther[primaryid %in% target_ids,
                             .(primaryid, dsg_drug_seq, start_dt, end_dt)]

    indi <- name_cols(read_faers_txt(f("INDI")), INDI_COLS)
    if (!is.null(indi))
      drug_indi[[q]] <- indi[primaryid %in% target_ids,
                             .(primaryid, indi_pt = toupper(str_trim(indi_pt)))]
  }

  quarter_log[[q]] <- data.table(quarter = q, n_reports = n_q,
                                 n_target_reports = length(target_ids))
  message(sprintf("  [%2d/%d] %s  报告 %d  目标药报告 %d",
                  i, length(quarters), q, n_q, length(target_ids)))

  rm(demo, reac, drug); invisible(gc(verbose = FALSE))
}

## --- 合并结果 ---------------------------------------------------------------
DEMO_T <- rbindlist(drug_demo, fill = TRUE)
REAC_T <- rbindlist(drug_reac, fill = TRUE)
DRUG_T <- rbindlist(drug_drug, fill = TRUE)
OUTC_T <- rbindlist(drug_outc, fill = TRUE)
THER_T <- rbindlist(drug_ther, fill = TRUE)
INDI_T <- rbindlist(drug_indi, fill = TRUE)
QLOG   <- rbindlist(quarter_log, fill = TRUE)

PT_TOTAL <- data.table(pt = names(pt_total_count),
                       n_total_pt = unlist(pt_total_count, use.names = FALSE))

message("[STEP 3 完成] 全库去重后总报告数 = ", format(total_reports, big.mark = ","))
message("             目标药物报告数 = ", format(uniqueN(DEMO_T$primaryid), big.mark = ","))
message("             PT 总条目数 = ", format(nrow(PT_TOTAL), big.mark = ","))

fwrite(QLOG,     file.path(OUT_DIR, "00_季度处理日志.csv"))
fwrite(PT_TOTAL, file.path(OUT_DIR, "00_全库PT计数.csv"))
saveRDS(list(DEMO_T, REAC_T, DRUG_T, OUTC_T, THER_T, INDI_T,
             PT_TOTAL, total_reports),
        file.path(OUT_DIR, "faers_target.rds"))


## ============================================================================
## 4. 构建分析数据集
## ============================================================================

## 一个报告可能同时含多种目标药（如东莨菪碱+茶苯海明复方）。
## 做法：为每个"药物-报告"对建立记录，主分析按药物分组独立计算。
REP <- unique(DRUG_T, by = c("primaryid","drug_group"))[, .(primaryid, drug_group)]

## 报告层面信息
DEMO_T[, age_grp := cut(age_yr, breaks = c(-Inf, 17, 64, Inf),
                        labels = c("<18", "18-64", ">=65"))]
DEMO_T[, sex_f := factor(toupper(sex), levels = c("F","M","UNK","NS"))]

## 每个报告是否用于晕动症适应证
ms_indi <- unique(INDI_T[str_detect(indi_pt, regex(MS_INDI_PATTERN, ignore_case = TRUE)),
                         .(primaryid, ms_indi = TRUE)])

## 每个报告的严重结局
## 严重结局标志（每报告一行）
## 注意：data.table 的 j 混合"列引用 + 聚合标量"时，any() 会在全表求值
## 并按行回收（不分组）→ 所有报告的 death/hosp 全为 TRUE 的静默错误！
## 必须用 by=primaryid 分组聚合（已实测踩坑：真 DE=18 却算出 144）。
## 修复前（错误）：unique(OUTC_T[, .(primaryid, death=any(...)), by="primaryid"])
OUTC_SER <- OUTC_T[, .(death = any(outc_cod == "DE"),
                       hosp  = any(outc_cod == "HO"),
                       disab = any(outc_cod == "DS"),
                       life  = any(outc_cod == "LT")), by = primaryid]

## 每个报告的给药途径（取该药物对应的 route）+ 标准化分类
## 实测 2026Q1：途径缺失/UNKNOWN 约 78%；同义词必须合并，否则分层无效
std_route <- function(r) {
  r <- toupper(str_trim(as.character(r)))
  out <- rep("Unknown/Missing", length(r))
  out[grepl("TRANSDERM|PATCH", r)] <- "Transdermal"
  out[grepl("TOPICAL|CUTANEOUS", r) & out == "Unknown/Missing"] <- "Transdermal"
  out[grepl("ORAL|MOUTH|SUBLINGUAL|^PO$", r) & out == "Unknown/Missing"] <- "Oral"
  out[grepl("INTRAVEN|INTRAMUSC|SUBCUT|INJECT|INFUSION|^IV$|^IM$|^SC$", r) &
      out == "Unknown/Missing"] <- "Parenteral"
  out[r != "" & r != "UNKNOWN" & out == "Unknown/Missing"] <- "Other"
  out
}
RT <- unique(DRUG_T[, .(primaryid, drug_group, route)], by = c("primaryid","drug_group"))
RT[, route_std := std_route(route)]

## 途径分布汇总（论文 Table 1 素材；缺失率必须写进 Limitations，
## 并做"仅纳入可判定途径者"的敏感性分析——单季实测缺失约 78%）
route_summary <- RT[, .(n_reports = uniqueN(primaryid)), by = .(drug_group, route_std)][
  order(drug_group, -n_reports)]
fwrite(route_summary, file.path(OUT_DIR, "01b_给药途径分布.csv"))

## 报告-PT 明细
RD <- merge(REP, DEMO_T, by = "primaryid", all.x = TRUE, allow.cartesian = FALSE)
RD <- merge(RD, unique(REAC_T, by = c("primaryid","pt")), by = "primaryid", allow.cartesian = TRUE)
RD <- merge(RD, RT,  by = c("primaryid","drug_group"), all.x = TRUE)
RD <- merge(RD, ms_indi,  by = "primaryid", all.x = TRUE)
RD <- merge(RD, OUTC_SER, by = "primaryid", all.x = TRUE)
RD[is.na(ms_indi), ms_indi := FALSE]
RD[, is_MS := ms_indi]

fwrite(RD, file.path(OUT_DIR, "01_报告-PT明细.csv"))


## ============================================================================
## 5. 不成比例性分析（信号检测）
## ============================================================================
##
## 四格表（以药物 D 与不良事件 E 为例）：
##                  E       非E
##   药物 D          a        b
##   其他药物        c        d
##
##   c = 全库 E 的报告数 - a
##   d = (全库总报告数 - 药物 D 的报告数) - c
##
## ============================================================================

calc_signals <- function(dat, pt_total, n_total, drug_label) {
  ## dat: data.table(primaryid, pt)，该药物的报告-PT 明细
  n_drug <- uniqueN(dat$primaryid)
  cnt    <- dat[, .(a = uniqueN(primaryid)), by = pt]
  if (nrow(cnt) == 0) return(NULL)

  cnt <- merge(cnt, pt_total, by = "pt", all.x = TRUE)
  cnt[is.na(n_total_pt), n_total_pt := 0]
  cnt[, c := n_total_pt - a]
  cnt[c < 0, c := 0]

  b <- n_drug - cnt$a
  d <- (n_total - n_drug) - cnt$c
  d[d < 0] <- 0

  res <- cnt[, .(pt, a, b = n_drug - a, c, d = pmax(0, (n_total - n_drug) - c))]

  ## ROR（乘法必须转 double：千万级 a*d 会整数溢出 → NA，已实测踩坑）
  res[, ROR := (as.numeric(a) * as.numeric(d)) / (as.numeric(b) * as.numeric(c))]
  res[, se_ln_ror := sqrt(1/a + 1/b + 1/c + 1/d)]
  res[, ROR_lcl := exp(log(ROR) - 1.96 * se_ln_ror)]
  res[, ROR_ucl := exp(log(ROR) + 1.96 * se_ln_ror)]

  ## PRR
  res[, PRR := (a/(a+b)) / (c/(c+d))]
  res[, se_ln_prr := sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d))]
  res[, PRR_lcl := exp(log(PRR) - 1.96 * se_ln_prr)]
  res[, PRR_ucl := exp(log(PRR) + 1.96 * se_ln_prr)]
  ## 注意：全库千万级报告时 a*d 等乘积会超出 R 整数上限（int32 ≈ 21 亿），
  ## 结果为 NA → 所有卡方信号判定静默失效。必须先转 double。（已实测踩坑）
  res[, chisq := ((as.numeric(a)*as.numeric(d) - as.numeric(b)*as.numeric(c))^2) *
        as.numeric(a+b+c+d) /
        (as.numeric(a+b)*as.numeric(c+d)*as.numeric(a+c)*as.numeric(b+d))]

  ## BCPNN (IC / IC025)
  res[, n := a + b + c + d]
  res[, IC := log2(as.numeric(a) * as.numeric(n) /
        (as.numeric(a+b) * as.numeric(a+c)))]
  res[, IC025 := IC - 3.3 * n^(-0.5) - 2 * n^(-1.5)]

  ## 信号判定
  ## 信号判定策略（应对"论文工厂"污名化）：
  ##   主分析只用 ROR（Evans 2001 标准：a>=3 且 ROR 95%CI 下限 > 1）。
  ##   PRR + BCPNN 多法齐用已被点名为"论文工厂"识别特征（Khouri 2025：
  ##   "多种统计方法做同一分析，实属冗余"）。三法一致降级为敏感性分析。
  res[, sig_ROR   := a >= MIN_CASES & ROR_lcl > 1]
  res[, sig_PRR   := a >= MIN_CASES & PRR >= PRR_THRESH & chisq >= CHISQ_THRESH]
  res[, sig_BCPNN := a >= MIN_CASES & IC025 > 0]
  res[, signal := sig_ROR]                                # 主分析信号
  res[, sens_trio := sig_ROR & sig_PRR & sig_BCPNN]       # 敏感性：三法一致
  res[, drug := drug_label]

  setorder(res, -a)
  res[]
}

## 对每种药物分别计算
signal_all <- rbindlist(lapply(names(DRUG_LIST), function(nm) {
  dat <- RD[drug_group == nm, .(primaryid, pt)]
  message("  [信号检测] ", nm, " n = ", uniqueN(dat$primaryid))
  calc_signals(dat, PT_TOTAL, total_reports, nm)
}), fill = TRUE)

fwrite(signal_all, file.path(OUT_DIR, "02_信号检测_全部药物.csv"))

## 仅东莨菪碱、仅阳性信号
sig_scop <- signal_all[drug == "scopolamine" & signal == TRUE]
fwrite(sig_scop, file.path(OUT_DIR, "03_东莨菪碱_阳性信号.csv"))

## 敏感性分析：ROR+PRR+BCPNN 三法一致的信号（应对审稿质疑稳健性）
sig_trio <- signal_all[signal == TRUE & sens_trio == TRUE]
fwrite(sig_trio, file.path(OUT_DIR, "03b_东莨菪碱_三法一致敏感性.csv"))

## 适应证敏感性分析：只保留用于晕动症的报告
sens_all <- rbindlist(lapply(names(DRUG_LIST), function(nm) {
  dat <- RD[drug_group == nm & is_MS == TRUE, .(primaryid, pt)]
  if (nrow(dat) == 0) return(NULL)
  calc_signals(dat, PT_TOTAL, total_reports, paste0(nm, " [晕动症适应证]"))
}), fill = TRUE)
if (nrow(sens_all)) fwrite(sens_all, file.path(OUT_DIR, "04_敏感性分析_晕动症适应证.csv"))


## ============================================================================
## 6. 时间-发生分析（Weibull）
## ============================================================================

tto <- merge(
  THER_T[, .(primaryid, start_dt = parse_fda_date(start_dt))],
  unique(DEMO_T[, .(primaryid, event_dt)]), by = "primaryid", all.x = TRUE
)
tto <- tto[!is.na(start_dt) & !is.na(event_dt)]
tto[, tto_days := as.numeric(event_dt - start_dt)]
## 必须用 > 0 而非 >= 0：Weibull 密度在 0 处发散 → NaN → 拟合失败（已实测踩坑）
tto <- tto[tto_days > 0 & tto_days <= 3650]

weibull_res <- data.table()
for (nm in names(DRUG_LIST)) {
  ids <- RD[drug_group == nm, unique(primaryid)]
  sub <- unique(tto[primaryid %in% ids, .(primaryid, tto_days)], by = "primaryid")
  if (nrow(sub) < 20) {
    message("  [Weibull] ", nm, " 样本不足（n=", nrow(sub), "），跳过"); next
  }
  fit <- survreg(Surv(tto_days, rep(1, nrow(sub))) ~ 1, data = sub, dist = "weibull")
  sh  <- 1 / fit$scale          # shape
  sc  <- exp(coef(fit))         # scale
  se  <- sqrt(diag(vcov(fit)))[2] / (fit$scale^2)
  weibull_res <- rbind(weibull_res, data.table(
    drug = nm, n = nrow(sub),
    median_tto = median(sub$tto_days),
    shape = sh, shape_lcl = sh - 1.96 * se, shape_ucl = sh + 1.96 * se,
    scale = sc,
    pattern = ifelse(sh + 1.96*se < 1, "早期衰减型",
              ifelse(sh - 1.96*se > 1, "蓄积型", "随机型"))
  ))
}
if (nrow(weibull_res)) fwrite(weibull_res, file.path(OUT_DIR, "05_Weibull时间发生分析.csv"))


## ============================================================================
## 7. 亚组分析与严重结局
## ============================================================================

## 7.1 亚组：按药物 × 重点不良事件，计算各亚组的 ROR
focus_pts <- PT_TOTAL[str_detect(pt, regex(FOCUS_PT_PATTERN, ignore_case = TRUE)), pt]
focus_pts <- intersect(focus_pts, signal_all[a >= MIN_CASES, unique(pt)])

subgroup <- function(by_var, label) {
  out <- list()
  for (nm in names(DRUG_LIST)) {
    vals <- RD[[by_var]]
    for (lv in unique(vals[!is.na(vals)])) {   ## na.unique 不存在，改为显式去 NA（原代码 bug）
      if (is.na(lv)) next
      sel <- RD[drug_group == nm & get(by_var) == lv & pt %in% focus_pts]
      if (nrow(sel) == 0) next
      r <- calc_signals(sel[, .(primaryid, pt)], PT_TOTAL, total_reports, nm)
      if (is.null(r)) next
      out[[paste(nm, lv)]] <- r[, .(drug, pt, a, ROR, ROR_lcl, ROR_ucl, signal)][
        , `:=`(subgroup_var = label, subgroup_level = as.character(lv))]
    }
  }
  rbindlist(out, fill = TRUE)
}

sg_sex  <- subgroup("sex_f",      "性别")
sg_age  <- subgroup("age_grp",    "年龄组")
sg_ms   <- subgroup("is_MS",      "晕动症适应证")
sg_route<- subgroup("route_std",  "给药途径")

subgroup_all <- rbindlist(list(sg_sex, sg_age, sg_ms, sg_route), fill = TRUE)
if (nrow(subgroup_all)) fwrite(subgroup_all, file.path(OUT_DIR, "06_亚组分析.csv"))

## 7.2 严重结局构成
## 注意：RD 是"报告×PT"长表，同一报告会重复 PT 次！必须先按报告去重，
## 否则 n_death > n_reports 的荒谬结果（已实测踩坑：单季 death 1001 vs 报告 203）
REP_SER <- unique(RD[, .(primaryid, drug_group, death, hosp, disab, life)],
                  by = c("primaryid","drug_group"))
serious <- REP_SER[, .(n_reports = uniqueN(primaryid),
                  n_death   = sum(death, na.rm = TRUE),
                  n_hosp    = sum(hosp,  na.rm = TRUE),
                  n_disab   = sum(disab, na.rm = TRUE),
                  n_life    = sum(life,  na.rm = TRUE)), by = drug_group]
serious[, pct_death := round(100 * n_death / n_reports, 2)]
serious[, pct_hosp  := round(100 * n_hosp  / n_reports, 2)]
fwrite(serious, file.path(OUT_DIR, "07_严重结局构成.csv"))


## ============================================================================
## 8. 可视化
## ============================================================================

## 8.1 火山图：log2(ROR) vs -log10(校正 P)
if (nrow(signal_all)) {
  d <- signal_all[a >= MIN_CASES]
  d[, p_adj := p.adjust(pchisq(chisq, df = 1, lower.tail = FALSE), method = "BH")]
  d[, neglogp := -log10(p_adj)]
  d[, is_scop := drug == "scopolamine"]

  p1 <- ggplot(d, aes(x = log2(ROR), y = neglogp)) +
    geom_point(aes(color = is_scop, alpha = is_scop), size = 1.1) +
    scale_color_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#95A5A6")) +
    scale_alpha_manual(values = c("TRUE" = 0.85, "FALSE" = 0.28)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#7F8C8D") +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "#7F8C8D") +
    facet_wrap(~ drug, scales = "free_y") +
    labs(x = "log2(ROR)", y = "-log10(校正后 P)") +
    theme_minimal(base_size = 11) + theme(legend.position = "none")
  ggsave(file.path(OUT_DIR, "fig1_火山图.png"), p1, width = 10, height = 6, dpi = 300)
}

## 8.2 森林图：东莨菪碱 Top 20 信号
if (nrow(sig_scop) >= 5) {
  d2 <- head(sig_scop[order(-a)], 20)
  d2[, pt := factor(pt, levels = rev(pt))]
  p2 <- ggplot(d2, aes(x = ROR, y = pt)) +
    geom_point(size = 2.2, color = "#C0392B") +
    geom_errorbarh(aes(xmin = ROR_lcl, xmax = ROR_ucl), height = 0.2, color = "#C0392B") +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_log10() +
    labs(x = "ROR (95% CI, log scale)", y = "") +
    theme_minimal(base_size = 11)
  ggsave(file.path(OUT_DIR, "fig2_东莨菪碱Top20信号.png"), p2,
         width = 8, height = 0.35 * nrow(d2) + 2, dpi = 300)
}

## 8.3 五种药物信号谱热图（Top PT × 药物 的 ROR）
if (nrow(signal_all)) {
  top_pt <- signal_all[drug == "scopolamine" & a >= MIN_CASES][order(-a)][1:min(25, .N), pt]
  hm <- signal_all[pt %in% top_pt & a >= MIN_CASES, .(drug, pt, ROR)]
  hm <- dcast(hm, pt ~ drug, value.var = "ROR")
  hm[, pt := factor(pt, levels = rev(top_pt))]
  hm_l <- melt(hm, id.vars = "pt", variable.name = "drug", value.name = "ROR")
  p3 <- ggplot(hm_l, aes(x = drug, y = pt, fill = log2(ROR))) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#2980B9", mid = "#FDFEFE", high = "#C0392B",
                        midpoint = 0, na.value = "#ECF0F1") +
    labs(x = "", y = "", fill = "log2(ROR)") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(OUT_DIR, "fig3_信号谱比较热图.png"), p3, width = 8, height = 9, dpi = 300)
}


## ============================================================================
## 9. 运行摘要
## ============================================================================

cat("\n==================== 运行摘要 ====================\n")
cat("覆盖季度      : ", length(quarters), " 个（", min(quarters), " - ", max(quarters), "）\n", sep = "")
cat("全库总报告数  : ", format(total_reports, big.mark = ","), "\n", sep = "")
cat("\n各药物报告数：\n")
print(unique(RD[, .(primaryid, drug_group, is_MS, death, hosp)],
             by = c("primaryid","drug_group"))[, .(n_reports = uniqueN(primaryid),
             n_ms_indi = sum(is_MS, na.rm = TRUE),
             n_death   = sum(death, na.rm = TRUE),
             n_hosp    = sum(hosp,  na.rm = TRUE)), by = drug_group][order(-n_reports)])
cat("\n东莨菪碱阳性信号数：", nrow(sig_scop), "\n", sep = "")
if (nrow(sig_scop)) {
  cat("\nTop 15 信号（按报告数）：\n")
  print(head(sig_scop[order(-a), .(pt, a, ROR, ROR_lcl, ROR_ucl, PRR, IC025)], 15))
}
if (nrow(weibull_res)) {
  cat("\nWeibull 时间-发生分析：\n"); print(weibull_res)
}
cat("\n结果已写入：", OUT_DIR, "\n", sep = "")
cat("==================================================\n")

if (TEST_MODE) {
  cat("\n[提示] 当前为测试模式。确认流程无误后，把脚本顶部 TEST_MODE 改为 FALSE，\n",
      "      并下载全部季度数据，再完整运行一次。\n", sep = "")
}
