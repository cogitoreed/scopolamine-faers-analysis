###############################################################################
## Canada Vigilance 外部验证分析（东莨菪碱）
##
## 目的：对 FAERS 主分析检出的信号，在 Canada Vigilance 独立国家级数据库
##       中做方向性外部验证（directional external validation）。
##
## 数据：data/canada/extracted/cvponline_extract_20241130/
##       reports.txt      全库报告主表（1,154,017 报告，1965–2024-11）
##       report_drug.txt  药物表（列: REPORT_DRUG_ID, REPORT_ID, DRUG_PRODUCT_ID,
##                                DRUG_NAME, ROLE_EN, ...）
##       reactions.txt    反应表（列: REPORT_DRUG_ID, REPORT_ID, ..., PT_EN, PT_FR, SOC...）
##
## 方法：
##   1) 暴露定义：DRUG_NAME 匹配 scopolamine/hyoscine，排除 methscopolamine、
##      丁溴（butylscopolamine/buscopan）及复方（含 atropine/hyoscyamine/phenobarbital）
##   2) 主分析: ROLE_EN == "Suspect"；敏感性: 任意角色
##   3) 输出：scop 报告的 PT 计数（a 侧）+ 全库 PT 报告数（c 侧）
##      → 与 FAERS 信号表 join 后计算 Canada ROR（在 05_join_validation.R 完成）
##
## 运行：LC_ALL=en_US.UTF-8 R_LIBS_USER=~/Rlibs Rscript 04_canada_validation.R
## 实测：已于 2026-09-04 全流程跑通（约 3 分钟）
###############################################################################

suppressMessages(library(data.table))
setDTthreads(0)

CV_DIR <- "/Users/xiaocaixu/WorkBuddy/学位自救计划/FAERS分析/data/canada/extracted/cvponline_extract_20241130"
OUT    <- "/Users/xiaocaixu/WorkBuddy/学位自救计划/FAERS分析/output"
dir.create(OUT, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## 1. 读取药物表，定义东莨菪碱暴露
## ---------------------------------------------------------------------------
message("[1/4] 读取 report_drug.txt (909MB) ...")
## 注意：字段值内嵌双引号（如 VOLTAREN "NOVARTIS"）+ 法文 Latin-1 字节，
## 必须 quote=""（把引号当字面量，仅用 $ 分界）且不能用 trimws（坏字节致
## sub() 报 invalid UTF-8，已实测踩坑），引号用 gsub 手工去。
rd <- fread(file.path(CV_DIR, "report_drug.txt"), sep = "$", quote = "",
            header = FALSE, encoding = "Latin-1", fill = TRUE,
            colClasses = "character", showProgress = FALSE,
            select = c(1, 2, 4, 5))
setnames(rd, c("report_drug_id", "report_id", "drug_name", "role"))

rd[, drug_u := toupper(gsub('"', '', iconv(drug_name, "UTF-8", "ASCII", sub = "")))]
rd[, role    := gsub('"', '', role)]          ## 角色列同样要清洗引号（漏了会 Suspect=0，已踩坑）

## 诊断：行数一致性与匹配分布（校验 fread 与 shell grep 谁对）
message("    report_drug.txt 行数: ", nrow(rd))
message("    drug_name 含 scopolamine 行: ", sum(grepl("SCOPOLAMIN|HYOSCIN", rd$drug_u)))

## 暴露：含 scopolamine/hyoscine
is_scop <- grepl("SCOPOLAMIN|HYOSCIN", rd$drug_u)
## 排除：季铵类（不透过 BBB）、复方（含阿托品/莨菪碱/苯巴比妥，归因混淆）
is_excl <- grepl("METHSCOPOLAMINE|BUTYLBROMIDE|BUTYLCOPOLAMINE|BUSCOPAN", rd$drug_u) |
           grepl("ATROPINE|HYOSCYAMINE|PHENOBARBITAL", rd$drug_u)

scop <- rd[is_scop & !is_excl]
message("    东莨菪碱药物行: ", nrow(scop),
        " | 唯一报告: ", uniqueN(scop$report_id),
        " | Suspect 行: ", sum(scop$role == "Suspect"))

## ---------------------------------------------------------------------------
## 2. 暴露报告集合（主分析 Suspect；敏感性 任意角色）
## ---------------------------------------------------------------------------
rep_suspect <- scop[role == "Suspect", unique(report_id)]
rep_anyrole <- scop[, unique(report_id)]
message("[2/4] Suspect 报告: ", length(rep_suspect),
        " | 任意角色报告: ", length(rep_anyrole))

## ---------------------------------------------------------------------------
## 3. 读取反应表，提取暴露报告的 PT（a 侧）
## ---------------------------------------------------------------------------
message("[3/4] 读取 reactions.txt (708MB) ...")
reac <- fread(file.path(CV_DIR, "reactions.txt"), sep = "$", quote = "",
              header = FALSE, encoding = "Latin-1", fill = TRUE,
              colClasses = "character", showProgress = FALSE,
              select = c(2, 6))
setnames(reac, c("report_id", "pt"))
reac <- reac[pt != "" & !is.na(pt)]
reac[, pt := toupper(gsub('"', '', iconv(pt, "UTF-8", "ASCII", sub = "")))]
reac <- unique(reac, by = c("report_id", "pt"))   # 报告×PT 去重

## 全库 PT 计数（c 侧分母）
pt_all <- reac[, .(n_total = uniqueN(report_id)), by = pt]
fwrite(pt_all, file.path(OUT, "10_canada_全库PT计数.csv"))

## 东莨菪碱报告的 PT 计数（主分析 = Suspect）
scop_pt_susp <- reac[report_id %in% rep_suspect, .(a_suspect = uniqueN(report_id)), by = pt]
scop_pt_any  <- reac[report_id %in% rep_anyrole, .(a_anyrole = uniqueN(report_id)), by = pt]

## ---------------------------------------------------------------------------
## 4. 汇总输出（供与 FAERS 信号 join）
## ---------------------------------------------------------------------------
message("[4/4] 汇总输出 ...")
canada_scop <- merge(scop_pt_any, scop_pt_susp, by = "pt", all = TRUE)
canada_scop[is.na(a_suspect), a_suspect := 0]
canada_scop[is.na(a_anyrole),  a_anyrole  := 0]
canada_scop <- merge(canada_scop, pt_all, by = "pt")
canada_scop[order(-a_suspect)]
fwrite(canada_scop, file.path(OUT, "11_canada_东莨菪碱PT计数.csv"))

## 运行摘要
cat("\n================ Canada Vigilance 验证摘要 ================\n")
cat("全库报告数        : ", format(uniqueN(reac$report_id), big.mark = ","), "\n", sep = "")
cat("东莨菪碱报告      : Suspect ", length(rep_suspect),
    " | 任意角色 ", length(rep_anyrole), "\n", sep = "")
cat("Suspect 报告 PT 谱 Top 12：\n")
print(head(canada_scop[order(-a_suspect), .(pt, a_suspect, a_anyrole, n_total)], 12))
cat("\n输出: 10_canada_全库PT计数.csv / 11_canada_东莨菪碱PT计数.csv\n")
cat("（待 FAERS 全量信号表出来后，由 05_join_validation.R 计算 Canada ROR 与同向性）\n")
