###############################################################################
## 06_final_tables.R — 论文正式主表生成（剔除非事件 PT + 双库对比终版）
##
## 目的：
##   1) 剔除"报告行为类"非临床事件 PT（OFF LABEL USE、DRUG INEFFECTIVE 等），
##      它们不是不良反应——高质量 FAERS 论文标配，且回应 Khouri 批评
##   2) 生成论文 Table 2（英文列名，Top 25 信号）
##   3) 更新两库对比（复现率按剔除非事件后重算）
##   4) 双向核对 Canada 新信号（AKATHISIA / SUICIDAL IDEATION）在 FAERS 全量的数字
##
## 运行：LC_ALL=en_US.UTF-8 R_LIBS_USER=~/Rlibs Rscript 06_final_tables.R
###############################################################################

suppressMessages(library(data.table))
OUT <- "output_v2"   # run from the project root

faers <- fread(file.path(OUT, "03_东莨菪碱_阳性信号.csv"))
cmp   <- fread(file.path(OUT, "12_外部验证对比.csv"))

## ---------------------------------------------------------------------------
## 1. 非临床事件 PT（报告行为/管理类，非不良反应）——移入补充材料
##    保留 PRODUCT ADHESION ISSUE：透皮贴特有的事件，是递送系统叙事素材
## ---------------------------------------------------------------------------
non_event_pt <- c(
  "OFF LABEL USE",
  "PRODUCT USE IN UNAPPROVED INDICATION",
  "DRUG INEFFECTIVE",
  "DRUG INEFFECTIVE FOR UNAPPROVED INDICATION",
  "WRONG PATIENT RECEIVED PRODUCT",
  "PRODUCT STORAGE ERROR",
  "PRODUCT PREPARATION ISSUE",
  "PRESCRIBING ERROR",
  "DISPENSING ERROR",
  "MEDICATION ERROR"
)
faers[, is_non_event := pt %in% non_event_pt]

main <- faers[is_non_event == FALSE]
supp <- faers[is_non_event == TRUE]

cat("=== PT 分类 ===\n")
cat("总阳性信号:", nrow(faers),
    "| 主分析(临床事件):", nrow(main),
    "| 非事件(补充材料):", nrow(supp), "\n\n")

## ---------------------------------------------------------------------------
## 2. 论文 Table 2：Top 25 主分析信号（英文列名）
## ---------------------------------------------------------------------------
t2 <- head(main[order(-a)],
           min(25, nrow(main)))[, .(
  `Preferred term (PT)`            = pt,
  `Cases (a)`                      = a,
  `ROR`                            = round(ROR, 2),
  `ROR 95% LCL`                    = round(ROR_lcl, 2),
  `ROR 95% UCL`                    = round(ROR_ucl, 2),
  `PRR`                            = round(PRR, 2),
  `IC025`                          = round(IC025, 2),
  `Replicated in Canada Vigilance` = fifelse(pt %in% cmp[verdict == "replicated", pt],
                                             "Yes",
                                             fifelse(pt %in% cmp[verdict == "direction-consistent", pt],
                                                     "Directional", "Not evaluable")))]
fwrite(t2, file.path(OUT, "20_Table2_主分析信号_Top25.csv"))
cat("=== 论文 Table 2（Top 25 临床事件信号）===\n")
print(t2, row.names = FALSE)

## ---------------------------------------------------------------------------
## 3. 剔除非事件后，两库复现统计
## ---------------------------------------------------------------------------
cm <- cmp[pt %in% main$pt]
cat("\n=== 两库对比（剔除非事件 PT 后）===\n")
cat("FAERS 主分析信号:", nrow(main),
    "| Canada 可评:", sum(cm$verdict != "insufficient (a=0)"), "\n")
print(table(Verdict = cm$verdict))
fwrite(cm, file.path(OUT, "21_两库对比_剔除非事件.csv"))

## ---------------------------------------------------------------------------
## 4. 双向核对：Canada 新信号在 FAERS 全量的数字
## ---------------------------------------------------------------------------
cat("\n=== 双向核对：Canada 新信号（FAERS 全量侧）===\n")
chk <- faers[pt %in% c("AKATHISIA", "SUICIDAL IDEATION", "CONFUSIONAL STATE",
                        "HALLUCINATION", "TREMOR", "RESTLESSNESS"),
             .(pt, a, ROR = round(ROR, 1), lcl = round(ROR_lcl, 1),
               signal, in_main = !is_non_event)]
if (nrow(chk)) print(chk) else cat("（FAERS 中未达信号阈值）\n")

## Canada 侧这些 PT 的数字（来自 04 输出）
cv <- fread(file.path(gsub("output_v2$","output",OUT), "11_canada_东莨菪碱PT计数.csv"))
cv2 <- cv[pt %in% c("AKATHISIA", "SUICIDAL IDEATION", "CONFUSIONAL STATE", "HALLUCINATION")]
N_all <- 1153422; n_s <- 114
cv2[, c_ := n_total - a_suspect]
cv2[, d_ := N_all - n_s - c_]
cv2[, ROR_ca := round((as.numeric(a_suspect) * as.numeric(d_)) /
                        (as.numeric(n_s - a_suspect) * as.numeric(c_)), 1)]
cat("\nCanada 侧（Suspect n=114）：\n")
print(cv2[, .(pt, a_suspect, ROR_ca)], row.names = FALSE)

cat("\n输出: 20_Table2_主分析信号_Top25.csv / 21_两库对比_剔除非事件.csv\n")
