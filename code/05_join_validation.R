###############################################################################
## 05_join_validation.R — FAERS × Canada Vigilance 外部验证对比表
##
## 输入：
##   output/03_东莨菪碱_阳性信号.csv        FAERS 信号（pt, a, ROR, ROR_lcl, signal, sens_trio）
##   output/10_canada_全库PT计数.csv        Canada 每 PT 全库报告数
##   output/11_canada_东莨菪碱PT计数.csv    Canada 东莨菪碱 PT 计数（a_suspect/a_anyrole/n_total）
## 输出：
##   output/12_外部验证对比.csv             两库并列 ROR + 复现判定
##
## 判定标准：
##   复现   = FAERS ROR_lcl>1 且 Canada ROR_lcl>1（两库独立均为阳性信号）
##   部分复现 = 两库点估计同向（ROR 同 >1 或同 <1）但一侧 CI 含 1
##   未复现   = 方向相反或 Canada a=0（样本不足，标注 insufficient）
##
## 运行：LC_ALL=en_US.UTF-8 R_LIBS_USER=~/Rlibs Rscript 05_join_validation.R
###############################################################################

suppressMessages(library(data.table))
OUT <- "output_v2"   # run from the project root
OUT_V1 <- "output"   ## Canada-side tables deposited under output/

faers <- fread(file.path(OUT, "03_东莨菪碱_阳性信号.csv"))
cv    <- fread(file.path(OUT_V1, "11_canada_东莨菪碱PT计数.csv"))

## Canada 侧参数（与 04 脚本一致；如重跑 04 后数字变化请同步）
n_scop_canada <- 114      # Canada Suspect 报告数
N_all_canada  <- 1153422  # Canada 全库报告数

## ---------------------------------------------------------------------------
## 1. Canada 侧 ROR（四格表：a=scop Suspect 中该PT；b=n_scop-a；
##    c=该PT全库报告数-a；d=N_all-n_scop-c）
## ---------------------------------------------------------------------------
setnames(cv, "n_total", "n_pt_all")
cv[, a_c := a_suspect]
cv[, c_ := n_pt_all - a_c]
cv[, d_ := N_all_canada - n_scop_canada - c_]
cv[a_c > 0,
   `:=`(ror_c = (as.numeric(a_c) * as.numeric(d_)) /
          (as.numeric(n_scop_canada - a_c) * as.numeric(c_)),
        se_c  = sqrt(1/a_c + 1/(n_scop_canada-a_c) + 1/c_ + 1/d_))]
cv[a_c > 0, ror_c_lcl := exp(log(ror_c) - 1.96 * se_c)]
## a=0 的 PT：Canada 无报告，无法计算（insufficient）

## ---------------------------------------------------------------------------
## 2. join + 复现判定
## ---------------------------------------------------------------------------
j <- merge(faers[, .(pt, faers_a = a, faers_ror = ROR, faers_lcl = ROR_lcl,
                     faers_signal = signal, faers_trio = sens_trio)],
           cv[, .(pt, canada_a = a_c, canada_n_pt_all = n_pt_all,
                  canada_ror = ror_c, canada_lcl = ror_c_lcl)],
           by = "pt", all.x = TRUE)

j[, canada_a := fifelse(is.na(canada_a), 0L, canada_a)]
j[is.na(canada_ror) & canada_a == 0, verdict := "insufficient (a=0)"]
j[!is.na(canada_ror),
  verdict := fcase(
    faers_lcl > 1 & canada_lcl > 1,                         "replicated",
    (faers_ror > 1) == (canada_ror > 1),                    "direction-consistent",
    default = "discordant")]
j[faers_lcl <= 1 & !is.na(canada_ror),
  verdict := fifelse((faers_ror > 1) == (canada_ror > 1),
                     "direction-consistent (FAERS lcl<=1)", "discordant")]

setorder(j, -faers_a)
fwrite(j, file.path(OUT, "12_外部验证对比.csv"))

## ---------------------------------------------------------------------------
## 3. 摘要
## ---------------------------------------------------------------------------
cat("\n=============== FAERS × Canada 外部验证摘要 ===============\n")
cat("FAERS 信号数:", nrow(j),
    " | Canada 可评:", sum(j$verdict != "insufficient (a=0)"), "\n")
print(table(Verdict = j$verdict))
cat("\n复现信号明细（Canada ROR 排序）：\n")
print(head(j[verdict == "replicated",
            .(pt, faers_a, faers_ror = round(faers_ror, 1),
              canada_a, canada_ror = round(canada_ror, 1),
              canada_lcl = round(canada_lcl, 1))], 20))
cat("\n输出: 12_外部验证对比.csv\n")
