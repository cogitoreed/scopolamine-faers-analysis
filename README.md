# Route-stratified post-marketing safety analysis of scopolamine in FAERS with external validation in Canada Vigilance

**Repository: https://github.com/cogitoreed/scopolamine-faers-analysis (v1.0.0)**

Code and key output tables for:

> Xu Cai, Jiahui Chen, Ying Yuan, Xin Wang, Wei Gu, Yanli You.
> Post-marketing safety profile of scopolamine: a disproportionality analysis of the FDA Adverse Event Reporting System stratified by route of administration.

## Data sources (all publicly available; no raw data is deposited here)

- **FDA FAERS** quarterly ASCII extract files, 2012Q4–2026Q2 (55 quarters):
  https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html
- **Canada Vigilance** adverse reaction database (full extract):
  https://www.canada.ca/en/health-canada/services/drugs-health-products/medeffect-canada/adverse-reaction-database.html
- **DailyMed** — TRANSDERM SCŌP prescribing information (SPL setid `b877a694-a1d0-4280-937a-a06820b12a88`):
  https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/b877a694-a1d0-4280-937a-a06820b12a88.xml

## File name mapping

The analysis scripts emit result files with Chinese names (the working language of the
project); in this repository the deposited copies use the English names below:

| Script-emitted name (Chinese) | Deposited copy (English) |
|---|---|
| 00_季度处理日志.csv | 00_quarter_processing_log.csv |
| 00_全局去重日志.csv | 00_global_dedup_log.csv |
| 00_全库PT计数.csv | 00_library_PT_counts.csv |
| 01b_给药途径分布.csv | 01b_route_distribution.csv |
| 02_信号检测_全部药物.csv | 02_signal_detection_all_drugs.csv |
| 03_东莨菪碱_阳性信号.csv | 03_scopolamine_positive_signals.csv |
| 03b_东莨菪碱_三法一致敏感性.csv | 03b_scopolamine_trio_agreement.csv |
| 04_敏感性分析_晕动症适应证.csv | 04_sensitivity_motion_sickness_subgroup.csv |
| 05_Weibull时间发生分析.csv | 05_weibull_tto.csv |
| 06_亚组分析.csv | 06_subgroup_analysis.csv |
| 07_严重结局构成.csv | 07_serious_outcomes.csv |
| 12_外部验证对比.csv | 12_external_validation_comparison.csv |
| 20_Table2_主分析信号_Top25.csv | 20_Table2_top25_signals.csv |
| 21_两库对比_剔除非事件.csv | 21_two_database_comparison.csv |
| 22_Table1_基线特征.csv | 22_Table1_baseline.csv |
| 22b_核心PT_总体.csv | 22b_core_PT_overall.csv |
| 23_TTO_重拟_CI.csv | 23_TTO_refit_CI.csv |
| 24_heat_PT_ROR.csv | 24_heat_PT_ROR.csv |
| 24b_heat_PT_TTO汇总.csv | 24b_heat_PT_TTO_summary.csv |
| 24c_heat_TTO_拟合.csv | 24c_heat_TTO_fit.csv |
| 25_途径分层_核心PT.csv | 25_route_stratified_core_PT.csv |
| 26_头对头_tab_headtohead.csv | 26_head_to_head.csv |
| 27_途径xPT类别交叉表.csv | 27_route_x_PT_category.csv |
| 30_敏感性_PSound核心PT.csv | 30_sensitivity_PSonly_core_PT.csv |

## Repository layout

```
code/
  01_download_faers.sh, 02_parallel_download.sh, 02a_dl_one.sh  # FAERS quarterly download (8-way parallel, resumable)
  faers_scopolamine.R        # Main pipeline: read → cross-quarter dedup (FDA algorithm) → drug ID →
                             # route standardization → disproportionality (ROR/PRR/IC) → Weibull TTO →
                             # subgroup/stratified analyses → 5-drug head-to-head
  06b_global_dedup.R         # Cross-quarter deduplication pass (one record per caseid, highest caseversion)
  04_canada_validation.R     # Canada Vigilance external validation
  05_join_validation.R       # FAERS × Canada signal-level join and agreement classification
  06_final_tables.R          # Manuscript Table 2 and two-database comparison
  07_w1prime_tables.R        # Table 1, TTO dual-fit (fitdistr + survreg), route-stratified core PTs,
                             # route × PT-category cross-table, head-to-head table, label-coding inputs
  08_figures.R               # Figures 1–6
  09_label_coding.py         # Structured label review coding (Yes/No/Unclear) against the DailyMed SPL
  10_label_merge.py          # Merge FAERS a/ROR into the label-comparison table
output/
  Key result tables (CSV) and figures (PNG) as reported in the manuscript and supplements.
  Large intermediate files (full report–PT long table, .rds) are omitted; they are regenerated
  by the pipeline from the public sources above.
```

## How to run

> **Note**: `09_label_coding.py` / `10_label_merge.py` additionally require the TRANSDERM SCŌP
> SPL text (download from the DailyMed link above, place as `材料包_EN/label/transderm_scop_plain.txt`)
> and the extracted PLAIN-text conversion of the SPL. All other scripts are self-contained.

```bash
# 0) download FAERS quarters (edit quarter list inside as needed)
bash code/01_download_faers.sh && bash code/02_parallel_download.sh

# 1) cross-quarter deduplication (writes output_v2/global_surviving_primaryids.rds)
Rscript code/06b_global_dedup.R

# 2) main pipeline (55 quarters, ≈20 min on a desktop machine)
Rscript code/faers_scopolamine.R

# 3) external validation, manuscript tables, figures
Rscript code/05_join_validation.R
Rscript code/06_final_tables.R
Rscript code/07_w1prime_tables.R
Rscript code/08_figures.R
```

**Environment**: R 4.4.2 (packages: data.table, stringr, ggplot2, MASS, survival, writexl);
Python 3 (pypdf for SPL text extraction). Run R with UTF-8 locale:
`LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 Rscript ...`.

## Notes

- FAERS quarters before 2012Q4 (legacy AERS format) are excluded because the FDA-recommended
  cross-quarter deduplication algorithm (one record per caseid, highest caseversion) is not
  applicable to them.
- Non-event preferred terms (e.g., OFF LABEL USE, DRUG INEFFECTIVE) are excluded from the
  primary signal set and reported separately.
- This deposit contains analysis code and derived summary tables only; it contains no
  personally identifiable information.

## License

MIT License — see `LICENSE`.
