# Changelog

## v1.0.2 (2026-09-05)
- Synchronized output/figures with the manuscript (Figure 4 redesigned: clinical-cluster events per 1,000
  route-specific reports + PT x route ROR landscape heatmap; Figure 5 delta-method confidence bands;
  Figures 6-8 legend placement and CI refinements).
- No changes to analysis code or result tables relative to v1.0.1.

## v1.0.1 (2026-09-05)
- Fixed: survreg Weibull shape confidence intervals were computed from the wrong variance element
  (Var(intercept) instead of Var(Log(scale))). All survreg CIs corrected; the pooled heat-related TTO
  shape 95% CI changes from 0.249-1.039 (v1.0.0, incorrect) to 0.405-0.638. Independent large-sample
  check (SE(log shape) ~ 1.06/sqrt(3(n-2))) reproduces the corrected intervals.
- Figure 5 confidence bands recomputed via the delta method on (mu, log sigma).
- Regenerated tables: 23_TTO_refit_CI.csv, 24c_heat_TTO_fit.csv.

## v1.0.0 (2026-09-05)
- Initial archived release. Note: the v1.0.0 copies of 23_TTO_refit_CI.csv and 24c_heat_TTO_fit.csv
  contain the incorrect survreg CIs described above; cite v1.0.1 or later.
