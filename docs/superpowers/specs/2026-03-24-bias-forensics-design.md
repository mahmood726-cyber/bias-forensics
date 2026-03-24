# Publication Bias Forensics Engine — Design Spec

## 1. Problem

Publication bias is the most consequential methodological concern in meta-analysis — the Fragility Atlas showed it drives 10x more disagreement than estimator choice (η²=0.374). Yet each bias detection method has different assumptions, different power, and different false positive rates. No tool runs ALL methods simultaneously and compares them. Reviewers pick one method and hope for the best.

## 2. Goal

Run 8 publication bias detection/correction methods on the same 403 Cochrane reviews from the Fragility Atlas, producing a "bias fingerprint" for each review: which methods detect bias, which don't, and what drives disagreement.

## 3. Target

- **Journal:** Research Synthesis Methods (companion paper to the Fragility Atlas BMJ paper)
- **Location:** `C:\BiasForensics\`
- **Architecture:** Python pipeline (reuses Fragility Atlas data) + single-file HTML dashboard

## 4. The 8 Methods

### Detection tests (binary: significant/not)
1. **Egger's regression** — regress standardized effect on precision (1/SE). Significant intercept → asymmetry.
2. **Begg-Mazumdar rank correlation** — Kendall's tau between effect size and variance. Significant → asymmetry.
3. **P-curve** — distribution of p-values from significant studies. Right-skewed = evidential value, flat/left-skewed = bias or p-hacking.
4. **Z-curve** — observed vs expected discovery rate (EDR). Low EDR → selection bias.

### Correction methods (adjusted effect estimate)
5. **Trim-and-fill** — Duval-Tweedie L0 (already in Fragility Atlas). Adjusted theta.
6. **PET-PEESE** — conditional precision-effect test (already in Fragility Atlas). Adjusted theta.
7. **Selection model (Vevea-Hedges 3PSM)** — weight function model assuming step-function selection at p=0.025. Adjusted theta + likelihood ratio test.
8. **Limit meta-analysis (Rucker)** — regress effect on SE, extrapolate to SE=0. The "bias-free" intercept.

## 5. Pipeline Architecture

```
src/
  loader.py          — Load Fragility Atlas results + Pairwise70 data (reuse)
  egger.py           — Egger's regression test
  begg.py            — Begg-Mazumdar rank correlation
  pcurve.py          — P-curve analysis
  zcurve.py          — Z-curve (observed/expected discovery rate)
  selection_model.py — 3-parameter selection model (Vevea-Hedges)
  limit_meta.py      — Rucker's limit meta-analysis
  forensics.py       — Run all 8 methods, produce fingerprint
  pipeline.py        — Orchestrate: load → run all → classify → export
```

Trim-and-fill and PET-PEESE are already implemented in the Fragility Atlas — import from there or copy.

## 6. Output: The Bias Fingerprint

For each of 403 reviews, a row with:

| Field | Description |
|-------|-------------|
| review_id | Cochrane review ID |
| k | Number of studies |
| egger_p | Egger's test p-value |
| egger_sig | 1 if p < 0.10 (standard threshold) |
| begg_p | Begg-Mazumdar p-value |
| begg_sig | 1 if p < 0.10 |
| pcurve_z | P-curve Z-statistic (right-skew test) |
| pcurve_evidential | 1 if right-skew test p < 0.05 |
| pcurve_flat | 1 if flatness test p < 0.05 (lack of evidence) |
| zcurve_edr | Expected discovery rate (0-1) |
| zcurve_oir | Observed inclusion rate |
| tf_theta | Trim-and-fill adjusted effect |
| tf_k0 | Number of imputed studies |
| petpeese_theta | PET-PEESE adjusted effect |
| sel3psm_theta | 3PSM adjusted effect |
| sel3psm_p | Selection model LR test p-value |
| limit_theta | Limit meta-analysis intercept |
| n_methods_detect | Count of methods detecting bias (0-4 for tests) |
| concordance | % of correction methods agreeing on direction+significance |
| bias_class | Clean / Suspected / Confirmed / Discordant |

### Classification
- **Clean**: 0-1 tests detect bias, correction methods agree with unadjusted
- **Suspected**: 2 tests detect bias, or corrections shift effect >20%
- **Confirmed**: 3-4 tests detect bias AND corrections substantially change effect
- **Discordant**: Tests disagree (some detect, some don't) AND corrections point in different directions

## 7. Dashboard

Single-file HTML (`dashboard/index.html`) with embedded results:
- Overview: pie chart of Clean/Suspected/Confirmed/Discordant
- Heatmap: 403 reviews × 8 methods, colored by detection/direction
- Method agreement matrix: which methods tend to agree?
- Predictors: what review characteristics predict bias detection?
- Drill-down: click any review → full fingerprint

## 8. Validation

1. Egger's and Begg's verified against R `metafor::regtest()` and `metafor::ranktest()` on 10 reviews
2. P-curve verified against p-curve.com reference implementation
3. Selection model verified against `weightr::weightfunct()` on 5 reviews
4. All 403 reviews should complete without errors

## 9. Success Criteria

1. Pipeline processes all 403 reviews without crashes
2. R cross-validation matches on 10/10 test reviews
3. Classification produces non-trivial distribution (not all one category)
4. Dashboard loads in <3 seconds with all 403 reviews
5. Manuscript-ready tables and figures
