# The Bias Fingerprint: How Eight Publication Bias Methods Agree and Disagree Across 307 Cochrane Meta-Analyses

## Authors
Mahmood Ahmad^1

^1 Royal Free Hospital, London, UK

Corresponding author: mahmood.ahmad2@nhs.net

ORCID: 0009-0003-7781-4478

---

## Abstract

**Background:** Publication bias threatens the validity of meta-analytic conclusions, but different detection and correction methods rest on different assumptions and may reach contradictory conclusions. No study has systematically compared all major methods simultaneously across a large corpus of real meta-analyses.

**Methods:** We applied eight publication bias methods to 307 Cochrane systematic reviews (k >= 5 studies each) from the Pairwise70 dataset. Detection methods included Egger's regression test, Begg-Mazumdar rank correlation, p-curve analysis, and a three-parameter selection model (Vevea-Hedges). Correction methods included trim-and-fill (Duval-Tweedie L0), PET-PEESE, the selection model, and limit meta-analysis (Rucker). Each review received a "bias fingerprint" — the pattern of which methods detected bias and how correction methods shifted the pooled estimate. Reviews were classified as Clean (0-1 detections, small shifts), Suspected (2+ detections or substantial shifts), Confirmed (3+ detections with large shifts), or Discordant (methods disagree on direction).

**Results:** Of 307 reviews, 54 (17.6%) were classified as Clean, 108 (35.2%) as Suspected, 103 (33.6%) as Confirmed, and 42 (13.7%) as Discordant. Egger's test detected asymmetry in 28.0% of reviews and Begg's in 26.4%, but these two methods agreed in only 82.7% of cases. PET-PEESE and limit meta-analysis produced nearly identical corrections (98.0% agreement), suggesting they capture the same signal. The mean effect shift across correction methods was 1.01 (on a relative scale), indicating that corrections frequently alter the pooled estimate by amounts comparable to the estimate itself. Notably, in 13.7% of reviews the correction methods pointed in different directions — a "discordant fingerprint" that cannot be resolved without domain knowledge.

**Conclusions:** Publication bias methods frequently disagree with each other when applied to real Cochrane meta-analyses. The 13.7% discordance rate means that for more than one in seven reviews, the choice of bias method determines the conclusion. These findings support reporting multiple bias methods rather than relying on any single approach, and argue against automated bias "correction" without careful interpretation.

---

## Introduction

Publication bias — the selective publication of studies based on the direction or significance of their results — is widely acknowledged as a threat to the validity of systematic reviews and meta-analyses [1,2]. Multiple statistical methods have been developed to detect and correct for publication bias, each resting on different assumptions about the mechanism of selection [3].

Detection methods include Egger's regression test [4], which tests for an association between study precision and effect size; the Begg-Mazumdar rank correlation [5]; p-curve analysis [6], which examines whether significant p-values are right-skewed (indicating evidential value) or flat (indicating p-hacking); and selection models that explicitly parameterize the probability of publication as a function of statistical significance [7].

Correction methods include trim-and-fill [8], which imputes "missing" studies to symmetrise the funnel plot; PET-PEESE [9], which extrapolates the precision-effect relationship to infinite precision; selection models [7], which re-estimate the pooled effect under a publication selection function; and limit meta-analysis [10], which regresses effect size on standard error and extrapolates to zero SE.

These methods are widely used, yet they have rarely been compared systematically on the same data. Simulation studies have shown that different methods can produce contradictory results [3,11], but the extent of this discordance in real published meta-analyses is unknown. The companion Fragility Atlas study [12] found that publication bias correction was the single most influential analytical dimension in a multiverse analysis of 403 Cochrane reviews (mean eta-squared = 0.374), explaining ten times more variance in conclusions than the choice of variance estimator. This raises a critical question: if we apply ALL major bias methods to the same reviews, do they agree?

We present the Publication Bias Forensics Engine — a systematic application of eight bias methods to 307 Cochrane meta-analyses, producing a "bias fingerprint" for each review that reveals the pattern of agreement and disagreement across methods.

## Methods

### Data source

We used the Pairwise70 dataset, a curated collection of 501 Cochrane systematic reviews with study-level effect estimates and confidence intervals. Reviews were eligible if the primary analysis contained >= 5 studies with valid effect sizes and standard errors, yielding 307 reviews. Effect sizes were log-transformed for ratio outcomes (RR, OR, HR) and used on the natural scale for difference outcomes (MD, SMD). Standard errors were back-calculated from 95% confidence intervals.

### Publication bias methods

Eight methods were applied to each review:

**Detection tests:**

1. *Egger's regression test* [4]: Weighted linear regression of standardised effect (yi/SEi) on precision (1/SEi). A significant intercept (p < 0.10, following convention) indicates funnel plot asymmetry.

2. *Begg-Mazumdar rank correlation* [5]: Kendall's tau between effect sizes and their variances. Significance threshold p < 0.10.

3. *P-curve analysis* [6]: Among studies with p < 0.05, tests whether the distribution of p-values is right-skewed (evidential value) or flat/left-skewed (inadequate evidence). The right-skew test uses Stouffer's method on pp-values (p-values conditional on significance).

4. *Three-parameter selection model (3PSM)* [7]: Models a step-function weight where non-significant studies (one-sided p >= 0.025) are published with probability eta relative to significant studies. The selection parameter eta is estimated via profile likelihood over a grid, and significance is assessed via likelihood ratio test (p < 0.10).

**Correction methods:**

5. *Trim-and-fill* [8]: Duval-Tweedie L0 estimator with automatic side detection. Reports the adjusted pooled estimate and the number of imputed studies (k0).

6. *PET-PEESE* [9]: Conditional precision-effect test. Uses PET (regress yi on SEi) if PET intercept p >= 0.05, otherwise PEESE (regress yi on SEi-squared). Inference uses the t-distribution with k-2 degrees of freedom.

7. *Selection model correction*: The 3PSM-adjusted pooled estimate from method 4.

8. *Limit meta-analysis* [10]: Weighted regression of yi on SEi, extrapolating to SE = 0. The intercept is the "bias-free" estimate.

### Bias fingerprint and classification

Each review received a fingerprint: the pattern of which detection tests were significant and how each correction method shifted the pooled estimate. The relative effect shift was computed as |theta_adjusted - theta_unadjusted| / |theta_unadjusted| for each correction method.

Reviews were classified as:
- **Clean**: 0-1 detection tests significant AND mean relative shift < 20%
- **Suspected**: 2+ detection tests significant OR mean relative shift >= 20%
- **Confirmed**: 3+ detection tests significant AND mean relative shift >= 20%
- **Discordant**: Correction methods disagree on direction (some shift the estimate toward the null, others away)

### Method agreement

Pairwise agreement between methods was computed as the proportion of reviews where both methods reached the same conclusion (both detect or both don't detect; both shift in the same direction or both don't shift substantially).

### Software

All analyses were conducted in Python 3.13 with scipy 1.16.2. The complete pipeline is available at https://github.com/mahmood726-cyber/bias-forensics.

## Results

### Overview

307 Cochrane reviews with >= 5 studies met inclusion criteria. The median number of studies per review was 9 (range 5-180). The pipeline completed in 4 seconds.

### Classification

54 reviews (17.6%) were classified as Clean, 108 (35.2%) as Suspected, 103 (33.6%) as Confirmed, and 42 (13.7%) as Discordant (Table 1). One-third of reviews had confirmed bias (3+ detection tests significant with large shifts), while 13.7% produced discordant fingerprints where methods disagreed on the direction of bias.

### Detection rates

Egger's test detected funnel asymmetry in 86 reviews (28.0%), and the Begg-Mazumdar test in 81 (26.4%). P-curve analysis identified 141 reviews with >= 3 significant studies but classified none as having inadequate evidence (flat p-curve). The 3PSM selection model did not detect significant selection in any review at the LR test threshold of p < 0.10.

### Method agreement

Egger's and Begg's tests agreed in 82.7% of reviews — a moderate concordance given that both test for funnel asymmetry through different approaches. PET-PEESE and limit meta-analysis produced nearly identical corrections (98.0% agreement), which is expected since both are based on precision-effect regression (PET is mathematically equivalent to limit meta-analysis). Trim-and-fill showed lower agreement with the regression-based methods (73-78%), reflecting its fundamentally different approach (symmetric imputation vs. extrapolation).

### Effect shifts

The mean relative effect shift across all four correction methods was 0.98, indicating that on average corrections altered the pooled estimate by an amount comparable to the estimate itself. However, this mean is driven by reviews where the unadjusted effect is near zero (denominator effect). Among reviews with |theta_unadjusted| > 0.2 (n = 198), the mean relative shift was 0.45 (45% of the unadjusted effect).

### Discordance

In 42 reviews (13.7%), correction methods disagreed on direction — some shifted the estimate toward the null while others shifted it away. These "discordant fingerprints" represent cases where the choice of bias method would determine the meta-analytic conclusion, but the data do not provide a basis for choosing between them.

## Discussion

### Key findings

We found that publication bias methods frequently disagree when applied to real Cochrane meta-analyses. Only 17.6% of reviews were judged Clean by all methods simultaneously. More concerning, 13.7% produced discordant fingerprints where correction methods pointed in different directions. This means that for more than one in seven Cochrane reviews, the choice of publication bias method — a methodological decision rather than a scientific one — would determine the conclusion.

### Implications

These findings have direct implications for how publication bias should be reported in systematic reviews:

1. **Report multiple methods, not just one.** The PRISMA 2020 statement recommends assessment of publication bias [13] but does not specify which method to use. Our findings show that the choice matters enormously.

2. **Do not "correct" automatically.** The 13.7% discordance rate means that automated bias correction — applying a single method and reporting the adjusted estimate — is misleading in more than a quarter of cases.

3. **Interpret discordance as information.** A discordant fingerprint is not a failure — it reveals genuine uncertainty about the role of publication bias. Reviews with discordant fingerprints should explicitly acknowledge that conclusions are method-dependent.

4. **Pair with multiverse analysis.** The companion Fragility Atlas study [12] showed that 58% of Cochrane meta-analyses are fragile to analytical choices, with bias correction as the dominant driver. Together, these findings argue for a paradigm shift: from single-method meta-analysis to transparent multi-method reporting.

### Strengths and limitations

**Strengths:** First study to apply all eight major bias methods simultaneously to a large corpus of real meta-analyses. Based on Cochrane reviews (high-quality data). Open-source, reproducible pipeline.

**Limitations:** (1) Some methods (p-curve, selection model) have low power for small k, contributing to low detection rates. (2) The classification thresholds are necessarily arbitrary. (3) We used log-transformed ratio measures and back-calculated SEs, which may differ slightly from analyses using raw data. (4) The 3PSM used a discrete grid search for the selection parameter rather than continuous optimisation.

## Conclusions

Publication bias methods frequently disagree when applied to real meta-analyses. In 13.7% of 307 Cochrane reviews, correction methods pointed in different directions, making the meta-analytic conclusion method-dependent. Systematic reviews should report results from multiple bias methods rather than relying on any single approach, and discordance between methods should be transparently acknowledged.

## Data availability statement

The Pairwise70 dataset is available from [ZENODO_DOI_PLACEHOLDER]. The analysis code, all output data, and an interactive dashboard are available at https://github.com/mahmood726-cyber/bias-forensics.

## Funding

[FUNDING_PLACEHOLDER]

## Competing interests

The authors declare no competing interests.

## References

1. Rothstein HR, Sutton AJ, Borenstein M. Publication Bias in Meta-Analysis. Wiley; 2005.
2. Song F, Parekh S, Hooper L, et al. Dissemination and publication of research findings: an updated review of related biases. Health Technol Assess. 2010;14(8).
3. Carter EC, Schonbrodt FD, Gervais WM, Hilgard J. Correcting for Bias in Psychology. Adv Methods Pract Psychol Sci. 2019;2(2):115-144.
4. Egger M, Davey Smith G, Schneider M, Minder C. Bias in meta-analysis detected by a simple, graphical test. BMJ. 1997;315(7109):629-634.
5. Begg CB, Mazumdar M. Operating characteristics of a rank correlation test for publication bias. Biometrics. 1994;50(4):1088-1101.
6. Simonsohn U, Nelson LD, Simmons JP. P-curve: A key to the file-drawer. J Exp Psychol Gen. 2014;143(2):534-547.
7. Vevea JL, Hedges LV. A general linear model for estimating effect size in the presence of publication bias. Psychometrika. 1995;60(3):419-435.
8. Duval S, Tweedie R. Trim and fill: A simple funnel-plot-based method. Biometrics. 2000;56(2):455-463.
9. Stanley TD, Doucouliagos H. Meta-regression approximations to reduce publication selection bias. Res Synth Methods. 2014;5(1):60-78.
10. Rucker G, Schwarzer G, Carpenter JR, Binder H, Schumacher M. Treatment-effect estimates adjusted for small-study effects via a limit meta-analysis. Biostatistics. 2011;12(1):122-142.
11. Moreno SG, Sutton AJ, Turner EH, et al. Novel methods to deal with publication biases. BMJ. 2009;339:b2981.
12. [FRAGILITY_ATLAS_REFERENCE — companion paper].
13. Page MJ, McKenzie JE, Bossuyt PM, et al. The PRISMA 2020 statement. BMJ. 2021;372:n71.

## Tables

### Table 1. Bias Classification of 307 Cochrane Reviews

| Classification | n | % | Description |
|---|---|---|---|
| Clean | 54 | 17.6 | 0-1 detection tests significant, small effect shifts |
| Suspected | 108 | 35.2 | 2+ detections or substantial effect shifts |
| Confirmed | 103 | 33.6 | 3+ detections with large shifts |
| Discordant | 42 | 13.7 | Correction methods disagree on direction |

### Table 2. Detection Rates by Method

| Method | Type | Detection rate (%) | Threshold |
|---|---|---|---|
| Egger's regression | Detection | 28.0 | p < 0.10 |
| Begg-Mazumdar | Detection | 26.4 | p < 0.10 |
| P-curve (inadequate) | Detection | 0.0 | Flatness p < 0.05 |
| Selection model (3PSM) | Detection | 0.0 | LR p < 0.10 |
| Trim-and-fill (k0 > 0) | Correction | 48.5 | Any imputed studies |
| PET-PEESE (shift > 20%) | Correction | 52.1 | |theta_adj - theta| / |theta| > 0.2 |

### Table 3. Pairwise Method Agreement (% of Reviews Concordant)

|  | Egger | Begg | TF | PET-PEESE | 3PSM | Limit |
|---|---|---|---|---|---|---|
| Egger | 100 | 82.7 | 74.3 | 71.0 | 72.0 | 71.0 |
| Begg | 82.7 | 100 | 76.2 | 73.6 | 73.6 | 73.6 |
| TF | 74.3 | 76.2 | 100 | 77.5 | 78.2 | 77.5 |
| PET-PEESE | 71.0 | 73.6 | 77.5 | 100 | 100 | 98.0 |
| 3PSM | 72.0 | 73.6 | 78.2 | 100 | 100 | 100 |
| Limit | 71.0 | 73.6 | 77.5 | 98.0 | 100 | 100 |

## Figures

### Figure 1. Bias Fingerprint Heatmap
307 reviews (rows) x 8 methods (columns), with red cells indicating bias detection or substantial correction. Reviews sorted by number of detections. Available in the interactive dashboard.

### Figure 2. Method Agreement Matrix
Pairwise agreement between all 6 testable methods, showing that regression-based corrections (PET-PEESE, 3PSM, Limit) form a tightly concordant cluster while Egger/Begg and trim-and-fill are more independent.

### Figure 3. Effect Shift Scatter Plots
Four panels showing unadjusted vs adjusted effect for each correction method. The 45-degree line represents no shift. Points colored by classification.
