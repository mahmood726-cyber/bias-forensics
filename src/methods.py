"""All 8 publication bias detection/correction methods.

Each method takes yi (effect sizes) and sei (standard errors) as numpy arrays
and returns a dict with method-specific results.
"""

import math
import numpy as np
from scipy import stats


# ─── 1. Egger's Regression Test ───

def egger_test(yi, sei):
    """Egger's regression test for funnel plot asymmetry.

    Regress standardized effect (yi/sei) on precision (1/sei).
    A significant intercept indicates asymmetry.
    """
    k = len(yi)
    if k < 3:
        return {'p': 1.0, 'intercept': 0, 'se': 0, 'significant': False}

    # Standard Egger: regress yi/sei on 1/sei (equivalent to WLS of yi on sei)
    x = 1.0 / sei  # precision
    y = yi / sei    # standardized effect

    # Guard: if all precisions are identical, linregress will fail
    if np.std(x) < 1e-10:
        return {'p': 1.0, 'intercept': 0, 'se': 0, 't': 0, 'significant': False}

    slope, intercept, r, p, se_slope = stats.linregress(x, y)

    # The intercept's significance test (the actual Egger test)
    # SE of intercept from linear regression
    n = len(x)
    x_mean = np.mean(x)
    ss_x = np.sum((x - x_mean) ** 2)
    residuals = y - (slope * x + intercept)
    mse = np.sum(residuals ** 2) / (n - 2) if n > 2 else 1
    se_intercept = math.sqrt(mse * (1.0 / n + x_mean ** 2 / ss_x)) if ss_x > 0 else 1

    t_stat = intercept / se_intercept if se_intercept > 0 else 0
    df = max(1, n - 2)
    p_intercept = 2 * (1 - stats.t.cdf(abs(t_stat), df))

    return {
        'p': float(p_intercept),
        'intercept': float(intercept),
        'se': float(se_intercept),
        't': float(t_stat),
        'significant': p_intercept < 0.10  # standard threshold for Egger
    }


# ─── 2. Begg-Mazumdar Rank Correlation ───

def begg_test(yi, sei):
    """Begg-Mazumdar rank correlation test.

    Kendall's tau between CENTERED effect size and variance (P1-1 fix).
    """
    k = len(yi)
    if k < 3:
        return {'p': 1.0, 'tau': 0, 'significant': False}

    vi = sei ** 2
    # P1-1 FIX: center effects before rank correlation (matches metafor::ranktest)
    wi = 1.0 / vi
    theta_fe = np.sum(wi * yi) / np.sum(wi)
    tau, p = stats.kendalltau(yi - theta_fe, vi)

    return {
        'p': float(p) if not math.isnan(p) else 1.0,
        'tau': float(tau) if not math.isnan(tau) else 0.0,
        'significant': (float(p) if not math.isnan(p) else 1.0) < 0.10
    }


# ─── 3. P-curve Analysis ───

def p_curve(yi, sei):
    """P-curve analysis: tests whether significant p-values are right-skewed.

    Right-skewed (clustered near 0) → evidential value (real effect).
    Flat/left-skewed → no evidential value (p-hacking or null).
    """
    k = len(yi)
    if k < 3:
        return {'z_right': 0, 'p_right': 1.0, 'z_flat': 0, 'p_flat': 1.0,
                'evidential': False, 'inadequate': False, 'n_sig': 0}

    # Compute two-sided p-values for each study
    z_vals = yi / sei
    p_vals = 2 * (stats.norm.sf(np.abs(z_vals)))

    # Select only significant studies (p < 0.05)
    sig_mask = p_vals < 0.05
    sig_p = p_vals[sig_mask].copy()
    n_sig = len(sig_p)

    if n_sig < 3:
        return {'z_right': 0, 'p_right': 1.0, 'z_flat': 0, 'p_flat': 1.0,
                'evidential': False, 'inadequate': False, 'n_sig': n_sig}

    # Transform p-values to pp-values (conditional on p < 0.05)
    # Under H0: pp = p/0.05 ~ Uniform(0,1)
    # Right-skew test: are pp-values smaller than expected under uniform?
    pp_vals = sig_p / 0.05

    # Stouffer's method for right-skew (are pp-values below 0.5?)
    # Each pp-value transformed: z_i = qnorm(pp_i)
    z_pp = stats.norm.ppf(pp_vals)
    z_right = float(np.sum(z_pp) / math.sqrt(n_sig))
    p_right = float(stats.norm.cdf(z_right))  # one-sided: small = right-skewed

    # Flatness test: are pp-values uniformly distributed? (lack of evidence)
    # Binomial: count how many pp > 0.5; under right-skew, most should be < 0.5
    n_above_half = np.sum(pp_vals > 0.5)
    # Binomial test: is proportion above 0.5 significantly high?
    p_flat = float(stats.binom_test(n_above_half, n_sig, 0.5, alternative='greater')
                   if hasattr(stats, 'binom_test')
                   else stats.binomtest(n_above_half, n_sig, 0.5, alternative='greater').pvalue)

    return {
        'z_right': z_right,
        'p_right': p_right,
        'z_flat': 0,  # simplified
        'p_flat': float(p_flat),
        'evidential': p_right < 0.05,
        'inadequate': p_flat < 0.05,
        'n_sig': n_sig
    }


# ─── 4. Z-curve ───

def z_curve(yi, sei):
    """Z-curve analysis: observed vs expected discovery rate.

    Computes the observed inclusion rate (proportion significant) and
    estimates the expected discovery rate by modeling the z-value distribution.
    """
    k = len(yi)
    if k < 3:
        return {'oir': 0, 'edr': 0, 'ratio': 1.0}

    z_vals = np.abs(yi / sei)
    n_sig = np.sum(z_vals > 1.96)
    oir = float(n_sig / k)

    # Expected discovery rate: estimate the mean power of significant studies
    # Simple approach: use the mean z-value of significant studies to estimate power
    sig_z = z_vals[z_vals > 1.96].copy()
    if len(sig_z) == 0:
        return {'oir': oir, 'edr': 0, 'ratio': 0 if oir > 0 else 1.0}

    # Average power of significant studies (using observed z to estimate non-centrality)
    mean_z_sig = float(np.mean(sig_z))
    # Power at mean noncentrality: P(Z > 1.96 | delta = mean_z_sig)
    edr_est = float(stats.norm.sf(1.96 - mean_z_sig))

    # Adjusted EDR for the full set
    edr = edr_est * (n_sig / k) if k > 0 else 0

    ratio = edr / oir if oir > 0 else 1.0

    return {
        'oir': oir,
        'edr': float(edr),
        'ratio': float(min(ratio, 2.0)),
        'mean_z_sig': mean_z_sig,
        'n_sig': int(n_sig)
    }


# ─── 5. Trim-and-Fill ───

def trim_and_fill(yi, sei):
    """Duval-Tweedie trim-and-fill with auto side detection (L0 estimator)."""
    k = len(yi)
    if k < 3:
        return {'theta_adj': float(np.mean(yi)), 'k0': 0, 'theta_unadj': float(np.mean(yi))}

    # Unadjusted pooled (DL)
    wi = 1.0 / sei ** 2
    sum_w = np.sum(wi)
    theta_fe = float(np.sum(wi * yi) / sum_w)
    Q = float(np.sum(wi * (yi - theta_fe) ** 2))
    C = float(sum_w - np.sum(wi ** 2) / sum_w)
    tau2 = max(0, (Q - (k - 1)) / C) if C > 0 else 0
    wi_star = 1.0 / (sei ** 2 + tau2)
    theta0 = float(np.sum(wi_star * yi) / np.sum(wi_star))

    # Side detection
    di = yi - theta0
    n_right = np.sum(di > 0)
    n_left = np.sum(di < 0)
    side = 'right' if n_right > n_left else 'left'

    # L0 estimator
    ranks = np.argsort(np.abs(di))
    signed_ranks = np.zeros(k)
    for i, r in enumerate(ranks):
        signed_ranks[r] = (i + 1) * np.sign(di[r])

    if side == 'right':
        S = np.sum(signed_ranks[signed_ranks > 0])
    else:
        S = np.sum(np.abs(signed_ranks[signed_ranks < 0]))

    k0 = max(0, round((4 * S - k * (k + 1)) / (2 * k + 1)))

    if k0 == 0:
        return {'theta_adj': theta0, 'k0': 0, 'theta_unadj': theta0}

    # Impute
    order = np.argsort(di)
    if side == 'right':
        extreme_idx = order[-k0:].copy()
    else:
        extreme_idx = order[:k0].copy()

    yi_filled = np.concatenate([yi, 2 * theta0 - yi[extreme_idx]])
    sei_filled = np.concatenate([sei, sei[extreme_idx]])

    # Re-pool
    wi2 = 1.0 / sei_filled ** 2
    theta_fe2 = np.sum(wi2 * yi_filled) / np.sum(wi2)
    Q2 = np.sum(wi2 * (yi_filled - theta_fe2) ** 2)
    C2 = np.sum(wi2) - np.sum(wi2 ** 2) / np.sum(wi2)
    tau2_2 = max(0, (Q2 - len(yi_filled) + 1) / C2) if C2 > 0 else 0
    wi_star2 = 1.0 / (sei_filled ** 2 + tau2_2)
    theta_adj = float(np.sum(wi_star2 * yi_filled) / np.sum(wi_star2))

    return {'theta_adj': theta_adj, 'k0': int(k0), 'theta_unadj': theta0}


# ─── 6. PET-PEESE ───

def pet_peese(yi, sei):
    """Conditional PET-PEESE with t-distribution."""
    k = len(yi)
    if k < 3:
        return {'theta_adj': float(np.mean(yi)), 'method_used': 'none'}

    wi = 1.0 / sei ** 2
    df = max(1, k - 2)

    # PET: WLS of yi on sei
    pet_int, pet_se, pet_p = _wls_regression(yi, sei, wi, df, use_sq=False)

    if pet_p >= 0.05:
        intercept, se_int, p_val = pet_int, pet_se, pet_p
        method = 'PET'
    else:
        intercept, se_int, p_val = _wls_regression(yi, sei, wi, df, use_sq=True)
        method = 'PEESE'

    return {'theta_adj': float(intercept), 'method_used': method, 'p': float(p_val)}


def _wls_regression(yi, sei, wi, df, use_sq=False):
    x = sei ** 2 if use_sq else sei
    sw = np.sum(wi)
    sx = np.sum(wi * x)
    sy = np.sum(wi * yi)
    sxx = np.sum(wi * x ** 2)
    sxy = np.sum(wi * x * yi)
    denom = sw * sxx - sx ** 2
    scale_ref = max(abs(sw * sxx), abs(sx ** 2), 1e-30)

    if abs(denom) < 1e-10 * scale_ref:
        intercept = float(sy / sw) if sw > 0 else 0
        se = float(1.0 / math.sqrt(sw)) if sw > 0 else 1
        t_stat = intercept / se if se > 0 else 0
        p = 2 * (1 - stats.t.cdf(abs(t_stat), df))
        return intercept, se, p

    intercept = float((sxx * sy - sx * sxy) / denom)
    slope = float((sw * sxy - sx * sy) / denom)
    residuals = yi - intercept - slope * x
    sigma2 = float(np.sum(wi * residuals ** 2) / max(len(yi) - 2, 1))
    var_int = sigma2 * sxx / denom
    se = math.sqrt(max(0, float(var_int)))
    t_stat = intercept / se if se > 0 else 0
    p = 2 * (1 - stats.t.cdf(abs(t_stat), df))
    return intercept, se, p


# ─── 7. 3-Parameter Selection Model (Vevea-Hedges) ───

def selection_model_3psm(yi, sei):
    """Simplified 3-parameter selection model.

    Assumes step-function weight: w(p) = 1 if p < 0.025 (one-sided), else eta.
    Estimates theta, tau2, and eta via profile likelihood.
    """
    k = len(yi)
    if k < 5:
        return {'theta_adj': float(np.mean(yi)), 'eta': 1.0, 'lr_p': 1.0, 'significant': False}

    vi = sei ** 2
    z_vals = yi / sei
    p_onesided = stats.norm.sf(z_vals)  # one-sided p-values

    # Unadjusted DL
    wi = 1.0 / vi
    theta_fe = np.sum(wi * yi) / np.sum(wi)
    Q = np.sum(wi * (yi - theta_fe) ** 2)
    C = np.sum(wi) - np.sum(wi ** 2) / np.sum(wi)
    tau2_unadj = max(0, (Q - (k - 1)) / C) if C > 0 else 0

    # Profile likelihood over eta (selection weight for non-significant studies)
    best_ll = -np.inf
    best_eta = 1.0
    best_theta = theta_fe
    best_tau2 = tau2_unadj

    # P0-1 FIX: proper weighted likelihood WITH normalizing constant
    # L_i = w(p_i) * f(y_i|theta,tau2) / C_i
    # where C_i = integral w(p(y)) * f(y|theta,tau2) dy
    # For step-function at p=0.025: C_i = Phi(z_crit - theta/se_i) + eta*(1 - Phi(z_crit - theta/se_i))
    # where z_crit = 1.96 (one-sided 0.025)
    z_crit = stats.norm.ppf(0.975)  # 1.96

    for eta_try in [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]:
        wt = np.where(p_onesided < 0.025, 1.0, eta_try)

        # Estimate theta for this eta (weighted DL)
        wi_w = wt / (vi + tau2_unadj)
        if np.sum(wi_w) < 1e-10:
            continue
        theta_w = float(np.sum(wi_w * yi) / np.sum(wi_w))

        # Log-likelihood with normalizing constant
        ll = 0
        for i in range(k):
            var_i = vi[i].copy() + tau2_unadj
            se_i = math.sqrt(var_i)
            # f(y_i | theta, tau2): normal density
            log_f = -0.5 * math.log(2 * math.pi * var_i) - 0.5 * (yi[i] - theta_w) ** 2 / var_i
            # w(p_i): selection weight
            log_w = math.log(max(wt[i], 1e-30))
            # C_i: normalizing constant
            # P(significant | theta, se_i) = P(|Z| > z_crit) = 1 - Phi(z_crit - theta_w/se_i) + Phi(-z_crit - theta_w/se_i)
            p_sig = (stats.norm.sf(z_crit - theta_w / se_i) +
                     stats.norm.cdf(-z_crit - theta_w / se_i))
            p_nonsig = 1 - p_sig
            C_i = p_sig * 1.0 + p_nonsig * eta_try
            C_i = max(C_i, 1e-30)

            ll += log_f + log_w - math.log(C_i)

        if ll > best_ll:
            best_ll = ll
            best_eta = eta_try
            best_theta = theta_w
            best_tau2 = tau2_unadj

    # LR test: compare eta=1 (no selection) vs best eta
    # Null model: eta=1, same likelihood computation
    ll_null = 0
    wi_null = 1.0 / (vi + tau2_unadj)
    theta_null = float(np.sum(wi_null * yi) / np.sum(wi_null))
    for i in range(k):
        var_i = vi[i].copy() + tau2_unadj
        ll_null += -0.5 * math.log(2 * math.pi * var_i) - 0.5 * (yi[i] - theta_null) ** 2 / var_i
        # Under eta=1, C_i = 1 and log(w_i) = 0, so just the normal density

    lr_stat = 2 * (best_ll - ll_null)
    lr_p = 1 - stats.chi2.cdf(max(0, lr_stat), 1)

    return {
        'theta_adj': float(best_theta),
        'eta': float(best_eta),
        'lr_stat': float(lr_stat),
        'lr_p': float(lr_p),
        'significant': lr_p < 0.10
    }


# ─── 8. Limit Meta-Analysis (Rucker) ───

def limit_meta(yi, sei):
    """Rucker's limit meta-analysis.

    Regress effect size on SE, extrapolate to SE=0.
    The intercept is the "bias-free" estimate.
    """
    k = len(yi)
    if k < 3:
        return {'theta_limit': float(np.mean(yi)), 'se_limit': 0, 'p_limit': 1.0}

    # WLS regression of yi on sei, weighted by 1/sei^2
    wi = 1.0 / sei ** 2
    df = max(1, k - 2)
    intercept, se_int, p_val = _wls_regression(yi, sei, wi, df, use_sq=False)

    return {
        'theta_limit': float(intercept),
        'se_limit': float(se_int),
        'p_limit': float(p_val)
    }


# ─── Master forensics function ───

def run_all_methods(yi, sei):
    """Run all 8 publication bias methods on a single review."""
    return {
        'egger': egger_test(yi, sei),
        'begg': begg_test(yi, sei),
        'pcurve': p_curve(yi, sei),
        'zcurve': z_curve(yi, sei),
        'trimfill': trim_and_fill(yi, sei),
        'petpeese': pet_peese(yi, sei),
        'sel3psm': selection_model_3psm(yi, sei),
        'limit': limit_meta(yi, sei),
    }
