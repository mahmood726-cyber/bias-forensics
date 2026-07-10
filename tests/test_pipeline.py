"""Tests for src.pipeline.classify_bias (F2: previously zero coverage).

Builds synthetic fingerprint dicts matching the output shape of
src.methods.run_all_methods and asserts classify_bias returns the expected
bias_class for representative Clean / Confirmed / Discordant cases, plus the
theta_unadj~0 branch that switches shift from relative to absolute
(guards against division by ~0).
"""

import sys

sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent.parent))
from src.pipeline import classify_bias


def make_fingerprint(*, egger_sig=False, begg_sig=False, pcurve_inadequate=False,
                     sel3psm_sig=False, theta_unadj=0.5, corrections=None):
    """Build a minimal fingerprint dict of the shape classify_bias consumes.

    corrections = [trimfill.theta_adj, petpeese.theta_adj,
                   sel3psm.theta_adj, limit.theta_limit]
    """
    if corrections is None:
        corrections = [theta_unadj] * 4
    tf_adj, pp_adj, sel_adj, limit_adj = corrections
    return {
        'egger': {'significant': egger_sig},
        'begg': {'significant': begg_sig},
        'pcurve': {'inadequate': pcurve_inadequate},
        'sel3psm': {'significant': sel3psm_sig, 'theta_adj': sel_adj},
        'trimfill': {'theta_unadj': theta_unadj, 'theta_adj': tf_adj},
        'petpeese': {'theta_adj': pp_adj},
        'limit': {'theta_limit': limit_adj},
    }


def test_clean_no_evidence_no_shift():
    """No detection tests fire and corrections equal the unadjusted effect."""
    fp = make_fingerprint(theta_unadj=0.5, corrections=[0.5, 0.5, 0.5, 0.5])
    out = classify_bias(fp)
    assert out['bias_class'] == 'Clean'
    assert out['n_detect'] == 0
    assert out['mean_shift'] == 0.0
    assert out['concordance'] == 100.0


def test_confirmed_many_detections_and_large_shift():
    """>=3 detection tests significant AND mean shift >= 0.2 -> Confirmed."""
    # theta_unadj=0.5, corrections all 0.3 -> relative shift = 0.2/0.5 = 0.4 each.
    fp = make_fingerprint(egger_sig=True, begg_sig=True, sel3psm_sig=True,
                          theta_unadj=0.5, corrections=[0.3, 0.3, 0.3, 0.3])
    out = classify_bias(fp)
    assert out['n_detect'] == 3
    assert out['mean_shift'] >= 0.2
    assert out['bias_class'] == 'Confirmed'


def test_discordant_when_corrections_disagree_on_direction():
    """Corrections split in sign (n_agree_dir <= 2) -> Discordant."""
    # unadj +0.5; corrections split two positive / two negative.
    fp = make_fingerprint(theta_unadj=0.5, corrections=[0.5, 0.5, -0.5, -0.5])
    out = classify_bias(fp)
    assert out['bias_class'] == 'Discordant'
    assert out['concordance'] == 50.0


def test_theta_unadj_near_zero_uses_absolute_shift():
    """abs(theta_unadj) <= 0.01 switches to absolute shift (no div-by-~0)."""
    # theta_unadj tiny; corrections offset by 0.1 in absolute terms.
    fp = make_fingerprint(theta_unadj=0.005,
                          corrections=[0.105, 0.105, 0.105, 0.105])
    out = classify_bias(fp)
    # absolute shift = |0.105 - 0.005| = 0.1 for each correction
    assert abs(out['mean_shift'] - 0.1) < 1e-9
    assert out['max_shift'] == 0.1
    assert out['bias_class'] == 'Clean'


def test_theta_unadj_exactly_zero_no_crash():
    """theta_unadj == 0 must not raise (absolute-shift branch)."""
    fp = make_fingerprint(theta_unadj=0.0, corrections=[0.4, 0.4, 0.4, 0.4])
    out = classify_bias(fp)
    assert out['bias_class'] in ('Clean', 'Suspected', 'Confirmed', 'Discordant')
    assert abs(out['mean_shift'] - 0.4) < 1e-9


if __name__ == '__main__':
    import pytest
    pytest.main([__file__, '-v'])
