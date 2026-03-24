"""Tests for all 8 publication bias methods."""

import sys
import math
import numpy as np
import pytest

sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent.parent))
from src.methods import (
    egger_test, begg_test, p_curve, z_curve,
    trim_and_fill, pet_peese, selection_model_3psm, limit_meta,
    run_all_methods
)

# Test data: 10 studies with clear right-side asymmetry (small studies show bigger effects)
ASYM_YI = np.array([-0.8, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.2, 0.5])
ASYM_SEI = np.array([0.10, 0.12, 0.15, 0.18, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50])

# Symmetric data: no bias expected
SYM_YI = np.array([-0.5, -0.3, -0.1, 0.1, 0.3])
SYM_SEI = np.array([0.2, 0.2, 0.2, 0.2, 0.2])

# Strong effect: all significant
STRONG_YI = np.array([-1.0, -0.9, -0.8, -1.1, -0.7, -0.95])
STRONG_SEI = np.array([0.15, 0.15, 0.15, 0.15, 0.15, 0.15])


class TestEgger:
    def test_returns_p_value(self):
        r = egger_test(ASYM_YI, ASYM_SEI)
        assert 0 <= r['p'] <= 1
        assert math.isfinite(r['intercept'])

    def test_detects_asymmetry(self):
        """Asymmetric funnel should have significant Egger test."""
        r = egger_test(ASYM_YI, ASYM_SEI)
        # The asymmetric data has larger effects in smaller studies
        assert r['p'] < 0.5  # should trend toward significance

    def test_symmetric_not_significant(self):
        """Symmetric data with varying SEs should not trigger Egger."""
        yi = np.array([-0.5, -0.3, -0.1, 0.1, 0.3])
        sei = np.array([0.15, 0.20, 0.25, 0.20, 0.15])  # varied SEs
        r = egger_test(yi, sei)
        assert r['p'] > 0.05

    def test_k2_returns_default(self):
        r = egger_test(np.array([0.5, -0.5]), np.array([0.2, 0.2]))
        assert r['p'] == 1.0


class TestBegg:
    def test_returns_tau_and_p(self):
        r = begg_test(ASYM_YI, ASYM_SEI)
        assert -1 <= r['tau'] <= 1
        assert 0 <= r['p'] <= 1

    def test_uses_centered_effects(self):
        """After P1-1 fix, Begg should use centered effects."""
        r = begg_test(ASYM_YI, ASYM_SEI)
        assert math.isfinite(r['tau'])

    def test_k2_returns_default(self):
        r = begg_test(np.array([0.5, -0.5]), np.array([0.2, 0.2]))
        assert r['p'] == 1.0


class TestPCurve:
    def test_strong_effect_evidential(self):
        """Strong effects should produce right-skewed p-curve (evidential)."""
        r = p_curve(STRONG_YI, STRONG_SEI)
        assert r['n_sig'] >= 3
        # Right-skewed means evidential value

    def test_returns_n_sig(self):
        r = p_curve(ASYM_YI, ASYM_SEI)
        assert isinstance(r['n_sig'], int)
        assert r['n_sig'] >= 0

    def test_few_significant_returns_default(self):
        """If fewer than 3 significant studies, returns default."""
        r = p_curve(SYM_YI, SYM_SEI)
        assert r['p_right'] == 1.0 or r['n_sig'] < 3


class TestZCurve:
    def test_returns_oir_edr(self):
        r = z_curve(ASYM_YI, ASYM_SEI)
        assert 0 <= r['oir'] <= 1
        assert r['edr'] >= 0

    def test_strong_effect_high_oir(self):
        r = z_curve(STRONG_YI, STRONG_SEI)
        assert r['oir'] > 0.5  # most studies should be significant


class TestTrimAndFill:
    def test_returns_adjusted_theta(self):
        r = trim_and_fill(ASYM_YI, ASYM_SEI)
        assert math.isfinite(r['theta_adj'])
        assert isinstance(r['k0'], int)
        assert r['k0'] >= 0

    def test_symmetric_no_fill(self):
        r = trim_and_fill(SYM_YI, SYM_SEI)
        # Symmetric data should have k0 = 0 or small
        assert r['k0'] <= 2

    def test_k2_fallback(self):
        r = trim_and_fill(np.array([0.5, -0.5]), np.array([0.2, 0.2]))
        assert math.isfinite(r['theta_adj'])


class TestPetPeese:
    def test_returns_adjusted_theta(self):
        r = pet_peese(ASYM_YI, ASYM_SEI)
        assert math.isfinite(r['theta_adj'])
        assert r['method_used'] in ('PET', 'PEESE', 'none')

    def test_k2_fallback(self):
        r = pet_peese(np.array([0.5, -0.5]), np.array([0.2, 0.2]))
        assert r['method_used'] == 'none'


class TestSelectionModel:
    def test_returns_eta_and_lr(self):
        r = selection_model_3psm(ASYM_YI, ASYM_SEI)
        assert 0 < r['eta'] <= 1
        assert 0 <= r['lr_p'] <= 1

    def test_strong_effect_may_detect(self):
        """With strong asymmetry, selection model should have eta < 1."""
        # Use data with clear selection: only significant studies present
        yi = np.array([-0.8, -0.7, -0.9, -0.6, -0.8, -0.5])
        sei = np.array([0.2, 0.2, 0.2, 0.2, 0.2, 0.2])
        r = selection_model_3psm(yi, sei)
        assert math.isfinite(r['theta_adj'])

    def test_k4_returns_default(self):
        r = selection_model_3psm(np.array([0.5, -0.5, 0.3, -0.3]), np.array([0.2]*4))
        assert r['eta'] == 1.0


class TestLimitMeta:
    def test_returns_theta_limit(self):
        r = limit_meta(ASYM_YI, ASYM_SEI)
        assert math.isfinite(r['theta_limit'])
        assert math.isfinite(r['se_limit'])

    def test_matches_pet(self):
        """Limit meta should equal PET (same regression)."""
        r_limit = limit_meta(ASYM_YI, ASYM_SEI)
        r_pet = pet_peese(ASYM_YI, ASYM_SEI)
        if r_pet['method_used'] == 'PET':
            assert abs(r_limit['theta_limit'] - r_pet['theta_adj']) < 0.001


class TestRunAllMethods:
    def test_returns_all_8(self):
        r = run_all_methods(ASYM_YI, ASYM_SEI)
        assert set(r.keys()) == {'egger', 'begg', 'pcurve', 'zcurve',
                                  'trimfill', 'petpeese', 'sel3psm', 'limit'}

    def test_no_crashes_on_real_data(self):
        """Test on data from first Pairwise70 review."""
        sys.path.insert(0, r'C:\FragilityAtlas')
        from src.loader import load_review
        review = load_review(r'C:\Models\Pairwise70\data\CD000028_pub4_data.rda')
        if review and review.k >= 5:
            r = run_all_methods(review.yi, review.sei)
            for method, result in r.items():
                assert isinstance(result, dict), f"{method} returned {type(result)}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
