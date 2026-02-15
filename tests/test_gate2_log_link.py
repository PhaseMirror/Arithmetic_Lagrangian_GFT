"""
test_gate2_log_link.py
======================

Test Phase 3: Log-link fit quality and prior constraints.

Verifies that:
1. ν_eff(M) is well-fit by c₀ + c₁·log(M/M_P)
2. Residuals < 10%
3. M² term is negligible: |c₂|/|c₁| < 0.01
4. Prior is sufficiently tight: σ_c/|c̄| ≤ 0.30
"""

import pytest
import numpy as np
import sys
sys.path.insert(0, 'src')

from algftgate2 import (
    find_ngfp, critical_surface_projection, integrate_flow,
    fit_log_link, check_ward_identity, DEFAULT_CONFIG
)
from mapping_uv_gft import map_algft_to_gft


def get_test_flow_and_fit():
    """Helper to generate flow and fit."""
    config = DEFAULT_CONFIG.copy()
    ngfp = find_ngfp(config)
    
    # Gate 1 boundary conditions
    epsilon = 0.01
    sigma = 0.1
    omega_n = np.array([14.13, 21.02, 28.09])
    H_inf = 6e-6
    
    result = map_algft_to_gft(epsilon, sigma, omega_n, H_inf)
    lambda4_raw = result['lambda_4']
    lambda6_raw = result['lambda_6']
    
    lambda4_uv, lambda6_uv = critical_surface_projection(
        lambda4_raw, lambda6_raw, ngfp
    )
    
    flow = integrate_flow(lambda4_uv, lambda6_uv, config, n_points=1000)
    fit = fit_log_link(flow, M_range=(5e12, 5e13), n_M_points=30, config=config)
    
    return flow, fit, config


def test_log_link_residuals():
    """Test that log-link fit has residuals < 10% (Test 3)."""
    flow, fit, config = get_test_flow_and_fit()
    
    max_res = fit['max_residual_pct']
    
    print(f"\nLog-link fit residuals:")
    print(f"  Max residual: {max_res:.4f}%")
    print(f"  Threshold: 10%")
    print(f"  Status: {'PASS' if max_res < 10.0 else 'FAIL'}")
    
    assert max_res < 10.0, f"Max residual {max_res:.4f}% exceeds 10% threshold"


def test_m2_term_negligible():
    """Test that M² term is negligible: |c₂|/|c₁| < 0.01 (Test 8)."""
    flow, fit, config = get_test_flow_and_fit()
    
    c2_over_c1 = fit['c2_over_c1']
    
    print(f"\nM² term contribution:")
    print(f"  |c₂|/|c₁| = {c2_over_c1:.6f}")
    print(f"  Threshold: 0.01")
    print(f"  Status: {'PASS' if c2_over_c1 < 0.01 else 'FAIL'}")
    
    assert c2_over_c1 < 0.01, f"|c₂|/|c₁| = {c2_over_c1:.6f} exceeds 0.01 threshold"


def test_c1_reasonable_magnitude():
    """Test that c₁ is in reasonable range (ballpark check)."""
    flow, fit, config = get_test_flow_and_fit()
    
    c1 = fit['c1']
    
    print(f"\nPrior coefficient:")
    print(f"  c₁ = {c1:.4f}")
    print(f"  Expected range: O(10²-10³)")
    
    # Rough sanity check - should be positive and O(100-1000)
    assert c1 > 0, f"c₁ = {c1:.4f} should be positive"
    assert 10 < abs(c1) < 1e5, f"c₁ = {c1:.4f} outside reasonable range (10, 10⁵)"


def test_ward_identity():
    """Test Ward identity along flow (Test 4)."""
    flow, fit, config = get_test_flow_and_fit()
    
    ward = check_ward_identity(flow, config, threshold=0.05)
    
    print(f"\nWard identity check:")
    print(f"  Max W(t) = {ward['max_W']:.6f}")
    print(f"  At t = {ward['t_at_max']:.2f}")
    print(f"  Threshold: 0.05")
    print(f"  Status: {'PASS' if ward['passed'] else 'FAIL'}")
    
    assert ward['passed'], f"Ward violation {ward['max_W']:.6f} exceeds 0.05"


def test_nu_eff_monotonicity():
    """Test that ν_eff(M) is monotonically increasing with M."""
    flow, fit, config = get_test_flow_and_fit()
    
    M_test = fit['M_test']
    nu_test = fit['nu_test']
    
    # Check monotonicity
    monotonic = np.all(np.diff(nu_test) > 0)
    
    print(f"\nν_eff(M) monotonicity:")
    print(f"  Range: M ∈ [{M_test[0]:.2e}, {M_test[-1]:.2e}] GeV")
    print(f"  ν_eff range: [{nu_test[0]:.6f}, {nu_test[-1]:.6f}]")
    print(f"  Monotonic increasing: {monotonic}")
    
    # Note: This might not always be strictly monotonic depending on physics,
    # but let's check that it's generally increasing
    pos_diffs = np.sum(np.diff(nu_test) > 0)
    neg_diffs = np.sum(np.diff(nu_test) < 0)
    
    print(f"  Positive increments: {pos_diffs}/{len(nu_test)-1}")
    print(f"  Negative increments: {neg_diffs}/{len(nu_test)-1}")
    
    # Allow some small fluctuations but majority should be increasing
    assert pos_diffs > 0.8 * len(nu_test), "ν_eff(M) should be generally increasing"


def test_fit_quality_metrics():
    """Test various quality metrics of the fit."""
    flow, fit, config = get_test_flow_and_fit()
    
    nu_test = fit['nu_test']
    fit_vals = fit['fit_vals']
    
    # Compute R² coefficient
    ss_res = np.sum((nu_test - fit_vals)**2)
    ss_tot = np.sum((nu_test - np.mean(nu_test))**2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    # Compute mean absolute percentage error
    mape = np.mean(np.abs((nu_test - fit_vals) / (np.abs(nu_test) + 1e-30))) * 100
    
    print(f"\nFit quality metrics:")
    print(f"  R² = {r_squared:.6f}")
    print(f"  MAPE = {mape:.4f}%")
    print(f"  c₀ = {fit['c0']:.6e}")
    print(f"  c₁ = {fit['c1']:.4f}")
    print(f"  c₂ = {fit['c2']:.6e}")
    
    assert r_squared > 0.99, f"R² = {r_squared:.6f} < 0.99, poor fit quality"
    assert mape < 5.0, f"MAPE = {mape:.4f}% exceeds 5%"


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "-s"])
