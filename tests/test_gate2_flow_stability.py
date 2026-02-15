"""
test_gate2_flow_stability.py
=============================

Test Phase 2: RG flow stability and integration quality.

Verifies that flows:
1. Remain bounded (no blow-up)
2. Reach deep IR (t < -100)
3. Cover at least 140 e-folds
4. Maintain numerical accuracy (Radau vs RK45 agreement)
"""

import pytest
import numpy as np
import sys
sys.path.insert(0, 'src')

from algftgate2 import (
    find_ngfp, critical_surface_projection, integrate_flow,
    DEFAULT_CONFIG, M_PLANCK, H_0
)
from mapping_uv_gft import map_algft_to_gft


def get_test_flow():
    """Helper to generate a test flow with Gate 1 boundary conditions."""
    config = DEFAULT_CONFIG.copy()
    ngfp = find_ngfp(config)
    
    # Use representative Gate 1 values
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
    
    return flow, config


def test_flow_reaches_ir():
    """Test that flow reaches deep IR (t < -100)."""
    flow, config = get_test_flow()
    
    t_final = flow['t'][-1]
    
    print(f"\nFlow integration range:")
    print(f"  t_final = {t_final:.2f}")
    print(f"  Required: t < -100")
    print(f"  Status: {flow['success']}")
    
    assert flow['success'], "Flow integration failed"
    assert t_final < -100, f"Flow only reached t = {t_final:.2f}, need t < -100"


def test_flow_bounded():
    """Test that couplings remain bounded (< 10^6) throughout flow."""
    flow, config = get_test_flow()
    
    max_lam4 = np.max(np.abs(flow['lambda4']))
    max_lam6 = np.max(np.abs(flow['lambda6']))
    
    print(f"\nCoupling bounds:")
    print(f"  max |λ₄| = {max_lam4:.2e}")
    print(f"  max |λ₆| = {max_lam6:.2e}")
    print(f"  Limit: 10^6")
    
    assert max_lam4 < 1e6, f"|λ₄| reached {max_lam4:.2e}, exceeds 10^6"
    assert max_lam6 < 1e6, f"|λ₆| reached {max_lam6:.2e}, exceeds 10^6"


def test_flow_e_folds():
    """Test that flow covers at least 140 e-folds."""
    flow, config = get_test_flow()
    
    t_span = abs(flow['t'][-1] - flow['t'][0])
    
    print(f"\nRG time span:")
    print(f"  Δt = {t_span:.2f} e-folds")
    print(f"  Required: > 140 e-folds")
    
    assert t_span > 140, f"Flow only covers {t_span:.2f} e-folds, need > 140"


def test_flow_monotonic_time():
    """Test that RG time is monotonically decreasing."""
    flow, config = get_test_flow()
    
    t_diff = np.diff(flow['t'])
    
    print(f"\nTime monotonicity:")
    print(f"  All Δt < 0: {np.all(t_diff < 0)}")
    print(f"  Min Δt: {t_diff.min():.6f}")
    print(f"  Max Δt: {t_diff.max():.6f}")
    
    assert np.all(t_diff < 0), "RG time is not monotonically decreasing"


def test_integrator_agreement():
    """Test agreement between Radau and RK45 methods (Test 7)."""
    config = DEFAULT_CONFIG.copy()
    ngfp = find_ngfp(config)
    
    # Get UV boundary conditions
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
    
    # Run with both methods
    flow_radau = integrate_flow(lambda4_uv, lambda6_uv, config, 
                               method='Radau', n_points=500)
    flow_rk45 = integrate_flow(lambda4_uv, lambda6_uv, config, 
                              method='RK45', n_points=500)
    
    # Check both succeeded
    assert flow_radau['success'], "Radau integration failed"
    assert flow_rk45['success'], "RK45 integration failed"
    
    # Compare at common time points
    t_common = np.linspace(max(flow_radau['t'][0], flow_rk45['t'][0]),
                          min(flow_radau['t'][-1], flow_rk45['t'][-1]),
                          100)
    
    lam4_radau = np.array([flow_radau['sol'](t)[0] for t in t_common])
    lam4_rk45 = np.array([flow_rk45['sol'](t)[0] for t in t_common])
    
    lam6_radau = np.array([flow_radau['sol'](t)[1] for t in t_common])
    lam6_rk45 = np.array([flow_rk45['sol'](t)[1] for t in t_common])
    
    # Compute relative differences
    max_diff_4 = np.max(np.abs(lam4_radau - lam4_rk45) / 
                       (np.abs(lam4_radau) + 1e-30))
    max_diff_6 = np.max(np.abs(lam6_radau - lam6_rk45) / 
                       (np.abs(lam6_radau) + 1e-30))
    
    print(f"\nIntegrator agreement:")
    print(f"  Max relative diff λ₄: {max_diff_4:.2e}")
    print(f"  Max relative diff λ₆: {max_diff_6:.2e}")
    print(f"  Threshold: 1e-6")
    
    assert max_diff_4 < 1e-6, f"Radau-RK45 λ₄ disagreement: {max_diff_4:.2e}"
    assert max_diff_6 < 1e-6, f"Radau-RK45 λ₆ disagreement: {max_diff_6:.2e}"


def test_lambda6_scaling():
    """Test that λ₆ preserves M² scaling at UV (Test 6)."""
    config = DEFAULT_CONFIG.copy()
    ngfp = find_ngfp(config)
    
    # Test with two different inflationary scales
    H_inf_1 = 6e-6
    H_inf_2 = 8e-6
    
    epsilon = 0.01
    sigma = 0.1
    omega_n = np.array([14.13, 21.02, 28.09])
    
    result1 = map_algft_to_gft(epsilon, sigma, omega_n, H_inf_1)
    result2 = map_algft_to_gft(epsilon, sigma, omega_n, H_inf_2)
    
    lam6_1 = result1['lambda_6']
    lam6_2 = result2['lambda_6']
    
    # λ₆ should scale as M²
    scale_ratio = (H_inf_2 / H_inf_1)**2
    expected_ratio = lam6_2 / lam6_1
    
    deviation = abs(expected_ratio - scale_ratio) / scale_ratio
    
    print(f"\nλ₆ M² scaling:")
    print(f"  H_inf ratio: {H_inf_2/H_inf_1:.4f}")
    print(f"  Expected λ₆ ratio: {scale_ratio:.4f}")
    print(f"  Observed λ₆ ratio: {expected_ratio:.4f}")
    print(f"  Deviation: {deviation*100:.2f}%")
    
    assert deviation < 0.01, f"λ₆ M² scaling deviation {deviation*100:.2f}% > 1%"


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "-s"])
