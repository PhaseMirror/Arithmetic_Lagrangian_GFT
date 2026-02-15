"""
test_eps_zero_limit.py
======================

Phase 2 Test: Verify ε → 0 recovers ΛCDM spectrum.

Pass criterion: |M(k) - 1| < 10⁻⁶ for all k when ε = 0

Author: PhaseMirror/Arithmetic_Lagrangian_GFT
Date: February 2026
"""

import numpy as np
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from algftgate1 import zeta_comb_modulation, epsilon_zero_limit_test
from algft_sk import compute_modulation, verify_lambda_cdm_limit


def test_baseline_eps_zero(
    tolerance: float = 1e-6,
    verbose: bool = True
) -> bool:
    """
    Test that baseline implementation recovers ΛCDM for ε = 0.
    
    Parameters
    ----------
    tolerance : float
        Maximum allowed |M(k) - 1|
    verbose : bool
        Print detailed results
        
    Returns
    -------
    passes : bool
        True if test passes
    """
    k_array = np.logspace(-4, 0, 100)
    
    # Use the built-in test from algftgate1
    passes, max_dev = epsilon_zero_limit_test(k_array, tolerance)
    
    if verbose:
        print("Test: Baseline ε → 0 (ΛCDM recovery)")
        print("-" * 70)
        print(f"k range: [{k_array.min():.2e}, {k_array.max():.2e}] Mpc⁻¹")
        print(f"Number of points: {len(k_array)}")
        print(f"Max |M(k) - 1|: {max_dev:.2e}")
        print(f"Required: < {tolerance:.2e}")
        
        if passes:
            print("✓ Baseline correctly recovers ΛCDM")
        else:
            print("✗ Baseline fails ΛCDM recovery")
        print()
    
    return passes


def test_sk_eps_zero(
    tolerance: float = 1e-6,
    verbose: bool = True
) -> bool:
    """
    Test that SK implementation recovers ΛCDM for ε = 0.
    
    Parameters
    ----------
    tolerance : float
        Maximum allowed |M(k) - 1|
    verbose : bool
        Print detailed results
        
    Returns
    -------
    passes : bool
        True if test passes
    """
    k_array = np.logspace(-4, 0, 100)
    omega_n = np.array([14.13, 21.02, 28.09])
    phi_n = np.array([0.0, np.pi/4, np.pi/3])
    
    # Use the built-in test from algft_sk
    passes, max_dev = verify_lambda_cdm_limit(k_array, omega_n, phi_n, tolerance)
    
    if verbose:
        print("Test: SK ε → 0 (ΛCDM recovery)")
        print("-" * 70)
        print(f"k range: [{k_array.min():.2e}, {k_array.max():.2e}] Mpc⁻¹")
        print(f"Number of points: {len(k_array)}")
        print(f"Max |M(k) - 1|: {max_dev:.2e}")
        print(f"Required: < {tolerance:.2e}")
        
        if passes:
            print("✓ SK correctly recovers ΛCDM")
        else:
            print("✗ SK fails ΛCDM recovery")
        print()
    
    return passes


def test_small_epsilon_continuity(verbose: bool = True) -> bool:
    """
    Test that M(k) → 1 continuously as ε → 0.
    
    This verifies there are no discontinuities or numerical artifacts
    in the small-ε regime.
    
    Returns
    -------
    passes : bool
        True if test passes
    """
    k_array = np.logspace(-3, -1, 50)  # Intermediate scales
    omega_n = np.array([14.13, 21.02])
    phi_n = np.array([0.0, np.pi/4])
    sigma = 0.1
    k_star = 0.05
    
    # Test decreasing values of epsilon
    epsilon_values = np.array([1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6])
    max_deviations_baseline = []
    max_deviations_sk = []
    
    for eps in epsilon_values:
        # Baseline
        M_base = zeta_comb_modulation(k_array, eps, sigma, k_star, omega_n, phi_n)
        max_deviations_baseline.append(np.max(np.abs(M_base - 1)))
        
        # SK
        M_sk = compute_modulation(k_array, eps, sigma, omega_n, phi_n, k_star)
        max_deviations_sk.append(np.max(np.abs(M_sk - 1)))
    
    max_deviations_baseline = np.array(max_deviations_baseline)
    max_deviations_sk = np.array(max_deviations_sk)
    
    # Check that deviations decrease monotonically (approximately)
    # Allow for some numerical noise at very small ε
    baseline_monotonic = True
    sk_monotonic = True
    
    for i in range(len(epsilon_values) - 1):
        if max_deviations_baseline[i+1] > 2 * max_deviations_baseline[i]:
            baseline_monotonic = False
        if max_deviations_sk[i+1] > 2 * max_deviations_sk[i]:
            sk_monotonic = False
    
    # Check final values are small
    baseline_small = max_deviations_baseline[-1] < 1e-5
    sk_small = max_deviations_sk[-1] < 1e-5
    
    passes = baseline_monotonic and sk_monotonic and baseline_small and sk_small
    
    if verbose:
        print("Test: Continuity as ε → 0")
        print("-" * 70)
        print("ε values and max |M(k) - 1|:")
        print()
        print("  ε          Baseline      SK")
        print("  " + "-" * 35)
        for i, eps in enumerate(epsilon_values):
            print(f"  {eps:.0e}     {max_deviations_baseline[i]:.2e}    "
                  f"{max_deviations_sk[i]:.2e}")
        print()
        
        if passes:
            print("✓ Both implementations show smooth ε → 0 limit")
        else:
            print("✗ Non-monotonic or non-vanishing behavior detected")
        print()
    
    return passes


def test_power_spectrum_ratio(verbose: bool = True) -> bool:
    """
    Test that P_ζ(k) / P_0(k) = M(k) correctly.
    
    This verifies the relationship between the modulation and the
    power spectrum enhancement.
    
    Returns
    -------
    passes : bool
        True if test passes
    """
    from algftgate1 import power_spectrum_primordial, power_spectrum_algft
    
    k_array = np.logspace(-3, -1, 50)
    epsilon = 0.02  # Moderate value for clear signal
    sigma = 0.1
    k_star = 0.05
    omega_n = np.array([14.13, 21.02, 28.09])
    phi_n = np.array([0.0, np.pi/4, np.pi/3])
    
    # Compute power spectra
    P_0 = power_spectrum_primordial(k_array)
    P_zeta = power_spectrum_algft(k_array, epsilon, sigma, k_star, 
                                   omega_n, phi_n)
    
    # Compute ratio
    ratio = P_zeta / P_0
    
    # Compute M(k) from baseline
    M_k = zeta_comb_modulation(k_array, epsilon, sigma, k_star, omega_n, phi_n)
    
    # They should match exactly (within numerical precision)
    residuals = np.abs(ratio - M_k) / M_k
    max_residual = np.max(residuals)
    
    passes = max_residual < 1e-10
    
    if verbose:
        print("Test: P_ζ(k) / P_0(k) = M(k) identity")
        print("-" * 70)
        print(f"Max fractional difference: {max_residual:.2e}")
        print(f"Required: < 1e-10")
        
        if passes:
            print("✓ Power spectrum ratio correctly equals M(k)")
        else:
            print("✗ Mismatch between power spectrum ratio and M(k)")
        print()
    
    return passes


# ============================================================================
# Main test runner
# ============================================================================

if __name__ == "__main__":
    print()
    print("=" * 70)
    print("Phase 2: Testing ε → 0 Limit (ΛCDM Recovery)")
    print("=" * 70)
    print()
    
    # Run all tests
    test1_pass = test_baseline_eps_zero(verbose=True)
    test2_pass = test_sk_eps_zero(verbose=True)
    test3_pass = test_small_epsilon_continuity(verbose=True)
    test4_pass = test_power_spectrum_ratio(verbose=True)
    
    print("=" * 70)
    print("Test Summary")
    print("=" * 70)
    print(f"Baseline ε=0 test:        {'PASS ✓' if test1_pass else 'FAIL ✗'}")
    print(f"SK ε=0 test:              {'PASS ✓' if test2_pass else 'FAIL ✗'}")
    print(f"Continuity test:          {'PASS ✓' if test3_pass else 'FAIL ✗'}")
    print(f"P_ζ/P_0 = M(k) test:      {'PASS ✓' if test4_pass else 'FAIL ✗'}")
    print()
    
    all_pass = test1_pass and test2_pass and test3_pass and test4_pass
    
    if all_pass:
        print("✓✓✓ ALL TESTS PASSED ✓✓✓")
        print()
        print("Phase 2 Criterion Met: ε → 0 correctly recovers ΛCDM")
        sys.exit(0)
    else:
        print("✗✗✗ SOME TESTS FAILED ✗✗✗")
        sys.exit(1)
