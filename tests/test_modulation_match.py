"""
test_modulation_match.py
========================

Phase 2 Test: Verify SK-derived modulation matches baseline implementation.

This test compares:
    algft_sk.compute_modulation(k) vs algftgate1.zeta_comb_modulation(k)

Pass criterion: Max fractional residual < 2% over k ∈ [10⁻⁴, 1] Mpc⁻¹

Author: PhaseMirror/Arithmetic_Lagrangian_GFT
Date: February 2026
"""

import numpy as np
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from algftgate1 import zeta_comb_modulation
from algft_sk import compute_modulation


def test_modulation_match(
    k_min: float = 1e-4,
    k_max: float = 1.0,
    n_k: int = 100,
    max_residual_allowed: float = 0.02,
    verbose: bool = True
) -> bool:
    """
    Test that SK-derived M(k) matches baseline M(k) within 2%.
    
    Parameters
    ----------
    k_min : float
        Minimum wavenumber (Mpc⁻¹)
    k_max : float
        Maximum wavenumber (Mpc⁻¹)
    n_k : int
        Number of k points
    max_residual_allowed : float
        Maximum fractional residual (default 0.02 = 2%)
    verbose : bool
        Print detailed results
        
    Returns
    -------
    passes : bool
        True if test passes
    """
    # Test parameters (must match between implementations)
    epsilon = 1e-2
    sigma = 0.1
    k_star = 0.05
    omega_n = np.array([14.13, 21.02, 28.09, 35.44, 42.89])
    phi_n = np.array([0.0, np.pi/4, np.pi/3, -np.pi/6, np.pi/2])
    
    # Generate k array
    k_array = np.logspace(np.log10(k_min), np.log10(k_max), n_k)
    
    # Compute modulation from baseline
    M_baseline = zeta_comb_modulation(
        k_array, epsilon, sigma, k_star, omega_n, phi_n
    )
    
    # Compute modulation from SK derivation
    M_sk = compute_modulation(
        k_array, epsilon, sigma, omega_n, phi_n, k_star
    )
    
    # Compute fractional residuals
    # Use M_baseline as reference (avoid division by zero)
    residuals = np.abs(M_sk - M_baseline) / np.abs(M_baseline)
    
    max_residual = np.max(residuals)
    mean_residual = np.mean(residuals)
    
    # Test passes if max residual < 2%
    passes = max_residual < max_residual_allowed
    
    if verbose:
        print("=" * 70)
        print("Test: SK Modulation vs Baseline Modulation")
        print("=" * 70)
        print(f"k range: [{k_min:.2e}, {k_max:.2e}] Mpc⁻¹")
        print(f"Number of points: {n_k}")
        print(f"Parameters: ε={epsilon}, σ={sigma}, k*={k_star}")
        print()
        print("Results:")
        print("-" * 70)
        print(f"M_baseline range: [{M_baseline.min():.6f}, {M_baseline.max():.6f}]")
        print(f"M_SK range:       [{M_sk.min():.6f}, {M_sk.max():.6f}]")
        print()
        print("Residual Analysis:")
        print(f"  Max fractional residual:  {max_residual:.4%}")
        print(f"  Mean fractional residual: {mean_residual:.4%}")
        print(f"  Allowed threshold:        {max_residual_allowed:.4%}")
        print()
        
        if passes:
            print("✓ TEST PASSED: SK modulation matches baseline within 2%")
        else:
            print("✗ TEST FAILED: SK modulation deviates by more than 2%")
            print()
            print("Top 5 largest residuals:")
            idx_sorted = np.argsort(residuals)[::-1]
            for i in range(min(5, len(residuals))):
                idx = idx_sorted[i]
                print(f"  k={k_array[idx]:.4e}: M_base={M_baseline[idx]:.6f}, "
                      f"M_SK={M_sk[idx]:.6f}, residual={residuals[idx]:.4%}")
        
        print("=" * 70)
    
    return passes


def test_modulation_match_zero_epsilon(verbose: bool = True) -> bool:
    """
    Test that both implementations give M(k) = 1 when ε = 0.
    
    Returns
    -------
    passes : bool
        True if test passes
    """
    k_array = np.logspace(-4, 0, 50)
    epsilon = 0.0
    sigma = 0.1
    k_star = 0.05
    omega_n = np.array([14.13, 21.02])
    phi_n = np.array([0.0, np.pi/4])
    
    M_baseline = zeta_comb_modulation(
        k_array, epsilon, sigma, k_star, omega_n, phi_n
    )
    
    M_sk = compute_modulation(
        k_array, epsilon, sigma, omega_n, phi_n, k_star
    )
    
    # Both should be exactly 1.0
    baseline_correct = np.allclose(M_baseline, 1.0, atol=1e-10)
    sk_correct = np.allclose(M_sk, 1.0, atol=1e-10)
    
    passes = baseline_correct and sk_correct
    
    if verbose:
        print()
        print("Test: ε = 0 limit (both implementations)")
        print("-" * 70)
        print(f"Baseline: max |M-1| = {np.max(np.abs(M_baseline - 1)):.2e}")
        print(f"SK:       max |M-1| = {np.max(np.abs(M_sk - 1)):.2e}")
        
        if passes:
            print("✓ Both implementations correctly return M(k) = 1 for ε = 0")
        else:
            print("✗ Error: One or both implementations fail ε = 0 test")
    
    return passes


def test_single_mode(verbose: bool = True) -> bool:
    """
    Test with a single mode to verify the oscillation pattern.
    
    Returns
    -------
    passes : bool
        True if test passes
    """
    k_array = np.logspace(-4, 0, 200)
    epsilon = 0.05  # Larger for visibility
    sigma = 0.1
    k_star = 0.05
    omega_n = np.array([14.13])  # Single mode
    phi_n = np.array([0.0])
    
    M_baseline = zeta_comb_modulation(
        k_array, epsilon, sigma, k_star, omega_n, phi_n
    )
    
    M_sk = compute_modulation(
        k_array, epsilon, sigma, omega_n, phi_n, k_star
    )
    
    residuals = np.abs(M_sk - M_baseline) / np.abs(M_baseline)
    max_residual = np.max(residuals)
    
    passes = max_residual < 0.02
    
    if verbose:
        print()
        print("Test: Single-mode oscillation")
        print("-" * 70)
        print(f"Mode: ω = {omega_n[0]:.2f}, φ = {phi_n[0]:.2f}")
        print(f"Amplitude: ε = {epsilon}")
        print(f"Max residual: {max_residual:.4%}")
        
        if passes:
            print("✓ Single-mode test passed")
        else:
            print("✗ Single-mode test failed")
    
    return passes


# ============================================================================
# Main test runner
# ============================================================================

if __name__ == "__main__":
    print()
    print("=" * 70)
    print("Phase 2: Testing SK Modulation vs Baseline")
    print("=" * 70)
    print()
    
    # Run all tests
    test1_pass = test_modulation_match(verbose=True)
    test2_pass = test_modulation_match_zero_epsilon(verbose=True)
    test3_pass = test_single_mode(verbose=True)
    
    print()
    print("=" * 70)
    print("Test Summary")
    print("=" * 70)
    print(f"Main comparison test:     {'PASS ✓' if test1_pass else 'FAIL ✗'}")
    print(f"Epsilon = 0 test:         {'PASS ✓' if test2_pass else 'FAIL ✗'}")
    print(f"Single-mode test:         {'PASS ✓' if test3_pass else 'FAIL ✗'}")
    print()
    
    all_pass = test1_pass and test2_pass and test3_pass
    
    if all_pass:
        print("✓✓✓ ALL TESTS PASSED ✓✓✓")
        print()
        print("Phase 2 Criterion Met: SK modulation matches baseline within 2%")
        sys.exit(0)
    else:
        print("✗✗✗ SOME TESTS FAILED ✗✗✗")
        print()
        print("Phase 2 requires debugging before proceeding to Phase 3")
        sys.exit(1)
