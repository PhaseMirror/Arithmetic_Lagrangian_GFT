"""
test_gate2_fixed_point.py
==========================

Test Phase 0.5 (BLOCKING): Literature cross-check for NGFP.

The fixed point finder MUST reproduce known results from
Benedetti et al. 2015 within 20% tolerance before any
AL-GFT data enters the pipeline.
"""

import pytest
import numpy as np
import sys
sys.path.insert(0, 'src')

from algftgate2 import find_ngfp, LITERATURE_NGFP, DEFAULT_CONFIG


def test_ngfp_critical_exponents():
    """Test that θ₁ matches literature within 20%."""
    config = DEFAULT_CONFIG.copy()
    ngfp = find_ngfp(config)
    
    theta1_lit = LITERATURE_NGFP['theta1']
    theta1_computed = ngfp['theta1']
    
    deviation = abs(theta1_computed - theta1_lit) / theta1_lit
    
    print(f"\nCritical exponent θ₁:")
    print(f"  Literature: {theta1_lit:.4f}")
    print(f"  Computed:   {theta1_computed:.4f}")
    print(f"  Deviation:  {deviation*100:.2f}%")
    
    assert deviation < 0.20, f"θ₁ deviation {deviation*100:.2f}% exceeds 20% tolerance"


def test_ngfp_is_saddle():
    """Test that NGFP has one positive and one negative critical exponent."""
    config = DEFAULT_CONFIG.copy()
    ngfp = find_ngfp(config)
    
    theta1 = ngfp['theta1']
    theta2 = ngfp['theta2']
    
    print(f"\nNGFP stability:")
    print(f"  θ₁ = {theta1:.4f} (should be > 0, UV attractive)")
    print(f"  θ₂ = {theta2:.4f} (should be < 0, UV repulsive)")
    
    assert theta1 > 0, f"θ₁ = {theta1:.4f} should be positive (UV attractive)"
    assert theta2 < 0, f"θ₂ = {theta2:.4f} should be negative (UV repulsive)"


def test_ngfp_beta_vanishes():
    """Test that beta functions vanish at the fixed point."""
    from algftgate2 import beta_system
    
    config = DEFAULT_CONFIG.copy()
    ngfp = find_ngfp(config)
    
    y_fp = [ngfp['lambda4_star'], ngfp['lambda6_star']]
    beta = beta_system(0.0, y_fp, config)
    
    print(f"\nBeta functions at fixed point:")
    print(f"  β₄ = {beta[0]:.6e}")
    print(f"  β₆ = {beta[1]:.6e}")
    
    assert abs(beta[0]) < 1e-6, f"β₄ = {beta[0]:.6e} should vanish at FP"
    assert abs(beta[1]) < 1e-6, f"β₆ = {beta[1]:.6e} should vanish at FP"


def test_ngfp_eigenvectors_orthogonal():
    """Test that critical eigenvectors are orthogonal."""
    config = DEFAULT_CONFIG.copy()
    ngfp = find_ngfp(config)
    
    v_att = ngfp['v_attractive']
    v_rep = ngfp['v_repulsive']
    
    dot_product = np.dot(v_att, v_rep)
    
    print(f"\nEigenvector orthogonality:")
    print(f"  v_att · v_rep = {dot_product:.6e}")
    
    assert abs(dot_product) < 1e-6, f"Eigenvectors not orthogonal: dot = {dot_product:.6e}"


def test_ngfp_stability_across_regulators():
    """Test that NGFP exists and is qualitatively similar for different regulators."""
    regulators = ['litim', 'exponential']
    results = {}
    
    for reg in regulators:
        config = DEFAULT_CONFIG.copy()
        config['regulator'] = reg
        results[reg] = find_ngfp(config)
    
    print(f"\nNGFP across regulators:")
    for reg in regulators:
        r = results[reg]
        print(f"  {reg:12s}: λ₄* = {r['lambda4_star']:.6f}, "
              f"λ₆* = {r['lambda6_star']:.6f}, θ₁ = {r['theta1']:.4f}")
    
    # Both should have positive first critical exponent
    for reg in regulators:
        assert results[reg]['theta1'] > 0, f"{reg} gives θ₁ < 0"


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "-s"])
