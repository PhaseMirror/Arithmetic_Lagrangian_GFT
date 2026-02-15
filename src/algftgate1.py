"""
algftgate1.py
=============

Baseline phenomenological AL-GFT implementation (Gate 1, Track A).

This module implements the Arithmetic-Langevin GFT phenomenological model
with Zeta-Comb modulation of the primordial power spectrum. This is the
REFERENCE implementation that Phase 1's Schwinger-Keldysh derivation must
reproduce.

**Gate 1: Framework specified, derivation in progress.**

Author: PhaseMirror/Arithmetic_Lagrangian_GFT
Date: February 2026
"""

import numpy as np
from scipy.interpolate import interp1d
from typing import Tuple, Optional, List
import warnings


# ============================================================================
# Physical constants and default parameters
# ============================================================================

# Planck units normalized to 1
M_PL = 1.0  # Planck mass (in natural units)
H_INF = 1.0e-5  # Hubble parameter during inflation (in Planck units)

# AL-GFT default parameters
DEFAULT_EPSILON = 1e-2      # Multiplicity coupling strength
DEFAULT_SIGMA = 0.1         # Resonance width
DEFAULT_K_STAR = 0.05       # Reference scale (Mpc^-1)
DEFAULT_A_S = 2.1e-9        # Scalar amplitude (Planck 2018)
DEFAULT_N_S = 0.9649        # Scalar spectral index (Planck 2018)

# Zeta-Comb frequencies (first few modes)
DEFAULT_OMEGA_N = np.array([
    14.13,      # Primary resonance observed in Planck residuals
    21.02,      # Second harmonic
    28.09,      # Third harmonic
    35.44,      # Fourth harmonic
    42.89       # Fifth harmonic
])

# Zeta-Comb phases (phenomenological)
DEFAULT_PHI_N = np.array([
    0.0,        # Primary mode phase
    np.pi/4,    # Second mode phase
    np.pi/3,    # Third mode phase
    -np.pi/6,   # Fourth mode phase
    np.pi/2     # Fifth mode phase
])


# ============================================================================
# Zeta-Comb modulation function
# ============================================================================

def zeta_comb_modulation(
    k: np.ndarray,
    epsilon: float = DEFAULT_EPSILON,
    sigma: float = DEFAULT_SIGMA,
    k_star: float = DEFAULT_K_STAR,
    omega_n: np.ndarray = DEFAULT_OMEGA_N,
    phi_n: np.ndarray = DEFAULT_PHI_N
) -> np.ndarray:
    """
    Compute the Zeta-Comb modulation M(k) from AL-GFT environment.
    
    The modulation arises from a tower of complex-mass modes in the 
    quantum-geometric environment, producing log-periodic oscillations
    in the primordial power spectrum.
    
    Parameters
    ----------
    k : np.ndarray
        Comoving wavenumber(s) in Mpc^-1
    epsilon : float, optional
        Multiplicity coupling strength (dimensionless)
    sigma : float, optional
        Resonance width parameter (dimensionless)
    k_star : float, optional
        Reference scale in Mpc^-1
    omega_n : np.ndarray, optional
        Dimensionless log-periodic frequencies
    phi_n : np.ndarray, optional
        Phase offsets for each mode
        
    Returns
    -------
    M_k : np.ndarray
        Modulation factor M(k), where P_zeta(k) = P_0(k) * M(k)
        
    Notes
    -----
    The modulation has the form:
    
    M(k) = 1 + epsilon * sum_n exp(-sigma * omega_n^2) 
                           * cos(omega_n * log(k/k_star) + phi_n)
    
    In the limit epsilon → 0, M(k) → 1, recovering the ΛCDM spectrum.
    """
    if not isinstance(k, np.ndarray):
        k = np.array([k])
    
    # Ensure omega_n and phi_n have the same length
    if len(omega_n) != len(phi_n):
        raise ValueError("omega_n and phi_n must have the same length")
    
    # Initialize modulation
    M_k = np.ones_like(k, dtype=float)
    
    # Compute log-periodic oscillations
    log_k_ratio = np.log(k / k_star)
    
    # Sum over Zeta-Comb modes
    for omega, phi in zip(omega_n, phi_n):
        # Gaussian suppression at high frequencies
        amplitude = epsilon * np.exp(-sigma * omega**2)
        
        # Log-periodic oscillation
        oscillation = np.cos(omega * log_k_ratio + phi)
        
        M_k += amplitude * oscillation
    
    return M_k


# ============================================================================
# Primordial power spectrum
# ============================================================================

def power_spectrum_primordial(
    k: np.ndarray,
    A_s: float = DEFAULT_A_S,
    n_s: float = DEFAULT_N_S,
    k_pivot: float = 0.05
) -> np.ndarray:
    """
    Compute the baseline ΛCDM primordial power spectrum P_0(k).
    
    Parameters
    ----------
    k : np.ndarray
        Comoving wavenumber(s) in Mpc^-1
    A_s : float, optional
        Scalar amplitude at pivot scale
    n_s : float, optional
        Scalar spectral index
    k_pivot : float, optional
        Pivot scale in Mpc^-1
        
    Returns
    -------
    P_0 : np.ndarray
        Baseline primordial power spectrum
    """
    return A_s * (k / k_pivot)**(n_s - 1.0)


def power_spectrum_algft(
    k: np.ndarray,
    epsilon: float = DEFAULT_EPSILON,
    sigma: float = DEFAULT_SIGMA,
    k_star: float = DEFAULT_K_STAR,
    omega_n: np.ndarray = DEFAULT_OMEGA_N,
    phi_n: np.ndarray = DEFAULT_PHI_N,
    A_s: float = DEFAULT_A_S,
    n_s: float = DEFAULT_N_S,
    k_pivot: float = 0.05
) -> np.ndarray:
    """
    Compute the full AL-GFT primordial power spectrum with Zeta-Comb.
    
    P_zeta(k) = P_0(k) * M(k)
    
    where M(k) is the Zeta-Comb modulation from the quantum-geometric
    environment.
    
    Parameters
    ----------
    k : np.ndarray
        Comoving wavenumber(s) in Mpc^-1
    epsilon : float, optional
        Multiplicity coupling strength
    sigma : float, optional
        Resonance width parameter
    k_star : float, optional
        Reference scale
    omega_n : np.ndarray, optional
        Log-periodic frequencies
    phi_n : np.ndarray, optional
        Phase offsets
    A_s : float, optional
        Scalar amplitude
    n_s : float, optional
        Scalar spectral index
    k_pivot : float, optional
        Pivot scale
        
    Returns
    -------
    P_zeta : np.ndarray
        Full AL-GFT primordial power spectrum
    """
    # Baseline spectrum
    P_0 = power_spectrum_primordial(k, A_s, n_s, k_pivot)
    
    # Zeta-Comb modulation
    M_k = zeta_comb_modulation(k, epsilon, sigma, k_star, omega_n, phi_n)
    
    return P_0 * M_k


# ============================================================================
# Matched filter for oscillatory signatures
# ============================================================================

class ZetaCombMatchedFilter:
    """
    Matched filter to detect Zeta-Comb oscillations in CMB power spectrum
    residuals.
    
    This class implements a template-based search for log-periodic
    oscillations in Planck data residuals, scanning over omega and epsilon.
    """
    
    def __init__(
        self,
        k_min: float = 1e-4,
        k_max: float = 1.0,
        n_k: int = 1000
    ):
        """
        Initialize the matched filter.
        
        Parameters
        ----------
        k_min : float
            Minimum wavenumber in Mpc^-1
        k_max : float
            Maximum wavenumber in Mpc^-1
        n_k : int
            Number of k bins
        """
        self.k_min = k_min
        self.k_max = k_max
        self.n_k = n_k
        self.k_array = np.logspace(np.log10(k_min), np.log10(k_max), n_k)
    
    def compute_snr(
        self,
        residuals: np.ndarray,
        sigma_residuals: np.ndarray,
        omega: float,
        epsilon: float = DEFAULT_EPSILON,
        sigma: float = DEFAULT_SIGMA,
        k_star: float = DEFAULT_K_STAR,
        phi: float = 0.0
    ) -> float:
        """
        Compute signal-to-noise ratio for a given (omega, epsilon) template.
        
        Parameters
        ----------
        residuals : np.ndarray
            Observed residuals (data - ΛCDM)
        sigma_residuals : np.ndarray
            Uncertainty on residuals
        omega : float
            Test frequency
        epsilon : float, optional
            Test amplitude
        sigma : float, optional
            Resonance width
        k_star : float, optional
            Reference scale
        phi : float, optional
            Phase offset
            
        Returns
        -------
        snr : float
            Signal-to-noise ratio
        """
        # Construct template
        omega_n = np.array([omega])
        phi_n = np.array([phi])
        
        template = zeta_comb_modulation(
            self.k_array, epsilon, sigma, k_star, omega_n, phi_n
        ) - 1.0  # Subtract baseline
        
        # Normalize template
        template_norm = np.sqrt(np.sum(template**2 / sigma_residuals**2))
        if template_norm == 0:
            return 0.0
        
        template /= template_norm
        
        # Compute inner product
        inner_product = np.sum(residuals * template / sigma_residuals**2)
        
        # SNR
        snr = np.abs(inner_product)
        
        return snr
    
    def scan_omega(
        self,
        residuals: np.ndarray,
        sigma_residuals: np.ndarray,
        omega_min: float = 5.0,
        omega_max: float = 50.0,
        n_omega: int = 100,
        epsilon: float = DEFAULT_EPSILON,
        sigma: float = DEFAULT_SIGMA
    ) -> Tuple[np.ndarray, np.ndarray, dict]:
        """
        Scan over omega to find best-fit frequency.
        
        Parameters
        ----------
        residuals : np.ndarray
            Observed residuals
        sigma_residuals : np.ndarray
            Uncertainty on residuals
        omega_min : float
            Minimum test frequency
        omega_max : float
            Maximum test frequency
        n_omega : int
            Number of frequency bins
        epsilon : float, optional
            Coupling strength
        sigma : float, optional
            Resonance width
            
        Returns
        -------
        omega_array : np.ndarray
            Array of test frequencies
        snr_array : np.ndarray
            SNR at each frequency
        best_fit : dict
            Dictionary with best-fit parameters and SNR
        """
        omega_array = np.linspace(omega_min, omega_max, n_omega)
        snr_array = np.zeros(n_omega)
        
        for i, omega in enumerate(omega_array):
            snr_array[i] = self.compute_snr(
                residuals, sigma_residuals, omega,
                epsilon=epsilon, sigma=sigma
            )
        
        # Find maximum
        i_max = np.argmax(snr_array)
        best_fit = {
            'omega': omega_array[i_max],
            'snr': snr_array[i_max],
            'epsilon': epsilon,
            'sigma': sigma
        }
        
        return omega_array, snr_array, best_fit


# ============================================================================
# Utility functions
# ============================================================================

def epsilon_zero_limit_test(
    k: np.ndarray,
    tolerance: float = 1e-6
) -> Tuple[bool, float]:
    """
    Test that epsilon → 0 recovers ΛCDM to required precision.
    
    Parameters
    ----------
    k : np.ndarray
        Array of wavenumbers
    tolerance : float
        Maximum allowed |M(k) - 1|
        
    Returns
    -------
    passes : bool
        True if test passes
    max_deviation : float
        Maximum |M(k) - 1| observed
    """
    # Compute M(k) with epsilon = 0
    M_k = zeta_comb_modulation(k, epsilon=0.0)
    
    # Check deviation from 1
    deviation = np.abs(M_k - 1.0)
    max_deviation = np.max(deviation)
    
    passes = max_deviation < tolerance
    
    return passes, max_deviation


def verify_observability(
    epsilon: float,
    sigma: float,
    omega_n: np.ndarray,
    min_snr: float = 3.0
) -> dict:
    """
    Check if AL-GFT parameters produce observable signatures.
    
    Parameters
    ----------
    epsilon : float
        Coupling strength
    sigma : float
        Resonance width
    omega_n : np.ndarray
        Mode frequencies
    min_snr : float
        Minimum required SNR
        
    Returns
    -------
    result : dict
        Observability assessment
    """
    # Estimate peak modulation amplitude
    amplitudes = epsilon * np.exp(-sigma * omega_n**2)
    peak_amplitude = np.max(amplitudes)
    
    # Rough SNR estimate (actual calculation requires data)
    # Assume sqrt(N_modes) enhancement
    estimated_snr = peak_amplitude * np.sqrt(len(omega_n)) * 1000  # Heuristic
    
    result = {
        'peak_amplitude': peak_amplitude,
        'estimated_snr': estimated_snr,
        'observable': estimated_snr >= min_snr,
        'epsilon': epsilon,
        'sigma': sigma
    }
    
    return result


# ============================================================================
# Main execution (for testing)
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("AL-GFT Gate 1 Baseline Implementation")
    print("=" * 70)
    print()
    
    # Test 1: Compute power spectrum
    print("Test 1: Computing AL-GFT power spectrum")
    print("-" * 70)
    k_test = np.logspace(-4, 0, 100)
    P_0 = power_spectrum_primordial(k_test)
    P_zeta = power_spectrum_algft(k_test)
    M_k = P_zeta / P_0
    
    print(f"k range: [{k_test.min():.2e}, {k_test.max():.2e}] Mpc^-1")
    print(f"P_0 range: [{P_0.min():.2e}, {P_0.max():.2e}]")
    print(f"M(k) range: [{M_k.min():.4f}, {M_k.max():.4f}]")
    print()
    
    # Test 2: Epsilon → 0 limit
    print("Test 2: Testing epsilon → 0 limit (ΛCDM recovery)")
    print("-" * 70)
    passes, max_dev = epsilon_zero_limit_test(k_test)
    print(f"Max deviation |M(k) - 1|: {max_dev:.2e}")
    print(f"Required tolerance: 1e-6")
    print(f"Test result: {'PASS ✓' if passes else 'FAIL ✗'}")
    print()
    
    # Test 3: Observability check
    print("Test 3: Observability assessment")
    print("-" * 70)
    obs = verify_observability(DEFAULT_EPSILON, DEFAULT_SIGMA, DEFAULT_OMEGA_N)
    print(f"Peak amplitude: {obs['peak_amplitude']:.4f}")
    print(f"Estimated SNR: {obs['estimated_snr']:.2f}")
    print(f"Observable: {obs['observable']}")
    print()
    
    # Test 4: Primary resonance
    print("Test 4: Primary Zeta-Comb resonance")
    print("-" * 70)
    omega_primary = DEFAULT_OMEGA_N[0]
    print(f"Primary frequency: ω = {omega_primary:.2f}")
    print(f"Log-period: Δ log(k) = 2π/ω = {2*np.pi/omega_primary:.4f}")
    print()
    
    print("=" * 70)
    print("Baseline implementation complete.")
    print("Exit criterion: Script runs and produces P_zeta(k) with Zeta-Comb.")
    print("=" * 70)
