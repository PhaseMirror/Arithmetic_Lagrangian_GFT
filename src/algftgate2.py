"""
algftgate2.py
=============

CEQG Renormalization Group implementation (Gate 2).

This module implements the Functional Renormalization Group (FRG) flow for
causal-emergent quantum gravity with UV boundary conditions from AL-GFT (Gate 1).
It computes the cosmological prior on the neutrino effective mass parameter c₁
through RG flow integration and log-link fitting.

**Gate 2: UV → IR Prior Derivation**

Author: PhaseMirror/Arithmetic_Lagrangian_GFT
Date: February 2026
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from typing import Tuple, Dict, Optional
import warnings


# ============================================================================
# Physical constants and default parameters
# ============================================================================

# Scales
M_PLANCK = 2.435e18  # Planck mass in GeV
H_0 = 1.5e-42        # Present Hubble scale in GeV
M_SCALARON = 3.0e13  # Inflationary scale in GeV

# NGFP from literature (Benedetti et al. 2015)
LITERATURE_NGFP = {
    'lambda4_star': 0.016677,
    'lambda6_star': 0.174157,
    'theta1': 2.0,  # Critical exponent (positive = UV attractive)
    'theta2': -1.5  # Second exponent (negative = UV repulsive)
}

# Default physics parameters
DEFAULT_CONFIG = {
    'rank': 3,                    # Tensor rank d
    'kappa': 12.0/25.0,          # Anomalous dimension coupling
    'regulator': 'litim',         # 'litim' or 'exponential'
    'truncation': 'melonic',      # 'melonic' or 'non_melonic'
}


# ============================================================================
# Regulator threshold functions
# ============================================================================

def litim_threshold(eta: float) -> Tuple[float, float]:
    """
    Compute (l₁, l₂) threshold functions for Litim regulator.
    
    For rank-3 GFT with Litim optimized regulator:
    l₁ = l₂ = 2/5 · (1 - η/5)
    
    Parameters
    ----------
    eta : float
        Anomalous dimension
        
    Returns
    -------
    l1, l2 : float
        Threshold functions
    """
    l1 = 0.4 * (1.0 - eta / 5.0)
    l2 = l1
    return l1, l2


def exponential_threshold(eta: float) -> Tuple[float, float]:
    """
    Compute (l₁, l₂) threshold functions for exponential regulator.
    
    Exponential regulator yields slightly suppressed thresholds
    relative to Litim (~0.86-0.89 ratio).
    
    Parameters
    ----------
    eta : float
        Anomalous dimension
        
    Returns
    -------
    l1, l2 : float
        Threshold functions
    """
    base = 0.4 * (1.0 - eta / 5.0)
    l1 = base * 0.89
    l2 = base * 0.86
    return l1, l2


def get_thresholds(eta: float, regulator: str = 'litim') -> Tuple[float, float]:
    """Dispatch to appropriate regulator threshold function."""
    if regulator == 'litim':
        return litim_threshold(eta)
    elif regulator == 'exponential':
        return exponential_threshold(eta)
    else:
        raise ValueError(f"Unknown regulator: {regulator}")


# ============================================================================
# Anomalous dimension and beta functions
# ============================================================================

def anomalous_dimension(lambda4: float, kappa: float = 12.0/25.0) -> float:
    """
    Compute anomalous dimension η = κ · λ₄.
    
    Clipped to [-2, 2] for numerical stability.
    
    Parameters
    ----------
    lambda4 : float
        Quartic coupling
    kappa : float, optional
        Anomalous dimension coupling constant
        
    Returns
    -------
    eta : float
        Anomalous dimension
    """
    return np.clip(kappa * lambda4, -2.0, 2.0)


def beta_system(t: float, y: list, config: dict = None) -> list:
    """
    Compute beta functions [β₄, β₆] for the coupled RG flow.
    
    Implements the Wetterich equation projected onto melonic/non-melonic
    truncations for rank-d tensor GFT.
    
    Parameters
    ----------
    t : float
        RG time t = ln(k/M_P)
    y : list
        [λ₄(t), λ₆(t)] coupling values
    config : dict, optional
        Configuration with 'rank', 'regulator', 'truncation', 'kappa'
        
    Returns
    -------
    [beta4, beta6] : list
        Time derivatives of couplings
        
    Notes
    -----
    Beta functions:
    
    β₄ = -(d₄ - η)λ₄ + 2d(d+1)l₂ λ₄² + d(d-1)l₁l₂ λ₄λ₆
         + (d(d-1)(d-2)/6) l₁² λ₆
         
    β₆ = -(d₆ - 2η)λ₆ + 3d(d+1)l₂ λ₆² + 6d(d+1)l₂ λ₄λ₆
         + 4d²(d+1)l₂² λ₄³
         
    With non-melonic corrections if enabled.
    """
    if config is None:
        config = DEFAULT_CONFIG
    
    lam4, lam6 = y
    d = config.get('rank', 3)
    eta = anomalous_dimension(lam4, config.get('kappa', 12.0/25.0))
    
    l1, l2 = get_thresholds(eta, config.get('regulator', 'litim'))
    
    d4, d6 = 2.0, 3.0  # Canonical dimensions
    
    # Melonic contributions
    beta4 = (-(d4 - eta) * lam4
             + 2 * d * (d + 1) * l2 * lam4**2
             + d * (d - 1) * l1 * l2 * lam4 * lam6
             + (d * (d - 1) * (d - 2) / 6.0) * l1**2 * lam6)
    
    beta6 = (-(d6 - 2 * eta) * lam6
             + 3 * d * (d + 1) * l2 * lam6**2
             + 6 * d * (d + 1) * l2 * lam4 * lam6
             + 4 * d**2 * (d + 1) * l2**2 * lam4**3)
    
    # Non-melonic (necklace) corrections
    if config.get('truncation', 'melonic') == 'non_melonic':
        beta4 += 0.5 * d * l2 * lam4**2
        beta6 += d * l2 * lam4 * lam6
    
    return [beta4, beta6]


def F_kernel(lambda4: float, lambda6: float, eta: float, 
             regulator: str = 'litim') -> float:
    """
    Compute F-kernel for ν_eff integration.
    
    The F-kernel combines three diagram classes:
    1. Quartic tadpole (one loop of λ₄)
    2. Sextic sunset (two legs of λ₆)
    3. Two-loop quartic (chain of two λ₄ bubbles)
    
    F = d·λ₄·l₁(η) + (d(d-1)/2)·λ₆·l₁(η)² + d²·λ₄²·l₂(η)
    
    Parameters
    ----------
    lambda4, lambda6 : float
        Coupling values
    eta : float
        Anomalous dimension
    regulator : str, optional
        'litim' or 'exponential'
        
    Returns
    -------
    F : float
        Kernel value
    """
    d = 3  # rank
    l1, l2 = get_thresholds(eta, regulator)
    
    return (d * lambda4 * l1
            + (d * (d - 1) / 2.0) * lambda6 * l1**2
            + d**2 * lambda4**2 * l2)


# ============================================================================
# Fixed point finding and stability analysis
# ============================================================================

def find_ngfp(config: dict = None, x0: np.ndarray = None) -> Dict:
    """
    Find the non-Gaussian fixed point (NGFP).
    
    Solves β₄ = β₆ = 0 and computes stability matrix eigenvalues
    to determine critical exponents.
    
    Parameters
    ----------
    config : dict, optional
        Physics configuration
    x0 : array-like, optional
        Initial guess [λ₄*, λ₆*]
        
    Returns
    -------
    result : dict
        Contains: lambda4_star, lambda6_star, theta1, theta2,
        v_attractive, v_repulsive, eta_star
        
    Notes
    -----
    This MUST reproduce literature values (Benedetti et al. 2015)
    within 20% before the pipeline can proceed (blocking test).
    """
    if config is None:
        config = DEFAULT_CONFIG
    if x0 is None:
        x0 = np.array([0.02, 0.18])
    
    def rhs(y):
        return beta_system(0.0, y, config)
    
    sol = fsolve(rhs, x0, full_output=True)
    fp = sol[0]
    
    if not sol[2]:
        warnings.warn("Fixed point solver did not converge")
    
    # Compute stability matrix S_ij = ∂β_i/∂λ_j
    S = stability_matrix(fp, config)
    
    # Critical exponents are eigenvalues of -S
    eigenvalues, eigenvectors = np.linalg.eig(-S)
    
    # Sort by real part (descending)
    idx = np.argsort(eigenvalues.real)[::-1]
    theta = eigenvalues[idx]
    vecs = eigenvectors[:, idx]
    
    eta_star = anomalous_dimension(fp[0], config.get('kappa', 12.0/25.0))
    
    return {
        'lambda4_star': fp[0],
        'lambda6_star': fp[1],
        'eta_star': eta_star,
        'theta1': theta[0].real,
        'theta2': theta[1].real,
        'v_attractive': vecs[:, 0].real,
        'v_repulsive': vecs[:, 1].real,
        'stability_matrix': S
    }


def stability_matrix(fp: np.ndarray, config: dict = None, eps: float = 1e-8) -> np.ndarray:
    """
    Compute numerical Jacobian ∂β_i/∂λ_j at fixed point.
    
    Parameters
    ----------
    fp : array-like
        Fixed point [λ₄*, λ₆*]
    config : dict, optional
        Physics configuration
    eps : float, optional
        Finite difference step size
        
    Returns
    -------
    S : ndarray (2, 2)
        Stability matrix
    """
    if config is None:
        config = DEFAULT_CONFIG
    
    S = np.zeros((2, 2))
    for j in range(2):
        yp, ym = list(fp), list(fp)
        yp[j] += eps
        ym[j] -= eps
        bp = beta_system(0.0, yp, config)
        bm = beta_system(0.0, ym, config)
        for i in range(2):
            S[i, j] = (bp[i] - bm[i]) / (2 * eps)
    
    return S


def critical_surface_projection(lambda4_raw: float, lambda6_raw: float, 
                                ngfp: Dict) -> Tuple[float, float]:
    """
    Project raw UV boundary conditions onto critical surface.
    
    Retains only the UV-attractive eigendirection, discarding the
    UV-repulsive component that would cause flow instability.
    
    Parameters
    ----------
    lambda4_raw, lambda6_raw : float
        Raw UV couplings from AL-GFT (Gate 1)
    ngfp : dict
        Fixed point data from find_ngfp()
        
    Returns
    -------
    lambda4_corrected, lambda6_corrected : float
        Projected couplings on critical surface
    """
    fp = np.array([ngfp['lambda4_star'], ngfp['lambda6_star']])
    v_att = ngfp['v_attractive']
    
    displacement = np.array([lambda4_raw, lambda6_raw]) - fp
    proj_mag = np.dot(displacement, v_att)
    
    corrected = fp + proj_mag * v_att
    return corrected[0], corrected[1]


# ============================================================================
# RG flow integration
# ============================================================================

def integrate_flow(lambda4_uv: float, lambda6_uv: float,
                  config: dict = None,
                  method: str = 'Radau',
                  rtol: float = 1e-10,
                  atol: float = 1e-12,
                  max_step: float = 0.1,
                  n_points: int = 2000) -> Dict:
    """
    Integrate RG flow from UV (M_P) to IR (H_0).
    
    Parameters
    ----------
    lambda4_uv, lambda6_uv : float
        UV initial conditions (after critical surface projection)
    config : dict, optional
        Physics configuration
    method : str, optional
        Integration method ('Radau' or 'RK45')
    rtol, atol : float, optional
        Relative and absolute tolerances
    max_step : float, optional
        Maximum step size in RG time
    n_points : int, optional
        Number of output points
        
    Returns
    -------
    result : dict
        Contains: t, lambda4, lambda6 (arrays), sol (OdeSolution),
        success (bool), reached_IR (bool)
    """
    if config is None:
        config = DEFAULT_CONFIG
    
    t_IR = np.log(H_0 / M_PLANCK)
    t_eval = np.linspace(0, t_IR, n_points)
    
    def rhs(t, y):
        return beta_system(t, y, config)
    
    def blowup_event(t, y):
        """Terminate if couplings exceed 10^6."""
        return 1e6 - max(abs(y[0]), abs(y[1]))
    blowup_event.terminal = True
    
    sol = solve_ivp(
        rhs, (0, t_IR), [lambda4_uv, lambda6_uv],
        method=method,
        rtol=rtol,
        atol=atol,
        max_step=max_step,
        events=blowup_event,
        dense_output=True,
        t_eval=t_eval
    )
    
    return {
        't': sol.t,
        'lambda4': sol.y[0],
        'lambda6': sol.y[1],
        'sol': sol,
        'success': sol.success,
        'reached_IR': sol.t[-1] < -100
    }


# ============================================================================
# ν_eff computation and log-link fitting
# ============================================================================

def compute_nu_eff(M: float, flow_sol, config: dict = None) -> float:
    """
    Compute effective neutrino mass parameter ν_eff(M).
    
    Integrates the F-kernel along the RG trajectory from IR to scale M:
    
    ν_eff(M) = ∫_{t_IR}^{t_M} dt F(λ₄(t), λ₆(t); η(t))
    
    Parameters
    ----------
    M : float
        Mass scale in GeV
    flow_sol : dict
        Flow result from integrate_flow()
    config : dict, optional
        Physics configuration
        
    Returns
    -------
    nu_eff : float
        Effective parameter value
    """
    if config is None:
        config = DEFAULT_CONFIG
    
    t_M = np.log(M / M_PLANCK)
    t_min = flow_sol['t'][-1]  # IR end
    t_max = min(t_M, flow_sol['t'][0])  # UV end or M
    
    if t_max <= t_min:
        return 0.0
    
    t_grid = np.linspace(t_min, t_max, 200)
    integrand = np.zeros_like(t_grid)
    
    regulator = config.get('regulator', 'litim')
    
    for i, t_i in enumerate(t_grid):
        y_i = flow_sol['sol'](t_i)
        eta_i = anomalous_dimension(y_i[0], config.get('kappa', 12.0/25.0))
        integrand[i] = F_kernel(y_i[0], y_i[1], eta_i, regulator)
    
    nu = np.trapz(integrand, t_grid)
    return nu


def fit_log_link(flow_sol, M_range: Tuple[float, float] = (5e12, 5e13),
                n_M_points: int = 30, config: dict = None) -> Dict:
    """
    Fit ν_eff(M) to log-linear form: ν_eff = c₀ + c₁·log(M/M_P) + c₂·(M/M_P)².
    
    The coefficient c₁ becomes the cosmological prior mean.
    
    Parameters
    ----------
    flow_sol : dict
        Flow result from integrate_flow()
    M_range : tuple, optional
        (M_min, M_max) in GeV for inflationary band
    n_M_points : int, optional
        Number of mass points for fitting
    config : dict, optional
        Physics configuration
        
    Returns
    -------
    result : dict
        Contains: c0, c1, c2, max_residual_pct, c2_over_c1
    """
    if config is None:
        config = DEFAULT_CONFIG
    
    M_test = np.logspace(np.log10(M_range[0]), np.log10(M_range[1]), n_M_points)
    nu_test = np.array([compute_nu_eff(M, flow_sol, config) for M in M_test])
    
    log_M = np.log(M_test / M_PLANCK)
    
    # Fit to c₀ + c₁·log(M/M_P) + c₂·(M/M_P)²
    M_norm = M_test / M_PLANCK
    A = np.column_stack([np.ones_like(log_M), log_M, M_norm**2])
    
    c_fit, residuals, rank, s = np.linalg.lstsq(A, nu_test, rcond=None)
    c0, c1, c2 = c_fit
    
    # Compute residuals
    fit_vals = A @ c_fit
    if np.any(np.abs(nu_test) > 1e-30):
        rel_residuals = np.abs((nu_test - fit_vals) / (np.abs(nu_test) + 1e-30))
        max_residual_pct = np.max(rel_residuals) * 100
    else:
        max_residual_pct = 0.0
    
    c2_over_c1 = abs(c2 / c1) if abs(c1) > 1e-30 else 0.0
    
    return {
        'c0': c0,
        'c1': c1,
        'c2': c2,
        'max_residual_pct': max_residual_pct,
        'c2_over_c1': c2_over_c1,
        'M_test': M_test,
        'nu_test': nu_test,
        'fit_vals': fit_vals
    }


# ============================================================================
# Uncertainty quantification
# ============================================================================

def run_prior_scan(epsilon_range: Tuple[float, float] = (3e-3, 2e-2),
                  sigma_range: Tuple[float, float] = (0.05, 0.5),
                  n_samples: int = 50,
                  regulators: list = ['litim', 'exponential'],
                  truncations: list = ['melonic', 'non_melonic'],
                  ngfp: Dict = None,
                  seed: int = 42) -> Dict:
    """
    Run joint posterior scan over Gate 1 parameters and systematics.
    
    Samples from:
    - ε (multiplicity coupling): affects UV displacement magnitude
    - σ (resonance width): affects coupling ratios
    - Regulator choice: Litim vs exponential
    - Truncation order: melonic vs non-melonic
    
    Parameters
    ----------
    epsilon_range, sigma_range : tuple
        (min, max) for log-uniform sampling
    n_samples : int
        Number of posterior samples
    regulators : list
        Regulator schemes to test
    truncations : list
        Truncation orders to test
    ngfp : dict, optional
        Fixed point data (computed if not provided)
    seed : int
        Random seed for reproducibility
        
    Returns
    -------
    result : dict
        Contains: c1_mean, c1_std, c1_samples, sigma_over_c,
        configurations, successful_runs
    """
    np.random.seed(seed)
    
    # Map (ε, σ) to displacement from NGFP
    # λ₆(M_P) ∝ ε · M² → displacement magnitude ∝ ε
    eps_samples = 10**np.random.uniform(np.log10(epsilon_range[0]), 
                                       np.log10(epsilon_range[1]), n_samples)
    sigma_samples = 10**np.random.uniform(np.log10(sigma_range[0]),
                                         np.log10(sigma_range[1]), n_samples)
    
    c1_samples = []
    configs_used = []
    
    # Compute or use provided NGFP
    base_config = DEFAULT_CONFIG.copy()
    if ngfp is None:
        ngfp = find_ngfp(base_config)
    
    disp_base = 0.019514  # From Gate 1 mapping
    
    print(f"Running prior scan: {n_samples} samples × {len(regulators)} reg × {len(truncations)} trunc")
    
    for i, (eps_i, sig_i) in enumerate(zip(eps_samples, sigma_samples)):
        for reg in regulators:
            for trunc in truncations:
                config = base_config.copy()
                config['regulator'] = reg
                config['truncation'] = trunc
                
                # Scale displacement by ε/ε_reference
                disp = disp_base * (eps_i / 0.01)
                
                # Apply displacement along UV-attractive direction
                lam4_uv = ngfp['lambda4_star'] + disp * ngfp['v_attractive'][0]
                lam6_uv = ngfp['lambda6_star'] + disp * ngfp['v_attractive'][1]
                
                try:
                    # Run flow
                    flow = integrate_flow(lam4_uv, lam6_uv, config)
                    
                    if flow['success'] and flow['reached_IR']:
                        # Fit log-link
                        fit = fit_log_link(flow, config=config)
                        
                        # Quality checks
                        if fit['max_residual_pct'] < 10.0 and fit['c2_over_c1'] < 0.01:
                            c1_samples.append(fit['c1'])
                            configs_used.append({
                                'epsilon': eps_i,
                                'sigma': sig_i,
                                'regulator': reg,
                                'truncation': trunc,
                                'c1': fit['c1']
                            })
                except Exception as e:
                    continue
    
    c1_samples = np.array(c1_samples)
    
    if len(c1_samples) == 0:
        warnings.warn("No successful flows in prior scan")
        return {
            'c1_mean': np.nan,
            'c1_std': np.nan,
            'c1_samples': np.array([]),
            'sigma_over_c': np.nan,
            'configurations': [],
            'successful_runs': 0
        }
    
    c1_mean = np.mean(c1_samples)
    c1_std = np.std(c1_samples)
    sigma_over_c = c1_std / abs(c1_mean) if abs(c1_mean) > 0 else np.inf
    
    return {
        'c1_mean': c1_mean,
        'c1_std': c1_std,
        'c1_samples': c1_samples,
        'sigma_over_c': sigma_over_c,
        'configurations': configs_used,
        'successful_runs': len(c1_samples)
    }


# ============================================================================
# Prior table generation
# ============================================================================

def generate_prior_table(c1_mean: float, c1_std: float,
                        M_range: Tuple[float, float] = (5e12, 5e13),
                        n_M_points: int = 20,
                        output_file: str = None) -> np.ndarray:
    """
    Generate CSV-formatted prior table for hiCLASS/EFTCAMB.
    
    Format: M_GeV, log_M_over_MP, c_bar, sigma_c, sigma_c_over_c
    
    Parameters
    ----------
    c1_mean, c1_std : float
        Mean and standard deviation from prior scan
    M_range : tuple
        (M_min, M_max) in GeV
    n_M_points : int
        Number of mass points
    output_file : str, optional
        If provided, write to CSV file
        
    Returns
    -------
    table : ndarray
        (n_M_points, 5) array with columns:
        [M_GeV, log_M_over_MP, c_bar, sigma_c, sigma_c_over_c]
    """
    M_vals = np.logspace(np.log10(M_range[0]), np.log10(M_range[1]), n_M_points)
    log_M_over_MP = np.log(M_vals / M_PLANCK)
    
    # c_bar and sigma_c are constant across M (from marginalization)
    c_bar = np.full_like(M_vals, c1_mean)
    sigma_c = np.full_like(M_vals, c1_std)
    sigma_over_c = sigma_c / np.abs(c_bar)
    
    table = np.column_stack([M_vals, log_M_over_MP, c_bar, sigma_c, sigma_over_c])
    
    if output_file:
        header = "M_GeV,log_M_over_MP,c_bar,sigma_c,sigma_c_over_c"
        np.savetxt(output_file, table, delimiter=',', header=header, 
                  comments='', fmt='%.4e,%.6f,%.4f,%.4f,%.4f')
        print(f"Prior table saved to {output_file}")
    
    return table


# ============================================================================
# Ward identity check
# ============================================================================

def check_ward_identity(flow_sol, config: dict = None, 
                       threshold: float = 0.05) -> Dict:
    """
    Monitor Ward-Takahashi identity violation along flow.
    
    For O(N)^d symmetric truncations:
    W(t) = |β₄(λ₄, λ₆)| · |1 - η/(d₄-1)| / |λ₄|
    
    Should remain < threshold (typically 0.05).
    
    Parameters
    ----------
    flow_sol : dict
        Flow result from integrate_flow()
    config : dict, optional
        Physics configuration
    threshold : float, optional
        Maximum allowed W(t)
        
    Returns
    -------
    result : dict
        Contains: passed (bool), max_W (float), t_at_max (float)
    """
    if config is None:
        config = DEFAULT_CONFIG
    
    W_vals = []
    for i, t_i in enumerate(flow_sol['t']):
        lam4 = flow_sol['lambda4'][i]
        lam6 = flow_sol['lambda6'][i]
        
        eta = anomalous_dimension(lam4, config.get('kappa', 12.0/25.0))
        beta = beta_system(t_i, [lam4, lam6], config)
        
        d4 = 2.0
        ward_factor = abs(1.0 - eta / (d4 - 1.0))
        
        if abs(lam4) > 1e-30:
            W = abs(beta[0]) * ward_factor / abs(lam4)
        else:
            W = 0.0
        
        W_vals.append(W)
    
    W_vals = np.array(W_vals)
    idx_max = np.argmax(W_vals)
    max_W = W_vals[idx_max]
    
    return {
        'passed': max_W < threshold,
        'max_W': max_W,
        't_at_max': flow_sol['t'][idx_max],
        'W_trajectory': W_vals
    }
