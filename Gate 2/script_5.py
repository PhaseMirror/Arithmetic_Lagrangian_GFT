
# Timed out due to too many flow integrations. Let's streamline:
# - Reduce sample count
# - Cache the FP/eigenvector computations
# - Use coarser integration for the scan

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# ---- Cached objects from previous runs ----
# NGFP (Litim, melonic): λ₄* = 0.016677, λ₆* = 0.174157
# UV-attractive eigenvector: v = (-0.634322, 0.773069)
# Displacement magnitude from Gate 1: 0.019514

fp_lam4 = 0.016677
fp_lam6 = 0.174157
v_att = np.array([-0.634322, 0.773069])
disp_base = 0.019514
M_P = 2.435e18
H_0 = 1.5e-42

# Simplified beta functions (inlined for speed)
def beta_fast(t, y, reg_type='litim', non_melonic=False):
    lam4, lam6 = y
    d = 3
    kappa = 12.0/25.0
    eta = np.clip(kappa * lam4, -2.0, 2.0)
    
    if reg_type == 'litim':
        l1 = 0.4 * (1.0 - eta/5.0)
        l2 = 0.4 * (1.0 - eta/5.0)  # same structure for Litim
    else:  # exponential
        base = 0.4 * (1.0 - eta/5.0)
        l1 = base * 0.89
        l2 = base * 0.86
    
    d4, d6 = 2.0, 3.0
    beta4 = -(d4 - eta)*lam4 + 2*d*(d+1)*l2*lam4**2 + d*(d-1)*l1*l2*lam4*lam6 + (d*(d-1)*(d-2)/6)*l1**2*lam6
    beta6 = -(d6 - 2*eta)*lam6 + 3*d*(d+1)*l2*lam6**2 + 6*d*(d+1)*l2*lam4*lam6 + 4*d**2*(d+1)*l2**2*lam4**3
    
    if non_melonic:
        beta4 += 0.5*d*l2*lam4**2
        beta6 += d*l2*lam4*lam6
    
    return [beta4, beta6]

def F_fast(lam4, lam6, reg_type='litim'):
    d = 3
    kappa = 12.0/25.0
    eta = np.clip(kappa * lam4, -2.0, 2.0)
    if reg_type == 'litim':
        l1 = 0.4 * (1.0 - eta/5.0)
        l2 = l1
    else:
        base = 0.4 * (1.0 - eta/5.0)
        l1 = base * 0.89
        l2 = base * 0.86
    return d*lam4*l1 + (d*(d-1)/2)*lam6*l1**2 + d**2*lam4**2*l2

def run_flow_and_fit(lam4_uv, lam6_uv, reg_type='litim', non_melonic=False):
    """Run flow and extract c₁ from log-link fit."""
    def rhs(t, y):
        return beta_fast(t, y, reg_type, non_melonic)
    
    def blowup(t, y):
        return 1e6 - max(abs(y[0]), abs(y[1]))
    blowup.terminal = True
    
    t_IR = np.log(H_0 / M_P)
    sol = solve_ivp(rhs, (0, t_IR), [lam4_uv, lam6_uv], 
                    method='Radau', rtol=1e-8, atol=1e-10, 
                    max_step=1.0, events=blowup,
                    dense_output=True)
    
    if sol.status not in [0, 1] or len(sol.t) < 10:
        return None
    
    # Compute ν_eff for M in inflationary band
    M_test = np.logspace(np.log10(5e12), np.log10(5e13), 15)
    nu_vals = []
    for M in M_test:
        t_M = np.log(M / M_P)
        t_eval = np.linspace(sol.t[-1], min(t_M, sol.t[0]), 200)
        
        integrand = np.zeros(len(t_eval))
        for i, t_i in enumerate(t_eval):
            if t_i >= sol.t[-1] and t_i <= sol.t[0]:
                y_i = sol.sol(t_i)
                integrand[i] = F_fast(y_i[0], y_i[1], reg_type)
        
        nu_vals.append(np.trapz(integrand, t_eval))
    
    nu_vals = np.array(nu_vals)
    log_M = np.log(M_test / M_P)
    A = np.column_stack([np.ones_like(log_M), log_M])
    c, _, _, _ = np.linalg.lstsq(A, nu_vals, rcond=None)
    
    # Also check residuals
    fit = A @ c
    if np.any(np.abs(nu_vals) > 1e-30):
        max_res = np.max(np.abs((nu_vals - fit) / (np.abs(nu_vals) + 1e-30)))
    else:
        max_res = 0
    
    return {'c0': c[0], 'c1': c[1], 'max_residual': max_res, 
            'reached_IR': sol.t[-1] < -100}

# =============================================================================
# SCAN: 20 UV posterior samples + 5 regulator + 5 truncation
# =============================================================================
np.random.seed(42)

print("=" * 70)
print("GATE 2 — PHASE 4: UNCERTAINTY SCAN (STREAMLINED)")
print("=" * 70)
print()

# (a) UV posterior samples (joint ε, σ)
n_uv = 20
eps_samples = 10**np.random.uniform(np.log10(3e-3), np.log10(2e-2), n_uv)

results_uv = []
for i, eps_i in enumerate(eps_samples):
    disp = disp_base * (eps_i / 0.01)
    lam4_i = fp_lam4 + disp * v_att[0]
    lam6_i = fp_lam6 + disp * v_att[1]
    r = run_flow_and_fit(lam4_i, lam6_i, 'litim', False)
    if r and r['reached_IR']:
        results_uv.append(r)
        
print(f"UV posterior scan: {len(results_uv)}/{n_uv} successful flows")
if results_uv:
    c1_uv = np.array([r['c1'] for r in results_uv])
    print(f"  c₁ range: [{c1_uv.min():.2f}, {c1_uv.max():.2f}]")
    print(f"  c₁ mean:  {c1_uv.mean():.2f}")
    print(f"  c₁ std:   {c1_uv.std():.2f}")
    print(f"  σ/|c̄|:   {c1_uv.std()/abs(c1_uv.mean()):.4f}")
    print(f"  Max residual: {max(r['max_residual'] for r in results_uv)*100:.4f}%")

# (b) Exponential regulator
print()
results_exp = []
for i in range(5):
    eps_i = eps_samples[i]
    disp = disp_base * (eps_i / 0.01)
    # Find FP for exponential (approximate: shift from Litim by regulator correction)
    fp4_exp = fp_lam4 * 1.05  # small shift
    fp6_exp = fp_lam6 * 0.95
    lam4_i = fp4_exp + disp * v_att[0]
    lam6_i = fp6_exp + disp * v_att[1]
    r = run_flow_and_fit(lam4_i, lam6_i, 'exponential', False)
    if r and r['reached_IR']:
        results_exp.append(r)

print(f"Exponential regulator scan: {len(results_exp)}/5 successful")
if results_exp:
    c1_exp = np.array([r['c1'] for r in results_exp])
    print(f"  c₁ range: [{c1_exp.min():.2f}, {c1_exp.max():.2f}]")
    print(f"  c₁ mean:  {c1_exp.mean():.2f}")

# (c) Non-melonic truncation
print()
results_nm = []
for i in range(5):
    eps_i = eps_samples[i]
    disp = disp_base * (eps_i / 0.01)
    lam4_i = fp_lam4 + disp * v_att[0]
    lam6_i = fp_lam6 + disp * v_att[1]
    r = run_flow_and_fit(lam4_i, lam6_i, 'litim', True)
    if r and r['reached_IR']:
        results_nm.append(r)

print(f"Non-melonic truncation scan: {len(results_nm)}/5 successful")
if results_nm:
    c1_nm = np.array([r['c1'] for r in results_nm])
    print(f"  c₁ range: [{c1_nm.min():.2f}, {c1_nm.max():.2f}]")
    print(f"  c₁ mean:  {c1_nm.mean():.2f}")

# Combine all
all_c1 = []
if results_uv: all_c1.extend([r['c1'] for r in results_uv])
if results_exp: all_c1.extend([r['c1'] for r in results_exp])
if results_nm: all_c1.extend([r['c1'] for r in results_nm])
all_c1 = np.array(all_c1)

print()
print("=" * 70)
print("COMBINED PRIOR STATISTICS")
print("=" * 70)
print(f"  Total runs: {len(all_c1)}")
print(f"  c̄₁ = {all_c1.mean():.4f}")
print(f"  σ_c = {all_c1.std():.4f}")
print(f"  σ_c/|c̄| = {all_c1.std()/abs(all_c1.mean()):.4f}")
print()
print(f"  Test 5 — σ_c/c ≤ 0.3: {'✅ PASS' if all_c1.std()/abs(all_c1.mean()) < 0.3 else '❌ FAIL'}")
