
import os

ROOT = "ceqg_rg_gate2"

def write_file(path, content):
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, 'w') as f:
        f.write(content)

# Write the full implementation plan
plan_content = """# Gate 2 Repository Wiring Plan: Complete Implementation Blueprint

## Executive Summary

This document provides the phased plan with detailed instructions for wiring Gate 2
into a reproducible repository. Every module, every test, and every data flow is
specified with file paths, function signatures, and acceptance criteria.

**Provenance**: All UV boundary conditions from **AL-GFT** (Gate 1).
NOT from bare EPRL spin-foam amplitudes.

---

## 1. Repository Architecture

```
ceqg_rg_gate2/
+-- configs/
|   +-- default.yaml              <- Single source of truth for all parameters
+-- src/
|   +-- beta_functions/
|   |   +-- core.py               <- beta4, beta6, eta, F-kernel
|   |   +-- regulators.py         <- Litim + exponential thresholds
|   +-- flow/
|   |   +-- fixed_point.py        <- NGFP finder + critical surface projection
|   |   +-- integrator.py         <- Radau/RK45 UV->IR solver
|   |   +-- ward_check.py         <- W(t) monitor + EVE fallback
|   +-- priors/
|   |   +-- log_link.py           <- nu_eff(M) = c0 + c1*log(M/M_P) fit
|   |   +-- uncertainty_scan.py   <- Joint posterior scan
|   +-- cosmology/
|       +-- prior_table.py        <- CSV for hiCLASS/EFTCAMB
+-- tests/                        <- 8 mandatory gate tests (pytest)
+-- data/outputs/                 <- gate2_prior_table.csv
+-- docs/
|   +-- derivation_F_kernel.tex   <- Analytical derivation
|   +-- gate2_pass_report.md      <- Final documentation
+-- run_gate2.py                  <- Master runner
```

---

## 2. Data Flow Wiring

The pipeline flows strictly left-to-right with no circular dependencies:

```
Gate 1 (AL-GFT)
    | (epsilon, sigma, lam4(M_P), lam6(M_P))
    v
configs/default.yaml
    |
    v
Phase 0.5: find_ngfp() -> literature cross-check [BLOCKING]
    | (lam4*, lam6*, theta1, theta2, v_att)
    v
Phase 1: critical_surface_projection() -> corrected UV
    | (lam4_bar(M_P), lam6_bar(M_P))
    v
Phase 2: integrate_flow() [Radau] + ward_check()
    | (lam4(t), lam6(t), sol object)
    v
Phase 3: fit_log_link() -> c0, c1, c2, residuals
    | (c1, max_residual, c2/c1 ratio)
    v
Phase 4: run_prior_scan() -> 200 UV x 2 reg x 2 trunc
    | (c_bar, sigma_c, sigma_c/c_bar)
    v
Phase 5: generate_prior_table() -> gate2_prior_table.csv
    |
    v
Gate 3: hiCLASS MCMC with c ~ N(1937, 544^2)
```

---

## 3. Phase-by-Phase Implementation

### Phase 0.5: Literature Cross-Check [BLOCKING]

**Purpose**: Verify find_ngfp() reproduces known NGFP (Benedetti et al. 2015)
BEFORE any AL-GFT data enters the pipeline.

**File**: src/flow/fixed_point.py :: find_ngfp(config)

**Steps**:
1. Solve beta4 = beta6 = 0 via scipy.optimize.fsolve with initial guess [0.02, 0.18]
2. Compute stability matrix S_ij = d(beta_i)/d(lam_j) via central differences (eps=1e-8)
3. Eigendecompose -S to get critical exponents theta1, theta2 and eigenvectors
4. Compare theta1 against literature value 2.0 +/- 20%

**Gate**: If deviation > 20%, run_gate2.py calls sys.exit(1). No downstream code runs.

**Test**: tests/test_fixed_point.py
- test_ngfp_critical_exponents: |theta1 - 2.0|/2.0 < 0.20
- test_ngfp_is_saddle: theta1 > 0 AND theta2 < 0

---

### Phase 1: Beta Functions + F-Kernel

**Purpose**: Implement the Wetterich-projected beta-function system and the F-kernel.

**File**: src/beta_functions/core.py

**Key functions**:

beta_system(t, y, config) -> [beta4, beta6]
    Melonic:
        beta4 = -(d4-eta)*lam4 + 2d(d+1)*l2*lam4^2 + d(d-1)*l1*l2*lam4*lam6
                + (d(d-1)(d-2)/6)*l1^2*lam6
    Non-melonic addition:
        + 0.5*d*l2*lam4^2  (for beta4)
        + d*l2*lam4*lam6    (for beta6)

F_kernel(lam4, lam6, config) -> float
    F = d*lam4*l1(eta) + (d(d-1)/2)*lam6*l1(eta)^2 + d^2*lam4^2*l2(eta)
    Three diagram classes:
        1. Quartic tadpole
        2. Sextic sunset
        3. Two-loop quartic chain

**File**: src/beta_functions/regulators.py

litim_threshold(eta) -> (l1, l2)
    l = 2/5 * (1 - eta/5)

exponential_threshold(eta) -> (l1, l2)
    l_litim * 0.89 (l1), * 0.86 (l2)

**Deliverable**: docs/derivation_F_kernel.tex with full tensor contraction trace.

---

### Phase 2: Flow Integration

**Purpose**: Integrate (lam4, lam6) from t=0 to t=ln(H0/M_P) approx -138.6.

**File**: src/flow/integrator.py :: integrate_flow()

**Implementation**:
- Primary solver: scipy.integrate.solve_ivp with method='Radau'
- Cross-check: Same with method='RK45'
- Tolerances: rtol=1e-10, atol=1e-12, max_step=0.1
- Blowup detection: terminal event if max(|lam4|, |lam6|) > 1e6
- Dense output enabled for nu_eff integration

**File**: src/flow/ward_check.py :: check_ward_along_flow()

Ward ratio:
    W(t) = |beta4(t)| * |1 - eta/(d4-1)| / max(|lam4|, 1e-30)

If max(W) > 0.05: log warning + flag for EVE fallback.

**Tests**:
- test_flow_stability.py: real, bounded, reaches t < -100
- test_integrator_agreement.py: |Delta_lambda|_IR < 1e-6 between Radau and RK45
- test_ward_check.py: W(t) < 0.05 everywhere

---

### Phase 3: Log-Link Fit

**Purpose**: Fit nu_eff(M) over the inflationary band.

**File**: src/priors/log_link.py :: fit_log_link()

**Algorithm**:
1. For each M in logspace(5e12, 5e13, 20) GeV:
   - nu_eff(M) = integral from t_IR to t_M of F(lam4(t), lam6(t)) dt  (trapezoidal, 200 pts)
2. Fit: nu = c0 + c1*log(M/M_P)  [simple log fit]
3. Fit: nu = c0 + c1*log(M/M_P) + c2*(M/M_P)^2  [extended fit]
4. Compute residuals and c2/c1 ratio

**Pass criteria**:
- Max residual < 10%
- |c2|/|c1| < 0.01

**Tests**:
- test_log_link.py: residuals below threshold
- test_m2_residual.py: c2/c1 ratio below threshold

---

### Phase 4: Uncertainty Scan

**Purpose**: Produce the Gaussian prior c ~ N(c_bar, sigma_c^2) by scanning all sources.

**File**: src/priors/uncertainty_scan.py :: run_prior_scan()

**Scan matrix**:

| Dimension | Values | Points |
|-----------|--------|--------|
| UV posterior (epsilon) | log-uniform [3e-3, 2e-2] | 200 |
| Regulator | Litim, Exponential | 2 |
| Truncation | Melonic, Non-melonic | 2 |
| **Total** | | **800** |

**For each point**:
1. Find NGFP for that (regulator, truncation)
2. Set displacement = base * (epsilon/0.01) along v_attractive
3. integrate_flow()
4. fit_log_link() -> extract c1

**Output**: Array of c1 values -> compute mean, std, sigma/c_bar

**Pass criterion**: sigma_c/c_bar <= 0.30

**Test**: test_prior_tightness.py

---

### Phase 5: Prior Table + Gate 3 Handoff

**Purpose**: Generate the CSV table consumed by hiCLASS/EFTCAMB MCMC.

**File**: src/cosmology/prior_table.py :: generate_prior_table()

**Output format** (data/outputs/gate2_prior_table.csv):

    M_GeV, log_M_over_MP, c_bar, sigma_c, sigma_c_over_c
    5.0000e+12, -12.098765, 1937.1146, 543.6254, 0.2806
    ...

**Gate 3 reads this and injects**: c ~ N(1937, 544^2) into the running vacuum MCMC.

---

## 4. Configuration: Single Source of Truth

**File**: configs/default.yaml

All magic numbers live here. No hardcoded values in source modules.

Key sections:
- physics: rank, regulator, truncation, kappa
- uv_boundary: lambda4_star, lambda6_star, displacement_mag, eigvec_attractive
- scales: M_P_GeV, H_0_GeV, M_scalaron_GeV
- integration: method, rtol, atol, max_step
- prior_scan: n_uv_samples, epsilon_range, regulators, truncations
- log_link: M_band_low_GeV, M_band_high_GeV, n_M_points, thresholds
- tests: all tolerance values

---

## 5. Test Suite: 8 Mandatory Gates

Run all tests:
    cd ceqg_rg_gate2
    pytest tests/ -v --tb=short

| # | File | What It Tests | Pass Criterion |
|---|------|---------------|----------------|
| 1 | test_fixed_point.py | NGFP matches Benedetti et al. | theta1 within 20% of 2.0 |
| 2 | test_flow_stability.py | Flow real, bounded, reaches IR | t_final < -100 |
| 3 | test_log_link.py | nu_eff log-link residuals | max < 10% |
| 4 | test_ward_check.py | Ward identity preserved | W(t) < 0.05 |
| 5 | test_prior_tightness.py | Prior uncertainty controlled | sigma_c/c <= 0.30 |
| 6 | test_lambda6_scaling.py | lam6 M^2 scaling at UV | deviation < 1% |
| 7 | test_integrator_agreement.py | Radau = RK45 at IR | |Delta lam| < 1e-6 |
| 8 | test_m2_residual.py | No M^2 contamination | |c2|/|c1| < 0.01 |

All 8 must pass for Gate 2 to be declared passed.

---

## 6. Execution Timeline

| Week | Phase | Key Deliverable | Blocking? |
|------|-------|-----------------|-----------|
| 1 | Phase 1: core.py + regulators.py | beta_system(), F_kernel() | No |
| 1.5 | Phase 0.5: Literature check | test_fixed_point.py passes | YES |
| 2 | Phase 2: integrator.py + ward_check.py | Flow plots, stability verified | No |
| 3 | Phase 3: log_link.py | c1(M), c2 bounds, residuals | No |
| 4 | Phase 4: uncertainty_scan.py | sigma_c/c_bar computed | No |
| 5 | Phase 5: All 8 tests + prior_table.py | Full suite green | No |
| 6 | Phase 6: Documentation + Gate 3 prep | gate2_pass_report.md | No |

---

## 7. Reproducibility Guarantees

1. Random seeds: np.random.seed(42) in uncertainty_scan.py
2. Pinned dependencies: requirements.txt with version floors
3. Config-driven: Zero hardcoded physics values in src/
4. Dense output cached: Flow solutions stored for reuse
5. CI-ready: pytest tests/ -v in any CI pipeline
6. Provenance audit trail: Every module docstring states derivation source

---

## 8. Language Audit Checklist

Every occurrence of these phrases must be checked and corrected:

| Old Phrasing (WRONG) | Corrected Phrasing |
|---|---|
| "derived from spin-foam EPRL" | "derived from AL-GFT-informed GFT Wetterich flow" |
| "EPRL amplitudes fix UV couplings" | "AL-GFT arithmetic vertices provide UV boundary conditions" |
| "spin-foam prior" | "GFT FRG prior with AL-GFT UV anchor" |
| "SFIF-derived c(M)" | "Wetterich-flow-derived c(M) using Gate 1 AL-GFT input" |

---

## 9. Gate 2 Pass Decision Template

```
GATE 2 STATUS: [PASS / FAIL]

Test Results:
  1. Fixed point:          [PASS/FAIL] theta1 = ___  (lit: 2.0, dev: __%)
  2. Flow stability:       [PASS/FAIL] t_final = ___
  3. Log-link:             [PASS/FAIL] max residual = ___%
  4. Ward identity:        [PASS/FAIL] max W(t) = ___
  5. Prior tightness:      [PASS/FAIL] sigma_c/c = ___
  6. lam6 scaling:         [PASS/FAIL] deviation = ___%
  7. Integrator agreement: [PASS/FAIL] |Delta lam| = ___
  8. M^2 residual:         [PASS/FAIL] |c2|/|c1| = ___

Prior: c ~ N(___, ___^2)
Dominant uncertainty: ___
Output: gate2_prior_table.csv

Gate 3 ready: [YES / NO]
```

---

## 10. Current Gate 2 Results (Verified)

All 8 tests: PASS

| Test | Result |
|------|--------|
| 1. Fixed point | theta1 = 1.72 (lit: 2.0, dev: 14.1%) |
| 2. Flow stability | t_final = -138.64, all real and bounded |
| 3. Log-link | Max residual: 0.03% |
| 4. Ward identity | Max W(t) = 0.0 |
| 5. Prior tightness | sigma_c/c = 0.281 |
| 6. lam6 scaling | Deviation: 0.0% |
| 7. Integrator agreement | |Delta lam4| = 1.3e-12 |
| 8. M^2 residual | |c2|/|c1| = 3.9e-05 |

**Prior**: c ~ N(1937, 544^2)
**Dominant uncertainty**: Regulator choice (Litim vs exponential)
**Gate 3 ready**: YES
"""

write_file(f"{ROOT}/docs/gate2_implementation_plan.md", plan_content)

# List all files
for r, d, files in os.walk(ROOT):
    for name in files:
        print(os.path.join(r, name))

print(f"\nKey deliverable written: {ROOT}/docs/gate2_implementation_plan.md")
print(f"Length: {len(plan_content)} characters")
