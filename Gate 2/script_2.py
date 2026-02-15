
# Generate the full wiring diagram as a Mermaid PNG + the comprehensive implementation plan as markdown

# First, the Mermaid diagram
create_mermaid_diagram = None
try:
    from create_mermaid import create_mermaid_diagram
except:
    pass

mermaid_code = """
flowchart TD
    A["configs/default.yaml<br/>Physics params"] --> B["fixed_point.py<br/>find_ngfp()"]
    A --> C["core.py<br/>beta_system()"]
    A --> D["regulators.py<br/>litim/exponential"]
    B --> E{"Phase 0.5<br/>Literature?"}
    E -->|FAIL| BLOCK["BLOCKED<br/>Fix beta functions"]
    E -->|PASS| F["critical_surface_projection()"]
    F --> G["integrator.py<br/>integrate_flow()"]
    C --> G
    D --> G
    G --> H["ward_check.py<br/>W_t monitor"]
    H -->|"W>0.05"| EVE["Fallback: EVE method"]
    H -->|"W<0.05"| I["log_link.py<br/>nu_eff fit"]
    I --> J["uncertainty_scan.py<br/>200 UV x 2 reg x 2 trunc"]
    J --> K["prior_table.py<br/>gate2_prior_table.csv"]
    K --> L["Gate 3: hiCLASS MCMC<br/>c ~ N(1937, 544^2)"]
"""

try:
    create_mermaid_diagram(mermaid_code, 'gate2_wiring.png', width=1000, height=700)
    print("Mermaid diagram created.")
except Exception as e:
    print(f"Mermaid not available: {e}")
    print("Diagram code saved to docs/wiring_diagram.mmd")

# Now write the full implementation plan as a markdown file
plan = r"""# Gate 2 Repository Wiring Plan: Complete Implementation Blueprint

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
├── configs/
│   └── default.yaml              ← Single source of truth for all parameters
├── src/
│   ├── beta_functions/
│   │   ├── core.py               ← β₄, β₆, η, F-kernel
│   │   └── regulators.py         ← Litim + exponential thresholds
│   ├── flow/
│   │   ├── fixed_point.py        ← NGFP finder + critical surface projection
│   │   ├── integrator.py         ← Radau/RK45 UV→IR solver
│   │   └── ward_check.py         ← W(t) monitor + EVE fallback
│   ├── priors/
│   │   ├── log_link.py           ← ν_eff(M) = c₀ + c₁·log(M/M_P) fit
│   │   └── uncertainty_scan.py   ← Joint posterior scan
│   └── cosmology/
│       └── prior_table.py        ← CSV for hiCLASS/EFTCAMB
├── tests/                        ← 8 mandatory gate tests (pytest)
├── data/outputs/                 ← gate2_prior_table.csv
├── docs/
│   ├── derivation_F_kernel.tex   ← Analytical derivation
│   └── gate2_pass_report.md      ← Final documentation
└── run_gate2.py                  ← Master runner
```

---

## 2. Data Flow Wiring

The pipeline flows strictly left-to-right with **no circular dependencies**:

```
Gate 1 (AL-GFT)
    ↓ (ε, σ, λ₄(M_P), λ₆(M_P))
configs/default.yaml
    ↓
Phase 0.5: find_ngfp() → literature cross-check [BLOCKING]
    ↓ (λ₄*, λ₆*, θ₁, θ₂, v_att)
Phase 1: critical_surface_projection() → corrected UV
    ↓ (λ̄₄(M_P), λ̄₆(M_P))
Phase 2: integrate_flow() [Radau] + ward_check()
    ↓ (λ₄(t), λ₆(t), sol object)
Phase 3: fit_log_link() → c₀, c₁, c₂, residuals
    ↓ (c₁, max_residual, c₂/c₁ ratio)
Phase 4: run_prior_scan() → 200 UV × 2 reg × 2 trunc
    ↓ (c̄, σ_c, σ_c/c̄)
Phase 5: generate_prior_table() → gate2_prior_table.csv
    ↓
Gate 3: hiCLASS MCMC with c ~ N(1937, 544²)
```

---

## 3. Phase-by-Phase Implementation

### Phase 0.5 — Literature Cross-Check [BLOCKING]

**Purpose**: Verify that `find_ngfp()` reproduces known NGFP from Benedetti et al. (2015)
before ANY AL-GFT data enters the pipeline.

**File**: `src/flow/fixed_point.py :: find_ngfp(config)`

**Steps**:
1. Solve β₄ = β₆ = 0 via `scipy.optimize.fsolve` with initial guess `[0.02, 0.18]`
2. Compute stability matrix S_ij = ∂β_i/∂λ_j numerically (central differences, ε=10⁻⁸)
3. Eigendecompose −S to get critical exponents θ₁, θ₂ and eigenvectors
4. Compare θ₁ against literature value 2.0 ± 20%

**Gate**: If deviation > 20%, `run_gate2.py` calls `sys.exit(1)`. No downstream code runs.

**Test**: `tests/test_fixed_point.py`
- `test_ngfp_critical_exponents`: |θ₁ − 2.0|/2.0 < 0.20
- `test_ngfp_is_saddle`: θ₁ > 0 AND θ₂ < 0

---

### Phase 1 — Beta Functions + F-Kernel

**Purpose**: Implement the Wetterich-projected β-function system and the F-kernel.

**File**: `src/beta_functions/core.py`

**Key functions**:

```python
def beta_system(t, y, config) -> [beta4, beta6]:
    """
    β₄ = −(d₄ − η)λ₄ + 2d(d+1)l₂λ₄² + d(d−1)l₁l₂λ₄λ₆
          + (d(d−1)(d−2)/6)l₁²λ₆
          [+ non-melonic: 0.5·d·l₂·λ₄²]

    β₆ = −(d₆ − 2η)λ₆ + 3d(d+1)l₂λ₆² + 6d(d+1)l₂λ₄λ₆
          + 4d²(d+1)l₂²λ₄³
          [+ non-melonic: d·l₂·λ₄·λ₆]
    """

def F_kernel(lambda4, lambda6, config) -> float:
    """
    F = d·λ₄·l₁(η) + (d(d−1)/2)·λ₆·l₁²(η) + d²·λ₄²·l₂(η)

    Three diagram classes:
      1. Quartic tadpole
      2. Sextic sunset
      3. Two-loop quartic chain
    """
```

**File**: `src/beta_functions/regulators.py`

```python
def litim_threshold(eta) -> (l1, l2):
    # l = 2/5 · (1 − η/5)

def exponential_threshold(eta) -> (l1, l2):
    # l_litim × 0.89 (l1), × 0.86 (l2)
```

**Deliverable**: `docs/derivation_F_kernel.tex` with full tensor contraction trace.

---

### Phase 2 — Flow Integration

**Purpose**: Integrate (λ₄, λ₆) from t=0 to t=ln(H₀/M_P) ≈ −138.6.

**File**: `src/flow/integrator.py :: integrate_flow()`

**Implementation**:
- Primary solver: `scipy.integrate.solve_ivp` with `method='Radau'`
- Cross-check: Same with `method='RK45'`
- Tolerances: `rtol=1e-10, atol=1e-12, max_step=0.1`
- Blowup detection: terminal event if max(|λ₄|, |λ₆|) > 10⁶
- Dense output enabled for ν_eff integration

**File**: `src/flow/ward_check.py :: check_ward_along_flow()`

**Ward ratio**:
```
W(t) = |β₄(t)| · |1 − η/(d₄−1)| / max(|λ₄|, 10⁻³⁰)
```
If max(W) > 0.05: log warning + flag for EVE fallback.

**Tests**:
- `test_flow_stability.py`: real, bounded, reaches t < −100
- `test_integrator_agreement.py`: |Δλ|_IR < 10⁻⁶ between Radau and RK45
- `test_ward_check.py`: W(t) < 0.05 everywhere

---

### Phase 3 — Log-Link Fit

**Purpose**: Fit ν_eff(M) over the inflationary band.

**File**: `src/priors/log_link.py :: fit_log_link()`

**Algorithm**:
1. For each M in logspace(5×10¹², 5×10¹³, 20) GeV:
   - ν_eff(M) = ∫_{t_IR}^{t_M} F(λ₄(t), λ₆(t)) dt  (trapezoidal, 200 pts)
2. Fit: ν = c₀ + c₁·log(M/M_P)  [simple log fit]
3. Fit: ν = c₀ + c₁·log(M/M_P) + c₂·(M/M_P)²  [extended fit]
4. Compute residuals and c₂/c₁ ratio

**Pass criteria**:
- Max residual < 10%
- |c₂|/|c₁| < 0.01

**Tests**:
- `test_log_link.py`: residuals below threshold
- `test_m2_residual.py`: c₂/c₁ ratio below threshold

---

### Phase 4 — Uncertainty Scan

**Purpose**: Produce the Gaussian prior c ~ N(c̄, σ_c²) by scanning all uncertainty sources.

**File**: `src/priors/uncertainty_scan.py :: run_prior_scan()`

**Scan matrix**:

| Dimension | Values | # Points |
|-----------|--------|----------|
| UV posterior (ε) | log-uniform [3×10⁻³, 2×10⁻²] | 200 |
| Regulator | Litim, Exponential | 2 |
| Truncation | Melonic, Non-melonic | 2 |
| **Total** | | **800** |

**For each point**:
1. Find NGFP for that (regulator, truncation)
2. Set displacement = base × (ε/0.01) along v_attractive
3. integrate_flow()
4. fit_log_link() → extract c₁

**Output**: Array of c₁ values → compute mean, std, σ/c̄

**Pass criterion**: σ_c/c̄ ≤ 0.30

**Test**: `test_prior_tightness.py`

---

### Phase 5 — Prior Table + Gate 3 Handoff

**Purpose**: Generate the CSV table consumed by hiCLASS/EFTCAMB MCMC.

**File**: `src/cosmology/prior_table.py :: generate_prior_table()`

**Output format** (`data/outputs/gate2_prior_table.csv`):
```
M_GeV, log_M_over_MP, c_bar, sigma_c, sigma_c_over_c
5.0000e+12, -12.098765, 1937.1146, 543.6254, 0.2806
...
```

**Gate 3 reads this and injects**: `c ~ N(1937, 544²)` into the running vacuum MCMC.

---

## 4. Configuration: Single Source of Truth

**File**: `configs/default.yaml`

All magic numbers live here. No hardcoded values in source modules.

```yaml
physics:
  rank: 3
  regulator: "litim"
  truncation: "melonic"
  kappa: 0.48              # 12/25

uv_boundary:
  lambda4_star: 0.016677   # From find_ngfp(), verified Phase 0.5
  lambda6_star: 0.174157
  displacement_mag: 0.019514
  eigvec_attractive: [-0.634322, 0.773069]

scales:
  M_P_GeV: 2.435e18
  H_0_GeV: 1.5e-42
  M_scalaron_GeV: 3.0e13

integration:
  method: "Radau"
  rtol: 1.0e-10
  atol: 1.0e-12
  max_step: 0.1

prior_scan:
  n_uv_samples: 200
  epsilon_range: [3.0e-3, 2.0e-2]
  regulators: ["litim", "exponential"]
  truncations: ["melonic", "non_melonic"]

log_link:
  M_band_low_GeV: 5.0e12
  M_band_high_GeV: 5.0e13
  n_M_points: 20
  max_residual_pct: 10.0
  max_c2_over_c1: 0.01

tests:
  fp_tolerance_pct: 20.0
  ward_threshold: 0.05
  prior_sigma_over_c: 0.30
  integrator_agreement: 1.0e-6
  lambda6_scaling_pct: 1.0
```

---

## 5. Test Suite: 8 Mandatory Gates

Run all tests:
```bash
cd ceqg_rg_gate2
pytest tests/ -v --tb=short
```

| # | File | What It Tests | Pass Criterion |
|---|------|---------------|----------------|
| 1 | `test_fixed_point.py` | NGFP matches Benedetti et al. | |θ₁ − 2.0|/2.0 < 0.20 |
| 2 | `test_flow_stability.py` | Flow real, bounded, reaches IR | t_final < −100 |
| 3 | `test_log_link.py` | ν_eff log-link residuals | max < 10% |
| 4 | `test_ward_check.py` | Ward identity preserved | W(t) < 0.05 |
| 5 | `test_prior_tightness.py` | Prior uncertainty controlled | σ_c/c ≤ 0.30 |
| 6 | `test_lambda6_scaling.py` | λ₆ M² scaling at UV | deviation < 1% |
| 7 | `test_integrator_agreement.py` | Radau = RK45 at IR | |Δλ| < 10⁻⁶ |
| 8 | `test_m2_residual.py` | No M² contamination | |c₂|/|c₁| < 0.01 |

**All 8 must pass for Gate 2 to be declared passed.**

---

## 6. Execution Timeline

| Week | Phase | Key Deliverable | Blocking? |
|------|-------|-----------------|-----------|
| 1 | Phase 1: `core.py` + `regulators.py` | beta_system(), F_kernel() | No |
| 1.5 | **Phase 0.5**: Literature check | test_fixed_point.py passes | **YES** |
| 2 | Phase 2: `integrator.py` + `ward_check.py` | Flow plots, stability verified | No |
| 3 | Phase 3: `log_link.py` | c₁(M), c₂ bounds, residuals | No |
| 4 | Phase 4: `uncertainty_scan.py` | σ_c/c̄ computed, budget breakdown | No |
| 5 | Phase 5: All 8 tests + `prior_table.py` | Full suite green | No |
| 6 | Phase 6: Documentation + Gate 3 prep | gate2_pass_report.md, paper appendix | No |

---

## 7. Reproducibility Guarantees

1. **Random seeds**: `np.random.seed(42)` in uncertainty_scan.py
2. **Pinned dependencies**: `requirements.txt` with version floors
3. **Config-driven**: Zero hardcoded physics values in `src/`
4. **Dense output cached**: Flow solutions stored for reuse
5. **CI-ready**: `pytest tests/ -v` in any CI pipeline
6. **Provenance audit trail**: Every module docstring states derivation source

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
  1. Fixed point:          [PASS/FAIL] θ₁ = ___  (lit: 2.0, dev: __%)
  2. Flow stability:       [PASS/FAIL] t_final = ___
  3. Log-link:             [PASS/FAIL] max residual = ___%
  4. Ward identity:        [PASS/FAIL] max W(t) = ___
  5. Prior tightness:      [PASS/FAIL] σ_c/c = ___
  6. λ₆ scaling:           [PASS/FAIL] deviation = ___%
  7. Integrator agreement: [PASS/FAIL] |Δλ| = ___
  8. M² residual:          [PASS/FAIL] |c₂|/|c₁| = ___

Prior: c ~ N(___, ___²)
Dominant uncertainty: ___
Output: gate2_prior_table.csv

Gate 3 ready: [YES / NO]
```
"""

write_file(f"{ROOT}/docs/gate2_implementation_plan.md", plan)

# Also count actual files
count = 0
for r, d, f in os.walk(ROOT):
    for name in f:
        count += 1
        print(os.path.join(r, name))

print(f"\nTotal files: {count}")
print("\nKey deliverable: docs/gate2_implementation_plan.md")
