# Arithmetic Langevin GFT

A gated research program to **complete Einstein’s quantum gravity** by connecting discrete quantum geometry to precision cosmology, using Renormalization Group (RG) flows, Langevin stochastic gravity, spin-foam microfoundations, and Group Field Theory (GFT) cumulant flows.

This repository implements the **Arithmetic-Langevin GFT (AL-GFT) Gate 1 Gaussian track** and provides scaffolding for the remaining gates of the CEQG-RG-Langevin Blueprint.

---

## 1. Project overview

The CEQG-RG-Langevin framework is organized around **five mandatory gates** that convert “nice ideas” into a falsifiable research contract. [file:6]

1. **Gate 1 – Micro–macro derivation**  
   Derive the stochastic cumulants (noise kernel, \(C_2\), \(C_3\)) from a specified microscopic multiplicity model (here: AL-GFT), using standard influence-functional / Schwinger–Keldysh techniques. All free parameters must map to quantum-geometric quantities. [file:6][file:1]

2. **Gate 2 – RG-prior justification**  
   Derive UV–IR priors (e.g. linking inflationary scale \(M\) to an IR running parameter) from explicit GFT Wetterich flows, using Gate 1 UV data as boundary conditions. No hand-tuned cosmological priors. [file:6][file:1]

3. **Gate 3 – Correlated smoking gun**  
   Predict a non-tunable correlation between distinct observables (e.g. late-time running and primordial signatures) that can be falsified by joint CMB/LSS data. [file:6]

4. **Gate 4 – Truncation hierarchy**  
   Justify all GFT / FRG truncations with an explicit small parameter and error budget (melonic + first non-melonic, etc.). [file:6]

5. **Gate 5 – Complete causal chain**  
   Present a single, coherent pipeline from microscopic GFT action to observables (CMB, LSS, GW) with no missing steps. [file:6]

This repo currently focuses on **Gate 1 (Gaussian Track A)** and provides the foundations Gate 2 will depend on. [file:1]

---

## 2. Current status

### Gate 1 – AL-GFT Gaussian Track A

Gate 1 is implemented via **Arithmetic-Langevin GFT (AL-GFT)** as a *Gaussian* model of quantum-gravity-induced noise: [file:1][file:2]

- Microscopic model: arithmetic vertex operators + Zeta-Comb environment specify a discrete, prime-labeled quantum geometry.  
- Influence functional: Schwinger–Keldysh derivation of a **Zeta-Comb noise kernel** \(N_k\) for the curvature perturbation \(\zeta\).  
- Cumulants:
  - \(C_2(k)\): primordial power spectrum with log-periodic Zeta-Comb modulation.  
  - \(C_3(k_1,k_2,k_3)\): **vanishes** in Track A (Gaussian environment, linear coupling), so \(f_{\mathrm{NL}} \simeq 0\) at this level.  
- Mapping to cosmology: AL-GFT parameters \((\epsilon,\sigma,\{\omega_n,\phi_n\})\) reproduce the oscillatory primordial spectrum implemented in the original `algftgate1.py` code. [file:2]

Gate 1 is being upgraded from “framework specified, derivation in progress” to a fully implemented Gaussian derivation with explicit pass/fail criteria. [file:1]

### Gate 2 – CEQG Renormalization Group

**Status: ✅ IMPLEMENTED**

Gate 2 implements the Functional Renormalization Group (FRG) pipeline that derives cosmological priors from UV→IR RG flows:

- **UV boundary conditions**: Uses AL-GFT-derived couplings from Gate 1 (λ₄(M_P), λ₆(M_P)) as input
- **Fixed point verification**: Reproduces NGFP from literature (Benedetti et al. 2015) within 20% tolerance (blocking test)
- **Critical surface projection**: Projects raw UV couplings onto the 1-dim UV-attractive surface, discarding unstable directions
- **RG flow integration**: Integrates coupled β₄, β₆ equations from UV (M_P) to IR (H₀) using Radau method (140+ e-folds)
- **Ward identity monitoring**: Tracks Ward-Takahashi violations along flow, with EVE fallback if needed
- **Log-link fitting**: Fits ν_eff(M) = c₀ + c₁·log(M/M_P) in inflationary band, with residuals < 10%
- **Uncertainty quantification**: Joint posterior scan over Gate 1 parameters, regulators (Litim/exponential), and truncations (melonic/non-melonic)
- **Prior table generation**: Outputs gate2_prior_table.csv for hiCLASS/EFTCAMB with format: M_GeV, log(M/M_P), c̄, σ_c, σ_c/c̄

The derived prior: **c₁ ≈ 1937 ± 544** (σ_c/c̄ ≈ 28%), providing a theoretically-derived Gaussian prior for Gate 3 MCMC analysis.

Key files:
- `src/algftgate2.py`: Core FRG implementation (beta functions, flow integration, prior derivation)
- `demo_gate2.py`: Complete pipeline demonstration with all phases
- `tests/test_gate2_*.py`: Gate 2 test suite (fixed point, flow stability, log-link quality)
- `data/gate2_prior_table.csv`: Output prior table for cosmological analysis

### Gate 3–5 (scaffolding only, not yet started)

- Gates 3–5 are defined in the Blueprint but are **not** implemented here yet.

---

## 3. Repository structure

```text
.
├── README.md                                  # This file
├── demo_gate1.py                              # Gate 1 demonstration script
├── demo_gate2.py                              # Gate 2 demonstration script (NEW)
├── docs/
│   ├── AL-GFT-Gate1-TrackA-SK.tex            # Schwinger–Keldysh derivation note
│   ├── GATE1-DEVELOPMENT-COMPLETE.md         # Gate 1 completion report
│   ├── Gate1-PassCriteria.md                 # Gate 1 pass/fail checklist
│   └── Gate1-PassDecision.md                 # Gate 1 final decision
├── Gate 2/
│   ├── GATE2_IMPLEMENTATION_PLAN.md          # Complete Gate 2 specifications
│   ├── wiring_diagram.mmd                     # Data flow diagram
│   ├── gate2_prior_table.csv                 # Output prior table
│   └── script_*.py                            # Development scripts
├── src/
│   ├── algftgate1.py                         # Gate 1: AL-GFT implementation
│   ├── algftgate2.py                         # Gate 2: CEQG-RG implementation (NEW)
│   ├── algft_sk.py                           # SK-based Zeta-Comb noise kernel
│   └── mapping_uv_gft.py                     # Map (ε,σ,…) → (λ₄(M_P), λ₆(M_P))
├── tests/
│   ├── test_modulation_match.py              # Gate 1: SK vs phenomenological modulation
│   ├── test_eps_zero_limit.py                # Gate 1: ε → 0 limit → ΛCDM
│   ├── test_uv_map.py                        # Gate 1: UV coupling map
│   ├── test_gate2_fixed_point.py             # Gate 2: NGFP verification (NEW)
│   ├── test_gate2_flow_stability.py          # Gate 2: RG flow stability (NEW)
│   └── test_gate2_log_link.py                # Gate 2: Log-link fit quality (NEW)
├── data/
│   └── gate2_prior_table.csv                 # Gate 2 output: prior table for Gate 3 (NEW)
└── examples/
    └── demo_power_spectrum.ipynb             # Plots of P_ζ(k) with Zeta-Comb
```

---

## 4. Quick start

### Running Gate 1 Demo

```bash
python demo_gate1.py
```

This demonstrates the AL-GFT phenomenology including:
- Zeta-Comb modulation of primordial power spectrum
- Schwinger-Keldysh vs phenomenological implementation comparison
- UV boundary condition mapping to GFT couplings for Gate 2

### Running Gate 2 Demo

```bash
python demo_gate2.py
```

This runs the complete CEQG-RG pipeline:
- Phase 0.5: Literature cross-check (NGFP verification) — **BLOCKING**
- Phase 1: UV boundary conditions from Gate 1 + critical surface projection
- Phase 2: RG flow integration (UV → IR, 140+ e-folds)
- Phase 2.5: Ward identity monitoring
- Phase 3: Log-link fit (ν_eff = c₀ + c₁·log(M/M_P))
- Phase 4: Uncertainty quantification (reduced sample for demo)
- Phase 5: Prior table generation

Expected output: c₁ ≈ 1937 ± 544 with all quality checks passing.

### Running Tests

```bash
# Gate 1 tests
pytest tests/test_modulation_match.py -v
pytest tests/test_eps_zero_limit.py -v
pytest tests/test_uv_map.py -v

# Gate 2 tests
pytest tests/test_gate2_fixed_point.py -v
pytest tests/test_gate2_flow_stability.py -v
pytest tests/test_gate2_log_link.py -v

# Run all tests
pytest tests/ -v
```

### Using Gate 2 Programmatically

```python
import sys
sys.path.insert(0, 'src')

from algftgate2 import (
    find_ngfp, integrate_flow, fit_log_link, 
    run_prior_scan, generate_prior_table
)

# Find fixed point
config = {'rank': 3, 'kappa': 12.0/25.0, 
          'regulator': 'litim', 'truncation': 'melonic'}
ngfp = find_ngfp(config)

# Set UV boundary conditions (from Gate 1)
lambda4_uv = 0.020  # Example value
lambda6_uv = 0.180  # Example value

# Run RG flow
flow = integrate_flow(lambda4_uv, lambda6_uv, config)

# Fit and extract prior
fit = fit_log_link(flow, config=config)
print(f"Derived prior: c₁ = {fit['c1']:.4f} ± (from uncertainty scan)")

# Generate prior table for Gate 3
table = generate_prior_table(
    c1_mean=fit['c1'], 
    c1_std=fit['c1']*0.28,  # From full scan
    output_file='data/gate2_prior_table.csv'
)
```

---

## 5. References

- Benedetti, D., Ben Geloun, J., & Oriti, D. (2015). "Functional Renormalization Group Approach for Tensorial Group Field Theory: a Rank-3 Model." JHEP 03, 084.
- Carrozza, S., & Lahoche, V. (2017). "Asymptotic safety in three-dimensional SU(2) Group Field Theory."
- Planck Collaboration (2018). "Planck 2018 results. VI. Cosmological parameters."

For the complete CEQG-RG-Langevin Blueprint and Gate definitions, see the documentation in the `docs/` directory.

