# Arithmetic Lagrangian GFT

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

### Gate 2+ (scaffolding only, not yet started)

- Gate 2 will use AL-GFT-derived UV couplings \(\lambda_4(M_P) = (1.73 \pm 0.24) \times 10^{-3}\), \(\lambda_6(M_P) = (1.08 \pm 0.21) \times 10^{-4}\) as boundary conditions for a sextic GFT Wetterich flow (Litim regulator, melonic + first non-melonic truncation), producing an RG-derived prior for an IR running parameter. [file:1][file:6]
- Gates 3–5 are defined in the Blueprint but are **not** implemented here yet.

---

## 3. Repository structure

Proposed layout (adjust as you implement):

```text
.
├── README.md
├── docs/
│   ├── CEQG-RG-Langevin-Blueprint.pdf        # High-level Blueprint (5 gates) [external]
│   ├── AL-GFT-Gate1-TrackA-SK.tex            # Schwinger–Keldysh derivation note
│   ├── AL-GFT-Gate1-TrackA-SK.pdf
│   └── Gate1-PassCriteria.md                 # Explicit Gate 1 pass/fail checklist
├── src/
│   ├── algftgate1.py                         # Original phenomenological AL-GFT implementation
│   ├── algft_sk.py                           # SK-based Zeta-Comb noise kernel implementation
│   └── mapping_uv_gft.py                     # (planned) map (ε,σ,…) → (λ4(MP), λ6(MP))
├── tests/
│   ├── test_modulation_match.py              # SK vs original modulation M(k)
│   └── test_eps_zero_limit.py                # ε → 0 limit → ΛCDM spectrum
└── examples/
    ├── demo_power_spectrum.ipynb             # Plots of Pζ(k) with Zeta-Comb
    └── demo_gate1_to_gate2.ipynb             # Sketch: AL-GFT UV → GFT couplings
