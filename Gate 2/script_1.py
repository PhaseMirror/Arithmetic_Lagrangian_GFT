
import os
import textwrap

ROOT = "ceqg_rg_gate2"

def write_file(path, content):
    os.makedirs(os.path.dirname(path) if os.path.dirname(path) else '.', exist_ok=True)
    with open(path, 'w') as f:
        f.write(content.strip() + '\n')

dirs = [
    "", "src", "src/beta_functions", "src/flow", "src/priors",
    "src/cosmology", "src/utils", "tests", "data", "data/literature",
    "data/outputs", "notebooks", "docs", "configs",
]
for d in dirs:
    os.makedirs(os.path.join(ROOT, d), exist_ok=True)

# ---- Tree listing ----
tree = """
ceqg_rg_gate2/
|-- configs/
|   |-- default.yaml            # ALL physics params, tolerances, scan ranges
|-- src/
|   |-- __init__.py
|   |-- beta_functions/
|   |   |-- __init__.py          # Public API: beta_system, F_kernel
|   |   |-- core.py              # beta4, beta6, eta, F-kernel (3 diagram classes)
|   |   |-- regulators.py        # Litim + exponential threshold functions
|   |-- flow/
|   |   |-- __init__.py
|   |   |-- fixed_point.py       # NGFP finder, stability matrix, critical surface projection
|   |   |-- integrator.py        # Radau/RK45 UV->IR integration
|   |   |-- ward_check.py        # Ward ratio monitor + EVE fallback flag
|   |-- priors/
|   |   |-- __init__.py
|   |   |-- log_link.py          # nu_eff(M) computation + c0,c1,c2 fit
|   |   |-- uncertainty_scan.py  # Joint (UV x regulator x truncation x xi) scan
|   |-- cosmology/
|   |   |-- __init__.py
|   |   |-- prior_table.py       # CSV generator for hiCLASS/EFTCAMB
|-- tests/
|   |-- conftest.py              # Shared fixtures: config, ngfp, reference_flow
|   |-- test_fixed_point.py      # Test 1: literature NGFP within 20%
|   |-- test_flow_stability.py   # Test 2: real, bounded, 140 e-folds
|   |-- test_log_link.py         # Test 3: residuals < 10%
|   |-- test_ward_check.py       # Test 4: W(t) < 0.05
|   |-- test_prior_tightness.py  # Test 5: sigma_c/c <= 0.30
|   |-- test_lambda6_scaling.py  # Test 6: M^2 preserved at UV
|   |-- test_integrator_agreement.py  # Test 7: Radau vs RK45 < 1e-6
|   |-- test_m2_residual.py      # Test 8: |c2|/|c1| < 0.01
|-- data/
|   |-- literature/              # Benedetti et al. reference values
|   |-- outputs/                 # gate2_prior_table.csv
|-- notebooks/
|   |-- 01_flow_visualization.ipynb
|   |-- 02_prior_analysis.ipynb
|-- docs/
|   |-- derivation_F_kernel.tex  # Full derivation of F(lam4,lam6)
|   |-- gate2_pass_report.md     # Final pass/fail documentation
|-- run_gate2.py                 # Master runner (all phases)
|-- requirements.txt
|-- README.md
"""

write_file(f"{ROOT}/REPO_TREE.txt", tree)

# ---- Generate the wiring diagram as Mermaid ----
mermaid_code = """
flowchart TD
    A[configs/default.yaml] --> B[src/flow/fixed_point.py<br>find_ngfp]
    A --> C[src/beta_functions/core.py<br>beta_system + F_kernel]
    A --> D[src/beta_functions/regulators.py<br>litim / exponential]
    
    B --> E[Phase 0.5: Literature Check<br>test_fixed_point.py]
    B --> F[critical_surface_projection]
    
    F --> G[src/flow/integrator.py<br>integrate_flow]
    C --> G
    D --> G
    
    G --> H[src/flow/ward_check.py<br>W_t < 0.05]
    G --> I[src/priors/log_link.py<br>nu_eff + c0,c1,c2 fit]
    
    I --> J[src/priors/uncertainty_scan.py<br>200 UV x 2 reg x 2 trunc]
    
    J --> K[src/cosmology/prior_table.py<br>gate2_prior_table.csv]
    
    K --> L[Gate 3: hiCLASS MCMC<br>c ~ N_1937_544sq]
    
    E --> M{Pass?}
    M -->|No| N[BLOCKED]
    M -->|Yes| F
    
    H --> O{W < 0.05?}
    O -->|No| P[Fallback: EVE method]
    O -->|Yes| I
"""

write_file(f"{ROOT}/docs/wiring_diagram.mmd", mermaid_code)

# Now let's count what we generated
file_count = 0
for root_dir, subdirs, files in os.walk(ROOT):
    for f in files:
        file_count += 1
        
print(f"Repository scaffold generated: {ROOT}/")
print(f"Total files created: {file_count}")
print()
print("DIRECTORY TREE:")
print(tree)
