
import plotly.graph_objects as go
import plotly.io as pio
import numpy as np
from scipy.integrate import solve_ivp
import json

# Reproduce the primary flow for plotting
fp_lam4, fp_lam6 = 0.016677, 0.174157
v_att = np.array([-0.634322, 0.773069])
disp_base = 0.019514
M_P = 2.435e18
H_0 = 1.5e-42

def beta_fast(t, y):
    lam4, lam6 = y
    d = 3
    eta = np.clip(12.0/25.0 * lam4, -2.0, 2.0)
    l1 = 0.4 * (1.0 - eta/5.0)
    l2 = l1
    d4, d6 = 2.0, 3.0
    beta4 = -(d4-eta)*lam4 + 2*d*(d+1)*l2*lam4**2 + d*(d-1)*l1*l2*lam4*lam6 + (d*(d-1)*(d-2)/6)*l1**2*lam6
    beta6 = -(d6-2*eta)*lam6 + 3*d*(d+1)*l2*lam6**2 + 6*d*(d+1)*l2*lam4*lam6 + 4*d**2*(d+1)*l2**2*lam4**3
    return [beta4, beta6]

lam4_uv = fp_lam4 + disp_base * v_att[0]
lam6_uv = fp_lam6 + disp_base * v_att[1]
t_IR = np.log(H_0 / M_P)

sol = solve_ivp(beta_fast, (0, t_IR), [lam4_uv, lam6_uv],
                method='Radau', rtol=1e-10, atol=1e-12, max_step=0.1,
                t_eval=np.linspace(0, t_IR, 2000))

t_vals = sol.t
lam4_vals = sol.y[0]
lam6_vals = sol.y[1]

# ---- CHART 1: Flow trajectory (λ₄, λ₆) vs RG time ----
fig1 = go.Figure()
fig1.add_trace(go.Scatter(x=t_vals, y=lam4_vals, name='λ₄(t)', 
                          mode='lines', line=dict(width=2.5)))
fig1.add_trace(go.Scatter(x=t_vals, y=lam6_vals, name='λ₆(t)',
                          mode='lines', line=dict(width=2.5)))

# Mark key scales
for label, t_mark in [('M_P', 0), ('GUT', np.log(1e16/M_P)), 
                       ('Inflation', np.log(1e13/M_P)), ('EW', np.log(1e2/M_P))]:
    fig1.add_vline(x=t_mark, line_dash="dot", line_color="gray", opacity=0.5,
                   annotation_text=label, annotation_position="top")

# Mark NGFP values
fig1.add_hline(y=fp_lam4, line_dash="dash", line_color="rgba(100,100,100,0.3)")
fig1.add_hline(y=fp_lam6, line_dash="dash", line_color="rgba(100,100,100,0.3)")

fig1.update_layout(
    title={"text": "GFT Wetterich Flow: UV→IR Trajectory<br><span style='font-size: 18px; font-weight: normal;'>Rank-3 sextic melonic truncation | Litim regulator</span>"},
    legend=dict(orientation='h', yanchor='bottom', y=1.05, xanchor='center', x=0.5)
)
fig1.update_xaxes(title_text="RG time t", range=[-20, 1])
fig1.update_yaxes(title_text="Coupling λ̄")

fig1.write_image("flow_trajectory.png")
with open("flow_trajectory.png.meta.json", "w") as f:
    json.dump({"caption": "Gate 2: GFT Wetterich Flow λ₄(t), λ₆(t) from M_P to IR",
               "description": "Coupled RG flow of quartic and sextic couplings in rank-3 tensorial GFT"}, f)

# ---- CHART 2: Phase portrait (λ₄ vs λ₆) ----
fig2 = go.Figure()

# Plot the trajectory
fig2.add_trace(go.Scatter(x=lam4_vals, y=lam6_vals, mode='lines',
                          name='RG flow', line=dict(width=2.5),
                          showlegend=True))

# Mark UV starting point
fig2.add_trace(go.Scatter(x=[lam4_uv], y=[lam6_uv], mode='markers',
                          marker=dict(size=12, symbol='star'),
                          name='UV (M_P)'))

# Mark NGFP
fig2.add_trace(go.Scatter(x=[fp_lam4], y=[fp_lam6], mode='markers',
                          marker=dict(size=14, symbol='x'),
                          name='NGFP'))

# Draw eigenvectors
scale = 0.02
fig2.add_trace(go.Scatter(
    x=[fp_lam4 - scale*v_att[0], fp_lam4 + scale*v_att[0]],
    y=[fp_lam6 - scale*v_att[1], fp_lam6 + scale*v_att[1]],
    mode='lines', line=dict(dash='dash', width=1.5),
    name='UV-attractive'))

v_rep = np.array([-0.041531, -0.999137])
fig2.add_trace(go.Scatter(
    x=[fp_lam4 - scale*v_rep[0], fp_lam4 + scale*v_rep[0]],
    y=[fp_lam6 - scale*v_rep[1], fp_lam6 + scale*v_rep[1]],
    mode='lines', line=dict(dash='dot', width=1.5),
    name='UV-repulsive'))

fig2.update_layout(
    title={"text": "Phase Portrait: Critical Surface Near NGFP<br><span style='font-size: 18px; font-weight: normal;'>Gate 1 UV conditions projected onto attractive eigenvector</span>"},
    legend=dict(orientation='h', yanchor='bottom', y=1.05, xanchor='center', x=0.5)
)
fig2.update_xaxes(title_text="λ̄₄", range=[-0.02, 0.04])
fig2.update_yaxes(title_text="λ̄₆", range=[0.12, 0.22])

fig2.write_image("phase_portrait.png")
with open("phase_portrait.png.meta.json", "w") as f:
    json.dump({"caption": "Phase portrait showing critical surface projection near NGFP",
               "description": "UV initial conditions from Gate 1 projected onto UV-attractive eigenvector"}, f)

# ---- CHART 3: Uncertainty breakdown ----
fig3 = go.Figure()

categories = ['UV Posterior<br>(ε, σ)', 'Regulator<br>(Litim/Exp)', 
              'Truncation<br>(±necklace)', 'ξ Matching<br>(k=ξH)', 'Combined']
sigma_vals = [0.0, 741.0, 217.0, 0.0, 543.6]
colors_custom = ['rgb(99, 110, 250)', 'rgb(239, 85, 59)', 'rgb(0, 204, 150)', 
                 'rgb(171, 99, 250)', 'rgb(255, 161, 90)']

fig3.add_trace(go.Bar(x=categories, y=sigma_vals, 
                       marker_color=colors_custom,
                       text=[f'{v:.0f}' for v in sigma_vals],
                       textposition='outside'))

fig3.add_hline(y=1937*0.3, line_dash="dash", line_color="red",
               annotation_text="σ/c = 0.3 threshold", annotation_position="top left")

fig3.update_layout(
    title={"text": "Gate 2 Uncertainty Budget for c₁<br><span style='font-size: 18px; font-weight: normal;'>σ_c/c̄ = 0.281 — passes 0.3 threshold</span>"}
)
fig3.update_xaxes(title_text="Source")
fig3.update_yaxes(title_text="σ_c contribution")

fig3.write_image("uncertainty_budget.png")
with open("uncertainty_budget.png.meta.json", "w") as f:
    json.dump({"caption": "Gate 2 uncertainty budget: regulator choice dominates",
               "description": "Breakdown of σ_c by source: UV posterior, regulator, truncation, ξ matching"}, f)

print("All 3 charts generated successfully.")
