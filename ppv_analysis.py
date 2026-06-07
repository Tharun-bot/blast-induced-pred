"""
PPV Attenuation Analysis — Blast Monitoring Dataset
====================================================
Columns in blast_monitoring_data.csv:
  Point   : monitoring point number
  D_m     : distance from blast (m)
  W_kg    : charge per delay (kg)
  SD      : scaled distance m/kg^0.5  (provided; we also recompute as SD_check)
  PPV_mms : resultant PPV (mm/s)
  PPV_V   : vertical PPV component
  PPV_L   : longitudinal PPV component
  PPV_T   : transverse PPV component

Day 1 — Clean + log-transform + scatter plot
Day 2 — Holmberg-Persson empirical fit  (PPV = K × SD^-alpha)
Bonus  — Component breakdown plot
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import json

# ─────────────────────────────────────────────────────────────────────────────
# DAY 1 — Load & Clean
# ─────────────────────────────────────────────────────────────────────────────

df = pd.read_csv("blast_monitoring_data.csv")

# Rename columns to standard names matching the actual file structure
df.columns = ['Point', 'D_m', 'W_kg', 'SD', 'PPV_mms', 'PPV_V', 'PPV_L', 'PPV_T']

# Recompute SD from raw columns (serves as a cross-check)
df['SD_check'] = df['D_m'] / np.sqrt(df['W_kg'])

# Quality filters
df = df[df['PPV_mms'] > 0.1]    # remove zero / noise readings
df = df[df['PPV_mms'] < 500]    # remove instrument spikes
df = df[df['D_m'] > 5]          # too-close readings are unreliable

# Log-transform columns for regression diagnostics
df['Log_PPV'] = np.log10(df['PPV_mms'])
df['Log_SD']  = np.log10(df['SD'])

df.to_csv("master_ppv_clean.csv", index=False)
print(f"Clean dataset: {df.shape[0]} rows × {df.shape[1]} cols")

# Day 1 plot — log-log scatter (should show clear downward trend)
plt.figure(figsize=(7, 5))
plt.scatter(df['Log_SD'], df['Log_PPV'],
            alpha=0.6, edgecolors='k', s=40, color='steelblue')
plt.xlabel('Log₁₀(Scaled Distance)  [m/kg⁰·⁵]', fontsize=12)
plt.ylabel('Log₁₀(PPV)  [mm/s]', fontsize=12)
plt.title('PPV Attenuation — Blast Monitoring Dataset', fontsize=13)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('ppv_scatter.png', dpi=300, bbox_inches='tight')
plt.show()
print("Saved: ppv_scatter.png")

# ─────────────────────────────────────────────────────────────────────────────
# DAY 2 — Holmberg-Persson Empirical Fit
# Formula: PPV = K × (D / W^0.5)^(-alpha)  ≡  K × SD^(-alpha)
# ─────────────────────────────────────────────────────────────────────────────

def holmberg_persson(SD, K, alpha):
    return K * (SD ** (-alpha))

popt, _ = curve_fit(
    holmberg_persson, df['SD'], df['PPV_mms'],
    p0=[1000, 1.5], maxfev=5000
)
K, alpha = popt

# Goodness-of-fit (R²)
PPV_pred = holmberg_persson(df['SD'], K, alpha)
ss_res   = np.sum((df['PPV_mms'] - PPV_pred) ** 2)
ss_tot   = np.sum((df['PPV_mms'] - df['PPV_mms'].mean()) ** 2)
r2       = 1 - ss_res / ss_tot

print(f"\nHolmberg-Persson fit:")
print(f"  K     = {K:.2f}")
print(f"  alpha = {alpha:.4f}")
print(f"  R²    = {r2:.4f}")

# Save parameters for use in UDEC / further modelling
with open("hp_params.json", "w") as f:
    json.dump({"K": round(K, 4), "alpha": round(alpha, 4), "R2": round(r2, 4)}, f, indent=2)
print("Saved: hp_params.json")

# H-P fit plot
SD_range  = np.linspace(df['SD'].min(), df['SD'].max(), 300)
PPV_curve = holmberg_persson(SD_range, K, alpha)

plt.figure(figsize=(7, 5))
plt.scatter(df['SD'], df['PPV_mms'],
            alpha=0.6, edgecolors='k', s=40, color='steelblue', label='Measured PPV')
plt.plot(SD_range, PPV_curve, 'r-', linewidth=2,
         label=f'H-P Fit:  K={K:.1f},  α={alpha:.3f}  (R²={r2:.3f})')
plt.xlabel('Scaled Distance  [m/kg⁰·⁵]', fontsize=12)
plt.ylabel('PPV  [mm/s]', fontsize=12)
plt.title('Holmberg-Persson Attenuation Fit', fontsize=13)
plt.legend(fontsize=10)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('hp_fit.png', dpi=300, bbox_inches='tight')
plt.show()
print("Saved: hp_fit.png")

# ─────────────────────────────────────────────────────────────────────────────
# BONUS — PPV Component Breakdown (Vertical / Longitudinal / Transverse)
# Useful for directional anisotropy analysis before UDEC input prep
# ─────────────────────────────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(8, 5))
ax.scatter(df['SD'], df['PPV_V'], s=25, alpha=0.5, label='Vertical',     color='royalblue')
ax.scatter(df['SD'], df['PPV_L'], s=25, alpha=0.5, label='Longitudinal', color='darkorange')
ax.scatter(df['SD'], df['PPV_T'], s=25, alpha=0.5, label='Transverse',   color='green')
ax.set_xlabel('Scaled Distance  [m/kg⁰·⁵]', fontsize=12)
ax.set_ylabel('PPV Component  [mm/s]', fontsize=12)
ax.set_title('PPV Component Breakdown vs Scaled Distance', fontsize=13)
ax.legend(fontsize=10)
ax.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('ppv_components.png', dpi=300, bbox_inches='tight')
plt.show()
print("Saved: ppv_components.png")

# =====================================================
# DAY 3 — DIRECTIONAL ATTENUATION ANALYSIS
# =====================================================

import os

os.makedirs("figures", exist_ok=True)

# -----------------------------------------------------
# Derived variables
# -----------------------------------------------------

df['logSD'] = np.log10(df['SD'])

df['ratio_VL'] = df['PPV_V'] / df['PPV_L']
df['ratio_VT'] = df['PPV_V'] / df['PPV_T']
df['ratio_LT'] = df['PPV_L'] / df['PPV_T']

# =====================================================
# FIGURE 1 — 4 PANEL ATTENUATION PLOT
# =====================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

targets = ['PPV_mms', 'PPV_V', 'PPV_L', 'PPV_T']

titles = [
    'Resultant PPV',
    'Vertical (V)',
    'Longitudinal (L)',
    'Transverse (T)'
]

print("\n========== FIGURE 1 SLOPES ==========")

for ax, col, title in zip(axes.flat, targets, titles):

    y = np.log10(df[col])

    ax.scatter(
        df['logSD'],
        y,
        alpha=0.6,
        s=30
    )

    m, c = np.polyfit(df['logSD'], y, 1)

    ax.plot(
        df['logSD'],
        m * df['logSD'] + c,
        'r-',
        lw=2
    )

    print(f"{title}: slope = {m:.4f}")

    ax.set_xlabel("log10(Scaled Distance)")
    ax.set_ylabel(f"log10({title})")
    ax.set_title(f"{title}\nSlope = {m:.3f}")
    ax.grid(True, alpha=0.3)

plt.tight_layout()

plt.savefig(
    "figures/fig1_attenuation_4panel.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Saved: figures/fig1_attenuation_4panel.png")

# =====================================================
# FIGURE 2 — COMPONENT RATIOS
# =====================================================

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

pairs = [
    ('ratio_VL', 'V/L'),
    ('ratio_VT', 'V/T'),
    ('ratio_LT', 'L/T')
]

print("\n========== FIGURE 2 SLOPES ==========")

for ax, (col, label) in zip(axes, pairs):

    ax.scatter(
        df['SD'],
        df[col],
        alpha=0.5,
        s=25
    )

    ax.axhline(
        1.0,
        color='red',
        linestyle='--',
        label='ratio = 1'
    )

    m, c = np.polyfit(df['SD'], df[col], 1)

    x = np.linspace(
        df['SD'].min(),
        df['SD'].max(),
        100
    )

    ax.plot(
        x,
        m * x + c,
        'g-',
        lw=2
    )

    print(f"{label}: slope = {m:.6f}")

    ax.set_xlabel("Scaled Distance")
    ax.set_ylabel(label)
    ax.set_title(f"{label} Ratio")
    ax.legend()

plt.tight_layout()

plt.savefig(
    "figures/fig2_component_ratios.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Saved: figures/fig2_component_ratios.png")

print("\nDay 3 Complete")

# =====================================================
# DAY 4 — HP fitting per component
# =====================================================

