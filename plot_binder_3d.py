#!/usr/bin/env python3
"""
3D Ising Model Finite-Size Scaling Analysis Visualization

This script visualizes the finite-size scaling (FSS) analysis of the 3D Ising model,
including Binder cumulant crossing to determine the critical temperature.

Usage:
    python3 plot_binder_3d.py [input_file]

If no input file is specified, it reads from 'fss_data_3d.txt'
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

# Read data file
input_file = sys.argv[1] if len(sys.argv) > 1 else 'fss_data_3d.txt'

try:
    data = np.loadtxt(input_file)
except FileNotFoundError:
    print(f"Error: File '{input_file}' not found.")
    print("Please run './ising3d_fss > fss_data_3d.txt' first")
    sys.exit(1)

# Extract columns
N_vals = data[:, 0]
T_vals = data[:, 1]
M_vals = data[:, 2]
E_vals = data[:, 3]
chi_vals = data[:, 4]
C_vals = data[:, 5]
U_vals = data[:, 6]

# Get unique system sizes
sizes = np.unique(N_vals)
print(f"System sizes: {sizes}")

# Organize data by system size
data_by_size = {}
for size in sizes:
    mask = N_vals == size
    data_by_size[size] = {
        'T': T_vals[mask],
        'M': M_vals[mask],
        'E': E_vals[mask],
        'chi': chi_vals[mask],
        'C': C_vals[mask],
        'U': U_vals[mask]
    }

# ============================================================================
# Figure 1: Binder Cumulant Crossing (Main result for Tc determination)
# ============================================================================

fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1a: Full Binder cumulant curves
colors = plt.cm.viridis(np.linspace(0, 0.9, len(sizes)))
for i, size in enumerate(sizes):
    ax1.plot(data_by_size[size]['T'], data_by_size[size]['U'],
             'o-', color=colors[i], label=f'L={int(size)}', markersize=4)

ax1.axhline(y=0.61, color='red', linestyle='--', alpha=0.3, label='Universal value ≈ 0.61')
ax1.set_xlabel('Temperature (T)', fontsize=12)
ax1.set_ylabel('Binder Cumulant (U)', fontsize=12)
ax1.set_title('3D Ising Model: Binder Cumulant vs Temperature', fontsize=13, fontweight='bold')
ax1.legend(loc='best', fontsize=9)
ax1.grid(True, alpha=0.3)

# Plot 1b: Zoomed view near critical point
T_min, T_max = 4.3, 4.7
for i, size in enumerate(sizes):
    mask = (data_by_size[size]['T'] >= T_min) & (data_by_size[size]['T'] <= T_max)
    ax2.plot(data_by_size[size]['T'][mask], data_by_size[size]['U'][mask],
             'o-', color=colors[i], label=f'L={int(size)}', markersize=5, linewidth=2)

ax2.axhline(y=0.61, color='red', linestyle='--', alpha=0.3, label='Universal value')
ax2.axvline(x=4.51, color='gray', linestyle=':', alpha=0.5, label='Known Tc ≈ 4.51')
ax2.set_xlabel('Temperature (T)', fontsize=12)
ax2.set_ylabel('Binder Cumulant (U)', fontsize=12)
ax2.set_title('Zoomed Near Critical Point', fontsize=13, fontweight='bold')
ax2.legend(loc='best', fontsize=9)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('binder_crossing_3d.png', dpi=300, bbox_inches='tight')
print("Saved: binder_crossing_3d.png")

# ============================================================================
# Figure 2: Complete FSS Analysis (All thermodynamic quantities)
# ============================================================================

fig2, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 2a: Magnetization
ax = axes[0, 0]
for i, size in enumerate(sizes):
    ax.plot(data_by_size[size]['T'], data_by_size[size]['M'],
            'o-', color=colors[i], label=f'L={int(size)}', markersize=3)
ax.axvline(x=4.51, color='gray', linestyle=':', alpha=0.5, label='Tc ≈ 4.51')
ax.set_xlabel('Temperature (T)', fontsize=11)
ax.set_ylabel('Magnetization per spin |M|', fontsize=11)
ax.set_title('(a) Magnetization', fontsize=12, fontweight='bold')
ax.legend(loc='best', fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 2b: Energy
ax = axes[0, 1]
for i, size in enumerate(sizes):
    ax.plot(data_by_size[size]['T'], data_by_size[size]['E'],
            'o-', color=colors[i], label=f'L={int(size)}', markersize=3)
ax.axvline(x=4.51, color='gray', linestyle=':', alpha=0.5, label='Tc ≈ 4.51')
ax.set_xlabel('Temperature (T)', fontsize=11)
ax.set_ylabel('Energy per spin (E)', fontsize=11)
ax.set_title('(b) Energy', fontsize=12, fontweight='bold')
ax.legend(loc='best', fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 2c: Susceptibility
ax = axes[1, 0]
for i, size in enumerate(sizes):
    ax.plot(data_by_size[size]['T'], data_by_size[size]['chi'],
            'o-', color=colors[i], label=f'L={int(size)}', markersize=3)
ax.axvline(x=4.51, color='gray', linestyle=':', alpha=0.5, label='Tc ≈ 4.51')
ax.set_xlabel('Temperature (T)', fontsize=11)
ax.set_ylabel('Susceptibility (χ)', fontsize=11)
ax.set_title('(c) Magnetic Susceptibility', fontsize=12, fontweight='bold')
ax.legend(loc='best', fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 2d: Heat Capacity
ax = axes[1, 1]
for i, size in enumerate(sizes):
    ax.plot(data_by_size[size]['T'], data_by_size[size]['C'],
            'o-', color=colors[i], label=f'L={int(size)}', markersize=3)
ax.axvline(x=4.51, color='gray', linestyle=':', alpha=0.5, label='Tc ≈ 4.51')
ax.set_xlabel('Temperature (T)', fontsize=11)
ax.set_ylabel('Heat Capacity (C)', fontsize=11)
ax.set_title('(d) Heat Capacity', fontsize=12, fontweight='bold')
ax.legend(loc='best', fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('fss_complete_analysis_3d.png', dpi=300, bbox_inches='tight')
print("Saved: fss_complete_analysis_3d.png")

# ============================================================================
# Analysis Summary
# ============================================================================

print("\n" + "="*60)
print("3D ISING MODEL FINITE-SIZE SCALING ANALYSIS")
print("="*60)
print(f"Known critical temperature: Tc ≈ 4.5107 (Exact)")
print(f"System sizes analyzed: {[int(s) for s in sizes]}")
print(f"Temperature range: {T_vals.min():.2f} - {T_vals.max():.2f}")

# Estimate Tc from Binder cumulant crossing
# Find temperature where U ≈ 0.61 for largest system
largest_size = sizes[-1]
U_data = data_by_size[largest_size]['U']
T_data = data_by_size[largest_size]['T']

# Find closest point to U = 0.61
idx = np.argmin(np.abs(U_data - 0.61))
estimated_Tc = T_data[idx]

print(f"\nEstimated Tc from Binder cumulant (L={int(largest_size)}): {estimated_Tc:.3f}")
print(f"Universal Binder cumulant value: U* ≈ 0.61")

# Show peak positions in susceptibility (another way to estimate Tc)
print("\nSusceptibility peak positions:")
for size in sizes:
    chi_data = data_by_size[size]['chi']
    T_data = data_by_size[size]['T']
    peak_idx = np.argmax(chi_data)
    T_peak = T_data[peak_idx]
    print(f"  L={int(size):2d}: T_peak = {T_peak:.3f}")

print("\n" + "="*60)
print("Physical Interpretation:")
print("="*60)
print("- The 3D Ising model exhibits second-order phase transition")
print("- At T < Tc: Ordered phase (ferromagnetic, M ≠ 0)")
print("- At T > Tc: Disordered phase (paramagnetic, M → 0)")
print("- Binder cumulant curves cross at Tc (size-independent)")
print("- Critical exponents: β ≈ 0.326, γ ≈ 1.237, ν ≈ 0.630")
print("- 3D Ising model is in the same universality class as")
print("  liquid-gas transitions and many magnetic systems")
print("="*60)

plt.show()
