#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

# Load data
data = np.loadtxt('fss_data.txt')

# Parse data by system size
sizes = np.unique(data[:, 0])
print(f"System sizes: {sizes}")

# Organize data by size
data_by_size = {}
for L in sizes:
    mask = data[:, 0] == L
    data_by_size[int(L)] = {
        'T': data[mask, 1],
        'mag': data[mask, 2],
        'energy': data[mask, 3],
        'chi': data[mask, 4],
        'C': data[mask, 5]
    }

# Find critical temperature from susceptibility maxima
Tc_estimates = []
chi_max_values = []

print("\nSusceptibility maxima:")
for L in sizes:
    L_int = int(L)
    chi = data_by_size[L_int]['chi']
    T = data_by_size[L_int]['T']

    # Find temperature where chi is maximum
    max_idx = np.argmax(chi)
    T_max = T[max_idx]
    chi_max = chi[max_idx]

    Tc_estimates.append(T_max)
    chi_max_values.append(chi_max)

    print(f"L={L_int:3d}: T_max = {T_max:.4f}, χ_max = {chi_max:.4f}")

# Finite-size scaling: chi_max ~ L^(gamma/nu) where gamma/nu = 7/4 for 2D Ising
# T_max(L) = Tc + a*L^(-1/nu) where 1/nu = 1 for 2D Ising
print(f"\nEstimated Tc from different system sizes: {Tc_estimates}")
print(f"Average Tc estimate: {np.mean(Tc_estimates):.4f} ± {np.std(Tc_estimates):.4f}")

# Linear fit of T_max vs L^(-1) to extrapolate Tc
L_inv = 1.0 / sizes
coeffs = np.polyfit(L_inv, Tc_estimates, 1)
Tc_extrapolated = coeffs[1]  # intercept at L -> infinity
print(f"\nExtrapolated Tc (linear fit): {Tc_extrapolated:.4f}")
print(f"Theoretical Tc (Onsager): {2/np.log(1+np.sqrt(2)):.4f}")

# Create comprehensive plots
fig = plt.figure(figsize=(15, 10))

# Plot 1: Susceptibility vs Temperature
ax1 = plt.subplot(2, 3, 1)
for L in sizes:
    L_int = int(L)
    T = data_by_size[L_int]['T']
    chi = data_by_size[L_int]['chi']
    ax1.plot(T, chi, 'o-', label=f'L={L_int}', markersize=3)
ax1.axvline(x=2/np.log(1+np.sqrt(2)), color='red', linestyle='--',
            label='Theoretical Tc', linewidth=2)
ax1.axvline(x=Tc_extrapolated, color='green', linestyle='--',
            label='Extrapolated Tc', linewidth=2)
ax1.set_xlabel('Temperature', fontsize=12)
ax1.set_ylabel('Susceptibility χ', fontsize=12)
ax1.set_title('Susceptibility vs Temperature', fontsize=14)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: Heat Capacity vs Temperature
ax2 = plt.subplot(2, 3, 2)
for L in sizes:
    L_int = int(L)
    T = data_by_size[L_int]['T']
    C = data_by_size[L_int]['C']
    ax2.plot(T, C, 'o-', label=f'L={L_int}', markersize=3)
ax2.axvline(x=2/np.log(1+np.sqrt(2)), color='red', linestyle='--',
            label='Theoretical Tc', linewidth=2)
ax2.axvline(x=Tc_extrapolated, color='green', linestyle='--',
            label='Extrapolated Tc', linewidth=2)
ax2.set_xlabel('Temperature', fontsize=12)
ax2.set_ylabel('Heat Capacity C', fontsize=12)
ax2.set_title('Heat Capacity vs Temperature', fontsize=14)
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: Magnetization vs Temperature
ax3 = plt.subplot(2, 3, 3)
for L in sizes:
    L_int = int(L)
    T = data_by_size[L_int]['T']
    mag = data_by_size[L_int]['mag']
    ax3.plot(T, mag, 'o-', label=f'L={L_int}', markersize=3)
ax3.axvline(x=2/np.log(1+np.sqrt(2)), color='red', linestyle='--',
            label='Theoretical Tc', linewidth=2)
ax3.axvline(x=Tc_extrapolated, color='green', linestyle='--',
            label='Extrapolated Tc', linewidth=2)
ax3.set_xlabel('Temperature', fontsize=12)
ax3.set_ylabel('Magnetization', fontsize=12)
ax3.set_title('Magnetization vs Temperature', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3)

# Plot 4: T_max vs 1/L (finite-size scaling for Tc)
ax4 = plt.subplot(2, 3, 4)
ax4.plot(L_inv, Tc_estimates, 'bo', markersize=8, label='Simulation data')
L_inv_fit = np.linspace(0, max(L_inv), 100)
T_fit = coeffs[0] * L_inv_fit + coeffs[1]
ax4.plot(L_inv_fit, T_fit, 'b--', label=f'Linear fit: Tc = {Tc_extrapolated:.4f}')
ax4.axhline(y=2/np.log(1+np.sqrt(2)), color='red', linestyle='--',
            label=f'Theoretical Tc = {2/np.log(1+np.sqrt(2)):.4f}', linewidth=2)
ax4.set_xlabel('1/L', fontsize=12)
ax4.set_ylabel('T_max (from χ peak)', fontsize=12)
ax4.set_title('Finite-Size Scaling: Critical Temperature', fontsize=14)
ax4.legend()
ax4.grid(True, alpha=0.3)

# Plot 5: chi_max vs L (scaling: chi_max ~ L^(gamma/nu))
ax5 = plt.subplot(2, 3, 5)
ax5.loglog(sizes, chi_max_values, 'ro', markersize=8, label='Simulation data')
# Fit: log(chi_max) = (gamma/nu) * log(L) + const
log_L = np.log(sizes)
log_chi_max = np.log(chi_max_values)
fit_coeffs = np.polyfit(log_L, log_chi_max, 1)
gamma_over_nu = fit_coeffs[0]
L_fit = np.logspace(np.log10(min(sizes)), np.log10(max(sizes)), 100)
chi_fit = np.exp(fit_coeffs[1]) * L_fit**gamma_over_nu
ax5.loglog(L_fit, chi_fit, 'r--',
          label=f'Power law: γ/ν = {gamma_over_nu:.3f} (theory: 1.75)')
ax5.set_xlabel('System size L', fontsize=12)
ax5.set_ylabel('χ_max', fontsize=12)
ax5.set_title('Finite-Size Scaling: Susceptibility', fontsize=14)
ax5.legend()
ax5.grid(True, alpha=0.3)

# Plot 6: Data collapse for susceptibility
ax6 = plt.subplot(2, 3, 6)
# Scale susceptibility by L^(-gamma/nu) and temperature by (T-Tc)*L^(1/nu)
# For 2D Ising: gamma/nu = 7/4, 1/nu = 1
gamma_over_nu_theory = 7.0/4.0
one_over_nu_theory = 1.0
Tc_use = Tc_extrapolated

for L in sizes:
    L_int = int(L)
    T = data_by_size[L_int]['T']
    chi = data_by_size[L_int]['chi']

    # Scaled quantities
    T_scaled = (T - Tc_use) * (L ** one_over_nu_theory)
    chi_scaled = chi / (L ** gamma_over_nu_theory)

    ax6.plot(T_scaled, chi_scaled, 'o', label=f'L={L_int}', markersize=4, alpha=0.7)

ax6.set_xlabel(r'$(T - T_c) L^{1/\nu}$', fontsize=12)
ax6.set_ylabel(r'$\chi L^{-\gamma/\nu}$', fontsize=12)
ax6.set_title('Data Collapse (Finite-Size Scaling)', fontsize=14)
ax6.legend()
ax6.grid(True, alpha=0.3)
ax6.set_xlim(-10, 10)

plt.tight_layout()
plt.savefig('fss_analysis.png', dpi=300, bbox_inches='tight')
print(f"\nPlot saved as 'fss_analysis.png'")

# Summary
print("\n" + "="*60)
print("FINITE-SIZE SCALING ANALYSIS SUMMARY")
print("="*60)
print(f"Extrapolated Tc:     {Tc_extrapolated:.4f}")
print(f"Theoretical Tc:      {2/np.log(1+np.sqrt(2)):.4f}")
print(f"Difference:          {abs(Tc_extrapolated - 2/np.log(1+np.sqrt(2))):.4f}")
print(f"γ/ν from fit:        {gamma_over_nu:.3f} (theory: 1.75)")
print("="*60)

plt.show()
