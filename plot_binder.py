#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# Read the data from fss_data.txt
# Format: N, T, mag, energy, susceptibility, heat_capacity, binder_cumulant
data = np.loadtxt('fss_data.txt')

N_values = data[:, 0]
T_values = data[:, 1]
binder_values = data[:, 6]

# Get unique system sizes
unique_N = np.unique(N_values)

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Binder cumulant crossing
colors = plt.cm.viridis(np.linspace(0, 1, len(unique_N)))
for i, N in enumerate(unique_N):
    mask = N_values == N
    T = T_values[mask]
    binder = binder_values[mask]
    ax1.plot(T, binder, 'o-', label=f'L={int(N)}',
             color=colors[i], markersize=3, linewidth=1.5, alpha=0.8)

ax1.set_xlabel('Temperature (T)', fontsize=12)
ax1.set_ylabel('Binder Cumulant $U_L$', fontsize=12)
ax1.set_title('Binder Cumulant Crossing at $T_c$', fontsize=13, fontweight='bold')
ax1.legend(loc='best', fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.axvline(x=2.269, color='red', linestyle='--', linewidth=1, alpha=0.5, label='$T_c$ (exact)')

# Plot 2: Zoomed in near Tc
ax2_zoom = 0.15  # zoom range
Tc_approx = 2.269
for i, N in enumerate(unique_N):
    mask = N_values == N
    T = T_values[mask]
    binder = binder_values[mask]
    # Only plot points near Tc
    zoom_mask = np.abs(T - Tc_approx) < ax2_zoom
    ax2.plot(T[zoom_mask], binder[zoom_mask], 'o-', label=f'L={int(N)}',
             color=colors[i], markersize=4, linewidth=2, alpha=0.8)

ax2.set_xlabel('Temperature (T)', fontsize=12)
ax2.set_ylabel('Binder Cumulant $U_L$', fontsize=12)
ax2.set_title(f'Zoomed Near $T_c \\approx {Tc_approx}$', fontsize=13, fontweight='bold')
ax2.legend(loc='best', fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.axvline(x=2.269, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='$T_c$ (exact)')
ax2.set_xlim(Tc_approx - ax2_zoom, Tc_approx + ax2_zoom)

plt.tight_layout()
plt.savefig('binder_crossing.png', dpi=300, bbox_inches='tight')
print("Binder cumulant crossing plot saved as 'binder_crossing.png'")

# Also create a combined FSS plot
fig2, ((ax3, ax4), (ax5, ax6)) = plt.subplots(2, 2, figsize=(14, 10))

# Magnetization
for i, N in enumerate(unique_N):
    mask = N_values == N
    T = T_values[mask]
    mag = data[mask, 2]
    ax3.plot(T, mag, 'o-', label=f'L={int(N)}', color=colors[i], markersize=3, linewidth=1.5, alpha=0.8)
ax3.set_xlabel('Temperature (T)', fontsize=11)
ax3.set_ylabel('Magnetization |<M>|', fontsize=11)
ax3.set_title('Magnetization vs Temperature', fontsize=12, fontweight='bold')
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3)
ax3.axvline(x=2.269, color='red', linestyle='--', linewidth=1, alpha=0.5)

# Energy
for i, N in enumerate(unique_N):
    mask = N_values == N
    T = T_values[mask]
    energy = data[mask, 3]
    ax4.plot(T, energy, 'o-', label=f'L={int(N)}', color=colors[i], markersize=3, linewidth=1.5, alpha=0.8)
ax4.set_xlabel('Temperature (T)', fontsize=11)
ax4.set_ylabel('Energy per spin', fontsize=11)
ax4.set_title('Energy vs Temperature', fontsize=12, fontweight='bold')
ax4.legend(fontsize=9)
ax4.grid(True, alpha=0.3)
ax4.axvline(x=2.269, color='red', linestyle='--', linewidth=1, alpha=0.5)

# Susceptibility
for i, N in enumerate(unique_N):
    mask = N_values == N
    T = T_values[mask]
    chi = data[mask, 4]
    ax5.plot(T, chi, 'o-', label=f'L={int(N)}', color=colors[i], markersize=3, linewidth=1.5, alpha=0.8)
ax5.set_xlabel('Temperature (T)', fontsize=11)
ax5.set_ylabel('Susceptibility $\\chi$', fontsize=11)
ax5.set_title('Susceptibility vs Temperature', fontsize=12, fontweight='bold')
ax5.legend(fontsize=9)
ax5.grid(True, alpha=0.3)
ax5.axvline(x=2.269, color='red', linestyle='--', linewidth=1, alpha=0.5)

# Heat Capacity
for i, N in enumerate(unique_N):
    mask = N_values == N
    T = T_values[mask]
    capacity = data[mask, 5]
    ax6.plot(T, capacity, 'o-', label=f'L={int(N)}', color=colors[i], markersize=3, linewidth=1.5, alpha=0.8)
ax6.set_xlabel('Temperature (T)', fontsize=11)
ax6.set_ylabel('Heat Capacity $C$', fontsize=11)
ax6.set_title('Heat Capacity vs Temperature', fontsize=12, fontweight='bold')
ax6.legend(fontsize=9)
ax6.grid(True, alpha=0.3)
ax6.axvline(x=2.269, color='red', linestyle='--', linewidth=1, alpha=0.5)

plt.tight_layout()
plt.savefig('fss_complete_analysis.png', dpi=300, bbox_inches='tight')
print("Complete FSS analysis plot saved as 'fss_complete_analysis.png'")

plt.show()
