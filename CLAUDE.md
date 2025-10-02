# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains C implementations of Ising model simulations for statistical physics research. The Ising model is a mathematical model of ferromagnetism in statistical mechanics, used to study phase transitions and critical phenomena.

## Architecture

The codebase consists of six main simulation programs:

- **ising1d.c**: One-dimensional Ising model simulation with periodic boundary conditions
- **ising2d.c**: Two-dimensional Ising model simulation with periodic boundary conditions
- **ising2d_fss.c**: Two-dimensional finite-size scaling (FSS) analysis with Binder cumulant calculation
- **ising3d.c**: Three-dimensional Ising model simulation with periodic boundary conditions
- **ising3d_fss.c**: Three-dimensional finite-size scaling (FSS) analysis with Binder cumulant calculation

All implementations use Monte Carlo methods with the Metropolis algorithm to simulate spin configurations and calculate thermodynamic properties.

### Core Data Structure

All programs use a similar `Ising` struct but with different neighbor connectivity:
- 1D: Each spin has `left` and `right` neighbors (2 neighbors)
- 2D: Each spin has `up`, `down`, `left`, and `right` neighbors (4 neighbors)
- 3D: Each spin has `up`, `down`, `left`, `right`, `front`, and `back` neighbors (6 neighbors)

### Simulation Flow

1. **Initialization**: Random spin configuration with periodic boundary conditions
2. **Thermalization**: 200,000 Monte Carlo steps to reach equilibrium
3. **Measurement**: 200,000 Monte Carlo steps collecting statistics
4. **Analysis**: Calculate magnetization, energy, susceptibility (χ), and heat capacity

## Build and Run

### Random Number Generator

The code now uses **PCG (Permuted Congruential Generator)**, a modern high-quality random number generator that replaces the original closed-source `twist.c`. PCG is:

- Faster than Mersenne Twister
- Better statistical properties
- Widely used in scientific computing (NumPy's default since 2019)
- Perfect for Monte Carlo simulations

The implementation is in `pcg_random.h` with function aliases matching the original `twist.c` interface:
- `init_rnd(seed)`: Initialize with seed
- `drnd()`: Generate random double in [0,1)
- `gus()`: Generate seed from current time

### Compilation

```bash
gcc -o ising1d ising1d.c -lm
gcc -o ising2d ising2d.c -lm
gcc -o ising2d_fss ising2d_fss.c -lm
gcc -o ising3d ising3d.c -lm
gcc -o ising3d_fss ising3d_fss.c -lm
```

### Execution

```bash
# Run 1D simulation with N=100 spins
./ising1d 100

# Run 2D simulation with N×N=100×100 lattice
./ising2d 100

# Run 3D simulation with N×N×N=20×20×20 lattice
./ising3d 20

# Run 2D FSS analysis (outputs to stdout, redirect to file)
./ising2d_fss > fss_data.txt

# Run 3D FSS analysis (outputs to stdout, redirect to file)
./ising3d_fss > fss_data_3d.txt
```

### Output Format

**ising1d.c, ising2d.c, and ising3d.c** output tab-separated values for each temperature:
```
Temperature  Magnetization  Energy  Susceptibility  Heat_Capacity
```

**ising2d_fss.c and ising3d_fss.c** output tab-separated values with additional Binder cumulant:
```
N  Temperature  Magnetization  Energy  Susceptibility  Heat_Capacity  Binder_Cumulant
```

## Key Functions

### Shared Functions (1D, 2D, and 3D)
- `initialize()`: Sets up random initial configuration and neighbor connectivity
- `mcstep()`: Performs one Monte Carlo sweep using Metropolis algorithm
- `magnetization()`: Calculates absolute magnetization per spin
- `energy()`: Calculates energy per spin
- `cleanup()`: Frees allocated memory (in FSS programs)

### Random Number Functions (from pcg_random.h)
- `init_rnd(seed)`: Initialize PCG with seed
- `drnd()`: Generate random double in [0,1) using PCG
- `gus()`: Generate seed from current time
- `pcg32_random()`: Core PCG function generating 32-bit integers
- `pcg32_random_double()`: Convert to double precision

## Temperature Range

- **1D and 2D**: T = 3.01 → 0.01 in steps of 0.01
- **3D**: T = 3.01 → 0.01 in steps of 0.01 (basic simulation)
- **2D FSS**: T = 2.6 → 1.9 in steps of 0.005 (focused near T_c ≈ 2.269)
- **3D FSS**: T = 5.5 → 3.5 in steps of 0.02 (focused near T_c ≈ 4.51)

## Physics Notes

- Uses natural units where coupling constant J=1 and Boltzmann constant k_B=1
- Implements periodic boundary conditions for finite-size scaling studies
- Susceptibility and heat capacity are extensive quantities (scaled by system size)

### Critical Temperatures

- **1D Ising**: T_c = 0 (no phase transition at finite temperature)
- **2D Ising**: T_c ≈ 2.269 (exact: 2/ln(1+√2) ≈ 2.26918531...)
- **3D Ising**: T_c ≈ 4.5107 (exact value from high-precision simulations)

## Finite-Size Scaling (FSS) Analysis

The **ising2d_fss.c** and **ising3d_fss.c** programs implement advanced finite-size scaling analysis using the Binder cumulant crossing method to precisely determine the critical temperature T_c.

### Binder Cumulant

The Binder cumulant is a dimensionless quantity defined as:
```
U_L = 1 - <M⁴> / (3<M²>²)
```

At the critical temperature T_c, this quantity becomes size-independent. When plotted for different system sizes L, all curves intersect at a single point—the critical temperature.

### FSS Implementation Details

#### 2D FSS (ising2d_fss.c)
- **System Sizes**: L = 8, 16, 24, 32, 48, 64
- **Temperature Range**: T = 2.6 → 1.9 (focused near T_c ≈ 2.269)
- **Temperature Resolution**: ΔT = 0.005 (finer than basic simulation)
- **Thermalization**: 200,000 Monte Carlo steps per temperature
- **Measurement**: 200,000 Monte Carlo steps per temperature

#### 3D FSS (ising3d_fss.c)
- **System Sizes**: L = 4, 6, 8, 10, 12, 16
- **Temperature Range**: T = 5.5 → 3.5 (focused near T_c ≈ 4.51)
- **Temperature Resolution**: ΔT = 0.02
- **Thermalization**: 200,000 Monte Carlo steps per temperature
- **Measurement**: 200,000 Monte Carlo steps per temperature
- **Note**: Smaller sizes used for 3D due to N³ memory/computation scaling

### Visualization

#### 2D FSS Plots (plot_binder.py)

Generates two plots for 2D analysis:

1. **binder_crossing.png**: Binder cumulant crossing plot
   - Shows all system sizes crossing at T_c
   - Includes zoomed view near the critical point

2. **fss_complete_analysis.png**: Complete FSS analysis
   - Magnetization vs temperature for all L
   - Energy vs temperature for all L
   - Susceptibility vs temperature for all L
   - Heat capacity vs temperature for all L

#### 3D FSS Plots (plot_binder_3d.py)

Generates two plots for 3D analysis:

1. **binder_crossing_3d.png**: Binder cumulant crossing plot
   - Shows all system sizes crossing at T_c ≈ 4.51
   - Includes zoomed view near the critical point

2. **fss_complete_analysis_3d.png**: Complete FSS analysis
   - All thermodynamic quantities vs temperature for all L
   - Shows critical exponents: β ≈ 0.326, γ ≈ 1.237, ν ≈ 0.630

### Usage

```bash
# Run 2D FSS simulation
./ising2d_fss > fss_data.txt
python3 plot_binder.py

# Run 3D FSS simulation
./ising3d_fss > fss_data_3d.txt
python3 plot_binder_3d.py
```

## Testing

The **test_suite.c** includes comprehensive tests:

- PCG random number generator quality and reproducibility
- 1D and 2D Ising model physics (high-T and low-T behavior)
- 3D Ising model physics (high-T, low-T, and near-critical behavior)
- Binder cumulant calculation accuracy (2D and 3D)
- Finite-size scaling behavior

Run tests with:
```bash
gcc -o test_suite test_suite.c -lm
./test_suite
```

## Dimensionality Comparison

| Property | 1D | 2D | 3D |
|----------|----|----|-----|
| **Critical Temperature** | T_c = 0 | T_c ≈ 2.269 | T_c ≈ 4.51 |
| **Phase Transition** | None at T > 0 | Second-order | Second-order |
| **Critical Exponents** | N/A | β=1/8, γ=7/4, ν=1 | β≈0.326, γ≈1.237, ν≈0.630 |
| **Neighbors per Spin** | 2 | 4 | 6 |
| **Memory Scaling** | O(N) | O(N²) | O(N³) |
| **Universal Class** | N/A | 2D Ising | 3D Ising |

### Physical Insights

- **1D**: No finite-temperature phase transition (thermal fluctuations destroy order)
- **2D**: Exactly solvable (Onsager 1944); exhibits continuous phase transition
- **3D**: Most realistic model for real ferromagnets; belongs to same universality class as liquid-gas transitions

The Binder cumulant method works particularly well for 2D and 3D systems, where the crossing point precisely identifies T_c independent of system size at the critical point.