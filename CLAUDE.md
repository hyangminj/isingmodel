# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains C implementations of Ising model simulations for statistical physics research. The Ising model is a mathematical model of ferromagnetism in statistical mechanics, used to study phase transitions and critical phenomena.

## Architecture

The codebase consists of three main simulation programs:

- **ising1d.c**: One-dimensional Ising model simulation with periodic boundary conditions
- **ising2d.c**: Two-dimensional Ising model simulation with periodic boundary conditions
- **ising2d_fss.c**: Two-dimensional finite-size scaling (FSS) analysis with Binder cumulant calculation

All implementations use Monte Carlo methods with the Metropolis algorithm to simulate spin configurations and calculate thermodynamic properties.

### Core Data Structure

Both programs use a similar `Ising` struct but with different neighbor connectivity:
- 1D: Each spin has `left` and `right` neighbors
- 2D: Each spin has `up`, `down`, `left`, and `right` neighbors

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
```

### Execution

```bash
# Run 1D simulation with N=100 spins
./ising1d 100

# Run 2D simulation with N×N=100×100 lattice
./ising2d 100

# Run FSS analysis (outputs to stdout, redirect to file)
./ising2d_fss > fss_data.txt
```

### Output Format

**ising1d.c and ising2d.c** output tab-separated values for each temperature:
```
Temperature  Magnetization  Energy  Susceptibility  Heat_Capacity
```

**ising2d_fss.c** outputs tab-separated values with additional Binder cumulant:
```
N  Temperature  Magnetization  Energy  Susceptibility  Heat_Capacity  Binder_Cumulant
```

## Key Functions

### Shared Functions (both 1D and 2D)
- `initialize()`: Sets up random initial configuration and neighbor connectivity
- `mcstep()`: Performs one Monte Carlo sweep using Metropolis algorithm
- `magnetization()`: Calculates absolute magnetization per spin
- `energy()`: Calculates energy per spin

### Random Number Functions (from pcg_random.h)
- `init_rnd(seed)`: Initialize PCG with seed
- `drnd()`: Generate random double in [0,1) using PCG
- `gus()`: Generate seed from current time
- `pcg32_random()`: Core PCG function generating 32-bit integers
- `pcg32_random_double()`: Convert to double precision

## Temperature Range

Both simulations sweep temperature from T=3.01 down to T=0.01 in steps of 0.01. This range spans across the critical temperature for phase transitions.

## Physics Notes

- Uses natural units where coupling constant J=1 and Boltzmann constant k_B=1
- Implements periodic boundary conditions for finite-size scaling studies
- Susceptibility and heat capacity are extensive quantities (scaled by system size)
- The 2D Ising model has a known critical temperature T_c ≈ 2.269 (exact: 2/ln(1+√2))

## Finite-Size Scaling (FSS) Analysis

The **ising2d_fss.c** program implements advanced finite-size scaling analysis using the Binder cumulant crossing method to precisely determine the critical temperature T_c.

### Binder Cumulant

The Binder cumulant is a dimensionless quantity defined as:
```
U_L = 1 - <M⁴> / (3<M²>²)
```

At the critical temperature T_c, this quantity becomes size-independent. When plotted for different system sizes L, all curves intersect at a single point—the critical temperature.

### FSS Implementation Details

- **System Sizes**: L = 8, 16, 24, 32, 48, 64
- **Temperature Range**: T = 2.6 → 1.9 (focused near T_c ≈ 2.269)
- **Temperature Resolution**: ΔT = 0.005 (finer than basic simulation)
- **Thermalization**: 200,000 Monte Carlo steps per temperature
- **Measurement**: 200,000 Monte Carlo steps per temperature

### Visualization

The **plot_binder.py** Python script generates two plots:

1. **binder_crossing.png**: Binder cumulant crossing plot
   - Shows all system sizes crossing at T_c
   - Includes zoomed view near the critical point

2. **fss_complete_analysis.png**: Complete FSS analysis
   - Magnetization vs temperature for all L
   - Energy vs temperature for all L
   - Susceptibility vs temperature for all L
   - Heat capacity vs temperature for all L

### Usage

```bash
# Run FSS simulation
./ising2d_fss > fss_data.txt

# Generate plots
python3 plot_binder.py
```

## Testing

The **test_suite.c** includes comprehensive tests:

- PCG random number generator quality and reproducibility
- Ising model physics (high-T and low-T behavior)
- Binder cumulant calculation accuracy
- Finite-size scaling behavior

Run tests with:
```bash
gcc -o test_suite test_suite.c -lm
./test_suite
```