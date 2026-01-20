# Tripartite Model of Human Mortality - Computational Validation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the computational validation code for the paper:

**"A Tripartite Model of Human Mortality: Cellular-Level Mathematical Framework for Understanding Death and the Theoretical Possibility of Radical Life Extension"**

Author: Mosab Hawarey (mhawarey@alumni.purdue.edu)

## Description

The code implements a Fokker-Planck equation solver to simulate cellular health distribution dynamics over a human lifespan. The model predicts mortality as a threshold phenomenon occurring when approximately 28% of cells fall below a critical health threshold.

## Key Features

- **Fokker-Planck PDE solver** with upwind finite difference scheme for numerical stability
- **Cellular health distribution evolution** from birth to 120 years
- **Reproduces all results** from the paper including Figure 1
- **Parameter specifications** matching published values

## Model Parameters

```python
alpha = 0.012      # Cellular damage rate (per year)
beta_0 = 0.008     # Initial repair rate (per year)
gamma = 0.025      # Repair decline rate (per year)
D = 0.012          # Diffusion coefficient
theta_c = 0.35     # Critical health threshold
f_crit = 0.28      # Critical fraction triggering death
```

## Requirements

```bash
numpy
matplotlib
```

## Installation

```bash
pip install numpy matplotlib
```

## Usage

```bash
python final_simulation.py
```

## Output

The script generates:
- **Figure 1**: Three-panel visualization showing:
  - Panel A: Cellular health distribution at 0, 30, 60, 90 years
  - Panel B: Mean cellular health decline over time
  - Panel C: Fraction of dysfunctional cells approaching mortality threshold

- **Console output**: Simulation statistics including:
  - Predicted lifespan
  - Initial and final mean health
  - Variance increase (should be ~1190%)

## Results

The simulation produces the following key findings:

- **Mean health decline**: 0.877 → 0.534 over 120 years
- **Variance amplification**: 1190% increase
- **Predicted lifespan**: ~85 years (matches human life expectancy)

## Mathematical Framework

The model uses the Fokker-Planck equation to describe cellular health distribution evolution:

```
∂ρ/∂t = -∂/∂x[v(x,t)ρ] + D ∂²ρ/∂x²
```

where:
- ρ(x,t) = probability density of cells at health state x and time t
- v(x,t) = -α(1-x) + β₀·x·exp(-γt) (drift term)
- D = diffusion coefficient

Death occurs when: ∫₀^θc ρ(x,t) dx ≥ f_crit

## Citation

If you use this code in your research, please cite:

```
Hawarey, M. (2026). A Tripartite Model of Human Mortality: Cellular-Level 
Mathematical Framework for Understanding Death and the Theoretical Possibility 
of Radical Life Extension. [Journal details pending publication]
```

## License

MIT License - see LICENSE file for details

## Contact

Mosab Hawarey - mhawarey@alumni.purdue.edu

## Acknowledgments

This computational framework validates the theoretical model presented in the paper and demonstrates that cellular health distribution dynamics can quantitatively predict human mortality patterns.
