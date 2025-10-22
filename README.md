# White Worms: Modeling Defenses for Networked Systems in IoT Security

This repository contains the code for reproducing the analyses presented in:

**"Modelling White Worms as Defenses for Networked Systems in IoT Security"**  
*Andreia Sofia Teixeira, Ignacio Echegoyen, Rasha Shanaz, and Alberto Aleta*  
arXiv: ...

This work investigates the spreading dynamics of two competing worms—a malicious black worm and a benign white worm—on networked IoT devices using different network topologies (Erdős-Rényi, Scale-Free, Complete Graph, and Stochastic Block Model).

This project was quickstarted during the Complexity72h workshop held at IFISC in Palma, Spain, June 26-30, 2023, with the participation of Francesco Bonacina, Diego Escribano, Marcus Krellner and Francesco Paolo Nerini who we want to thank for their contributions. More information: https://www.complexity72h.com

## Installation

1. Clone the repository:
```bash
git clone https://github.com/aaleta/whiteworms.git
cd whiteworms
```

2. Create a virtual environment (recommended):
```bash
python3.10 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Project Structure

- `main.py` - Main script for running simulations
- `parser.py` - Command-line argument parser
- `deterministic/` - Deterministic approximation methods
- `stochastic/` - Gillespie algorithm implementation
- `networks/` - Network topology files (edgelist format)
- `notebooks/` - Jupyter notebooks for analysis and visualization
- `results/` - Output directory for simulation results

## Usage

The main simulation script supports three types of analysis:

### 1. Protected Nodes Analysis
Estimates the fraction of nodes protected by the white worm:
```bash
python main.py -a protected -n ER_N10000_k10 -bb 0.5 -bw 0.5 -e 0.1 -g 0.1 -m 0.1 -i 100
```

### 2. Botnet Size Analysis
Estimates the final botnet size (nodes infected by black worm):
```bash
python main.py -a botnet -n SF_N10000_k10 -bb 0.5 -bw 0.5 -e 0.1 -g 0.1 -m 0.1 -i 100
```

### 3. Botnet Threshold Analysis
Estimates botnet size when a threshold fraction is infected:
```bash
python main.py -a threshold -n SBM_N10000_k10 -bb 0.5 -bw 0.5 -e 0.1 -g 0.1 -m 0.1 -th 0.05 -i 100
```

### Parameters

- `-a, --analysis-type`: Type of analysis (`protected`, `botnet`, or `threshold`)
- `-n, --network`: Network file name (without `.edgelist` extension)
- `-bb, --beta_b`: Black worm transmission rate
- `-bw, --beta_w`: White worm transmission rate
- `-e, --epsilon`: White worm spontaneous removal rate
- `-g, --gamma`: Prompted transition rate
- `-m, --mu`: Forced fix rate
- `-i, --iterations`: Number of simulation iterations (default: 100)
- `-nb, --n_black`: Initial number of black worm seeds (default: 1)
- `-nw, --n_white`: Initial number of white worm seeds (default: 1)
- `-th, --threshold`: Fraction threshold for botnet analysis (default: 0)

### Available Networks

- `ER_N10000_k10` - Erdős-Rényi random graph (N=10,000, ⟨k⟩=10)
- `SF_N10000_k10` - Scale-Free network (N=10,000, ⟨k⟩=10)
- `SBM_N10000_k10` - Stochastic Block Model (N=10,000, ⟨k⟩=10)
- `CG_N100` - Complete Graph (N=100)
- `CG_N1000` - Complete Graph (N=1,000)

**Note:** The complete graph for N=10,000 exceeds GitHub's file size limits and must be generated locally using the notebook `notebooks/generate_networks.ipynb`.

## Notebooks

The `notebooks/` directory contains Jupyter notebooks for:
- Generating network topologies
- Running deterministic approximations
- Analyzing and visualizing results

## Results

Simulation results are saved as pickle files in the `results/` directory with filenames encoding the parameters used.

## Citation

If you use this code in your research, please cite:
```bibtex
@article{teixeira2023modelling,
  title={Modelling White Worms as Defenses for Networked Systems in IoT Security},
  author={Teixeira, Andreia Sofia and Echegoyen, Ignacio and Shanaz, Rasha and Aleta, Alberto},
  journal={...},
  year={2025}
}
```