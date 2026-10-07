# exosome-retention-pk-model

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.xxxxxx.svg)](https://doi.org/10.5281/zenodo.xxxxxx)
[![Preprint](https://img.shields.io/badge/Research_Square-rs.3.rs--7939584/v1-blue)](https://doi.org/10.21203/rs.3.rs-7939584/v1)

A Quantitative Systems Pharmacology (QSP) and Non-Linear Phase-Space Control Framework for Engineered Exosome Biodistribution & Retention PK/PD Dynamics.

---

## 📌 Overview

This repository contains the numerical implementations, dynamic simulations, and optimal control optimization scripts for the **Exosome Retention PK Model** based on *Shen's Laws of Biological Phase-Space Control*.

Key capabilities include:
- **Phase-Space Trajectory Modeling**: Simulating dynamic transitions across biological basins under multi-channel therapies.
- **RGD/Targeting Ligand Saturation & Retention**: Modeling tissue-specific binding kinetics ($\sigma_{\text{target}}$) and MPS clearance ($\sigma_{\text{clear}}$).
- **Optimal Dosing Optimization**: Implementing Pontryagin's Minimum Principle and SLSQP algorithm to generate front-loaded exponential decay dosing profiles ($\mathbf{U}^*(t) = \mathbf{U}_0 e^{-\lambda t}$).

---

## 🛠️ Quick Start

### 1. Installation
```bash
git clone [https://github.com/alexshen2468-ai/exosome-retention-pk-model.git](https://github.com/alexshen2468-ai/exosome-retention-pk-model.git)
cd exosome-retention-pk-model
pip install -r requirements.txt

2. Run Simulations & Optimization
# Run baseline phase-space dynamic simulation
python run_phase_space_simulation.py

# Run SLSQP optimal dosing profile computation
python optimize_dosing_slsqp.py


📑 Citation
If you use this model or codebase in your research, please cite our preprint:

Shen, Z. Algorithmization of Mesoscopic Bio-Assembly Informed by Room-Temperature Superconductivity Physics / Shen's Three Laws of Biological Phase-Space Control. Research Square (2026). DOI: 10.21203/rs.3.rs-7939584/v1

