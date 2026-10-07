# exosome-retention-pk-model

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Preprint](https://img.shields.io/badge/Research_Square-rs.3.rs--7939584/v1-blue)](https://doi.org/10.21203/rs.3.rs-7939584/v1)

A Mechanistic Quantitative Systems Pharmacology (QSP) and Non-Linear Phase-Space Control Framework for Engineered Exosome Biodistribution & Retention PK/PD Dynamics.

---

## 📌 Overview

This repository contains the numerical implementations, dynamic simulations, and optimal control optimization scripts for the **Exosome Retention PK Model** based on *Shen's Laws of Biological Phase-Space Control*.

Key capabilities include:
- **Phase-Space Trajectory Modeling**: Simulating dynamic transitions across biological basins under multi-channel therapies.
- **RGD/Targeting Ligand Saturation & Retention**: Modeling tissue-specific binding kinetics ($\sigma_{\text{target}}$) and MPS clearance ($\sigma_{\text{clear}}$).
- **Optimal Dosing Optimization**: Implementing Pontryagin's Minimum Principle and SLSQP algorithm to generate front-loaded exponential decay dosing profiles ($\mathbf{U}^*(t) = \mathbf{U}_0 e^{-\lambda t}$).

---

## 🛠️ Installation & Requirements

Ensure you have Python 3.9+ installed. Clone the repository and install dependencies:

```bash
git clone [https://github.com/alexshen2468-ai/exosome-retention-pk-model.git](https://github.com/alexshen2468-ai/exosome-retention-pk-model.git)
cd exosome-retention-pk-model
pip install -r requirements.txt

Key Dependencies
 ⁠numpy >= 1.21.0⁠
 ⁠scipy >= 1.7.0⁠ (SLSQP, ⁠solve_ivp⁠)
 ⁠matplotlib >= 3.4.0⁠
 ⁠torch >= 1.9.0⁠ (Optional: for PINNs formulation)
🚀 Usage & Quick Start
1. Run Baseline Phase-Space Simulation
python scripts/run_phase_space_simulation.py --mode triple_channel
2. Dosing Schedule Optimization
python scripts/optimize_dosing_slsqp.py --target_efficiency 0.95

📑 Citation & References
If you use this model or codebase in your research, please cite our preprint:
Shen, Z. Algorithmization of Mesoscopic Bio-Assembly Informed by Room-Temperature Superconductivity Physics / Shen's Three Laws of Biological Phase-Space Control. Research Square (2026). DOI: 10.21203/rs.3.rs-7939584/v1

@article{shen2026exosome,
  title={Exosome retention PK model and non-linear phase-space control repository},
  author={Shen, Zhuofan},
  journal={Research Square Preprint},
  year={2026},
  doi={10.21203/rs.3.rs-7939584/v1}
}


