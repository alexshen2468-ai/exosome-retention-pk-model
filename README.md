README = r”””# Exosome Retention PK Model

[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.8+-blue)](https://www.python.org/)
[![Status](https://img.shields.io/badge/status-manuscript_submitted-orange)]()

Mechanistic pharmacokinetic model for exosome retention dynamics based on retention-ratio theory.

-----

## Overview

Extracellular vesicles (exosomes) are widely studied as drug delivery vehicles due to their biocompatibility, intrinsic targeting capability, and ability to transport therapeutic molecules. However, experimental studies consistently report a paradox:

- rapid systemic clearance (blood half-life ≈ 1–2 h)
- prolonged organ retention (liver accumulation persisting 24–168 h)

Traditional PK models capture clearance but cannot mechanistically account for this decoupling. This repository provides the full computational implementation of a **reversible organ reservoir model (M3)** that resolves this paradox through a single dimensionless parameter:

```
R = k_bind / k_rel
```

R represents the balance between tissue binding and release, and directly determines the bound fraction R/(1+R) at pseudo-steady state.

-----

## Key Results

|Metric                        |Value                   |
|------------------------------|------------------------|
|Fitted retention ratio R      |5.33                    |
|Hepatic bound fraction at PSS |~84%                    |
|Baseline plateau duration     |~60 h                   |
|AICc improvement over M2      |8.4 units (decisive)    |
|AICc improvement over M1      |19.6 units (decisive)   |
|Cross-validation MAE (blood)  |0.090–0.118             |
|Cross-validation MAE (liver)  |0.073–0.080             |
|Dominant sensitivity parameter|k_rel (Sobol S1 = 0.182)|

A novel **label persistence sensitivity analysis** (Fig 12) demonstrates that RMSE is minimised at zero label correction (α = 0), formally confirming that R correctly captures kinetic timescales independently of the ⁹⁹Zr radiolabel artefact.

-----

## Quick Start

```bash
pip install numpy scipy matplotlib SALib
python exosome_pk_model_v3_final.py
```

Generates all figures (Fig 1–12 and SFig 1–4) in the current directory. SALib is optional — if unavailable, Sobol analysis falls back to a built-in pure-numpy estimator automatically.

**Runtime:** approximately 5–10 minutes on a standard laptop (bootstrap n = 500 is the bottleneck).

-----

## Repository Structure

```
exosome_pk_model_v3_final.py          # Main simulation script (all analyses + figures)
mirzaaghasi_normalized.csv            # Cross-validation dataset 1 (Mirzaaghasi et al. 2021)
wiklander_normalized.csv              # Cross-validation dataset 2 (Wiklander et al. 2015)
requirements.txt                      # Python dependencies
README.md                             # This file
LICENSE                               # MIT licence
```

-----

## Figure Output

|File                       |Description                                 |Manuscript|
|---------------------------|--------------------------------------------|----------|
|Fig01.png                  |Blood exosome concentration fit             |Fig 1     |
|Fig02.png                  |Hepatic exosome concentration fit           |Fig 2     |
|Fig03.png                  |Intracellular uptake rate [illustrative]    |Fig 3     |
|Fig04.png                  |Local drug concentration [illustrative]     |Fig 4     |
|Fig05.png                  |Normalised therapeutic effect [illustrative]|Fig 5     |
|Fig06.png                  |Model comparison M1 / M2 / M3               |Fig 6     |
|Fig07.png                  |AICc bar chart                              |Fig 7     |
|Fig08.png                  |Bootstrap distribution of k_rel             |Fig 8     |
|Fig09.png                  |Sobol global sensitivity indices            |Fig 9     |
|Fig10.png                  |Model structure schematic                   |Fig 10    |
|Fig11.png                  |Plateau duration vs retention ratio R       |Fig 11    |
|Fig12.png                  |Label persistence sensitivity analysis      |Fig 12    |
|SFig01_hepatic_profiles.png|Hepatic profiles at R = 1, 3, 6, 10         |S Fig 1   |
|SFig02_2D_heatmap.png      |2D heatmap: plateau duration vs (R, k_rel)  |S Fig 2   |
|SFig03_crossval.png        |Independent cross-validation                |S Fig 3   |
|SFig04_krel_reduction.png  |k_rel reduction simulations                 |S Fig 4   |

-----

## Model Structure

```
M1  One-compartment (blood only, first-order clearance)
M2  Two-compartment (blood + peripheral); organ tracks blood
M3  Reservoir model (proposed):

  Blood (B) ──k_to──> Organ free (L_f) <──k_rel── Organ bound (L_b)
                                              ───k_bind──>
                        L_f ──Michaelis-Menten──> Intracellular (C_cell)
                        C_cell ──q──> Drug (C_drug)
                        Effect = C_drug / (C_drug + EC50)

  Retention ratio:  R = k_bind / k_rel
  Bound fraction:   R / (1 + R)
```

-----

## Datasets

|Dataset                |Source              |Exosome type|Label|Timepoints|Use             |
|-----------------------|--------------------|------------|-----|----------|----------------|
|Choi et al. 2022       |Pharmaceutics       |GMP-grade   |⁹⁹Zr |0.25–168 h|Primary fit     |
|Mirzaaghasi et al. 2021|Pharmaceutics       |HEK293T     |DiR  |1–8 h     |Cross-validation|
|Wiklander et al. 2015  |J Extracell Vesicles|C2C12       |DiI  |1–24 h    |Cross-validation|

Primary data (Choi et al.) are hard-coded in the script. Cross-validation data are provided as CSV files in the repository.

-----

## Computational Methods

- ODE integration: adaptive RK4(5), `scipy.integrate.solve_ivp` (rtol=1e⁻⁶, atol=1e⁻⁹)
- Parameter estimation: Levenberg-Marquardt, 10 random initialisations
- Model selection: bias-corrected AIC (AICc); ΔAICc > 10 = decisive evidence
- Uncertainty: non-parametric bootstrap resampling (n = 500)
- Sensitivity: Sobol variance decomposition via SALib (n = 1024 Saltelli samples, ±50% ranges)
- Label sensitivity: extended M3 with degraded-label pool; α scanned 0 → 1

-----

## Research Identity

This repository is maintained by **Zhuofan Shen**, an undergraduate researcher in biotechnology focusing on computational modeling of exosome-mediated drug delivery systems, with emphasis on retention, biodistribution, and pharmacokinetics.

- ORCID: https://orcid.org/0009-0005-4304-8391
- GitHub: https://github.com/alexshen2468-ai

-----

## Citation

If you use this code, please cite:

> Shen Z. Mechanistic modeling of exosome pharmacokinetics reveals a retention ratio governing sustained organ accumulation. *Manuscript submitted*, 2026.

-----

## Keywords

exosome · pharmacokinetics · drug delivery · nanomedicine · computational biology · PK/PD modeling · biodistribution · retention ratio · organ accumulation

-----

## License

MIT License — see <LICENSE> for details.
“””

with open(“README.md”, “w”, encoding=“utf-8”) as f:
f.write(README)

print(“README.md written successfully.”)
print(f”Total lines: {len(README.splitlines())}”)