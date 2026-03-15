# Exosome Retention PK Model

Mechanistic pharmacokinetic model for exosome retention dynamics based on **retention ratio theory**.

This repository provides the computational implementation, simulation results, and manuscript associated with the study:

**Simulation-Based Exosome Nanocarrier Model for Targeted Opioid Analgesia**

Preprint available on Research Square:  
https://doi.org/10.21203/rs.3.rs-7939584/v1

---

# Overview

Exosomes are increasingly investigated as nanocarriers for targeted drug delivery.  
However, quantitative modeling of **retention vs clearance dynamics** remains limited.

This project introduces a **retention ratio–based pharmacokinetic framework** that models how binding and release dynamics influence sustained accumulation of exosome carriers in biological tissues.

The computational model integrates:

- ODE-based pharmacokinetic modeling
- parameter estimation
- global sensitivity analysis
- uncertainty quantification
- cross-dataset validation

---

# Repository Structure
exosome-retention-pk-model
│
├ README.md
├ LICENSE
├ exosome-pk-model.py
├ Python Exosome Pk Model v4.pdf
└ Retention Ratio PLOSONE v4 revised.docx

**File descriptions**

| File | Description |
|----|----|
| exosome-pk-model.py | Python implementation of the PK simulation model |
| Python Exosome Pk Model v4.pdf | Technical description of the computational model |
| Retention Ratio PLOSONE v4 revised.docx | Manuscript describing the retention ratio framework |
| LICENSE | Open-source license |

---

# Methods

The computational workflow includes:

- ODE-based pharmacokinetic modeling
- parameter fitting
- Sobol global sensitivity analysis
- bootstrap uncertainty analysis
- cross-dataset validation

The retention ratio framework describes the balance between:

**binding dynamics vs release kinetics** governing sustained tissue accumulation of exosome carriers.

---

# Usage

To reproduce the pharmacokinetic simulation:
python3 exosome-pk-model.py

---

# Requirements

Recommended environment:

Python 3.9+
numpy
scipy
matplotlib
pandas
SALib

Install dependencies:

pip install numpy scipy matplotlib pandas SALib

---

# Scientific Context

This project contributes to research in:

- exosome engineering
- nanomedicine pharmacokinetics
- targeted drug delivery systems
- computational biotechnology

The retention ratio model aims to provide a **conceptual and computational framework** for analyzing the pharmacokinetics of engineered extracellular vesicle systems.

---

# Author

**Zhuofan Shen**

Research focus:

- exosome engineering
- computational biotechnology
- nanocarrier pharmacokinetics

ORCID  
https://orcid.org/0009-0005-4304-8391

---

# Citation

If you use this code or model, please cite:

Shen Z.  
Simulation-Based Exosome Nanocarrier Model for Targeted Opioid Analgesia.  
Research Square Preprint.  
https://doi.org/10.21203/rs.3.rs-7939584/v1

---

# License

This project is released under the MIT License.
