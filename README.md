# Exosome Retention PK Model

![DOI](https://img.shields.io/badge/DOI-10.21203%2Frs.3.rs--7939584%2Fv1-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Python](https://img.shields.io/badge/python-3.x-blue)
![Status](https://img.shields.io/badge/status-research%20prototype-orange)

Version: v1.0 (March 2026)

This repository contains the computational implementation of a mechanistic pharmacokinetic (PK) model designed to investigate **exosome biodistribution and tissue retention dynamics**.

The model introduces a **retention-ratio framework** describing the balance between tissue binding/uptake and release processes, offering a simplified mechanistic explanation for prolonged organ accumulation observed in experimental exosome studies.

---

# Overview

Extracellular vesicles (exosomes) are increasingly studied as drug delivery vehicles due to their biocompatibility, intrinsic targeting capability, and ability to transport therapeutic molecules.

However, experimental biodistribution studies frequently report a discrepancy between:

- rapid systemic clearance from circulation
- prolonged retention in organs such as the liver and spleen

Traditional pharmacokinetic models typically describe systemic decay but often lack explicit mechanisms explaining **tissue retention dynamics**.

To address this, the present model introduces a **retention ratio parameter**

\[
R = \frac{k_{bind}}{k_{rel}}
\]

which represents the effective balance between tissue binding/uptake kinetics and release or clearance dynamics.

This simplified framework enables exploration of how retention processes influence organ accumulation and pharmacological duration.

---

# Repository Structure
exosome-retention-pk-model
│
README.md
LICENSE
CITATION.cff
exosome_pk_model.py
Retention Ratio PLOSONE v4 revised.docx

**exosome_pk_model.py**

Main simulation script implementing the pharmacokinetic model.

**Retention Ratio PLOSONE v4 revised.docx**

Manuscript describing the modeling framework and computational analysis.

**CITATION.cff**

Machine-readable citation file enabling automatic citation generation on GitHub.

---

# Model Description

The model implements a mechanistic PK structure including:

- systemic circulation compartment
- organ distribution
- reversible tissue binding dynamics

Key kinetic processes include:

- systemic clearance
- blood–organ exchange
- reversible tissue binding

The **retention ratio**

\[
R = k_{bind}/k_{rel}
\]

acts as a simplified descriptor of tissue-level retention behavior.

---

# Computational Methods

Simulations include:

- ordinary differential equation (ODE) modeling
- nonlinear parameter estimation
- sensitivity analysis
- retention-dynamics exploration

All simulations were implemented in **Python**.

---

# Installation

Python 3.x is required.

Install required packages:
pip install numpy scipy matplotlib
---

# Running the Model

Run the simulation script:
python exosome_pk_model.py
The script performs pharmacokinetic simulations and generates retention-dynamics analysis.

---

# Reproducibility

The repository provides the full computational implementation used in the associated preprint.

Researchers can reproduce the simulations by running the provided Python script using the specified dependencies.

---

# Biological Interpretation

The retention-ratio parameter represents the effective balance between tissue binding/uptake processes and release or clearance dynamics.

Biologically, this parameter may capture mechanisms known to influence exosome biodistribution, including:

- membrane adhesion interactions
- receptor-mediated cellular uptake
- extracellular matrix trapping
- intracellular processing or degradation pathways

Experimental studies have shown that exosome accumulation in organs such as the liver and spleen is strongly influenced by these interactions.

The retention-ratio framework therefore provides a simplified mechanistic representation of these combined biological processes.

---

# Preprint

The full manuscript describing this model is available as a preprint:

Shen, Z. (2026)  
**Simulation-Based Exosome Nanocarrier Model for Targeted Opioid Analgesia**  
Research Square

DOI  
https://doi.org/10.21203/rs.3.rs-7939584/v1

---

# Citation

If you use this repository or model in academic work, please cite:

Shen, Z. (2026)  
Simulation-Based Exosome Nanocarrier Model for Targeted Opioid Analgesia  
Research Square preprint  
https://doi.org/10.21203/rs.3.rs-7939584/v1

This repository also includes a machine-readable citation file (`CITATION.cff`).

---

# Code Availability

The computational implementation used in this study is publicly available in this repository.

---

# Author

Zhuofan Shen  

ORCID  
https://orcid.org/0009-0005-4304-8391

---

# License

This project is released under the **MIT License**.

---

Maintained as part of ongoing research in exosome pharmacokinetics and computational nanomedicine.