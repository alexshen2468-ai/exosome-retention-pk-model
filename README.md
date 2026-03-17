# Exosome Retention PK Model

![DOI](https://img.shields.io/badge/DOI-10.21203%2Frs.3.rs--7939584%2Fv1-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Python](https://img.shields.io/badge/python-3.x-blue)
![Status](https://img.shields.io/badge/status-research_prototype-orange)

Mechanistic pharmacokinetic model for exosome retention dynamics based on retention-ratio theory.

---

## Research Identity

This repository is maintained by **Zhuofan Shen**, an undergraduate researcher in biotechnology focusing on **computational modeling of exosome-mediated drug delivery systems**, with an emphasis on retention, biodistribution, and pharmacokinetics.

This work represents ongoing research in computational nanomedicine and exosome pharmacokinetics.

- ORCID: https://orcid.org/0009-0005-4304-8391  
- Preprint (Research Square): https://doi.org/10.21203/rs.3.rs-7939584/v1  
- GitHub Profile: https://github.com/alexshen2468-ai  

---

## Overview

Extracellular vesicles (exosomes) are increasingly studied as drug delivery vehicles due to their biocompatibility, intrinsic targeting capability, and ability to transport therapeutic molecules.

However, experimental biodistribution studies frequently report a discrepancy between:

- rapid systemic clearance from circulation  
- prolonged retention in organs such as the liver and spleen  

Traditional pharmacokinetic models typically describe systemic decay but often lack explicit mechanisms explaining **tissue retention dynamics**.

To address this, the present model introduces a **retention ratio parameter**.

---

## Model Description

The model implements a mechanistic pharmacokinetic (PK) structure including:

- systemic circulation compartment  
- organ distribution  
- reversible tissue binding dynamics  

Key kinetic processes include:

- systemic clearance  
- blood–organ exchange  
- reversible tissue binding  

The **retention ratio** is defined as:
R = k_bind / k_rel
This parameter represents the effective balance between tissue binding/uptake kinetics and release or clearance dynamics.

The framework enables exploration of how retention processes influence:

- organ accumulation  
- pharmacological duration  
- drug delivery efficiency  

---

## Computational Methods

Simulations include:

- ordinary differential equation (ODE) modeling  
- nonlinear parameter estimation  
- sensitivity analysis  
- retention-dynamics exploration  

All simulations are implemented in Python.

---

## Installation

Python 3.x is required.

Install required packages:

```bash
pip install numpy scipy matplotlib

Running the Model

Run the simulation script:
python exosome_pk_model.py
The script performs:
	•	time-course simulation of exosome distribution
	•	retention vs clearance analysis
	•	parameter sensitivity exploration

Reproducibility

The repository provides the full computational implementation used in the associated preprint.

Researchers can reproduce the simulations by running the provided Python script with the specified dependencies.

⸻

Biological Interpretation

The retention-ratio parameter captures the effective balance between tissue binding and clearance processes.

Biologically, it may reflect mechanisms influencing exosome biodistribution, including:
	•	membrane adhesion interactions
	•	receptor-mediated cellular uptake
	•	extracellular matrix trapping
	•	intracellular processing or degradation pathways
	Preprint

The full manuscript describing this model is available as a preprint:

Shen, Z. (2026)
Simulation-Based Exosome Nanocarrier Model for Targeted Opioid Analgesia
Research Square

DOI: https://doi.org/10.21203/rs.3.rs-7939584/v1
Citation

If you use this repository or model in academic work, please cite:

Shen, Z. (2026)
Simulation-Based Exosome Nanocarrier Model for Targeted Opioid Analgesia
Research Square preprint

https://doi.org/10.21203/rs.3.rs-7939584/v1

This repository also includes a machine-readable citation file (CITATION.cff).

⸻

Code Availability

The computational implementation used in this study is publicly available in this repository.

Repository Structure
exosome_pk_model.py                  # Main simulation script
README.md                            # Project documentation
LICENSE                              # MIT license
CITATION.cff                         # Machine-readable citation metadata
Retention Ratio PLOSONE v4 revised.docx   # Manuscript draft

Author

Zhuofan Shen
Undergraduate researcher in biotechnology

Research focus:
	•	Computational modeling of exosome-mediated drug delivery
	•	Retention dynamics and biodistribution
	•	Pharmacokinetic (PK/PD) modeling

Affiliation: University of Debrecen
	•	ORCID: https://orcid.org/0009-0005-4304-8391
	•	GitHub: https://github.com/alexshen2468-ai

⸻

Keywords

exosome · pharmacokinetics · drug delivery · nanomedicine · computational biology · PK modeling · biodistribution

License

This project is released under the MIT License.

⸻

Notes

Maintained as part of ongoing research in exosome pharmacokinetics and computational nanomedicine.