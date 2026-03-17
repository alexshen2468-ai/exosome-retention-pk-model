# Exosome Retention PK Model

![DOI](https://img.shields.io/badge/DOI-10.21203%2Frs.3.rs--7939584%2Fv1-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Python](https://img.shields.io/badge/python-3.x-blue)
![Status](https://img.shields.io/badge/status-research_prototype-orange)

Mechanistic pharmacokinetic model for exosome retention dynamics based on retention-ratio theory.

---

## Research Identity

This repository is maintained by **Zhuofan Shen**, an undergraduate researcher in biotechnology focusing on **computational modeling of exosome-mediated drug delivery systems**, with emphasis on retention, biodistribution, and pharmacokinetics.

- ORCID: https://orcid.org/0009-0005-4304-8391  
- Preprint: https://doi.org/10.21203/rs.3.rs-7939584/v1  
- GitHub: https://github.com/alexshen2468-ai  

---

## Quick Start 🚀

Run the model in seconds:

```bash
pip install numpy scipy matplotlib
python exosome_pk_model.py
This will generate:
	•	time-course simulation of exosome distribution
	•	retention vs clearance dynamics
	•	parameter sensitivity analysis

What This Model Does

This model provides a simplified computational framework to study exosome retention behavior in biological systems.

It explains the commonly observed phenomenon:
	•	rapid clearance from circulation
	•	prolonged retention in organs (e.g., liver, spleen)

The model introduces a key parameter:
R = k_bind / k_rel
which represents the balance between:
	•	tissue binding / uptake
	•	release / clearance

Model Structure

The pharmacokinetic (PK) model includes:
	•	systemic circulation compartment
	•	organ distribution
	•	reversible tissue binding

Key processes:
	•	systemic clearance
	•	blood–organ exchange
	•	reversible binding

Example Use Cases
	•	studying exosome biodistribution
	•	testing retention hypotheses
	•	early-stage drug delivery modeling
	•	conceptual PK/PD exploration

Computational Approach
	•	ordinary differential equation (ODE) modeling
	•	nonlinear parameter estimation
	•	sensitivity analysis
	•	simulation-based exploration

Implemented entirely in Python.

Reproducibility

This repository contains the full computational implementation used in the associated preprint.

All results can be reproduced using the provided script and dependencies.

Preprint

Shen, Z. (2026)
Simulation-Based Exosome Nanocarrier Model for Targeted Opioid Analgesia
Research Square

https://doi.org/10.21203/rs.3.rs-7939584/v1

Citation

If you use this repository, please cite:

Shen, Z. (2026)
Simulation-Based Exosome Nanocarrier Model for Targeted Opioid Analgesia
Research Square

https://doi.org/10.21203/rs.3.rs-7939584/v1

Repository Structure
exosome_pk_model.py                  # Main simulation script
README.md                            # Documentation
LICENSE                              # MIT license
CITATION.cff                         # Citation metadata
Retention Ratio PLOSONE v4 revised.docx   # Manuscript draft

Author

Zhuofan Shen
Undergraduate researcher in biotechnology

Research focus:
	•	computational modeling of exosome-mediated drug delivery
	•	retention dynamics and biodistribution
	•	pharmacokinetic (PK/PD) modeling

Affiliation: University of Debrecen
	•	ORCID: https://orcid.org/0009-0005-4304-8391
	•	GitHub: https://github.com/alexshen2468-ai

Keywords

exosome · pharmacokinetics · drug delivery · nanomedicine · computational biology · PK modeling · biodistribution

License

MIT License