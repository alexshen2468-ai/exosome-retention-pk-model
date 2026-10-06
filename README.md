# Exosome Retention PK Model

Mechanistic multi‑scale pharmacokinetic model for exosome biodistribution and hepatic retention. Introduces the dimensionless **retention ratio $R = k_{\text{bind}} / k_{\text{rel}}$** as a compact, biologically interpretable design parameter for exosome‑based drug delivery.
Model developed following ISEV extracellular vesicle research standards.

> **Core thesis: Write the organ curve first. Then the vesicle.**
>
> Antibodies have $K_D$, epidemics have $R_0$. Exosome liver tropism can be quantified with retention ratio $R$. This framework treats $R = k_{\text{bind}}/k_{\text{rel}}$ as a design spec you put on the vial — and supports inverse design to solve for feasible $(R, k_{\text{rel}})$ pairs from target organ PK constraints.
>
> Interactive web demo: **R‑SPEC (Inverse design for exosome retention)**
>
> X(Twitter) thread: https://x.com/kangzhous/status/2097641515044499484?s=46

-----

## Quick Start

```bash
pip install numpy scipy matplotlib SALib
python "exosome retention model (New).py"
```

Reproduces all 12 main figures and 4 supplementary figures from the manuscript.

-----

## Background

Intravenously administered exosomes show a consistent paradox: blood concentrations fall rapidly (initial half‑life ≈ 1–2 h), yet hepatic levels remain substantially elevated for 24–168 h. Classical one‑ or two‑compartment PK models cannot explain this — they predict organ concentrations that track blood.

This model resolves the paradox by introducing an explicit **reversible hepatic reservoir** governed by $R = k_{\text{bind}} / k_{\text{rel}}$.

The R‑SPEC web interface builds directly on this ODE reservoir model, enabling **inverse design**: define your desired hepatic plateau, liver Cmax ceiling, and blood half‑life floor; the solver searches the $(R, k_{\text{rel}})$ feasible plane to output engineering specs for vesicle fabrication.

> Key design insight from Sobol sensitivity: $k_{\text{rel}}$ (off‑rate) dominates hepatic AUC ($S_1=0.182$ vs $k_{\text{bind}}$ $S_1=0.09$). Slowing release outperforms raising binding avidity. Raising $R$ only via higher $k_{\text{bind}}$ still allows washout if $k_{\text{rel}}$ remains fast. $R$ and $k_{\text{rel}}$ are a paired specification, not a single slogan.

-----

## Model Architecture

Five coupled ODEs:

```
dB/dt  = −(k_clear + k_to + k_bp)·B + k_pb·P                        (1)
dP/dt  =   k_bp·B − k_pb·P                                           (2)
dLf/dt =   k_to·B + k_rel·Lb − k_bind·Lf − Vmax·Lf/(Km + Lf)       (3)
dLb/dt =   k_bind·Lf − k_rel·Lb                                      (4)
R      =   k_bind / k_rel                                             (5)
```

|Variable|Description                  |
|--------|-----------------------------|
|B       |Blood exosome (%ID)          |
|P       |Peripheral tissue (%ID)      |
|Lf      |Free hepatic exosome (%ID)   |
|Lb      |Bound hepatic reservoir (%ID)|

At pseudo‑steady state, bound fraction = $R / (1 + R)$. Fitted $R ≈ 5.3$ → ~84 % of hepatic exosomes in reservoir.

A pharmacodynamic sub‑model (Michaelis‑Menten uptake → drug release → hyperbolic $E_\text{max}$) is included as an illustrative projection; PD parameters are fixed at literature‑informed values and not fitted to data.

An optional irreversible loss term $k_{\text{loss}}$ models Kupffer cell degradation / biliary clearance, recovering the late‑phase ⁹⁹Zr tracer falloff observed in Choi et al. 2022.

-----

## 🖥️ R‑SPEC Interactive Inverse Design Tool
R‑SPEC is the web frontend built on top of this reservoir PK model for **exosome retention inverse design**.

### Capabilities
1. **Simulation mode**: Tune $R$, $k_{\text{bind}}$, $k_{\text{rel}}$, $k_{\text{loss}}$ via sliders and visualize real‑time trajectories: blood PK, liver observed‑scale curve, reservoir split ($L_f$, $L_b$, peripheral pool), PD projection, and residual table vs Choi ⁹⁹Zr dataset.
2. **Design mode**: Define your therapeutic intent as hard constraints:
   - Hepatic plateau duration
   - Liver Cmax ceiling (%ID)
   - Minimum blood half‑life
3. Feasible plane search: solver scans $(R, k_{\text{rel}})$ parameter space, scores candidates, and outputs lab‑ready lot specifications.
4. Sensitivity & validation modules, export trajectory data.

> Screenshot preview of R‑SPEC UI (add image to repo):
> ![R‑SPEC UI](Fig_RSPEC_demo.png)

### Core design rule
> You specify the organ‑level PK target first, then derive the required vesicle biophysical parameters ($R$, $k_{\text{rel}}$). R‑SPEC does not invent PEG‑density calibration; it outputs $R$ and $k_{\text{rel}}$ specs that you measure on your final vesicle lot.

-----

## Key Results

### Model comparison (AICc)

|Model |Structure                    |ΔAICc|Support          |
|------|-----------------------------|-----|-----------------|
|M1    |1‑compartment, no reservoir  |19.6 |None             |
|M2    |2‑compartment, no reservoir  |8.4  |Considerably less|
|**M3**|**2‑compartment + reservoir**|**0**|**Best         |

$\Delta\text{AICc} > 10$ = decisive evidence against (Burnham–Anderson criteria).

### Fitted parameters

|Parameter|Value  |Units|Bootstrap 95% CI|
|---------|-------|-----|----------------|
|k_clear  |0.350  |h⁻¹  |[0.301, 0.408]  |
|k_to     |0.150  |h⁻¹  |[0.118, 0.192]  |
|k_bind   |0.020  |h⁻¹  |[0.014, 0.029]  |
|k_rel    |0.004  |h⁻¹  |[0.0028, 0.0061]|
|**R**    |**5.3**|—    |—               |

Bootstrap $n = 500$ resamples. RMSE = 3.06 %ID across both compartments.

### Sobol sensitivity (hepatic AUC)

$k_{\text{rel}}$ is the dominant parameter ($S_1 = 0.182$). Reducing $k_{\text{rel}}$ (slow off‑rate engineering) is more effective at prolonging hepatic retention than increasing $k_{\text{bind}}$ (binding avidity).

-----

## Data

### Primary dataset — hard‑coded in script

Choi et al. (2022), ⁹⁹Zr‑labelled GMP‑grade exosomes, ICR mice, $n = 4$/timepoint:

|time (h)|blood (%ID)|±SD|liver (%ID)|±SD|
|--------|-----------|---|-----------|---|
|0.25    |52.1       |6.8|8.3        |1.4|
|1       |28.4       |4.2|22.6       |3.1|
|2       |18.6       |3.8|28.1       |3.7|
|6       |8.3        |2.1|38.9       |4.2|
|24      |3.1        |0.9|45.2       |5.8|
|48      |1.4        |0.4|41.7       |6.1|
|72      |0.8        |0.3|36.4       |5.4|
|120     |0.4        |0.2|28.8       |4.9|
|168     |0.2        |0.1|21.3       |3.6|

### Cross‑validation datasets — loaded from CSV

Normalised to $t = 1\ \text{h}$ (fluorescence; no absolute %ID calibration).

**`mirzaaghasi_normalized.csv`** — DiR‑labelled HEK293T, ICR mice, 1–8 h

|time_h|blood_norm|liver_norm|
|------|----------|----------|
|1.0   |1.000     |0.190     |
|2.0   |0.682     |0.420     |
|3.0   |0.441     |0.601     |
|4.0   |0.312     |0.710     |
|8.0   |0.133     |0.882     |

**`wiklander_normalized.csv`** — DiI‑labelled C2C12, 1–24 h

|time_h|blood_norm|liver_norm|
|------|----------|----------|
|1.0   |1.000     |0.280     |
|3.0   |0.300     |0.640     |
|6.0   |0.095     |0.800     |
|24.0  |0.085     |0.960     |

No parameter refitting was performed for cross‑validation.

-----

## Repository Structure

```
exosome retention model (New).py   # Main simulation script
mirzaaghasi_normalized.csv         # Cross‑validation dataset 1
wiklander_normalized.csv           # Cross‑validation dataset 2
Fig_retention_curves.png           # Example output
Fig_RSPEC_demo.png                 # R‑SPEC web UI screenshot (add your image here)
requirements.txt                   # numpy scipy matplotlib SALib
LICENSE                            # MIT
README.md
.gitignore
```

-----

## Figures Generated

|Figure|Description                                     |
|------|------------------------------------------------|
|Fig 1 |Blood concentration — M3 vs Choi et al.         |
|Fig 2 |Hepatic concentration — M3 vs ⁹⁹Zr signal       |
|Fig 3 |Intracellular uptake rate (illustrative)        |
|Fig 4 |Local drug concentration (illustrative)         |
|Fig 5 |Normalised therapeutic effect (illustrative)    |
|Fig 6 |M1 / M2 / M3 comparison on blood data           |
|Fig 7 |AICc bar chart                                  |
|Fig 8 |Bootstrap distribution of k_rel (n = 500)       |
|Fig 9 |Sobol sensitivity indices                       |
|Fig 10|Classical vs reservoir model schematic          |
|Fig 11|Plateau duration vs retention ratio R           |
|Fig 12|Label persistence sensitivity analysis (4‑panel)|
|Fig S1|Hepatic profiles at R = 1, 3, 6, 10             |
|Fig S2|2D landscape: plateau duration vs (R, k_rel)    |
|Fig S3|Cross‑validation — Mirzaaghasi + Wiklander      |
|Fig S4|Progressive k_rel reduction scenarios           |

-----

## Limitations

- Absolute hepatic concentration is underestimated (~3–4×) due to ⁹⁹Zr radiolabel persistence; formally characterised in Fig 12
- Cross‑validation on normalised fluorescence data only — not absolute %ID
- Liver treated as a single well‑mixed compartment; Kupffer cell / LSEC heterogeneity not modelled
- PD parameters ($V_\text{max}$, $K_m$, $\text{EC}_{50}$) fixed at literature values; Figs 3–5 are illustrative only
- Protein corona dynamics not modelled ($k_{\text{bind}}$, $k_{\text{rel}}$ assumed time‑invariant)
- R‑SPEC inverse design is a computational prototype; outputs are **in silico specs**, not validated experimental predictions.

-----

## Preprint

Shen, Z. (2026). *Mechanistic modeling of exosome pharmacokinetics reveals a retention ratio governing sustained organ accumulation.* Research Square.
https://doi.org/10.21203/rs.3.rs-7939584/v1

> Related work: [A Simplified Computational Framework for Modeling Exosome‑Mediated Local Drug Delivery](https://doi.org/10.21203/rs.3.rs-8968566/v1)

-----

## Citation

```bibtex
@misc{shen2026exosome,
  author    = {Shen, Zhuofan},
  title     = {Mechanistic modeling of exosome pharmacokinetics reveals
               a retention ratio governing sustained organ accumulation},
  year      = {2026},
  publisher = {Research Square},
  doi       = {10.21203/rs.3.rs-7939584/v1}
}
```

-----

## Author

**Zhuofan Shen (Kangzhou Shen)**
> Junior Member, International Society for Extracellular Vesicles (ISEV, 2026‑2028)
> BSc Biotechnology, University of Debrecen | BSc Bioengineering, Universitat Internacional de Catalunya (UIC Barcelona)

- ORCID: [0009‑0005‑4304‑8391](https://orcid.org/0009-0005-4304-8391)
- Email: alexshen2468@gmail.com
- GitHub: [@alexshen2468‑ai](https://github.com/alexshen2468-ai)
- X(Twitter): [@kangzhous](https://twitter.com/kangzhous)

-----

## Related Resources
- Zenodo Dataset: [10.5281/zenodo.19193005](https://zenodo.org/record/19193005)
- X Thread for R‑SPEC demo: https://x.com/kangzhous/status/2097641515044499484?s=46
- Profile Repository: [github.com/alexshen2468‑ai](https://github.com/alexshen2468-ai)

## License

MIT © 2026 alexshen2468‑ai

-----

**Keywords:** exosome · pharmacokinetics · retention ratio · drug delivery · ODE modeling · biodistribution · reservoir model · nanomedicine · computational pharmacology · inverse design · R‑SPEC
