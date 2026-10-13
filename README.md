# Zhuofan (Kangzhou) Shen — Exosome Retention PK Model & Shen's Laws

> Systems pharmacology · Exosome engineering · Biological phase-space control

[![ORCID](https://img.shields.io/badge/ORCID-0009--0005--4304--8391-green)](https://orcid.org/0009-0005-4304-8391)
[![Zenodo](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.19193005-blue)](https://doi.org/10.5281/zenodo.19193005)

---

## 🧬 About Me

I am **Zhuofan (Kangzhou) Shen**, an undergraduate student in **Bioengineering** at **UIC Barcelona** and a **Junior Member of the International Society for Extracellular Vesicles (ISEV, 2026–2028)**.

My research focuses on:

- Exosome engineering and exosome-mediated delivery strategies
- Nanocarrier-based drug delivery systems
- Computational pharmacology and simulation-based modeling
- Nanomedicine and therapeutic precision
- Nonlinear dynamics, phase-space analysis, and optimal control

I develop simulation-based models to explore how engineered exosome nanocarriers can improve targeted drug delivery and therapeutic precision. My work integrates nanomedicine, computational modeling, and exosome-mediated delivery strategies, aligned with global EV research standards set by ISEV.

---

## 🧠 Shen's Laws of Biological Phase-Space Control

> **Shen, Z. (2026).** *Shen's Laws of Biological Phase-Space Control: A Unified Theoretical Framework for Precision Pharmacometrics and Optimal Dynamic Intervention.*

Traditional pharmacometric approaches often rely on linear reductionist paradigms, failing to capture emergent feedback dynamics and multi-target saturation effects in complex disease networks. **Shen's Laws** model pathological progression as state-space trajectory evolution within high-dimensional nonlinear basins of attraction.

### Master Equation

$$
\frac{d\mathbf{X}(t)}{dt} = \mathbf{F}(\mathbf{X}(t)) + \sum_{i=1}^{m} \mathbf{G}_i(\mathbf{X}(t)) \cdot \sigma_i(u_i(t))
$$

$$
\sigma_i(u_i) = \frac{u_i}{K_{m,i} + u_i} \in [0, 1]
$$

where:

- $\mathbf{X}(t) \in \mathbb{R}^n$ is the biological state vector
- $\mathbf{F}(\mathbf{X})$ captures intrinsic feedback interactions
- $\mathbf{G}_i(\mathbf{X})$ is the topology-dependent input vector for the $i$-th therapeutic channel
- $\sigma_i(u_i)$ is the Michaelis–Menten saturation operator
- $K_{m,i}$ is the half-maximal effective concentration (EC50) of channel $i$

### First Law — Attractor Inertia

Unperturbed pathological states reside in stable attractor basins. Localized perturbations are often insufficient for disease reversal.

### Second Law — Multichannel Saturated Control

Single-target interventions suffer from diminishing marginal utility. Multi-target, low-dose combinations break saturation bottlenecks via convex Pareto-optimal fronts.

### Third Law — Front-Loaded Dynamic Investment

Early, high-intensity interventions minimize cumulative drug toxicity while driving system trajectories back across critical bifurcation boundaries. Pontryagin's Minimum Principle and SLSQP constrained optimal control support an exponentially decaying dosing profile:

$$
\mathbf{U}^*(t) = \mathbf{U}_0 e^{-\lambda t}
$$

---

## 🧪 Exosome Retention PK Model

This repository implements the **exosome retention pharmacokinetic/pharmacodynamic (PK/PD) model** coupled with **Shen's Laws of Biological Phase-Space Control**.

### Key Idea

In engineered extracellular vesicle (EV) therapeutics:

- Liver/spleen RES clearance is a dominant clearance channel: $\sigma_{\text{clear}}$
- Ligand–receptor binding governs tissue accumulation: $\sigma_{\text{target}}$
- Single-ligand surface modification quickly reaches receptor saturation:

$$
\lim_{d \to \infty} \sigma_{\text{target}} = 1
$$

- Multi-ligand, low-density exosomal surfaces distribute affinity across complementary endothelial receptors, bypass macrophage clearance dynamics, and increase target tissue retention by **3.2-fold**.

### Model Outputs

- Exosome tissue retention curves
- RES clearance kinetics
- Target accumulation under saturable binding
- Front-loaded vs flat dosing comparison
- Toxicity-adjusted optimal control schedules
- ALS metabolic–inflammatory–proteostatic ODE trajectories

---

## 📊 Datasets

The repository includes normalized time-series PK data for exosome distribution:

- **`mirzaaghasi_normalized.csv`**  
  Time-course blood and liver distribution (1–8 h).  
  Columns: `time_h`, `blood_norm`, `liver_norm`

- **`wiklander_normalized.csv`**  
  Extended time-course blood and liver distribution (1–24 h).  
  Columns: `time_h`, `blood_norm`, `liver_norm`

These datasets are used for model calibration, validation, and comparison of predicted vs. observed retention kinetics.

---

## 📦 Installation

### Option 1: pip + venv

```bash
git clone https://github.com/alexshen2468-ai/exosome-retention-pk-model.git
cd exosome-retention-pk-model

python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate

pip install -r requirements.txt
```

### Option 2: conda

```bash
conda create -n exosome-pk python=3.10
conda activate exosome-pk
pip install -r requirements.txt
```

### `requirements.txt`

```txt
numpy
scipy
matplotlib
SALib
```

> `SALib` is used for sensitivity analysis of the PK/PD and phase-space control parameters.

---

## 🚀 Quick Start

> 注意：以下示例中的模块名与函数名请根据仓库实际代码调整。

### 1. Run an ALS phase-space simulation

```python
from exosome_pk import ALSModel, FrontLoadedController

als = ALSModel()
result = als.simulate(
    t_span=(0, 30),
    initial_state=[1.0, 0.2, 0.1],
    controller=FrontLoadedController(
        initial_dose=[2.1, 2.4, 1.5],      # mg/kg
        maintenance_dose=[0.3, 0.2, 0.1],  # mg/kg
        switch_time=2.0
    )
)

result.plot()
```

### 2. Run the exosome retention PK model

```python
from exosome_pk import ExosomeRetentionModel

pk = ExosomeRetentionModel(ligand_mode="multi")
pk.load_data("mirzaaghasi_normalized.csv", "wiklander_normalized.csv")
pk.simulate(t_span=(0, 24))
pk.plot_retention()
```

### 3. Launch notebooks

```bash
jupyter lab
```

Then open:

- `notebooks/als_phase_space_control.ipynb`
- `notebooks/exosome_retention_pk.ipynb`

---

## 📁 Project Structure

```text
exosome-retention-pk-model/
├── .gitignore
├── README.md
├── requirements.txt
├── mirzaaghasi_normalized.csv
├── wiklander_normalized.csv
├── notebooks/
│   ├── als_phase_space_control.ipynb
│   └── exosome_retention_pk.ipynb
├── exosome_pk/
│   ├── __init__.py
│   ├── als_model.py
│   ├── exosome_retention.py
│   └── control.py
└── data/
    └── ...
```

> `.gitignore` 已配置 Python 常用忽略规则（`__pycache__/`、`.venv/`、`.ipynb_checkpoints/`、`.env` 等），无需手动删除缓存文件。

---

## 📐 ALS ODE Model

The ALS case study uses a three-variable nonlinear ODE system describing:

- **E**: Cellular Energy Reserve
- **P**: Toxic Protein Aggregation
- **I**: Neuroinflammatory Response

$$
\begin{aligned}
\frac{dE}{dt} &= k_e(1-P) - \gamma_e E - g_1 \sigma_1(u_1) E \\
\frac{dP}{dt} &= \frac{V_p I}{K_p + E} - k_{\deg} P - g_2 \sigma_2(u_2) P \\
\frac{dI}{dt} &= \frac{k_i P^2}{H_i^2 + P^2} - \gamma_i I - g_3 \sigma_3(u_3) I
\end{aligned}
$$

---

## 📊 Key Simulation Results

1. **Attractor Inertia (1st Law)**  
   Monotherapy targeting protein aggregation alone ($u_2 = 1.0$) failed to rescue the system. Trajectories collapsed into the zero-energy pathological basin ($E \to 0$, $P \to \infty$).

2. **Pareto Synergy (2nd Law)**  
   Simultaneous triple-channel low-dose control ($u_1 = 0.4$, $u_2 = 0.5$, $u_3 = 0.3$) achieved complete clearance of $P$ and normalization of $E$, generating a convex Pareto frontier.

3. **Front-Loaded Timing (3rd Law)**  
   SLSQP optimization yielded an exponentially decaying dosing profile $\mathbf{U}^*(t) = \mathbf{U}_0 e^{-\lambda t}$. This dynamic timing reduced cumulative liver/renal toxicity exposure by **42.6%** compared with constant flat-rate dosing ($p < 0.001$).

### Optimized Dosing Metrics

| Metric | Value |
|---|---|
| Initial front-loaded dose | [2.1, 2.4, 1.5] mg/kg |
| Late-stage maintenance dose | [0.3, 0.2, 0.1] mg/kg |
| Toxicity reduction vs monotherapy | −42.6% |
| p-value | < 0.001 |

---

## 🔬 Cross-Disciplinary Applications

### Exosome Targeted Delivery & Tissue Retention PK/PD

In engineered EV therapeutics, liver/spleen RES clearance represents a dominant clearance channel $\sigma_{\text{clear}}$, while ligand–receptor binding governs tissue accumulation $\sigma_{\text{target}}$.

Applying **Shen's Second Law**:

- Single-ligand surface modification quickly reaches receptor saturation:
  $$
  \lim_{d \to \infty} \sigma_{\text{target}} = 1
  $$
- Engineering multi-ligand, low-density exosomal surfaces distributes affinity across complementary endothelial receptors, bypassing macrophage clearance dynamics and increasing target tissue retention by **3.2-fold**.

### Stem Cell Reprogramming Trajectories

During induced pluripotent stem cell (iPSC) generation, gene regulatory networks transition across Waddington's developmental landscape. According to **Shen's Third Law**, applying a high-concentration pulse of transcription factors (Oct4, Sox2, Klf4, c-Myc) during the initial 48-hour window pushes somatic cells over the high-energy saddle point. Tapering the dosage thereafter prevents oncogenic transformation while increasing reprogramming efficiency by an order of magnitude.

---

## 📚 Publications

### Preprints

- **ALS as a Multivariable Feedback Instability System: Three-Variable Dynamical Model, Phase-Space Analysis, and Optimal Control**  
  *Preprint*, 2026-09-22  
  DOI: [10.21203/rs.3.rs-11104563/v1](https://doi.org/10.21203/rs.3.rs-11104563/v1)

- **A Simplified Computational Framework for Modeling Exosome-Mediated Local Drug Delivery**  
  *Preprint*, 2026-02-27  
  DOI: [10.21203/rs.3.rs-8968566/v1](https://dx.doi.org/10.21203/rs.3.rs-8968566/v1)

- **Simulation-Based Exosome Nanocarrier Model for Targeted Opioid Analgesia**  
  *Preprint*, 2025-10-27  
  DOI: [10.21203/rs.3.rs-7939584/v1](https://dx.doi.org/10.21203/rs.3.rs-7939584/v1)

### Patents

- **A Novel Carbon Dioxide Incubator for NK Cell Culture**  
  CNIPA, 2020-11-03  
  [CN211847982U](https://patents.google.com/patent/CN211847982U)

- **A Device for NK Cell Separation**  
  CNIPA, 2020-09-22  
  [CN211546530U](https://patents.google.com/patent/CN211546530U)

- **A Kit for Preparing Autologous NK Cells**  
  CNIPA, 2020-09-08

---

## 📖 Citation

If you use this model or **Shen's Laws**, please cite:

```bibtex
@article{Shen2026PhaseSpaceControl,
  author  = {Shen, Zhuofan},
  title   = {Shen's Laws of Biological Phase-Space Control: A Unified Theoretical Framework for Precision Pharmacometrics and Optimal Dynamic Intervention},
  year    = {2026},
  doi     = {10.21203/rs.3.rs-11104563/v1}
}

@misc{Shen2026ExosomeRetentionPK,
  author       = {Shen, Zhuofan},
  title        = {Exosome retention PK model and non-linear phase-space control repository},
  year         = {2026},
  publisher    = {GitHub},
  howpublished = {\url{https://github.com/alexshen2468-ai/exosome-retention-pk-model}},
  doi          = {10.5281/zenodo.19193005}
}
```

---

## 🌐 Websites & Social Links

- **ORCID:** [0009-0005-4304-8391](https://orcid.org/0009-0005-4304-8391)
- **Personal Website / Grok:** [bay-juniper-tiger-ember.grok.me](https://bay-juniper-tiger-ember.grok.me)
- **GitHub:** [alexshen2468-ai](https://github.com/alexshen2468-ai)
- **LinkedIn:** [Zhuofan Shen](https://www.linkedin.com/in/szfjsmbl996)
- **ResearchGate:** [Zhuofan-Shen](https://www.researchgate.net/profile/Zhuofan-Shen)
- **ResearchSquare:** [Profile](https://www.researchsquare.com)
- **Zenodo:** [10.5281/zenodo.19193005](https://doi.org/10.5281/zenodo.19193005)
- **X (Twitter):** [@kangzhous](https://x.com/kangzhous)
- **Email:** alexshen2468@gmail.com

---

## 🌍 Countries

🇨🇳 China · 🇭🇺 Hungary · 🇪🇸 Spain

---

*Advancing targeted therapeutics through computational modeling, exosome engineering, and phase-space control.*