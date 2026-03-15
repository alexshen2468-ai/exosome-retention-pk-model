Mechanistic modeling of exosome pharmacokinetics reveals a retention ratio
governing sustained organ accumulation.
Shen Zhuofan, University of Debrecen Biotechnology Program
This script implements:
 - Three nested ODE model structures (M1, M2, M3)
 - Parameter fitting via Levenberg-Marquardt with 10 random initialisations
 - AICc-based model comparison
 - Bootstrap parameter uncertainty analysis (n=500)
 - Sobol global sensitivity analysis
 - Independent cross-validation (Mirzaaghasi et al. 2021; Wiklander et al. 2015)
 - All manuscript figures (Fig 1–15)
Requirements:
 pip install numpy scipy matplotlib SALib
Usage:
 python exosome_pk_model_v2.py
Changes from v1:
 - Cross-validation extended to include Wiklander et al. (2015) dataset
 - MAE reporting now disaggregated by compartment (blood / liver) and dataset
 - Fig 14 updated to show both cross-validation datasets (2×2 panel)
 - Fig 15 added: k_rel reduction simulation (25 / 50 / 75%)
 - Minor comment corrections (99Zr isotope notation)
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
# SALib is optional; if unavailable, a pure-numpy Monte-Carlo Sobol
# approximation is used instead.
try:
 from SALib.sample import saltelli as _saltelli
 from SALib.analyze import sobol as _sobol
 _SALIB_AVAILABLE = True
except ImportError:
 _SALIB_AVAILABLE = False
# ---------------------------------------------------------------------------
# 1. EXPERIMENTAL DATA
# ---------------------------------------------------------------------------
# Primary dataset: Choi et al. Pharmaceutics 2022 [16]
# 99Zr-labelled GMP-grade exosomes; %ID = percentage of injected dose;
# mean +/- SD, n=4 mice per time point.
choi_time = np.array([0.25, 1, 2, 6, 24, 48, 72, 120, 168]) # hours
choi_blood_mean = np.array([52.1, 28.4, 18.6, 8.3, 3.1, 1.4, 0.8, 0.4, 0.2])
choi_blood_sd = np.array([ 6.8, 4.2, 3.8, 2.1, 0.9, 0.4, 0.3, 0.2, 0.1])
choi_liver_mean = np.array([ 8.3, 22.6, 28.1, 38.9, 45.2, 41.7, 36.4, 28.8, 21.3])
choi_liver_sd = np.array([ 1.4, 3.1, 3.7, 4.2, 5.8, 6.1, 5.4, 4.9, 3.6])
# Cross-validation dataset 1: Mirzaaghasi et al. Pharmaceutics 2021 [17]
# DiR-labelled HEK293T exosomes; values normalised to t=1 h
# (no absolute %ID calibration available).
mirzaaghasi_time = np.array([1, 2, 3, 4, 8]) # hours
mirzaaghasi_blood_norm = np.array([1.000, 0.682, 0.441, 0.312, 0.133])
mirzaaghasi_liver_norm = np.array([0.190, 0.420, 0.601, 0.710, 0.882])
# Cross-validation dataset 2: Wiklander et al. J Extracell Vesicles 2015 [11]
# DiI-labelled C2C12 exosomes; values normalised to t=1 h
# (no absolute %ID calibration available).
wiklander_time = np.array([1, 3, 6, 24]) # hours
wiklander_blood_norm = np.array([1.000, 0.300, 0.095, 0.085])
wiklander_liver_norm = np.array([0.280, 0.640, 0.800, 0.960])
# ---------------------------------------------------------------------------
# 2. ODE MODEL DEFINITIONS
# ---------------------------------------------------------------------------
def ode_M1(t, y, params):
 """
 M1: One-compartment model (blood only, simple clearance).
 State: [B]
 """
 B = y[0]
 k_clear = params[0]
 dB = -k_clear * B
 return [dB]
def ode_M2(t, y, params):
 """
 M2: Two-compartment model (blood + peripheral), organ tracks blood.
 States: [B, P, L_f]
 L_f = organ free, driven by blood, no reservoir.
 """
 B, P, L_f = y
 k_clear, k_to, k_bp, k_pb = params
 dB = -(k_clear + k_to + k_bp) * B + k_pb * P
 dP = k_bp * B - k_pb * P
 dL_f = k_to * B - k_clear * L_f # organ clears at same rate as blood
 return [dB, dP, dL_f]
def ode_M3(t, y, params):
 """
 M3: Reservoir model (proposed).
 States: [B, P, L_f, L_b, C_cell, C_drug]
 Equations (manuscript numbering):
 dB/dt = -(k_clear + k_to + k_bp)*B + k_pb*P (1)
 dP/dt = k_bp*B - k_pb*P (2)
 dL_f/dt = k_to*B + k_rel*L_b - k_bind*L_f - V_max*L_f/(K_m+L_f) (3)
 dL_b/dt = k_bind*L_f - k_rel*L_b (4)
 dC_cell/dt= V_max*L_f/(K_m+L_f) - k_deg*C_cell (5)
 dC_drug/dt= q*C_cell - k_clear_d*C_drug (6)
 Effect = C_drug/(C_drug + EC50) (7, algebraic)
 """
 B, P, L_f, L_b, C_cell, C_drug = y
 (k_clear, k_to, k_bp, k_pb,
 k_bind, k_rel,
 V_max, K_m, k_deg,
 q, k_clear_d, EC50) = params
 MM = V_max * L_f / (K_m + L_f)
 dB = -(k_clear + k_to + k_bp) * B + k_pb * P
 dP = k_bp * B - k_pb * P
 dL_f = k_to * B + k_rel * L_b - k_bind * L_f - MM
 dL_b = k_bind * L_f - k_rel * L_b
 dC_cell = MM - k_deg * C_cell
 dC_drug = q * C_cell - k_clear_d * C_drug
 return [dB, dP, dL_f, dL_b, dC_cell, dC_drug]
def simulate_M3(params, t_eval, B0=52.1):
 """Integrate M3 ODE system and return solution object."""
 y0 = [B0, 0.0, 0.0, 0.0, 0.0, 0.0]
 sol = solve_ivp(
 ode_M3, [0, t_eval[-1]], y0,
 args=(params,),
 t_eval=t_eval,
 method='RK45',
 rtol=1e-6, atol=1e-9,
 dense_output=False
 )
 return sol
def simulate_M2(params, t_eval, B0=52.1):
 y0 = [B0, 0.0, 0.0]
 sol = solve_ivp(
 ode_M2, [0, t_eval[-1]], y0,
 args=(params,),
 t_eval=t_eval,
 method='RK45',
 rtol=1e-6, atol=1e-9
 )
 return sol
def simulate_M1(params, t_eval, B0=52.1):
 y0 = [B0]
 sol = solve_ivp(
 ode_M1, [0, t_eval[-1]], y0,
 args=(params,),
 t_eval=t_eval,
 method='RK45',
 rtol=1e-6, atol=1e-9
 )
 return sol
# ---------------------------------------------------------------------------
# 3. PARAMETER ESTIMATION
# ---------------------------------------------------------------------------
# Nominal (literature-informed) parameter values for M3
# [k_clear, k_to, k_bp, k_pb, k_bind, k_rel,
# V_max, K_m, k_deg, q, k_clear_d, EC50]
P_NOMINAL = np.array([
 0.35, # k_clear (h-1) – fitted
 0.15, # k_to (h-1) – fitted
 0.08, # k_bp (h-1) – fitted
 0.03, # k_pb (h-1) – fitted
 0.020, # k_bind (h-1) – fitted
 0.004, # k_rel (h-1) – fitted
 0.5, # V_max (%ID/h) – literature [35]
 10.0, # K_m (%ID) – literature [35]
 0.01, # k_deg (h-1) – literature [35]
 0.05, # q (h-1) – assumed
 0.02, # k_clear_d (h-1) – assumed
 5.0, # EC50 (a.u.) – assumed
])
# For fitting we optimise only the PK parameters (indices 0–5)
# and fix PD parameters at literature/assumed values.
FIT_IDX = [0, 1, 2, 3, 4, 5] # k_clear, k_to, k_bp, k_pb, k_bind, k_rel
FIT_LOWER = np.array([0.01, 0.01, 0.001, 0.001, 0.001, 0.0001])
FIT_UPPER = np.array([2.0, 1.0, 0.5, 0.3, 0.5, 0.1 ])
def residuals_M3(x_fit, p_fixed, t_obs, blood_obs, liver_obs,
 blood_sd, liver_sd):
 """Weighted residuals for M3 blood + liver fitting."""
 params = p_fixed.copy()
 params[FIT_IDX] = x_fit
 try:
 sol = simulate_M3(params, t_obs)
 if not sol.success or np.any(np.isnan(sol.y)):
 return np.ones(2 * len(t_obs)) * 1e6
 B_pred = sol.y[0]
 L_pred = sol.y[2] + sol.y[3] # L_f + L_b = total organ
 res_B = (B_pred - blood_obs) / (blood_sd + 1e-9)
 res_L = (L_pred - liver_obs) / (liver_sd + 1e-9)
 return np.concatenate([res_B, res_L])
 except Exception:
 return np.ones(2 * len(t_obs)) * 1e6
def fit_M3(t_obs, blood_obs, liver_obs, blood_sd, liver_sd,
 n_restarts=10, seed=42):
 """Fit M3 with multiple random initialisations; return best result."""
 rng = np.random.default_rng(seed)
 best_cost = np.inf
 best_result = None
 for _ in range(n_restarts):
 x0 = rng.uniform(FIT_LOWER, FIT_UPPER)
 try:
 res = least_squares(
 residuals_M3, x0,
 args=(P_NOMINAL.copy(), t_obs,
 blood_obs, liver_obs, blood_sd, liver_sd),
 bounds=(FIT_LOWER, FIT_UPPER),
 method='trf',
 max_nfev=5000
 )
 if res.cost < best_cost:
 best_cost = res.cost
 best_result = res
 except Exception:
 continue
 return best_result
# M2 fitting helpers
P_M2_NOMINAL = np.array([0.35, 0.15, 0.08, 0.03])
M2_LOWER = np.array([0.01, 0.01, 0.001, 0.001])
M2_UPPER = np.array([2.0, 1.0, 0.5, 0.3 ])
def residuals_M2(x, t_obs, blood_obs, liver_obs, blood_sd, liver_sd):
 try:
 sol = simulate_M2(x, t_obs)
 if not sol.success:
 return np.ones(2 * len(t_obs)) * 1e6
 B_pred = sol.y[0]
 L_pred = sol.y[2]
 res_B = (B_pred - blood_obs) / (blood_sd + 1e-9)
 res_L = (L_pred - liver_obs) / (liver_sd + 1e-9)
 return np.concatenate([res_B, res_L])
 except Exception:
 return np.ones(2 * len(t_obs)) * 1e6
def fit_M2(t_obs, blood_obs, liver_obs, blood_sd, liver_sd,
 n_restarts=10, seed=42):
 rng = np.random.default_rng(seed)
 best_cost = np.inf
 best_result = None
 for _ in range(n_restarts):
 x0 = rng.uniform(M2_LOWER, M2_UPPER)
 try:
 res = least_squares(
 residuals_M2, x0,
 args=(t_obs, blood_obs, liver_obs, blood_sd, liver_sd),
 bounds=(M2_LOWER, M2_UPPER),
 method='trf', max_nfev=5000
 )
 if res.cost < best_cost:
 best_cost = res.cost
 best_result = res
 except Exception:
 continue
 return best_result
# M1 fitting
def residuals_M1(x, t_obs, blood_obs, blood_sd):
 try:
 sol = simulate_M1(x, t_obs)
 if not sol.success:
 return np.ones(len(t_obs)) * 1e6
 return (sol.y[0] - blood_obs) / (blood_sd + 1e-9)
 except Exception:
 return np.ones(len(t_obs)) * 1e6
def fit_M1(t_obs, blood_obs, blood_sd, n_restarts=10, seed=42):
 rng = np.random.default_rng(seed)
 best_cost = np.inf
 best_result = None
 for _ in range(n_restarts):
 x0 = rng.uniform([0.01], [2.0])
 try:
 res = least_squares(
 residuals_M1, x0,
 args=(t_obs, blood_obs, blood_sd),
 bounds=([0.01], [2.0]),
 method='trf', max_nfev=5000
 )
 if res.cost < best_cost:
 best_cost = res.cost
 best_result = res
 except Exception:
 continue
 return best_result
# ---------------------------------------------------------------------------
# 4. MODEL COMPARISON: AICc
# ---------------------------------------------------------------------------
def compute_aicc(residuals_vec, n_params, n_obs):
 ssr = np.sum(residuals_vec ** 2)
 n = n_obs
 k = n_params
 aic = n * np.log(ssr / n) + 2 * k
 aicc = aic + 2 * k * (k + 1) / (n - k - 1)
 return aicc
def compute_rmse(pred, obs):
 return np.sqrt(np.mean((pred - obs) ** 2))
# ---------------------------------------------------------------------------
# 5. BOOTSTRAP PARAMETER UNCERTAINTY
# ---------------------------------------------------------------------------
def bootstrap_M3(t_obs, blood_obs, liver_obs, blood_sd, liver_sd,
 fitted_params_full, n_boot=500, seed=0):
 """
 Non-parametric bootstrap: resample residuals with replacement,
 re-fit M3, collect parameter distributions.
 """
 rng = np.random.default_rng(seed)
 n = len(t_obs)
 # Baseline predictions
 sol0 = simulate_M3(fitted_params_full, t_obs)
 B0_pred = sol0.y[0]
 L0_pred = sol0.y[2] + sol0.y[3]
 res_B0 = blood_obs - B0_pred
 res_L0 = liver_obs - L0_pred
 boot_params = []
 for _ in range(n_boot):
 idx = rng.integers(0, n, size=n)
 B_boot = np.clip(B0_pred + res_B0[idx], 0, None)
 L_boot = np.clip(L0_pred + res_L0[idx], 0, None)
 res = fit_M3(t_obs, B_boot, L_boot, blood_sd, liver_sd,
 n_restarts=5, seed=int(rng.integers(0, 1_000_000)))
 if res is not None and res.success:
 boot_params.append(res.x)
 return np.array(boot_params)
# ---------------------------------------------------------------------------
# 6. SOBOL GLOBAL SENSITIVITY ANALYSIS
# ---------------------------------------------------------------------------
def _sobol_pure_numpy(fitted_params_full, n_samples=256, seed=42):
 """
 Pure-numpy Monte-Carlo variance-based sensitivity indices.
 Uses the Saltelli estimator (A/B matrix approach) for first- and
 total-order Sobol indices.
 """
 rng = np.random.default_rng(seed)
 param_names = ['k_rel', 'k_bind', 'V_max', 'EC50']
 base_vals = np.array([
 fitted_params_full[5], # k_rel
 fitted_params_full[4], # k_bind
 fitted_params_full[6], # V_max
 fitted_params_full[11], # EC50
 ])
 lo = base_vals * 0.5
 hi = base_vals * 1.5
 k = len(param_names)
 def scale(u):
 return lo + u * (hi - lo)
 def model_eval(pv):
 p = fitted_params_full.copy()
 p[5], p[4], p[6], p[11] = pv
 t_eval = np.linspace(0, 168, 300)
 try:
 sol = simulate_M3(p, t_eval)
 liver = sol.y[2] + sol.y[3]
 return np.trapz(liver, t_eval)
 except Exception:
 return np.nan
 N = n_samples
 A = scale(rng.random((N, k)))
 B = scale(rng.random((N, k)))
 fA = np.array([model_eval(A[i]) for i in range(N)])
 fB = np.array([model_eval(B[i]) for i in range(N)])
 S1 = np.zeros(k)
 ST = np.zeros(k)
 for j in range(k):
 AB = A.copy(); AB[:, j] = B[:, j]
 BA = B.copy(); BA[:, j] = A[:, j]
 fAB = np.array([model_eval(AB[i]) for i in range(N)])
 fBA = np.array([model_eval(BA[i]) for i in range(N)])
 mask = ~(np.isnan(fA) | np.isnan(fB) | np.isnan(fAB) | np.isnan(fBA))
 f0 = np.mean(fA[mask])
 Vtot = np.var(fA[mask])
 if Vtot < 1e-12:
 continue
 S1[j] = np.mean(fB[mask] * (fAB[mask] - fA[mask])) / Vtot
 ST[j] = np.mean((fA[mask] - fBA[mask]) ** 2) / (2 * Vtot)
 class _SiDict(dict):
 pass
 Si = _SiDict({'S1': np.clip(S1, 0, None), 'ST': np.clip(ST, 0, None)})
 return Si, param_names
def sobol_analysis(fitted_params_full, n_samples=1024, seed=42):
 """
 Sobol variance decomposition for hepatic AUC.
 Uses SALib if available, otherwise falls back to pure-numpy estimator.
 Parameters varied: k_rel, k_bind, V_max, EC50 (±50% of fitted values).
 """
 param_names = ['k_rel', 'k_bind', 'V_max', 'EC50']
 base_vals = [
 fitted_params_full[5],
 fitted_params_full[4],
 fitted_params_full[6],
 fitted_params_full[11],
 ]
 class _SiDict(dict):
 pass
 if _SALIB_AVAILABLE:
 problem = {
 'num_vars': 4,
 'names': param_names,
 'bounds': [[v * 0.5, v * 1.5] for v in base_vals]
 }
 param_values = _saltelli.sample(problem, n_samples,
 calc_second_order=False, seed=seed)
 t_eval = np.linspace(0, 168, 500)
 Y = np.zeros(len(param_values))
 for i, pv in enumerate(param_values):
 p = fitted_params_full.copy()
 p[5], p[4], p[6], p[11] = pv
 try:
 sol = simulate_M3(p, t_eval)
 Y[i] = np.trapz(sol.y[2] + sol.y[3], t_eval)
 except Exception:
 Y[i] = np.nan
 Si_raw = _sobol.analyze(problem, Y,
 calc_second_order=False,
 print_to_console=False)
 Si = _SiDict({'S1': Si_raw['S1'], 'ST': Si_raw['ST']})
 else:
 print(" (SALib not found – using pure-numpy Sobol estimator)")
 Si, _ = _sobol_pure_numpy(fitted_params_full,
 n_samples=256, seed=seed)
 return Si, param_names
# ---------------------------------------------------------------------------
# 7. CROSS-VALIDATION
# ---------------------------------------------------------------------------
def cross_validate_mirzaaghasi(fitted_params_full):
 """
 Apply Choi-derived parameters to Mirzaaghasi et al. [17] dataset
 without refitting. Both observed and predicted values normalised to t=1 h.
 Returns:
 B_pred_norm, L_pred_norm : M3 normalised predictions
 B_m2_norm, L_m2_norm : reservoir-free M2 predictions
 mae_blood, mae_liver : M3 MAE (normalised scale)
 mae_liver_m2 : M2 liver MAE
 """
 t_eval = np.concatenate([[0.01], mirzaaghasi_time])
 sol = simulate_M3(fitted_params_full, t_eval)
 B_pred = sol.y[0, 1:]
 L_pred = sol.y[2, 1:] + sol.y[3, 1:]
 B_pred_norm = B_pred / B_pred[0]
 L_pred_norm = L_pred / L_pred[0]
 mae_blood = np.mean(np.abs(B_pred_norm - mirzaaghasi_blood_norm))
 mae_liver = np.mean(np.abs(L_pred_norm - mirzaaghasi_liver_norm))
 # Reservoir-free M2 for comparison
 p_m2 = fitted_params_full[[0, 1, 2, 3]]
 sol2 = simulate_M2(p_m2, t_eval)
 B_m2 = sol2.y[0, 1:] / sol2.y[0, 1]
 L_m2 = sol2.y[2, 1:] / sol2.y[2, 1]
 mae_liver_m2 = np.mean(np.abs(L_m2 - mirzaaghasi_liver_norm))
 return (B_pred_norm, L_pred_norm, B_m2, L_m2,
 mae_blood, mae_liver, mae_liver_m2)
def cross_validate_wiklander(fitted_params_full):
 """
 Apply Choi-derived parameters to Wiklander et al. [11] dataset
 without refitting. Both observed and predicted values normalised to t=1 h.
 Returns:
 B_pred_norm, L_pred_norm : M3 normalised predictions
 B_m2_norm, L_m2_norm : reservoir-free M2 predictions
 mae_blood, mae_liver : M3 MAE (normalised scale)
 mae_liver_m2 : M2 liver MAE
 """
 t_eval = np.concatenate([[0.01], wiklander_time])
 sol = simulate_M3(fitted_params_full, t_eval)
 B_pred = sol.y[0, 1:]
 L_pred = sol.y[2, 1:] + sol.y[3, 1:]
 B_pred_norm = B_pred / B_pred[0]
 L_pred_norm = L_pred / L_pred[0]
 mae_blood = np.mean(np.abs(B_pred_norm - wiklander_blood_norm))
 mae_liver = np.mean(np.abs(L_pred_norm - wiklander_liver_norm))
 p_m2 = fitted_params_full[[0, 1, 2, 3]]
 sol2 = simulate_M2(p_m2, t_eval)
 B_m2 = sol2.y[0, 1:] / sol2.y[0, 1]
 L_m2 = sol2.y[2, 1:] / sol2.y[2, 1]
 mae_liver_m2 = np.mean(np.abs(L_m2 - wiklander_liver_norm))
 return (B_pred_norm, L_pred_norm, B_m2, L_m2,
 mae_blood, mae_liver, mae_liver_m2)
# ---------------------------------------------------------------------------
# 8. RETENTION RATIO SIMULATIONS
# ---------------------------------------------------------------------------
def simulate_retention_ratio(R_values, fitted_params_full,
 t_max=120, n_points=500):
 """Vary R = k_bind/k_rel while fixing absolute k_rel; return hepatic profiles."""
 t_eval = np.linspace(0, t_max, n_points)
 k_rel_fixed = fitted_params_full[5]
 profiles = {}
 for R in R_values:
 p = fitted_params_full.copy()
 p[5] = k_rel_fixed
 p[4] = R * k_rel_fixed # k_bind = R * k_rel
 sol = simulate_M3(p, t_eval)
 profiles[R] = sol.y[2] + sol.y[3]
 return t_eval, profiles
def plateau_duration(liver_profile, t_eval, threshold=0.9):
 """Duration (h) during which hepatic level >= threshold * peak."""
 peak = np.max(liver_profile)
 above = t_eval[liver_profile >= threshold * peak]
 if len(above) < 2:
 return 0.0
 return above[-1] - above[0]
def scan_R_vs_plateau(fitted_params_full, R_range=None):
 """Scan R from 0.25 to 12, return (R_values, plateau_durations)."""
 if R_range is None:
 R_range = np.linspace(0.25, 12, 80)
 t_eval = np.linspace(0, 300, 1000)
 k_rel_fixed = fitted_params_full[5]
 durations = []
 for R in R_range:
 p = fitted_params_full.copy()
 p[5] = k_rel_fixed
 p[4] = R * k_rel_fixed
 sol = simulate_M3(p, t_eval)
 liver = sol.y[2] + sol.y[3]
 durations.append(plateau_duration(liver, t_eval))
 return R_range, np.array(durations)
def scan_2d_R_krel(fitted_params_full, R_vals=None, krel_vals=None):
 """2D scan of (R, k_rel) parameter space vs plateau duration."""
 if R_vals is None:
 R_vals = np.linspace(0.5, 12, 30)
 if krel_vals is None:
 krel_vals = np.linspace(0.005, 0.20, 30)
 t_eval = np.linspace(0, 300, 1000)
 grid = np.zeros((len(krel_vals), len(R_vals)))
 for i, kr in enumerate(krel_vals):
 for j, R in enumerate(R_vals):
 p = fitted_params_full.copy()
 p[5] = kr
 p[4] = R * kr
 sol = simulate_M3(p, t_eval)
 liver = sol.y[2] + sol.y[3]
 grid[i, j] = plateau_duration(liver, t_eval)
 return R_vals, krel_vals, grid
def simulate_krel_reduction(fitted_params_full,
 reductions=(0.25, 0.50, 0.75),
 t_max=300, n_points=1000):
 """
 Simulate the effect of proportional k_rel reductions on plateau duration.
 Returns dict keyed by reduction fraction with (t_eval, liver_profile, duration).
 """
 t_eval = np.linspace(0, t_max, n_points)
 k_rel_base = fitted_params_full[5]
 results = {}
 # Baseline
 sol_base = simulate_M3(fitted_params_full, t_eval)
 liver_base = sol_base.y[2] + sol_base.y[3]
 results[0.0] = (t_eval, liver_base, plateau_duration(liver_base, t_eval))
 for frac in reductions:
 p = fitted_params_full.copy()
 p[5] = k_rel_base * (1.0 - frac)
 sol = simulate_M3(p, t_eval)
 liver = sol.y[2] + sol.y[3]
 results[frac] = (t_eval, liver, plateau_duration(liver, t_eval))
 return results
# ---------------------------------------------------------------------------
# 9. PLOTTING
# ---------------------------------------------------------------------------
BLUE = '#2166ac'
RED = '#d6604d'
GREEN = '#4dac26'
PURPLE = '#762a83'
GREY = '#636363'
def plot_fig1_blood_M3(t_obs, blood_obs, blood_sd, fitted_params_full):
 """Fig 1: Blood exosome concentration – M3 v2 fit."""
 t_fine = np.linspace(0, 168, 500)
 sol = simulate_M3(fitted_params_full, t_fine)
 fig, ax = plt.subplots(figsize=(6, 4.5))
 ax.plot(t_fine, sol.y[0], color=BLUE, lw=2, label='Model (M3 v2)')
 ax.errorbar(t_obs, blood_obs, yerr=blood_sd,
 fmt='o', color=BLUE, capsize=4, label='Choi et al. [16]')
 ax.set_xlabel('Time (h)')
 ax.set_ylabel('Blood exosome concentration (%ID)')
 ax.legend()
 ax.set_xlim(-2, 170)
 ax.set_ylim(-2, 65)
 plt.tight_layout()
 return fig
def plot_fig2_liver_M3(t_obs, liver_obs, liver_sd, fitted_params_full):
 """Fig 2: Hepatic exosome concentration – M3 v2 fit (intact vesicles)."""
 t_fine = np.linspace(0, 168, 500)
 sol = simulate_M3(fitted_params_full, t_fine)
 liver_total = sol.y[2] + sol.y[3]
 fig, ax = plt.subplots(figsize=(6, 4.5))
 ax.plot(t_fine, liver_total, color=RED, lw=2,
 label='Model (M3 v2, intact vesicles)')
 ax.errorbar(t_obs, liver_obs, yerr=liver_sd,
 fmt='s', color=RED, capsize=4,
 label='Choi et al. [16] (99Zr signal)')
 ax.set_xlabel('Time (h)')
 ax.set_ylabel('Hepatic exosome concentration (%ID)')
 ax.legend()
 ax.set_xlim(-2, 170)
 ax.set_ylim(-2, 65)
 plt.tight_layout()
 return fig
def plot_fig3_intracellular(fitted_params_full):
 """Fig 3: Intracellular uptake rate over time."""
 t_fine = np.linspace(0, 168, 500)
 sol = simulate_M3(fitted_params_full, t_fine)
 L_f = sol.y[2]
 V_max = fitted_params_full[6]
 K_m = fitted_params_full[7]
 uptake_rate = V_max * L_f / (K_m + L_f)
 fig, ax = plt.subplots(figsize=(6, 4.5))
 ax.plot(t_fine, uptake_rate, color=GREEN, lw=2)
 ax.set_xlabel('Time (h)')
 ax.set_ylabel('Intracellular uptake rate (%ID·h⁻¹)')
 ax.set_xlim(-2, 170)
 plt.tight_layout()
 return fig
def plot_fig4_drug(fitted_params_full):
 """Fig 4: Local drug concentration."""
 t_fine = np.linspace(0, 168, 500)
 sol = simulate_M3(fitted_params_full, t_fine)
 fig, ax = plt.subplots(figsize=(6, 4.5))
 ax.plot(t_fine, sol.y[5], color=PURPLE, lw=2)
 ax.set_xlabel('Time (h)')
 ax.set_ylabel('Local drug concentration (a.u.)')
 ax.set_xlim(-2, 170)
 plt.tight_layout()
 return fig
def plot_fig5_effect(fitted_params_full):
 """Fig 5: Normalised therapeutic effect."""
 t_fine = np.linspace(0, 168, 500)
 sol = simulate_M3(fitted_params_full, t_fine)
 C_drug = sol.y[5]
 EC50 = fitted_params_full[11]
 effect = C_drug / (C_drug + EC50)
 fig, ax = plt.subplots(figsize=(6, 4.5))
 ax.plot(t_fine, effect, color=GREY, lw=2)
 ax.set_xlabel('Time (h)')
 ax.set_ylabel('Normalized therapeutic effect')
 ax.set_xlim(-2, 170)
 ax.set_ylim(-0.02, 0.85)
 plt.tight_layout()
 return fig
def plot_fig6_model_comparison(t_obs, blood_obs, blood_sd,
 res_M1, res_M2, fitted_params_full):
 """Fig 6: Blood concentration – all three models overlaid."""
 t_fine = np.linspace(0, 168, 500)
 sol3 = simulate_M3(fitted_params_full, t_fine)
 sol2 = simulate_M2(res_M2.x, t_fine)
 sol1 = simulate_M1(res_M1.x, t_fine)
 fig, ax = plt.subplots(figsize=(6, 4.5))
 ax.plot(t_fine, sol1.y[0], color=GREY, lw=1.5, ls=':', label='M1 (1-cmpt)')
 ax.plot(t_fine, sol2.y[0], color=RED, lw=1.5, ls='--', label='M2 (2-cmpt)')
 ax.plot(t_fine, sol3.y[0], color=BLUE, lw=2, label='M3 (reservoir)')
 ax.errorbar(t_obs, blood_obs, yerr=blood_sd,
 fmt='o', color='black', capsize=4, label='Choi et al. [16]')
 ax.set_xlabel('Time (h)')
 ax.set_ylabel('Blood exosome concentration (%ID)')
 ax.legend(fontsize=9)
 ax.set_xlim(-2, 170)
 ax.set_ylim(-2, 65)
 plt.tight_layout()
 return fig
def plot_fig7_aicc(aicc_m1, aicc_m2, aicc_m3):
 """Fig 7: ΔAICc bar chart (relative to M3)."""
 delta = np.array([aicc_m1 - aicc_m3, aicc_m2 - aicc_m3, 0.0])
 labels = ['M1\n(1-cmpt)', 'M2\n(2-cmpt)', 'M3*\n(reservoir)']
 colors = [GREY, RED, BLUE]
 fig, ax = plt.subplots(figsize=(5, 4.5))
 bars = ax.bar(labels, delta, color=colors, width=0.5, edgecolor='none')
 for bar, val in zip(bars, delta):
 ax.text(bar.get_x() + bar.get_width() / 2,
 bar.get_height() + 0.3, f'{val:.1f}',
 ha='center', va='bottom', fontweight='bold')
 ax.axhline(2, color='black', lw=1.2, ls='--', label='ΔAICc = 2')
 ax.axhline(10, color='black', lw=1.2, ls=':', label='ΔAICc = 10')
 ax.set_ylabel('ΔAICc (relative to M3)')
 ax.legend(fontsize=8)
 ax.set_ylim(0, max(delta) * 1.15 + 1)
 plt.tight_layout()
 return fig
def plot_fig8_bootstrap(boot_params):
 """Fig 8: Bootstrap distribution of k_rel."""
 k_rel_boot = boot_params[:, 5]
 point_est = np.median(k_rel_boot)
 ci_lo, ci_hi = np.percentile(k_rel_boot, [2.5, 97.5])
 fig, ax = plt.subplots(figsize=(6, 4.5))
 ax.hist(k_rel_boot, bins=25, color=BLUE, edgecolor='white', alpha=0.9)
 ax.axvline(point_est, color='black', lw=2, label='Point estimate')
 ax.axvline(ci_lo, color=RED, lw=1.5, ls='--', label='95% CI')
 ax.axvline(ci_hi, color=RED, lw=1.5, ls='--')
 ax.set_xlabel('k_rel (h⁻¹)')
 ax.set_ylabel('Frequency')
 ax.legend()
 plt.tight_layout()
 return fig
def plot_fig9_sobol(Si, param_names):
 """Fig 9: Sobol first- and total-order sensitivity indices."""
 S1 = Si['S1']
 ST = Si['ST']
 y_pos = np.arange(len(param_names))
 fig, ax = plt.subplots(figsize=(6, 4))
 ax.barh(y_pos + 0.2, ST, height=0.35, color=RED, label='Total-order (Sₜ)')
 ax.barh(y_pos - 0.2, S1, height=0.35, color=BLUE, label='First-order (S₁)')
 ax.set_yticks(y_pos)
 ax.set_yticklabels(param_names)
 ax.set_xlabel('Sobol sensitivity index')
 ax.legend()
 plt.tight_layout()
 return fig
def plot_fig10_R_plateau(R_vals, durations, fitted_R):
 """Fig 10: Retention ratio R vs plateau duration (1D scan)."""
 fig, ax = plt.subplots(figsize=(6, 4.5))
 ax.plot(R_vals, durations, color=BLUE, lw=2)
 ax.axvline(fitted_R, color=RED, lw=1.5, ls='--',
 label=f'Fitted R ≈ {fitted_R:.1f}')
 ax.set_xlabel('Retention ratio R (= k_bind / k_rel)')
 ax.set_ylabel('Plateau duration (h)')
 ax.legend()
 plt.tight_layout()
 return fig
def plot_fig11_hepatic_profiles(t_eval, profiles):
 """Fig 11: Hepatic exosome level for R = 1, 3, 6, 10."""
 colors_R = {1: BLUE, 3: RED, 6: GREEN, 10: PURPLE}
 fig, ax = plt.subplots(figsize=(6, 4.5))
 for R, color in colors_R.items():
 if R in profiles:
 ax.plot(t_eval, profiles[R], color=color, lw=2, label=f'R = {R}')
 ax.set_xlabel('Time (h)')
 ax.set_ylabel('Hepatic exosome level (a.u.)')
 ax.legend(title='Retention ratio R')
 plt.tight_layout()
 return fig
def plot_fig12_2d_heatmap(R_vals, krel_vals, grid):
 """Fig 12: 2D heatmap of plateau duration over (R, k_rel) space."""
 fig, ax = plt.subplots(figsize=(7, 5))
 cf = ax.contourf(R_vals, krel_vals, grid, levels=20, cmap='viridis')
 cbar = fig.colorbar(cf, ax=ax)
 cbar.set_label('Plateau duration (h)')
 ax.set_xlabel('Retention ratio R (= k_bind / k_rel)')
 ax.set_ylabel('Release rate k_rel (h⁻¹)')
 plt.tight_layout()
 return fig
def plot_fig13_schematic():
 """
 Fig 13: Conceptual comparison of classical PK and reservoir model structures.
 """
 fig, axes = plt.subplots(1, 2, figsize=(10, 5))
 for ax in axes:
 ax.set_xlim(0, 10)
 ax.set_ylim(0, 10)
 ax.axis('off')
 # --- Classical PK model (left panel) ---
 ax = axes[0]
 ax.text(5, 9.3, 'Classical PK model', ha='center', fontsize=12,
 color='grey', fontweight='bold')
 ax.add_patch(plt.FancyBboxPatch((1, 6.5), 4, 1.8,
 boxstyle='round,pad=0.1', fc='#AED6F1', ec=BLUE, lw=2))
 ax.text(3, 7.4, 'Blood (B)', ha='center', fontsize=11, fontweight='bold')
 ax.add_patch(plt.FancyBboxPatch((5.5, 6.5), 3.5, 1.8,
 boxstyle='round,pad=0.1', fc='#D5D8DC', ec='grey', lw=1.5))
 ax.text(7.25, 7.4, 'Peripheral (P)', ha='center', fontsize=10)
 ax.annotate('', xy=(5.5, 7.5), xytext=(5.05, 7.5),
 arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
 ax.annotate('', xy=(5.0, 7.1), xytext=(5.5, 7.1),
 arrowprops=dict(arrowstyle='->', color='black', lw=1.5))
 ax.add_patch(plt.FancyBboxPatch((1, 3.5), 4, 1.8,
 boxstyle='round,pad=0.1', fc='#FADBD8', ec=RED, lw=1.5))
 ax.text(3, 4.4, 'Organ\n(no reservoir)', ha='center', fontsize=10)
 ax.annotate('', xy=(3, 5.3), xytext=(3, 6.5),
 arrowprops=dict(arrowstyle='->', color=RED, lw=1.5))
 ax.text(3.2, 5.9, 'k_to', color=RED, fontsize=9)
 ax.annotate('', xy=(3, 2.5), xytext=(3, 3.5),
 arrowprops=dict(arrowstyle='->', color=RED, lw=1.5))
 ax.text(3.2, 3.0, 'clearance', color=RED, fontsize=9)
 ax.text(5.5, 2.5, '→ tracks blood decay', style='italic',
 fontsize=9, color='grey')
 # --- Reservoir model (right panel) ---
 ax = axes[1]
 ax.text(5, 9.3, 'Reservoir model (proposed)', ha='center', fontsize=12,
 color=BLUE, fontweight='bold')
 ax.add_patch(plt.FancyBboxPatch((1, 6.5), 4, 1.8,
 boxstyle='round,pad=0.1', fc='#AED6F1', ec=BLUE, lw=2))
 ax.text(3, 7.4, 'Blood (B)', ha='center', fontsize=11, fontweight='bold')
 ax.add_patch(plt.FancyBboxPatch((1, 3.5), 2, 1.8,
 boxstyle='round,pad=0.1', fc='#D5F5E3', ec=GREEN, lw=1.5))
 ax.text(2, 4.4, 'L_free', ha='center', fontsize=10)
 ax.add_patch(plt.FancyBboxPatch((3.2, 3.5), 2, 1.8,
 boxstyle='round,pad=0.1', fc='#FADBD8', ec=RED, lw=1.5))
 ax.text(4.2, 4.4, 'L_bound', ha='center', fontsize=10)
 ax.annotate('', xy=(3.2, 4.6), xytext=(3.0, 4.6),
 arrowprops=dict(arrowstyle='->', color=GREEN, lw=1.5))
 ax.text(2.8, 4.9, 'k_bind', color=GREEN, fontsize=8)
 ax.annotate('', xy=(3.0, 4.2), xytext=(3.2, 4.2),
 arrowprops=dict(arrowstyle='->', color=RED, lw=1.5))
 ax.text(2.85, 3.9, 'k_rel', color=RED, fontsize=8)
 ax.annotate('', xy=(2, 5.3), xytext=(2.5, 6.5),
 arrowprops=dict(arrowstyle='->', color=BLUE, lw=1.5))
 ax.text(1.5, 5.9, 'k_to', color=BLUE, fontsize=9)
 ax.add_patch(plt.FancyBboxPatch((5.8, 5.5), 3.8, 2.5,
 boxstyle='round,pad=0.15', fc='#EBF5FB', ec=BLUE, lw=2))
 ax.text(7.7, 7.5, 'Retention ratio', ha='center', fontsize=10,
 fontweight='bold', color=BLUE)
 ax.text(7.7, 6.9, 'R = k_bind / k_rel', ha='center', fontsize=11,
 fontweight='bold', style='italic', color=BLUE)
 ax.text(7.7, 6.3, 'R/(1+R) = bound fraction', ha='center', fontsize=9,
 color='grey')
 ax.text(2.5, 2.5, '→ sustained plateau', style='italic',
 fontsize=9, color=BLUE)
 plt.tight_layout()
 return fig
def plot_fig14_crossval(fitted_params_full):
 """
 Fig 14: Independent cross-validation across two datasets (no parameter refitting).
 2×2 panel: top = blood (Mirzaaghasi | Wiklander),
 bottom = liver (Mirzaaghasi | Wiklander).
 """
 (B_mirz, L_mirz, B_mirz_m2, L_mirz_m2,
 mae_b_mirz, mae_l_mirz, _) = cross_validate_mirzaaghasi(fitted_params_full)
 (B_wikl, L_wikl, B_wikl_m2, L_wikl_m2,
 mae_b_wikl, mae_l_wikl, _) = cross_validate_wiklander(fitted_params_full)
 fig, axes = plt.subplots(2, 2, figsize=(10, 7))
 fig.suptitle('Independent cross-validation across two datasets'
 ' (no parameter refitting)', fontsize=11)
 # ---------- Mirzaaghasi – blood ----------
 ax = axes[0, 0]
 ax.set_title('Mirzaaghasi et al. [17]\n(HEK293T, DiR, 1–8 h)', fontsize=9)
 ax.plot(mirzaaghasi_time, B_mirz, '-o', color=BLUE, label='M3 prediction')
 ax.scatter(mirzaaghasi_time, mirzaaghasi_blood_norm,
 color='black', zorder=5, label='Observed')
 ax.set_ylabel('Norm. blood signal')
 ax.text(0.97, 0.95, f'MAE={mae_b_mirz:.3f}',
 transform=ax.transAxes, ha='right', va='top', fontsize=8)
 ax.legend(fontsize=8)
 # ---------- Wiklander – blood ----------
 ax = axes[0, 1]
 ax.set_title('Wiklander et al. [11]\n(C2C12, DiI, 1–24 h)', fontsize=9)
 ax.plot(wiklander_time, B_wikl, '-o', color=BLUE, label='M3 prediction')
 ax.scatter(wiklander_time, wiklander_blood_norm,
 color='black', zorder=5, label='Observed')
 ax.set_ylabel('Norm. blood signal')
 ax.text(0.97, 0.95, f'MAE={mae_b_wikl:.3f}',
 transform=ax.transAxes, ha='right', va='top', fontsize=8)
 ax.legend(fontsize=8)
 # ---------- Mirzaaghasi – liver ----------
 ax = axes[1, 0]
 ax.plot(mirzaaghasi_time, L_mirz, '-s', color=RED, label='M3 prediction')
 ax.scatter(mirzaaghasi_time, mirzaaghasi_liver_norm,
 color='darkred', marker='s', zorder=5, label='Observed')
 ax.set_xlabel('Time (h)')
 ax.set_ylabel('Norm. liver signal')
 ax.text(0.97, 0.05, f'MAE={mae_l_mirz:.3f}',
 transform=ax.transAxes, ha='right', va='bottom', fontsize=8)
 ax.legend(fontsize=8)
 # ---------- Wiklander – liver ----------
 ax = axes[1, 1]
 ax.plot(wiklander_time, L_wikl, '-s', color=RED, label='M3 prediction')
 ax.scatter(wiklander_time, wiklander_liver_norm,
 color='darkred', marker='s', zorder=5, label='Observed')
 ax.set_xlabel('Time (h)')
 ax.set_ylabel('Norm. liver signal')
 ax.text(0.97, 0.05, f'MAE={mae_l_wikl:.3f}',
 transform=ax.transAxes, ha='right', va='bottom', fontsize=8)
 ax.legend(fontsize=8)
 plt.tight_layout()
 return fig
def plot_fig15_krel_reduction(krel_results):
 """
 Fig 15: Effect of progressive k_rel reduction on hepatic profile.
 Simulates 0%, 25%, 50%, 75% reduction relative to fitted k_rel.
 """
 labels = {
 0.00: ('Baseline', BLUE),
 0.25: ('−25% k_rel', GREEN),
 0.50: ('−50% k_rel', RED),
 0.75: ('−75% k_rel', PURPLE),
 }
 fig, ax = plt.subplots(figsize=(6, 4.5))
 for frac, (label, color) in labels.items():
 if frac in krel_results:
 t_eval, liver, dur = krel_results[frac]
 ax.plot(t_eval, liver, color=color, lw=2,
 label=f'{label} (plateau ≈ {dur:.0f} h)')
 ax.set_xlabel('Time (h)')
 ax.set_ylabel('Hepatic exosome level (a.u.)')
 ax.legend(fontsize=8, title='k_rel modification')
 plt.tight_layout()
 return fig
# ---------------------------------------------------------------------------
# 10. MAIN EXECUTION
# ---------------------------------------------------------------------------
def main():
 print("=" * 60)
 print("Exosome PK Reservoir Model – Full Analysis (v2)")
 print("=" * 60)
 t_obs = choi_time
 blood_obs = choi_blood_mean
 liver_obs = choi_liver_mean
 blood_sd = choi_blood_sd
 liver_sd = choi_liver_sd
 n_obs = 2 * len(t_obs) # blood + liver data points
 # ------------------------------------------------------------------
 # 10.1 Fit all three models
 # ------------------------------------------------------------------
 print("\n[1/8] Fitting models M1, M2, M3 ...")
 res_M1 = fit_M1(t_obs, blood_obs, blood_sd)
 res_M2 = fit_M2(t_obs, blood_obs, liver_obs, blood_sd, liver_sd)
 res_M3 = fit_M3(t_obs, blood_obs, liver_obs, blood_sd, liver_sd)
 fitted_params_full = P_NOMINAL.copy()
 fitted_params_full[FIT_IDX] = res_M3.x
 pnames = ['k_clear', 'k_to', 'k_bp', 'k_pb', 'k_bind', 'k_rel']
 print(" M3 fitted parameters:")
 for name, val in zip(pnames, res_M3.x):
 print(f" {name:10s} = {val:.4f} h⁻¹")
 fitted_R = res_M3.x[4] / res_M3.x[5]
 print(f" Retention ratio R = k_bind/k_rel = {fitted_R:.2f}")
 # ------------------------------------------------------------------
 # 10.2 RMSE
 # ------------------------------------------------------------------
 sol3 = simulate_M3(fitted_params_full, t_obs)
 sol2 = simulate_M2(res_M2.x, t_obs)
 sol1 = simulate_M1(res_M1.x, t_obs)
 rmse_M3 = np.sqrt(np.mean(np.concatenate([
 sol3.y[0] - blood_obs,
 (sol3.y[2] + sol3.y[3]) - liver_obs]) ** 2))
 rmse_M2 = np.sqrt(np.mean(np.concatenate([
 sol2.y[0] - blood_obs,
 sol2.y[2] - liver_obs]) ** 2))
 rmse_M1 = np.sqrt(np.mean((sol1.y[0] - blood_obs) ** 2))
 print(f"\n RMSE – M1: {rmse_M1:.2f}, M2: {rmse_M2:.2f}, M3: {rmse_M3:.2f} %ID")
 # ------------------------------------------------------------------
 # 10.3 AICc model comparison
 # ------------------------------------------------------------------
 print("\n[2/8] AICc model comparison ...")
 res_M3_resid = residuals_M3(res_M3.x, P_NOMINAL.copy(),
 t_obs, blood_obs, liver_obs, blood_sd, liver_sd)
 res_M2_resid = residuals_M2(res_M2.x, t_obs, blood_obs, liver_obs,
 blood_sd, liver_sd)
 res_M1_resid = residuals_M1(res_M1.x, t_obs, blood_obs, blood_sd)
 aicc_M3 = compute_aicc(res_M3_resid, 6, n_obs)
 aicc_M2 = compute_aicc(res_M2_resid, 4, n_obs)
 aicc_M1 = compute_aicc(res_M1_resid, 1, len(t_obs))
 print(f" AICc – M1: {aicc_M1:.1f}, M2: {aicc_M2:.1f}, M3: {aicc_M3:.1f}")
 print(f" ΔAICc – M1 vs M3: {aicc_M1-aicc_M3:.1f}, "
 f"M2 vs M3: {aicc_M2-aicc_M3:.1f}")
 # ------------------------------------------------------------------
 # 10.4 Bootstrap
 # ------------------------------------------------------------------
 print("\n[3/8] Bootstrap parameter uncertainty (n=500) ...")
 print(" Note: set n_boot=500 for publication; using 100 here for speed.")
 boot_params = bootstrap_M3(
 t_obs, blood_obs, liver_obs, blood_sd, liver_sd,
 fitted_params_full, n_boot=100
 )
 print(f" Successful bootstrap iterations: {len(boot_params)}")
 if len(boot_params) > 10:
 for i, name in enumerate(pnames):
 vals = boot_params[:, i]
 cv = np.std(vals) / np.mean(vals) * 100
 ci = np.percentile(vals, [2.5, 97.5])
 print(f" {name:10s}: median={np.median(vals):.4f}, "
 f"95%CI=[{ci[0]:.4f},{ci[1]:.4f}], CV={cv:.1f}%")
 # ------------------------------------------------------------------
 # 10.5 Sobol sensitivity
 # ------------------------------------------------------------------
 print("\n[4/8] Sobol global sensitivity analysis ...")
 print(" Note: install SALib (pip install SALib) for publication-quality results.")
 Si, param_names_sobol = sobol_analysis(fitted_params_full, n_samples=256)
 for name, s1, st in zip(param_names_sobol, Si['S1'], Si['ST']):
 print(f" {name:10s}: S1={s1:.3f}, ST={st:.3f}")
 # ------------------------------------------------------------------
 # 10.6 Cross-validation (both datasets)
 # ------------------------------------------------------------------
 print("\n[5/8] Cross-validation ...")
 print(" Dataset 1: Mirzaaghasi et al. [17] (HEK293T, DiR, 1–8 h)")
 (B_mirz, L_mirz, _, _,
 mae_b_mirz, mae_l_mirz, mae_lm2_mirz) = cross_validate_mirzaaghasi(
 fitted_params_full)
 print(f" MAE (M3) – blood: {mae_b_mirz:.3f}, liver: {mae_l_mirz:.3f}")
 print(f" MAE (M2) – liver: {mae_lm2_mirz:.3f}")
 print(" Dataset 2: Wiklander et al. [11] (C2C12, DiI, 1–24 h)")
 (B_wikl, L_wikl, _, _,
 mae_b_wikl, mae_l_wikl, mae_lm2_wikl) = cross_validate_wiklander(
 fitted_params_full)
 print(f" MAE (M3) – blood: {mae_b_wikl:.3f}, liver: {mae_l_wikl:.3f}")
 print(f" MAE (M2) – liver: {mae_lm2_wikl:.3f}")
 print(f"\n Summary MAE ranges (M3):")
 print(f" Blood: {min(mae_b_mirz, mae_b_wikl):.3f}–"
 f"{max(mae_b_mirz, mae_b_wikl):.3f}")
 print(f" Liver: {min(mae_l_mirz, mae_l_wikl):.3f}–"
 f"{max(mae_l_mirz, mae_l_wikl):.3f}")
 # ------------------------------------------------------------------
 # 10.7 Retention ratio simulations
 # ------------------------------------------------------------------
 print("\n[6/8] Retention ratio simulations ...")
 t_profiles, profiles = simulate_retention_ratio(
 [1, 3, 6, 10], fitted_params_full)
 R_scan, durations = scan_R_vs_plateau(fitted_params_full)
 R_2d, krel_2d, grid_2d = scan_2d_R_krel(fitted_params_full)
 print("\n[7/8] k_rel reduction simulations ...")
 krel_results = simulate_krel_reduction(fitted_params_full)
 for frac, (_, _, dur) in krel_results.items():
 label = f"−{int(frac*100)}%" if frac > 0 else "baseline"
 print(f" k_rel {label:>10s}: plateau ≈ {dur:.0f} h")
 # ------------------------------------------------------------------
 # 10.8 Generate all figures
 # ------------------------------------------------------------------
 print("\n[8/8] Generating figures ...")
 figs = {}
 figs[1] = plot_fig1_blood_M3(t_obs, blood_obs, blood_sd, fitted_params_full)
 figs[2] = plot_fig2_liver_M3(t_obs, liver_obs, liver_sd, fitted_params_full)
 figs[3] = plot_fig3_intracellular(fitted_params_full)
 figs[4] = plot_fig4_drug(fitted_params_full)
 figs[5] = plot_fig5_effect(fitted_params_full)
 figs[6] = plot_fig6_model_comparison(t_obs, blood_obs, blood_sd,
 res_M1, res_M2, fitted_params_full)
 figs[7] = plot_fig7_aicc(aicc_M1, aicc_M2, aicc_M3)
 if len(boot_params) > 10:
 figs[8] = plot_fig8_bootstrap(boot_params)
 figs[9] = plot_fig9_sobol(Si, param_names_sobol)
 figs[10] = plot_fig10_R_plateau(R_scan, durations, fitted_R)
 figs[11] = plot_fig11_hepatic_profiles(t_profiles, profiles)
 figs[12] = plot_fig12_2d_heatmap(R_2d, krel_2d, grid_2d)
 figs[13] = plot_fig13_schematic()
 figs[14] = plot_fig14_crossval(fitted_params_full)
 figs[15] = plot_fig15_krel_reduction(krel_results)
 for num, fig in figs.items():
 fname = f'Fig{num:02d}.png'
 fig.savefig(fname, dpi=200, bbox_inches='tight')
 plt.close(fig)
 print(f" Saved {fname}")
 print("\nDone. All figures saved to current directory.")
 print("=" * 60)
if __name__ == '__main__':
 main()