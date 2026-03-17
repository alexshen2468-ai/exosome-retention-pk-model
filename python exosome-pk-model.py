import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# -----------------------------
# Mechanistic PK model (ODE)
# -----------------------------
def pk_model(y, t, k_in, k_out, k_bind, k_rel):
    C_blood, C_tissue, C_bound = y

    # Blood compartment
    dC_blood_dt = -k_in * C_blood - k_out * C_blood

    # Free tissue compartment
    dC_tissue_dt = k_in * C_blood - k_bind * C_tissue + k_rel * C_bound

    # Bound tissue compartment (retention)
    dC_bound_dt = k_bind * C_tissue - k_rel * C_bound

    return [dC_blood_dt, dC_tissue_dt, dC_bound_dt]


# -----------------------------
# Simulation function
# -----------------------------
def simulate(R):
    # Base parameters
    k_in = 0.1        # blood → tissue
    k_out = 0.05      # clearance from blood
    
    # Retention ratio
    k_rel = 0.01
    k_bind = R * k_rel   # R = k_bind / k_rel

    # Initial conditions
    y0 = [100, 0, 0]  # initial dose in blood

    # Time span
    t = np.linspace(0, 120, 500)

    # Solve ODE
    sol = odeint(pk_model, y0, t, args=(k_in, k_out, k_bind, k_rel))

    C_blood = sol[:, 0]
    C_tissue = sol[:, 1]
    C_bound = sol[:, 2]

    # Total hepatic level = free + bound
    hepatic_level = C_tissue + C_bound

    return t, hepatic_level


# -----------------------------
# Plotting
# -----------------------------
R_values = [1, 3, 6, 10]

plt.figure(figsize=(8, 5))

for R in R_values:
    t, y = simulate(R)
    plt.plot(t, y, linewidth=2, label=f"R = {R}")

# Labels
plt.xlabel("Time (h)")
plt.ylabel("Hepatic exosome level (a.u.)")
plt.title("Effect of Retention Ratio on Exosome Dynamics")

plt.legend(title="Retention ratio R")

plt.xlim(0, 120)

plt.tight_layout()

# Save figure for GitHub README
plt.savefig("Fig_retention_curves.png", dpi=300)

plt.show()