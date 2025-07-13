import numpy as np

G = 6.67430e-11  # m^3 kg^-1 s^-2
kpc_to_m = 3.086e19
solar_mass = 1.989e30
km_s = 1e3

def beta_function(r, beta_inner=3.0, beta_outer=2.0, R_damp=5.0):
    return beta_outer + (beta_inner - beta_outer) * np.exp(-r / R_damp)

def v_rotation_UDT(r, M_star, R0=20.0, alpha=0.125):
    r_m = r * kpc_to_m
    R0_m = R0 * kpc_to_m
    M = M_star * solar_mass
    v_star_squared = alpha * G * M / r_m
    beta_r = beta_function(r, R_damp=R0 / np.e)
    D_factor = np.sqrt(1 + (r / R0)**beta_r)
    v = np.sqrt(v_star_squared) * D_factor / km_s  # km/s
    return v

# Example: UGCA442
r = np.array([0.42, 1.26, 2.11, 2.96, 3.79, 4.65, 5.48, 6.33])  # kpc
M_star = 5e7  # solar masses
v_pred = v_rotation_UDT(r, M_star)
print(v_pred)  # Compare to observed [14.2, 28.6, 41.0, 49.0, 54.8, 56.4, 57.8, 56.5]