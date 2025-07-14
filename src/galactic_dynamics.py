import numpy as np

G = 6.67430e-11  # m^3 kg^-1 s^-2
kpc_to_m = 3.086e19
solar_mass = 1.989e30
km_s = 1e3

def beta_function(r, beta_inner=3.0, beta_outer=2.0, R_damp=5.0):
    return beta_outer + (beta_inner - beta_outer) * np.exp(-r / R_damp)

def v_rotation_UDT(r, M_star, R0=20.0, alpha=0.125):
    r = np.array(r)
    valid = r > 0
    if not np.any(valid):
        return np.zeros(len(r))
    r_valid = r[valid]
    beta_r = 2.0 + 1.0 * np.exp(-r_valid / (R0 / np.e))
    G_astro = 4.302e-6  # kpc (km/s)^2 / M_sun
    v_star_sq = alpha * G_astro * M_star / r_valid
    D_factor = np.sqrt(1 + (r_valid / R0)**beta_r)
    v_total = np.sqrt(v_star_sq) * D_factor
    v_out = np.zeros(len(r))
    v_out[valid] = v_total
    return v_out

# Example: UGCA442
r = np.array([0.42, 1.26, 2.11, 2.96, 3.79, 4.65, 5.48, 6.33])  # kpc
M_star = 5e7  # solar masses
v_pred = v_rotation_UDT(r, M_star)
print(v_pred)  # Compare to observed [14.2, 28.6, 41.0, 49.0, 54.8, 56.4, 57.8, 56.5]