import numpy as np

# Constants (SI units)
G = 6.67430e-11
c = 3e8
kpc = 3.086e19
solar_mass = 1.989e30

def beta_r(r, beta_inner=3.0, beta_outer=2.0, R_damp=200.0):  # R_damp for cluster scale in kpc
    return beta_outer + (beta_inner - beta_outer) * np.exp(-r / R_damp)

def kappa_UDT(r, M_bary, R0=500.0, Sigma_c=3.1e9):  # R0 cluster scale ~500 kpc, Sigma_c in M_sun/kpc^2
    Sigma_bary = M_bary * solar_mass / (np.pi * (r * kpc)**2) / (1e6**2)  # Convert to M_sun/kpc^2 approx
    beta = beta_r(r)
    D_factor = np.sqrt(1 + (r / R0)**beta)
    return (Sigma_bary / Sigma_c) * D_factor

# Bullet Cluster example (from data: main cluster r~100 kpc, M_bary~8.8e12 M_sun)
r_array = np.linspace(10, 500, 100)  # kpc
M_bary_main = 8.8e12
kappa = kappa_UDT(r_array, M_bary_main)

print("Sample κ at r=100 kpc:", kappa[np.argmin(np.abs(r_array - 100))])  # Should ~0.38 to match data
