import numpy as np
from scipy.integrate import quad

def beta_function(z, beta_inner=3.0, beta_outer=1.0, damp=2.0):
    return beta_outer + (beta_inner - beta_outer) * np.exp(-z / damp)

def c_eff(z, beta_z):
    return 3e5 * (1 + z**beta_z)**(-0.5)  # km/s

def H_z(z, H0=70, Omega_b=0.05):
    beta_z = beta_function(z)
    dilation = 1 + z**beta_z
    c_ratio = 3e5 / c_eff(z, beta_z)
    return H0 * np.sqrt(Omega_b * (1+z)**3 + (1-Omega_b) * dilation) * c_ratio

def D_L(z, H0=70):
    integral = quad(lambda zp: 1 / H_z(zp, H0), 0, z)[0]
    return (1+z) * integral * (3e5 / H0)  # Mpc

print("D_L at z=1:", D_L(1))  # ~44.1 expected for μ