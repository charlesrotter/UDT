import sympy as sp
import numpy as np

# Symbols
r, M, G, c, R0, beta = sp.symbols('r M G c R0 beta')

# Metric g_rr inverse for precession approximation
g_rr_inv = 1 - 2*G*M/(c**2 * r) - (r/R0)**beta

# Precession δφ per orbit ~ 6π G M / (c^2 a (1-e^2)) + correction from β term
def precession_UDT(a, e, M_sun=1.989e30, G=6.6743e-11, c=3e8, R0=1e12, beta=2.5):
    GR_term = 6 * np.pi * G * M_sun / (c**2 * a * (1 - e**2))
    beta_term = (a / R0)**beta * np.pi  # Approximate dilation correction
    return GR_term + beta_term  # arcsec/century

# Mercury: a=5.79e10 m, e=0.2056
delta_phi = precession_UDT(5.79e10, 0.2056)
print("UDT precession for Mercury:", delta_phi, "arcsec/century")  # Should ~43 to match GR/data
