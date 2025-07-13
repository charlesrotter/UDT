import sympy as sp

# Define symbols
r, t, theta, phi, M, G, c, R0, beta = sp.symbols('r t theta phi M G c R0 beta')
dt, dr, dtheta, dphi = sp.symbols('dt dr dtheta dphi')

# Emergent c_eff (simplified for galactic scales)
c_eff = sp.sqrt(G * M / r) * (1 + (r / R0)**beta)**(-sp.Rational(1,2))

# Metric components
g_tt = -c_eff**2
g_rr = 1 / (1 - 2 * G * M / (c_eff**2 * r) - (r / R0)**beta)
g_theta_theta = r**2
g_phi_phi = r**2 * sp.sin(theta)**2

# Line element ds²
ds2 = g_tt * dt**2 + g_rr * dr**2 + g_theta_theta * dtheta**2 + g_phi_phi * dphi**2

# Print for verification
sp.pprint(ds2)
