import numpy as np
import os
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl

# Constants
a = b = 0.5
h = 0.004
E = 200e9
nu = 0.3
q0 = 1000
yield_stress = 450e6

# Calculations
def compute_coefficients(N):
    coeffs = np.zeros((N, N))
    for m in range(1, N + 1):
        for n in range(1, N + 1):
            if m % 2 == 1 and n % 2 == 1:
                coeffs[m-1, n-1] = (16 * q0 * a**4) / (D * np.pi**6 * m * n * (m**2 + n**2)**2)
    return coeffs

def compute_fields(x, y, coeffs):
    N = coeffs.shape[0]
    w, Mx, My, Mxy = 0, 0, 0, 0
    for m in range(1, N + 1):
        for n in range(1, N + 1):
            c = coeffs[m-1, n-1]
            if c == 0:
                continue
            sin_mx = np.sin(m * np.pi * x / a)
            sin_ny = np.sin(n * np.pi * y / b)
            cos_mx = np.cos(m * np.pi * x / a)
            cos_ny = np.cos(n * np.pi * y / b)
            w += c * sin_mx * sin_ny
            Mx += -D * ((m * np.pi / a)**2 + nu * (n * np.pi / b)**2) * c * sin_mx * sin_ny
            My += -D * (nu * (m * np.pi / a)**2 + (n * np.pi / b)**2) * c * sin_mx * sin_ny
            Mxy += (1 - nu) * D * (m * n * np.pi**2 / (a * b)) * c * cos_mx * cos_ny
    return w, Mx, My, Mxy

x = np.linspace(0, a, 50)
y = np.linspace(0, b, 50)
X, Y = np.meshgrid(x, y)
z = np.linspace(-h / 2, h / 2, 100)
center_x, center_y = a / 2, b / 2

plt.figure(figsize=(10, 6))

strain_energy = []

for N in [2, 4, 6]:
    coeffs = compute_coefficients(N)
    W, VM = np.zeros_like(X), np.zeros_like(X)
    for i in range(len(x)):
        for j in range(len(y)):
            w, Mx, My, Mxy = compute_fields(X[i, j], Y[i, j], coeffs)
            W[i, j] = w
            z_top = h / 2
            sigma_xx = -12 * Mx * z_top / h**3
            sigma_yy = -12 * My * z_top / h**3
            sigma_xy = -12 * Mxy * z_top / h**3
            VM[i, j] = np.sqrt(sigma_xx**2 + sigma_yy**2 - sigma_xx * sigma_yy + 3 * sigma_xy**2)

    plt.figure(figsize=(6, 5))
    plt.contourf(X, Y, W, levels=20, cmap='viridis')
    plt.colorbar()
    plt.title(f'Displacement (N={N})')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()
    
    plt.figure(figsize=(6, 5))
    plt.contourf(X, Y, VM, levels=20, cmap='plasma')
    plt.colorbar()
    plt.title(f'Von Mises Stress (N={N})')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.show()

    energy = 0
    for m in range(1, N + 1):
        for n in range(1, N + 1):
            if m % 2 == 0 or n % 2 == 0:
                continue
            denom = ((m * np.pi / a)**2 + (n * np.pi / b)**2)**2
            Wmn = (16 * q0) / (D * np.pi**6 * m * n * denom)
            energy += Wmn**2 * denom
    strain_energy.append(D / 2 * energy * a * b / 4)

plt.figure(figsize=(8, 6))
for N in [2, 4, 6]:
    coeffs = compute_coefficients(N)
    w, Mx, My, Mxy = compute_fields(center_x, center_y, coeffs)
    sigma_xx = -12 * Mx * z / h**3
    sigma_yy = -12 * My * z / h**3
    sigma_xy = -12 * Mxy * z / h**3
    plt.plot(sigma_xx, z, label=f'N={N} σ_xx', linestyle='-')
    plt.plot(sigma_yy, z, label=f'N={N} σ_yy', linestyle='--')
    plt.plot(sigma_xy, z, label=f'N={N} σ_xy', linestyle=':')

plt.legend()
plt.xlabel('Stress (Pa)')
plt.ylabel('z (m)')
plt.title('Through-Thickness Stresses')

plt.figure()
plt.plot([2, 4, 6], strain_energy, marker='o', linestyle='-', color='b')
plt.xlabel('N')
plt.ylabel('Strain Energy (J)')
plt.title('Strain Energy vs N')
plt.grid()

plt.show()