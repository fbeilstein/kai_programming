"""
Demo: Stationary States of Infinite and Finite Square Wells
============================================================
Visualizes the first N eigenstates and probability densities
for both infinite and finite 1D potential wells.

Usage:
    python demo_well_wavefunctions.py

Requires: numpy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt

# --- Physical constants (atomic units for clarity) ---
hbar = 1.0545718e-34   # J·s
m_e  = 9.1093837e-31   # kg (electron mass)
eV   = 1.602176634e-19  # J per eV

# --- Parameters ---
L     = 1e-9      # well width: 1 nm
N_pts = 500       # grid points
N_states = 5      # number of eigenstates to show


def infinite_well_analytical(x, n, L):
    """Analytical wavefunctions for infinite square well."""
    return np.sqrt(2 / L) * np.sin(n * np.pi * x / L)


def energy_analytical(n, L, m=m_e):
    """Analytical energy levels for infinite square well."""
    return (n**2 * np.pi**2 * hbar**2) / (2 * m * L**2)


def solve_well_numerically(x, V, m=m_e):
    """
    Solve time-independent Schrödinger equation via matrix diagonalization.
    Returns eigenvalues (energies) and eigenvectors (wavefunctions).
    """
    dx = x[1] - x[0]
    N = len(x)

    # Kinetic energy matrix (finite difference Laplacian)
    diag = np.ones(N) * (-2.0)
    off  = np.ones(N - 1)
    T = (-hbar**2 / (2 * m * dx**2)) * (
        np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
    )

    # Hamiltonian = T + V
    H = T + np.diag(V)

    # Solve eigenvalue problem
    energies, states = np.linalg.eigh(H)
    return energies, states


def plot_wavefunctions(ax, x, energies, states, n_show, title, E_scale=1.0):
    """Plot wavefunctions offset by their energy levels."""
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, n_show))
    x_nm = x * 1e9  # convert to nm

    for i in range(n_show):
        E_eV = energies[i] / eV
        psi = states[:, i]
        # Normalize for display
        psi_display = psi / np.max(np.abs(psi)) * E_scale * 0.3
        ax.fill_between(x_nm, E_eV, E_eV + psi_display, alpha=0.3, color=colors[i])
        ax.plot(x_nm, E_eV + psi_display, color=colors[i], lw=1.5,
                label=f'n={i+1}, E={E_eV:.3f} eV')
        ax.axhline(E_eV, color=colors[i], ls='--', alpha=0.4, lw=0.8)

    ax.set_xlabel('x (nm)')
    ax.set_ylabel('Energy (eV)')
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.legend(fontsize=8, loc='upper right')


# === Main ===
fig, axes = plt.subplots(1, 3, figsize=(18, 7))

# --- 1. Infinite well: analytical ---
x_inf = np.linspace(0, L, N_pts)
E_analytical = np.array([energy_analytical(n, L) for n in range(1, N_states + 1)])
psi_analytical = np.array([infinite_well_analytical(x_inf, n, L) for n in range(1, N_states + 1)])

# Create a fake "states" array for plotting
states_ana = psi_analytical.T
plot_wavefunctions(axes[0], x_inf, E_analytical, states_ana, N_states,
                   'Infinite Well (Analytical)', E_scale=E_analytical[-1])

# --- 2. Infinite well: numerical ---
# Use hard wall BCs by making potential huge outside
x_num = np.linspace(-0.2 * L, 1.2 * L, N_pts)
V_inf = np.where((x_num >= 0) & (x_num <= L), 0.0, 1e6 * eV)
E_num, states_num = solve_well_numerically(x_num, V_inf)
plot_wavefunctions(axes[1], x_num, E_num, states_num, N_states,
                   'Infinite Well (Numerical)', E_scale=E_num[N_states - 1])

# --- 3. Finite well ---
V0 = 5.0 * eV  # barrier height: 5 eV
a = L / 2       # half-width
x_fin = np.linspace(-1.5 * L, 1.5 * L, N_pts)
V_fin = np.where(np.abs(x_fin - L/2) <= a, 0.0, V0)
E_fin, states_fin = solve_well_numerically(x_fin, V_fin)

# Only show bound states (E < V0)
n_bound = min(N_states, np.sum(E_fin < V0))
plot_wavefunctions(axes[2], x_fin, E_fin, states_fin, n_bound,
                   f'Finite Well (V₀ = {V0/eV:.1f} eV)', E_scale=V0 / eV)
axes[2].axhline(V0 / eV, color='red', ls='-', lw=2, alpha=0.5, label=f'V₀ = {V0/eV:.0f} eV')
axes[2].legend(fontsize=8, loc='upper right')

plt.tight_layout()
plt.suptitle('Quantum Particle in a Box — Stationary States', fontsize=15, fontweight='bold', y=1.02)
plt.savefig('demo_well_wavefunctions.png', dpi=150, bbox_inches='tight')
plt.show()
print("✓ Saved: demo_well_wavefunctions.png")
