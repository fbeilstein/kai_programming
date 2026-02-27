import cupy as cp
import numpy as np

def solve_magnetic_field(X, mu_grid, M_grid, max_iter=3000, tol=1e-6):
    print("\n--- SOLVER: Harmonic Mean & Neumann BCs ---")

    dx = float(X[1,0,0] - X[0,0,0])
    dy = dz = dx

    mu = cp.array(mu_grid, dtype=cp.float32)
    M = cp.array(M_grid, dtype=cp.float32)

    # --- Compute divergence of M ---
    div_M = cp.zeros_like(mu)
    div_M[1:,:,:] += (M[1:,:,:,0] - M[:-1,:,:,0]) / dx
    div_M[:,1:,:] += (M[:,1:,:,1] - M[:,:-1,:,1]) / dy
    div_M[:,:,1:] += (M[:,:,1:,2] - M[:,:,:-1,2]) / dz
    
    # Zero out global charge to prevent runaway monopole lines
    div_M -= cp.mean(div_M)

    psi = cp.zeros_like(mu)

    # --- HARMONIC averaging for μ at faces ---
    # This keeps the magnet boundaries mathematically sharp
    def harmonic_mean(m1, m2):
        return (2.0 * m1 * m2) / (m1 + m2)

    mu_xp = harmonic_mean(mu[1:-1,1:-1,1:-1], mu[2:,1:-1,1:-1])
    mu_xm = harmonic_mean(mu[1:-1,1:-1,1:-1], mu[:-2,1:-1,1:-1])
    mu_yp = harmonic_mean(mu[1:-1,1:-1,1:-1], mu[1:-1,2:,1:-1])
    mu_ym = harmonic_mean(mu[1:-1,1:-1,1:-1], mu[1:-1,:-2,1:-1])
    mu_zp = harmonic_mean(mu[1:-1,1:-1,1:-1], mu[1:-1,1:-1,2:])
    mu_zm = harmonic_mean(mu[1:-1,1:-1,1:-1], mu[1:-1,1:-1,:-2])

    B_denom = mu_xp + mu_xm + mu_yp + mu_ym + mu_zp + mu_zm

    for i in range(max_iter):
        old = psi.copy()

        A = (
            mu_xp * psi[2:,1:-1,1:-1] +
            mu_xm * psi[:-2,1:-1,1:-1] +
            mu_yp * psi[1:-1,2:,1:-1] +
            mu_ym * psi[1:-1,:-2,1:-1] +
            mu_zp * psi[1:-1,1:-1,2:] +
            mu_zm * psi[1:-1,1:-1,:-2]
        )

        psi[1:-1,1:-1,1:-1] = (A - div_M[1:-1,1:-1,1:-1] * dx**2) / B_denom

        # --- NEUMANN BOUNDARY CONDITIONS ---
        # Allow the field lines to "breathe" beyond the grid edges
        psi[0, :, :] = psi[1, :, :]
        psi[-1, :, :] = psi[-2, :, :]
        psi[:, 0, :] = psi[:, 1, :]
        psi[:, -1, :] = psi[:, -2, :]
        psi[:, :, 0] = psi[:, :, 1]
        psi[:, :, -1] = psi[:, :, -2]

        if i % 100 == 0:
            err = cp.max(cp.abs(psi[1:-1,1:-1,1:-1] - old[1:-1,1:-1,1:-1]))
            if i > 0 and err < tol:
                print(f"Converged at {i}, err={err:.2e}")
                break


    Hx = cp.zeros_like(mu)
    Hy = cp.zeros_like(mu)
    Hz = cp.zeros_like(mu)

    # This ensures H aligns perfectly with the cell-centered 'mu'
    Hx[1:-1, :, :] = -(psi[2:, :, :] - psi[:-2, :, :]) / (2 * dx)
    Hy[:, 1:-1, :] = -(psi[:, 2:, :] - psi[:, :-2, :]) / (2 * dy)
    Hz[:, :, 1:-1] = -(psi[:, :, 2:] - psi[:, :, :-2]) / (2 * dz)

    # --- B = μH + M ---
    Bx = mu * Hx + M[..., 0]
    By = mu * Hy + M[..., 1]
    Bz = mu * Hz + M[..., 2]

    return cp.asnumpy(cp.stack((Bx, By, Bz), axis=-1))

