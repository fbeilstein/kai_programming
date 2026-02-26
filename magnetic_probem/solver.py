import numpy as np
import cupy as cp

def solve_magnetic_field(X, mu_grid, M_grid, max_iter=5000, tol=1e-5):
    print(f"\n--- SOLVER INITIALIZED ---")
    print(f"Hardware Acceleration: CUPY (GTX 950M)")
    
    # 1. Grid Spacing
    dx = float(X[1, 0, 0] - X[0, 0, 0])
    dy, dz = dx, dx 

    # 2. Upload to GPU
    mu = cp.array(mu_grid, dtype=cp.float32)
    M = cp.array(M_grid, dtype=cp.float32)
    
    # 3. Compute Divergence of M (RHS of our equation)
    div_M = cp.zeros_like(mu)
    div_M[1:-1, 1:-1, 1:-1] = (
        (M[2:, 1:-1, 1:-1, 0] - M[:-2, 1:-1, 1:-1, 0]) / (2 * dx) +
        (M[1:-1, 2:, 1:-1, 1] - M[1:-1, :-2, 1:-1, 1]) / (2 * dy) +
        (M[1:-1, 1:-1, 2:, 2] - M[1:-1, 1:-1, :-2, 2]) / (2 * dz)
    )

    # 4. Define Coefficients (Used inside the loop)
    cx_p = 0.5 * (mu[2:, 1:-1, 1:-1] + mu[1:-1, 1:-1, 1:-1]) / dx**2
    cx_m = 0.5 * (mu[:-2, 1:-1, 1:-1] + mu[1:-1, 1:-1, 1:-1]) / dx**2
    cy_p = 0.5 * (mu[1:-1, 2:, 1:-1] + mu[1:-1, 1:-1, 1:-1]) / dy**2
    cy_m = 0.5 * (mu[1:-1, :-2, 1:-1] + mu[1:-1, 1:-1, 1:-1]) / dy**2
    cz_p = 0.5 * (mu[1:-1, 1:-1, 2:] + mu[1:-1, 1:-1, 1:-1]) / dz**2
    cz_m = 0.5 * (mu[1:-1, 1:-1, :-2] + mu[1:-1, 1:-1, 1:-1]) / dz**2

    c_sum = cx_p + cx_m + cy_p + cy_m + cz_p + cz_m
    rho_interior = div_M[1:-1, 1:-1, 1:-1]

    # 5. Initialize Potential
    psi = cp.zeros_like(mu)
    
    print("Running Jacobi Relaxation...")

    # 6. Relaxation Loop
    for i in range(max_iter):
        p_old = psi.copy()
        
        # Jacobi Update
        psi[1:-1, 1:-1, 1:-1] = (
            cx_p * p_old[2:, 1:-1, 1:-1] + cx_m * p_old[:-2, 1:-1, 1:-1] +
            cy_p * p_old[1:-1, 2:, 1:-1] + cy_m * p_old[1:-1, :-2, 1:-1] +
            cz_p * p_old[1:-1, 1:-1, 2:] + cz_m * p_old[1:-1, 1:-1, :-2] - rho_interior
        ) / c_sum

        # Check convergence every 100 iterations
        if i % 100 == 0:
            diff = cp.max(cp.abs(psi[1:-1, 1:-1, 1:-1] - p_old[1:-1, 1:-1, 1:-1]))
            if i > 0 and diff < tol: 
                print(f"Converged at iter {i} (Err: {diff:.6e})")
                break
                
    if i == max_iter - 1:
        print("Reached max iterations without full convergence.")

    # 7. Calculate B-Field = mu_r * (-grad psi) + M
    Hx = cp.zeros_like(mu)
    Hy = cp.zeros_like(mu)
    Hz = cp.zeros_like(mu)

    Hx[1:-1, :, :] = -(psi[2:, :, :] - psi[:-2, :, :]) / (2 * dx)
    Hy[:, 1:-1, :] = -(psi[:, 2:, :] - psi[:, :-2, :]) / (2 * dy)
    Hz[:, :, 1:-1] = -(psi[:, :, 2:] - psi[:, :, :-2]) / (2 * dz)

    Bx = mu * Hx + M[..., 0]
    By = mu * Hy + M[..., 1]
    Bz = mu * Hz + M[..., 2]
    
    # 8. Stack and Return to CPU
    B_grid = cp.stack((Bx, By, Bz), axis=-1)
    return cp.asnumpy(B_grid)
