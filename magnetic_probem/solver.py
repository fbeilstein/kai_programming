import numpy as np
import torch

def solve_magnetic_field(X, mu_grid, M_grid, max_iter=5000, tol=1e-5):
    """
    Solves for the magnetic scalar potential using a PyTorch vectorized Jacobi iteration.
    Returns the B-field vector grid (Bx, By, Bz).
    """
    # 1. Hardware Engineering: Auto-detect CUDA
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"\n--- SOLVER INITIALIZED ---")
    print(f"Hardware Acceleration: {device.type.upper()}")
    
    # Extract uniform grid spacing
    dx = X[1, 0, 0] - X[0, 0, 0]
    dy = dx # Assuming uniform grid based on our linspace
    dz = dx 

    # 2. Upload arrays to GPU (or CPU fallback)
    mu = torch.tensor(mu_grid, dtype=torch.float32, device=device)
    M = torch.tensor(M_grid, dtype=torch.float32, device=device)
    
    # 3. Compute Divergence of M (The "Magnetic Charge Density")
    # Using central differences for the interior of the grid
    div_M = torch.zeros_like(mu)
    div_M[1:-1, :, :] += (M[2:, :, :, 0] - M[:-2, :, :, 0]) / (2 * dx)
    div_M[:, 1:-1, :] += (M[:, 2:, :, 1] - M[:, :-2, :, 1]) / (2 * dy)
    div_M[:, :, 1:-1] += (M[:, :, 2:, 2] - M[:, :, :-2, 2]) / (2 * dz)

    # 4. Setup Averaged Permeability Coefficients
    # We average mu between adjacent nodes to handle the boundaries between air and iron smoothly
    cx_p = 0.5 * (mu[2:, 1:-1, 1:-1] + mu[1:-1, 1:-1, 1:-1]) / dx**2
    cx_m = 0.5 * (mu[:-2, 1:-1, 1:-1] + mu[1:-1, 1:-1, 1:-1]) / dx**2
    cy_p = 0.5 * (mu[1:-1, 2:, 1:-1] + mu[1:-1, 1:-1, 1:-1]) / dy**2
    cy_m = 0.5 * (mu[1:-1, :-2, 1:-1] + mu[1:-1, 1:-1, 1:-1]) / dy**2
    cz_p = 0.5 * (mu[1:-1, 1:-1, 2:] + mu[1:-1, 1:-1, 1:-1]) / dz**2
    cz_m = 0.5 * (mu[1:-1, 1:-1, :-2] + mu[1:-1, 1:-1, 1:-1]) / dz**2

    c_sum = cx_p + cx_m + cy_p + cy_m + cz_p + cz_m
    rho_interior = div_M[1:-1, 1:-1, 1:-1]

    # Initialize potential grid
    psi = torch.zeros_like(mu)

    print("Running Jacobi Relaxation Loop...")
    
    # 5. The GPU Relaxation Loop
    for i in range(max_iter):
        psi_interior = psi[1:-1, 1:-1, 1:-1]
        
        # Grab neighbor potentials
        p_xp = psi[2:, 1:-1, 1:-1]
        p_xm = psi[:-2, 1:-1, 1:-1]
        p_yp = psi[1:-1, 2:, 1:-1]
        p_ym = psi[1:-1, :-2, 1:-1]
        p_zp = psi[1:-1, 1:-1, 2:]
        p_zm = psi[1:-1, 1:-1, :-2]

        # Calculate the new potential based on neighbors and magnetic charge
        psi_new = (
            cx_p * p_xp + cx_m * p_xm +
            cy_p * p_yp + cy_m * p_ym +
            cz_p * p_zp + cz_m * p_zm - rho_interior
        ) / c_sum

        # Check for convergence every 100 iterations to save processing time
        if i % 100 == 0:
            diff = torch.max(torch.abs(psi_new - psi_interior))
            if diff < tol:
                print(f"Converged at iteration {i} (Error: {diff:.6f})")
                break
        
        # Apply the update to the interior
        psi[1:-1, 1:-1, 1:-1] = psi_new

    if i == max_iter - 1:
        print(f"Warning: Reached max iterations ({max_iter}) before strict convergence.")

    # 6. Calculate H-Field = -gradient(psi)
    Hx = torch.zeros_like(mu)
    Hy = torch.zeros_like(mu)
    Hz = torch.zeros_like(mu)

    Hx[1:-1, :, :] = -(psi[2:, :, :] - psi[:-2, :, :]) / (2 * dx)
    Hy[:, 1:-1, :] = -(psi[:, 2:, :] - psi[:, :-2, :]) / (2 * dy)
    Hz[:, :, 1:-1] = -(psi[:, :, 2:] - psi[:, :, :-2]) / (2 * dz)

    # 7. Calculate final B-Field = mu_r * H + M
    Bx = mu * Hx + M[..., 0]
    By = mu * Hy + M[..., 1]
    Bz = mu * Hz + M[..., 2]
    
    # Stack vectors back together and download to CPU RAM
    B_grid = torch.stack((Bx, By, Bz), dim=-1)
    
    return B_grid.cpu().numpy()
