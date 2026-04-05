import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

# 1. Define the spatial grid
N = 1000                 # Number of grid points
L_domain = 10.0          # Total size of the simulation space
x = np.linspace(-L_domain/2, L_domain/2, N)
dx = x[1] - x[0]

# 2. Define the Finite Potential Well V(x)
well_width = 4.0
V0 = 15.0                # Depth of the well (Finite wall height)
V = np.where(np.abs(x) <= well_width/2, 0, V0)

# 3. Build the Kinetic Energy Matrix (K) using Finite Differences
# Formula: d^2/dx^2 f(x) ≈ (f(x+dx) - 2f(x) + f(x-dx)) / dx^2
# We use natural units: hbar^2 / 2m = 1
main_diag = 2.0 / dx**2 * np.ones(N)
off_diag = -1.0 / dx**2 * np.ones(N-1)
K = diags([off_diag, main_diag, off_diag], [-1, 0, 1])

# 4. Build the total Hamiltonian matrix
# H = K + V (V is added to the main diagonal)
H = K + diags(V)

# 5. Solve the eigenvalue problem
# We ask for the 'k' lowest energy states ('SM' = Smallest Magnitude)
num_states = 10
evals, evecs = eigsh(H, k=num_states, which='SM')

# --- Plotting the Results ---
plt.figure(figsize=(10, 6))
plt.plot(x, V, label="Potential V(x)", color="black", linewidth=2)

for i in range(num_states):
    #if evals[i] < V0: # Only plot bound states!
    # Scale and shift the eigenvector to float at its energy level
    wave_amplitude = 2.0  # Just for visual scaling
    psi = evecs[:, i] * wave_amplitude + evals[i]
    
    plt.plot(x, psi, label=f"E_{i} = {evals[i]:.2f}")
    plt.axhline(evals[i], linestyle="--", alpha=0.5)

plt.ylim(-1, V0 + 5)
plt.xlabel("Position (x)")
plt.ylabel("Energy")
plt.legend()
plt.title("Numerical Solution of the Finite Potential Well")
plt.show()
