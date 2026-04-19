import numpy as np
from numpy.fft import fft, ifft, fftfreq
from scipy.fft import dst, idst

# --- Natural units: hbar = m = 1 ---
hbar = 1.0
m    = 1.0

# =============================================================================
#  STUDENT IMPLEMENTATION (ASSIGNMENTS 1-3)
# =============================================================================

def gaussian_packet(x, x0, sigma, k0):
    """
    Create a normalized Gaussian wave packet.

        psi(x) = A * exp(-(x - x0)^2 / (4*sigma^2)) * exp(i*k0*x)

    where A is chosen so that the integral of |psi|^2 dx over the grid = 1.

    Args:
        x     : np.ndarray, position grid
        x0    : float, center of the packet
        sigma : float, spatial width
        k0    : float, central wavenumber (momentum in natural units)

    Returns:
        psi : complex np.ndarray, normalized wavefunction
    """
    pass


def momentum_wavefunction(psi, dx):
    """
    Compute the momentum-space wavefunction phi(k) from psi(x) using the FFT.

    The continuous Fourier transform is approximated as:
        phi(k) ≈ FFT(psi) * dx

    Normalize phi so that the integral of |phi|^2 dk = 1.

    Args:
        psi : complex np.ndarray, position-space wavefunction
        dx  : float, grid spacing

    Returns:
        k   : np.ndarray, wavenumber grid (use 2*pi*fftfreq(N, d=dx))
        phi : complex np.ndarray, normalized momentum wavefunction

    Hint: dk = 2*pi / (N*dx)
    """
    pass


def evolve_free_particle(psi, k, dt):
    """
    Evolve a free-particle wavefunction by one time step dt.

    Free particle has Hamiltonian H = hbar^2 k^2 / (2m).
    The exact propagator in momentum space is:
        phi(k, t+dt) = phi(k, t) * exp(-i * hbar * k^2 * dt / (2*m))

    Args:
        psi : complex np.ndarray, current position-space wavefunction
        k   : np.ndarray, wavenumber grid
        dt  : float, time step

    Returns:
        psi_new : complex np.ndarray, evolved wavefunction

    Hint: FFT → apply phase factor → IFFT back.
    """
    pass


# =============================================================================
#  STUDENT IMPLEMENTATION (ASSIGNMENT 4)
# =============================================================================

def well_eigenfunction(x, n, L):
    """
    Compute the n-th eigenfunction of a particle in an infinite square well [0, L].

        psi_n(x) = sqrt(2/L) * sin(n * pi * x / L)

    For integer n >= 1 this satisfies the boundary conditions psi(0) = psi(L) = 0.

    Args:
        x : np.ndarray, position values (should span [0, L])
        n : float, quantum number (try non-integer values to see what breaks!)
        L : float, well length

    Returns:
        psi : real np.ndarray, eigenfunction (normalized if n is a positive integer)

    Experiment: move the slider for n to a non-integer and watch psi(L) ≠ 0.
    """
    pass


# =============================================================================
#  STUDENT IMPLEMENTATION (ASSIGNMENTS 5-6)
# =============================================================================

def split_operator_step(psi, k, V, dt):
    """
    Advance psi by one time step using the split-operator (Trotter) method.

    The Hamiltonian H = T + V is split as:
        exp(-i H dt/hbar) ≈ exp(-i V dt/2hbar) · exp(-i T dt/hbar) · exp(-i V dt/2hbar)

    Steps:
        1. Half potential kick:    psi *= exp(-i V dt / (2*hbar))
        2. FFT to k-space
        3. Kinetic phase factor:   phi *= exp(-i hbar k^2 dt / (2*m))
        4. IFFT back to x-space
        5. Half potential kick again

    Args:
        psi : complex np.ndarray, current wavefunction
        k   : np.ndarray, wavenumber grid (from fftfreq)
        V   : np.ndarray, potential at each grid point
        dt  : float, time step

    Returns:
        psi_new : complex np.ndarray, evolved wavefunction
    """
    pass


def dst_energy_levels(N, L):
    """
    Compute the DST energy eigenvalues for an infinite square well of length L.

    These are the energies associated with the N sine modes used by the
    Discrete Sine Transform propagator.  In natural units (hbar = m = 1):

        E_n = (n * pi / L)^2 / 2,    n = 1, 2, ..., N

    Args:
        N : int, number of grid points
        L : float, well length

    Returns:
        E_k : np.ndarray of length N, energy eigenvalues

    Hint: use np.arange(1, N+1) for the mode indices.
    """
    pass


# =============================================================================
#  STUDENT IMPLEMENTATION (ASSIGNMENT 7)
# =============================================================================

def absorbing_mask(N, gobble_frac=0.1):
    """
    Build an absorbing boundary mask (the "Gobbler") to prevent wrap-around.

    The mask equals 1 in the interior and tapers smoothly to 0 at both edges
    using a sine-squared envelope:

        taper = sin(linspace(0, pi/2, gobble_width))^2

    Args:
        N           : int, number of grid points
        gobble_frac : float, fraction of grid absorbed at each edge (default 0.1)

    Returns:
        mask : real np.ndarray of length N, values in [0, 1]
               mask ≈ 1 in the interior, smoothly → 0 at both ends
    """
    pass


# =============================================================================
#  BONUS (EXTRA CREDIT)
# =============================================================================

def compute_tunneling_probability(V0, k0, a):
    """
    Analytical quantum tunneling transmission probability through a rectangular barrier.

    In natural units (hbar = m = 1), the particle energy is E = k0^2 / 2.

    For E < V0 (classically forbidden region):
        kappa = sqrt(2*(V0 - E))
        T = 1 / (1 + V0^2 * sinh(kappa*a)^2 / (4*E*(V0-E)))

    For E >= V0 (classically allowed): T = 1.0  (simplified assumption).

    Args:
        V0 : float, barrier height
        k0 : float, central wavenumber of the packet
        a  : float, barrier full width

    Returns:
        T : float, transmission probability in [0, 1]
    """
    pass
