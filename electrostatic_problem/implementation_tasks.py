import numpy as np
from numba import njit


# =========================================================
# PHYSICAL & ALGORITHMIC CONSTANTS
# =========================================================

# The squared radius of the charge "hitbox" (0.3^2).
# Used to stop lines from piercing the charge and to prevent divide-by-zero singularities.
CORE_RADIUS_SQ = 0.09 

# Minimum allowed electric field magnitude to prevent division by zero during normalization.
MIN_E_MAG = 1e-9      

# =========================================================
# PART 1: THE STUDENT ASSIGNMENT
# =========================================================

@njit(cache=True)
def calc_E_dir(pt, charges_positions, charges_values):
    """Calculate the normalized E-field direction, use NumPy broadcasting where possible."""
    r = pt - charges_positions 
    r2 = np.sum(r**2, axis=1)
    
    # Prevent singularities using the physical core boundary
    r2[r2 < CORE_RADIUS_SQ] = CORE_RADIUS_SQ
        
    r_mag = np.sqrt(r2)
    factor = charges_values / (r_mag * r2)
    
    E = np.sum(r * factor.reshape(-1, 1), axis=0)
    E_mag = np.linalg.norm(E)
    
    # Prevent division by zero if the field is perfectly canceled out
    if E_mag < MIN_E_MAG: 
        E_mag = MIN_E_MAG
    
    return E / E_mag

@njit(cache=True)
def propagate_trajectory(trajectory, current_charge_sign, charges_positions, charges_values, max_steps, step_size):
    """
    STUDENT TODO: 
    Implement the RK2 (Midpoint) numerical integration loop to trace a single 
    electric field line. Modify the 'trajectory' array in-place.
    Stop early if the line hits a charge and return the number of steps taken.
    """
    current_point = trajectory[0]
    
    for step in range(1, max_steps):
        # 1. E-field at current point
        E1_dir = calc_E_dir(current_point, charges_positions, charges_values)
        
        # 2. Scout the midpoint
        mid_pt = current_point + E1_dir * step_size * 0.5 * current_charge_sign
        
        # 3. E-field at midpoint
        E2_dir = calc_E_dir(mid_pt, charges_positions, charges_values)
        
        # 4. Take final step
        current_point = current_point + E2_dir * step_size * current_charge_sign
        
        trajectory[step] = current_point
        
        # 5. Vectorized Collision Check: stop if we hit the core of any charge
        dist_sq = np.sum((current_point - charges_positions)**2, axis=1)
        if np.any(dist_sq < CORE_RADIUS_SQ):  
            return step + 1 
            
    return max_steps
