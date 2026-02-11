import numpy as np

# =============================================================================
#  STUDENT IMPLEMENTATION (ASSIGNMENTS 1-4)
# =============================================================================

CONST_EPSILON = 0#1e-5

def z_cross(v,u):
    return v[0] * u[1] - u[0] * v[1]


def intersect_line_infinite(ray_origin, ray_dir, p1, p2):
    v = p2 - p1
    denom = z_cross(v, ray_dir)
    if abs(denom) <= CONST_EPSILON: return float('inf'), float('inf')
    u = z_cross(ray_origin - p1, ray_dir) / denom
    t = z_cross(ray_origin - p1, v) / denom
    return t, u

#def intersect_line_infinite(ray_origin, ray_dir, p1, p2):
#    pass
    

def intersect_segment(ray_origin, ray_dir, p1, p2):
    # LEVEL 2: Implement boundary check
    t, u = intersect_line_infinite(ray_origin, ray_dir, p1, p2)
    if 0 <= u <= 1 and t > CONST_EPSILON: return t
    return float('inf')

#def intersect_segment(ray_origin, ray_dir, p1, p2):
#    pass


def intersect_circle_infinite(ray_origin, ray_dir, center, radius):
    # LEVEL 3: Standard quadratic intersection
    OC = ray_origin - center
    b = 2 * np.dot(ray_dir, OC)
    c = np.dot(OC, OC) - radius**2
    discriminant = b**2 - 4*c
    if discriminant < 0: return []
    sqrt_d = np.sqrt(discriminant)
    return [(-b - sqrt_d) / 2, (-b + sqrt_d) / 2]

#def intersect_circle_infinite(ray_origin, ray_dir, center, radius):
#    pass


def intersect_arc(ray_origin, ray_dir, center, radius, axis, cos_half_angle):
    # LEVEL 4: Angular sector check
    ts = intersect_circle_infinite(ray_origin, ray_dir, center, radius)
    best_t = float('inf')
    for t in ts:
        if t > CONST_EPSILON:
            P = ray_origin + t * ray_dir
            vec_CP = (P - center) / radius
            if np.dot(vec_CP, axis) >= cos_half_angle - CONST_EPSILON:
                if t < best_t: best_t = t
    return best_t

#def intersect_arc(ray_origin, ray_dir, center, radius, axis, cos_half_angle):
#    pass

# =============================================================================
#  INTERSECTION DISPATCHER, IMPLEMENTATION PROVIDED
# =============================================================================
def intersect_curve(ray_origin, ray_dir, curve):
    if curve['type'] == 'line':
        return intersect_segment(ray_origin, ray_dir, curve['p1'], curve['p2'])
    elif curve['type'] == 'arc':
        return intersect_arc(ray_origin, ray_dir, curve['center'], 
                             curve['radius'], curve['axis'], curve['cos_half_angle'])
    return float('inf')
# =============================================================================


# =============================================================================
#  STUDENT IMPLEMENTATION (ASSIGNMENTS 5-6)
# =============================================================================


def calculate_normal_segment(ray_dir, p1, p2):
    # LEVEL 5: Line Normal + Flip Logic
    tangent = p2 - p1
    normal = np.array([-tangent[1], tangent[0]], dtype=float)
    normal /= np.linalg.norm(normal)
    if np.dot(ray_dir, normal) > 0: normal = -normal
    return normal

#def calculate_normal_segment(ray_dir, p1, p2):
#    pass


def calculate_normal_arc(hit_point, ray_dir, center):
    # LEVEL 6: Arc Normal + Flip Logic
    normal = (hit_point - center).astype(float)
    normal /= np.linalg.norm(normal)
    if np.dot(ray_dir, normal) > 0: normal = -normal
    return normal

#def calculate_normal_arc(hit_point, ray_dir, center):
#    pass

# =============================================================================
#  NORMAL DISPATCHER, IMPLEMENTATION PROVIDED
# =============================================================================
def calculate_normal(hit_point, ray_dir, curve):
    if curve['type'] == 'line':
        return calculate_normal_segment(ray_dir, curve['p1'], curve['p2'])
    elif curve['type'] == 'arc':
        return calculate_normal_arc(hit_point, ray_dir, curve['center'])
    return None
# =============================================================================


# =============================================================================
#  STUDENT IMPLEMENTATION (ASSIGNMENT 7)
# =============================================================================


def refract_vector(ray_dir, normal, n1, n2):
    # LEVEL 7: Snell's Law
    eta = n1 / n2
    cos_theta1 = -np.dot(ray_dir, normal)
    sin2_theta1 = 1 - cos_theta1**2
    sin2_theta2 = eta**2 * sin2_theta1
    if sin2_theta2 > 1.0: return None 
    cos_theta2 = np.sqrt(1 - sin2_theta2)
    return eta * ray_dir + (eta * cos_theta1 - cos_theta2) * normal

#def refract_vector(ray_dir, normal, n1, n2):
#    pass

# =============================================================================



# =============================================================================
#  RAY TRACING, IMPLEMENTATION PROVIDED
# =============================================================================

def trace_ray_step(ray_origin, ray_dir, current_medium_id, curves, media):
    # --- 1. GEOMETRY: FIND NEAREST INTERSECTION ---
    # Create list of (t, curve) tuples
    hits = [(intersect_curve(ray_origin, ray_dir, c), c) for c in curves]
    
    # Filter for valid hits and find the one with the minimum t
    min_t, hit_curve = min(hits, key=lambda x: x[0], default=(float('inf'), None))

    if hit_curve is None or min_t == float('inf'):
        return None, None, None 

    best_hit = ray_origin + min_t * ray_dir

    # --- 2. PHYSICS & 3. REFRACTION ---
    # (The rest of your logic remains the same)
    n1 = media[current_medium_id]['n']
    new_medium_id = 0 if hit_curve['medium_id'] == current_medium_id else hit_curve['medium_id']
    n2 = media[new_medium_id]['n']

    normal = calculate_normal(best_hit, ray_dir, hit_curve)
    new_dir = refract_vector(ray_dir, normal, n1, n2)

    return best_hit, new_dir, new_medium_id
    
# =============================================================================


# =============================================================================
#  STUDENT IMPLEMENTATION (ADDITIONAL ASSIGNMENT)
# =============================================================================


def calculate_focal_length(R1, R2, d, n):
    power = (n - 1) * (1/R1 - 1/R2 + ((n - 1) * d * 1/R1 * 1/R2) / n)
    return 1.0 / power
    

def calculate_slab_displacement(d, n, theta):
    theta_prime = np.arcsin(np.sin(theta) / n)
    h = d * np.sin(theta - theta_prime) / np.cos(theta_prime)
    return h
    
    
