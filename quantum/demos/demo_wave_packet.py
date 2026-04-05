"""
Demo: Wave Packet Propagation via Split-Operator Method
========================================================
Animates a Gaussian wave packet:
  1. Bouncing in a mathematically perfect infinite well (DST method)
  2. Tunneling through a rectangular barrier with UI controls and ABCs

Usage:
    python demo_wave_packet.py [--mode well|barrier]

Requires: numpy, scipy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.widgets as widgets
from numpy.fft import fft, ifft, fftfreq
from scipy.fft import dst, idst
import argparse

# --- Natural units for simplicity: hbar = m = 1 ---
hbar = 1.0
m    = 1.0

def gaussian_packet(x, x0, sigma, k0):
    """Create a normalized Gaussian wave packet."""
    psi = (1 / (2 * np.pi * sigma**2))**0.25 * \
          np.exp(-(x - x0)**2 / (4 * sigma**2)) * \
          np.exp(1j * k0 * x)
    dx = x[1] - x[0]
    psi /= np.sqrt(np.sum(np.abs(psi)**2) * dx)
    return psi

def split_operator_step_fft(psi, k, V, dt):
    """Standard FFT split-operator for finite potentials (Barrier)."""
    psi = psi * np.exp(-1j * V * dt / (2 * hbar))
    psi_k = fft(psi)
    psi_k = psi_k * np.exp(-1j * hbar * k**2 * dt / (2 * m))
    psi = ifft(psi_k)
    psi = psi * np.exp(-1j * V * dt / (2 * hbar))
    return psi

def split_operator_step_dst(psi, E_k, dt):
    """DST split-operator for mathematically perfect infinite wells."""
    psi_k = dst(psi, type=1, norm='ortho')
    psi_k = psi_k * np.exp(-1j * E_k * dt / hbar)
    psi = idst(psi_k, type=1, norm='ortho')
    return psi

def run_simulation(mode='well'):
    N = 1024
    
    # --- Shared Animation State ---
    sim_time = 0.0
    norms = []
    
    if mode == 'well':
        L = 20.0
        x = np.linspace(0, L, N, endpoint=False)
        dx = x[1] - x[0]

        # DST Energy Levels (replaces V array for infinite walls)
        n_modes = np.arange(1, N + 1)
        E_k = (hbar**2 * (n_modes * np.pi / L)**2) / (2 * m)

        x0 = L * 0.3
        sigma_init = L * 0.05
        k0_init = 15.0
        dt = 0.005
        n_steps_per_frame = 10
        title = "Wave Packet in a True Infinite Well (V=∞)"

    elif mode == 'barrier':
        L = 40.0
        x = np.linspace(-L/2, L/2, N, endpoint=False)
        dx = x[1] - x[0]

        barrier_width = 0.5
        v0_init = 25.0
        
        # Absorbing Boundary Condition (Gobbler)
        mask = np.ones(N)
        gobble_width = int(0.1 * N)
        taper = np.sin(np.linspace(0, np.pi/2, gobble_width))**2
        mask[:gobble_width] = taper
        mask[-gobble_width:] = taper[::-1]

        k = 2 * np.pi * fftfreq(N, d=dx)

        x0 = -L * 0.25
        sigma_init = 1.0
        k0_init = 10.0
        dt = 0.002
        n_steps_per_frame = 15
        title = "Quantum Tunneling Explorer"

    else:
        raise ValueError(f"Unknown mode: {mode}")

    # --- Setup Plot Layout ---
    # Made bottom margin larger (0.30) to fit 4 sliders
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), gridspec_kw={'height_ratios': [3, 1]})
    fig.subplots_adjust(bottom=0.30) 
    fig.suptitle(title, fontsize=14, fontweight='bold')

    # Initial wave state
    psi = gaussian_packet(x, x0, sigma_init, k0_init)

    # Plot wave lines
    prob_line, = ax1.plot(x, np.abs(psi)**2, color='#2196F3', lw=2, label='|ψ(x)|²')
    real_line, = ax1.plot(x, psi.real, color='#4CAF50', lw=0.8, alpha=0.5, label='Re[ψ(x)]')

    y_max = np.max(np.abs(psi)**2) * 2.5
    ax1.set_ylim(-y_max * 0.3, y_max)
    ax1.set_ylabel('|ψ(x)|²')
    ax1.legend(loc='upper right', fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Define dynamic potential polygon for barrier mode
    v_poly = None

    if mode == 'barrier':
        V = np.where(np.abs(x) < barrier_width / 2, v0_init, 0.0)
        V_display = np.clip(V / max(np.max(np.abs(psi)**2) * 5, 1), 0, y_max * 0.8)
        v_poly = ax1.fill_between(x, 0, V_display, color='#FF5722', alpha=0.3, label='V(x)')
        ax1.set_xlim(x[0], x[-1])
        
        # Visually shade the Absorbing Zones (the Gobbler)
        ax1.axvspan(x[0], x[gobble_width], color='#607D8B', alpha=0.2)
        ax1.axvspan(x[-gobble_width], x[-1], color='#607D8B', alpha=0.2)
        ax1.text(x[gobble_width//2], y_max*0.75, 'Absorbing\nZone', ha='center', va='center', rotation=90, color='#455A64', alpha=0.7)
        ax1.text(x[-gobble_width//2], y_max*0.75, 'Absorbing\nZone', ha='center', va='center', rotation=90, color='#455A64', alpha=0.7)
    
    elif mode == 'well':
        # Draw infinite walls
        ax1.axvspan(x[0] - L*0.1, x[0], color='gray', alpha=0.3, hatch='//')
        ax1.axvspan(x[-1], x[-1] + L*0.1, color='gray', alpha=0.3, hatch='//')
        ax1.axvline(x[0], color='black', lw=4)
        ax1.axvline(x[-1], color='black', lw=4)
        ax1.text(x[0] - L*0.05, y_max*0.8, r'$V=\infty$', fontsize=16, ha='center', weight='bold')
        ax1.text(x[-1] + L*0.05, y_max*0.8, r'$V=\infty$', fontsize=16, ha='center', weight='bold')
        ax1.set_xlim(x[0] - L*0.1, x[-1] + L*0.1)

    # Setup Norm Plot
    norm_line, = ax2.plot([], [], color='#9C27B0', lw=1.5)
    ax2.set_xlim(0, 400)
    ax2.set_ylim(0.0, 1.1)
    ax2.set_xlabel('Calculation Frames')
    ax2.set_ylabel('∫|ψ|² dx')
    ax2.axhline(1.0, color='gray', ls='--', alpha=0.5)
    ax2.grid(True, alpha=0.3)

    time_text = ax1.text(0.02, 0.95, '', transform=ax1.transAxes, fontsize=11,
                         verticalalignment='top', fontfamily='monospace',
                         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # --- Setup UI Controls ---
    # Shifted Y-coordinates up to fit 4 sliders cleanly
    ax_k0    = plt.axes([0.15, 0.20, 0.65, 0.03])
    ax_sigma = plt.axes([0.15, 0.15, 0.65, 0.03])
    
    slider_k0 = widgets.Slider(ax_k0, 'Momentum ($k_0$)', 2.0, 30.0, valinit=k0_init)
    slider_sigma = widgets.Slider(ax_sigma, 'Width ($\sigma$)', 0.5, 3.0, valinit=sigma_init)
    
    slider_v0 = None
    slider_a = None
    if mode == 'barrier':
        ax_v0 = plt.axes([0.15, 0.10, 0.65, 0.03])
        ax_a  = plt.axes([0.15, 0.05, 0.65, 0.03])
        slider_v0 = widgets.Slider(ax_v0, 'Barrier Height ($V_0$)', 0.0, 50.0, valinit=v0_init)
        slider_a  = widgets.Slider(ax_a, 'Barrier Width ($a$)', 0.1, 2.0, valinit=barrier_width)

    # Fire Button placed slightly off to the side
    ax_btn = plt.axes([0.82, 0.10, 0.1, 0.08])
    btn_fire = widgets.Button(ax_btn, 'Fire Packet', color='#81C784', hovercolor='#66BB6A')

    # Update logic when sliders change
    def reset_sim(event=None):
        # Added barrier_width to nonlocal so it persists into the animate loop
        nonlocal psi, sim_time, norms, V, v_poly, barrier_width 
        sim_time = 0.0
        norms.clear()
        
        # Build new initial packet
        psi = gaussian_packet(x, x0, slider_sigma.val, slider_k0.val)
        
        # Update barrier if applicable
        if mode == 'barrier' and slider_v0 and slider_a:
            barrier_width = slider_a.val
            V = np.where(np.abs(x) < barrier_width / 2, slider_v0.val, 0.0)
            V_display = np.clip(V / max(np.max(np.abs(psi)**2) * 5, 1), 0, y_max * 0.8)
            v_poly.remove()
            v_poly = ax1.fill_between(x, 0, V_display, color='#FF5722', alpha=0.3)

    # Attach UI triggers
    slider_k0.on_changed(reset_sim)
    slider_sigma.on_changed(reset_sim)
    btn_fire.on_clicked(reset_sim)
    if slider_v0:
        slider_v0.on_changed(reset_sim)
        slider_a.on_changed(reset_sim)

    # --- Animation Loop ---
    def animate(frame):
        nonlocal psi, sim_time
        
        for _ in range(n_steps_per_frame):
            if mode == 'well':
                psi = split_operator_step_dst(psi, E_k, dt)
            elif mode == 'barrier':
                psi = split_operator_step_fft(psi, k, V, dt)
                psi *= mask
        
        sim_time += n_steps_per_frame * dt

        prob = np.abs(psi)**2
        prob_line.set_ydata(prob)
        real_line.set_ydata(psi.real)

        # Update Norm Tracking
        norm = np.sum(prob) * dx
        norms.append(norm)
        if len(norms) > ax2.get_xlim()[1]:
            ax2.set_xlim(0, len(norms) + 100)
            
        norm_line.set_data(range(len(norms)), norms)

        # Calculate dynamic Transmission and Reflection proportions
        if mode == 'barrier':
            # Dynamically uses the updated barrier_width
            R = np.sum(prob[x < -barrier_width/2]) * dx
            T = np.sum(prob[x > barrier_width/2]) * dx
            status = f't = {sim_time:.2f} | R: {R*100:04.1f}% | T: {T*100:04.1f}%'
        else:
            status = f't = {sim_time:.3f} | norm = {norm:.6f}'

        time_text.set_text(status)
        return prob_line, real_line, norm_line, time_text

    anim = animation.FuncAnimation(fig, animate, interval=30, blit=False, cache_frame_data=False)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Wave packet propagation demo')
    parser.add_argument('--mode', choices=['well', 'barrier'], default='barrier',
                        help='Simulation mode: "well" (infinite well) or "barrier" (tunneling)')
    args = parser.parse_args()
    run_simulation(args.mode)
