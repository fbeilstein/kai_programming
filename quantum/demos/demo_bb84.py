"""
Demo: BB84 Quantum Key Distribution Protocol
==============================================
Simulates the BB84 protocol showing:
  - Key generation with random bits and bases
  - Key sifting (matching bases)
  - Effect of an eavesdropper on error rate
  - Side-by-side histogram comparison

Usage:
    python demo_bb84.py

Requires: numpy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt


def simulate_bb84(n_bits=1000, eve_present=False, eve_intercept_prob=1.0):
    """
    Simulate one round of BB84 protocol.

    Parameters:
        n_bits: number of qubits Alice sends
        eve_present: whether Eve intercepts
        eve_intercept_prob: fraction of qubits Eve intercepts

    Returns:
        dict with results
    """
    # Alice: random bits and random bases
    alice_bits  = np.random.randint(0, 2, n_bits)
    alice_bases = np.random.randint(0, 2, n_bits)  # 0 = rectilinear (+), 1 = diagonal (×)

    # Bob: random measurement bases
    bob_bases = np.random.randint(0, 2, n_bits)

    # Transmitted bits (what Bob measures)
    bob_bits = np.copy(alice_bits)

    eve_bits = None
    if eve_present:
        eve_bases = np.random.randint(0, 2, n_bits)
        eve_bits = np.zeros(n_bits, dtype=int)
        intercepts = np.random.random(n_bits) < eve_intercept_prob

        for i in range(n_bits):
            if intercepts[i]:
                # Eve measures
                if eve_bases[i] == alice_bases[i]:
                    # Correct basis: Eve gets the right bit
                    eve_bits[i] = alice_bits[i]
                else:
                    # Wrong basis: Eve gets random result
                    eve_bits[i] = np.random.randint(0, 2)
                # Eve re-sends in her measured basis
                # If Eve's basis differs from Alice's, the state is disturbed
                if eve_bases[i] != alice_bases[i]:
                    # Bob receives a potentially disturbed qubit
                    bob_bits[i] = np.random.randint(0, 2)
            # else: Eve doesn't intercept, qubit passes through unchanged

    # For qubits where Bob uses wrong basis: random result
    wrong_basis = bob_bases != alice_bases
    bob_bits[wrong_basis] = np.random.randint(0, 2, np.sum(wrong_basis))

    # Key sifting: keep only matching bases
    matching = alice_bases == bob_bases
    sifted_alice = alice_bits[matching]
    sifted_bob   = bob_bits[matching]

    # Error rate in sifted key
    errors = sifted_alice != sifted_bob
    error_rate = np.mean(errors) if len(errors) > 0 else 0.0

    return {
        'n_bits': n_bits,
        'n_sifted': len(sifted_alice),
        'error_rate': error_rate,
        'sifted_alice': sifted_alice,
        'sifted_bob': sifted_bob,
        'eve_present': eve_present,
    }


def run_demo():
    """Run multiple BB84 simulations and display results."""
    n_trials = 200
    n_bits = 1000

    # Collect error rates
    errors_no_eve = []
    errors_with_eve = []

    for _ in range(n_trials):
        result = simulate_bb84(n_bits, eve_present=False)
        errors_no_eve.append(result['error_rate'])

        result = simulate_bb84(n_bits, eve_present=True)
        errors_with_eve.append(result['error_rate'])

    errors_no_eve = np.array(errors_no_eve)
    errors_with_eve = np.array(errors_with_eve)

    # --- Plotting ---
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle('BB84 Quantum Key Distribution — Eavesdropper Detection',
                 fontsize=15, fontweight='bold')

    # 1. Histogram comparison
    ax = axes[0]
    bins = np.linspace(0, 0.4, 40)
    ax.hist(errors_no_eve, bins=bins, alpha=0.7, color='#4CAF50', label='No Eve', edgecolor='white')
    ax.hist(errors_with_eve, bins=bins, alpha=0.7, color='#F44336', label='Eve present', edgecolor='white')
    ax.axvline(0.11, color='black', ls='--', lw=2, label='Security threshold (11%)')
    ax.set_xlabel('Error Rate in Sifted Key')
    ax.set_ylabel('Count')
    ax.set_title('Error Rate Distribution')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 2. Single round visualization
    ax = axes[1]
    result_safe = simulate_bb84(100, eve_present=False)
    result_eve  = simulate_bb84(100, eve_present=True)

    n_show = min(30, result_safe['n_sifted'], result_eve['n_sifted'])
    x = np.arange(n_show)

    safe_match = result_safe['sifted_alice'][:n_show] == result_safe['sifted_bob'][:n_show]
    eve_match  = result_eve['sifted_alice'][:n_show] == result_eve['sifted_bob'][:n_show]

    colors_safe = ['#4CAF50' if m else '#F44336' for m in safe_match]
    colors_eve  = ['#4CAF50' if m else '#F44336' for m in eve_match]

    ax.bar(x - 0.2, result_safe['sifted_alice'][:n_show], 0.35, color=colors_safe, alpha=0.8, label='No Eve')
    ax.bar(x + 0.2, result_eve['sifted_alice'][:n_show], 0.35, color=colors_eve, alpha=0.5, label='With Eve')
    ax.set_xlabel('Sifted Key Bit Index')
    ax.set_ylabel('Bit Value')
    ax.set_title('Key Comparison (green=match, red=error)')
    ax.set_yticks([0, 1])
    ax.grid(True, alpha=0.3, axis='y')

    # 3. Eve intercept probability vs error rate
    ax = axes[2]
    intercept_probs = np.linspace(0, 1, 20)
    mean_errors = []
    for p in intercept_probs:
        errs = []
        for _ in range(100):
            r = simulate_bb84(n_bits, eve_present=True, eve_intercept_prob=p)
            errs.append(r['error_rate'])
        mean_errors.append(np.mean(errs))

    ax.plot(intercept_probs * 100, np.array(mean_errors) * 100, 'o-', color='#2196F3', lw=2, markersize=5)
    ax.axhline(11, color='black', ls='--', lw=2, label='Security threshold')
    ax.axhline(25, color='gray', ls=':', lw=1, label='Theoretical max (25%)')
    ax.set_xlabel("Eve's Intercept Probability (%)")
    ax.set_ylabel('Mean Error Rate (%)')
    ax.set_title('Error Rate vs Eavesdropping Intensity')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('demo_bb84.png', dpi=150, bbox_inches='tight')
    plt.show()

    # Print summary
    print(f"\n{'='*60}")
    print(f"BB84 Simulation Results ({n_trials} trials, {n_bits} qubits each)")
    print(f"{'='*60}")
    print(f"Without Eve:  mean error = {errors_no_eve.mean()*100:.2f}% ± {errors_no_eve.std()*100:.2f}%")
    print(f"With Eve:     mean error = {errors_with_eve.mean()*100:.2f}% ± {errors_with_eve.std()*100:.2f}%")
    print(f"Security threshold: 11%")
    print(f"Detection: {'YES ✓' if errors_with_eve.mean() > 0.11 else 'NO ✗'} — Eve is detectable!")
    print(f"{'='*60}")


if __name__ == '__main__':
    run_demo()
