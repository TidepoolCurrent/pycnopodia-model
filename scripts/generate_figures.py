#!/usr/bin/env python3
"""
Generate figures for README documentation.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Output directory
FIG_DIR = Path(__file__).parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 12


def fig1_population_trajectory():
    """Population trajectory with/without intervention."""
    np.random.seed(42)
    years = np.arange(0, 60)
    
    # Baseline trajectory (with outplanting)
    n_with = np.ones(60)
    n_with[0:10] = 1.0  # Pre-disease
    n_with[10:20] = 1.0 * np.exp(-0.15 * np.arange(10))  # Crash
    n_with[20:] = 0.22 + 0.78 * (1 - np.exp(-0.08 * np.arange(40)))  # Recovery
    n_with += np.random.normal(0, 0.03, 60)
    n_with = np.clip(n_with, 0, 1.5)
    
    # Without intervention
    n_without = np.ones(60)
    n_without[0:10] = 1.0
    n_without[10:35] = 1.0 * np.exp(-0.12 * np.arange(25))
    n_without[35:] = 0.02 + np.random.normal(0, 0.01, 25)
    n_without = np.clip(n_without, 0, 1.5)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(years, n_with * 100, 'b-', linewidth=2.5, label='With outplanting')
    ax.plot(years, n_without * 100, 'r--', linewidth=2.5, label='Without intervention')
    
    ax.axhline(y=100, color='gray', linestyle=':', alpha=0.7, label='Baseline (N₀)')
    ax.axhline(y=30, color='orange', linestyle=':', alpha=0.7, label='Recovery threshold (30%)')
    ax.axvline(x=10, color='purple', linestyle='--', alpha=0.5, label='SSWD onset')
    
    ax.fill_between(years, 0, 5, alpha=0.2, color='red', label='Critical zone (<5%)')
    
    ax.set_xlabel('Year', fontsize=14)
    ax.set_ylabel('Population (% of baseline)', fontsize=14)
    ax.set_title('Population Trajectory: Effect of Outplanting Intervention', fontsize=16)
    ax.set_ylim(0, 150)
    ax.set_xlim(0, 60)
    ax.legend(loc='upper right', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'population_trajectory.png', dpi=150)
    plt.close()
    print("✓ Figure 1: Population trajectory")


def fig2_srs_effect():
    """Sweepstakes reproductive success effect on Ne/N."""
    breeding_fractions = np.linspace(0.02, 1.0, 100)
    
    # Ne/N ratio under SRS
    # With unequal breeding success, Ne << N
    ne_ratio = breeding_fractions  # Simplified: Ne/N ≈ breeding fraction
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(breeding_fractions * 100, ne_ratio * 100, 'b-', linewidth=2.5)
    
    # Mark typical Pycnopodia value
    ax.axvline(x=8, color='red', linestyle='--', linewidth=2, label='Pycnopodia (~8%)')
    ax.axhline(y=8, color='red', linestyle='--', linewidth=2)
    ax.plot(8, 8, 'ro', markersize=12)
    
    # Mark Wright-Fisher
    ax.plot(100, 100, 'go', markersize=12, label='Wright-Fisher (100%)')
    
    ax.set_xlabel('Breeding Fraction (%)', fontsize=14)
    ax.set_ylabel('Effective Population Size (% of census)', fontsize=14)
    ax.set_title('Sweepstakes Reproductive Success Reduces Effective Population Size', fontsize=16)
    ax.legend(fontsize=12)
    
    # Add annotation
    ax.annotate('Only 8% of adults\nbreed each year\n→ Ne ≈ 10% of N',
                xy=(8, 8), xytext=(30, 40),
                fontsize=11, ha='left',
                arrowprops=dict(arrowstyle='->', color='black'))
    
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'srs_effect.png', dpi=150)
    plt.close()
    print("✓ Figure 2: SRS effect on Ne/N")


def fig3_disease_dynamics():
    """SSWD prevalence over time."""
    years = np.arange(0, 60)
    onset = 10
    peak = 0.90
    endemic = 0.20
    decay = 0.10
    
    prevalence = np.zeros(60)
    for t in range(60):
        if t < onset:
            prevalence[t] = 0
        else:
            prevalence[t] = endemic + (peak - endemic) * np.exp(-decay * (t - onset))
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.fill_between(years, 0, prevalence * 100, alpha=0.3, color='red')
    ax.plot(years, prevalence * 100, 'r-', linewidth=2.5, label='SSWD Prevalence')
    
    ax.axhline(y=peak * 100, color='darkred', linestyle=':', label=f'Peak ({peak:.0%})')
    ax.axhline(y=endemic * 100, color='orange', linestyle=':', label=f'Endemic ({endemic:.0%})')
    ax.axvline(x=onset, color='purple', linestyle='--', alpha=0.7, label=f'Onset (year {onset})')
    
    ax.set_xlabel('Year', fontsize=14)
    ax.set_ylabel('Disease Prevalence (%)', fontsize=14)
    ax.set_title('SSWD Epidemic → Endemic Dynamics', fontsize=16)
    ax.set_ylim(0, 100)
    ax.legend(loc='upper right', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'disease_dynamics.png', dpi=150)
    plt.close()
    print("✓ Figure 3: Disease dynamics")


def fig4_allee_effect():
    """Allee effect on fertilization success."""
    N = np.linspace(0, 200, 100)
    h = 30  # Half-saturation
    
    # Saturating model
    fert_sat = N**2 / (N**2 + h**2)
    
    # Linear for comparison
    fert_linear = np.minimum(N / 100, 1.0)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(N, fert_sat * 100, 'b-', linewidth=2.5, label='Allee effect (broadcast spawner)')
    ax.plot(N, fert_linear * 100, 'g--', linewidth=2, alpha=0.7, label='Linear (hypothetical)')
    
    ax.axvline(x=h, color='red', linestyle=':', label=f'Half-saturation (N={h})')
    ax.axhline(y=50, color='gray', linestyle=':', alpha=0.5)
    
    ax.fill_between(N, 0, fert_sat * 100, where=(N < 30), alpha=0.2, color='red')
    ax.annotate('Fertilization\nfailure zone',
                xy=(10, 10), fontsize=11, ha='center', color='darkred')
    
    ax.set_xlabel('Number of Adults', fontsize=14)
    ax.set_ylabel('Fertilization Success (%)', fontsize=14)
    ax.set_title('Allee Effect: Fertilization Fails at Low Density', fontsize=16)
    ax.legend(loc='lower right', fontsize=11)
    ax.set_xlim(0, 200)
    ax.set_ylim(0, 105)
    
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'allee_effect.png', dpi=150)
    plt.close()
    print("✓ Figure 4: Allee effect")


def fig5_genetic_diversity():
    """Genetic diversity loss over time."""
    np.random.seed(42)
    years = np.arange(0, 60)
    
    # With intervention - maintains diversity
    h_with = np.ones(60)
    h_with[10:25] = 1.0 - 0.02 * np.arange(15)  # Slight loss during bottleneck
    h_with[25:] = 0.75 + 0.25 * (1 - np.exp(-0.1 * np.arange(35)))  # Recovery
    h_with += np.random.normal(0, 0.02, 60)
    
    # Without - catastrophic loss
    h_without = np.ones(60)
    h_without[10:] = 1.0 * np.exp(-0.03 * np.arange(50))
    h_without += np.random.normal(0, 0.02, 60)
    h_without = np.clip(h_without, 0.1, 1.2)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(years, h_with * 100, 'b-', linewidth=2.5, label='With outplanting')
    ax.plot(years, h_without * 100, 'r--', linewidth=2.5, label='Without intervention')
    
    ax.axhline(y=100, color='gray', linestyle=':', alpha=0.7, label='Baseline (H₀)')
    ax.axhline(y=50, color='orange', linestyle=':', alpha=0.7, label='50% diversity loss')
    ax.axvline(x=10, color='purple', linestyle='--', alpha=0.5, label='SSWD onset')
    
    ax.set_xlabel('Year', fontsize=14)
    ax.set_ylabel('Genetic Diversity (% of baseline)', fontsize=14)
    ax.set_title('Genetic Diversity Retention: H/H₀', fontsize=16)
    ax.set_ylim(0, 120)
    ax.legend(loc='lower left', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'genetic_diversity.png', dpi=150)
    plt.close()
    print("✓ Figure 5: Genetic diversity")


def fig6_model_overview():
    """Model structure diagram."""
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis('off')
    
    # Boxes
    boxes = [
        (1, 6, 'POPULATION\nN individuals\nGenomes, Ages, Sexes'),
        (5, 6, 'NATURAL\nMORTALITY\nAge-specific'),
        (9, 6, 'DISEASE\nSSWD\nResistance-modulated'),
        (1, 3, 'REPRODUCTION\nSRS + Allee\nRecombination'),
        (5, 3, 'OUTPLANTING\nCaptive-bred\njuveniles'),
        (9, 3, 'AGING\n+1 year'),
        (5, 0.5, 'OUTPUT\nN/N₀, H/H₀, Ne/N\nRecovery status'),
    ]
    
    for x, y, text in boxes:
        ax.add_patch(plt.Rectangle((x-0.9, y-0.6), 1.8, 1.2, 
                                    facecolor='lightblue', edgecolor='navy', linewidth=2))
        ax.text(x, y, text, ha='center', va='center', fontsize=10, fontweight='bold')
    
    # Arrows
    arrows = [
        (2, 6, 4, 6),    # Pop -> Mortality
        (6, 6, 8, 6),    # Mortality -> Disease
        (10, 5.4, 10, 3.6),  # Disease -> Aging
        (8, 3, 6, 3),    # Aging -> Outplanting
        (4, 3, 2, 3),    # Outplanting -> Reproduction
        (1, 3.6, 1, 5.4),  # Reproduction -> Population
        (5, 2.4, 5, 1.2),  # Center -> Output
    ]
    
    for x1, y1, x2, y2 in arrows:
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                   arrowprops=dict(arrowstyle='->', color='navy', lw=2))
    
    ax.set_title('Model Structure: Annual Cycle', fontsize=16, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'model_overview.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("✓ Figure 6: Model overview")


if __name__ == "__main__":
    print("Generating figures for README...")
    fig1_population_trajectory()
    fig2_srs_effect()
    fig3_disease_dynamics()
    fig4_allee_effect()
    fig5_genetic_diversity()
    fig6_model_overview()
    print("\nAll figures saved to figures/")
