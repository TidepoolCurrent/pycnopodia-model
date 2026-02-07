#!/usr/bin/env python3
"""
Visualize disease spread as a spatial wave through the network.

Creates a heatmap showing:
- X-axis: Site index (spatial position along coast)
- Y-axis: Year
- Color: Disease prevalence (0 = uninfected, 1 = fully infected)

Demonstrates connectivity-based disease spread pattern.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add parent to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from pycnopodia.network import (
    NetworkConfig, NetworkSimulation, ConnectivityType
)


def generate_disease_spread_visualization():
    """Generate disease spread wave visualization."""
    
    # Configure for clear wave visualization
    config = NetworkConfig(
        n_sites=200,  # More sites for smoother visualization
        n_years=25,
        disease_onset_year=5,
        disease_onset_site=100,  # Start in middle
        disease_connectivity_type=ConnectivityType.DISTANCE_DECAY,
        disease_dispersal_scale=30,  # Reasonable spread distance
        disease_asymmetry=0.2,  # Slight southward bias (typical for SSWD)
        disease_transmission_prob=0.7,
        disease_mortality=0.90,  # High but not 99% to see spread
        refugia_fraction=0.02,  # 2% refugia
    )
    
    print("Running simulation...")
    sim = NetworkSimulation(config, seed=42)
    result = sim.run()
    
    # Extract disease prevalence over time
    n_years = len(result.states)
    n_sites = config.n_sites
    prevalence_matrix = np.zeros((n_years, n_sites))
    
    for year, state in enumerate(result.states):
        prevalence_matrix[year, :] = state.disease_prevalence
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Full heatmap of disease spread
    ax1 = axes[0, 0]
    im1 = ax1.imshow(
        prevalence_matrix,
        aspect='auto',
        cmap='YlOrRd',
        origin='lower',
        extent=[0, n_sites, 0, n_years],
        vmin=0, vmax=1
    )
    ax1.set_xlabel('Site Index (spatial position)')
    ax1.set_ylabel('Year')
    ax1.set_title('Disease Spread as Spatial Wave')
    ax1.axhline(y=config.disease_onset_year, color='white', linestyle='--', 
                label='Disease onset', alpha=0.7)
    ax1.axvline(x=config.disease_onset_site, color='white', linestyle=':', 
                label='Origin site', alpha=0.7)
    plt.colorbar(im1, ax=ax1, label='Disease Prevalence')
    
    # Plot 2: Wave front progression
    ax2 = axes[0, 1]
    # Find wave front (first site with prevalence > 0.5 at each distance from origin)
    infected_fraction = (prevalence_matrix > 0.1).sum(axis=1) / n_sites
    ax2.plot(range(n_years), infected_fraction * 100, 'b-', linewidth=2)
    ax2.axhline(y=50, color='red', linestyle='--', alpha=0.5, label='50% infected')
    ax2.axhline(y=90, color='darkred', linestyle='--', alpha=0.5, label='90% infected')
    ax2.axvline(x=config.disease_onset_year, color='gray', linestyle=':', alpha=0.5)
    ax2.set_xlabel('Year')
    ax2.set_ylabel('Percent of Sites Infected')
    ax2.set_title('Disease Wave Front Progression')
    ax2.legend()
    ax2.set_xlim(0, n_years)
    ax2.set_ylim(0, 100)
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Cross-sections at different times
    ax3 = axes[1, 0]
    years_to_show = [config.disease_onset_year, 
                     config.disease_onset_year + 1,
                     config.disease_onset_year + 2,
                     config.disease_onset_year + 4,
                     min(config.disease_onset_year + 8, n_years - 1)]
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(years_to_show)))
    
    for year, color in zip(years_to_show, colors):
        ax3.plot(range(n_sites), prevalence_matrix[year, :], 
                color=color, label=f'Year {year}', linewidth=1.5)
    ax3.axvline(x=config.disease_onset_site, color='gray', linestyle=':', alpha=0.5)
    ax3.set_xlabel('Site Index')
    ax3.set_ylabel('Disease Prevalence')
    ax3.set_title('Disease Wave Cross-Sections')
    ax3.legend(loc='upper right')
    ax3.set_xlim(0, n_sites)
    ax3.set_ylim(0, 1.05)
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Population decline with disease spread
    ax4 = axes[1, 1]
    years = result.years
    n_ratio = result.n_ratio * 100
    occupied = result.occupied_ratio * 100
    infected = np.array([s.infected_sites_ratio for s in result.states]) * 100
    
    ax4.plot(years, n_ratio, 'b-', linewidth=2, label='Population %')
    ax4.plot(years, occupied, 'g--', linewidth=1.5, label='Sites occupied %')
    ax4.plot(years, infected, 'r:', linewidth=2, label='Sites infected %')
    ax4.axvline(x=config.disease_onset_year, color='gray', linestyle=':', alpha=0.5)
    ax4.set_xlabel('Year')
    ax4.set_ylabel('Percentage')
    ax4.set_title('Population Response to Disease Wave')
    ax4.legend(loc='right')
    ax4.set_xlim(0, n_years)
    ax4.set_ylim(0, 105)
    ax4.grid(True, alpha=0.3)
    
    # Overall title
    fig.suptitle(
        'Disease Connectivity-Based Spread Pattern\n'
        f'(Disease: {config.disease_connectivity_type.value}, '
        f'Dispersal scale: {config.disease_dispersal_scale}, '
        f'Transmission prob: {config.disease_transmission_prob})',
        fontsize=12, y=1.02
    )
    
    plt.tight_layout()
    
    # Save figure
    output_path = Path(__file__).parent.parent / 'figures' / 'investigation' / 'disease_spread_wave.png'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Saved to: {output_path}")
    
    # Also save a comparison figure showing different connectivity types
    generate_connectivity_comparison()
    
    return result


def generate_connectivity_comparison():
    """Compare disease spread under different connectivity types."""
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    connectivity_types = [
        ConnectivityType.STEPPING_STONE,
        ConnectivityType.DISTANCE_DECAY,
        ConnectivityType.ASYMMETRIC_FLOW,
        ConnectivityType.UNIFORM,
        ConnectivityType.HUB_NETWORK,
        ConnectivityType.MODULAR,
    ]
    
    for ax, conn_type in zip(axes.flat, connectivity_types):
        config = NetworkConfig(
            n_sites=100,
            n_years=20,
            disease_onset_year=5,
            disease_onset_site=50,
            disease_connectivity_type=conn_type,
            disease_dispersal_scale=15,
            disease_transmission_prob=0.7,
            disease_mortality=0.85,
        )
        
        sim = NetworkSimulation(config, seed=42)
        result = sim.run()
        
        # Build prevalence matrix
        prevalence_matrix = np.array([s.disease_prevalence for s in result.states])
        
        im = ax.imshow(
            prevalence_matrix,
            aspect='auto',
            cmap='YlOrRd',
            origin='lower',
            extent=[0, 100, 0, 20],
            vmin=0, vmax=1
        )
        ax.set_title(f'{conn_type.value.title()}')
        ax.set_xlabel('Site')
        ax.set_ylabel('Year')
        ax.axhline(y=5, color='white', linestyle='--', alpha=0.5)
        ax.axvline(x=50, color='white', linestyle=':', alpha=0.5)
    
    fig.suptitle('Disease Spread Under Different Connectivity Types', fontsize=14)
    fig.colorbar(im, ax=axes.ravel().tolist(), label='Disease Prevalence', 
                 shrink=0.6, pad=0.02)
    
    plt.tight_layout()
    
    output_path = Path(__file__).parent.parent / 'figures' / 'investigation' / 'disease_connectivity_comparison.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Saved to: {output_path}")


if __name__ == '__main__':
    generate_disease_spread_visualization()
