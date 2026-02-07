#!/usr/bin/env python3
"""
Comprehensive investigation of the pycnopodia network model.
Generates publication-quality figures for all model outputs.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import networkx as nx
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Tuple
import sys

# Add repo to path
sys.path.insert(0, str(Path(__file__).parent))

from pycnopodia.network import (
    NetworkConfig, NetworkSimulation, ConnectivityType,
    build_connectivity_matrix
)

# Output directory
FIGURES_DIR = Path(__file__).parent / "figures" / "investigation"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Plotting style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 10

# Color schemes
TOPOLOGY_COLORS = {
    'STEPPING_STONE': '#2ecc71',  # green
    'DISTANCE_DECAY': '#3498db',  # blue
    'UNIFORM': '#9b59b6',         # purple
    'HUB_NETWORK': '#e74c3c',     # red
}


def save_figure(fig, name: str):
    """Save figure with consistent settings."""
    path = FIGURES_DIR / f"{name}.png"
    fig.savefig(path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: {path}")


# ==============================================================================
# PART 1: NETWORK TOPOLOGIES
# ==============================================================================

def visualize_network_structure(connectivity_type: ConnectivityType, n_sites: int = 100):
    """Create graph visualization of network topology."""
    config = NetworkConfig(
        n_sites=n_sites,
        connectivity_type=connectivity_type,
        n_hubs=5 if connectivity_type == ConnectivityType.HUB_NETWORK else 10,
    )
    C = build_connectivity_matrix(config)
    
    # Create networkx graph
    G = nx.DiGraph()
    G.add_nodes_from(range(n_sites))
    
    # Add edges (threshold at 0.01 for visibility)
    for i in range(n_sites):
        for j in range(n_sites):
            if i != j and C[i, j] > 0.01:
                G.add_edge(i, j, weight=C[i, j])
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Position nodes in a circle
    pos = nx.circular_layout(G)
    
    # Draw edges with varying alpha based on weight
    edges = G.edges(data=True)
    weights = [d['weight'] for u, v, d in edges]
    max_weight = max(weights) if weights else 1
    
    # Sample edges for large graphs
    if len(edges) > 1000:
        sampled_edges = list(edges)[:1000]
    else:
        sampled_edges = list(edges)
    
    for u, v, d in sampled_edges:
        alpha = min(d['weight'] / max_weight * 2, 0.8)
        ax.annotate("", xy=pos[v], xytext=pos[u],
                   arrowprops=dict(arrowstyle="->", alpha=alpha, 
                                   color=TOPOLOGY_COLORS.get(connectivity_type.name, 'gray'),
                                   lw=0.5))
    
    # Draw nodes
    node_colors = TOPOLOGY_COLORS.get(connectivity_type.name, 'gray')
    nx.draw_networkx_nodes(G, pos, node_size=30, node_color=node_colors, 
                          alpha=0.8, ax=ax)
    
    ax.set_title(f"Network Topology: {connectivity_type.value.upper()}\n({n_sites} sites)", 
                fontsize=14, fontweight='bold')
    ax.axis('off')
    
    return fig


def run_topology_simulations():
    """Run simulations for all topology types."""
    print("\n" + "="*60)
    print("PART 1: NETWORK TOPOLOGY ANALYSIS")
    print("="*60)
    
    topologies = [
        ConnectivityType.STEPPING_STONE,
        ConnectivityType.DISTANCE_DECAY,
        ConnectivityType.UNIFORM,
        ConnectivityType.HUB_NETWORK,
    ]
    
    results = {}
    
    for conn_type in topologies:
        print(f"\n  Running {conn_type.name}...")
        
        # Network structure visualization (smaller for speed)
        fig = visualize_network_structure(conn_type, n_sites=50)
        save_figure(fig, f"topology_{conn_type.value}_network")
        
        # Run simulation
        config = NetworkConfig(
            connectivity_type=conn_type,
            n_sites=200,  # Smaller for faster runs
            n_years=100,
            n_hubs=10 if conn_type == ConnectivityType.HUB_NETWORK else 10,
        )
        sim = NetworkSimulation(config, seed=42)
        result = sim.run()
        results[conn_type.name] = result
        
        # Population trajectory
        fig, ax = plt.subplots(figsize=(10, 6))
        years = [s.year for s in result.states]
        n_ratios = [s.n_ratio for s in result.states]
        ax.plot(years, n_ratios, color=TOPOLOGY_COLORS[conn_type.name], lw=2)
        ax.axhline(y=0.3, color='gray', ls='--', alpha=0.7, label='Recovery threshold (30%)')
        ax.axvline(x=config.disease_onset_year, color='red', ls=':', alpha=0.7, label='Disease onset')
        ax.set_xlabel('Year')
        ax.set_ylabel('Population (N/N₀)')
        ax.set_title(f'Population Trajectory: {conn_type.value.upper()}')
        ax.legend()
        ax.set_ylim(0, 1.1)
        save_figure(fig, f"topology_{conn_type.value}_trajectory")
        
        # Spatial distribution at timepoints
        timepoints = [0, 25, 50, 100]
        fig, axes = plt.subplots(2, 4, figsize=(16, 8))
        
        for idx, year in enumerate(timepoints):
            state = result.states[min(year, len(result.states)-1)]
            
            # Resistance allele frequency (mean across loci for polygenic model)
            ax = axes[0, idx]
            # Handle both 1D (legacy) and 2D (polygenic) resistance_freqs
            if state.resistance_freqs.ndim == 2:
                site_resistance = state.resistance_freqs.mean(axis=1)
            else:
                site_resistance = state.resistance_freqs
            sites = np.arange(len(site_resistance))
            ax.bar(sites, site_resistance, color=TOPOLOGY_COLORS[conn_type.name], alpha=0.7)
            ax.set_xlabel('Site')
            ax.set_ylabel('Resistance Freq')
            ax.set_title(f'Year {year}')
            ax.set_ylim(0, 1)
            
            # Population density
            ax = axes[1, idx]
            density = state.populations / config.carrying_capacity_per_site
            ax.bar(sites, density, color=TOPOLOGY_COLORS[conn_type.name], alpha=0.7)
            ax.set_xlabel('Site')
            ax.set_ylabel('Pop Density (N/K)')
            ax.set_ylim(0, 1.2)
        
        fig.suptitle(f'Spatial Distribution: {conn_type.value.upper()}', fontsize=14, fontweight='bold')
        axes[0, 0].set_ylabel('Resistance Allele\nFrequency')
        axes[1, 0].set_ylabel('Population Density\n(N/K)')
        plt.tight_layout()
        save_figure(fig, f"topology_{conn_type.value}_spatial")
    
    # Combined trajectory comparison
    fig, ax = plt.subplots(figsize=(12, 7))
    for conn_type, result in results.items():
        years = [s.year for s in result.states]
        n_ratios = [s.n_ratio for s in result.states]
        ax.plot(years, n_ratios, color=TOPOLOGY_COLORS[conn_type], lw=2.5, label=conn_type)
    
    ax.axhline(y=0.3, color='gray', ls='--', alpha=0.7, label='Recovery threshold')
    ax.axvline(x=10, color='red', ls=':', alpha=0.7, label='Disease onset')
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Population (N/N₀)', fontsize=12)
    ax.set_title('Population Trajectories by Network Topology', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    ax.set_ylim(0, 1.1)
    save_figure(fig, "topology_comparison_trajectories")
    
    return results


# ==============================================================================
# PART 2: PARAMETER SENSITIVITY
# ==============================================================================

def run_disease_mortality_sweep():
    """Sweep disease mortality rates."""
    print("\n  Disease mortality sweep...")
    mortalities = [0.80, 0.90, 0.95, 0.99]
    results = {}
    
    for mort in mortalities:
        config = NetworkConfig(
            n_sites=200,
            n_years=100,
            disease_mortality=mort,
        )
        sim = NetworkSimulation(config, seed=42)
        results[mort] = sim.run()
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 7))
    colors = plt.cm.Reds(np.linspace(0.3, 0.9, len(mortalities)))
    
    for (mort, result), color in zip(results.items(), colors):
        years = [s.year for s in result.states]
        n_ratios = [s.n_ratio for s in result.states]
        ax.plot(years, n_ratios, color=color, lw=2.5, label=f'{mort:.0%} mortality')
    
    ax.axhline(y=0.3, color='gray', ls='--', alpha=0.7)
    ax.axvline(x=10, color='red', ls=':', alpha=0.5)
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Population (N/N₀)', fontsize=12)
    ax.set_title('Disease Mortality Impact on Population', fontsize=14, fontweight='bold')
    ax.legend()
    ax.set_ylim(0, 1.1)
    save_figure(fig, "sensitivity_disease_mortality")
    
    return results


def run_self_recruitment_sweep():
    """Sweep self-recruitment rates."""
    print("\n  Self-recruitment sweep...")
    self_recs = [0.1, 0.3, 0.5, 0.7, 0.9]
    results = {}
    
    for sr in self_recs:
        config = NetworkConfig(
            n_sites=200,
            n_years=100,
            self_recruitment=sr,
        )
        sim = NetworkSimulation(config, seed=42)
        results[sr] = sim.run()
    
    # Plot trajectories
    fig, ax = plt.subplots(figsize=(12, 7))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(self_recs)))
    
    for (sr, result), color in zip(results.items(), colors):
        years = [s.year for s in result.states]
        n_ratios = [s.n_ratio for s in result.states]
        ax.plot(years, n_ratios, color=color, lw=2.5, label=f'{sr:.0%} self-recruitment')
    
    ax.axhline(y=0.3, color='gray', ls='--', alpha=0.7)
    ax.axvline(x=10, color='red', ls=':', alpha=0.5)
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Population (N/N₀)', fontsize=12)
    ax.set_title('Self-Recruitment Impact on Recovery', fontsize=14, fontweight='bold')
    ax.legend()
    ax.set_ylim(0, 1.1)
    save_figure(fig, "sensitivity_self_recruitment")
    
    # Final population bar chart
    fig, ax = plt.subplots(figsize=(10, 6))
    final_pops = [results[sr].final_n_ratio for sr in self_recs]
    bars = ax.bar([f'{sr:.0%}' for sr in self_recs], final_pops, color=colors)
    ax.axhline(y=0.3, color='red', ls='--', alpha=0.7, label='Recovery threshold')
    ax.set_xlabel('Self-Recruitment Rate', fontsize=12)
    ax.set_ylabel('Final Population (N/N₀)', fontsize=12)
    ax.set_title('Final Population by Self-Recruitment Rate', fontsize=14, fontweight='bold')
    ax.legend()
    save_figure(fig, "sensitivity_self_recruitment_final")
    
    return results


def run_outplanting_sweep():
    """Sweep outplanting intensity."""
    print("\n  Outplanting intensity sweep...")
    intensities = [0, 100, 500, 1000, 5000]
    results = {}
    
    # Define outplanting sites (every 10th site)
    outplanting_sites = list(range(0, 200, 10))
    
    for intensity in intensities:
        config = NetworkConfig(
            n_sites=200,
            n_years=100,
            outplanting_n=intensity,
            outplanting_sites=outplanting_sites if intensity > 0 else [],
            outplanting_start=15,
        )
        sim = NetworkSimulation(config, seed=42)
        results[intensity] = sim.run()
    
    # Plot trajectories
    fig, ax = plt.subplots(figsize=(12, 7))
    colors = plt.cm.Greens(np.linspace(0.2, 0.9, len(intensities)))
    
    for (intensity, result), color in zip(results.items(), colors):
        years = [s.year for s in result.states]
        n_ratios = [s.n_ratio for s in result.states]
        label = 'No outplanting' if intensity == 0 else f'{intensity} per event'
        ax.plot(years, n_ratios, color=color, lw=2.5, label=label)
    
    ax.axhline(y=0.3, color='gray', ls='--', alpha=0.7)
    ax.axvline(x=15, color='green', ls=':', alpha=0.5, label='Outplanting starts')
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Population (N/N₀)', fontsize=12)
    ax.set_title('Outplanting Intensity Impact', fontsize=14, fontweight='bold')
    ax.legend()
    ax.set_ylim(0, 1.1)
    save_figure(fig, "sensitivity_outplanting_intensity")
    
    # Final population bar chart
    fig, ax = plt.subplots(figsize=(10, 6))
    final_pops = [results[i].final_n_ratio for i in intensities]
    labels = ['None'] + [str(i) for i in intensities[1:]]
    bars = ax.bar(labels, final_pops, color=colors)
    ax.axhline(y=0.3, color='red', ls='--', alpha=0.7, label='Recovery threshold')
    ax.set_xlabel('Outplanting Intensity (per event)', fontsize=12)
    ax.set_ylabel('Final Population (N/N₀)', fontsize=12)
    ax.set_title('Final Population by Outplanting Intensity', fontsize=14, fontweight='bold')
    ax.legend()
    save_figure(fig, "sensitivity_outplanting_final")
    
    return results


def run_resistance_freq_sweep():
    """Sweep initial resistance allele frequency."""
    print("\n  Resistance allele frequency sweep...")
    freqs = [0.01, 0.03, 0.05, 0.10]
    results = {}
    extinction_counts = {}
    
    n_replicates = 20
    
    for freq in freqs:
        replicate_results = []
        for rep in range(n_replicates):
            config = NetworkConfig(
                n_sites=200,
                n_years=100,
                initial_resistance_freq=freq,
            )
            sim = NetworkSimulation(config, seed=42 + rep)
            replicate_results.append(sim.run())
        
        results[freq] = replicate_results
        extinction_counts[freq] = sum(1 for r in replicate_results if r.extinct)
    
    # Plot extinction risk
    fig, ax = plt.subplots(figsize=(10, 6))
    freqs_pct = [f'{f:.0%}' for f in freqs]
    extinctions = [extinction_counts[f] / n_replicates * 100 for f in freqs]
    colors = plt.cm.RdYlGn_r(np.linspace(0.2, 0.8, len(freqs)))
    
    bars = ax.bar(freqs_pct, extinctions, color=colors)
    ax.set_xlabel('Initial Resistance Allele Frequency', fontsize=12)
    ax.set_ylabel('Extinction Risk (%)', fontsize=12)
    ax.set_title(f'Extinction Risk by Initial Resistance Frequency\n(n={n_replicates} replicates)', 
                fontsize=14, fontweight='bold')
    ax.set_ylim(0, 100)
    
    # Add value labels
    for bar, val in zip(bars, extinctions):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2, 
               f'{val:.0f}%', ha='center', va='bottom', fontsize=11)
    
    save_figure(fig, "sensitivity_resistance_freq_extinction")
    
    # Mean trajectory plot
    fig, ax = plt.subplots(figsize=(12, 7))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(freqs)))
    
    for (freq, reps), color in zip(results.items(), colors):
        # Average across replicates
        all_trajectories = []
        for r in reps:
            all_trajectories.append([s.n_ratio for s in r.states])
        mean_traj = np.mean(all_trajectories, axis=0)
        std_traj = np.std(all_trajectories, axis=0)
        years = list(range(len(mean_traj)))
        
        ax.plot(years, mean_traj, color=color, lw=2.5, label=f'{freq:.0%} initial freq')
        ax.fill_between(years, mean_traj - std_traj, mean_traj + std_traj, 
                       color=color, alpha=0.2)
    
    ax.axhline(y=0.3, color='gray', ls='--', alpha=0.7)
    ax.axvline(x=10, color='red', ls=':', alpha=0.5)
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Population (N/N₀)', fontsize=12)
    ax.set_title('Mean Population Trajectory by Resistance Frequency\n(±1 SD)', 
                fontsize=14, fontweight='bold')
    ax.legend()
    ax.set_ylim(0, 1.1)
    save_figure(fig, "sensitivity_resistance_freq_trajectories")
    
    return results


def run_sensitivity_analysis():
    """Run all parameter sensitivity analyses."""
    print("\n" + "="*60)
    print("PART 2: PARAMETER SENSITIVITY ANALYSIS")
    print("="*60)
    
    run_disease_mortality_sweep()
    run_self_recruitment_sweep()
    run_outplanting_sweep()
    run_resistance_freq_sweep()


# ==============================================================================
# PART 3: INTERVENTION COMPARISONS
# ==============================================================================

def run_intervention_comparison():
    """Compare different intervention strategies."""
    print("\n" + "="*60)
    print("PART 3: INTERVENTION COMPARISONS")
    print("="*60)
    
    n_sites = 200
    outplanting_sites = list(range(0, n_sites, 10))  # Every 10th site
    
    interventions = {
        'No Intervention': {
            'outplanting_n': 0,
            'outplanting_sites': [],
        },
        'Wild-caught': {
            'outplanting_n': 500,
            'outplanting_sites': outplanting_sites,
            'outplant_resistance_mode': 'wild',
        },
        'Enhanced (50% resistant)': {
            'outplanting_n': 500,
            'outplanting_sites': outplanting_sites,
            'outplant_resistance_mode': 'enhanced',
            'outplant_enhanced_resistance': 0.5,
        },
        'Enhanced (95% resistant)': {
            'outplanting_n': 500,
            'outplanting_sites': outplanting_sites,
            'outplant_resistance_mode': 'enhanced',
            'outplant_enhanced_resistance': 0.95,
        },
    }
    
    results = {}
    
    print("\n  Running intervention comparisons...")
    for name, params in interventions.items():
        config = NetworkConfig(
            n_sites=n_sites,
            n_years=100,
            outplanting_start=15,
            **params
        )
        sim = NetworkSimulation(config, seed=42)
        results[name] = sim.run()
    
    # Plot comparison - population trajectory
    fig, ax = plt.subplots(figsize=(12, 7))
    colors = ['#e74c3c', '#3498db', '#2ecc71', '#9b59b6']
    
    for (name, result), color in zip(results.items(), colors):
        years = [s.year for s in result.states]
        n_ratios = [s.n_ratio for s in result.states]
        ax.plot(years, n_ratios, color=color, lw=2.5, label=name)
    
    ax.axhline(y=0.3, color='gray', ls='--', alpha=0.7, label='Recovery threshold')
    ax.axvline(x=15, color='green', ls=':', alpha=0.5, label='Outplanting starts')
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Population (N/N₀)', fontsize=12)
    ax.set_title('Intervention Strategy Comparison', fontsize=14, fontweight='bold')
    ax.legend()
    ax.set_ylim(0, 1.1)
    save_figure(fig, "intervention_comparison")
    
    # NEW: Resistance evolution comparison figure
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Population trajectory (log scale for bottleneck visibility)
    ax = axes[0]
    for (name, result), color in zip(results.items(), colors):
        years = [s.year for s in result.states]
        n_ratios = [max(s.n_ratio, 1e-4) for s in result.states]  # Floor for log
        ax.plot(years, n_ratios, color=color, lw=2.5, label=name)
    
    ax.axhline(y=0.3, color='gray', ls='--', alpha=0.7)
    ax.axvline(x=10, color='red', ls=':', alpha=0.5, label='Disease onset')
    ax.axvline(x=15, color='green', ls=':', alpha=0.5, label='Outplanting starts')
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Population (N/N₀, log scale)', fontsize=12)
    ax.set_title('Population Trajectory', fontsize=12, fontweight='bold')
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 1.5)
    ax.legend(loc='lower right')
    
    # Resistance evolution
    ax = axes[1]
    for (name, result), color in zip(results.items(), colors):
        years = [s.year for s in result.states]
        resistance = [s.mean_resistance_freq * 100 for s in result.states]
        ax.plot(years, resistance, color=color, lw=2.5, label=name)
    
    ax.axhline(y=95, color='gray', ls='--', alpha=0.7, label='95% ceiling')
    ax.axvline(x=10, color='red', ls=':', alpha=0.5)
    ax.axvline(x=15, color='green', ls=':', alpha=0.5)
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Resistance Allele Frequency (%)', fontsize=12)
    ax.set_title('Resistance Evolution', fontsize=12, fontweight='bold')
    ax.set_ylim(0, 100)
    ax.legend(loc='lower right')
    
    plt.tight_layout()
    save_figure(fig, "intervention_resistance_comparison")
    
    # Print summary statistics
    print("\n  Intervention Summary:")
    print("  " + "-"*60)
    print(f"  {'Strategy':<30} {'Final N/N₀':>12} {'Final R%':>10}")
    print("  " + "-"*60)
    for name, result in results.items():
        final_n = result.final_n_ratio
        final_r = result.states[-1].mean_resistance_freq * 100
        print(f"  {name:<30} {final_n:>11.2%} {final_r:>9.1f}%")
    print("  " + "-"*60)
    
    # Timing comparison
    print("\n  Running timing comparison...")
    timing_results = {}
    
    for start_year in [5, 10, 15, 20]:
        config = NetworkConfig(
            n_sites=n_sites,
            n_years=100,
            outplanting_n=500,
            outplanting_sites=outplanting_sites,
            outplanting_start=start_year,
            outplant_resistance_mode='enhanced',
            outplant_enhanced_resistance=0.5,
        )
        sim = NetworkSimulation(config, seed=42)
        timing_results[start_year] = sim.run()
    
    fig, ax = plt.subplots(figsize=(12, 7))
    colors = plt.cm.plasma(np.linspace(0.2, 0.8, len(timing_results)))
    
    for (year, result), color in zip(timing_results.items(), colors):
        years = [s.year for s in result.states]
        n_ratios = [s.n_ratio for s in result.states]
        ax.plot(years, n_ratios, color=color, lw=2.5, label=f'Start year {year}')
        ax.axvline(x=year, color=color, ls=':', alpha=0.5)
    
    ax.axhline(y=0.3, color='gray', ls='--', alpha=0.7)
    ax.axvline(x=10, color='red', ls=':', alpha=0.7, label='Disease onset')
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Population (N/N₀)', fontsize=12)
    ax.set_title('Outplanting Timing: Early vs Late Start', fontsize=14, fontweight='bold')
    ax.legend()
    ax.set_ylim(0, 1.1)
    save_figure(fig, "intervention_timing")
    
    # Site selection comparison
    print("\n  Running site selection comparison...")
    site_strategies = {
        'Every 10th site': list(range(0, n_sites, 10)),
        'Clustered (0-30)': list(range(0, 30)),
        'Random': list(np.random.default_rng(42).choice(n_sites, 20, replace=False)),
        'Distributed (edges+center)': list(range(0, 10)) + list(range(95, 105)) + list(range(190, 200)),
    }
    
    site_results = {}
    for name, sites in site_strategies.items():
        config = NetworkConfig(
            n_sites=n_sites,
            n_years=100,
            outplanting_n=500,
            outplanting_sites=sites,
            outplanting_start=15,
            outplant_resistance_mode='enhanced',
            outplant_enhanced_resistance=0.5,
        )
        sim = NetworkSimulation(config, seed=42)
        site_results[name] = sim.run()
    
    fig, ax = plt.subplots(figsize=(12, 7))
    colors = ['#2ecc71', '#e74c3c', '#9b59b6', '#f39c12']
    
    for (name, result), color in zip(site_results.items(), colors):
        years = [s.year for s in result.states]
        n_ratios = [s.n_ratio for s in result.states]
        ax.plot(years, n_ratios, color=color, lw=2.5, label=name)
    
    ax.axhline(y=0.3, color='gray', ls='--', alpha=0.7)
    ax.axvline(x=15, color='green', ls=':', alpha=0.5)
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Population (N/N₀)', fontsize=12)
    ax.set_title('Outplanting Site Selection Strategy Comparison', fontsize=14, fontweight='bold')
    ax.legend()
    ax.set_ylim(0, 1.1)
    save_figure(fig, "intervention_site_selection")
    
    return results


# ==============================================================================
# PART 4: INDIVIDUAL-BASED MODEL
# ==============================================================================

def run_individual_based_model():
    """Run the individual-based model and compare with network model."""
    print("\n" + "="*60)
    print("PART 4: INDIVIDUAL-BASED MODEL")
    print("="*60)
    
    # Try to import and run the individual-based model
    try:
        from pycnopodia.config import PRESETS
        from pycnopodia.simulation import run_replicates, summarize_replicates
        
        print("\n  Running individual-based model (baseline)...")
        config = PRESETS['baseline']
        config.n_replicates = 10
        config.n_years = 80
        
        ibm_results = run_replicates(config, progress=False)
        summary = summarize_replicates(ibm_results, recovery_threshold=0.3)
        
        # Plot IBM trajectories
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Population trajectories
        ax = axes[0]
        for i, result in enumerate(ibm_results[:5]):
            ax.plot(result.n_ratio, alpha=0.7, lw=1.5, label=f'Rep {i+1}' if i < 3 else None)
        
        mean_traj = np.mean([r.n_ratio for r in ibm_results], axis=0)
        ax.plot(mean_traj, color='black', lw=3, label='Mean')
        ax.axhline(y=0.3, color='red', ls='--', alpha=0.7, label='Recovery threshold')
        ax.set_xlabel('Year', fontsize=12)
        ax.set_ylabel('Population (N/N₀)', fontsize=12)
        ax.set_title('Individual-Based Model: Population Trajectories', fontsize=12, fontweight='bold')
        ax.legend()
        ax.set_ylim(0, 1.5)
        
        # Genetic diversity
        ax = axes[1]
        for i, result in enumerate(ibm_results[:5]):
            ax.plot(result.h_ratio, alpha=0.7, lw=1.5)
        
        mean_h = np.mean([r.h_ratio for r in ibm_results], axis=0)
        ax.plot(mean_h, color='black', lw=3, label='Mean')
        ax.set_xlabel('Year', fontsize=12)
        ax.set_ylabel('Genetic Diversity (H/H₀)', fontsize=12)
        ax.set_title('Individual-Based Model: Genetic Diversity', fontsize=12, fontweight='bold')
        ax.legend()
        ax.set_ylim(0, 1.5)
        
        plt.tight_layout()
        save_figure(fig, "ibm_trajectories")
        
        # Summary statistics
        fig, ax = plt.subplots(figsize=(10, 6))
        
        metrics = ['Extinction\nProb', 'Recovery\nProb', 'Final\nN/N₀', 'Final\nH/H₀']
        values = [
            summary['extinction_probability'] * 100,
            summary['recovery_probability'] * 100,
            summary['mean_final_n_ratio'] * 100,
            summary['mean_final_h_ratio'] * 100,
        ]
        colors = ['#e74c3c', '#2ecc71', '#3498db', '#9b59b6']
        
        bars = ax.bar(metrics, values, color=colors)
        ax.set_ylabel('Percentage (%)', fontsize=12)
        ax.set_title('Individual-Based Model: Summary Statistics\n(Baseline preset, n=10 replicates)', 
                    fontsize=14, fontweight='bold')
        ax.set_ylim(0, 110)
        
        for bar, val in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2, 
                   f'{val:.1f}%', ha='center', va='bottom', fontsize=11)
        
        save_figure(fig, "ibm_summary")
        
        print(f"\n  IBM Summary:")
        print(f"    Extinction probability: {summary['extinction_probability']:.1%}")
        print(f"    Recovery probability: {summary['recovery_probability']:.1%}")
        print(f"    Final N/N₀: {summary['mean_final_n_ratio']:.2%}")
        print(f"    Final H/H₀: {summary['mean_final_h_ratio']:.2%}")
        
        return ibm_results, summary
        
    except Exception as e:
        print(f"  Warning: Could not run individual-based model: {e}")
        return None, None


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("="*60)
    print("PYCNOPODIA MODEL INVESTIGATION")
    print("="*60)
    print(f"\nOutput directory: {FIGURES_DIR}")
    
    # Part 1: Network topologies
    topology_results = run_topology_simulations()
    
    # Part 2: Parameter sensitivity
    run_sensitivity_analysis()
    
    # Part 3: Intervention comparisons
    run_intervention_comparison()
    
    # Part 4: Individual-based model
    ibm_results, ibm_summary = run_individual_based_model()
    
    print("\n" + "="*60)
    print("INVESTIGATION COMPLETE")
    print("="*60)
    print(f"\nAll figures saved to: {FIGURES_DIR}")
    
    # List generated figures
    figures = list(FIGURES_DIR.glob("*.png"))
    print(f"\nGenerated {len(figures)} figures:")
    for fig in sorted(figures):
        print(f"  - {fig.name}")


if __name__ == "__main__":
    main()
