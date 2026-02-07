#!/usr/bin/env python3
"""
Exhaustive interrogation and visualization of the Pycnopodia network model.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from pycnopodia.network import (
    NetworkSimulation, NetworkConfig, ConnectivityType,
    run_scenario, build_connectivity_matrix
)

# Output directory
OUT_DIR = Path(__file__).parent.parent / "figures"
OUT_DIR.mkdir(exist_ok=True)

plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11


def fig1_single_trajectory():
    """Figure 1: Single trajectory showing crash and recovery dynamics."""
    print("Figure 1: Single trajectory...")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Run with and without outplanting
    scenarios = [
        ("No intervention", NetworkConfig(n_sites=500, n_years=100)),
        ("5k/yr × 50 sites (wild)", NetworkConfig(
            n_sites=500, n_years=100,
            outplanting_n=5000, outplanting_sites=list(range(0, 500, 10)),
            outplanting_start=15, outplant_resistance_mode='wild'
        )),
        ("5k/yr × 50 sites (enhanced 95%)", NetworkConfig(
            n_sites=500, n_years=100,
            outplanting_n=5000, outplanting_sites=list(range(0, 500, 10)),
            outplanting_start=15, outplant_resistance_mode='enhanced',
            outplant_enhanced_resistance=0.95
        )),
    ]
    
    colors = ['#d62728', '#2ca02c', '#1f77b4']
    
    for (name, cfg), color in zip(scenarios, colors):
        sim = NetworkSimulation(cfg, seed=42)
        result = sim.run()
        
        years = [s.year for s in result.states]
        n_ratio = [s.n_ratio * 100 for s in result.states]
        resistance = [s.mean_resistance_freq * 100 for s in result.states]
        infected = [s.infected_sites_ratio * 100 for s in result.states]
        occupied = [s.occupied_ratio * 100 for s in result.states]
        
        axes[0, 0].plot(years, n_ratio, label=name, color=color, linewidth=2)
        axes[0, 1].plot(years, resistance, label=name, color=color, linewidth=2)
        axes[1, 0].plot(years, infected, label=name, color=color, linewidth=2)
        axes[1, 1].plot(years, occupied, label=name, color=color, linewidth=2)
    
    axes[0, 0].set_ylabel('Population (% of baseline)')
    axes[0, 0].set_title('Population Trajectory')
    axes[0, 0].axvline(x=10, color='gray', linestyle='--', alpha=0.5, label='Disease onset')
    axes[0, 0].axvline(x=15, color='gray', linestyle=':', alpha=0.5, label='Outplanting start')
    axes[0, 0].set_yscale('log')
    axes[0, 0].set_ylim(0.001, 100)
    axes[0, 0].legend(loc='upper right')
    
    axes[0, 1].set_ylabel('Resistance allele frequency (%)')
    axes[0, 1].set_title('Resistance Evolution')
    axes[0, 1].axhline(y=95, color='gray', linestyle='--', alpha=0.5, label='95% ceiling')
    axes[0, 1].legend(loc='lower right')
    
    axes[1, 0].set_ylabel('Sites infected (%)')
    axes[1, 0].set_xlabel('Year')
    axes[1, 0].set_title('Disease Spread')
    
    axes[1, 1].set_ylabel('Sites occupied (%)')
    axes[1, 1].set_xlabel('Year')
    axes[1, 1].set_title('Range Occupancy')
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "fig1_trajectory.png", dpi=150)
    plt.close()
    print(f"  Saved: {OUT_DIR / 'fig1_trajectory.png'}")


def fig2_outplanting_intensity():
    """Figure 2: Effect of outplanting intensity."""
    print("Figure 2: Outplanting intensity sweep...")
    
    intensities = [0, 100, 500, 1000, 2000, 5000, 10000]
    n_sites_target = [10, 25, 50, 100]
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    for n_target in n_sites_target:
        final_pops = []
        final_resist = []
        
        for intensity in intensities:
            if intensity == 0:
                sites = []
            else:
                sites = list(range(0, min(n_target * 10, 500), 10))[:n_target]
            
            cfg = NetworkConfig(
                n_sites=500, n_years=100,
                outplanting_n=intensity,
                outplanting_sites=sites,
                outplanting_start=15,
            )
            
            # Average over 3 replicates
            pops, resists = [], []
            for seed in range(3):
                sim = NetworkSimulation(cfg, seed=seed)
                result = sim.run()
                pops.append(result.final_n_ratio * 100)
                resists.append(result.states[-1].mean_resistance_freq * 100)
            
            final_pops.append(np.mean(pops))
            final_resist.append(np.mean(resists))
        
        axes[0].plot(intensities, final_pops, 'o-', label=f'{n_target} sites', linewidth=2, markersize=6)
        axes[1].plot(intensities, final_resist, 'o-', label=f'{n_target} sites', linewidth=2, markersize=6)
    
    axes[0].set_xlabel('Outplanting intensity (individuals/year)')
    axes[0].set_ylabel('Final population (% of baseline)')
    axes[0].set_title('Recovery vs Outplanting Intensity')
    axes[0].legend()
    axes[0].set_xscale('symlog', linthresh=100)
    
    axes[1].set_xlabel('Outplanting intensity (individuals/year)')
    axes[1].set_ylabel('Final resistance (%)')
    axes[1].set_title('Resistance Evolution vs Intensity')
    axes[1].axhline(y=95, color='gray', linestyle='--', alpha=0.5)
    axes[1].legend()
    axes[1].set_xscale('symlog', linthresh=100)
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "fig2_intensity.png", dpi=150)
    plt.close()
    print(f"  Saved: {OUT_DIR / 'fig2_intensity.png'}")


def fig3_enhanced_vs_wild():
    """Figure 3: Enhanced resistance levels vs wild-caught."""
    print("Figure 3: Enhanced vs wild comparison...")
    
    enhanced_levels = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95]
    
    wild_pops = []
    wild_resist = []
    enhanced_pops = []
    enhanced_resist = []
    
    # Wild baseline
    cfg_wild = NetworkConfig(
        n_sites=500, n_years=100,
        outplanting_n=5000,
        outplanting_sites=list(range(0, 500, 10)),
        outplanting_start=15,
        outplant_resistance_mode='wild',
    )
    
    for seed in range(5):
        sim = NetworkSimulation(cfg_wild, seed=seed)
        result = sim.run()
        wild_pops.append(result.final_n_ratio * 100)
        wild_resist.append(result.states[-1].mean_resistance_freq * 100)
    
    wild_mean_pop = np.mean(wild_pops)
    wild_mean_resist = np.mean(wild_resist)
    
    # Enhanced at different levels
    for level in enhanced_levels:
        cfg = NetworkConfig(
            n_sites=500, n_years=100,
            outplanting_n=5000,
            outplanting_sites=list(range(0, 500, 10)),
            outplanting_start=15,
            outplant_resistance_mode='enhanced',
            outplant_enhanced_resistance=level,
        )
        
        pops, resists = [], []
        for seed in range(5):
            sim = NetworkSimulation(cfg, seed=seed)
            result = sim.run()
            pops.append(result.final_n_ratio * 100)
            resists.append(result.states[-1].mean_resistance_freq * 100)
        
        enhanced_pops.append(np.mean(pops))
        enhanced_resist.append(np.mean(resists))
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    enhanced_pct = [x * 100 for x in enhanced_levels]
    
    axes[0].plot(enhanced_pct, enhanced_pops, 'o-', color='#1f77b4', linewidth=2, markersize=8, label='Enhanced')
    axes[0].axhline(y=wild_mean_pop, color='#2ca02c', linestyle='--', linewidth=2, label=f'Wild ({wild_mean_pop:.1f}%)')
    axes[0].fill_between([10, 95], 0, wild_mean_pop, alpha=0.1, color='red')
    axes[0].fill_between([10, 95], wild_mean_pop, 20, alpha=0.1, color='green')
    axes[0].set_xlabel('Enhanced resistance level (%)')
    axes[0].set_ylabel('Final population (% of baseline)')
    axes[0].set_title('Population Recovery: Enhanced vs Wild')
    axes[0].legend()
    axes[0].text(50, wild_mean_pop - 2, 'Underperforms wild', ha='center', fontsize=10, color='red')
    axes[0].text(85, wild_mean_pop + 1, 'Matches/exceeds wild', ha='center', fontsize=10, color='green')
    
    axes[1].plot(enhanced_pct, enhanced_resist, 'o-', color='#1f77b4', linewidth=2, markersize=8, label='Enhanced')
    axes[1].axhline(y=wild_mean_resist, color='#2ca02c', linestyle='--', linewidth=2, label=f'Wild ({wild_mean_resist:.1f}%)')
    axes[1].axhline(y=95, color='gray', linestyle=':', alpha=0.5, label='95% ceiling')
    axes[1].plot([enhanced_pct[i] for i in range(len(enhanced_pct))], 
                 [enhanced_pct[i] for i in range(len(enhanced_pct))], 
                 'k:', alpha=0.3, label='1:1 line')
    axes[1].set_xlabel('Enhanced resistance level (%)')
    axes[1].set_ylabel('Final resistance (%)')
    axes[1].set_title('Resistance: Input vs Output')
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "fig3_enhanced_vs_wild.png", dpi=150)
    plt.close()
    print(f"  Saved: {OUT_DIR / 'fig3_enhanced_vs_wild.png'}")


def fig4_refugia_sensitivity():
    """Figure 4: Sensitivity to refugia fraction."""
    print("Figure 4: Refugia sensitivity...")
    
    refugia_fractions = [0.0, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20]
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    for intervention in [False, True]:
        final_pops = []
        extinctions = []
        
        for frac in refugia_fractions:
            if intervention:
                cfg = NetworkConfig(
                    n_sites=500, n_years=100,
                    refugia_fraction=frac,
                    outplanting_n=5000,
                    outplanting_sites=list(range(0, 500, 10)),
                    outplanting_start=15,
                )
            else:
                cfg = NetworkConfig(
                    n_sites=500, n_years=100,
                    refugia_fraction=frac,
                )
            
            pops = []
            ext_count = 0
            for seed in range(10):
                sim = NetworkSimulation(cfg, seed=seed)
                result = sim.run()
                pops.append(result.final_n_ratio * 100)
                if result.extinct:
                    ext_count += 1
            
            final_pops.append(np.mean(pops))
            extinctions.append(ext_count / 10 * 100)
        
        label = "With outplanting" if intervention else "No intervention"
        color = '#2ca02c' if intervention else '#d62728'
        
        axes[0].plot([f * 100 for f in refugia_fractions], final_pops, 'o-', 
                     label=label, color=color, linewidth=2, markersize=8)
        axes[1].plot([f * 100 for f in refugia_fractions], extinctions, 'o-', 
                     label=label, color=color, linewidth=2, markersize=8)
    
    axes[0].set_xlabel('Refugia fraction (%)')
    axes[0].set_ylabel('Final population (% of baseline)')
    axes[0].set_title('Recovery vs Refugia Availability')
    axes[0].legend()
    axes[0].set_yscale('log')
    axes[0].set_ylim(0.001, 30)
    
    axes[1].set_xlabel('Refugia fraction (%)')
    axes[1].set_ylabel('Extinction probability (%)')
    axes[1].set_title('Extinction Risk vs Refugia')
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "fig4_refugia.png", dpi=150)
    plt.close()
    print(f"  Saved: {OUT_DIR / 'fig4_refugia.png'}")


def fig5_mortality_sensitivity():
    """Figure 5: Sensitivity to disease mortality rate."""
    print("Figure 5: Mortality sensitivity...")
    
    mortality_rates = [0.50, 0.70, 0.80, 0.90, 0.95, 0.99]
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    for intervention in [False, True]:
        final_pops = []
        final_resist = []
        
        for mort in mortality_rates:
            if intervention:
                cfg = NetworkConfig(
                    n_sites=500, n_years=100,
                    disease_mortality=mort,
                    outplanting_n=5000,
                    outplanting_sites=list(range(0, 500, 10)),
                    outplanting_start=15,
                )
            else:
                cfg = NetworkConfig(
                    n_sites=500, n_years=100,
                    disease_mortality=mort,
                )
            
            pops, resists = [], []
            for seed in range(5):
                sim = NetworkSimulation(cfg, seed=seed)
                result = sim.run()
                pops.append(result.final_n_ratio * 100)
                if not result.extinct:
                    resists.append(result.states[-1].mean_resistance_freq * 100)
            
            final_pops.append(np.mean(pops))
            final_resist.append(np.mean(resists) if resists else 0)
        
        label = "With outplanting" if intervention else "No intervention"
        color = '#2ca02c' if intervention else '#d62728'
        
        axes[0].plot([m * 100 for m in mortality_rates], final_pops, 'o-', 
                     label=label, color=color, linewidth=2, markersize=8)
        axes[1].plot([m * 100 for m in mortality_rates], final_resist, 'o-', 
                     label=label, color=color, linewidth=2, markersize=8)
    
    axes[0].axvline(x=99, color='gray', linestyle='--', alpha=0.5, label='Calibrated (99%)')
    axes[0].set_xlabel('Disease mortality (%)')
    axes[0].set_ylabel('Final population (% of baseline)')
    axes[0].set_title('Recovery vs Disease Severity')
    axes[0].legend()
    axes[0].set_yscale('log')
    axes[0].set_ylim(0.001, 50)
    
    axes[1].axvline(x=99, color='gray', linestyle='--', alpha=0.5)
    axes[1].axhline(y=95, color='gray', linestyle=':', alpha=0.5, label='95% ceiling')
    axes[1].set_xlabel('Disease mortality (%)')
    axes[1].set_ylabel('Final resistance (%)')
    axes[1].set_title('Selection Strength vs Mortality')
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "fig5_mortality.png", dpi=150)
    plt.close()
    print(f"  Saved: {OUT_DIR / 'fig5_mortality.png'}")


def fig6_connectivity_comparison():
    """Figure 6: Effect of connectivity type."""
    print("Figure 6: Connectivity comparison...")
    
    conn_types = [
        ConnectivityType.UNIFORM,
        ConnectivityType.STEPPING_STONE,
        ConnectivityType.DISTANCE_DECAY,
        ConnectivityType.ASYMMETRIC_FLOW,
        ConnectivityType.HUB_NETWORK,
        ConnectivityType.MODULAR,
    ]
    
    results = []
    
    for conn in conn_types:
        cfg = NetworkConfig(
            n_sites=500, n_years=100,
            connectivity_type=conn,
            outplanting_n=5000,
            outplanting_sites=list(range(0, 500, 10)),
            outplanting_start=15,
        )
        
        pops, resists, occ = [], [], []
        for seed in range(5):
            sim = NetworkSimulation(cfg, seed=seed)
            result = sim.run()
            pops.append(result.final_n_ratio * 100)
            resists.append(result.states[-1].mean_resistance_freq * 100)
            occ.append(result.states[-1].occupied_ratio * 100)
        
        results.append({
            'type': conn.value,
            'pop_mean': np.mean(pops),
            'pop_std': np.std(pops),
            'resist': np.mean(resists),
            'occupied': np.mean(occ),
        })
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    types = [r['type'] for r in results]
    x = np.arange(len(types))
    
    axes[0].bar(x, [r['pop_mean'] for r in results], yerr=[r['pop_std'] for r in results],
                capsize=5, color='steelblue', edgecolor='black')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(types, rotation=45, ha='right')
    axes[0].set_ylabel('Final population (%)')
    axes[0].set_title('Recovery by Connectivity Type')
    
    axes[1].bar(x, [r['resist'] for r in results], color='forestgreen', edgecolor='black')
    axes[1].axhline(y=95, color='gray', linestyle='--', alpha=0.5)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(types, rotation=45, ha='right')
    axes[1].set_ylabel('Final resistance (%)')
    axes[1].set_title('Resistance by Connectivity')
    
    axes[2].bar(x, [r['occupied'] for r in results], color='coral', edgecolor='black')
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(types, rotation=45, ha='right')
    axes[2].set_ylabel('Sites occupied (%)')
    axes[2].set_title('Range Occupancy by Connectivity')
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "fig6_connectivity.png", dpi=150)
    plt.close()
    print(f"  Saved: {OUT_DIR / 'fig6_connectivity.png'}")


def fig7_timing_sensitivity():
    """Figure 7: Effect of outplanting timing."""
    print("Figure 7: Timing sensitivity...")
    
    start_years = [11, 15, 20, 25, 30, 40, 50]  # Years after year 0
    
    final_pops = []
    final_resist = []
    
    for start in start_years:
        cfg = NetworkConfig(
            n_sites=500, n_years=100,
            outplanting_n=5000,
            outplanting_sites=list(range(0, 500, 10)),
            outplanting_start=start,
        )
        
        pops, resists = [], []
        for seed in range(5):
            sim = NetworkSimulation(cfg, seed=seed)
            result = sim.run()
            pops.append(result.final_n_ratio * 100)
            resists.append(result.states[-1].mean_resistance_freq * 100)
        
        final_pops.append(np.mean(pops))
        final_resist.append(np.mean(resists))
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Convert to years after disease (disease at year 10)
    years_after = [s - 10 for s in start_years]
    
    axes[0].plot(years_after, final_pops, 'o-', color='steelblue', linewidth=2, markersize=8)
    axes[0].axvline(x=5, color='gray', linestyle='--', alpha=0.5, label='Baseline (5 yr)')
    axes[0].set_xlabel('Outplanting start (years after disease)')
    axes[0].set_ylabel('Final population (%)')
    axes[0].set_title('Recovery vs Intervention Timing')
    axes[0].legend()
    
    axes[1].plot(years_after, final_resist, 'o-', color='forestgreen', linewidth=2, markersize=8)
    axes[1].axhline(y=95, color='gray', linestyle=':', alpha=0.5)
    axes[1].set_xlabel('Outplanting start (years after disease)')
    axes[1].set_ylabel('Final resistance (%)')
    axes[1].set_title('Resistance vs Timing')
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "fig7_timing.png", dpi=150)
    plt.close()
    print(f"  Saved: {OUT_DIR / 'fig7_timing.png'}")


def fig8_stochasticity():
    """Figure 8: Stochastic variability across replicates."""
    print("Figure 8: Stochastic variability...")
    
    cfg = NetworkConfig(
        n_sites=500, n_years=100,
        outplanting_n=5000,
        outplanting_sites=list(range(0, 500, 10)),
        outplanting_start=15,
    )
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    all_trajectories = []
    final_pops = []
    
    for seed in range(20):
        sim = NetworkSimulation(cfg, seed=seed)
        result = sim.run()
        
        years = [s.year for s in result.states]
        n_ratio = [s.n_ratio * 100 for s in result.states]
        
        all_trajectories.append(n_ratio)
        final_pops.append(result.final_n_ratio * 100)
        
        axes[0].plot(years, n_ratio, alpha=0.3, color='steelblue')
    
    # Mean trajectory
    mean_traj = np.mean(all_trajectories, axis=0)
    axes[0].plot(years, mean_traj, color='darkblue', linewidth=3, label='Mean')
    axes[0].set_xlabel('Year')
    axes[0].set_ylabel('Population (% of baseline)')
    axes[0].set_title('Trajectory Variability (20 replicates)')
    axes[0].set_yscale('log')
    axes[0].set_ylim(0.001, 100)
    axes[0].legend()
    
    axes[1].hist(final_pops, bins=10, color='steelblue', edgecolor='black', alpha=0.7)
    axes[1].axvline(x=np.mean(final_pops), color='red', linestyle='--', 
                    label=f'Mean: {np.mean(final_pops):.1f}%')
    axes[1].axvline(x=np.median(final_pops), color='orange', linestyle=':', 
                    label=f'Median: {np.median(final_pops):.1f}%')
    axes[1].set_xlabel('Final population (%)')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Distribution of Outcomes')
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "fig8_stochasticity.png", dpi=150)
    plt.close()
    print(f"  Saved: {OUT_DIR / 'fig8_stochasticity.png'}")


def fig9_heatmap_intensity_sites():
    """Figure 9: Heatmap of intensity × sites."""
    print("Figure 9: Intensity × sites heatmap...")
    
    intensities = [100, 500, 1000, 2000, 5000, 10000]
    n_sites_list = [5, 10, 25, 50, 100]
    
    pop_matrix = np.zeros((len(n_sites_list), len(intensities)))
    
    for i, n_target in enumerate(n_sites_list):
        for j, intensity in enumerate(intensities):
            sites = list(range(0, min(n_target * 10, 500), 10))[:n_target]
            
            cfg = NetworkConfig(
                n_sites=500, n_years=100,
                outplanting_n=intensity,
                outplanting_sites=sites,
                outplanting_start=15,
            )
            
            pops = []
            for seed in range(3):
                sim = NetworkSimulation(cfg, seed=seed)
                result = sim.run()
                pops.append(result.final_n_ratio * 100)
            
            pop_matrix[i, j] = np.mean(pops)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    im = ax.imshow(pop_matrix, cmap='YlGn', aspect='auto')
    
    ax.set_xticks(np.arange(len(intensities)))
    ax.set_yticks(np.arange(len(n_sites_list)))
    ax.set_xticklabels(intensities)
    ax.set_yticklabels(n_sites_list)
    
    ax.set_xlabel('Outplanting intensity (individuals/year)')
    ax.set_ylabel('Number of target sites')
    ax.set_title('Final Population (%) by Intensity × Sites')
    
    # Add text annotations
    for i in range(len(n_sites_list)):
        for j in range(len(intensities)):
            text = ax.text(j, i, f'{pop_matrix[i, j]:.1f}',
                          ha='center', va='center', color='black', fontsize=10)
    
    cbar = plt.colorbar(im)
    cbar.set_label('Final population (%)')
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "fig9_heatmap.png", dpi=150)
    plt.close()
    print(f"  Saved: {OUT_DIR / 'fig9_heatmap.png'}")


def fig10_connectivity_matrix():
    """Figure 10: Visualization of connectivity matrices."""
    print("Figure 10: Connectivity matrices...")
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    conn_types = [
        ConnectivityType.UNIFORM,
        ConnectivityType.STEPPING_STONE,
        ConnectivityType.DISTANCE_DECAY,
        ConnectivityType.ASYMMETRIC_FLOW,
        ConnectivityType.HUB_NETWORK,
        ConnectivityType.MODULAR,
    ]
    
    for ax, conn in zip(axes, conn_types):
        cfg = NetworkConfig(n_sites=50, connectivity_type=conn)
        C = build_connectivity_matrix(cfg)
        
        im = ax.imshow(C, cmap='viridis', aspect='auto')
        ax.set_title(conn.value.replace('_', ' ').title())
        ax.set_xlabel('Destination site')
        ax.set_ylabel('Source site')
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    
    plt.suptitle('Connectivity Matrices (50 sites)', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "fig10_connectivity_matrices.png", dpi=150)
    plt.close()
    print(f"  Saved: {OUT_DIR / 'fig10_connectivity_matrices.png'}")


def main():
    print("=" * 60)
    print("EXHAUSTIVE MODEL INTERROGATION")
    print("=" * 60)
    print()
    
    fig1_single_trajectory()
    fig2_outplanting_intensity()
    fig3_enhanced_vs_wild()
    fig4_refugia_sensitivity()
    fig5_mortality_sensitivity()
    fig6_connectivity_comparison()
    fig7_timing_sensitivity()
    fig8_stochasticity()
    fig9_heatmap_intensity_sites()
    fig10_connectivity_matrix()
    
    print()
    print("=" * 60)
    print(f"All figures saved to: {OUT_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
