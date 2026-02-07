#!/usr/bin/env python3
"""Quick visualization of key model dynamics."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from pycnopodia.network import NetworkSimulation, NetworkConfig

OUT_DIR = Path(__file__).parent.parent / "figures"
OUT_DIR.mkdir(exist_ok=True)

plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11

def main():
    print("Generating figures...")
    
    # Figure 1: Trajectory comparison
    print("  Fig 1: Trajectories...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    scenarios = [
        ("No intervention", NetworkConfig(n_sites=200, n_years=80)),
        ("5k/yr wild", NetworkConfig(n_sites=200, n_years=80,
            outplanting_n=5000, outplanting_sites=list(range(0, 200, 10)),
            outplanting_start=15)),
        ("5k/yr enhanced 95%", NetworkConfig(n_sites=200, n_years=80,
            outplanting_n=5000, outplanting_sites=list(range(0, 200, 10)),
            outplanting_start=15, outplant_resistance_mode='enhanced',
            outplant_enhanced_resistance=0.95)),
    ]
    
    colors = ['#d62728', '#2ca02c', '#1f77b4']
    
    for (name, cfg), color in zip(scenarios, colors):
        sim = NetworkSimulation(cfg, seed=42)
        result = sim.run()
        years = [s.year for s in result.states]
        
        axes[0,0].plot(years, [s.n_ratio*100 for s in result.states], label=name, color=color, lw=2)
        axes[0,1].plot(years, [s.mean_resistance_freq*100 for s in result.states], label=name, color=color, lw=2)
        axes[1,0].plot(years, [s.infected_sites_ratio*100 for s in result.states], label=name, color=color, lw=2)
        axes[1,1].plot(years, [s.occupied_ratio*100 for s in result.states], label=name, color=color, lw=2)
    
    axes[0,0].set_ylabel('Population (% baseline)'); axes[0,0].set_title('Population'); axes[0,0].legend()
    axes[0,0].set_yscale('log'); axes[0,0].set_ylim(0.001, 100)
    axes[0,0].axvline(10, color='gray', ls='--', alpha=0.5)
    axes[0,1].set_ylabel('Resistance (%)'); axes[0,1].set_title('Resistance Evolution')
    axes[0,1].axhline(95, color='gray', ls='--', alpha=0.5); axes[0,1].legend()
    axes[1,0].set_ylabel('Sites infected (%)'); axes[1,0].set_xlabel('Year'); axes[1,0].set_title('Disease')
    axes[1,1].set_ylabel('Sites occupied (%)'); axes[1,1].set_xlabel('Year'); axes[1,1].set_title('Range')
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "fig1_trajectories.png", dpi=150)
    plt.close()
    
    # Figure 2: Enhanced vs Wild
    print("  Fig 2: Enhanced vs Wild...")
    enhanced_levels = [0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95]
    
    # Wild baseline
    cfg = NetworkConfig(n_sites=200, n_years=80,
        outplanting_n=5000, outplanting_sites=list(range(0, 200, 10)),
        outplanting_start=15)
    wild_pops = [NetworkSimulation(cfg, seed=s).run().final_n_ratio*100 for s in range(3)]
    wild_mean = np.mean(wild_pops)
    
    enhanced_pops = []
    for level in enhanced_levels:
        cfg = NetworkConfig(n_sites=200, n_years=80,
            outplanting_n=5000, outplanting_sites=list(range(0, 200, 10)),
            outplanting_start=15, outplant_resistance_mode='enhanced',
            outplant_enhanced_resistance=level)
        pops = [NetworkSimulation(cfg, seed=s).run().final_n_ratio*100 for s in range(3)]
        enhanced_pops.append(np.mean(pops))
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot([l*100 for l in enhanced_levels], enhanced_pops, 'o-', color='#1f77b4', lw=2, ms=8, label='Enhanced')
    ax.axhline(wild_mean, color='#2ca02c', ls='--', lw=2, label=f'Wild ({wild_mean:.1f}%)')
    ax.fill_between([10, 95], 0, wild_mean, alpha=0.1, color='red')
    ax.fill_between([10, 95], wild_mean, max(enhanced_pops)+2, alpha=0.1, color='green')
    ax.set_xlabel('Enhanced resistance level (%)')
    ax.set_ylabel('Final population (% baseline)')
    ax.set_title('Enhanced Breeding vs Wild-Caught Outplanting')
    ax.legend()
    ax.text(40, wild_mean-1.5, 'Underperforms wild', ha='center', color='red')
    ax.text(87, wild_mean+0.8, 'Matches wild', ha='center', color='green')
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "fig2_enhanced_vs_wild.png", dpi=150)
    plt.close()
    
    # Figure 3: Refugia sensitivity
    print("  Fig 3: Refugia...")
    refugia = [0.0, 0.01, 0.02, 0.05, 0.10, 0.15]
    
    pops_no_int = []
    pops_with_int = []
    for r in refugia:
        cfg1 = NetworkConfig(n_sites=200, n_years=80, refugia_fraction=r)
        cfg2 = NetworkConfig(n_sites=200, n_years=80, refugia_fraction=r,
            outplanting_n=5000, outplanting_sites=list(range(0, 200, 10)), outplanting_start=15)
        pops_no_int.append(np.mean([NetworkSimulation(cfg1, seed=s).run().final_n_ratio*100 for s in range(3)]))
        pops_with_int.append(np.mean([NetworkSimulation(cfg2, seed=s).run().final_n_ratio*100 for s in range(3)]))
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot([r*100 for r in refugia], pops_no_int, 'o-', color='#d62728', lw=2, ms=8, label='No intervention')
    ax.plot([r*100 for r in refugia], pops_with_int, 'o-', color='#2ca02c', lw=2, ms=8, label='With outplanting')
    ax.set_xlabel('Refugia fraction (%)')
    ax.set_ylabel('Final population (% baseline)')
    ax.set_title('Recovery vs Refugia Availability')
    ax.legend()
    ax.set_yscale('log')
    ax.set_ylim(0.001, 30)
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "fig3_refugia.png", dpi=150)
    plt.close()
    
    # Figure 4: Intensity heatmap
    print("  Fig 4: Intensity heatmap...")
    intensities = [100, 500, 1000, 2000, 5000]
    n_sites_list = [5, 10, 25, 50]
    
    matrix = np.zeros((len(n_sites_list), len(intensities)))
    for i, ns in enumerate(n_sites_list):
        for j, intensity in enumerate(intensities):
            sites = list(range(0, min(ns*10, 200), 10))[:ns]
            cfg = NetworkConfig(n_sites=200, n_years=80,
                outplanting_n=intensity, outplanting_sites=sites, outplanting_start=15)
            matrix[i,j] = np.mean([NetworkSimulation(cfg, seed=s).run().final_n_ratio*100 for s in range(2)])
    
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(matrix, cmap='YlGn', aspect='auto')
    ax.set_xticks(range(len(intensities))); ax.set_xticklabels(intensities)
    ax.set_yticks(range(len(n_sites_list))); ax.set_yticklabels(n_sites_list)
    ax.set_xlabel('Outplanting intensity (ind/yr)')
    ax.set_ylabel('Number of sites')
    ax.set_title('Final Population (%) by Intensity × Sites')
    for i in range(len(n_sites_list)):
        for j in range(len(intensities)):
            ax.text(j, i, f'{matrix[i,j]:.1f}', ha='center', va='center')
    plt.colorbar(im, label='Final pop (%)')
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "fig4_heatmap.png", dpi=150)
    plt.close()
    
    # Figure 5: Stochasticity
    print("  Fig 5: Stochasticity...")
    cfg = NetworkConfig(n_sites=200, n_years=80,
        outplanting_n=5000, outplanting_sites=list(range(0, 200, 10)), outplanting_start=15)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    all_traj = []
    final_pops = []
    for seed in range(15):
        result = NetworkSimulation(cfg, seed=seed).run()
        years = [s.year for s in result.states]
        traj = [s.n_ratio*100 for s in result.states]
        all_traj.append(traj)
        final_pops.append(result.final_n_ratio*100)
        axes[0].plot(years, traj, alpha=0.3, color='steelblue')
    
    axes[0].plot(years, np.mean(all_traj, axis=0), color='darkblue', lw=3, label='Mean')
    axes[0].set_xlabel('Year'); axes[0].set_ylabel('Population (%)'); axes[0].set_title('15 Replicates')
    axes[0].set_yscale('log'); axes[0].set_ylim(0.001, 100); axes[0].legend()
    
    axes[1].hist(final_pops, bins=8, color='steelblue', edgecolor='black')
    axes[1].axvline(np.mean(final_pops), color='red', ls='--', label=f'Mean: {np.mean(final_pops):.1f}%')
    axes[1].set_xlabel('Final pop (%)'); axes[1].set_ylabel('Count'); axes[1].set_title('Outcome Distribution')
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "fig5_stochasticity.png", dpi=150)
    plt.close()
    
    # Figure 6: Mortality sensitivity  
    print("  Fig 6: Mortality...")
    mortalities = [0.7, 0.8, 0.9, 0.95, 0.99]
    
    pops_no = []
    pops_with = []
    for m in mortalities:
        cfg1 = NetworkConfig(n_sites=200, n_years=80, disease_mortality=m)
        cfg2 = NetworkConfig(n_sites=200, n_years=80, disease_mortality=m,
            outplanting_n=5000, outplanting_sites=list(range(0, 200, 10)), outplanting_start=15)
        pops_no.append(np.mean([NetworkSimulation(cfg1, seed=s).run().final_n_ratio*100 for s in range(3)]))
        pops_with.append(np.mean([NetworkSimulation(cfg2, seed=s).run().final_n_ratio*100 for s in range(3)]))
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot([m*100 for m in mortalities], pops_no, 'o-', color='#d62728', lw=2, ms=8, label='No intervention')
    ax.plot([m*100 for m in mortalities], pops_with, 'o-', color='#2ca02c', lw=2, ms=8, label='With outplanting')
    ax.axvline(99, color='gray', ls='--', alpha=0.5, label='Calibrated (99%)')
    ax.set_xlabel('Disease mortality (%)')
    ax.set_ylabel('Final population (%)')
    ax.set_title('Recovery vs Disease Severity')
    ax.legend()
    ax.set_yscale('log')
    ax.set_ylim(0.001, 50)
    
    plt.tight_layout()
    plt.savefig(OUT_DIR / "fig6_mortality.png", dpi=150)
    plt.close()
    
    print(f"\nAll figures saved to: {OUT_DIR}")
    print("Done!")

if __name__ == "__main__":
    main()
