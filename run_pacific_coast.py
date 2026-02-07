#!/usr/bin/env python3
"""
Run Pacific Coast Pycnopodia simulation with extensive visualizations.

Generates figures to validate:
- Geographic structure (8 regions, temperature gradient)
- Connectivity patterns (larval vs disease, asymmetric currents)
- Disease spread (starting in Southern CA, moving north)
- Population trajectories by region
- Resistance evolution
- BC Fjords as refugia

All figures saved to figures/pacific_coast/
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from pycnopodia.pacific_coast import (
    PacificCoastConfig,
    PacificCoastSimulation,
    PacificCoastResult,
    PACIFIC_COAST_REGIONS,
    REGION_ORDER,
    RegionType,
    get_site_region_boundaries,
    get_temperature_disease_modifier,
)


# Create output directory
FIGURE_DIR = Path(__file__).parent / "figures" / "pacific_coast"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)


def save_figure(fig, name: str, dpi: int = 150):
    """Save figure with consistent settings."""
    path = FIGURE_DIR / f"{name}.png"
    fig.savefig(path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: {path}")


def plot_geographic_structure(result: PacificCoastResult):
    """
    Plot map-like visualization of the 8 regions.
    
    Shows sites as 1D coastline (south to north on y-axis) with:
    - Initial population density
    - Temperature gradient
    - Region boundaries labeled
    """
    print("\n📍 Generating geographic structure plots...")
    
    sites = result.sites
    config = result.config
    boundaries = get_site_region_boundaries(sites)
    
    n_sites = len(sites)
    latitudes = np.array([s.latitude for s in sites])
    temperatures = np.array([s.temperature for s in sites])
    initial_pops = np.array([s.initial_population for s in sites])
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 10))
    
    # --- Panel 1: Initial Population by Region ---
    ax1 = axes[0]
    for region_id in REGION_ORDER:
        start, end = boundaries[region_id]
        region = config.regions[region_id]
        ax1.barh(
            range(start, end), 
            initial_pops[start:end],
            color=region.color,
            label=region.short_name,
            height=1.0,
        )
    
    ax1.set_xlabel("Initial Relative Density", fontsize=12)
    ax1.set_ylabel("Site Index (South → North)", fontsize=12)
    ax1.set_title("Initial Population Density by Region", fontsize=14)
    ax1.legend(loc='upper right', fontsize=9)
    
    # Add region boundary lines
    for region_id in REGION_ORDER:
        start, end = boundaries[region_id]
        ax1.axhline(y=start, color='black', linewidth=0.5, linestyle='--', alpha=0.5)
    
    # --- Panel 2: Temperature Gradient ---
    ax2 = axes[1]
    colors = plt.cm.coolwarm(np.linspace(0, 1, len(temperatures)))
    for i, (lat, temp) in enumerate(zip(range(n_sites), temperatures)):
        ax2.barh(i, temp, color=plt.cm.coolwarm((temp - 7) / 10), height=1.0)
    
    ax2.set_xlabel("Temperature (°C)", fontsize=12)
    ax2.set_ylabel("Site Index (South → North)", fontsize=12)
    ax2.set_title("Temperature Gradient", fontsize=14)
    ax2.axvline(x=12, color='red', linestyle='--', linewidth=2, label='Disease threshold')
    ax2.legend()
    
    # Add region labels on right
    for region_id in REGION_ORDER:
        start, end = boundaries[region_id]
        mid = (start + end) / 2
        region = config.regions[region_id]
        ax2.text(
            ax2.get_xlim()[1] + 0.3, mid, 
            region.short_name, 
            fontsize=9, 
            va='center',
            color=region.color,
            fontweight='bold'
        )
    
    # --- Panel 3: Region Type and Refugia ---
    ax3 = axes[2]
    
    type_colors = {
        RegionType.OUTER_COAST: '#4ECDC4',
        RegionType.INLAND_SEA: '#FF6B6B',
        RegionType.FJORD: '#45B7D1',
    }
    
    for i, site in enumerate(sites):
        region = config.regions[site.region_id]
        color = type_colors[region.region_type]
        if site.is_refugia:
            ax3.barh(i, 1, color='gold', height=1.0)
        else:
            ax3.barh(i, 1, color=color, height=1.0, alpha=0.7)
    
    ax3.set_xlabel("Region Type", fontsize=12)
    ax3.set_ylabel("Site Index (South → North)", fontsize=12)
    ax3.set_title("Region Types & Refugia", fontsize=14)
    ax3.set_xticks([])
    
    # Legend for region types
    patches = [
        mpatches.Patch(color='#4ECDC4', label='Outer Coast'),
        mpatches.Patch(color='#FF6B6B', label='Inland Sea'),
        mpatches.Patch(color='gold', label='Fjord (Refugia)'),
    ]
    ax3.legend(handles=patches, loc='upper right')
    
    plt.tight_layout()
    save_figure(fig, "01_geographic_structure")
    
    # --- Additional: Latitude profile ---
    fig2, ax = plt.subplots(figsize=(10, 8))
    
    for region_id in REGION_ORDER:
        start, end = boundaries[region_id]
        region = config.regions[region_id]
        lats = latitudes[start:end]
        ax.scatter(
            [region.short_name] * len(lats), 
            lats,
            c=[region.color] * len(lats),
            s=30,
            alpha=0.6
        )
    
    ax.set_xlabel("Region", fontsize=12)
    ax.set_ylabel("Latitude (°N)", fontsize=12)
    ax.set_title("Site Distribution by Latitude", fontsize=14)
    plt.xticks(rotation=45, ha='right')
    
    plt.tight_layout()
    save_figure(fig2, "01b_latitude_distribution")


def plot_connectivity_matrices(result: PacificCoastResult):
    """
    Plot connectivity matrix heatmaps.
    
    Shows:
    - Larval connectivity (California Current pattern)
    - Disease connectivity
    - Side-by-side comparison
    """
    print("\n🔗 Generating connectivity matrix plots...")
    
    larval_C = result.larval_connectivity
    disease_D = result.disease_connectivity
    config = result.config
    boundaries = get_site_region_boundaries(result.sites)
    
    # --- Larval Connectivity ---
    fig1, ax = plt.subplots(figsize=(12, 10))
    
    im = ax.imshow(larval_C, cmap='Blues', aspect='auto', origin='lower')
    plt.colorbar(im, ax=ax, label='Connection Probability')
    
    # Add region boundary lines
    for region_id in REGION_ORDER:
        start, end = boundaries[region_id]
        ax.axhline(y=start, color='red', linewidth=0.5, alpha=0.5)
        ax.axvline(x=start, color='red', linewidth=0.5, alpha=0.5)
    
    # Add region labels
    for region_id in REGION_ORDER:
        start, end = boundaries[region_id]
        mid = (start + end) / 2
        region = config.regions[region_id]
        ax.text(-15, mid, region.short_name, fontsize=8, va='center', ha='right')
        ax.text(mid, -15, region.short_name, fontsize=8, va='top', ha='center', rotation=45)
    
    ax.set_xlabel("Target Site (receiving larvae)", fontsize=12)
    ax.set_ylabel("Source Site (sending larvae)", fontsize=12)
    ax.set_title("Larval Connectivity Matrix\n(California Current: North→South dominant)", fontsize=14)
    
    plt.tight_layout()
    save_figure(fig1, "02_larval_connectivity")
    
    # --- Disease Connectivity ---
    fig2, ax = plt.subplots(figsize=(12, 10))
    
    im = ax.imshow(disease_D, cmap='Reds', aspect='auto', origin='lower')
    plt.colorbar(im, ax=ax, label='Transmission Weight')
    
    for region_id in REGION_ORDER:
        start, end = boundaries[region_id]
        ax.axhline(y=start, color='black', linewidth=0.5, alpha=0.5)
        ax.axvline(x=start, color='black', linewidth=0.5, alpha=0.5)
    
    ax.set_xlabel("Target Site", fontsize=12)
    ax.set_ylabel("Source Site", fontsize=12)
    ax.set_title("Disease Connectivity Matrix\n(Temperature-modulated spread)", fontsize=14)
    
    plt.tight_layout()
    save_figure(fig2, "02_disease_connectivity")
    
    # --- Side-by-side comparison ---
    fig3, axes = plt.subplots(1, 2, figsize=(18, 8))
    
    im1 = axes[0].imshow(larval_C, cmap='Blues', aspect='auto', origin='lower')
    axes[0].set_title("Larval Connectivity\n(Oceanographic currents)", fontsize=14)
    plt.colorbar(im1, ax=axes[0], label='Probability', shrink=0.8)
    
    im2 = axes[1].imshow(disease_D, cmap='Reds', aspect='auto', origin='lower')
    axes[1].set_title("Disease Connectivity\n(Waterborne pathogen)", fontsize=14)
    plt.colorbar(im2, ax=axes[1], label='Weight', shrink=0.8)
    
    for ax in axes:
        for region_id in REGION_ORDER:
            start, _ = boundaries[region_id]
            ax.axhline(y=start, color='white', linewidth=0.5, alpha=0.8)
            ax.axvline(x=start, color='white', linewidth=0.5, alpha=0.8)
        ax.set_xlabel("Target Site", fontsize=11)
        ax.set_ylabel("Source Site", fontsize=11)
    
    plt.suptitle("Connectivity Comparison: Larvae vs Disease", fontsize=16, y=1.02)
    plt.tight_layout()
    save_figure(fig3, "02c_connectivity_comparison")


def plot_disease_spread_spacetime(result: PacificCoastResult):
    """
    Plot disease spread as space-time heatmap.
    
    x = site (south to north)
    y = year
    color = disease prevalence
    
    Should show wave starting in S. CA and moving north.
    BC Fjords should stay disease-free.
    """
    print("\n🦠 Generating disease spread space-time plot...")
    
    n_sites = len(result.sites)
    n_years = len(result.states)
    
    # Build space-time matrix
    disease_matrix = np.zeros((n_years, n_sites))
    for t, state in enumerate(result.states):
        disease_matrix[t, :] = state.disease_prevalence
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    im = ax.imshow(
        disease_matrix, 
        aspect='auto', 
        cmap='YlOrRd',
        origin='lower',
        extent=[0, n_sites, 0, n_years],
        vmin=0, vmax=1,
    )
    plt.colorbar(im, ax=ax, label='Disease Prevalence')
    
    # Add region boundaries
    boundaries = get_site_region_boundaries(result.sites)
    for region_id in REGION_ORDER:
        start, end = boundaries[region_id]
        ax.axvline(x=start, color='white', linewidth=1, linestyle='--', alpha=0.7)
        
        # Label at top
        mid = (start + end) / 2
        region = result.config.regions[region_id]
        ax.text(
            mid, n_years + 1, 
            region.short_name, 
            fontsize=9, 
            ha='center', 
            va='bottom',
            color=region.color,
            fontweight='bold'
        )
    
    # Mark disease onset
    onset = result.config.disease_onset_year
    ax.axhline(y=onset, color='red', linewidth=2, linestyle=':', label=f'Disease onset (year {onset})')
    
    ax.set_xlabel("Site Index (South → North)", fontsize=12)
    ax.set_ylabel("Year", fontsize=12)
    ax.set_title("Disease Spread Across Pacific Coast\n(Wave from Southern CA northward, Fjords protected)", fontsize=14)
    ax.legend(loc='lower right')
    
    plt.tight_layout()
    save_figure(fig, "03_disease_spread_spacetime")


def plot_population_trajectories_by_region(result: PacificCoastResult):
    """
    Plot population trajectories for each region.
    
    Shows:
    - Southern regions crash first and hardest
    - BC Fjords maintain population
    - Different recovery patterns
    """
    print("\n📊 Generating population trajectory plots...")
    
    trajectories = result.get_region_trajectories()
    years = result.years
    config = result.config
    
    # --- All regions on one plot ---
    fig1, ax = plt.subplots(figsize=(12, 8))
    
    for region_id in REGION_ORDER:
        region = config.regions[region_id]
        pop_ratio = trajectories[region_id]["population_ratio"]
        ax.plot(
            years, pop_ratio,
            label=region.short_name,
            color=region.color,
            linewidth=2,
        )
    
    ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Baseline')
    ax.axhline(y=0.1, color='gray', linestyle=':', alpha=0.5, label='10% threshold')
    ax.axvline(x=config.disease_onset_year, color='red', linestyle=':', alpha=0.5, label='Disease onset')
    
    ax.set_xlabel("Year", fontsize=12)
    ax.set_ylabel("Population Ratio (relative to initial)", fontsize=12)
    ax.set_title("Population Trajectories by Region", fontsize=14)
    ax.legend(loc='upper right', fontsize=9)
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 2.0)
    ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.3)
    ax.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    save_figure(fig1, "04_population_trajectories")
    
    # --- Separate subplots per region ---
    n_regions = len(REGION_ORDER)
    ncols = min(4, n_regions)
    nrows = (n_regions + ncols - 1) // ncols
    fig2, axes = plt.subplots(nrows, ncols, figsize=(16, 4 * nrows))
    axes = axes.flatten()
    
    for i, region_id in enumerate(REGION_ORDER):
        ax = axes[i]
        region = config.regions[region_id]
        
        pop_ratio = trajectories[region_id]["population_ratio"]
        disease = trajectories[region_id]["disease_prevalence"]
        resistance = trajectories[region_id]["resistance"]
        
        ax.plot(years, pop_ratio, color=region.color, linewidth=2, label='Population')
        ax.plot(years, disease, color='red', linewidth=1.5, linestyle='--', label='Disease', alpha=0.7)
        ax.fill_between(years, 0, pop_ratio, color=region.color, alpha=0.2)
        
        ax.set_title(f"{region.name}\n(Observed survival: {region.post_sswd_survival:.0%})", fontsize=10)
        ax.set_xlabel("Year", fontsize=9)
        ax.set_ylabel("Ratio", fontsize=9)
        ax.set_yscale('log')
        ax.set_ylim(1e-4, 2.0)
        ax.grid(True, alpha=0.3)
        
        if i == 0:
            ax.legend(fontsize=8)
    
    plt.suptitle("Regional Population & Disease Dynamics", fontsize=14, y=1.02)
    plt.tight_layout()
    save_figure(fig2, "04b_regional_trajectories_detail")


def plot_temperature_effects(result: PacificCoastResult):
    """
    Plot temperature effects on disease.
    
    Shows:
    - Temperature gradient across sites
    - Temperature-mortality relationship
    - Comparison with/without temperature effect
    """
    print("\n🌡️ Generating temperature effect plots...")
    
    sites = result.sites
    config = result.config
    n_sites = len(sites)
    
    temperatures = np.array([s.temperature for s in sites])
    
    # Get mortality after disease for each site
    final_state = result.states[-1]
    final_pop_ratios = final_state.populations / result.states[0].populations
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # --- Panel 1: Temperature gradient ---
    ax1 = axes[0, 0]
    ax1.scatter(range(n_sites), temperatures, c=temperatures, cmap='coolwarm', s=20)
    ax1.set_xlabel("Site Index (South → North)", fontsize=11)
    ax1.set_ylabel("Temperature (°C)", fontsize=11)
    ax1.set_title("Temperature Gradient Across Coast", fontsize=12)
    ax1.axhline(y=config.disease_temp_threshold, color='red', linestyle='--', 
                label=f'Threshold ({config.disease_temp_threshold}°C)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # --- Panel 2: Temperature vs final population ---
    ax2 = axes[0, 1]
    ax2.scatter(temperatures, final_pop_ratios, c=temperatures, cmap='coolwarm', s=30, alpha=0.6)
    ax2.set_xlabel("Site Temperature (°C)", fontsize=11)
    ax2.set_ylabel("Final Population Ratio", fontsize=11)
    ax2.set_title("Temperature vs Population Survival", fontsize=12)
    ax2.axvline(x=config.disease_temp_threshold, color='red', linestyle='--', alpha=0.7)
    ax2.grid(True, alpha=0.3)
    
    # --- Panel 3: Temperature-mortality modifier curve ---
    ax3 = axes[1, 0]
    temp_range = np.linspace(7, 18, 100)
    modifiers = [get_temperature_disease_modifier(t, config) for t in temp_range]
    ax3.plot(temp_range, modifiers, linewidth=2, color='red')
    ax3.fill_between(temp_range, 1, modifiers, alpha=0.2, color='red')
    ax3.set_xlabel("Temperature (°C)", fontsize=11)
    ax3.set_ylabel("Mortality Modifier", fontsize=11)
    ax3.set_title("Temperature-Dependent Mortality Modifier", fontsize=12)
    ax3.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)
    ax3.axvline(x=config.disease_temp_threshold, color='blue', linestyle=':', 
                label=f'Threshold ({config.disease_temp_threshold}°C)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # --- Panel 4: Region temperature vs survival scatter ---
    ax4 = axes[1, 1]
    
    for region_id in REGION_ORDER:
        region = config.regions[region_id]
        ax4.scatter(
            region.base_temperature,
            region.post_sswd_survival,
            s=200,
            c=[region.color],
            label=region.short_name,
            marker='o',
            edgecolors='black',
        )
    
    ax4.set_xlabel("Mean Region Temperature (°C)", fontsize=11)
    ax4.set_ylabel("Observed Post-SSWD Survival", fontsize=11)
    ax4.set_title("Temperature vs Observed Survival by Region", fontsize=12)
    ax4.legend(loc='upper right', fontsize=8)
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    save_figure(fig, "05_temperature_effects")


def plot_resistance_evolution(result: PacificCoastResult):
    """
    Plot resistance allele frequency evolution by region.
    
    Northern cold regions may evolve resistance slower
    (less selection pressure if less disease).
    """
    print("\n🧬 Generating resistance evolution plots...")
    
    trajectories = result.get_region_trajectories()
    years = result.years
    config = result.config
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # --- Panel 1: Resistance by region ---
    ax1 = axes[0]
    
    for region_id in REGION_ORDER:
        region = config.regions[region_id]
        resistance = trajectories[region_id]["resistance"]
        ax1.plot(
            years, resistance,
            label=region.short_name,
            color=region.color,
            linewidth=2,
        )
    
    ax1.axhline(y=config.initial_resistance_freq, color='gray', linestyle='--', 
                alpha=0.5, label='Initial')
    ax1.axhline(y=config.max_resistance_freq, color='gray', linestyle=':', 
                alpha=0.5, label='Ceiling')
    ax1.axvline(x=config.disease_onset_year, color='red', linestyle=':', 
                alpha=0.5, label='Disease onset')
    
    ax1.set_xlabel("Year", fontsize=12)
    ax1.set_ylabel("Mean Resistance Frequency", fontsize=12)
    ax1.set_title("Resistance Evolution by Region", fontsize=14)
    ax1.legend(loc='upper left', fontsize=8, ncol=2)
    ax1.set_ylim(0, 1.0)
    ax1.grid(True, alpha=0.3)
    
    # --- Panel 2: Resistance vs Temperature ---
    ax2 = axes[1]
    
    # Get final resistance by region
    final_resistance = {}
    for region_id in REGION_ORDER:
        final_resistance[region_id] = trajectories[region_id]["resistance"][-1]
    
    temps = [config.regions[r].base_temperature for r in REGION_ORDER]
    resists = [final_resistance[r] for r in REGION_ORDER]
    colors = [config.regions[r].color for r in REGION_ORDER]
    
    ax2.scatter(temps, resists, c=colors, s=200, edgecolors='black')
    
    for i, region_id in enumerate(REGION_ORDER):
        ax2.annotate(
            config.regions[region_id].short_name,
            (temps[i], resists[i]),
            xytext=(5, 5),
            textcoords='offset points',
            fontsize=9,
        )
    
    ax2.set_xlabel("Region Temperature (°C)", fontsize=12)
    ax2.set_ylabel("Final Resistance Frequency", fontsize=12)
    ax2.set_title("Final Resistance vs Temperature", fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    save_figure(fig, "06_resistance_evolution")
    
    # --- Detailed: resistance + disease + population overlay ---
    n_regions = len(REGION_ORDER)
    ncols = min(4, n_regions)
    nrows = (n_regions + ncols - 1) // ncols
    fig2, axes2 = plt.subplots(nrows, ncols, figsize=(16, 4 * nrows))
    axes2 = axes2.flatten()
    
    for i, region_id in enumerate(REGION_ORDER):
        ax = axes2[i]
        region = config.regions[region_id]
        
        resistance = trajectories[region_id]["resistance"]
        disease = trajectories[region_id]["disease_prevalence"]
        pop_ratio = trajectories[region_id]["population_ratio"]
        
        ax.plot(years, resistance, color='green', linewidth=2, label='Resistance')
        ax.plot(years, disease, color='red', linewidth=1.5, linestyle='--', label='Disease')
        ax.plot(years, pop_ratio, color=region.color, linewidth=1.5, linestyle=':', label='Population')
        
        ax.set_title(f"{region.short_name}", fontsize=11)
        ax.set_xlabel("Year", fontsize=9)
        ax.set_ylim(0, 1.0)
        ax.grid(True, alpha=0.3)
        
        if i == 0:
            ax.legend(fontsize=8)
    
    plt.suptitle("Resistance, Disease, and Population by Region", fontsize=14, y=1.02)
    plt.tight_layout()
    save_figure(fig2, "06b_resistance_disease_population")


def plot_regional_summary_dashboard(result: PacificCoastResult):
    """
    Single figure dashboard with population, resistance, disease for each region.
    """
    print("\n📈 Generating regional summary dashboard...")
    
    trajectories = result.get_region_trajectories()
    years = result.years
    config = result.config
    
    n_regions = len(REGION_ORDER)
    fig, axes = plt.subplots(3, n_regions, figsize=(3 * n_regions, 10))
    
    metrics = ['population_ratio', 'resistance', 'disease_prevalence']
    titles = ['Population Ratio', 'Resistance Freq', 'Disease Prevalence']
    colors_map = {'population_ratio': 'blue', 'resistance': 'green', 'disease_prevalence': 'red'}
    
    for col, region_id in enumerate(REGION_ORDER):
        region = config.regions[region_id]
        
        for row, metric in enumerate(metrics):
            ax = axes[row, col]
            data = trajectories[region_id][metric]
            
            ax.plot(years, data, color=colors_map[metric], linewidth=1.5)
            ax.fill_between(years, 0, data, color=colors_map[metric], alpha=0.2)
            
            if row == 0:
                ax.set_title(region.short_name, fontsize=10, fontweight='bold', color=region.color)
            
            if col == 0:
                ax.set_ylabel(titles[row], fontsize=9)
            
            if row == 2:
                ax.set_xlabel("Year", fontsize=9)
            
            ax.set_ylim(0, 1.0)
            ax.set_xlim(0, config.n_years)
            ax.tick_params(labelsize=7)
            ax.grid(True, alpha=0.2)
    
    plt.suptitle("Pacific Coast Regional Dashboard: Population, Resistance, Disease", 
                 fontsize=16, y=1.02)
    plt.tight_layout()
    save_figure(fig, "07_regional_dashboard")


def _run_outplanting_scenario(config, seed, outplanting_resistance):
    """
    Run a simulation with outplanting at year 15 in Salish Sea.
    
    Outplants 500 individuals per site with given resistance level
    (all loci set to that frequency).
    """
    sim = PacificCoastSimulation(config, seed=seed)
    boundaries = get_site_region_boundaries(sim.sites)
    salish_start, salish_end = boundaries["salish_sea"]
    outplant_year = 15
    outplant_n = 500

    result = PacificCoastResult(
        config=config,
        sites=sim.sites,
        states=[],
        larval_connectivity=sim.larval_connectivity,
        disease_connectivity=sim.disease_connectivity,
    )

    for year in range(config.n_years):
        # Outplant at specified year
        if year == outplant_year:
            for i in range(salish_start, salish_end):
                old_pop = sim.populations[i]
                new_pop = old_pop + outplant_n
                if new_pop > 0:
                    weight_old = old_pop / new_pop
                    weight_new = outplant_n / new_pop
                    for locus in range(config.n_loci):
                        sim.resistance_freqs[i, locus] = (
                            weight_old * sim.resistance_freqs[i, locus]
                            + weight_new * outplanting_resistance
                        )
                sim.adults[i] += outplant_n * 0.6
                sim.juveniles[i] += outplant_n * 0.4
                sim.populations[i] = sim.adults[i] + sim.juveniles[i]

        state = sim._simulate_year(year)
        result.states.append(state)

        if sim.populations.sum() < 1:
            from pycnopodia.pacific_coast import PacificCoastState
            for y in range(year + 1, config.n_years):
                result.states.append(PacificCoastState(
                    year=y,
                    populations=np.zeros(sim.n_sites),
                    resistance_freqs=np.zeros((sim.n_sites, config.n_loci)),
                    disease_prevalence=np.zeros(sim.n_sites),
                    temperatures=sim.temperatures.copy(),
                    initial_populations=sim.initial_populations,
                    locus_effects=sim.locus_effects,
                ))
            break

    return result


def run_intervention_comparison(seed: int = 42):
    """
    Run multiple intervention scenarios and compare.
    
    Scenarios:
    1. No intervention (baseline)
    2. Outplanting 50% resistance (all loci at 0.50)
    3. Outplanting 95% resistance (all loci at 0.95)
    """
    print("\n🔬 Running intervention comparison scenarios...")
    
    config = PacificCoastConfig()
    
    # Run baseline
    print("  Running baseline (no intervention)...")
    sim_base = PacificCoastSimulation(config, seed=seed)
    result_baseline = sim_base.run()
    
    # Run 50% resistance outplanting
    print("  Running 50% resistance outplanting...")
    result_50 = _run_outplanting_scenario(config, seed, outplanting_resistance=0.50)
    
    # Run 95% resistance outplanting
    print("  Running 95% resistance outplanting...")
    result_95 = _run_outplanting_scenario(config, seed, outplanting_resistance=0.95)
    
    years = result_baseline.years
    traj_base = result_baseline.get_region_trajectories()
    traj_50 = result_50.get_region_trajectories()
    traj_95 = result_95.get_region_trajectories()
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Helper: mean pop ratio across regions
    def mean_pop(traj):
        total = np.zeros(len(years))
        for r in REGION_ORDER:
            total += traj[r]["population_ratio"]
        return total / len(REGION_ORDER)
    
    pop_base = mean_pop(traj_base)
    pop_50 = mean_pop(traj_50)
    pop_95 = mean_pop(traj_95)
    
    # --- Panel 1: Overall population comparison ---
    ax1 = axes[0, 0]
    ax1.plot(years, pop_base, label='No intervention', linewidth=2, color='black')
    ax1.plot(years, pop_50, label='50% resistance outplant', linewidth=2, color='orange', linestyle='--')
    ax1.plot(years, pop_95, label='95% resistance outplant', linewidth=2, color='green', linestyle='-.')
    ax1.axvline(x=15, color='blue', linestyle=':', alpha=0.5, label='Outplant year')
    ax1.set_xlabel("Year", fontsize=11)
    ax1.set_ylabel("Mean Population Ratio", fontsize=11)
    ax1.set_title("Population Under Outplanting Scenarios", fontsize=12)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 1.0)
    
    # --- Panel 2: Salish Sea focus ---
    ax2 = axes[0, 1]
    ax2.plot(years, traj_base["salish_sea"]["population_ratio"], label='Baseline', linewidth=2, color='black')
    ax2.plot(years, traj_50["salish_sea"]["population_ratio"], label='50% resist', linewidth=2, color='orange', linestyle='--')
    ax2.plot(years, traj_95["salish_sea"]["population_ratio"], label='95% resist', linewidth=2, color='green', linestyle='-.')
    ax2.set_xlabel("Year", fontsize=11)
    ax2.set_ylabel("Population Ratio", fontsize=11)
    ax2.set_title("Salish Sea - Outplanting Target", fontsize=12)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0.1, color='red', linestyle=':', alpha=0.5)
    
    # --- Panel 3: Resistance evolution in Salish Sea ---
    ax3 = axes[1, 0]
    ax3.plot(years, traj_base["salish_sea"]["resistance"], label='Baseline', linewidth=2, color='black')
    ax3.plot(years, traj_50["salish_sea"]["resistance"], label='50% resist', linewidth=2, color='orange', linestyle='--')
    ax3.plot(years, traj_95["salish_sea"]["resistance"], label='95% resist', linewidth=2, color='green', linestyle='-.')
    ax3.set_xlabel("Year", fontsize=11)
    ax3.set_ylabel("Mean Resistance Frequency", fontsize=11)
    ax3.set_title("Resistance Evolution - Salish Sea", fontsize=12)
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(0, 1.0)
    
    # --- Panel 4: Summary bar chart ---
    ax4 = axes[1, 1]
    scenarios = ['No Intervention', '50% Resist', '95% Resist']
    final_pops = [pop_base[-1], pop_50[-1], pop_95[-1]]
    colors = ['gray', 'orange', 'green']
    
    bars = ax4.bar(scenarios, final_pops, color=colors, edgecolor='black')
    ax4.set_ylabel("Final Mean Population Ratio", fontsize=11)
    ax4.set_title("Final Outcomes by Outplanting Resistance", fontsize=12)
    ax4.set_ylim(0, max(final_pops) * 1.3 + 0.05)
    ax4.axhline(y=0.1, color='red', linestyle='--', label='10% threshold')
    ax4.legend()
    
    plt.suptitle("Outplanting Scenario Comparison: 50% vs 95% Resistance\n(Pacific Coast Pycnopodia)", fontsize=14, y=1.02)
    plt.tight_layout()
    save_figure(fig, "08_intervention_comparison")


def print_summary_statistics(result: PacificCoastResult):
    """Print summary statistics by region."""
    print("\n" + "="*60)
    print("PACIFIC COAST SIMULATION SUMMARY")
    print("="*60)
    
    config = result.config
    
    initial_state = result.states[0]
    final_state = result.states[-1]
    
    # Get state at peak disease
    peak_disease_year = config.disease_onset_year + 5
    if peak_disease_year < len(result.states):
        peak_state = result.states[peak_disease_year]
    else:
        peak_state = final_state
    
    print(f"\nSimulation: {config.n_years} years, {len(result.sites)} sites, 8 regions")
    print(f"Disease onset: Year {config.disease_onset_year} in {config.disease_origin_region}")
    
    print("\n" + "-"*60)
    print("REGIONAL OUTCOMES")
    print("-"*60)
    print(f"{'Region':<20} {'Initial':<10} {'Final':<10} {'Decline':<10} {'Resist':<10} {'Disease':<10}")
    print("-"*60)
    
    initial_summary = initial_state.get_region_summary(result.sites, config)
    final_summary = final_state.get_region_summary(result.sites, config)
    
    total_initial = 0
    total_final = 0
    
    for region_id in REGION_ORDER:
        region = config.regions[region_id]
        init = initial_summary[region_id]
        fin = final_summary[region_id]
        
        total_initial += init["population"]
        total_final += fin["population"]
        
        decline = 1 - fin["population_ratio"]
        
        print(f"{region.short_name:<20} {init['population']:>10.0f} {fin['population']:>10.0f} "
              f"{decline:>9.1%} {fin['mean_resistance']:>9.2%} {fin['mean_disease_prevalence']:>9.2%}")
    
    print("-"*60)
    total_decline = 1 - (total_final / total_initial)
    print(f"{'TOTAL':<20} {total_initial:>10.0f} {total_final:>10.0f} {total_decline:>9.1%}")
    
    print("\n" + "-"*60)
    print("CALIBRATION CHECK (vs observed SSWD patterns)")
    print("-"*60)
    
    # Check if results match observed patterns
    checks = [
        ("SE Alaska North (fjords) ~60% survival", final_summary.get("se_alaska_north", {}).get("population_ratio", 0), 0.3, 0.7),
        ("SE Alaska South ~10% survival", final_summary.get("se_alaska_south", {}).get("population_ratio", 0), 0.0, 0.20),
        ("BC Fjords ~50% survival (refugia)", final_summary["bc_fjords"]["population_ratio"], 0.3, 0.7),
        ("Salish Sea <5% survival", final_summary["salish_sea"]["population_ratio"], 0.0, 0.10),
        ("S. California ~0% survival", final_summary["s_california"]["population_ratio"], 0.0, 0.05),
        ("Overall ~90% decline", total_decline, 0.80, 0.98),
    ]
    
    for name, value, low, high in checks:
        status = "✓" if low <= value <= high else "✗"
        print(f"  {status} {name}: {value:.1%} (target: {low:.0%}-{high:.0%})")


def _run_moss_landing_outplanting(config, seed, broodstock_resistance):
    """
    Run simulation with outplanting at year 17 (2027) in Central California.
    
    Adds 1000 one-year-old juveniles per site in c_california (2 years to maturity)
    with resistance allele frequencies set to broodstock_resistance at each locus.
    Uses weighted blending with existing population.
    """
    sim = PacificCoastSimulation(config, seed=seed)
    boundaries = get_site_region_boundaries(sim.sites)
    cc_start, cc_end = boundaries["c_california"]
    outplant_year = 17  # 2027 = 2010 + 17
    outplant_n = 1000   # per site

    result = PacificCoastResult(
        config=config,
        sites=sim.sites,
        states=[],
        larval_connectivity=sim.larval_connectivity,
        disease_connectivity=sim.disease_connectivity,
    )

    for year in range(config.n_years):
        if year == outplant_year:
            for i in range(cc_start, cc_end):
                old_pop = sim.populations[i]
                new_pop = old_pop + outplant_n
                if new_pop > 0:
                    weight_old = old_pop / new_pop
                    weight_new = outplant_n / new_pop
                    for locus in range(config.n_loci):
                        sim.resistance_freqs[i, locus] = (
                            weight_old * sim.resistance_freqs[i, locus]
                            + weight_new * broodstock_resistance
                        )
                sim.juveniles[i] += 1000  # 1000 one-year-olds (2 yrs to maturity)
                sim.populations[i] = sim.adults[i] + sim.juveniles[i]

        state = sim._simulate_year(year)
        result.states.append(state)

        if sim.populations.sum() < 1:
            from pycnopodia.pacific_coast import PacificCoastState
            for y in range(year + 1, config.n_years):
                result.states.append(PacificCoastState(
                    year=y,
                    populations=np.zeros(sim.n_sites),
                    resistance_freqs=np.zeros((sim.n_sites, config.n_loci)),
                    disease_prevalence=np.zeros(sim.n_sites),
                    temperatures=sim.temperatures.copy(),
                    initial_populations=sim.initial_populations,
                    locus_effects=sim.locus_effects,
                ))
            break

    return result


def plot_moss_landing_comparison(result_2019, result_2026, baseline):
    """
    Plot comparison of Moss Landing outplanting scenarios.
    
    6-panel figure comparing baseline, 2019 broodstock, and 2026 broodstock.
    """
    print("\n🌟 Generating Moss Landing outplanting comparison plots...")
    
    years = baseline.years
    config = baseline.config
    
    traj_base = baseline.get_region_trajectories()
    traj_2019 = result_2019.get_region_trajectories()
    traj_2026 = result_2026.get_region_trajectories()
    
    # Helper: mean pop ratio across all regions
    def mean_pop(traj):
        total = np.zeros(len(years))
        for r in REGION_ORDER:
            total += traj[r]["population_ratio"]
        return total / len(REGION_ORDER)
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    
    line_styles = [
        ('Baseline (no outplanting)', 'black', '-'),
        ('2019 broodstock (resist=0.15)', '#E67E22', '--'),
        ('2026 broodstock (resist=0.35)', '#27AE60', '-.'),
    ]
    
    # Convert model years to calendar years
    cal_years = np.array(years) + 2010
    
    # --- Panel 1: C. California population over time ---
    ax = axes[0, 0]
    for label, color, ls, traj in zip(
        [s[0] for s in line_styles], [s[1] for s in line_styles], [s[2] for s in line_styles],
        [traj_base, traj_2019, traj_2026]
    ):
        ax.plot(cal_years, traj["c_california"]["population_ratio"], label=label, color=color, linestyle=ls, linewidth=2)
    ax.axvline(x=2027, color='blue', linestyle=':', alpha=0.5, label='Release year (2027)')
    ax.axvline(x=2013, color='red', linestyle=':', alpha=0.5, label='SSWD onset (2013)')
    ax.set_xlabel("Year"); ax.set_ylabel("Population Ratio")
    ax.set_title("Central California Population", fontsize=12, fontweight='bold')
    ax.legend(fontsize=7, loc='upper right'); ax.grid(True, alpha=0.3); ax.set_ylim(0, None)
    
    # --- Panel 2: Overall Pacific Coast population ---
    ax = axes[0, 1]
    for label, color, ls, traj in zip(
        [s[0] for s in line_styles], [s[1] for s in line_styles], [s[2] for s in line_styles],
        [traj_base, traj_2019, traj_2026]
    ):
        ax.plot(cal_years, mean_pop(traj), label=label, color=color, linestyle=ls, linewidth=2)
    ax.axvline(x=2027, color='blue', linestyle=':', alpha=0.5)
    ax.set_xlabel("Year"); ax.set_ylabel("Mean Population Ratio")
    ax.set_title("Overall Pacific Coast Population", fontsize=12, fontweight='bold')
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3); ax.set_ylim(0, None)
    
    # --- Panel 3: C. California resistance evolution ---
    ax = axes[0, 2]
    for label, color, ls, traj in zip(
        [s[0] for s in line_styles], [s[1] for s in line_styles], [s[2] for s in line_styles],
        [traj_base, traj_2019, traj_2026]
    ):
        ax.plot(cal_years, traj["c_california"]["resistance"], label=label, color=color, linestyle=ls, linewidth=2)
    ax.axvline(x=2027, color='blue', linestyle=':', alpha=0.5)
    ax.set_xlabel("Year"); ax.set_ylabel("Mean Resistance Frequency")
    ax.set_title("C. California Resistance Evolution", fontsize=12, fontweight='bold')
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3); ax.set_ylim(0, 1.0)
    
    # --- Panel 4: Disease prevalence in C. California ---
    ax = axes[1, 0]
    for label, color, ls, traj in zip(
        [s[0] for s in line_styles], [s[1] for s in line_styles], [s[2] for s in line_styles],
        [traj_base, traj_2019, traj_2026]
    ):
        ax.plot(cal_years, traj["c_california"]["disease_prevalence"], label=label, color=color, linestyle=ls, linewidth=2)
    ax.axvline(x=2027, color='blue', linestyle=':', alpha=0.5)
    ax.set_xlabel("Year"); ax.set_ylabel("Disease Prevalence")
    ax.set_title("C. California Disease Prevalence", fontsize=12, fontweight='bold')
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3); ax.set_ylim(0, 1.0)
    
    # --- Panel 5: Regional comparison bar chart at year 2060 (year 50) ---
    ax = axes[1, 1]
    year_idx_2060 = min(50, len(years) - 1)
    
    x_pos = np.arange(len(REGION_ORDER))
    width = 0.25
    
    for offset, (label, color, _), traj in zip(
        [-width, 0, width],
        line_styles,
        [traj_base, traj_2019, traj_2026]
    ):
        vals = [traj[r]["population_ratio"][year_idx_2060] for r in REGION_ORDER]
        ax.bar(x_pos + offset, vals, width, label=label, color=color, edgecolor='black', alpha=0.8)
    
    region_labels = [config.regions[r].short_name for r in REGION_ORDER]
    ax.set_xticks(x_pos); ax.set_xticklabels(region_labels, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel("Population Ratio"); ax.set_title("Regional Populations at 2060", fontsize=12, fontweight='bold')
    ax.legend(fontsize=6); ax.grid(True, alpha=0.3, axis='y')
    
    # --- Panel 6: C. California zoomed 2025-2060 ---
    ax = axes[1, 2]
    zoom_start = max(0, 15)   # year 15 = 2025
    zoom_end = min(50, len(years) - 1)  # year 50 = 2060
    zoom_cal = cal_years[zoom_start:zoom_end]
    
    for label, color, ls, traj in zip(
        [s[0] for s in line_styles], [s[1] for s in line_styles], [s[2] for s in line_styles],
        [traj_base, traj_2019, traj_2026]
    ):
        ax.plot(zoom_cal, traj["c_california"]["population_ratio"][zoom_start:zoom_end],
                label=label, color=color, linestyle=ls, linewidth=2)
    ax.axvline(x=2027, color='blue', linestyle=':', alpha=0.5, label='Release year')
    ax.set_xlabel("Year"); ax.set_ylabel("Population Ratio")
    ax.set_title("C. California Zoomed (2025-2060)", fontsize=12, fontweight='bold')
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3); ax.set_ylim(0, None)
    
    plt.suptitle(
        "Moss Landing Outplanting: 2019 vs 2026 Broodstock\n"
        "(1000 individuals/site in Central CA, released 2027)",
        fontsize=15, fontweight='bold', y=1.02
    )
    plt.tight_layout()
    save_figure(fig, "09_moss_landing_outplanting")


def run_moss_landing_scenarios(seed: int = 42):
    """Run Moss Landing outplanting comparison scenarios."""
    print("\n🌊 Running Moss Landing outplanting scenarios...")
    
    config = PacificCoastConfig(n_years=100)
    
    print("  Running baseline...")
    sim_base = PacificCoastSimulation(config, seed=seed)
    baseline = sim_base.run()
    
    print("  Running 2019 broodstock (resistance=0.15)...")
    result_2019 = _run_moss_landing_outplanting(config, seed, broodstock_resistance=0.15)
    
    print("  Running 2026 broodstock (resistance=0.35)...")
    result_2026 = _run_moss_landing_outplanting(config, seed, broodstock_resistance=0.35)
    
    plot_moss_landing_comparison(result_2019, result_2026, baseline)
    
    # Print summary
    traj_base = baseline.get_region_trajectories()
    traj_2019 = result_2019.get_region_trajectories()
    traj_2026 = result_2026.get_region_trajectories()
    
    print("\n" + "="*60)
    print("MOSS LANDING OUTPLANTING SUMMARY")
    print("="*60)
    
    for name, traj in [("Baseline", traj_base), ("2019 Broodstock", traj_2019), ("2026 Broodstock", traj_2026)]:
        cc_pop_2060 = traj["c_california"]["population_ratio"][min(50, len(baseline.years)-1)]
        cc_resist_final = traj["c_california"]["resistance"][-1]
        cc_disease_final = traj["c_california"]["disease_prevalence"][-1]
        print(f"\n  {name}:")
        print(f"    C. California pop ratio @ 2060: {cc_pop_2060:.3f}")
        print(f"    C. California final resistance:  {cc_resist_final:.3f}")
        print(f"    C. California final disease:     {cc_disease_final:.3f}")


def export_json(result: PacificCoastResult, output_path: str = "dashboard/data.json"):
    """
    Export simulation results to JSON format for dashboard visualization.
    
    Args:
        result: PacificCoastResult from simulation
        output_path: Where to save the JSON file
    """
    import json
    from pathlib import Path
    
    print(f"\n💾 Exporting data to {output_path}...")
    
    config = result.config
    
    # Build sites data with coordinates (using latitude as y-coordinate)
    # For x-coordinate, we'll use longitude proxy: alternate slightly east/west
    sites_data = []
    for i, site in enumerate(result.sites):
        # Simple longitude assignment: oscillate around -125
        lon = -125.0 + (i % 10 - 5) * 0.5
        sites_data.append({
            "idx": site.idx,
            "lat": float(site.latitude),
            "lon": float(lon),
            "region": site.region_id,
            "is_refugia": bool(site.is_refugia),
            "temperature": float(site.temperature),
        })
    
    # Build sparse connectivity matrix (only edges > 0.001)
    connectivity_sparse = []
    threshold = 0.001
    for i in range(result.larval_connectivity.shape[0]):
        row = []
        for j in range(result.larval_connectivity.shape[1]):
            val = float(result.larval_connectivity[i, j])
            if val > threshold:
                row.append(val)
            else:
                row.append(0)
        connectivity_sparse.append(row)
    
    # Extract timeseries data
    n_sites = len(result.sites)
    n_years = len(result.states)
    
    populations = []
    disease = []
    resistance = []
    
    for site_idx in range(n_sites):
        site_pops = []
        site_disease = []
        site_resist = []
        
        for state in result.states:
            site_pops.append(float(state.populations[site_idx]))
            site_disease.append(float(state.disease_prevalence[site_idx]))
            # Mean resistance across loci
            site_resist.append(float(state.resistance_freqs[site_idx].mean()))
        
        populations.append(site_pops)
        disease.append(site_disease)
        resistance.append(site_resist)
    
    # Build region info with colors and names
    regions_info = {}
    for region_id, region in config.regions.items():
        regions_info[region_id] = {
            "name": region.name,
            "short_name": region.short_name,
            "color": region.color,
            "base_temperature": float(region.base_temperature),
        }
    
    # Assemble final data structure
    data = {
        "sites": sites_data,
        "connectivity": connectivity_sparse,
        "years": [int(s.year) for s in result.states],
        "timeseries": {
            "populations": populations,
            "disease": disease,
            "resistance": resistance,
        },
        "regions": regions_info,
        "config": {
            "disease_onset_year": config.disease_onset_year,
            "disease_origin_region": config.disease_origin_region,
        }
    }
    
    # Write to file
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2)
    
    print(f"  ✓ Exported {n_sites} sites × {n_years} years")
    print(f"  ✓ Connectivity matrix: {len(connectivity_sparse)}×{len(connectivity_sparse[0])}")
    print(f"  ✓ File size: {output_path.stat().st_size / 1024:.1f} KB")


class EnsembleResult:
    """Container for ensemble simulation results."""
    
    def __init__(self, mean_trajectories, ci_low, ci_high, all_runs, config):
        """
        Args:
            mean_trajectories: dict[region_id -> dict[metric -> array]]
            ci_low: dict[region_id -> dict[metric -> array]] (2.5th percentile)
            ci_high: dict[region_id -> dict[metric -> array]] (97.5th percentile)
            all_runs: list of PacificCoastResult objects
            config: PacificCoastConfig
        """
        self.mean_trajectories = mean_trajectories
        self.ci_low = ci_low
        self.ci_high = ci_high
        self.all_runs = all_runs
        self.config = config
        self.n_runs = len(all_runs)
        self.years = all_runs[0].years if all_runs else []


def run_ensemble(n_runs: int = 20, config: PacificCoastConfig = None, verbose: bool = True) -> EnsembleResult:
    """
    Run ensemble of simulations with different random seeds.
    
    Args:
        n_runs: Number of replicate simulations to run
        config: PacificCoastConfig (defaults to standard config)
        verbose: Print progress updates
    
    Returns:
        EnsembleResult with mean trajectories and 95% confidence intervals
    """
    if config is None:
        config = PacificCoastConfig(n_years=100)
    
    if verbose:
        print(f"\n🔄 Running ensemble with {n_runs} replicates...")
    
    all_results = []
    
    for seed in range(n_runs):
        if verbose and (seed + 1) % 5 == 0:
            print(f"  Completed {seed + 1}/{n_runs} runs...")
        
        sim = PacificCoastSimulation(config, seed=seed)
        result = sim.run()
        all_results.append(result)
    
    if verbose:
        print(f"  ✓ All {n_runs} runs complete!")
    
    # Collect trajectories from all runs
    # Structure: runs[run_idx][region_id][metric][year]
    n_years = len(all_results[0].years)
    
    # Initialize storage for all runs
    all_trajectories = {region_id: {
        "population_ratio": np.zeros((n_runs, n_years)),
        "resistance": np.zeros((n_runs, n_years)),
        "disease_prevalence": np.zeros((n_runs, n_years)),
    } for region_id in REGION_ORDER}
    
    # Collect data from each run
    for run_idx, result in enumerate(all_results):
        trajectories = result.get_region_trajectories()
        for region_id in REGION_ORDER:
            for metric in ["population_ratio", "resistance", "disease_prevalence"]:
                all_trajectories[region_id][metric][run_idx, :] = trajectories[region_id][metric]
    
    # Compute mean and confidence intervals
    mean_trajectories = {}
    ci_low = {}
    ci_high = {}
    
    for region_id in REGION_ORDER:
        mean_trajectories[region_id] = {}
        ci_low[region_id] = {}
        ci_high[region_id] = {}
        
        for metric in ["population_ratio", "resistance", "disease_prevalence"]:
            data = all_trajectories[region_id][metric]
            mean_trajectories[region_id][metric] = np.mean(data, axis=0)
            ci_low[region_id][metric] = np.percentile(data, 2.5, axis=0)
            ci_high[region_id][metric] = np.percentile(data, 97.5, axis=0)
    
    return EnsembleResult(mean_trajectories, ci_low, ci_high, all_results, config)


def plot_ensemble_trajectories(ensemble: EnsembleResult):
    """
    Plot ensemble mean trajectories with 95% confidence intervals.
    
    Shows mean population ratio per region as solid line with shaded CI band.
    Matches style of plot_population_trajectories_by_region.
    """
    print("\n📊 Generating ensemble trajectory plot...")
    
    years = ensemble.years
    config = ensemble.config
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    for region_id in REGION_ORDER:
        region = config.regions[region_id]
        
        mean = ensemble.mean_trajectories[region_id]["population_ratio"]
        ci_low = ensemble.ci_low[region_id]["population_ratio"]
        ci_high = ensemble.ci_high[region_id]["population_ratio"]
        
        # Plot mean line
        ax.plot(
            years, mean,
            label=region.short_name,
            color=region.color,
            linewidth=2,
        )
        
        # Plot 95% CI as shaded band
        ax.fill_between(
            years, ci_low, ci_high,
            color=region.color,
            alpha=0.2,
        )
    
    ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Baseline')
    ax.axhline(y=0.1, color='gray', linestyle=':', alpha=0.5, label='10% threshold')
    ax.axvline(x=config.disease_onset_year, color='red', linestyle=':', alpha=0.5, label='Disease onset')
    
    ax.set_xlabel("Year", fontsize=12)
    ax.set_ylabel("Population Ratio (relative to initial)", fontsize=12)
    ax.set_title(f"Ensemble Population Trajectories (n={ensemble.n_runs}, mean ± 95% CI)", fontsize=14)
    ax.legend(loc='upper right', fontsize=9)
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 2.0)
    ax.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    save_figure(fig, "10_ensemble_trajectories")


def plot_ensemble_summary(ensemble: EnsembleResult):
    """
    Plot ensemble summary: mean final survival per region with 95% CI error bars.
    
    Bar chart showing final population ratio by region.
    """
    print("\n📊 Generating ensemble summary plot...")
    
    config = ensemble.config
    
    # Get final year index
    final_year_idx = -1
    
    # Extract final values for each region
    regions = []
    means = []
    ci_lows = []
    ci_highs = []
    colors = []
    
    for region_id in REGION_ORDER:
        region = config.regions[region_id]
        regions.append(region.short_name)
        
        mean = ensemble.mean_trajectories[region_id]["population_ratio"][final_year_idx]
        ci_low = ensemble.ci_low[region_id]["population_ratio"][final_year_idx]
        ci_high = ensemble.ci_high[region_id]["population_ratio"][final_year_idx]
        
        means.append(mean)
        ci_lows.append(mean - ci_low)  # Error bar is distance from mean
        ci_highs.append(ci_high - mean)
        colors.append(region.color)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    x_pos = np.arange(len(regions))
    
    # Plot bars with error bars
    ax.bar(
        x_pos, means,
        color=colors,
        edgecolor='black',
        alpha=0.8,
    )
    
    ax.errorbar(
        x_pos, means,
        yerr=[ci_lows, ci_highs],
        fmt='none',
        ecolor='black',
        capsize=5,
        capthick=2,
        alpha=0.7,
    )
    
    ax.set_xticks(x_pos)
    ax.set_xticklabels(regions, rotation=45, ha='right')
    ax.set_ylabel("Final Population Ratio (mean ± 95% CI)", fontsize=12)
    ax.set_title(f"Ensemble Final Survival by Region (n={ensemble.n_runs} replicates)", fontsize=14)
    ax.axhline(y=0.1, color='red', linestyle='--', alpha=0.5, label='10% threshold')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    save_figure(fig, "10b_ensemble_summary")


def print_ensemble_statistics(ensemble: EnsembleResult):
    """Print ensemble statistics: mean ± std for each region."""
    print("\n" + "="*70)
    print(f"ENSEMBLE STATISTICS (n={ensemble.n_runs} replicates)")
    print("="*70)
    
    config = ensemble.config
    final_year_idx = -1
    
    print(f"\n{'Region':<20} {'Mean':<12} {'Std':<12} {'95% CI Range':<20} {'Obs Target':<15}")
    print("-"*70)
    
    for region_id in REGION_ORDER:
        region = config.regions[region_id]
        
        # Get final population ratio across all runs
        final_pops = np.array([
            run.get_region_trajectories()[region_id]["population_ratio"][final_year_idx]
            for run in ensemble.all_runs
        ])
        
        mean = np.mean(final_pops)
        std = np.std(final_pops)
        ci_low = np.percentile(final_pops, 2.5)
        ci_high = np.percentile(final_pops, 97.5)
        
        ci_range = f"[{ci_low:.3f}, {ci_high:.3f}]"
        obs_survival = f"{region.post_sswd_survival:.1%}"
        
        print(f"{region.short_name:<20} {mean:>6.3f} ± {std:<5.3f} {ci_range:<20} {obs_survival:<15}")
    
    print("-"*70)
    
    # Overall statistics
    all_final_pops = []
    for run in ensemble.all_runs:
        traj = run.get_region_trajectories()
        total = sum(traj[r]["population_ratio"][final_year_idx] for r in REGION_ORDER) / len(REGION_ORDER)
        all_final_pops.append(total)
    
    overall_mean = np.mean(all_final_pops)
    overall_std = np.std(all_final_pops)
    overall_ci_low = np.percentile(all_final_pops, 2.5)
    overall_ci_high = np.percentile(all_final_pops, 97.5)
    
    print(f"\n{'OVERALL MEAN':<20} {overall_mean:>6.3f} ± {overall_std:<5.3f} [{overall_ci_low:.3f}, {overall_ci_high:.3f}]")
    
    print("\n" + "="*70)
    print("VARIANCE CHECK")
    print("="*70)
    
    # Show which regions have highest variance
    variances = []
    for region_id in REGION_ORDER:
        final_pops = np.array([
            run.get_region_trajectories()[region_id]["population_ratio"][final_year_idx]
            for run in ensemble.all_runs
        ])
        cv = np.std(final_pops) / (np.mean(final_pops) + 1e-10)  # Coefficient of variation
        variances.append((region_id, cv, np.std(final_pops)))
    
    variances.sort(key=lambda x: x[2], reverse=True)
    
    print(f"\n{'Region':<20} {'Std Dev':<12} {'Coeff. of Var.':<15}")
    print("-"*50)
    for region_id, cv, std in variances:
        region = config.regions[region_id]
        print(f"{region.short_name:<20} {std:>8.4f}      {cv:>8.2f}")


def main():
    """Run Pacific Coast simulation with all visualizations."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Run Pacific Coast Pycnopodia simulation")
    parser.add_argument('--export-json', action='store_true', 
                       help='Export simulation data to dashboard/data.json')
    parser.add_argument('--json-only', action='store_true',
                       help='Only export JSON, skip visualization plots')
    parser.add_argument('--moss-landing', action='store_true',
                       help='Run Moss Landing outplanting scenarios (2019 vs 2026 broodstock)')
    parser.add_argument('--ensemble', type=int, nargs='?', const=20, default=0,
                       help='Run ensemble averaging with N replicates (default: 20)')
    args = parser.parse_args()
    
    print("="*60)
    print("PACIFIC COAST PYCNOPODIA SIMULATION")
    print("="*60)
    
    # Check if ensemble mode is requested
    if args.ensemble > 0:
        # Run ensemble
        config = PacificCoastConfig(n_years=100)
        ensemble = run_ensemble(n_runs=args.ensemble, config=config, verbose=True)
        
        # Generate ensemble plots
        plot_ensemble_trajectories(ensemble)
        plot_ensemble_summary(ensemble)
        
        # Print ensemble statistics
        print_ensemble_statistics(ensemble)
        
        print("\n" + "="*60)
        print(f"Ensemble figures saved to: {FIGURE_DIR}")
        print("="*60)
        
        return ensemble
    
    # Otherwise run single simulation
    print("\n🌊 Running simulation...")
    config = PacificCoastConfig(n_years=100)
    sim = PacificCoastSimulation(config, seed=43)
    result = sim.run()
    
    print(f"  Completed: {len(result.sites)} sites, {len(result.states)} years")
    
    # Export JSON if requested
    if args.export_json or args.json_only:
        export_json(result)
    
    # Generate visualizations unless json-only
    if not args.json_only:
        print("\n📊 Generating visualizations...")
        
        plot_geographic_structure(result)
        plot_connectivity_matrices(result)
        plot_disease_spread_spacetime(result)
        plot_population_trajectories_by_region(result)
        plot_temperature_effects(result)
        plot_resistance_evolution(result)
        plot_regional_summary_dashboard(result)
        run_intervention_comparison(seed=42)
        
        # Moss Landing scenarios
        if args.moss_landing:
            run_moss_landing_scenarios(seed=42)
        
        # Print summary
        print_summary_statistics(result)
        
        print("\n" + "="*60)
        print(f"All figures saved to: {FIGURE_DIR}")
        print("="*60)
    
    return result


if __name__ == "__main__":
    main()
