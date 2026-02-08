#!/usr/bin/env python3
"""
Outplanting intervention experiments.

Adds 95% resistant 1-year-old juvenile stars near Monterey Bay at year 2027
(model year 24) and compares outcomes across release sizes.

Release sizes: 100, 10000, 100000
Each run as 10-seed ensemble for statistical power.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pycnopodia.geo_model import GeoSimulation, GeoConfig, GeoState
from data.real_sites import ALL_SITES
import json


MONTEREY_IDX = next(i for i, s in enumerate(ALL_SITES) if 'Monterey' in s.name)
OUTPLANT_YEAR = 24  # 2003 + 24 = 2027
OUTPLANT_SEASON = 1  # Spring (after winter spawn, before summer disease)
OUTPLANT_STEP = OUTPLANT_YEAR * 4 + OUTPLANT_SEASON

RELEASE_SIZES = [0, 100, 10_000, 100_000]
N_SEEDS = 10
N_YEARS = 80

REGION_ORDER = [
    "se_alaska_north", "se_alaska_south", "bc_fjords", "bc_outer",
    "salish_sea", "wa_or_outer", "n_california", "c_california", "s_california"
]

REGION_LABELS = {
    "se_alaska_north": "SE AK N", "se_alaska_south": "SE AK S",
    "bc_fjords": "BC Fjords", "bc_outer": "BC Coast",
    "salish_sea": "Salish Sea", "wa_or_outer": "WA/OR Coast",
    "n_california": "N. CA", "c_california": "C. CA", "s_california": "S. CA"
}


class OutplantingSimulation(GeoSimulation):
    """GeoSimulation with outplanting intervention."""
    
    def __init__(self, config, sites=None, seed=0, 
                 release_size=0, release_site=MONTEREY_IDX,
                 release_step=OUTPLANT_STEP, resistance_freq=0.95):
        super().__init__(config, sites, seed)
        self.release_size = release_size
        self.release_site = release_site
        self.release_step = release_step
        self.release_resistance_freq = resistance_freq
        self._released = False
    
    def _simulate_season(self, step, year, season):
        """Override to inject outplanted individuals."""
        # Run normal simulation step
        state = super()._simulate_season(step, year, season)
        
        # Outplant at the specified step
        if step == self.release_step and not self._released and self.release_size > 0:
            self._released = True
            site = self.release_site
            
            # Add juveniles (1-year-old = 4 seasonal steps old, not yet mature)
            self.juveniles[site] += self.release_size
            self.populations[site] += self.release_size
            
            # Set resistance allele frequencies high for this site
            # Weighted average: existing pop freqs + new resistant individuals
            existing_pop = state.populations[site]  # Before outplanting
            new_pop = existing_pop + self.release_size
            
            if new_pop > 0:
                weight_existing = existing_pop / new_pop
                weight_new = self.release_size / new_pop
                
                for l in range(self.config.n_loci):
                    self.resistance_freqs[site, l] = (
                        weight_existing * self.resistance_freqs[site, l] +
                        weight_new * self.release_resistance_freq
                    )
                    self.resistance_freqs[site, l] = np.clip(
                        self.resistance_freqs[site, l], 0.001, 
                        self.config.max_resistance_freq
                    )
            
            # Update state to reflect outplanting
            state = GeoState(
                step=step, year=year, season=season,
                populations=self.populations.copy(),
                adults=self.adults.copy(),
                juveniles=self.juveniles.copy(),
                disease_prevalence=self.disease_prevalence.copy(),
                resistance_freqs=self.resistance_freqs.copy(),
                temperatures=self.temperatures.copy(),
                locus_effects=self.locus_effects.copy(),
            )
        
        return state


def run_experiment():
    """Run all outplanting scenarios."""
    config = GeoConfig(n_years=N_YEARS)
    sites = list(ALL_SITES)
    
    results = {}  # release_size -> list of GeoResult
    
    for size in RELEASE_SIZES:
        label = f"{size:,}" if size > 0 else "Control"
        print(f"\n{'='*60}")
        print(f"Release size: {label} (95% resistant, 1yr juveniles at Monterey Bay)")
        print(f"Release time: Year {OUTPLANT_YEAR} ({2003+OUTPLANT_YEAR}), Spring")
        print(f"{'='*60}")
        
        runs = []
        for seed in range(N_SEEDS):
            sim = OutplantingSimulation(
                config=config, sites=sites, seed=seed,
                release_size=size
            )
            result = sim.run()
            runs.append(result)
            print(f"  Seed {seed+1}/{N_SEEDS} complete")
        
        results[size] = runs
    
    return results, sites, config


def plot_monterey_comparison(results, sites, config, outdir):
    """Plot Monterey Bay + C. California trajectories across scenarios."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    colors = {0: '#888888', 100: '#3498db', 10_000: '#e67e22', 100_000: '#e74c3c'}
    labels = {0: 'Control', 100: '100 released', 10_000: '10,000 released', 100_000: '100,000 released'}
    
    n_steps = config.n_years * config.seasons_per_year
    years = [2003 + s/4 for s in range(n_steps)]
    
    # Panel 1: Monterey Bay population
    ax = axes[0, 0]
    for size in RELEASE_SIZES:
        trajectories = []
        for result in results[size]:
            traj = [result.states[s].populations[MONTEREY_IDX] for s in range(len(result.states))]
            trajectories.append(traj)
        mean_traj = np.mean(trajectories, axis=0)
        std_traj = np.std(trajectories, axis=0)
        ax.semilogy(years[:len(mean_traj)], mean_traj + 1, color=colors[size], label=labels[size], linewidth=2)
        ax.fill_between(years[:len(mean_traj)], 
                        np.maximum(mean_traj - std_traj, 1),
                        mean_traj + std_traj + 1,
                        color=colors[size], alpha=0.15)
    ax.axvline(x=2027, color='green', linestyle='--', alpha=0.5, label='Release (2027)')
    ax.axvline(x=2013, color='red', linestyle='--', alpha=0.3, label='SSWD onset')
    ax.set_title('Monterey Bay Population', fontsize=13)
    ax.set_ylabel('Population (log scale)')
    ax.legend(fontsize=9)
    ax.set_xlim(2003, 2083)
    
    # Panel 2: C. California total
    ax = axes[0, 1]
    ca_idx = [i for i, s in enumerate(sites) if s.region == 'c_california']
    for size in RELEASE_SIZES:
        trajectories = []
        for result in results[size]:
            traj = [sum(result.states[s].populations[j] for j in ca_idx) for s in range(len(result.states))]
            trajectories.append(traj)
        mean_traj = np.mean(trajectories, axis=0)
        ax.semilogy(years[:len(mean_traj)], mean_traj + 1, color=colors[size], label=labels[size], linewidth=2)
    ax.axvline(x=2027, color='green', linestyle='--', alpha=0.5)
    ax.axvline(x=2013, color='red', linestyle='--', alpha=0.3)
    ax.set_title('Central California Total Population', fontsize=13)
    ax.set_ylabel('Population (log scale)')
    ax.legend(fontsize=9)
    ax.set_xlim(2003, 2083)
    
    # Panel 3: Resistance at Monterey
    ax = axes[1, 0]
    for size in RELEASE_SIZES:
        trajectories = []
        for result in results[size]:
            traj = [result.states[s].resistance_freqs[MONTEREY_IDX].mean() for s in range(len(result.states))]
            trajectories.append(traj)
        mean_traj = np.mean(trajectories, axis=0)
        ax.plot(years[:len(mean_traj)], mean_traj, color=colors[size], label=labels[size], linewidth=2)
    ax.axvline(x=2027, color='green', linestyle='--', alpha=0.5)
    ax.set_title('Mean Resistance Allele Frequency (Monterey)', fontsize=13)
    ax.set_ylabel('Frequency')
    ax.legend(fontsize=9)
    ax.set_xlim(2003, 2083)
    
    # Panel 4: Range-wide total
    ax = axes[1, 1]
    for size in RELEASE_SIZES:
        trajectories = []
        for result in results[size]:
            traj = [result.states[s].populations.sum() for s in range(len(result.states))]
            trajectories.append(traj)
        mean_traj = np.mean(trajectories, axis=0)
        ax.semilogy(years[:len(mean_traj)], mean_traj + 1, color=colors[size], label=labels[size], linewidth=2)
    ax.axvline(x=2027, color='green', linestyle='--', alpha=0.5)
    ax.axvline(x=2013, color='red', linestyle='--', alpha=0.3)
    ax.set_title('Range-wide Total Population', fontsize=13)
    ax.set_ylabel('Population (log scale)')
    ax.legend(fontsize=9)
    ax.set_xlim(2003, 2083)
    
    fig.suptitle('Outplanting Intervention: 95% Resistant Juveniles at Monterey Bay (2027)', fontsize=15, y=1.01)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'outplanting_comparison.png'), dpi=150, bbox_inches='tight')
    print(f"Saved: {outdir}/outplanting_comparison.png")
    plt.close()


def print_summary(results, sites, config):
    """Print numerical summary of all scenarios."""
    print(f"\n{'='*80}")
    print("OUTPLANTING EXPERIMENT SUMMARY")
    print(f"{'='*80}")
    
    ca_idx = [i for i, s in enumerate(sites) if s.region == 'c_california']
    
    for size in RELEASE_SIZES:
        label = f"{size:>7,}" if size > 0 else "Control"
        
        monterey_y50 = []
        monterey_y80 = []
        ca_y50 = []
        ca_y80 = []
        total_y80 = []
        resist_y80 = []
        
        for result in results[size]:
            n_states = len(result.states)
            s50 = min(50 * 4, n_states - 1)
            s80 = n_states - 1
            
            monterey_y50.append(result.states[s50].populations[MONTEREY_IDX])
            monterey_y80.append(result.states[s80].populations[MONTEREY_IDX])
            ca_y50.append(sum(result.states[s50].populations[j] for j in ca_idx))
            ca_y80.append(sum(result.states[s80].populations[j] for j in ca_idx))
            total_y80.append(result.states[s80].populations.sum())
            resist_y80.append(result.states[s80].resistance_freqs[MONTEREY_IDX].mean())
        
        print(f"\n  Release: {label}")
        print(f"    Monterey Y50: {np.mean(monterey_y50):>12,.0f} ± {np.std(monterey_y50):>10,.0f}")
        print(f"    Monterey Y80: {np.mean(monterey_y80):>12,.0f} ± {np.std(monterey_y80):>10,.0f}")
        print(f"    C. CA    Y50: {np.mean(ca_y50):>12,.0f} ± {np.std(ca_y50):>10,.0f}")
        print(f"    C. CA    Y80: {np.mean(ca_y80):>12,.0f} ± {np.std(ca_y80):>10,.0f}")
        print(f"    Total    Y80: {np.mean(total_y80):>12,.0f} ± {np.std(total_y80):>10,.0f}")
        print(f"    Resist   Y80: {np.mean(resist_y80):>8.4f} ± {np.std(resist_y80):.4f}")
    
    print(f"\n{'='*80}")


if __name__ == "__main__":
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures', 'outplanting')
    os.makedirs(outdir, exist_ok=True)
    
    results, sites, config = run_experiment()
    plot_monterey_comparison(results, sites, config, outdir)
    print_summary(results, sites, config)
