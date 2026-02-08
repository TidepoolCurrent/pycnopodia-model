#!/usr/bin/env python3
"""Sweep fecundity values to find stable pre-disease equilibrium.

Fecundity = settled recruits per breeding pair.
Breeders = N_adults × SRS (0.001).
Recruits = n_breeders/2 (pairs) × fecundity × Allee_fertilization × stochastic_breeding

Tests: 100, 1000, 10000, 100000
Also tests wild WA broodstock variant (outplanting from Washington).
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from copy import deepcopy
from pycnopodia.geo_model import GeoSimulation, GeoConfig, GeoState
from data.real_sites import ALL_SITES

REGION_ORDER = [
    "se_alaska_north", "se_alaska_south", "bc_fjords", "bc_outer",
    "salish_sea", "wa_or_outer", "n_california", "c_california", "s_california"
]
REGION_LABELS = {
    "se_alaska_north": "SE AK N", "se_alaska_south": "SE AK S",
    "bc_fjords": "BC Fjords", "bc_outer": "BC Coast",
    "salish_sea": "Salish Sea", "wa_or_outer": "WA/OR Coast",
    "n_california": "N. CA", "c_california": "C. CA", "s_california": "S. CA",
}
REGION_COLORS = {
    "s_california": "#e74c3c", "c_california": "#e67e22", "n_california": "#e91e95",
    "wa_or_outer": "#795548", "salish_sea": "#9b59b6", "bc_outer": "#27ae60",
    "bc_fjords": "#00bcd4", "se_alaska_south": "#2196f3", "se_alaska_north": "#00e5ff",
}

# Monkey-patch _reproduce to use fecundity parameter
_original_reproduce = GeoSimulation._reproduce

def _patched_reproduce(self):
    """Reproduction with fecundity = settled recruits per breeding pair."""
    fecundity = getattr(self.config, '_fecundity', None)
    if fecundity is None:
        return _original_reproduce(self)
    
    recruits = np.zeros(self.n_sites)
    srs = self.config.srs_breeding_fraction
    
    for i in range(self.n_sites):
        n_adults = self.adults[i]
        if n_adults < self.config.allee_threshold:
            fertilization = 0.01
        else:
            h = self.config.allee_half_sat
            fertilization = n_adults**2 / (n_adults**2 + h**2)
        
        # Number of real breeders (SRS)
        n_real_adults = n_adults  # already in units of thousands
        lam = n_real_adults * srs
        if lam > 1e15:
            n_breeders = int(lam)  # Deterministic for huge populations
        else:
            n_breeders = max(0, self.rng.poisson(lam))
        n_pairs = n_breeders // 2
        
        if n_pairs == 0:
            recruits[i] = 0
            continue
        
        # Each pair produces fecundity settled recruits (with stochastic variance)
        breeding_success = self.rng.beta(2, 20) / 0.091  # normalize by beta mean
        recruits[i] = min(n_pairs * fecundity * fertilization * breeding_success, 1e12)
    
    return recruits

GeoSimulation._reproduce = _patched_reproduce


def get_total_pop(result, sites):
    """Get total population trajectory."""
    n_steps = len(result.states)
    total = np.zeros(n_steps)
    for t in range(n_steps):
        total[t] = sum(result.states[t].populations)
    return total


def get_region_pops(result, sites):
    """Get population trajectories by region."""
    regions = {}
    for i, site in enumerate(sites):
        if site.region not in regions:
            regions[site.region] = []
        regions[site.region].append(i)
    
    n_steps = len(result.states)
    trajs = {}
    for region, indices in regions.items():
        traj = np.zeros(n_steps)
        for t in range(n_steps):
            traj[t] = sum(result.states[t].populations[i] for i in indices)
        trajs[region] = traj
    return trajs


def run_sweep():
    fecundity_values = [100, 1000, 10000, 100000]
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures', 'fecundity_sweep')
    os.makedirs(outdir, exist_ok=True)
    
    sites = list(ALL_SITES)
    
    # --- Part 1: Fecundity sweep (no disease, 20 years to check equilibrium) ---
    print("=" * 70)
    print("FECUNDITY SWEEP — Pre-disease stability check (20 years, no disease)")
    print("=" * 70)
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes_flat = axes.flatten()
    
    for idx, fec in enumerate(fecundity_values):
        config = GeoConfig(n_years=20, disease_onset_year=999)  # No disease
        config._fecundity = fec
        
        sim = GeoSimulation(config=config, sites=sites, seed=42)
        result = sim.run()
        
        total = get_total_pop(result, sites)
        # Convert from ×1000 to billions
        total_B = total / 1e6
        
        ax = axes_flat[idx]
        years = np.arange(len(total)) / 4  # quarterly steps to years
        ax.plot(years, total_B, 'b-', linewidth=2)
        ax.axhline(y=6.1, color='red', linestyle='--', alpha=0.5, label='Target: 6.1B')
        ax.set_xlabel('Year')
        ax.set_ylabel('Total Population (billions)')
        ax.set_title(f'Fecundity = {fec:,}')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        print(f"\nFecundity = {fec:,}")
        print(f"  Y0: {total_B[0]:.2f}B")
        print(f"  Y5: {total_B[20]:.2f}B")
        print(f"  Y10: {total_B[40]:.2f}B")
        print(f"  Y20: {total_B[-1]:.2f}B")
        print(f"  Change Y0→Y10: {(total_B[40]/total_B[0] - 1)*100:+.1f}%")
    
    fig.suptitle('Pre-Disease Population Stability by Fecundity\n(recruits per breeding pair, SRS=0.001)', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, '01_fecundity_stability.png'), dpi=150)
    plt.close()
    print(f"\nSaved: {outdir}/01_fecundity_stability.png")
    
    # --- Part 2: Full 80-year runs with disease ---
    print("\n" + "=" * 70)
    print("FULL 80-YEAR RUNS WITH DISEASE")
    print("=" * 70)
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes_flat = axes.flatten()
    
    for idx, fec in enumerate(fecundity_values):
        config = GeoConfig(n_years=80)
        config._fecundity = fec
        
        sim = GeoSimulation(config=config, sites=sites, seed=42)
        result = sim.run()
        
        # Region trajectories
        region_pops = get_region_pops(result, sites)
        ax = axes_flat[idx]
        years = np.arange(len(result.states)) / 4
        
        for region in REGION_ORDER:
            if region not in region_pops:
                continue
            traj = region_pops[region]
            init = traj[0]
            if init == 0:
                continue
            ratio = np.maximum(traj / init, 1e-5)
            color = REGION_COLORS.get(region, '#333')
            label = REGION_LABELS.get(region, region)
            ax.semilogy(years, ratio, color=color, linewidth=1.5, label=label)
        
        ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)
        ax.axhline(y=0.1, color='gray', linestyle=':', alpha=0.5)
        ax.axvline(x=10, color='red', linestyle=':', alpha=0.3)
        ax.set_xlabel('Year')
        ax.set_ylabel('Population Ratio')
        ax.set_title(f'Fecundity = {fec:,}')
        ax.legend(fontsize=7, loc='lower left')
        ax.set_ylim(1e-5, 10)
        ax.grid(True, alpha=0.3)
        
        # Print Y17 stats
        y17_step = 17 * 4
        print(f"\nFecundity = {fec:,} — Y17 (2020) decline:")
        for region in REGION_ORDER:
            if region not in region_pops:
                continue
            init = region_pops[region][0]
            y17 = region_pops[region][y17_step]
            if init > 0:
                decline = (1 - y17/init) * 100
                label = REGION_LABELS.get(region, region)
                print(f"  {label:<12}: {decline:.1f}%")
    
    fig.suptitle('80-Year Regional Trajectories by Fecundity\n(with SSWD, log scale)', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, '02_fecundity_disease.png'), dpi=150)
    plt.close()
    print(f"\nSaved: {outdir}/02_fecundity_disease.png")
    
    # --- Part 3: Wild WA broodstock ---
    # Approach: Monkey-patch _simulate_season to inject outplanting
    print("\n" + "=" * 70)
    print("WILD WASHINGTON BROODSTOCK SCENARIOS")
    print("=" * 70)
    
    wa_or_indices = [i for i, s in enumerate(sites) if s.region == "wa_or_outer"]
    salish_indices = [i for i, s in enumerate(sites) if s.region == "salish_sea"]
    
    _original_simulate = GeoSimulation._simulate_season
    
    def _outplant_simulate(self, step, year, season):
        """Simulate + inject outplanting."""
        state = _original_simulate(self, step, year, season)
        outplant_amount = getattr(self.config, '_outplant_amount', 0)
        outplant_indices = getattr(self.config, '_outplant_indices', [])
        outplant_start = getattr(self.config, '_outplant_start_year', 999)
        
        # Outplant every spring after start year
        if year >= outplant_start and season == 1 and outplant_amount > 0:
            for i in outplant_indices:
                self.populations[i] += outplant_amount
                self.juveniles[i] += outplant_amount
        return state
    
    # Test fecundity=1000 (closest to Hamilton), with WA broodstock
    fec = 1000
    outplant_amounts = [10, 100, 1000]  # ×1000 = 10K, 100K, 1M per site per year
    
    # Baseline (no outplanting)
    config_base = GeoConfig(n_years=80)
    config_base._fecundity = fec
    sim_base = GeoSimulation(config=config_base, sites=sites, seed=42)
    result_base = sim_base.run()
    
    # Outplanting runs
    GeoSimulation._simulate_season = _outplant_simulate
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    for ax_idx, outplant_amt in enumerate(outplant_amounts):
        config_op = GeoConfig(n_years=80)
        config_op._fecundity = fec
        config_op._outplant_amount = outplant_amt
        config_op._outplant_indices = wa_or_indices + salish_indices
        config_op._outplant_start_year = 17
        
        sim_op = GeoSimulation(config=config_op, sites=sites, seed=42)
        result_op = sim_op.run()
        
        ax = axes[ax_idx]
        years = np.arange(len(result_base.states)) / 4
        
        for region, indices, color, label in [
            ("wa_or_outer", wa_or_indices, "#795548", "WA/OR"),
            ("salish_sea", salish_indices, "#9b59b6", "Salish"),
        ]:
            # Baseline
            base_traj = np.array([sum(result_base.states[t].populations[i] for i in indices) for t in range(len(result_base.states))])
            init = base_traj[0]
            if init == 0: continue
            
            # Outplant
            op_traj = np.array([sum(result_op.states[t].populations[i] for i in indices) for t in range(len(result_op.states))])
            
            ax.semilogy(years, np.maximum(base_traj/init, 1e-6), color=color, linestyle='--', alpha=0.5, label=f'{label} baseline')
            ax.semilogy(years, np.maximum(op_traj/init, 1e-6), color=color, linewidth=2, label=f'{label} + outplant')
        
        real_amt = outplant_amt * 1000
        ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.3)
        ax.axvline(x=10, color='red', linestyle=':', alpha=0.3)
        ax.axvline(x=17, color='green', linestyle=':', alpha=0.3)
        ax.set_xlabel('Year')
        ax.set_ylabel('Population Ratio')
        ax.set_title(f'Outplant {real_amt:,}/site/yr')
        ax.legend(fontsize=8)
        ax.set_ylim(1e-6, 10)
        ax.grid(True, alpha=0.3)
        
        # Print stats
        y40 = min(40*4, len(result_op.states)-1)
        for region, indices, label in [("wa_or_outer", wa_or_indices, "WA/OR"), ("salish_sea", salish_indices, "Salish")]:
            init = sum(result_base.states[0].populations[i] for i in indices)
            base_y40 = sum(result_base.states[y40].populations[i] for i in indices)
            op_y40 = sum(result_op.states[y40].populations[i] for i in indices)
            print(f"  {label} @ Y40, outplant={real_amt:,}: baseline={base_y40/init:.4f}, outplant={op_y40/init:.4f}")
    
    fig.suptitle(f'Wild WA Broodstock Outplanting (fecundity={fec:,}, from Y17)', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, '03_wa_broodstock.png'), dpi=150)
    plt.close()
    print(f"Saved: {outdir}/03_wa_broodstock.png")
    
    # Restore
    GeoSimulation._simulate_season = _original_simulate


if __name__ == "__main__":
    run_sweep()
