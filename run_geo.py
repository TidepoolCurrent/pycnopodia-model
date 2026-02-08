#!/usr/bin/env python3
"""Run geography-based Pycnopodia simulation with ensemble and figures."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pycnopodia.geo_model import GeoSimulation, GeoConfig, run_geo_ensemble
from data.real_sites import ALL_SITES
import argparse


REGION_ORDER = [
    "se_alaska_north", "se_alaska_south", "bc_fjords", "bc_outer",
    "salish_sea", "wa_or_outer", "n_california", "c_california", "s_california"
]

REGION_LABELS = {
    "se_alaska_north": "SE AK N",
    "se_alaska_south": "SE AK S", 
    "bc_fjords": "BC Fjords",
    "bc_outer": "BC Coast",
    "salish_sea": "Salish Sea",
    "wa_or_outer": "WA/OR Coast",
    "n_california": "N. CA",
    "c_california": "C. CA",
    "s_california": "S. CA",
}

REGION_COLORS = {
    "s_california": "#e74c3c",
    "c_california": "#e67e22",
    "n_california": "#e91e95",
    "wa_or_outer": "#795548",
    "salish_sea": "#9b59b6",
    "bc_outer": "#27ae60",
    "bc_fjords": "#00bcd4",
    "se_alaska_south": "#2196f3",
    "se_alaska_north": "#00e5ff",
}


def get_region_trajectories(result, sites):
    """Extract population trajectories by region."""
    regions = {}
    for i, site in enumerate(sites):
        if site.region not in regions:
            regions[site.region] = {"indices": [], "init_total": 0}
        regions[site.region]["indices"].append(i)
        regions[site.region]["init_total"] += result.states[0].populations[i]
    
    trajectories = {}
    for region, info in regions.items():
        traj = []
        for state in result.states:
            total = sum(state.populations[i] for i in info["indices"])
            traj.append(total / info["init_total"] if info["init_total"] > 0 else 0)
        trajectories[region] = np.array(traj)
    
    return trajectories


def plot_ensemble_trajectories(results, sites, outdir):
    """Plot ensemble trajectories with confidence bands."""
    n_years = len(results[0].states)
    years = np.arange(n_years)
    
    # Collect trajectories per region
    region_trajs = {r: [] for r in REGION_ORDER}
    for result in results:
        trajs = get_region_trajectories(result, sites)
        for r in REGION_ORDER:
            if r in trajs:
                region_trajs[r].append(trajs[r])
    
    fig, ax = plt.subplots(1, 1, figsize=(14, 8))
    
    for region in REGION_ORDER:
        if not region_trajs[region]:
            continue
        data = np.array(region_trajs[region])
        mean = data.mean(axis=0)
        lo = np.percentile(data, 2.5, axis=0)
        hi = np.percentile(data, 97.5, axis=0)
        
        color = REGION_COLORS.get(region, "#333333")
        label = REGION_LABELS.get(region, region)
        
        # Clip for log scale
        mean_plot = np.maximum(mean, 1e-5)
        lo_plot = np.maximum(lo, 1e-5)
        hi_plot = np.maximum(hi, 1e-5)
        
        ax.semilogy(years, mean_plot, color=color, linewidth=2, label=label)
        ax.fill_between(years, lo_plot, hi_plot, color=color, alpha=0.15)
    
    ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='Baseline')
    ax.axhline(y=0.1, color='gray', linestyle=':', alpha=0.5, label='10% threshold')
    ax.axvline(x=10, color='red', linestyle=':', alpha=0.3, label='Disease onset')
    
    ax.set_xlabel('Year', fontsize=12)
    ax.set_ylabel('Population Ratio (relative to initial)', fontsize=12)
    ax.set_title(f'Geo-Model Ensemble Trajectories (n={len(results)}, mean ± 95% CI)', fontsize=14)
    ax.legend(loc='upper right', fontsize=9)
    ax.set_ylim(1e-5, 3)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'geo_ensemble_trajectories.png'), dpi=150)
    plt.close()
    print(f"Saved: {outdir}/geo_ensemble_trajectories.png")


def plot_fjord_details(results, sites, outdir):
    """Plot individual fjord site trajectories."""
    fjord_sites = [(i, s) for i, s in enumerate(sites) 
                   if s.site_type == "fjord" and s.sill_depth_m and s.sill_depth_m < 50]
    
    if not fjord_sites:
        return
    
    n_years = len(results[0].states)
    years = np.arange(n_years)
    
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    # Population trajectories
    ax = axes[0]
    for idx, site in fjord_sites:
        trajs = []
        for result in results:
            traj = [result.states[y].populations[idx] / result.states[0].populations[idx] 
                    for y in range(n_years)]
            trajs.append(traj)
        data = np.array(trajs)
        mean = data.mean(axis=0)
        mean_plot = np.maximum(mean, 1e-5)
        
        marker = '🌊' if site.has_freshwater_lens else ''
        ax.semilogy(years, mean_plot, linewidth=1.5, 
                     label=f"{site.name} (sill={site.sill_depth_m}m) {marker}")
    
    ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(y=0.1, color='gray', linestyle=':', alpha=0.5)
    ax.axvline(x=10, color='red', linestyle=':', alpha=0.3)
    ax.set_ylabel('Population Ratio')
    ax.set_title(f'Shallow-Sill Fjord Sites (n={len(results)} ensemble)', fontsize=13)
    ax.legend(fontsize=7, loc='lower right', ncol=2)
    ax.set_ylim(1e-3, 3)
    ax.grid(True, alpha=0.3)
    
    # Resistance evolution
    ax = axes[1]
    for idx, site in fjord_sites:
        trajs = []
        for result in results:
            traj = [result.states[y].resistance_freqs[idx].mean() for y in range(n_years)]
            trajs.append(traj)
        data = np.array(trajs)
        mean = data.mean(axis=0)
        
        ax.plot(years, mean, linewidth=1.5, label=f"{site.name}")
    
    ax.axvline(x=10, color='red', linestyle=':', alpha=0.3)
    ax.set_xlabel('Year')
    ax.set_ylabel('Mean Resistance Allele Frequency')
    ax.set_title('Resistance Evolution in Fjord Refugia', fontsize=13)
    ax.legend(fontsize=7, loc='upper left', ncol=2)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'geo_fjord_details.png'), dpi=150)
    plt.close()
    print(f"Saved: {outdir}/geo_fjord_details.png")


def plot_map(results, sites, outdir):
    """Plot geographic map of final population ratios."""
    n_years = len(results[0].states)
    
    # Average final state across ensemble
    final_ratios = np.zeros(len(sites))
    for result in results:
        for i in range(len(sites)):
            init = result.states[0].populations[i]
            final = result.states[-1].populations[i]
            final_ratios[i] += (final / init if init > 0 else 0)
    final_ratios /= len(results)
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 14))
    
    lats = [s.lat for s in sites]
    lons = [s.lon for s in sites]
    
    # Size proportional to population ratio
    sizes = np.maximum(final_ratios * 200, 5)
    
    # Color by ratio
    scatter = ax.scatter(lons, lats, c=final_ratios, s=sizes, 
                         cmap='RdYlGn', vmin=0, vmax=1, 
                         edgecolors='black', linewidth=0.5, alpha=0.8)
    
    # Label fjord refugia
    for i, site in enumerate(sites):
        if site.site_type == "fjord" and final_ratios[i] > 0.3:
            ax.annotate(site.name, (lons[i], lats[i]), fontsize=6,
                       xytext=(5, 5), textcoords='offset points')
    
    plt.colorbar(scatter, ax=ax, label='Final Population Ratio (Year 79)', shrink=0.5)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(f'Pycnopodia Recovery Map (Year {n_years-1}, n={len(results)} ensemble)')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'geo_map.png'), dpi=150)
    plt.close()
    print(f"Saved: {outdir}/geo_map.png")


def export_dashboard_json(result, sites, config, output_path="dashboard/data.json"):
    """Export simulation result to dashboard JSON format."""
    import json
    n_sites = len(sites)
    n_steps = len(result.states)
    seasons = ["Winter", "Spring", "Summer", "Fall"]
    start_year = 2003  # Y0=2003, so Y10=2013 (SSWD onset)
    
    # Step labels
    step_labels = []
    for step in range(n_steps):
        year = start_year + step // 4
        season = seasons[step % 4]
        step_labels.append(f"{year} {season}")
    
    # Site info
    site_data = []
    for s in sites:
        site_data.append({
            "name": s.name, "lat": s.lat, "lon": s.lon,
            "region": s.region, "site_type": s.site_type,
            "sill_depth": s.sill_depth_m, "has_freshwater_lens": s.has_freshwater_lens,
            "base_temp": s.base_temp_C
        })
    
    # Site timeseries
    site_timeseries = []
    for i in range(n_sites):
        pops = [float(result.states[t].populations[i]) for t in range(n_steps)]
        prevs = [float(result.states[t].disease_prevalence[i]) for t in range(n_steps)]
        res = [float(result.states[t].resistance_freqs[i].mean()) for t in range(n_steps)]
        site_timeseries.append({"population": pops, "prevalence": prevs, "resistance": res})
    
    # Region timeseries
    regions = sorted(set(s.region for s in sites))
    region_timeseries = {}
    for region in regions:
        ridx = [i for i, s in enumerate(sites) if s.region == region]
        pops = [sum(result.states[t].populations[i] for i in ridx) for t in range(n_steps)]
        prevs = [float(np.mean([result.states[t].disease_prevalence[i] for i in ridx])) for t in range(n_steps)]
        res_vals = [float(np.mean([result.states[t].resistance_freqs[i].mean() for i in ridx])) for t in range(n_steps)]
        region_timeseries[region] = {"population": pops, "prevalence": prevs, "resistance": res_vals}
    
    # Connectivity (sparse - only links > 0.001)
    from pycnopodia.geo_model import GeoSimulation
    sim = GeoSimulation(config=config, sites=sites, seed=0)
    larval = sim.larval_connectivity.tolist()
    disease = sim.disease_connectivity.tolist()
    
    data = {
        "model": "geo_model",
        "n_years": config.n_years,
        "n_sites": n_sites,
        "seasons_per_year": 4,
        "total_steps": n_steps,
        "sites": site_data,
        "regions": regions,
        "region_timeseries": region_timeseries,
        "site_timeseries": site_timeseries,
        "larval_connectivity": larval,
        "disease_connectivity": disease,
        "step_labels": step_labels,
        "years": list(range(n_steps))
    }
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(data, f)
    print(f"Dashboard data exported to {output_path} ({os.path.getsize(output_path)/1e6:.1f} MB)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ensemble', type=int, default=5)
    parser.add_argument('--years', type=int, default=80)
    parser.add_argument('--dense', action='store_true', help='Use 200-site dense network')
    parser.add_argument('--dashboard', action='store_true', help='Export dashboard data.json')
    args = parser.parse_args()
    
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'figures', 'geo')
    os.makedirs(outdir, exist_ok=True)
    
    config = GeoConfig(n_years=args.years)
    if args.dense:
        from data.coastline import generate_full_site_network
        sites = generate_full_site_network()
    else:
        sites = list(ALL_SITES)
    
    print(f"Running {args.ensemble}-seed ensemble ({args.years} years, {len(sites)} sites)...")
    from pycnopodia.geo_model import GeoSimulation
    results = []
    for seed in range(args.ensemble):
        sim = GeoSimulation(config=config, sites=sites, seed=seed)
        results.append(sim.run())
        print(f"  Run {seed+1}/{args.ensemble} complete")
    
    print("\nGenerating figures...")
    plot_ensemble_trajectories(results, sites, outdir)
    plot_fjord_details(results, sites, outdir)
    plot_map(results, sites, outdir)
    
    # Print summary statistics
    print("\n" + "=" * 70)
    print(f"GEO-MODEL ENSEMBLE STATISTICS (n={args.ensemble})")
    print("=" * 70)
    
    for region in REGION_ORDER:
        region_indices = [i for i, s in enumerate(sites) if s.region == region]
        if not region_indices:
            continue
        
        y17_ratios = []
        final_ratios = []
        for result in results:
            init = sum(result.states[0].populations[i] for i in region_indices)
            y17_step = 17 * config.seasons_per_year  # Year 17, not step 17
            y17 = sum(result.states[y17_step].populations[i] for i in region_indices)
            final = sum(result.states[-1].populations[i] for i in region_indices)
            y17_ratios.append(1 - y17/init if init > 0 else 1)
            final_ratios.append(final/init if init > 0 else 0)
        
        label = REGION_LABELS.get(region, region)
        n_fjords = sum(1 for i in region_indices if sites[i].site_type == "fjord")
        print(f"  {label:<12} ({n_fjords}f): Y17 decline={np.mean(y17_ratios)*100:.1f}%, "
              f"Final ratio={np.mean(final_ratios):.3f}±{np.std(final_ratios):.3f}")
    
    print("=" * 70)
    
    # Export dashboard data (use first run)
    if args.dashboard or True:  # Always export for now
        export_dashboard_json(results[0], sites, config)


if __name__ == "__main__":
    main()
