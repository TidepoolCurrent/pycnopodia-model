#!/usr/bin/env python3
"""
Run ensemble simulation with retuned parameters and check against Hamilton et al. targets.
"""

import sys
sys.path.insert(0, '.')

from run_pacific_coast import run_ensemble, plot_ensemble_trajectories, plot_ensemble_summary
from pycnopodia.pacific_coast import PacificCoastConfig, REGION_ORDER

def main():
    print("Running ensemble with retuned parameters (n=20)...")
    print("=" * 70)
    
    config = PacificCoastConfig(n_years=100)
    ensemble = run_ensemble(n_runs=20, config=config, verbose=True)
    
    print("\n" + "=" * 70)
    print("ENSEMBLE RESULTS VS HAMILTON ET AL. 2021 TARGETS")
    print("=" * 70)
    
    # Get final state (year ~10 = 2020 in model time)
    year_2020_idx = 10  # Year 10 = 2020 (onset was 2013, year 10 in model)
    
    print(f"\n{'Region':<25} {'Mean Survival':<15} {'Target':<20} {'Status':<10}")
    print("-" * 70)
    
    targets = {
        "se_alaska_north": (0.02, 0.08, "SE AK: 96% decline (~4% remain)"),
        "se_alaska_south": (0.01, 0.04, "SE AK: 96% decline (~4% remain)"),
        "bc_outer": (0.08, 0.15, "BC: 87.9% decline (~12% remain)"),
        "bc_fjords": (0.08, 0.18, "BC: 87.9% decline (~12% remain)"),
        "salish_sea": (0.05, 0.10, "Salish: 92.4% decline (~7.6%)"),
        "wa_or_outer": (0.00, 0.02, ">99.2% decline"),
        "n_california": (0.00, 0.02, ">99.2% decline"),
        "c_california": (0.00, 0.02, ">99.2% decline"),
        "s_california": (0.00, 0.01, ">99.2% decline"),
    }
    
    for region_id in REGION_ORDER:
        mean_surv = ensemble.mean_trajectories[region_id]["population_ratio"][year_2020_idx]
        region_name = config.regions[region_id].short_name
        
        if region_id in targets:
            low, high, desc = targets[region_id]
            in_range = low <= mean_surv <= high
            status = "✓" if in_range else "✗"
            target_str = f"{low:.1%}-{high:.1%}"
        else:
            target_str = "N/A"
            status = ""
        
        print(f"{region_name:<25} {mean_surv:>13.1%}  {target_str:<20} {status:<10}")
    
    # Overall Pacific Coast
    print("-" * 70)
    overall_initial = sum(ensemble.all_runs[0].states[0].get_region_summary(
        ensemble.all_runs[0].sites, config)[r]["population"] for r in REGION_ORDER)
    overall_2020 = sum(ensemble.mean_trajectories[r]["population_ratio"][year_2020_idx] * 
                      ensemble.all_runs[0].states[0].get_region_summary(
                          ensemble.all_runs[0].sites, config)[r]["population"] 
                      for r in REGION_ORDER)
    overall_surv = overall_2020 / overall_initial
    print(f"{'OVERALL PACIFIC COAST':<25} {overall_surv:>13.1%}  {'5-10% (94% decline)':<20} {'✓' if 0.05 <= overall_surv <= 0.12 else '✗':<10}")
    
    print("\n" + "=" * 70)
    print("GENERATING PLOTS...")
    print("=" * 70)
    
    plot_ensemble_trajectories(ensemble)
    plot_ensemble_summary(ensemble)
    
    print("\n✓ Ensemble analysis complete!")
    print("  - Figures saved to figures/pacific_coast/")

if __name__ == "__main__":
    main()
