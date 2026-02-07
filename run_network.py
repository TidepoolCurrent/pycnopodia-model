#!/usr/bin/env python3
"""
Run network-based metapopulation simulations.

1000 sites × 1000 individuals = 1 million total population.
Test different connectivity scenarios.
"""

import argparse
import numpy as np
from pycnopodia.network import (
    NetworkConfig, NetworkSimulation, ConnectivityType,
    run_scenario, compare_scenarios
)


def main():
    parser = argparse.ArgumentParser(description="Run network metapopulation model")
    
    parser.add_argument("--scenario", choices=[t.value for t in ConnectivityType],
                       default="stepping",
                       help="Connectivity scenario")
    parser.add_argument("--asymmetry", type=float, default=0.0,
                       help="Flow asymmetry (0=symmetric, 1=one-way)")
    parser.add_argument("--self-recruitment", type=float, default=0.5,
                       help="Fraction staying at natal site")
    parser.add_argument("--n-replicates", type=int, default=10,
                       help="Number of replicates")
    parser.add_argument("--n-years", type=int, default=100,
                       help="Simulation years")
    parser.add_argument("--seed", type=int, default=42,
                       help="Random seed")
    parser.add_argument("--compare", action="store_true",
                       help="Compare all scenarios")
    parser.add_argument("--single", action="store_true",
                       help="Run single simulation with detailed output")
    
    args = parser.parse_args()
    
    print("╔══════════════════════════════════════════════════════════╗")
    print("║      Network Metapopulation Model (1000 sites)           ║")
    print("╚══════════════════════════════════════════════════════════╝")
    print()
    
    if args.compare:
        print("Comparing all connectivity scenarios...")
        print()
        results = compare_scenarios(n_replicates=args.n_replicates, seed=args.seed)
        
        print()
        print("═" * 70)
        print("SUMMARY: Best scenarios by final population ratio")
        print("═" * 70)
        sorted_results = sorted(results, key=lambda x: x['mean_final_n_ratio'], reverse=True)
        for r in sorted_results[:10]:
            print(f"  {r['connectivity_type']:12} asym={r['asymmetry']:.1f} self={r['self_recruitment']:.1f} → "
                  f"N/N₀={r['mean_final_n_ratio']:6.1%} extinct={r['extinction_probability']:.0%}")
        
    elif args.single:
        # Single detailed run
        conn_type = ConnectivityType(args.scenario)
        config = NetworkConfig(
            connectivity_type=conn_type,
            asymmetry=args.asymmetry,
            self_recruitment=args.self_recruitment,
            n_years=args.n_years,
        )
        
        print(f"  Scenario:         {conn_type.value}")
        print(f"  Asymmetry:        {args.asymmetry}")
        print(f"  Self-recruitment: {args.self_recruitment}")
        print(f"  Sites:            {config.n_sites}")
        print(f"  Initial N:        {config.total_initial_population:,}")
        print()
        
        sim = NetworkSimulation(config, seed=args.seed)
        result = sim.run()
        
        print("TRAJECTORY (every 10 years):")
        print("-" * 60)
        print(f"{'Year':>6} {'N/N₀':>10} {'Sites':>10} {'Mean Site':>12}")
        print("-" * 60)
        
        for state in result.states[::10]:
            print(f"{state.year:>6} {state.n_ratio:>10.1%} "
                  f"{state.occupied_ratio:>10.1%} {state.mean_site_n_ratio:>12.1%}")
        
        print("-" * 60)
        print()
        print("FINAL STATE:")
        print(f"  Population (N/N₀):     {result.final_n_ratio:.1%}")
        print(f"  Bottleneck (min N/N₀): {result.min_n_ratio:.1%} (year {result.bottleneck_year})")
        print(f"  Extinct:               {result.extinct}")
        print(f"  Recovered (>30%):      {result.recovered(0.3)}")
        
    else:
        # Run scenario with replicates
        conn_type = ConnectivityType(args.scenario)
        
        print(f"  Scenario:         {conn_type.value}")
        print(f"  Asymmetry:        {args.asymmetry}")
        print(f"  Self-recruitment: {args.self_recruitment}")
        print(f"  Replicates:       {args.n_replicates}")
        print()
        
        result = run_scenario(
            connectivity_type=conn_type,
            asymmetry=args.asymmetry,
            self_recruitment=args.self_recruitment,
            n_replicates=args.n_replicates,
            n_years=args.n_years,
            seed=args.seed,
        )
        
        print("RESULTS:")
        print(f"  Extinction probability: {result['extinction_probability']:.1%}")
        print(f"  Final N/N₀:            {result['mean_final_n_ratio']:.1%} ± {result['std_final_n_ratio']:.1%}")
        print(f"  Bottleneck (min N/N₀): {result['mean_min_n_ratio']:.1%}")
        print(f"  Bottleneck year:       {result['mean_bottleneck_year']:.1f}")


if __name__ == "__main__":
    main()
