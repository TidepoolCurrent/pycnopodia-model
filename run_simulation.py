#!/usr/bin/env python3
"""
Run Pycnopodia population simulations.

Examples:
    # Run baseline scenario
    python run_simulation.py
    
    # Run with specific preset
    python run_simulation.py --preset metapopulation
    
    # Custom parameters
    python run_simulation.py --n-replicates 100 --n-years 100
"""

import argparse
import json
import numpy as np
from pathlib import Path
from datetime import datetime

from pycnopodia import Simulation, Config
from pycnopodia.config import PRESETS
from pycnopodia.simulation import run_replicates, summarize_replicates


def main():
    parser = argparse.ArgumentParser(description="Run Pycnopodia population model")
    
    parser.add_argument("--preset", choices=list(PRESETS.keys()), default="baseline",
                       help="Configuration preset to use")
    parser.add_argument("--n-replicates", type=int, default=50,
                       help="Number of replicate simulations")
    parser.add_argument("--n-years", type=int, default=80,
                       help="Simulation duration in years")
    parser.add_argument("--seed", type=int, default=42,
                       help="Random seed for reproducibility")
    parser.add_argument("--output", type=str, default="results",
                       help="Output directory")
    parser.add_argument("--no-outplanting", action="store_true",
                       help="Disable outplanting intervention")
    parser.add_argument("--verbose", action="store_true",
                       help="Print detailed progress")
    
    args = parser.parse_args()
    
    # Load preset and override with command-line args
    config = PRESETS[args.preset]
    config.n_replicates = args.n_replicates
    config.n_years = args.n_years
    config.random_seed = args.seed
    
    if args.no_outplanting:
        config.outplanting_n_per_year = 0
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"=== Pycnopodia Population Model ===")
    print(f"Preset: {args.preset}")
    print(f"Replicates: {args.n_replicates}")
    print(f"Years: {args.n_years}")
    print(f"Outplanting: {config.outplanting_n_per_year}/year")
    print()
    
    # Run simulations
    print("Running simulations...")
    results = run_replicates(config, progress=args.verbose)
    
    # Summarize
    summary = summarize_replicates(results)
    
    print()
    print("=== Results ===")
    print(f"Extinction probability: {summary['extinction_probability']:.1%}")
    print(f"Mean final population: {summary['mean_final_population']:.0f} ± {summary['std_final_population']:.0f}")
    print(f"Mean final resistance freq: {summary['mean_final_resistance_freq']:.3f}")
    
    # Save results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Save summary
    summary_file = output_dir / f"summary_{args.preset}_{timestamp}.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSummary saved to {summary_file}")
    
    # Save trajectories
    trajectories_file = output_dir / f"trajectories_{args.preset}_{timestamp}.npz"
    np.savez(
        trajectories_file,
        n_total=np.array([r.n_total for r in results]),
        n_adults=np.array([r.n_adults for r in results]),
        resistance_freq=np.array([[rec.resistance_allele_freq for rec in r.records] for r in results]),
    )
    print(f"Trajectories saved to {trajectories_file}")


if __name__ == "__main__":
    main()
