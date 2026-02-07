#!/usr/bin/env python3
"""
Run Pycnopodia population simulations.

All outputs are expressed as RATIOS relative to baseline:
- N/N₀: Population relative to initial
- Ne/N: Effective vs census size
- H/H₀: Genetic diversity retention

Examples:
    # Run baseline scenario
    python run_simulation.py
    
    # Run with specific preset
    python run_simulation.py --preset metapopulation
    
    # Custom parameters
    python run_simulation.py --n-replicates 100 --n-years 100
    
    # Set recovery threshold
    python run_simulation.py --recovery-threshold 0.50
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
    parser = argparse.ArgumentParser(
        description="Run Pycnopodia population model (ratio-based outputs)"
    )
    
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
    parser.add_argument("--recovery-threshold", type=float, default=0.30,
                       help="Recovery threshold as fraction of baseline (default 0.30)")
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
    
    print(f"╔══════════════════════════════════════════════════════════╗")
    print(f"║         Pycnopodia Population Model (Ratio-Based)        ║")
    print(f"╚══════════════════════════════════════════════════════════╝")
    print()
    print(f"  Preset:             {args.preset}")
    print(f"  Replicates:         {args.n_replicates}")
    print(f"  Years:              {args.n_years}")
    print(f"  Outplanting:        {config.outplanting_n_per_year}/year")
    print(f"  Recovery threshold: {args.recovery_threshold:.0%} of baseline")
    print()
    
    # Run simulations
    print("Running simulations...")
    results = run_replicates(config, progress=args.verbose)
    
    # Summarize with ratio-based metrics
    summary = summarize_replicates(results, recovery_threshold=args.recovery_threshold)
    
    # Print results (RATIO-BASED)
    print()
    print("╔══════════════════════════════════════════════════════════╗")
    print("║                       RESULTS                            ║")
    print("╚══════════════════════════════════════════════════════════╝")
    print()
    print("  OUTCOMES:")
    print(f"    Extinction probability:  {summary['extinction_probability']:.1%}")
    print(f"    Recovery probability:    {summary['recovery_probability']:.1%} (to {args.recovery_threshold:.0%} of baseline)")
    print()
    print("  FINAL STATE (as fraction of baseline):")
    print(f"    Population (N/N₀):       {summary['mean_final_n_ratio']:.2%} ± {summary['std_final_n_ratio']:.2%}")
    print(f"    Diversity (H/H₀):        {summary['mean_final_h_ratio']:.2%} ± {summary['std_final_h_ratio']:.2%}")
    print(f"    Resistance allele freq:  {summary['mean_final_resistance_freq']:.3f} ± {summary['std_final_resistance_freq']:.3f}")
    print()
    print("  BOTTLENECK:")
    print(f"    Minimum (N/N₀):          {summary['mean_min_n_ratio']:.2%}")
    print(f"    Year of minimum:         {summary['mean_bottleneck_year']:.1f}")
    print()
    if summary['mean_years_to_recovery'] is not None:
        print("  RECOVERY:")
        print(f"    Mean years to {args.recovery_threshold:.0%}:       {summary['mean_years_to_recovery']:.1f}")
    print()
    
    # Save results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Save summary
    summary_file = output_dir / f"summary_{args.preset}_{timestamp}.json"
    with open(summary_file, "w") as f:
        # Convert numpy types for JSON serialization
        summary_json = {k: float(v) if isinstance(v, (np.floating, np.integer)) else v 
                       for k, v in summary.items()}
        json.dump(summary_json, f, indent=2)
    print(f"Summary saved to {summary_file}")
    
    # Save ratio trajectories
    trajectories_file = output_dir / f"trajectories_{args.preset}_{timestamp}.npz"
    np.savez(
        trajectories_file,
        # PRIMARY: Ratio trajectories
        n_ratio=np.array([r.n_ratio for r in results]),
        h_ratio=np.array([r.h_ratio for r in results]),
        ne_ratio=np.array([r.ne_ratio for r in results]),
        resistance_freq=np.array([r.resistance_freq for r in results]),
        breeding_fraction=np.array([r.breeding_fraction for r in results]),
        # Baselines
        N0=np.array([r.N0 for r in results]),
        H0=np.array([r.H0 for r in results]),
    )
    print(f"Trajectories saved to {trajectories_file}")
    
    print()
    print("Note: All population metrics are RATIOS relative to baseline.")
    print("      N/N₀=0.50 means 50% of initial population size.")


if __name__ == "__main__":
    main()
