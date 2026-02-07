"""
Main simulation class orchestrating all model components.
"""

import numpy as np
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field
from .config import Config
from .population import Population
from .reproduction import reproduce
from .disease import apply_sswd, calculate_prevalence
from .intervention import Broodstock, outplant, should_outplant


@dataclass
class YearlyRecord:
    """
    Record of one year's dynamics.
    
    Primary outputs are RATIOS (meaningful across scenarios):
    - n_ratio: Population as fraction of baseline (N/N0)
    - ne_ratio: Effective/census ratio (Ne/N) 
    - h_ratio: Heterozygosity retention (H/H0)
    
    Absolute counts prefixed with _ are for internal use only.
    """
    year: int
    
    # PRIMARY: Ratios (comparable across scenarios)
    n_ratio: float           # N/N0 - population relative to baseline
    ne_ratio: float          # Ne/N - effective vs census size
    h_ratio: float           # H/H0 - genetic diversity retention
    breeding_fraction: float # Fraction of adults that bred
    resistance_allele_freq: float
    mean_resistance: float
    sswd_prevalence: float
    
    # SECONDARY: Absolute counts (internal use, prefixed)
    _n_total: int = 0
    _n_adults: int = 0
    _n_recruits: int = 0
    _n_sswd_deaths: int = 0
    _n_outplanted: int = 0
    _n_breeders_female: int = 0
    _n_breeders_male: int = 0
    _ne: float = 0.0
    
    broodstock_inbreeding: float = 0.0
    extinct: bool = False


@dataclass
class SimulationResult:
    """
    Results from a single simulation run.
    
    Primary outputs are ratio trajectories (N/N0, Ne/N, H/H0).
    These are comparable across scenarios regardless of absolute numbers.
    """
    config: Config
    records: List[YearlyRecord] = field(default_factory=list)
    
    # Baselines stored for reference
    N0: int = 0
    H0: float = 0.0
    
    @property
    def years(self) -> np.ndarray:
        return np.array([r.year for r in self.records])
    
    # PRIMARY: Ratio trajectories
    @property
    def n_ratio(self) -> np.ndarray:
        """Population trajectory as fraction of baseline."""
        return np.array([r.n_ratio for r in self.records])
    
    @property
    def ne_ratio(self) -> np.ndarray:
        """Effective/census ratio trajectory."""
        return np.array([r.ne_ratio for r in self.records])
    
    @property
    def h_ratio(self) -> np.ndarray:
        """Genetic diversity retention trajectory."""
        return np.array([r.h_ratio for r in self.records])
    
    @property
    def resistance_freq(self) -> np.ndarray:
        """Resistance allele frequency trajectory."""
        return np.array([r.resistance_allele_freq for r in self.records])
    
    @property
    def breeding_fraction(self) -> np.ndarray:
        """Breeding fraction trajectory."""
        return np.array([r.breeding_fraction for r in self.records])
    
    # Recovery metrics (ratio-based)
    def years_to_recovery(self, threshold: float = 0.30) -> Optional[int]:
        """Years to reach threshold fraction of baseline (default 30%)."""
        for r in self.records:
            if r.n_ratio >= threshold:
                return r.year
        return None
    
    def recovered(self, threshold: float = 0.30) -> bool:
        """Whether population reached recovery threshold."""
        return self.years_to_recovery(threshold) is not None
    
    # Status
    @property
    def extinct(self) -> bool:
        return len(self.records) > 0 and self.records[-1].extinct
    
    @property
    def extinction_year(self) -> Optional[int]:
        for r in self.records:
            if r.extinct:
                return r.year
        return None
    
    @property
    def final_n_ratio(self) -> float:
        """Final population as fraction of baseline."""
        return self.records[-1].n_ratio if self.records else 0.0
    
    @property
    def final_h_ratio(self) -> float:
        """Final heterozygosity as fraction of baseline."""
        return self.records[-1].h_ratio if self.records else 0.0
    
    @property
    def final_resistance_freq(self) -> float:
        return self.records[-1].resistance_allele_freq if self.records else 0.0
    
    @property
    def min_n_ratio(self) -> float:
        """Minimum population ratio (bottleneck depth)."""
        if not self.records:
            return 0.0
        return min(r.n_ratio for r in self.records)
    
    @property 
    def bottleneck_year(self) -> int:
        """Year of minimum population."""
        if not self.records:
            return 0
        return min(self.records, key=lambda r: r.n_ratio).year


class Simulation:
    """
    Main simulation class.
    
    Orchestrates population dynamics, disease, reproduction, and intervention.
    """
    
    def __init__(
        self,
        config: Config = None,
        seed: int = None
    ):
        """
        Initialize simulation.
        
        Parameters
        ----------
        config : Config, optional
            Model configuration. Uses defaults if not provided.
        seed : int, optional
            Random seed for reproducibility.
        """
        self.config = config or Config()
        self.config.validate()
        
        # Set up RNG
        if seed is not None:
            self.rng = np.random.default_rng(seed)
        elif self.config.random_seed is not None:
            self.rng = np.random.default_rng(self.config.random_seed)
        else:
            self.rng = np.random.default_rng()
        
        # Initialize population
        self.population = Population(self.config)
        
        # Initialize broodstock if outplanting is configured
        self.broodstock = None
        if self.config.outplanting_n_per_year > 0:
            self.broodstock = Broodstock(self.config)
    
    def run(self) -> SimulationResult:
        """
        Run full simulation.
        
        Returns
        -------
        result : SimulationResult
            Complete simulation results with ratio-based outputs
        """
        # Store baselines in result
        result = SimulationResult(
            config=self.config,
            N0=self.population.N0,
            H0=self.population.H0
        )
        
        for year in range(self.config.n_years):
            record = self._simulate_year(year)
            result.records.append(record)
            
            # Check for extinction
            if self.population.n == 0:
                record.extinct = True
                # Pad remaining years with extinction records (all ratios = 0)
                for future_year in range(year + 1, self.config.n_years):
                    result.records.append(YearlyRecord(
                        year=future_year,
                        n_ratio=0.0,
                        ne_ratio=0.0,
                        h_ratio=0.0,
                        breeding_fraction=0.0,
                        resistance_allele_freq=0.0,
                        mean_resistance=0.0,
                        sswd_prevalence=calculate_prevalence(future_year, self.config),
                        _n_total=0, _n_adults=0, _n_recruits=0,
                        _n_sswd_deaths=0, _n_outplanted=0,
                        _n_breeders_female=0, _n_breeders_male=0,
                        _ne=0.0,
                        extinct=True
                    ))
                break
            
            # Compact population periodically to save memory
            if year % 10 == 0:
                self.population.compact()
        
        return result
    
    def _simulate_year(self, year: int) -> YearlyRecord:
        """
        Simulate one year of population dynamics.
        
        Order of operations:
        1. Natural mortality (age-specific)
        2. Disease mortality (SSWD)
        3. Reproduction (SRS, Allee effects)
        4. Outplanting (if configured)
        5. Aging
        
        Returns YearlyRecord with ratio-based primary outputs.
        """
        pop = self.population
        
        # 1. Natural mortality
        pop.apply_natural_mortality()
        
        # 2. Disease mortality
        n_sswd_deaths, sswd_prevalence = apply_sswd(pop, year, self.rng)
        
        # 3. Reproduction
        n_recruits, n_f_breeders, n_m_breeders, breeding_frac = reproduce(pop, self.rng)
        
        # Calculate Ne from breeders
        ne = pop.effective_population_size(n_f_breeders, n_m_breeders)
        
        # Update cumulative Ne for long-term calculation
        pop.update_cumulative_ne(ne)
        
        # 4. Outplanting
        n_outplanted = 0
        broodstock_inbreeding = 0.0
        if self.broodstock is not None:
            n_outplanted = outplant(pop, self.broodstock, year, self.rng)
            broodstock_inbreeding = self.broodstock.mean_inbreeding()
        
        # 5. Aging
        pop.age_population()
        pop.year = year
        
        # Record state with RATIOS as primary outputs
        return YearlyRecord(
            year=year,
            # Primary: Ratios
            n_ratio=pop.n_ratio,
            ne_ratio=pop.ne_ratio(ne),
            h_ratio=pop.h_ratio,
            breeding_fraction=breeding_frac,
            resistance_allele_freq=pop.resistance_allele_frequency(),
            mean_resistance=pop.mean_resistance(),
            sswd_prevalence=sswd_prevalence,
            # Secondary: Counts (prefixed)
            _n_total=pop.n,
            _n_adults=pop.n_adults,
            _n_recruits=n_recruits,
            _n_sswd_deaths=n_sswd_deaths,
            _n_outplanted=n_outplanted,
            _n_breeders_female=n_f_breeders,
            _n_breeders_male=n_m_breeders,
            _ne=ne,
            broodstock_inbreeding=broodstock_inbreeding,
        )


def run_replicates(
    config: Config,
    n_replicates: int = None,
    seeds: List[int] = None,
    progress: bool = True
) -> List[SimulationResult]:
    """
    Run multiple replicate simulations.
    
    Parameters
    ----------
    config : Config
        Model configuration
    n_replicates : int, optional
        Number of replicates. Uses config value if not provided.
    seeds : list of int, optional
        Random seeds for each replicate. Auto-generated if not provided.
    progress : bool
        Whether to print progress
        
    Returns
    -------
    results : list of SimulationResult
        Results from all replicates
    """
    if n_replicates is None:
        n_replicates = config.n_replicates
    
    if seeds is None:
        base_seed = config.random_seed or 42
        seeds = [base_seed + i for i in range(n_replicates)]
    
    results = []
    for i, seed in enumerate(seeds):
        if progress and (i + 1) % 10 == 0:
            print(f"  Replicate {i + 1}/{n_replicates}")
        
        sim = Simulation(config, seed=seed)
        result = sim.run()
        results.append(result)
    
    return results


def summarize_replicates(results: List[SimulationResult], recovery_threshold: float = 0.30) -> Dict:
    """
    Summarize results across multiple replicates.
    
    All population metrics are expressed as RATIOS relative to baseline.
    
    Parameters
    ----------
    results : list of SimulationResult
    recovery_threshold : float
        Fraction of baseline considered "recovered" (default 0.30 = 30%)
    
    Returns
    -------
    summary : dict
        Summary statistics with ratio-based metrics
    """
    n = len(results)
    
    # Extinction
    extinction_count = sum(1 for r in results if r.extinct)
    extinction_prob = extinction_count / n
    
    # Recovery (ratio-based)
    recovery_count = sum(1 for r in results if r.recovered(recovery_threshold))
    recovery_prob = recovery_count / n
    
    # Final state (ratios)
    final_n_ratios = [r.final_n_ratio for r in results]
    final_h_ratios = [r.final_h_ratio for r in results]
    final_freqs = [r.final_resistance_freq for r in results]
    
    # Bottleneck (minimum ratio reached)
    min_ratios = [r.min_n_ratio for r in results]
    bottleneck_years = [r.bottleneck_year for r in results]
    
    # Filter surviving for conditional stats
    surviving = [r for r in results if not r.extinct]
    recovered = [r for r in results if r.recovered(recovery_threshold)]
    
    return {
        "n_replicates": n,
        
        # Outcomes (probabilities)
        "extinction_probability": extinction_prob,
        "recovery_probability": recovery_prob,
        "recovery_threshold": recovery_threshold,
        
        # Final state (RATIOS)
        "mean_final_n_ratio": np.mean(final_n_ratios),
        "std_final_n_ratio": np.std(final_n_ratios),
        "mean_final_h_ratio": np.mean(final_h_ratios),
        "std_final_h_ratio": np.std(final_h_ratios),
        "mean_final_resistance_freq": np.mean(final_freqs),
        "std_final_resistance_freq": np.std(final_freqs),
        
        # Bottleneck severity
        "mean_min_n_ratio": np.mean(min_ratios),
        "mean_bottleneck_year": np.mean(bottleneck_years),
        
        # Conditional on survival
        "mean_final_n_ratio_surviving": np.mean([r.final_n_ratio for r in surviving]) if surviving else 0.0,
        
        # Time to recovery
        "mean_years_to_recovery": np.mean([r.years_to_recovery(recovery_threshold) for r in recovered]) if recovered else None,
    }
