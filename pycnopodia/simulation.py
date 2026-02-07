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
    """Record of one year's dynamics."""
    year: int
    n_total: int
    n_adults: int
    n_recruits: int
    n_sswd_deaths: int
    n_outplanted: int
    n_breeders_female: int
    n_breeders_male: int
    ne: float
    breeding_fraction: float
    sswd_prevalence: float
    mean_resistance: float
    resistance_allele_freq: float
    broodstock_inbreeding: float = 0.0
    extinct: bool = False


@dataclass
class SimulationResult:
    """Results from a single simulation run."""
    config: Config
    records: List[YearlyRecord] = field(default_factory=list)
    
    @property
    def years(self) -> np.ndarray:
        return np.array([r.year for r in self.records])
    
    @property
    def n_total(self) -> np.ndarray:
        return np.array([r.n_total for r in self.records])
    
    @property
    def n_adults(self) -> np.ndarray:
        return np.array([r.n_adults for r in self.records])
    
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
    def final_population(self) -> int:
        return self.records[-1].n_total if self.records else 0
    
    @property
    def final_resistance_freq(self) -> float:
        return self.records[-1].resistance_allele_freq if self.records else 0.0


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
            Complete simulation results
        """
        result = SimulationResult(config=self.config)
        
        for year in range(self.config.n_years):
            record = self._simulate_year(year)
            result.records.append(record)
            
            # Check for extinction
            if self.population.n == 0:
                record.extinct = True
                # Pad remaining years with extinction records
                for future_year in range(year + 1, self.config.n_years):
                    result.records.append(YearlyRecord(
                        year=future_year,
                        n_total=0, n_adults=0, n_recruits=0,
                        n_sswd_deaths=0, n_outplanted=0,
                        n_breeders_female=0, n_breeders_male=0,
                        ne=0, breeding_fraction=0,
                        sswd_prevalence=calculate_prevalence(future_year, self.config),
                        mean_resistance=0, resistance_allele_freq=0,
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
        
        # 4. Outplanting
        n_outplanted = 0
        broodstock_inbreeding = 0.0
        if self.broodstock is not None:
            n_outplanted = outplant(pop, self.broodstock, year, self.rng)
            broodstock_inbreeding = self.broodstock.mean_inbreeding()
        
        # 5. Aging
        pop.age_population()
        pop.year = year
        
        # Record state
        return YearlyRecord(
            year=year,
            n_total=pop.n,
            n_adults=pop.n_adults,
            n_recruits=n_recruits,
            n_sswd_deaths=n_sswd_deaths,
            n_outplanted=n_outplanted,
            n_breeders_female=n_f_breeders,
            n_breeders_male=n_m_breeders,
            ne=ne,
            breeding_fraction=breeding_frac,
            sswd_prevalence=sswd_prevalence,
            mean_resistance=pop.mean_resistance(),
            resistance_allele_freq=pop.resistance_allele_frequency(),
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


def summarize_replicates(results: List[SimulationResult]) -> Dict:
    """
    Summarize results across multiple replicates.
    
    Returns
    -------
    summary : dict
        Summary statistics
    """
    n = len(results)
    
    extinction_count = sum(1 for r in results if r.extinct)
    extinction_prob = extinction_count / n
    
    final_pops = [r.final_population for r in results]
    final_freqs = [r.final_resistance_freq for r in results]
    
    # Filter out extinct runs for population stats
    surviving = [r for r in results if not r.extinct]
    
    return {
        "n_replicates": n,
        "extinction_probability": extinction_prob,
        "extinction_count": extinction_count,
        "mean_final_population": np.mean(final_pops),
        "std_final_population": np.std(final_pops),
        "mean_final_resistance_freq": np.mean(final_freqs),
        "std_final_resistance_freq": np.std(final_freqs),
        "mean_final_pop_surviving": np.mean([r.final_population for r in surviving]) if surviving else 0,
    }
