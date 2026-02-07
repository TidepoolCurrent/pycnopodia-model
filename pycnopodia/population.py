"""
Population class representing a Pycnopodia population.

Individuals are represented as structured arrays for computational efficiency.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Tuple, Dict
from .config import Config


@dataclass
class Population:
    """
    A population of Pycnopodia individuals.
    
    Attributes
    ----------
    genomes : ndarray, shape (N, n_loci), dtype uint8
        Genotype at each locus (0/1/2 = copies of resistance allele)
    ages : ndarray, shape (N,), dtype int8
        Age in years
    sexes : ndarray, shape (N,), dtype bool
        True = female, False = male
    alive : ndarray, shape (N,), dtype bool
        Whether individual is alive (for vectorized operations)
        
    Baseline tracking (for ratio calculations):
    - N0: Initial population size at simulation start
    - H0: Initial expected heterozygosity at simulation start
    - All outputs expressed as ratios relative to these baselines
    """
    
    config: Config
    genomes: np.ndarray = field(default=None, repr=False)
    ages: np.ndarray = field(default=None, repr=False)
    sexes: np.ndarray = field(default=None, repr=False)
    alive: np.ndarray = field(default=None, repr=False)
    
    # Effect sizes per locus (allows variable effects)
    locus_effects: np.ndarray = field(default=None, repr=False)
    
    # Baseline tracking for ratio calculations
    N0: int = 0  # Initial population size
    H0: float = 0.0  # Initial expected heterozygosity
    Ne_cumulative: float = 0.0  # Harmonic mean Ne tracker
    Ne_years: int = 0  # Years of Ne accumulation
    
    # Tracking
    year: int = 0
    
    def __post_init__(self):
        """Initialize population if arrays not provided."""
        if self.genomes is None:
            self._initialize_population()
        if self.locus_effects is None:
            self._initialize_locus_effects()
        # Set baselines
        if self.N0 == 0:
            self.N0 = self.n
        if self.H0 == 0.0:
            self.H0 = self.expected_heterozygosity()
    
    def _initialize_locus_effects(self):
        """Set up per-locus effect sizes."""
        if self.config.variable_effect_sizes:
            # Sample from gamma distribution with mean = per_allele_effect
            mean = self.config.per_allele_effect
            cv = self.config.effect_size_cv
            shape = 1 / (cv ** 2)
            scale = mean / shape
            self.locus_effects = np.random.gamma(shape, scale, self.config.n_loci)
        else:
            self.locus_effects = np.full(
                self.config.n_loci, 
                self.config.per_allele_effect
            )
    
    def _initialize_population(self):
        """Create initial population with given parameters."""
        n = self.config.initial_population
        n_loci = self.config.n_loci
        q = self.config.initial_resistance_freq
        
        # Initialize genomes: each locus independently in HW equilibrium
        # P(0) = (1-q)^2, P(1) = 2q(1-q), P(2) = q^2
        self.genomes = np.zeros((n, n_loci), dtype=np.uint8)
        for locus in range(n_loci):
            # Sample two alleles per individual
            allele1 = np.random.random(n) < q
            allele2 = np.random.random(n) < q
            self.genomes[:, locus] = allele1.astype(np.uint8) + allele2.astype(np.uint8)
        
        # Age distribution: approximate stable age distribution
        # More juveniles, fewer old adults
        age_probs = self._stable_age_distribution()
        self.ages = np.random.choice(
            len(age_probs), size=n, p=age_probs
        ).astype(np.int8)
        
        # 50/50 sex ratio
        self.sexes = np.random.random(n) < 0.5
        
        # All initially alive
        self.alive = np.ones(n, dtype=bool)
    
    def _stable_age_distribution(self) -> np.ndarray:
        """Calculate stable age distribution from survival rates."""
        max_age = self.config.max_age
        probs = np.zeros(max_age + 1)
        
        # Survivorship to each age
        probs[0] = 1.0  # All start at age 0
        for age in range(1, max_age + 1):
            probs[age] = probs[age - 1] * self._survival_rate(age - 1)
        
        # Normalize
        probs /= probs.sum()
        return probs
    
    def _survival_rate(self, age: int) -> float:
        """Get annual survival probability for given age."""
        if age == 0:
            return self.config.survival_juvenile_y0
        elif age < 5:
            return self.config.survival_juvenile_y1_4
        elif age < 19:
            return self.config.survival_adult
        else:
            return self.config.survival_senescent
    
    @property
    def n(self) -> int:
        """Total number of living individuals."""
        return int(self.alive.sum())
    
    @property
    def n_ratio(self) -> float:
        """Population size as ratio of baseline (N/N0)."""
        if self.N0 == 0:
            return 0.0
        return self.n / self.N0
    
    @property
    def n_adults(self) -> int:
        """Number of reproductive adults (age >= maturation_age)."""
        return int((self.alive & (self.ages >= self.config.maturation_age)).sum())
    
    @property
    def n_females(self) -> int:
        """Number of adult females."""
        return int((self.alive & self.sexes & (self.ages >= self.config.maturation_age)).sum())
    
    @property
    def n_males(self) -> int:
        """Number of adult males."""
        return int((self.alive & ~self.sexes & (self.ages >= self.config.maturation_age)).sum())
    
    def resistance_scores(self) -> np.ndarray:
        """
        Calculate resistance score for each individual.
        
        Returns array of mortality reduction factors (0 to 1).
        Higher = more resistant.
        """
        # Sum of (genotype * effect) across loci
        scores = (self.genomes * self.locus_effects).sum(axis=1)
        return np.clip(scores, 0, 1)  # Cap at 100% protection
    
    def mean_resistance(self) -> float:
        """Mean resistance score among living individuals."""
        if self.n == 0:
            return 0.0
        return float(self.resistance_scores()[self.alive].mean())
    
    def resistance_allele_frequency(self) -> float:
        """Mean resistance allele frequency across all loci."""
        if self.n == 0:
            return 0.0
        # Mean genotype / 2 = allele frequency
        return float(self.genomes[self.alive].mean() / 2)
    
    def expected_heterozygosity(self) -> float:
        """
        Calculate expected heterozygosity (gene diversity).
        
        He = 2pq averaged across loci, where p = allele frequency.
        This is the standard measure of genetic diversity.
        """
        if self.n == 0:
            return 0.0
        
        # Calculate allele frequency at each locus
        genomes_alive = self.genomes[self.alive]
        p = genomes_alive.mean(axis=0) / 2  # Divide by 2 since genotype is 0/1/2
        
        # He = 2pq = 2p(1-p) at each locus, then average
        he_per_locus = 2 * p * (1 - p)
        return float(he_per_locus.mean())
    
    @property
    def h_ratio(self) -> float:
        """Heterozygosity as ratio of baseline (H/H0)."""
        if self.H0 == 0.0:
            return 0.0
        return self.expected_heterozygosity() / self.H0
    
    def effective_population_size(self, n_breeders_f: int, n_breeders_m: int) -> float:
        """
        Calculate effective population size from breeder counts.
        
        Ne = 4 * Nf * Nm / (Nf + Nm)
        """
        if n_breeders_f + n_breeders_m == 0:
            return 0.0
        return 4 * n_breeders_f * n_breeders_m / (n_breeders_f + n_breeders_m)
    
    def ne_ratio(self, ne: float) -> float:
        """
        Calculate Ne/N ratio (genetic effective vs census).
        
        This ratio captures the SRS effect - lower ratios mean
        more reproductive skew and faster genetic drift.
        """
        if self.n == 0:
            return 0.0
        return ne / self.n
    
    def update_cumulative_ne(self, ne: float):
        """
        Update cumulative Ne for long-term effective size calculation.
        
        Long-term Ne is the harmonic mean across generations.
        """
        if ne > 0:
            self.Ne_cumulative += 1.0 / ne
            self.Ne_years += 1
    
    @property
    def long_term_ne(self) -> float:
        """
        Harmonic mean Ne across all years.
        
        This is the appropriate measure for cumulative genetic drift.
        """
        if self.Ne_years == 0 or self.Ne_cumulative == 0:
            return 0.0
        return self.Ne_years / self.Ne_cumulative
    
    def apply_mortality(self, mortality_probs: np.ndarray):
        """
        Apply mortality to population.
        
        Parameters
        ----------
        mortality_probs : ndarray, shape (N,)
            Per-individual mortality probability
        """
        deaths = np.random.random(len(self.alive)) < mortality_probs
        self.alive &= ~deaths
    
    def age_population(self):
        """Increment ages by one year."""
        self.ages[self.alive] += 1
        # Kill individuals exceeding max age
        too_old = self.ages > self.config.max_age
        self.alive &= ~too_old
    
    def apply_natural_mortality(self):
        """Apply age-specific natural mortality."""
        mortality = np.zeros(len(self.alive))
        
        for age in range(self.config.max_age + 1):
            mask = (self.ages == age) & self.alive
            survival = self._survival_rate(age)
            mortality[mask] = 1 - survival
        
        self.apply_mortality(mortality)
    
    def compact(self):
        """Remove dead individuals to free memory."""
        alive_idx = self.alive
        self.genomes = self.genomes[alive_idx]
        self.ages = self.ages[alive_idx]
        self.sexes = self.sexes[alive_idx]
        self.alive = np.ones(self.n, dtype=bool)
    
    def add_individuals(
        self, 
        genomes: np.ndarray, 
        ages: np.ndarray, 
        sexes: np.ndarray
    ):
        """Add new individuals to population."""
        self.genomes = np.vstack([self.genomes, genomes])
        self.ages = np.concatenate([self.ages, ages])
        self.sexes = np.concatenate([self.sexes, sexes])
        self.alive = np.concatenate([self.alive, np.ones(len(ages), dtype=bool)])
    
    def get_summary(self) -> Dict:
        """Get summary statistics (ratios are primary, counts secondary)."""
        return {
            # Primary outputs: RATIOS
            "n_ratio": self.n_ratio,
            "h_ratio": self.h_ratio,
            "resistance_allele_freq": self.resistance_allele_frequency(),
            "mean_resistance": self.mean_resistance(),
            
            # Secondary: raw counts (for internal use)
            "_n_total": self.n,
            "_n_adults": self.n_adults,
            "_n_females": self.n_females,
            "_n_males": self.n_males,
            "_mean_age": float(self.ages[self.alive].mean()) if self.n > 0 else 0,
            
            # Baselines (for reference)
            "_N0": self.N0,
            "_H0": self.H0,
        }
