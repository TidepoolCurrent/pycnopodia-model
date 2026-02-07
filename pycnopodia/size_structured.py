"""
Size-structured population model.

Pycnopodia biology is fundamentally size-dependent:
- Fecundity scales with size (~r^2.5)
- SSWD mortality is size-dependent (larger more susceptible)
- Growth rate varies with food availability
- Maturation is size-based, not age-based
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Tuple, Dict
from .config import Config


@dataclass
class SizeStructuredConfig:
    """Configuration for size-structured model."""
    
    # Size parameters (arm radius in cm)
    min_size: float = 1.0  # Post-settlement juvenile
    max_size: float = 60.0  # Maximum observed size
    maturation_size: float = 20.0  # Size at first reproduction
    
    # Growth
    von_bert_linf: float = 55.0  # Asymptotic size (cm)
    von_bert_k: float = 0.15  # Growth rate (per year)
    growth_cv: float = 0.2  # Individual variation in growth
    
    # Size-dependent survival
    juvenile_survival_intercept: float = 0.3  # Survival at min size
    adult_survival_max: float = 0.92  # Maximum survival (at optimal size)
    senescence_size: float = 50.0  # Size where senescence begins
    senescence_rate: float = 0.02  # Mortality increase per cm above senescence
    
    # Size-fecundity relationship
    fecundity_exponent: float = 2.5  # Eggs ~ size^exponent
    reference_size: float = 40.0  # Size at which fecundity = reference_eggs
    reference_eggs: float = 2e6  # Eggs at reference size
    
    # Size-dependent SSWD mortality
    sswd_size_effect: float = 0.3  # % mortality increase per 10cm above reference


@dataclass
class SizeStructuredPopulation:
    """
    Size-structured population of Pycnopodia.
    
    Individuals tracked by size rather than age.
    """
    
    config: SizeStructuredConfig
    base_config: Config  # For genetic parameters
    
    # Individual attributes
    sizes: np.ndarray = field(default=None, repr=False)
    genomes: np.ndarray = field(default=None, repr=False)
    sexes: np.ndarray = field(default=None, repr=False)
    alive: np.ndarray = field(default=None, repr=False)
    
    # Per-locus effect sizes
    locus_effects: np.ndarray = field(default=None, repr=False)
    
    year: int = 0
    
    def __post_init__(self):
        if self.sizes is None:
            self._initialize_population()
        if self.locus_effects is None:
            self._initialize_locus_effects()
    
    def _initialize_locus_effects(self):
        """Set per-locus effect sizes."""
        n_loci = self.base_config.n_loci
        effect = self.base_config.per_allele_effect
        self.locus_effects = np.full(n_loci, effect)
    
    def _initialize_population(self):
        """Create initial population with size distribution."""
        n = self.base_config.initial_population
        n_loci = self.base_config.n_loci
        q = self.base_config.initial_resistance_freq
        cfg = self.config
        
        # Size distribution: approximate stable size distribution
        # More small individuals, fewer large
        self.sizes = self._sample_stable_size_distribution(n)
        
        # Genomes in HW equilibrium
        self.genomes = np.zeros((n, n_loci), dtype=np.uint8)
        for locus in range(n_loci):
            allele1 = np.random.random(n) < q
            allele2 = np.random.random(n) < q
            self.genomes[:, locus] = allele1.astype(np.uint8) + allele2.astype(np.uint8)
        
        # 50/50 sex ratio
        self.sexes = np.random.random(n) < 0.5
        
        # All alive
        self.alive = np.ones(n, dtype=bool)
    
    def _sample_stable_size_distribution(self, n: int) -> np.ndarray:
        """
        Sample from approximate stable size distribution.
        
        Uses truncated exponential reflecting high juvenile mortality.
        """
        cfg = self.config
        
        # Exponential decay from small to large sizes
        # Rate parameter based on survival schedule
        rate = 2.0 / (cfg.von_bert_linf - cfg.min_size)
        
        sizes = np.random.exponential(1/rate, n) + cfg.min_size
        sizes = np.clip(sizes, cfg.min_size, cfg.max_size)
        
        return sizes
    
    @property
    def n(self) -> int:
        """Total living individuals."""
        return int(self.alive.sum())
    
    @property
    def n_adults(self) -> int:
        """Number of reproductive adults (size >= maturation_size)."""
        return int((self.alive & (self.sizes >= self.config.maturation_size)).sum())
    
    @property
    def n_females(self) -> int:
        """Number of adult females."""
        return int((self.alive & self.sexes & 
                   (self.sizes >= self.config.maturation_size)).sum())
    
    @property
    def n_males(self) -> int:
        """Number of adult males."""
        return int((self.alive & ~self.sexes & 
                   (self.sizes >= self.config.maturation_size)).sum())
    
    def mean_size(self) -> float:
        """Mean size of living individuals."""
        if self.n == 0:
            return 0.0
        return float(self.sizes[self.alive].mean())
    
    def size_distribution(self, bins: int = 10) -> Tuple[np.ndarray, np.ndarray]:
        """Get size distribution histogram."""
        sizes = self.sizes[self.alive]
        counts, edges = np.histogram(sizes, bins=bins, 
                                     range=(self.config.min_size, self.config.max_size))
        return counts, edges
    
    def resistance_scores(self) -> np.ndarray:
        """Calculate resistance score for each individual."""
        scores = (self.genomes * self.locus_effects).sum(axis=1)
        return np.clip(scores, 0, 1)
    
    def mean_resistance(self) -> float:
        """Mean resistance score among living."""
        if self.n == 0:
            return 0.0
        return float(self.resistance_scores()[self.alive].mean())
    
    def resistance_allele_frequency(self) -> float:
        """Mean resistance allele frequency."""
        if self.n == 0:
            return 0.0
        return float(self.genomes[self.alive].mean() / 2)
    
    def survival_probability(self) -> np.ndarray:
        """
        Calculate size-dependent annual survival probability.
        """
        cfg = self.config
        sizes = self.sizes
        
        # Juvenile survival increases with size
        juvenile_effect = (sizes - cfg.min_size) / (cfg.maturation_size - cfg.min_size)
        juvenile_effect = np.clip(juvenile_effect, 0, 1)
        
        # Adult survival is high but decreases with senescence
        senescence_effect = np.maximum(0, sizes - cfg.senescence_size) * cfg.senescence_rate
        
        # Combine
        survival = np.where(
            sizes < cfg.maturation_size,
            cfg.juvenile_survival_intercept + juvenile_effect * (cfg.adult_survival_max - cfg.juvenile_survival_intercept),
            cfg.adult_survival_max - senescence_effect
        )
        
        return np.clip(survival, 0.1, 0.99)
    
    def fecundity(self) -> np.ndarray:
        """
        Calculate size-dependent fecundity.
        
        Only females above maturation size reproduce.
        """
        cfg = self.config
        
        fec = np.zeros(len(self.sizes))
        
        # Reproductive females
        repro_mask = self.alive & self.sexes & (self.sizes >= cfg.maturation_size)
        
        # Fecundity scales with size^exponent
        size_ratio = self.sizes[repro_mask] / cfg.reference_size
        fec[repro_mask] = cfg.reference_eggs * (size_ratio ** cfg.fecundity_exponent)
        
        return fec
    
    def sswd_susceptibility(self) -> np.ndarray:
        """
        Calculate size-dependent SSWD susceptibility.
        
        Larger individuals are more susceptible.
        """
        cfg = self.config
        
        # Size effect: larger = more susceptible
        size_effect = 1 + cfg.sswd_size_effect * (self.sizes - cfg.reference_size) / 10
        size_effect = np.clip(size_effect, 0.5, 2.0)
        
        # Resistance effect
        resistance = self.resistance_scores()
        
        # Combined susceptibility
        susceptibility = size_effect * (1 - resistance)
        
        return np.clip(susceptibility, 0, 2)
    
    def grow(self, rng: np.random.Generator):
        """
        Apply one year of growth.
        
        Uses von Bertalanffy growth with individual variation.
        """
        cfg = self.config
        
        alive_mask = self.alive
        current_sizes = self.sizes[alive_mask]
        
        # von Bertalanffy: dL/dt = k * (Linf - L)
        expected_growth = cfg.von_bert_k * (cfg.von_bert_linf - current_sizes)
        expected_growth = np.maximum(0, expected_growth)
        
        # Individual variation
        growth_sd = expected_growth * cfg.growth_cv
        actual_growth = rng.normal(expected_growth, growth_sd)
        actual_growth = np.maximum(0, actual_growth)
        
        # Update sizes
        self.sizes[alive_mask] = np.clip(
            current_sizes + actual_growth,
            cfg.min_size,
            cfg.max_size
        )
    
    def apply_natural_mortality(self, rng: np.random.Generator):
        """Apply size-dependent natural mortality."""
        survival = self.survival_probability()
        mortality = 1 - survival
        
        deaths = rng.random(len(self.alive)) < mortality
        self.alive &= ~deaths
    
    def apply_sswd_mortality(
        self, 
        prevalence: float, 
        base_mortality: float,
        rng: np.random.Generator
    ) -> int:
        """
        Apply SSWD mortality with size dependence.
        
        Returns number of deaths.
        """
        if prevalence == 0:
            return 0
        
        susceptibility = self.sswd_susceptibility()
        mortality = prevalence * base_mortality * susceptibility
        mortality = np.clip(mortality, 0, 0.99)
        
        # Apply only to living
        mortality[~self.alive] = 0
        
        n_before = self.n
        deaths = rng.random(len(self.alive)) < mortality
        self.alive &= ~deaths
        
        return n_before - self.n
    
    def compact(self):
        """Remove dead individuals."""
        alive_idx = self.alive
        self.sizes = self.sizes[alive_idx]
        self.genomes = self.genomes[alive_idx]
        self.sexes = self.sexes[alive_idx]
        self.alive = np.ones(self.n, dtype=bool)
    
    def add_individuals(
        self,
        sizes: np.ndarray,
        genomes: np.ndarray,
        sexes: np.ndarray
    ):
        """Add new individuals."""
        self.sizes = np.concatenate([self.sizes, sizes])
        self.genomes = np.vstack([self.genomes, genomes])
        self.sexes = np.concatenate([self.sexes, sexes])
        self.alive = np.concatenate([self.alive, np.ones(len(sizes), dtype=bool)])
    
    def get_summary(self) -> Dict:
        """Get summary statistics."""
        return {
            "n_total": self.n,
            "n_adults": self.n_adults,
            "n_females": self.n_females,
            "n_males": self.n_males,
            "mean_size": self.mean_size(),
            "mean_resistance": self.mean_resistance(),
            "resistance_allele_freq": self.resistance_allele_frequency(),
        }
