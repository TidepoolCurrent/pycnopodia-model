"""
Captive breeding and outplanting interventions.

Implements:
- Broodstock maintenance
- Captive breeding with Mendelian inheritance
- Outplanting of captive-bred juveniles
- Optional inbreeding tracking
"""

import numpy as np
from typing import Tuple
from dataclasses import dataclass, field
from .config import Config


@dataclass
class Broodstock:
    """
    Captive broodstock population for breeding.
    
    Simplified representation focusing on genetics.
    """
    config: Config
    genomes: np.ndarray = field(default=None, repr=False)
    sexes: np.ndarray = field(default=None, repr=False)
    
    # Inbreeding tracking
    inbreeding_coefficients: np.ndarray = field(default=None, repr=False)
    generation: int = 0
    
    def __post_init__(self):
        if self.genomes is None:
            self._initialize_broodstock()
    
    def _initialize_broodstock(self):
        """Create initial broodstock population."""
        n = self.config.broodstock_n_parents
        n_loci = self.config.n_loci
        q = self.config.broodstock_resistance_freq
        
        # Initialize genomes in HW equilibrium
        self.genomes = np.zeros((n, n_loci), dtype=np.uint8)
        for locus in range(n_loci):
            allele1 = np.random.random(n) < q
            allele2 = np.random.random(n) < q
            self.genomes[:, locus] = allele1.astype(np.uint8) + allele2.astype(np.uint8)
        
        # 50/50 sex ratio
        self.sexes = np.random.random(n) < 0.5
        
        # Initial inbreeding = 0
        self.inbreeding_coefficients = np.zeros(n)
    
    @property
    def n_females(self) -> int:
        return int(self.sexes.sum())
    
    @property
    def n_males(self) -> int:
        return int((~self.sexes).sum())
    
    def resistance_allele_frequency(self) -> float:
        """Mean resistance allele frequency in broodstock."""
        return float(self.genomes.mean() / 2)
    
    def mean_inbreeding(self) -> float:
        """Mean inbreeding coefficient."""
        if self.inbreeding_coefficients is None:
            return 0.0
        return float(self.inbreeding_coefficients.mean())
    
    def breed_offspring(
        self,
        n_offspring: int,
        rng: np.random.Generator
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Breed offspring from broodstock.
        
        Parameters
        ----------
        n_offspring : int
            Number of offspring to produce
        rng : Generator
            Random number generator
            
        Returns
        -------
        genomes : ndarray
            Offspring genotypes
        ages : ndarray
            All zeros (juveniles)
        sexes : ndarray
            Random sex assignment
        """
        if n_offspring == 0:
            n_loci = self.config.n_loci
            return (
                np.zeros((0, n_loci), dtype=np.uint8),
                np.zeros(0, dtype=np.int8),
                np.zeros(0, dtype=bool)
            )
        
        female_idx = np.where(self.sexes)[0]
        male_idx = np.where(~self.sexes)[0]
        
        if len(female_idx) == 0 or len(male_idx) == 0:
            return (
                np.zeros((0, self.config.n_loci), dtype=np.uint8),
                np.zeros(0, dtype=np.int8),
                np.zeros(0, dtype=bool)
            )
        
        n_loci = self.config.n_loci
        
        # Sample parents
        mother_idx = rng.choice(female_idx, size=n_offspring)
        father_idx = rng.choice(male_idx, size=n_offspring)
        
        mother_genomes = self.genomes[mother_idx]
        father_genomes = self.genomes[father_idx]
        
        # Mendelian inheritance
        offspring_genomes = np.zeros((n_offspring, n_loci), dtype=np.uint8)
        
        for locus in range(n_loci):
            mat_geno = mother_genomes[:, locus]
            mat_allele = np.where(
                mat_geno == 0, 0,
                np.where(mat_geno == 2, 1, rng.integers(0, 2, n_offspring))
            )
            
            pat_geno = father_genomes[:, locus]
            pat_allele = np.where(
                pat_geno == 0, 0,
                np.where(pat_geno == 2, 1, rng.integers(0, 2, n_offspring))
            )
            
            offspring_genomes[:, locus] = mat_allele + pat_allele
        
        ages = np.zeros(n_offspring, dtype=np.int8)
        sexes = rng.random(n_offspring) < 0.5
        
        # Update inbreeding if tracking
        if self.config.track_broodstock_inbreeding:
            self._update_inbreeding()
        
        return offspring_genomes, ages, sexes
    
    def _update_inbreeding(self):
        """
        Update inbreeding coefficients after a generation.
        
        Simplified model: F increases by 1/(2*Ne) per generation
        """
        n = len(self.genomes)
        ne = 4 * self.n_females * self.n_males / (self.n_females + self.n_males + 0.001)
        
        delta_f = 1 / (2 * ne + 0.001)
        self.inbreeding_coefficients = self.inbreeding_coefficients * (1 - delta_f) + delta_f
        self.generation += 1
    
    def replace_individuals(
        self,
        n_replace: int,
        source_genomes: np.ndarray,
        source_sexes: np.ndarray,
        rng: np.random.Generator
    ):
        """
        Replace some broodstock individuals with wild-caught or new individuals.
        
        Parameters
        ----------
        n_replace : int
            Number to replace
        source_genomes : ndarray
            Genomes of potential replacements
        source_sexes : ndarray
            Sexes of potential replacements
        rng : Generator
            Random number generator
        """
        if n_replace == 0 or len(source_genomes) == 0:
            return
        
        n_replace = min(n_replace, len(self.genomes), len(source_genomes))
        
        # Select which to replace (random)
        replace_idx = rng.choice(len(self.genomes), size=n_replace, replace=False)
        
        # Select replacements from source
        source_idx = rng.choice(len(source_genomes), size=n_replace, replace=False)
        
        self.genomes[replace_idx] = source_genomes[source_idx]
        self.sexes[replace_idx] = source_sexes[source_idx]
        
        # Reset inbreeding for replaced individuals
        if self.inbreeding_coefficients is not None:
            self.inbreeding_coefficients[replace_idx] = 0


def should_outplant(year: int, config: Config) -> bool:
    """Check if outplanting should occur this year."""
    return config.outplanting_start_year <= year <= config.outplanting_end_year


def outplant(
    population,  # Population - avoid circular import
    broodstock: Broodstock,
    year: int,
    rng: np.random.Generator
) -> int:
    """
    Outplant captive-bred juveniles into wild population.
    
    Parameters
    ----------
    population : Population
        Wild population to outplant into
    broodstock : Broodstock
        Captive breeding population
    year : int
        Current year
    rng : Generator
        Random number generator
        
    Returns
    -------
    n_outplanted : int
        Number of individuals outplanted
    """
    config = population.config
    
    if not should_outplant(year, config):
        return 0
    
    n_to_outplant = config.outplanting_n_per_year
    
    # Breed offspring
    genomes, ages, sexes = broodstock.breed_offspring(n_to_outplant, rng)
    
    if len(genomes) == 0:
        return 0
    
    # Set age to 1 (post-settlement juveniles)
    ages = np.ones(len(ages), dtype=np.int8)
    
    # Add to population
    population.add_individuals(genomes, ages, sexes)
    
    # Optional: refresh broodstock from wild population
    if config.track_broodstock_inbreeding and population.n_adults > 10:
        n_refresh = max(1, int(config.broodstock_n_parents * config.broodstock_replacement_rate))
        adult_mask = (population.ages >= config.maturation_age) & population.alive
        adult_idx = np.where(adult_mask)[0]
        if len(adult_idx) > 0:
            broodstock.replace_individuals(
                n_refresh,
                population.genomes[adult_idx],
                population.sexes[adult_idx],
                rng
            )
    
    return len(genomes)
