"""
Reproduction module implementing:
- Sweepstakes reproductive success (SRS)
- Broadcast spawner Allee effects
- Mendelian inheritance
"""

import numpy as np
from typing import Tuple
from .config import Config
from .population import Population


def sample_breeding_fraction(config: Config, rng: np.random.Generator) -> float:
    """
    Sample the fraction of adults that successfully breed this year.
    
    Uses a Beta distribution to model high variance in reproductive success.
    Parameterized by mean and concentration (shape) parameter.
    
    NOTE: For a peaked (unimodal) distribution, need α > 1 and β > 1.
    With mean=0.08, this requires shape > 12.5. Lower shape values
    produce U-shaped distributions (bimodal at 0 and 1).
    
    Parameters
    ----------
    config : Config
        Model configuration
    rng : Generator
        Random number generator
        
    Returns
    -------
    breeding_fraction : float
        Fraction of adults that breed (0 to 1)
    """
    mean = config.srs_mean
    shape = config.srs_shape
    
    # Beta distribution parameterized by mean and concentration
    # Higher shape = less variance, more peaked around mean
    # For mean=0.08, shape=15 gives CV≈0.8, shape=5 gives CV≈1.4
    alpha = shape * mean
    beta = shape * (1 - mean)
    
    # Ensure valid parameters (α, β > 0)
    alpha = max(0.01, alpha)
    beta = max(0.01, beta)
    
    fraction = rng.beta(alpha, beta)
    
    # Clip to configured range
    fraction = np.clip(fraction, config.srs_min, config.srs_max)
    
    return float(fraction)


def select_breeders(
    population: Population,
    breeding_fraction: float,
    rng: np.random.Generator
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Select which adults breed this year (sweepstakes).
    
    Parameters
    ----------
    population : Population
        Current population
    breeding_fraction : float
        Fraction of adults that breed
    rng : Generator
        Random number generator
        
    Returns
    -------
    female_indices : ndarray
        Indices of breeding females
    male_indices : ndarray
        Indices of breeding males
    """
    config = population.config
    
    # Get adult indices
    adult_mask = (population.ages >= config.maturation_age) & population.alive
    adult_indices = np.where(adult_mask)[0]
    
    if len(adult_indices) == 0:
        return np.array([]), np.array([])
    
    # Determine which adults breed
    n_breeders = max(1, int(len(adult_indices) * breeding_fraction))
    breeder_indices = rng.choice(adult_indices, size=n_breeders, replace=False)
    
    # Split by sex
    female_mask = population.sexes[breeder_indices]
    female_indices = breeder_indices[female_mask]
    male_indices = breeder_indices[~female_mask]
    
    return female_indices, male_indices


def fertilization_success_saturating(
    n_adults: int,
    half_saturation: float
) -> float:
    """
    Saturating Allee effect on fertilization success.
    
    f = N^2 / (N^2 + h^2)
    
    At low density, gamete dilution reduces fertilization.
    """
    if n_adults == 0:
        return 0.0
    return n_adults**2 / (n_adults**2 + half_saturation**2)


def fertilization_success_levitan(
    n_females: int,
    n_males: int,
    config: Config
) -> float:
    """
    Mechanistic Levitan-style fertilization model.
    
    Based on sperm-egg collision kinetics in turbulent flow.
    See Levitan 1991, Yund 2000.
    
    Parameters
    ----------
    n_females, n_males : int
        Number of spawning adults
    config : Config
        Model configuration with sperm parameters
        
    Returns
    -------
    fertilization_rate : float
        Fraction of eggs fertilized (0 to 1)
    """
    if n_females == 0 or n_males == 0:
        return 0.0
    
    # Sperm concentration (sperm per m^3)
    # Assuming spawning in ~1000 m^3 water volume
    volume = 1000  # m^3
    total_sperm = n_males * config.sperm_release_rate
    sperm_conc = total_sperm / volume
    
    # Egg radius and contact rate
    r = config.egg_radius_m
    u = config.current_speed_m_s
    
    # Collision rate constant (simplified)
    beta = 4 * np.pi * r**2 * u
    
    # Fertilization probability per egg (Poisson process)
    # P(fertilized) = 1 - exp(-beta * S * t)
    contact_time = 3600  # seconds (1 hour spawning window)
    p_fert = 1 - np.exp(-beta * sperm_conc * contact_time)
    
    return float(np.clip(p_fert, 0, 1))


def calculate_recruitment(
    population: Population,
    female_indices: np.ndarray,
    male_indices: np.ndarray,
    config: Config,
    rng: np.random.Generator
) -> int:
    """
    Calculate number of recruits (settled juveniles).
    
    Incorporates:
    - Fecundity (eggs per female)
    - Larval survival
    - Allee effect on fertilization
    - Density-dependent ceiling
    - Fecundity cost of resistance
    """
    n_females = len(female_indices)
    n_males = len(male_indices)
    
    if n_females == 0 or n_males == 0:
        return 0
    
    # Fertilization success (Allee effect)
    if config.allee_model == "levitan":
        fert_success = fertilization_success_levitan(n_females, n_males, config)
    else:
        n_adults = population.n_adults
        fert_success = fertilization_success_saturating(n_adults, config.allee_half_saturation)
    
    # Calculate fecundity with resistance cost
    # Cost is per resistance ALLELE (not per effect unit)
    # Count total resistance alleles per female
    n_resistance_alleles = population.genomes[female_indices].sum(axis=1)
    fecundity_costs = n_resistance_alleles * config.fecundity_cost_per_allele
    effective_fecundity = config.fecundity * np.maximum(0, 1 - fecundity_costs)
    
    # Total potential settlers
    eggs_per_female = effective_fecundity * fert_success
    larvae_per_female = eggs_per_female * config.larval_survival
    total_potential = larvae_per_female.sum()
    
    # Density-dependent ceiling
    available_space = max(0, config.carrying_capacity - population.n)
    
    # Stochastic recruitment (Poisson)
    expected_recruits = min(total_potential, available_space)
    if expected_recruits <= 0:
        return 0
    
    n_recruits = rng.poisson(expected_recruits)
    return int(min(n_recruits, available_space))


def generate_offspring(
    population: Population,
    female_indices: np.ndarray,
    male_indices: np.ndarray,
    n_offspring: int,
    rng: np.random.Generator
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Generate offspring genomes via Mendelian inheritance.
    
    Parameters
    ----------
    population : Population
        Parent population
    female_indices, male_indices : ndarray
        Indices of breeding adults
    n_offspring : int
        Number of offspring to generate
    rng : Generator
        Random number generator
        
    Returns
    -------
    genomes : ndarray, shape (n_offspring, n_loci)
        Offspring genotypes
    ages : ndarray
        All zeros (newborns)
    sexes : ndarray
        Random sex assignment
    """
    if n_offspring == 0 or len(female_indices) == 0 or len(male_indices) == 0:
        return (
            np.zeros((0, population.config.n_loci), dtype=np.uint8),
            np.zeros(0, dtype=np.int8),
            np.zeros(0, dtype=bool)
        )
    
    n_loci = population.config.n_loci
    
    # Sample parents for each offspring (with replacement)
    mother_idx = rng.choice(female_indices, size=n_offspring)
    father_idx = rng.choice(male_indices, size=n_offspring)
    
    # Get parent genomes
    mother_genomes = population.genomes[mother_idx]
    father_genomes = population.genomes[father_idx]
    
    # Mendelian segregation - fully vectorized
    # For each parent genotype:
    #   0 -> always transmit 0
    #   1 -> 50% chance of 0 or 1  
    #   2 -> always transmit 1
    
    # Generate random draws for heterozygotes only (more efficient)
    rand_mat = rng.random((n_offspring, n_loci))
    rand_pat = rng.random((n_offspring, n_loci))
    
    # Maternal alleles: 0 if geno=0, 1 if geno=2, random if geno=1
    mat_allele = np.where(
        mother_genomes == 0, 0,
        np.where(mother_genomes == 2, 1, (rand_mat < 0.5).astype(np.uint8))
    )
    
    # Paternal alleles
    pat_allele = np.where(
        father_genomes == 0, 0,
        np.where(father_genomes == 2, 1, (rand_pat < 0.5).astype(np.uint8))
    )
    
    offspring_genomes = mat_allele + pat_allele
    
    # All newborns are age 0
    ages = np.zeros(n_offspring, dtype=np.int8)
    
    # 50/50 sex ratio
    sexes = rng.random(n_offspring) < 0.5
    
    return offspring_genomes, ages, sexes


def reproduce(
    population: Population,
    rng: np.random.Generator
) -> Tuple[int, int, int, float]:
    """
    Execute one year of reproduction.
    
    Returns
    -------
    n_recruits : int
        Number of new recruits added
    n_female_breeders : int
        Number of females that bred
    n_male_breeders : int
        Number of males that bred
    breeding_fraction : float
        Fraction of adults that bred
    """
    config = population.config
    
    # Sample SRS breeding fraction
    breeding_fraction = sample_breeding_fraction(config, rng)
    
    # Select breeders
    female_idx, male_idx = select_breeders(population, breeding_fraction, rng)
    
    if len(female_idx) == 0 or len(male_idx) == 0:
        return 0, 0, 0, breeding_fraction
    
    # Calculate recruitment
    n_recruits = calculate_recruitment(population, female_idx, male_idx, config, rng)
    
    if n_recruits > 0:
        # Generate offspring
        genomes, ages, sexes = generate_offspring(
            population, female_idx, male_idx, n_recruits, rng
        )
        # Add to population
        population.add_individuals(genomes, ages, sexes)
    
    return n_recruits, len(female_idx), len(male_idx), breeding_fraction
