"""
SSWD (Sea Star Wasting Disease) dynamics.

Implements:
- Time-varying prevalence (epidemic -> endemic)
- Age-dependent susceptibility
- Resistance-modified mortality
- Temperature-dependent dynamics (optional)
- Density-dependent transmission (optional)
"""

import numpy as np
from typing import Tuple
from .config import Config
from .population import Population


def calculate_prevalence(
    year: int,
    config: Config,
    temperature: float = None
) -> float:
    """
    Calculate SSWD prevalence at given year.
    
    Prevalence follows exponential decay from peak to endemic:
    P(t) = P_endemic + (P_peak - P_endemic) * exp(-δ * (t - t_onset))
    
    Parameters
    ----------
    year : int
        Current simulation year
    config : Config
        Model configuration
    temperature : float, optional
        Local temperature for temperature-dependent model
        
    Returns
    -------
    prevalence : float
        Disease prevalence (0 to 1)
    """
    onset = config.sswd_onset_year
    
    if year < onset:
        return 0.0
    
    peak = config.sswd_peak_prevalence
    endemic = config.sswd_endemic_prevalence
    decay = config.sswd_decay_rate
    
    t = year - onset
    prevalence = endemic + (peak - endemic) * np.exp(-decay * t)
    
    # Temperature modification (optional)
    if config.temperature_dependent_sswd and temperature is not None:
        temp_effect = config.sswd_temperature_sensitivity * (temperature - config.base_temperature_c)
        # Higher temperature -> higher prevalence
        prevalence = prevalence * (1 + temp_effect)
    
    return float(np.clip(prevalence, 0, 1))


def calculate_sswd_mortality(
    population: Population,
    prevalence: float,
    rng: np.random.Generator
) -> np.ndarray:
    """
    Calculate per-individual SSWD mortality probability.
    
    m_sswd(i) = m_base * P(t) * max(0, 1 - resistance_score(i))
    
    Parameters
    ----------
    population : Population
        Current population
    prevalence : float
        Current disease prevalence
    rng : Generator
        Random number generator
        
    Returns
    -------
    mortality_probs : ndarray
        Per-individual mortality probability
    """
    config = population.config
    n = len(population.alive)
    
    if prevalence == 0:
        return np.zeros(n)
    
    # Base mortality (age-dependent)
    base_mortality = np.zeros(n)
    adult_mask = population.ages >= config.maturation_age
    base_mortality[adult_mask] = config.sswd_adult_mortality
    base_mortality[~adult_mask] = config.sswd_juvenile_mortality
    
    # Resistance modification
    resistance = population.resistance_scores()
    susceptibility = np.maximum(0, 1 - resistance)
    
    # Final mortality = base * prevalence * susceptibility
    mortality_probs = base_mortality * prevalence * susceptibility
    
    # Only apply to living individuals
    mortality_probs[~population.alive] = 0
    
    return mortality_probs


def density_dependent_transmission(
    population: Population,
    base_prevalence: float,
    config: Config
) -> float:
    """
    Modify prevalence based on population density.
    
    At low density, transmission is reduced (fewer contacts).
    At high density, transmission is enhanced.
    
    Uses simple frequency-dependent transmission:
    P_effective = P_base * (N / K)
    """
    if not config.density_dependent_transmission:
        return base_prevalence
    
    density_ratio = population.n / config.carrying_capacity
    
    # Transmission scales with density
    effective_prevalence = base_prevalence * density_ratio
    
    return float(np.clip(effective_prevalence, 0, 1))


def apply_sswd(
    population: Population,
    year: int,
    rng: np.random.Generator,
    temperature: float = None
) -> Tuple[int, float]:
    """
    Apply SSWD mortality to population.
    
    Parameters
    ----------
    population : Population
        Current population
    year : int
        Current simulation year
    rng : Generator
        Random number generator
    temperature : float, optional
        Local temperature (for temperature-dependent model)
        
    Returns
    -------
    n_deaths : int
        Number of SSWD deaths
    prevalence : float
        Disease prevalence this year
    """
    config = population.config
    
    # Calculate prevalence
    prevalence = calculate_prevalence(year, config, temperature)
    
    # Modify by density if enabled
    prevalence = density_dependent_transmission(population, prevalence, config)
    
    if prevalence == 0:
        return 0, 0.0
    
    # Calculate mortality
    mortality_probs = calculate_sswd_mortality(population, prevalence, rng)
    
    # Count deaths before applying
    n_alive_before = population.n
    
    # Apply mortality
    population.apply_mortality(mortality_probs)
    
    n_deaths = n_alive_before - population.n
    
    return n_deaths, prevalence
