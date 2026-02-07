"""
Model configuration and default parameters.

All parameters are documented with sources from the literature.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional
import numpy as np


@dataclass
class Config:
    """Model configuration parameters."""
    
    # --- Demography ---
    carrying_capacity: int = 1500  # Representative local population
    initial_population: int = 1200  # Starting adult population
    
    # Age-specific survival (annual)
    survival_juvenile_y0: float = 0.15  # Post-settlement, extremely high mortality
    survival_juvenile_y1_4: float = 0.70  # Pre-reproductive juveniles
    survival_adult: float = 0.88  # Reproductive adults (5-18 years)
    survival_senescent: float = 0.26  # Sharp decline 19+ years
    
    maturation_age: int = 5  # Years to reproductive maturity (Hodin et al. 2021)
    max_age: int = 25  # Maximum lifespan
    
    # --- Reproduction ---
    fecundity: float = 2e6  # Eggs per female per spawning (Menge 1974; NOAA 2022)
    larval_survival: float = 5e-5  # Planktonic larval survival probability
    
    # Sweepstakes reproductive success (Hedgecock & Pudovkin 2011)
    srs_mean: float = 0.08  # Mean fraction of adults breeding per year
    srs_shape: float = 5.0  # Beta concentration param (higher = less variance)
    # NOTE: shape=5 with mean=0.08 gives α=0.4, β=4.6, CV≈1.3 (high variance)
    # For peaked distribution around mean, need shape > 12 (gives α > 1)
    
    # IMPROVEMENT: Allow exploration of more extreme SRS
    srs_min: float = 0.02  # Minimum breeding fraction
    srs_max: float = 0.30  # Maximum breeding fraction
    
    # Allee effect on fertilization (broadcast spawner gamete dilution)
    allee_half_saturation: float = 30  # Adults needed for 50% fertilization
    # IMPROVEMENT: Mechanistic Allee model option
    allee_model: str = "saturating"  # "saturating" or "levitan" (mechanistic)
    sperm_release_rate: float = 1e9  # Sperm per male per spawn (for Levitan model)
    egg_radius_m: float = 150e-6  # Egg radius in meters
    current_speed_m_s: float = 0.1  # Water current speed
    
    # --- Genetics ---
    n_loci: int = 20  # Resistance loci (subset of 51 from Schiebelhut et al. 2024)
    initial_resistance_freq: float = 0.08  # Starting resistance allele frequency
    per_allele_effect: float = 0.035  # 3.5% mortality reduction per resistance allele
    fecundity_cost_per_allele: float = 0.002  # 0.2% fecundity cost per resistance allele
    
    # IMPROVEMENT: Variable effect sizes across loci
    variable_effect_sizes: bool = False  # If True, sample effect sizes from distribution
    effect_size_cv: float = 0.5  # Coefficient of variation for effect sizes
    
    # --- Disease (SSWD) ---
    sswd_onset_year: int = 10  # Year of epidemic onset
    sswd_peak_prevalence: float = 0.90  # Peak disease prevalence (Hewson et al. 2014)
    sswd_endemic_prevalence: float = 0.20  # Long-term endemic level
    sswd_decay_rate: float = 0.10  # Exponential decay from peak to endemic
    sswd_adult_mortality: float = 0.60  # Base SSWD mortality for adults
    sswd_juvenile_mortality: float = 0.40  # Base SSWD mortality for juveniles
    
    # IMPROVEMENT: Temperature-dependent disease dynamics
    temperature_dependent_sswd: bool = False
    base_temperature_c: float = 12.0  # Reference temperature
    sswd_temperature_sensitivity: float = 0.05  # Mortality increase per degree C
    
    # IMPROVEMENT: Density-dependent transmission
    density_dependent_transmission: bool = False
    transmission_rate: float = 0.001  # Per-contact transmission probability
    
    # --- Intervention (Outplanting) ---
    outplanting_start_year: int = 15  # When outplanting begins
    outplanting_end_year: int = 50  # When outplanting ends
    outplanting_n_per_year: int = 200  # Juveniles outplanted per year
    broodstock_n_parents: int = 30  # Number of breeding adults in captive program
    broodstock_resistance_freq: float = 0.15  # Resistance allele freq in broodstock
    
    # IMPROVEMENT: Broodstock inbreeding accumulation
    track_broodstock_inbreeding: bool = False
    broodstock_replacement_rate: float = 0.1  # Fraction of broodstock replaced/year
    
    # --- Metapopulation ---
    n_subpopulations: int = 1  # 1 = single population, 5 = full metapopulation
    migration_rate: float = 0.005  # Adult migration between adjacent subpops
    
    # Subpopulation-specific parameters (when n_subpopulations > 1)
    subpop_names: List[str] = field(default_factory=lambda: [
        "Alaska", "British Columbia", "Washington", "Oregon", "California"
    ])
    subpop_carrying_capacities: List[int] = field(default_factory=lambda: [
        1800, 1500, 1200, 1000, 700  # Reflects habitat area
    ])
    subpop_sswd_onset_years: List[int] = field(default_factory=lambda: [
        14, 12, 11, 10, 10  # South-first epidemic spread
    ])
    subpop_sswd_peak_prevalence: List[float] = field(default_factory=lambda: [
        0.70, 0.80, 0.85, 0.90, 0.95  # Higher in warmer southern waters
    ])
    subpop_sswd_endemic_prevalence: List[float] = field(default_factory=lambda: [
        0.10, 0.15, 0.20, 0.25, 0.30  # Temperature gradient
    ])
    # IMPROVEMENT: Subpop-specific temperatures for thermal refugia modeling
    subpop_temperatures_c: List[float] = field(default_factory=lambda: [
        8.0, 10.0, 11.5, 13.0, 15.0
    ])
    
    # --- Simulation ---
    n_years: int = 80  # Simulation duration
    n_replicates: int = 200  # Number of stochastic replicates
    random_seed: Optional[int] = None  # For reproducibility
    
    def validate(self):
        """Validate parameter combinations."""
        assert 0 < self.srs_mean <= 1, "SRS mean must be in (0, 1]"
        assert self.n_loci > 0, "Must have at least 1 locus"
        assert self.carrying_capacity > 0, "K must be positive"
        if self.n_subpopulations > 1:
            assert len(self.subpop_names) >= self.n_subpopulations
        return True


# Default configuration
DEFAULT_CONFIG = Config()

# Preset configurations for different scenarios
PRESETS = {
    "baseline": Config(),
    
    "no_disease": Config(
        sswd_peak_prevalence=0.0,
        sswd_endemic_prevalence=0.0,
    ),
    
    "extreme_srs": Config(
        srs_mean=0.02,  # Only 2% breeding
        srs_shape=1.2,  # Even higher variance
    ),
    
    "wright_fisher": Config(
        srs_mean=1.0,  # All adults breed
        srs_shape=100,  # Very low variance
    ),
    
    "high_outplanting": Config(
        outplanting_n_per_year=500,
        broodstock_n_parents=50,
    ),
    
    "limited_broodstock": Config(
        outplanting_n_per_year=200,
        broodstock_n_parents=6,
        broodstock_resistance_freq=0.08,
    ),
    
    "metapopulation": Config(
        n_subpopulations=5,
    ),
    
    # IMPROVEMENT: Climate warming scenario
    "climate_warming": Config(
        temperature_dependent_sswd=True,
        # Assume 2C warming over simulation
    ),
}
