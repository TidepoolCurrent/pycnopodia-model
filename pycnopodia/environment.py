"""
Environmental stochasticity module.

Implements:
- Annual recruitment variability (good/bad years)
- Marine heatwave events
- Temperature time series from climate projections
- Oceanographic connectivity for larval dispersal
"""

import numpy as np
from typing import Tuple, Optional, List
from dataclasses import dataclass, field


@dataclass
class EnvironmentConfig:
    """Environmental variability parameters."""
    
    # Recruitment stochasticity
    recruitment_cv: float = 0.8  # Coefficient of variation (empirical: 0.5-1.5)
    recruitment_autocorr: float = 0.3  # Year-to-year correlation
    
    # Marine heatwaves
    heatwave_probability: float = 0.05  # Annual probability of MHW
    heatwave_intensity: float = 2.0  # Temperature anomaly (°C)
    heatwave_duration: int = 1  # Years (can be multi-year)
    heatwave_sswd_multiplier: float = 1.5  # SSWD mortality increase during MHW
    
    # Climate trend
    warming_rate: float = 0.03  # °C per year (roughly 2.4°C over 80 years)
    
    # Larval connectivity (for metapopulation)
    use_connectivity_matrix: bool = False
    connectivity_matrix: np.ndarray = None  # If provided, overrides simple migration


class EnvironmentState:
    """
    Tracks environmental state across simulation.
    
    Generates correlated recruitment multipliers and temperature anomalies.
    """
    
    def __init__(self, config: EnvironmentConfig, n_years: int, seed: int = None):
        self.config = config
        self.n_years = n_years
        self.rng = np.random.default_rng(seed)
        
        # Pre-generate environmental time series
        self._generate_recruitment_series()
        self._generate_temperature_series()
        self._generate_heatwaves()
    
    def _generate_recruitment_series(self):
        """
        Generate autocorrelated recruitment multipliers.
        
        Uses AR(1) process on log scale to create lognormal variation.
        """
        cv = self.config.recruitment_cv
        rho = self.config.recruitment_autocorr
        
        # Log-scale standard deviation
        sigma = np.sqrt(np.log(1 + cv**2))
        
        # AR(1) process
        innovations = self.rng.normal(0, sigma * np.sqrt(1 - rho**2), self.n_years)
        log_mult = np.zeros(self.n_years)
        log_mult[0] = innovations[0]
        
        for t in range(1, self.n_years):
            log_mult[t] = rho * log_mult[t-1] + innovations[t]
        
        # Convert to multipliers (median = 1)
        self.recruitment_multipliers = np.exp(log_mult - sigma**2/2)
    
    def _generate_temperature_series(self):
        """
        Generate temperature time series with trend and variability.
        """
        # Linear warming trend
        trend = np.arange(self.n_years) * self.config.warming_rate
        
        # Interannual variability (AR(1))
        variability = np.zeros(self.n_years)
        for t in range(1, self.n_years):
            variability[t] = 0.5 * variability[t-1] + self.rng.normal(0, 0.3)
        
        self.temperature_anomaly = trend + variability
    
    def _generate_heatwaves(self):
        """
        Generate discrete marine heatwave events.
        """
        self.heatwave_years = set()
        
        for t in range(self.n_years):
            if self.rng.random() < self.config.heatwave_probability:
                # MHW occurs
                for d in range(self.config.heatwave_duration):
                    if t + d < self.n_years:
                        self.heatwave_years.add(t + d)
    
    def get_recruitment_multiplier(self, year: int) -> float:
        """Get recruitment multiplier for given year."""
        if 0 <= year < self.n_years:
            return float(self.recruitment_multipliers[year])
        return 1.0
    
    def get_temperature_anomaly(self, year: int) -> float:
        """Get temperature anomaly (°C above baseline) for given year."""
        if 0 <= year < self.n_years:
            base = self.temperature_anomaly[year]
            # Add heatwave intensity if active
            if year in self.heatwave_years:
                base += self.config.heatwave_intensity
            return float(base)
        return 0.0
    
    def is_heatwave(self, year: int) -> bool:
        """Check if year is during a marine heatwave."""
        return year in self.heatwave_years
    
    def get_sswd_multiplier(self, year: int) -> float:
        """Get SSWD mortality multiplier based on temperature."""
        if year in self.heatwave_years:
            return self.config.heatwave_sswd_multiplier
        # Gradual increase with warming
        temp_effect = 1 + 0.1 * max(0, self.get_temperature_anomaly(year))
        return temp_effect


def generate_connectivity_matrix(
    n_subpops: int,
    dispersal_kernel: str = "stepping_stone",
    mean_dispersal: float = 0.1
) -> np.ndarray:
    """
    Generate larval connectivity matrix.
    
    Parameters
    ----------
    n_subpops : int
        Number of subpopulations
    dispersal_kernel : str
        "stepping_stone" - adjacent only
        "island" - equal to all
        "distance_decay" - decay with distance
    mean_dispersal : float
        Mean fraction of larvae dispersing
        
    Returns
    -------
    connectivity : ndarray, shape (n_subpops, n_subpops)
        connectivity[i,j] = fraction of larvae from i settling in j
    """
    C = np.zeros((n_subpops, n_subpops))
    
    if dispersal_kernel == "stepping_stone":
        # Only adjacent subpopulations exchange larvae
        for i in range(n_subpops):
            C[i, i] = 1 - mean_dispersal  # Self-recruitment
            if i > 0:
                C[i, i-1] = mean_dispersal / 2
            if i < n_subpops - 1:
                C[i, i+1] = mean_dispersal / 2
            # Normalize row
            C[i] /= C[i].sum()
    
    elif dispersal_kernel == "island":
        # Equal dispersal to all
        C = np.full((n_subpops, n_subpops), mean_dispersal / (n_subpops - 1))
        np.fill_diagonal(C, 1 - mean_dispersal)
    
    elif dispersal_kernel == "distance_decay":
        # Exponential decay with distance
        for i in range(n_subpops):
            for j in range(n_subpops):
                dist = abs(i - j)
                C[i, j] = np.exp(-dist / 2)
            C[i] /= C[i].sum()
    
    return C


# Default PNW connectivity based on simplified oceanography
# Columns: AK, BC, WA, OR, CA (source → destination)
PNW_CONNECTIVITY = np.array([
    [0.85, 0.10, 0.04, 0.01, 0.00],  # From Alaska
    [0.05, 0.80, 0.10, 0.04, 0.01],  # From BC
    [0.01, 0.08, 0.75, 0.12, 0.04],  # From WA
    [0.00, 0.02, 0.08, 0.70, 0.20],  # From OR
    [0.00, 0.00, 0.02, 0.18, 0.80],  # From CA
])
