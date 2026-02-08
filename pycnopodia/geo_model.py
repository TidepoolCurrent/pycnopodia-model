"""
Geography-based Pycnopodia simulation using real site coordinates.

SEASONAL VERSION: Uses quarterly time steps (4 seasons per year) to model:
- Winter spawning (broadcast spawning Dec-Feb)
- Spring larval settlement
- Summer disease peak (warmest water)
- Fall disease persistence and pre-winter mortality

Uses named sites from data/real_sites.py with:
- Distance-based larval and disease connectivity
- Sill-depth modulated disease transmission for fjords
- Freshwater lens disease reduction
- Site-specific temperatures with seasonal cycles
"""

import numpy as np
import math
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
import sys
import os

# Add parent dir to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from data.real_sites import ALL_SITES, RealSite, haversine_km
from data.coastline import generate_full_site_network

# Season constants
WINTER = 0
SPRING = 1
SUMMER = 2
FALL = 3

SEASON_NAMES = ["Winter", "Spring", "Summer", "Fall"]


@dataclass
class GeoConfig:
    """Configuration for geography-based simulation."""
    n_years: int = 100
    seasons_per_year: int = 4  # Could theoretically go monthly (12) later
    
    # Population
    base_density_per_site: int = 43000  # Base population per site (×1000; Gravem 2021: 6.1B total pre-SSWD)
    carrying_capacity_multiplier: float = 2.0
    survival_adult: float = 0.95  # Annual survival (converted to seasonal in code)
    survival_juvenile: float = 0.60  # Annual survival (converted to seasonal)
    maturation_years: int = 3  # Years to maturity (= 12 seasonal steps)
    recruitment_ratio: float = 0.35  # Per spawning event
    
    # Allee effect
    allee_threshold: int = 50
    allee_half_sat: int = 150
    
    # Disease
    disease_onset_year: int = 10  # Year 10 = 2013
    disease_base_mortality: float = 0.90  # Per year (converted to seasonal)
    disease_acute_years: int = 3
    disease_transmission_rate: float = 1.50
    
    # Larval dispersal
    larval_dispersal_scale_km: float = 50.0  # e-folding distance in km
    fjord_self_recruitment: float = 0.85  # High retention in fjords
    
    # Disease dispersal  
    disease_dispersal_scale_km: float = 200.0  # Disease spreads further than larvae
    
    # Sill effects (fjord-specific)
    sill_disease_reduction: float = 0.3  # Shallow sills reduce disease transmission
    sill_depth_threshold: float = 50.0  # Sills shallower than this reduce disease
    
    # Freshwater lens effect
    freshwater_lens_disease_reduction: float = 0.5  # 50% less disease in lens sites
    freshwater_lens_mortality_reduction: float = 0.3  # 30% less mortality (pushed to cold water)
    
    # Genetics
    n_loci: int = 10
    initial_resistance_freq: float = 0.02
    max_resistance_freq: float = 0.95
    resistance_effect: float = 0.70  # Max resistance at all loci fixed
    
    # Temperature effects
    disease_temp_threshold: float = 10.0  # °C, below this disease severity reduced
    disease_temp_optimum: float = 15.0  # °C, peak disease severity
    seasonal_temp_amplitude: float = 2.5  # °C, amplitude of seasonal variation


@dataclass
class GeoState:
    """State snapshot for one timestep (seasonal)."""
    step: int  # Absolute step number (0-319 for 80 years)
    year: int  # Calendar year
    season: int  # 0=Winter, 1=Spring, 2=Summer, 3=Fall
    populations: np.ndarray
    adults: np.ndarray
    juveniles: np.ndarray
    disease_prevalence: np.ndarray
    resistance_freqs: np.ndarray  # (n_sites, n_loci)
    temperatures: np.ndarray
    locus_effects: np.ndarray
    
    def season_name(self) -> str:
        """Return human-readable season name."""
        return SEASON_NAMES[self.season]


@dataclass 
class GeoResult:
    """Full simulation result."""
    config: GeoConfig
    sites: List[RealSite]
    states: List[GeoState]
    larval_connectivity: np.ndarray
    disease_connectivity: np.ndarray


# Latitude-based seasonal SST offsets from annual mean (°C)
# Based on NOAA buoy data and DFO station records
# [Winter, Spring, Summer, Fall]
_SEASONAL_SST_OFFSETS = [
    (61, 65, [-2.0, -0.5, 2.5, 0.5]),   # Aleutians (small range, maritime)
    (59, 61, [-2.5, -0.5, 3.0, 0.5]),   # PWS / Gulf of Alaska
    (55, 59, [-2.5, -0.5, 3.5, 0.5]),   # SE Alaska
    (52, 55, [-2.0,  0.0, 3.5, 0.5]),   # BC North Coast
    (49, 52, [-2.5,  0.0, 4.0, 0.5]),   # BC Central/South
    (47, 49, [-3.0, -0.5, 4.5, 1.0]),   # Salish Sea (semi-enclosed, larger range)
    (43, 47, [-2.5, -0.5, 3.0, 1.0]),   # WA/OR (upwelling suppresses summer)
    (39, 43, [-2.0, -0.5, 2.5, 1.0]),   # N California (strong upwelling)
    (35, 39, [-2.0, -0.5, 3.0, 0.5]),   # C California
    (29, 35, [-1.5, -0.5, 2.5, 0.5]),   # S California + Baja
    (27, 29, [-1.0, -0.5, 2.0, 0.5]),   # Southern Baja (warm, small range)
]


def _get_seasonal_offset(lat: float, season: int) -> float:
    """Get temperature offset from annual mean for latitude and season."""
    for lat_min, lat_max, offsets in _SEASONAL_SST_OFFSETS:
        if lat_min <= lat <= lat_max:
            return offsets[season]
    # Extrapolate for edge cases
    if lat < 29:
        return _SEASONAL_SST_OFFSETS[-1][2][season]
    if lat > 61:
        return _SEASONAL_SST_OFFSETS[0][2][season]
    # Interpolate between bands
    for i in range(len(_SEASONAL_SST_OFFSETS) - 1):
        lo_min, lo_max, lo_offs = _SEASONAL_SST_OFFSETS[i]
        hi_min, hi_max, hi_offs = _SEASONAL_SST_OFFSETS[i + 1]
        if hi_max <= lat < lo_min:
            frac = (lat - hi_max) / (lo_min - hi_max)
            return hi_offs[season] + frac * (lo_offs[season] - hi_offs[season])
    return 0.0


def _seasonal_temperature(base_temp: float, step: int, config: GeoConfig, 
                           year_offset: float = 0.0, lat: float = 50.0) -> float:
    """
    Calculate temperature for a given season using latitude-based SST offsets.
    
    Args:
        base_temp: Base annual mean temperature
        step: Current timestep
        config: Configuration
        year_offset: Climate warming + Blob offset (°C)
        lat: Site latitude for seasonal offset lookup
    
    Returns:
        Temperature in °C
    """
    season = step % config.seasons_per_year
    seasonal_offset = _get_seasonal_offset(lat, season)
    return base_temp + seasonal_offset + year_offset


def _temp_disease_modifier(temp: float, config: GeoConfig) -> float:
    """Temperature modifier for disease severity. Warmer = worse.
    Floor is HIGH — Hamilton shows 96% decline even in cold SE Alaska."""
    if temp >= config.disease_temp_optimum:
        return 1.0
    elif temp <= config.disease_temp_threshold - 5:
        return 0.80  # SSWD lethal at all temperatures
    else:
        return 0.80 + 0.20 * (temp - (config.disease_temp_threshold - 5)) / (config.disease_temp_optimum - config.disease_temp_threshold + 5)


def _temp_spread_modifier(temp: float, config: GeoConfig) -> float:
    """Temperature modifier for disease spread rate."""
    if temp >= config.disease_temp_optimum:
        return 1.0
    elif temp <= config.disease_temp_threshold - 5:
        return 0.4
    else:
        return 0.4 + 0.6 * (temp - (config.disease_temp_threshold - 5)) / (config.disease_temp_optimum - config.disease_temp_threshold + 5)


def build_larval_connectivity(sites: List[RealSite], config: GeoConfig) -> np.ndarray:
    """
    Build larval connectivity matrix from real geography.
    
    Key principles:
    - Distance decay (exponential kernel)
    - Fjord sites have high self-recruitment
    - Fjord→coast export reduced but nonzero (larvae are planktonic)
    - California Current creates north→south bias
    - Sill depth modulates exchange
    """
    n = len(sites)
    C = np.zeros((n, n))
    
    for i in range(n):
        src = sites[i]
        
        # First pass: compute total neighbor weights to scale self-recruitment
        neighbor_sum = 0.0
        for j in range(n):
            if i == j:
                continue
            dst = sites[j]
            dist = haversine_km(src.lat, src.lon, dst.lat, dst.lon)
            weight = math.exp(-dist / config.larval_dispersal_scale_km)
            if src.region == dst.region:
                weight *= 2.0
            neighbor_sum += weight
        
        # Set self-recruitment as ratio of neighbor sum to achieve target after normalization
        if src.site_type == "fjord":
            # Target ~85% self after normalization: self / (self + neighbors) = 0.85
            # So self = 0.85/(1-0.85) * neighbors = 5.67 * neighbors
            target_self = 0.85
        elif src.site_type == "inland_sea":
            target_self = 0.50
        else:
            target_self = 0.20
        
        raw_self = neighbor_sum * target_self / (1 - target_self)
        
        for j in range(n):
            if i == j:
                C[i, j] = raw_self
                continue
            
            dst = sites[j]
            dist = haversine_km(src.lat, src.lon, dst.lat, dst.lon)
            
            # Base dispersal kernel
            weight = math.exp(-dist / config.larval_dispersal_scale_km)
            
            # California Current bias: southward dispersal slightly favored
            lat_diff = src.lat - dst.lat  # positive = source is north of target
            if lat_diff > 0:
                # Southward (with current): slight boost
                weight *= 1.2
            elif lat_diff < -2:
                # Northward (against current): penalty for long distances
                weight *= 0.7
            
            # Fjord export reduction (sill blocks some larvae)
            if src.site_type == "fjord" and dst.site_type != "fjord":
                if src.sill_depth_m and src.sill_depth_m < config.sill_depth_threshold:
                    # Shallow sill reduces export
                    sill_factor = src.sill_depth_m / config.sill_depth_threshold
                    weight *= max(0.2, sill_factor)
                else:
                    weight *= 0.8  # Deep sill, moderate export
            
            # Fjord import (larvae entering fjord must cross sill)
            if dst.site_type == "fjord" and src.site_type != "fjord":
                if dst.sill_depth_m and dst.sill_depth_m < config.sill_depth_threshold:
                    sill_factor = dst.sill_depth_m / config.sill_depth_threshold
                    weight *= max(0.3, sill_factor)
            
            # Same region bonus
            if src.region == dst.region:
                weight *= 2.0
            
            C[i, j] = weight
    
    # Row-normalize
    for i in range(n):
        row_sum = C[i, :].sum()
        if row_sum > 0:
            C[i, :] /= row_sum
    
    return C


def build_disease_connectivity(sites: List[RealSite], config: GeoConfig) -> np.ndarray:
    """
    Build disease connectivity matrix.
    
    Disease spreads further than larvae (waterborne pathogen) but
    fjord sills and freshwater lenses impede transmission.
    """
    n = len(sites)
    D = np.zeros((n, n))
    
    for i in range(n):
        src = sites[i]
        for j in range(n):
            if i == j:
                D[i, j] = 1.0  # Self-transmission
                continue
            
            dst = sites[j]
            dist = haversine_km(src.lat, src.lon, dst.lat, dst.lon)
            
            # Wider disease dispersal kernel
            weight = math.exp(-dist / config.disease_dispersal_scale_km)
            
            # Sill blocks disease much more than larvae
            # (bacteria have 2-week transmission cycle — sill breaks the chain)
            if dst.site_type == "fjord":
                if dst.sill_depth_m and dst.sill_depth_m < config.sill_depth_threshold:
                    # Shallow sill: major disease barrier
                    sill_factor = dst.sill_depth_m / config.sill_depth_threshold
                    weight *= max(0.05, sill_factor * config.sill_disease_reduction)
                else:
                    weight *= 0.3  # Deep sill still impedes somewhat
            
            if src.site_type == "fjord":
                if src.sill_depth_m and src.sill_depth_m < config.sill_depth_threshold:
                    sill_factor = src.sill_depth_m / config.sill_depth_threshold
                    weight *= max(0.05, sill_factor * config.sill_disease_reduction)
            
            # Freshwater lens blocks waterborne pathogen
            if dst.has_freshwater_lens:
                weight *= (1 - config.freshwater_lens_disease_reduction)
            
            # Same region higher connectivity
            if src.region == dst.region:
                weight *= 1.5
            
            D[i, j] = weight
    
    # Row-normalize
    for i in range(n):
        row_sum = D[i, :].sum()
        if row_sum > 0:
            D[i, :] /= row_sum
    
    return D


class GeoSimulation:
    """Geography-based Pycnopodia population simulation with seasonal time steps."""
    
    def __init__(self, config: GeoConfig = None, sites: List[RealSite] = None, 
                 seed: int = 42, use_dense: bool = False):
        self.config = config or GeoConfig()
        if sites is not None:
            self.sites = sites
        elif use_dense:
            self.sites = generate_full_site_network()
        else:
            self.sites = list(ALL_SITES)
        self.n_sites = len(self.sites)
        self.rng = np.random.default_rng(seed)
        
        # Compute seasonal survival rates from annual rates
        # Annual survival S_annual → seasonal S_season = S_annual^(1/seasons_per_year)
        self.survival_adult_seasonal = self.config.survival_adult ** (1.0 / self.config.seasons_per_year)
        self.survival_juvenile_seasonal = self.config.survival_juvenile ** (1.0 / self.config.seasons_per_year)
        
        # Maturation: spread over seasons_per_year * maturation_years steps
        self.maturation_steps = self.config.maturation_years * self.config.seasons_per_year
        
        # Initialize populations
        self._init_populations()
        
        # Build connectivity matrices
        self.larval_connectivity = build_larval_connectivity(self.sites, self.config)
        self.disease_connectivity = build_disease_connectivity(self.sites, self.config)
        
        # Initialize disease
        self.disease_prevalence = np.zeros(self.n_sites)
        
        # Initialize genetics
        self.resistance_freqs = np.full(
            (self.n_sites, self.config.n_loci),
            self.config.initial_resistance_freq
        )
        # Add site-level variation
        self.resistance_freqs += self.rng.uniform(-0.01, 0.01, self.resistance_freqs.shape)
        self.resistance_freqs = np.clip(self.resistance_freqs, 0.001, 0.5)
        
        # Random locus effects
        raw = self.rng.dirichlet(np.ones(self.config.n_loci))
        self.locus_effects = raw * self.config.resistance_effect
        
        # Temperatures (from site data)
        self.base_temperatures = np.array([s.base_temp_C for s in self.sites])
        self.temperatures = self.base_temperatures.copy()
        
        # Larval pool for spring settlement (spawned in winter)
        self.larval_pool = np.zeros(self.n_sites)
    
    def _init_populations(self):
        """Initialize populations based on site characteristics."""
        self.populations = np.zeros(self.n_sites)
        
        for i, site in enumerate(self.sites):
            # Base density scaled by habitat quality
            base = self.config.base_density_per_site
            
            # Fjords have slightly lower density (less area)
            if site.site_type == "fjord":
                base *= 0.7
            elif site.site_type == "island":
                base *= 0.8
            
            # Temperature preference (Pycnopodia prefers cooler water)
            if site.base_temp_C < 8:
                base *= 1.2  # Cold = good
            elif site.base_temp_C > 13:
                base *= 0.6  # Warm = marginal
            
            self.populations[i] = base * 1000  # Scale to realistic numbers
        
        self.initial_populations = self.populations.copy()
        self.carrying_capacity = self.populations * self.config.carrying_capacity_multiplier
        self.adults = self.populations * 0.6
        self.juveniles = self.populations * 0.4
    
    def run(self) -> GeoResult:
        """Run full simulation."""
        result = GeoResult(
            config=self.config,
            sites=self.sites,
            states=[],
            larval_connectivity=self.larval_connectivity,
            disease_connectivity=self.disease_connectivity,
        )
        
        total_steps = self.config.n_years * self.config.seasons_per_year
        
        for step in range(total_steps):
            year = step // self.config.seasons_per_year
            season = step % self.config.seasons_per_year
            
            state = self._simulate_season(step, year, season)
            result.states.append(state)
            
            if self.populations.sum() < 1:
                # Pad remaining steps with zeros
                for s in range(step + 1, total_steps):
                    y = s // self.config.seasons_per_year
                    seas = s % self.config.seasons_per_year
                    result.states.append(GeoState(
                        step=s,
                        year=y,
                        season=seas,
                        populations=np.zeros(self.n_sites),
                        adults=np.zeros(self.n_sites),
                        juveniles=np.zeros(self.n_sites),
                        disease_prevalence=np.zeros(self.n_sites),
                        resistance_freqs=np.zeros((self.n_sites, self.config.n_loci)),
                        temperatures=self.temperatures.copy(),
                        locus_effects=self.locus_effects,
                    ))
                break
        
        return result
    
    def _simulate_season(self, step: int, year: int, season: int) -> GeoState:
        """Simulate one seasonal timestep."""
        self._current_step = step
        self._current_year = year
        self._current_season = season
        
        # 1. Natural mortality (all seasons)
        self.adults *= self.survival_adult_seasonal
        self.juveniles *= self.survival_juvenile_seasonal
        
        # 2. Temperature update (seasonal cycle + climate warming)
        self._update_temperatures(step, year, season)
        
        # 3. Disease (active all seasons, but severity varies with temperature)
        if year >= self.config.disease_onset_year:
            self._update_disease(step, year, season)
            self._apply_disease_mortality(season)
        
        # 4. Selection (stronger in summer when disease is worst)
        self._apply_selection(season)
        
        # 5. Maturation (gradual, all seasons)
        # Juveniles mature gradually: 1/maturation_steps per season
        maturing = self.juveniles / self.maturation_steps
        self.adults += maturing
        self.juveniles -= maturing
        
        # 6. Reproduction (WINTER ONLY - broadcast spawning)
        if season == WINTER:
            self.larval_pool = self._reproduce()
        
        # 7. Larval settlement (SPRING ONLY - larvae settle after winter spawn)
        if season == SPRING:
            settlers = self.larval_connectivity.T @ self.larval_pool
            settlers *= self.rng.uniform(0.8, 1.2, self.n_sites)
            
            # Density dependence on settlement
            density_effect = 1 - (self.populations / self.carrying_capacity)
            density_effect = np.clip(density_effect, 0.1, 1.0)
            settlers *= density_effect
            settlers = np.maximum(settlers, 0)
            
            self.juveniles += settlers
            
            # Update genetics from gene flow
            self._update_genetics(self.larval_pool, settlers)
            
            # Clear larval pool
            self.larval_pool = np.zeros(self.n_sites)
        
        # 8. Genetic drift
        self._apply_drift()
        
        # 9. Update totals
        self.populations = self.adults + self.juveniles
        noise = self.rng.normal(1.0, 0.03, self.n_sites)
        self.populations *= np.clip(noise, 0.9, 1.1)
        self.populations = np.maximum(self.populations, 0)
        self.populations[self.populations < 1] = 0
        self.adults[self.populations == 0] = 0
        self.juveniles[self.populations == 0] = 0
        
        total = self.populations.sum()
        if total > 0:
            adult_ratio = self.adults.sum() / (total + 1e-10)
            self.adults = self.populations * adult_ratio
            self.juveniles = self.populations * (1 - adult_ratio)
        
        return GeoState(
            step=step,
            year=year,
            season=season,
            populations=self.populations.copy(),
            adults=self.adults.copy(),
            juveniles=self.juveniles.copy(),
            disease_prevalence=self.disease_prevalence.copy(),
            resistance_freqs=self.resistance_freqs.copy(),
            temperatures=self.temperatures.copy(),
            locus_effects=self.locus_effects.copy(),
        )
    
    def _update_temperatures(self, step: int, year: int, season: int):
        """Apply seasonal cycle + climate warming + The Blob anomaly."""
        warming_rate = 0.02  # °C per year
        climate_offset = warming_rate * year
        
        # The Blob (years 10-13, strongest in summer/fall)
        blob_anomaly = 0.0
        if 10 <= year <= 12:
            if season in [SUMMER, FALL]:
                blob_anomaly = 2.5  # Peak warming in summer/fall
            else:
                blob_anomaly = 1.5  # Weaker in winter/spring
        elif year == 13:
            blob_anomaly = 0.5 if season in [SUMMER, FALL] else 0.3  # Fading
        
        # Calculate seasonal temperatures using latitude-based SST offsets
        for i in range(self.n_sites):
            self.temperatures[i] = _seasonal_temperature(
                self.base_temperatures[i], 
                step, 
                self.config, 
                climate_offset + blob_anomaly,
                lat=self.sites[i].lat
            )
    
    def _update_disease(self, step: int, year: int, season: int):
        """Update disease prevalence with seasonal dynamics."""
        years_since_onset = year - self.config.disease_onset_year
        new_prev = self.disease_prevalence.copy()
        is_acute = years_since_onset <= self.config.disease_acute_years
        
        # Seasonal disease factor
        # Winter: disease suppressed by cold
        # Spring: moderate
        # Summer: PEAK transmission
        # Fall: still high
        seasonal_disease_factor = {
            WINTER: 0.5,   # Cold suppression
            SPRING: 0.8,   # Moderate
            SUMMER: 1.3,   # PEAK
            FALL: 1.1      # Still active
        }[season]
        
        # Geographic disease spread: SSWD first detected ~47°N (WA/BC) in summer 2013
        # Spread south to California by fall 2013, north to Alaska by 2014-2015
        # Epicenter latitude
        epicenter_lat = 47.5  # Washington coast
        
        for i, site in enumerate(self.sites):
            lat = site.lat
            dist_from_epicenter = abs(lat - epicenter_lat)
            
            # Arrival delay: ~1 season per 5° latitude from epicenter
            # WA/BC (45-50°N): immediate; CA (33-42°N): 1-2 seasons; AK (55-61°N): 2-4 seasons
            arrival_delay_seasons = dist_from_epicenter / 5.0  # seasons
            # Fjords get hit later (isolation)
            if site.site_type == "fjord":
                arrival_delay_seasons += 2.0
                if site.has_freshwater_lens:
                    arrival_delay_seasons += 1.0
            
            seasons_since_onset = years_since_onset * 4 + season - SUMMER  # SUMMER of onset year = 0
            if seasons_since_onset < 0:
                continue
            
            # Has disease arrived at this site yet?
            seasons_exposed = seasons_since_onset - arrival_delay_seasons
            if seasons_exposed < 0:
                continue  # Disease hasn't reached here yet
            
            if self.disease_prevalence[i] < 0.1:
                # Initial infection at this site
                temp = self.temperatures[i]
                if temp >= 9.0:
                    base_prev = 0.80 + self.rng.uniform(0, 0.15)
                elif temp >= 7.0:
                    base_prev = 0.65 + self.rng.uniform(0, 0.15)
                else:
                    base_prev = 0.45 + self.rng.uniform(0, 0.20)
                
                # Ramp up over first 2 seasons of exposure
                if seasons_exposed < 2:
                    base_prev *= (0.5 + 0.25 * seasons_exposed)
                
                # Freshwater lens reduces initial disease establishment
                if site.has_freshwater_lens:
                    base_prev *= 0.6  # Stars pushed deeper, less exposure
                
                new_prev[i] = base_prev
                continue
        
        # Now handle ongoing disease dynamics for already-infected sites
        if True:
            for i, site in enumerate(self.sites):
                current = self.disease_prevalence[i]
                density_ratio = self.populations[i] / max(self.initial_populations[i], 1)
                
                # Fjord-specific disease dynamics
                is_fjord = site.site_type == "fjord"
                has_lens = site.has_freshwater_lens
                has_shallow_sill = (site.sill_depth_m is not None and 
                                    site.sill_depth_m < self.config.sill_depth_threshold)
                
                # Winter clearance in cold water (especially fjords)
                if season == WINTER and current > 0.1:
                    if is_fjord and has_shallow_sill and self.temperatures[i] < 6.0:
                        # Cold winter in protected fjord: disease clears significantly
                        new_prev[i] = current * 0.4
                        if new_prev[i] < 0.05:
                            new_prev[i] = 0.0
                        continue
                    elif has_lens and self.temperatures[i] < 7.0:
                        # Freshwater lens + cold: disease suppression
                        new_prev[i] = current * 0.6
                        continue
                
                if current < 0.1:
                    # Potential new infection
                    if density_ratio < 0.10:
                        new_prev[i] = 0.0
                        continue
                    
                    # Transmission pressure from neighbors
                    pressure = 0.0
                    for j in range(self.n_sites):
                        if self.disease_prevalence[j] > 0.1:
                            temp_mod = _temp_spread_modifier(self.temperatures[j], self.config)
                            pressure += (
                                self.disease_connectivity[j, i] *
                                self.disease_prevalence[j] *
                                self.config.disease_transmission_rate *
                                temp_mod *
                                seasonal_disease_factor
                            )
                    
                    infection_prob = 1.0 - math.exp(-pressure * 10.0)
                    
                    # Fjord protection during post-acute
                    if not is_acute and is_fjord and has_shallow_sill:
                        infection_prob *= 0.1  # Sill blocks reinfection
                    if has_lens:
                        infection_prob *= (1 - self.config.freshwater_lens_disease_reduction)
                    
                    if self.rng.random() < infection_prob:
                        new_prev[i] = 0.60 + self.rng.uniform(0, 0.20)
                else:
                    # Existing infection dynamics
                    if is_acute:
                        if is_fjord and has_shallow_sill:
                            # Fjord: disease decays during acute (isolation)
                            new_prev[i] = max(current * 0.90, 0.20)
                        else:
                            # Open coast: stays high
                            new_prev[i] = max(current * 0.97, 0.70)
                    else:
                        # Post-acute
                        if density_ratio < 0.10:
                            new_prev[i] = 0.0
                            continue
                        
                        if is_fjord and has_shallow_sill:
                            # Fjord: disease clears more rapidly post-acute
                            decay = 0.50 * seasonal_disease_factor
                            new_prev[i] = current * (1 - decay)
                            if new_prev[i] < 0.01:
                                new_prev[i] = 0.0
                        else:
                            # Coast: density-dependent endemic
                            neighbor_pressure = sum(
                                self.disease_connectivity[j, i] * self.disease_prevalence[j] * 0.3
                                for j in range(self.n_sites)
                                if j != i and self.disease_prevalence[j] > 0.01
                            )
                            temp_mod = _temp_spread_modifier(self.temperatures[i], self.config)
                            sustained = (density_ratio * temp_mod * 0.4 + neighbor_pressure * 0.3) * seasonal_disease_factor
                            sustained = min(sustained, 0.90)
                            
                            decay_rate = 0.15  # Per season (faster than annual 0.7)
                            target = sustained + (current - sustained) * (1 - decay_rate)
                            new_prev[i] = np.clip(
                                target + self.rng.normal(0, 0.01), 0.0, 0.95
                            )
        
        self.disease_prevalence = new_prev
    
    def _apply_disease_mortality(self, season: int):
        """Apply disease mortality with site-specific and seasonal modifiers."""
        site_resistance = self._compute_site_resistance()
        years_since = self._current_year - self.config.disease_onset_year
        is_acute = 0 <= years_since <= self.config.disease_acute_years
        
        # Seasonal mortality scaling (disease worse in summer/fall)
        seasonal_mort_factor = {
            WINTER: 0.7,   # Lower mortality in cold winter
            SPRING: 0.9,   # Moderate
            SUMMER: 1.3,   # PEAK mortality
            FALL: 1.1      # Still elevated
        }[season]
        
        # Acute phase multiplier
        acute_mult = 1.3 if is_acute else 1.0
        
        # Convert annual mortality to seasonal
        # If annual mortality is M_year, seasonal is approximately M_season = 1 - (1 - M_year)^(1/4)
        # For high mortality (0.90), this gives ~0.44 per season
        base_seasonal_mort = 1.0 - (1.0 - self.config.disease_base_mortality) ** (1.0 / self.config.seasons_per_year)
        
        for i in range(self.n_sites):
            if self.disease_prevalence[i] > 0.01:
                site = self.sites[i]
                temp_mod = _temp_disease_modifier(self.temperatures[i], self.config)
                
                effective_mort = (
                    base_seasonal_mort *
                    temp_mod *
                    self.disease_prevalence[i] *
                    (1 - site_resistance[i]) *
                    acute_mult *
                    seasonal_mort_factor
                )
                
                # Freshwater lens reduces mortality (stars go deeper into cold water)
                # Effect is weaker during acute phase (less snowmelt, Blob conditions)
                if site.has_freshwater_lens:
                    if is_acute:
                        # Gehman 2025: freshwater lens pushes stars into cold deep water
                        # Still protective during Blob, just less so
                        effective_mort *= (1 - self.config.freshwater_lens_mortality_reduction * 0.6)
                    else:
                        effective_mort *= (1 - self.config.freshwater_lens_mortality_reduction)
                
                survival = np.clip(1 - effective_mort, 0.01, 1.0)
                self.adults[i] *= survival
                self.juveniles[i] *= survival
    
    def _compute_site_resistance(self) -> np.ndarray:
        """Compute effective resistance per site from allele frequencies."""
        resistance = np.zeros(self.n_sites)
        for i in range(self.n_sites):
            for l in range(self.config.n_loci):
                resistance[i] += self.resistance_freqs[i, l] * self.locus_effects[l]
        # Cap at max possible resistance (sum of locus_effects = resistance_effect)
        np.clip(resistance, 0.0, self.config.resistance_effect, out=resistance)
        return resistance
    
    def _apply_selection(self, season: int):
        """Frequency-dependent selection on resistance loci (stronger in summer)."""
        # Selection strongest when disease is most active (summer)
        seasonal_selection_factor = {
            WINTER: 0.5,
            SPRING: 0.8,
            SUMMER: 1.5,  # Strong selection in summer
            FALL: 1.2
        }[season]
        
        for i in range(self.n_sites):
            if self.disease_prevalence[i] > 0.1:
                temp_mod = _temp_disease_modifier(self.temperatures[i], self.config)
                
                # Seasonal mortality for selection calculation
                base_seasonal_mort = 1.0 - (1.0 - self.config.disease_base_mortality) ** (1.0 / self.config.seasons_per_year)
                mortality = base_seasonal_mort * temp_mod * seasonal_selection_factor
                prevalence = self.disease_prevalence[i]
                
                for l in range(self.config.n_loci):
                    freq = self.resistance_freqs[i, l]
                    effect = self.locus_effects[l]
                    
                    w_S = 1.0 - mortality * prevalence
                    w_R = 1.0 - mortality * prevalence * (1.0 - effect)
                    
                    mean_w = freq * w_R + (1 - freq) * w_S
                    if mean_w > 0:
                        new_freq = freq * w_R / mean_w
                        self.resistance_freqs[i, l] = np.clip(
                            new_freq, 0.001, self.config.max_resistance_freq
                        )
    
    def _reproduce(self) -> np.ndarray:
        """
        Reproduction with Allee effect (WINTER ONLY - broadcast spawning).
        
        Returns larval pool that will settle in spring.
        """
        recruits = np.zeros(self.n_sites)
        for i in range(self.n_sites):
            n_adults = self.adults[i]
            if n_adults < self.config.allee_threshold:
                fertilization = 0.01
            else:
                h = self.config.allee_half_sat
                fertilization = n_adults**2 / (n_adults**2 + h**2)
            
            breeding_success = self.rng.beta(2, 20)
            recruits[i] = (
                n_adults * self.config.recruitment_ratio *
                fertilization * breeding_success / 0.091  # normalize by beta mean
            )
        return recruits
    
    def _update_genetics(self, recruits: np.ndarray, settlers: np.ndarray):
        """Update resistance frequencies from gene flow."""
        for locus in range(self.config.n_loci):
            gene_flow = self.larval_connectivity.T @ (recruits * self.resistance_freqs[:, locus])
            for j in range(self.n_sites):
                if settlers[j] > 0 and self.populations[j] > 0:
                    old_w = self.populations[j]
                    new_w = settlers[j]
                    settler_freq = gene_flow[j] / settlers[j] if settlers[j] > 0 else self.resistance_freqs[j, locus]
                    self.resistance_freqs[j, locus] = (
                        (old_w * self.resistance_freqs[j, locus] + new_w * settler_freq) /
                        (old_w + new_w)
                    )
        # Clip all frequencies to valid range
        np.clip(self.resistance_freqs, 0.001, self.config.max_resistance_freq, out=self.resistance_freqs)
    
    def _apply_drift(self):
        """Genetic drift in small populations."""
        for i in range(self.n_sites):
            if 0 < self.populations[i] < 500:
                drift_strength = 0.01 * (500 / max(self.populations[i], 1))
                for l in range(self.config.n_loci):
                    self.resistance_freqs[i, l] += self.rng.normal(0, drift_strength)
                    self.resistance_freqs[i, l] = np.clip(
                        self.resistance_freqs[i, l], 0.001, self.config.max_resistance_freq
                    )


def run_geo_ensemble(n_runs: int = 10, config: GeoConfig = None) -> List[GeoResult]:
    """Run ensemble of geography-based simulations."""
    config = config or GeoConfig()
    results = []
    for seed in range(n_runs):
        sim = GeoSimulation(config=config, seed=seed)
        results.append(sim.run())
        total_steps = config.n_years * config.seasons_per_year
        print(f"  Run {seed+1}/{n_runs} complete ({total_steps} seasonal steps)")
    return results


if __name__ == "__main__":
    print("Running geography-based Pycnopodia simulation (SEASONAL)...")
    print(f"Sites: {len(ALL_SITES)}")
    
    config = GeoConfig(n_years=80)
    total_steps = config.n_years * config.seasons_per_year
    print(f"Total timesteps: {total_steps} ({config.n_years} years × {config.seasons_per_year} seasons)")
    
    sim = GeoSimulation(config=config, seed=42)
    
    # Print connectivity check
    print("\nLarval connectivity (self-recruitment):")
    for i, site in enumerate(sim.sites):
        if site.site_type == "fjord" and site.sill_depth_m and site.sill_depth_m < 50:
            self_r = sim.larval_connectivity[i, i]
            print(f"  {site.name}: self={self_r:.3f} (sill={site.sill_depth_m}m)")
    
    result = sim.run()
    
    # Print results by region (comparing year 0 vs year 17 vs year 79)
    # Year indices in seasonal model: year N → steps N*4 through N*4+3
    print("\nResults by region:")
    regions = {}
    for i, site in enumerate(result.sites):
        if site.region not in regions:
            regions[site.region] = {"init": [], "y17": [], "y79": [], "sites": []}
        
        # Year 0 winter (step 0), Year 17 fall (step 17*4+3=71), Year 79 fall (step 319)
        step_y0 = 0
        step_y17 = 17 * 4 + 3  # End of year 17
        step_y79 = min(79 * 4 + 3, len(result.states) - 1)  # End of year 79
        
        regions[site.region]["init"].append(result.states[step_y0].populations[i])
        regions[site.region]["y17"].append(result.states[step_y17].populations[i])
        regions[site.region]["y79"].append(result.states[step_y79].populations[i])
        regions[site.region]["sites"].append(site)
    
    print(f"\n{'Region':<22} {'Sites':>5} {'Y0 Pop':>10} {'Y17 Decline':>12} {'Y79 Ratio':>10}")
    print("-" * 65)
    for region in ["se_alaska_north", "se_alaska_south", "bc_fjords", "bc_outer",
                    "salish_sea", "wa_or_outer", "n_california", "c_california", "s_california"]:
        if region in regions:
            r = regions[region]
            init = sum(r["init"])
            y17 = sum(r["y17"])
            y79 = sum(r["y79"])
            decline = (1 - y17/init) * 100 if init > 0 else 0
            ratio = y79/init if init > 0 else 0
            n_sites = len(r["sites"])
            fjords = sum(1 for s in r["sites"] if s.site_type == "fjord")
            print(f"{region:<22} {n_sites:>3}({fjords}f) {init:>10.0f} {decline:>10.1f}% {ratio:>10.4f}")
    
    # Print individual fjord sites
    print("\nFjord site details (year 79, fall):")
    step_y79 = min(79 * 4 + 3, len(result.states) - 1)
    for i, site in enumerate(result.sites):
        if site.site_type == "fjord" and site.sill_depth_m and site.sill_depth_m < 50:
            pop = result.states[step_y79].populations[i]
            init = result.states[0].populations[i]
            ratio = pop/init if init > 0 else 0
            prev = result.states[step_y79].disease_prevalence[i]
            res = result.states[step_y79].resistance_freqs[i].mean()
            print(f"  {site.name:<30} pop={pop:>8.0f} ({ratio:.3f}) prev={prev:.3f} resist={res:.4f} sill={site.sill_depth_m}m {'🌊' if site.has_freshwater_lens else ''}")
