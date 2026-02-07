"""
Pacific Coast Geographic Module for Pycnopodia Network Model.

Implements realistic geography with 8 regions from SE Alaska to Southern California.
Features:
- Temperature-dependent disease dynamics
- Asymmetric connectivity (California Current dominant)
- BC Fjords as refugia (cold, isolated)
- Calibrated to observed SSWD decline patterns

Author: Subagent for Friday Harbor Labs project
Date: 2026-02-07
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
from enum import Enum


class RegionType(Enum):
    """Type of marine habitat region."""
    OUTER_COAST = "outer_coast"  # California Current influenced
    INLAND_SEA = "inland_sea"    # Semi-enclosed (Salish Sea)
    FJORD = "fjord"              # Deep, cold, isolated (refugia)


@dataclass
class RegionConfig:
    """Configuration for a single geographic region."""
    name: str
    short_name: str
    latitude_range: Tuple[float, float]  # (south, north)
    region_type: RegionType
    base_temperature: float  # Mean annual SST (°C)
    temperature_variance: float  # Seasonal/spatial variance
    warming_rate: float  # °C per decade
    historical_density: float  # Pre-SSWD relative density (0-1)
    post_sswd_survival: float  # Observed survival fraction
    n_sites: int  # Number of sites in this region
    color: str  # For plotting


# Define the 8 Pacific Coast regions
PACIFIC_COAST_REGIONS = {
    "se_alaska_south": RegionConfig(
        name="SE Alaska (Southern)",
        short_name="SE AK S",
        latitude_range=(55.0, 57.5),
        region_type=RegionType.OUTER_COAST,
        base_temperature=8.5,  # Warmer southern range
        temperature_variance=2.0,
        warming_rate=0.2,
        historical_density=0.7,
        post_sswd_survival=0.10,  # Southern range decimated
        n_sites=20,
        color='#1f77b4',  # Blue
    ),
    "se_alaska_north": RegionConfig(
        name="SE Alaska (Northern Fjords)",
        short_name="SE AK N",
        latitude_range=(57.5, 60.0),
        region_type=RegionType.FJORD,  # Northern fjords = refugia
        base_temperature=7.0,  # Coldest region
        temperature_variance=1.5,
        warming_rate=0.15,
        historical_density=0.5,
        post_sswd_survival=0.60,  # Where survivors retreated
        n_sites=20,
        color='#17becf',  # Cyan
    ),
    "bc_outer": RegionConfig(
        name="BC Outer Coast",
        short_name="BC Coast",
        latitude_range=(49.5, 52.0),  # Non-overlapping with fjords
        region_type=RegionType.OUTER_COAST,
        base_temperature=9.5,  # 8-11°C range
        temperature_variance=2.5,
        warming_rate=0.3,
        historical_density=0.8,
        post_sswd_survival=0.20,  # 80% decline
        n_sites=50,
        color='#2ca02c',  # Green
    ),
    "bc_fjords": RegionConfig(
        name="BC Fjords (Refugia)",
        short_name="BC Fjords",
        latitude_range=(52.0, 55.0),  # Northernmost BC
        region_type=RegionType.FJORD,
        base_temperature=8.0,  # 7-9°C, cold deep water
        temperature_variance=1.5,
        warming_rate=0.2,
        historical_density=0.5,
        post_sswd_survival=0.50,  # ~50% survival - REFUGIA
        n_sites=25,
        color='#17becf',  # Cyan
    ),
    "salish_sea": RegionConfig(
        name="Salish Sea / Puget Sound",
        short_name="Salish Sea",
        latitude_range=(48.5, 49.5),  # Narrow range for site ordering
        region_type=RegionType.INLAND_SEA,
        base_temperature=11.0,  # 9-13°C range, warmer
        temperature_variance=3.0,
        warming_rate=0.4,
        historical_density=1.0,  # Highest historical density
        post_sswd_survival=0.05,  # 95%+ decline
        n_sites=50,
        color='#9467bd',  # Purple
    ),
    "wa_or_outer": RegionConfig(
        name="WA/OR Outer Coast",
        short_name="WA/OR Coast",
        latitude_range=(42.0, 48.5),
        region_type=RegionType.OUTER_COAST,
        base_temperature=10.5,  # 9-12°C
        temperature_variance=2.5,
        warming_rate=0.3,
        historical_density=0.6,
        post_sswd_survival=0.10,  # 90%+ decline
        n_sites=55,
        color='#8c564b',  # Brown
    ),
    "n_california": RegionConfig(
        name="Northern California",
        short_name="N. CA",
        latitude_range=(37.0, 42.0),
        region_type=RegionType.OUTER_COAST,
        base_temperature=12.0,  # 10-14°C
        temperature_variance=3.0,
        warming_rate=0.3,
        historical_density=0.4,
        post_sswd_survival=0.01,  # ~99% decline
        n_sites=45,
        color='#e377c2',  # Pink
    ),
    "c_california": RegionConfig(
        name="Central California",
        short_name="C. CA",
        latitude_range=(34.0, 37.0),
        region_type=RegionType.OUTER_COAST,
        base_temperature=13.0,  # 11-15°C
        temperature_variance=3.0,
        warming_rate=0.35,
        historical_density=0.3,
        post_sswd_survival=0.01,  # ~99% decline
        n_sites=35,
        color='#ff7f0e',  # Orange
    ),
    "s_california": RegionConfig(
        name="Southern CA / Channel Islands",
        short_name="S. CA",
        latitude_range=(32.0, 34.0),
        region_type=RegionType.OUTER_COAST,
        base_temperature=15.5,  # 13-18°C - warmest
        temperature_variance=3.5,
        warming_rate=0.4,
        historical_density=0.2,  # Range edge, always lower
        post_sswd_survival=0.00,  # ~100% decline, disease origin
        n_sites=30,
        color='#d62728',  # Red
    ),
}

# Region order for indexing (ordered by approximate latitude, 
# with inland seas grouped after their adjacent outer coast)
# Note: Salish Sea overlaps latitudinally with WA/OR - we place it after for convenience
REGION_ORDER = [
    "s_california",
    "c_california", 
    "n_california",
    "wa_or_outer",
    "salish_sea",  # Overlaps with WA/OR latitudinally but is inland
    "bc_outer",
    "bc_fjords",   # Part of BC coast
    "se_alaska_south",
    "se_alaska_north",  # Northern fjords = refugia
]


@dataclass
class Site:
    """A single site in the network."""
    idx: int
    region_id: str
    latitude: float
    temperature: float
    initial_population: float
    is_refugia: bool = False


@dataclass
class PacificCoastConfig:
    """
    Configuration for Pacific Coast Pycnopodia network model.
    
    Generates ~330 sites across 8 regions with realistic parameters.
    """
    
    regions: Dict[str, RegionConfig] = field(
        default_factory=lambda: PACIFIC_COAST_REGIONS.copy()
    )
    
    # Temperature-disease relationship
    # Mortality increases above temperature threshold
    disease_temp_threshold: float = 12.0  # °C
    disease_temp_coefficient: float = 0.06  # Mortality increase per °C above threshold
    
    # Disease spread parameters  
    # Tuned to produce ~90% range-wide decline with temperature-dependent mortality
    # Alaska (cold) keeps ~40%, BC Fjords (refugia) ~50%, south ~0%
    disease_base_mortality: float = 0.85  # Base per-year mortality at threshold temp (acute phase)
    disease_transmission_rate: float = 1.20  # Fast coastal spread
    disease_endemic_prevalence: float = 0.03  # Low background after acute phase
    disease_acute_years: int = 3  # Acute outbreak lasts ~3 years
    
    # Disease origin - Southern CA (warmest, southernmost)
    disease_origin_region: str = "s_california"
    disease_onset_year: int = 10
    
    # Larval dispersal parameters
    larval_duration_days: Tuple[int, int] = (14, 70)
    mean_larval_duration: int = 45
    
    # California Current parameters
    california_current_strength: float = 0.6  # North→South flow dominance
    davidson_current_strength: float = 0.15   # South→North (winter, weaker)
    
    # Connectivity decay
    larval_dispersal_scale: float = 8.0  # Sites (e-folding distance)
    disease_dispersal_scale: float = 30.0  # Very wide — SSWD spread entire coast in ~2yr
    
    # Within-region vs between-region connectivity
    within_region_weight: float = 0.6
    between_region_weight: float = 0.4
    
    # Inland sea isolation
    inland_sea_isolation: float = 0.85  # 85% of larvae stay within inland sea
    fjord_isolation: float = 0.95  # 95% stay within fjords (high isolation)
    
    # Genetics
    n_loci: int = 10
    initial_resistance_freq: float = 0.02
    resistance_effect: float = 0.70
    max_resistance_freq: float = 0.95
    
    # Demographics (matching network.py)
    survival_adult: float = 0.95
    survival_juvenile: float = 0.60
    maturation_years: int = 3
    carrying_capacity_multiplier: float = 5.0  # K relative to initial pop
    recruitment_ratio: float = 0.35  # Must offset 5% adult mortality + juvenile loss
    breeding_success_ratio: float = 0.08
    allee_threshold: int = 50
    allee_half_sat: int = 100
    
    # Simulation
    n_years: int = 100
    
    @property
    def n_sites(self) -> int:
        """Total number of sites across all regions."""
        return sum(r.n_sites for r in self.regions.values())
    
    @property
    def region_order(self) -> List[str]:
        """Region IDs ordered south to north."""
        return REGION_ORDER


def build_sites(config: PacificCoastConfig, rng: np.random.Generator = None) -> List[Site]:
    """
    Build list of sites with geographic parameters.
    
    Sites are ordered south to north for easy visualization.
    """
    if rng is None:
        rng = np.random.default_rng()
    
    sites = []
    site_idx = 0
    
    for region_id in config.region_order:
        region = config.regions[region_id]
        lat_min, lat_max = region.latitude_range
        
        for i in range(region.n_sites):
            # Distribute sites evenly across latitude range
            lat = lat_min + (i + 0.5) * (lat_max - lat_min) / region.n_sites
            
            # Temperature with some variation
            temp = region.base_temperature + rng.normal(0, region.temperature_variance / 3)
            
            # Initial population based on historical density with variation
            pop = region.historical_density * rng.uniform(0.7, 1.3)
            
            # Refugia probability depends on region type
            # Fjords: many deep-water refugia; Alaska: some; outer coast: rare
            if region.region_type == RegionType.FJORD:
                refugia_prob = 0.70  # Many deep-water refugia in fjords
            elif region.short_name == "SE AK N":
                refugia_prob = 0.35  # Northern fjord refugia
            elif region.short_name == "SE AK S":
                refugia_prob = 0.15  # Some deep-water refugia
            elif region.region_type == RegionType.INLAND_SEA:
                refugia_prob = 0.06  # Rare in inland seas
            else:
                refugia_prob = 0.03  # Very rare on outer coast
            is_refugia = rng.random() < refugia_prob
            
            sites.append(Site(
                idx=site_idx,
                region_id=region_id,
                latitude=lat,
                temperature=temp,
                initial_population=pop,
                is_refugia=is_refugia,
            ))
            site_idx += 1
    
    return sites


def get_site_region_mapping(sites: List[Site]) -> Dict[str, List[int]]:
    """Get mapping from region ID to list of site indices."""
    mapping = {r: [] for r in REGION_ORDER}
    for site in sites:
        mapping[site.region_id].append(site.idx)
    return mapping


def get_site_region_boundaries(sites: List[Site]) -> Dict[str, Tuple[int, int]]:
    """Get (start_idx, end_idx) for each region."""
    boundaries = {}
    current_region = None
    start_idx = 0
    
    for i, site in enumerate(sites):
        if site.region_id != current_region:
            if current_region is not None:
                boundaries[current_region] = (start_idx, i)
            current_region = site.region_id
            start_idx = i
    
    # Don't forget the last region
    if current_region is not None:
        boundaries[current_region] = (start_idx, len(sites))
    
    return boundaries


def get_temperature_disease_modifier(temperature: float, config: PacificCoastConfig) -> float:
    """
    Calculate disease mortality modifier based on temperature.
    
    Warmer water → higher disease mortality.
    From Eisenlord et al. 2016: mortality scales with temp above threshold.
    
    Returns modifier in [1.0, ~2.0] range.
    """
    # Cold water REDUCES lethality (e.g., SE Alaska ~90% vs S. California ~99%)
    # Warm water INCREASES lethality
    temp_diff = temperature - config.disease_temp_threshold
    modifier = 1.0 + config.disease_temp_coefficient * temp_diff
    return max(0.5, min(modifier, 2.0))  # Range: 0.5x to 2.0x (cold water protective)


def get_temperature_spread_modifier(temperature: float, config: PacificCoastConfig) -> float:
    """
    Calculate disease spread rate modifier based on temperature.
    
    Warmer water → faster spread.
    
    Returns modifier in [0.5, 1.5] range.
    """
    # Cold water dramatically slows spread
    # Center around 12°C; below 9°C disease spreads very slowly
    temp_deviation = (temperature - 12.0) / 4.0  # Steeper scaling
    modifier = 1.0 + 0.6 * temp_deviation
    return np.clip(modifier, 0.2, 1.8)  # Cold = 0.2x, warm = 1.8x


def build_larval_connectivity_matrix(
    sites: List[Site],
    config: PacificCoastConfig,
    rng: np.random.Generator = None,
) -> np.ndarray:
    """
    Build larval connectivity matrix based on Pacific Coast oceanography.
    
    Features:
    - California Current: Strong north→south flow on outer coast
    - Davidson Current: Weaker south→north (winter)
    - High within-region connectivity
    - Inland seas (Salish Sea) weakly connected to outer coast
    - Fjords highly isolated (refugia)
    
    Returns C where C[i,j] = probability larvae from i settle at j.
    Rows sum to 1.
    """
    if rng is None:
        rng = np.random.default_rng()
    
    n = len(sites)
    C = np.zeros((n, n))
    
    # Get region info
    region_map = get_site_region_mapping(sites)
    
    for i, source in enumerate(sites):
        source_region = config.regions[source.region_id]
        
        # Determine isolation factor based on region type
        if source_region.region_type == RegionType.FJORD:
            self_retention = config.fjord_isolation
        elif source_region.region_type == RegionType.INLAND_SEA:
            self_retention = config.inland_sea_isolation
        else:
            self_retention = 0.5  # Outer coast has more dispersal
        
        for j, target in enumerate(sites):
            target_region = config.regions[target.region_id]
            
            # Distance in site indices (proxy for geographic distance)
            dist = abs(j - i)
            
            # Same site - self-recruitment (much higher in fjords)
            if i == j:
                if source_region.region_type == RegionType.FJORD:
                    C[i, j] = 0.92  # Fjords trap larvae — semi-enclosed, high retention
                elif source_region.region_type == RegionType.INLAND_SEA:
                    C[i, j] = 0.50  # Inland seas retain more than open coast
                else:
                    C[i, j] = 0.15  # Open coast - most larvae disperse
                continue
            
            # Same region - high connectivity
            if source.region_id == target.region_id:
                # Distance decay within region
                weight = np.exp(-dist / (config.larval_dispersal_scale * 0.5))
                C[i, j] = weight * config.within_region_weight
                continue
            
            # Different regions - apply current patterns
            
            # Check if crossing inland/outer boundary
            inland_types = {RegionType.INLAND_SEA, RegionType.FJORD}
            source_inland = source_region.region_type in inland_types
            target_inland = target_region.region_type in inland_types
            
            # Inland <-> Outer coast crossing - very low connectivity
            if source_inland != target_inland:
                # Special case: Salish Sea connects to WA/OR outer coast only
                if source.region_id == "salish_sea" and target.region_id == "wa_or_outer":
                    C[i, j] = 0.03  # Weak connection via Strait of Juan de Fuca
                elif source.region_id == "wa_or_outer" and target.region_id == "salish_sea":
                    C[i, j] = 0.02  # Even weaker inbound
                elif source.region_id == "bc_fjords" and target.region_id == "bc_outer":
                    C[i, j] = 0.02  # Fjords weakly connect to outer
                elif source.region_id == "bc_outer" and target.region_id == "bc_fjords":
                    C[i, j] = 0.01  # Very weak inbound to fjords
                else:
                    C[i, j] = 0.001  # Negligible other crossings
                continue
            
            # Outer coast to outer coast - California Current patterns
            if not source_inland and not target_inland:
                # Distance decay
                base_weight = np.exp(-dist / config.larval_dispersal_scale)
                
                # Directional bias (California Current south, Davidson north)
                if j > i:  # Northward (against California Current)
                    direction_weight = config.davidson_current_strength
                else:  # Southward (with California Current)
                    direction_weight = config.california_current_strength
                
                C[i, j] = base_weight * direction_weight * config.between_region_weight
                continue
            
            # Within inland seas
            C[i, j] = 0.01  # Low connectivity between different inland areas
        
        # Normalize row to sum to 1
        row_sum = C[i, :].sum()
        if row_sum > 0:
            C[i, :] /= row_sum
    
    return C


def build_disease_connectivity_matrix(
    sites: List[Site],
    config: PacificCoastConfig,
    rng: np.random.Generator = None,
) -> np.ndarray:
    """
    Build disease connectivity matrix.
    
    Disease spreads via waterborne pathogen - different from larval dispersal:
    - Shorter range but can spread faster
    - Less affected by currents
    - Temperature modulates spread rate
    
    Returns D where D[i,j] = transmission weight from i to j.
    Rows sum to 1.
    """
    if rng is None:
        rng = np.random.default_rng()
    
    n = len(sites)
    D = np.zeros((n, n))
    
    for i, source in enumerate(sites):
        source_region = config.regions[source.region_id]
        
        # Temperature affects spread FROM this site
        temp_mod = get_temperature_spread_modifier(source.temperature, config)
        
        for j, target in enumerate(sites):
            target_region = config.regions[target.region_id]
            
            if i == j:
                D[i, j] = 0.2  # Local persistence
                continue
            
            # Distance decay
            dist = abs(j - i)
            base_weight = np.exp(-dist / config.disease_dispersal_scale)
            
            # Same region - higher transmission
            if source.region_id == target.region_id:
                D[i, j] = base_weight * 0.8 * temp_mod
                continue
            
            # Cross-region transmission
            # Inland seas are more isolated from disease too
            inland_types = {RegionType.INLAND_SEA, RegionType.FJORD}
            source_inland = source_region.region_type in inland_types
            target_inland = target_region.region_type in inland_types
            
            # Crossing inland/outer boundary - strongly reduced transmission
            if source_inland != target_inland:
                D[i, j] = base_weight * 0.02 * temp_mod
            else:
                # Within outer coast - moderate spread
                D[i, j] = base_weight * 0.4 * temp_mod
        
        # Normalize row
        row_sum = D[i, :].sum()
        if row_sum > 0:
            D[i, :] /= row_sum
    
    return D


@dataclass
class PacificCoastState:
    """State of Pacific Coast network at one time point."""
    year: int
    populations: np.ndarray  # Shape: (n_sites,)
    resistance_freqs: np.ndarray  # Shape: (n_sites, n_loci)
    disease_prevalence: np.ndarray  # Shape: (n_sites,)
    temperatures: np.ndarray  # Shape: (n_sites,) - current temps (may vary)
    
    # Reference values
    initial_populations: np.ndarray = None
    locus_effects: np.ndarray = None
    
    def get_region_summary(
        self, 
        sites: List[Site], 
        config: PacificCoastConfig
    ) -> Dict[str, Dict]:
        """Get summary statistics by region."""
        boundaries = get_site_region_boundaries(sites)
        
        summaries = {}
        for region_id in config.region_order:
            start, end = boundaries[region_id]
            pop = self.populations[start:end]
            resist = self.resistance_freqs[start:end]
            disease = self.disease_prevalence[start:end]
            
            # Calculate mean resistance across loci, weighted by locus effects
            if self.locus_effects is not None:
                site_resistance = (resist * self.locus_effects).sum(axis=1)
                max_resist = self.locus_effects.sum()
                mean_resist = site_resistance.mean() / max_resist if max_resist > 0 else 0
            else:
                mean_resist = resist.mean()
            
            summaries[region_id] = {
                "population": pop.sum(),
                "population_ratio": pop.sum() / (self.initial_populations[start:end].sum() + 1e-10)
                    if self.initial_populations is not None else 0,
                "mean_resistance": mean_resist,
                "mean_disease_prevalence": disease.mean(),
                "sites_with_population": (pop > 0).sum(),
                "n_sites": end - start,
            }
        
        return summaries


@dataclass
class PacificCoastResult:
    """Results from Pacific Coast simulation."""
    config: PacificCoastConfig
    sites: List[Site]
    states: List[PacificCoastState]
    larval_connectivity: np.ndarray
    disease_connectivity: np.ndarray
    
    @property
    def years(self) -> np.ndarray:
        return np.array([s.year for s in self.states])
    
    @property
    def n_ratio(self) -> np.ndarray:
        """Total population ratio over time."""
        initial = self.states[0].populations.sum()
        return np.array([s.populations.sum() / initial for s in self.states])
    
    def get_region_trajectories(self) -> Dict[str, Dict[str, np.ndarray]]:
        """Get population, resistance, disease trajectories by region."""
        boundaries = get_site_region_boundaries(self.sites)
        
        trajectories = {region_id: {
            "population": [],
            "population_ratio": [],
            "resistance": [],
            "disease_prevalence": [],
        } for region_id in self.config.region_order}
        
        for state in self.states:
            summaries = state.get_region_summary(self.sites, self.config)
            for region_id, summary in summaries.items():
                trajectories[region_id]["population"].append(summary["population"])
                trajectories[region_id]["population_ratio"].append(summary["population_ratio"])
                trajectories[region_id]["resistance"].append(summary["mean_resistance"])
                trajectories[region_id]["disease_prevalence"].append(summary["mean_disease_prevalence"])
        
        # Convert to arrays
        for region_id in trajectories:
            for key in trajectories[region_id]:
                trajectories[region_id][key] = np.array(trajectories[region_id][key])
        
        return trajectories


class PacificCoastSimulation:
    """
    Pacific Coast network simulation with realistic geography.
    
    Based on NetworkSimulation but with:
    - Temperature-dependent disease mortality
    - Geographic connectivity patterns
    - Disease origin in Southern CA
    - BC Fjords as refugia
    """
    
    def __init__(
        self, 
        config: PacificCoastConfig = None, 
        seed: int = None,
    ):
        self.config = config or PacificCoastConfig()
        self.rng = np.random.default_rng(seed)
        
        # Build geographic structure
        self.sites = build_sites(self.config, self.rng)
        self.n_sites = len(self.sites)
        self.region_boundaries = get_site_region_boundaries(self.sites)
        
        # Build connectivity matrices
        self.larval_connectivity = build_larval_connectivity_matrix(
            self.sites, self.config, self.rng
        )
        self.disease_connectivity = build_disease_connectivity_matrix(
            self.sites, self.config, self.rng
        )
        
        # Initialize populations
        self.populations = np.array([s.initial_population * 1000 for s in self.sites])
        self.initial_populations = self.populations.copy()
        
        # Carrying capacity per site
        self.carrying_capacity = self.populations * self.config.carrying_capacity_multiplier
        
        # Track adults vs juveniles
        self.adults = self.populations * 0.6
        self.juveniles = self.populations * 0.4
        
        # Site temperatures (can vary seasonally)
        self.temperatures = np.array([s.temperature for s in self.sites])
        
        # Initialize polygenic resistance
        self.locus_effects = self._sample_locus_effects()
        self.resistance_freqs = np.full(
            (self.n_sites, self.config.n_loci),
            self.config.initial_resistance_freq
        )
        # Add variation
        noise = self.rng.normal(0, 0.02, self.resistance_freqs.shape)
        self.resistance_freqs += noise
        self.resistance_freqs = np.clip(self.resistance_freqs, 0.01, 0.99)
        
        # Disease state
        self.disease_prevalence = np.zeros(self.n_sites)
        
        # Mark refugia sites (fjords don't get infected)
        self.refugia_sites = set()
        for site in self.sites:
            if site.is_refugia:
                self.refugia_sites.add(site.idx)
    
    def _sample_locus_effects(self) -> np.ndarray:
        """Sample per-locus effect sizes from gamma distribution."""
        raw = self.rng.gamma(2.0, 1.0, self.config.n_loci)
        normalized = raw / raw.sum() * self.config.resistance_effect
        return normalized
    
    def _compute_site_resistance(self) -> np.ndarray:
        """Compute overall resistance per site from polygenic loci."""
        raw = (self.resistance_freqs * self.locus_effects).sum(axis=1)
        max_resist = self.locus_effects.sum()
        if max_resist > 0:
            return raw / max_resist
        return raw
    
    def run(self) -> PacificCoastResult:
        """Run full simulation."""
        result = PacificCoastResult(
            config=self.config,
            sites=self.sites,
            states=[],
            larval_connectivity=self.larval_connectivity,
            disease_connectivity=self.disease_connectivity,
        )
        
        for year in range(self.config.n_years):
            state = self._simulate_year(year)
            result.states.append(state)
            
            if self.populations.sum() < 1:
                # Extinction - pad remaining years
                for y in range(year + 1, self.config.n_years):
                    result.states.append(PacificCoastState(
                        year=y,
                        populations=np.zeros(self.n_sites),
                        resistance_freqs=np.zeros((self.n_sites, self.config.n_loci)),
                        disease_prevalence=np.zeros(self.n_sites),
                        temperatures=self.temperatures.copy(),
                        initial_populations=self.initial_populations,
                        locus_effects=self.locus_effects,
                    ))
                break
        
        return result
    
    def _simulate_year(self, year: int) -> PacificCoastState:
        """Simulate one year of dynamics."""
        
        # 1. Natural mortality
        self.adults *= self.config.survival_adult
        self.juveniles *= self.config.survival_juvenile
        
        # 2. Disease dynamics
        if year >= self.config.disease_onset_year:
            self._update_disease_spread(year)
            self._apply_disease_mortality()
        
        # 3. Selection on resistance
        self._apply_selection()
        
        # 4. Maturation
        maturing = self.juveniles / self.config.maturation_years
        self.adults += maturing
        self.juveniles -= maturing
        
        # 5. Reproduction
        recruits = self._reproduce()
        
        # 6. Larval dispersal
        settlers, settler_genes = self._disperse_larvae(recruits)
        
        # 7. Add recruits
        self.juveniles += settlers
        
        # 8. Update genetics
        self._update_genetics_from_settlers(settlers, settler_genes)
        
        # 9. Genetic drift
        self._apply_genetic_drift()
        
        # 10. Update total populations
        self.populations = self.adults + self.juveniles
        
        # Add stochasticity
        noise = self.rng.normal(1.0, 0.05, self.n_sites)
        self.populations *= np.clip(noise, 0.9, 1.1)
        self.populations = np.maximum(self.populations, 0)
        
        # Discrete extinction threshold
        self.populations[self.populations < 1] = 0
        self.adults[self.populations == 0] = 0
        self.juveniles[self.populations == 0] = 0
        
        # Maintain adult/juvenile ratio
        total = self.populations.sum()
        if total > 0:
            adult_ratio = self.adults.sum() / (self.adults.sum() + self.juveniles.sum() + 1e-10)
            self.adults = self.populations * adult_ratio
            self.juveniles = self.populations * (1 - adult_ratio)
        
        return PacificCoastState(
            year=year,
            populations=self.populations.copy(),
            resistance_freqs=self.resistance_freqs.copy(),
            disease_prevalence=self.disease_prevalence.copy(),
            temperatures=self.temperatures.copy(),
            initial_populations=self.initial_populations,
            locus_effects=self.locus_effects,
        )
    
    def _update_disease_spread(self, year: int):
        """Disease spreads from origin, temperature-modulated."""
        years_since_onset = year - self.config.disease_onset_year
        
        if years_since_onset < 0:
            return
        
        new_prevalence = self.disease_prevalence.copy()
        
        if years_since_onset == 0:
            # Disease originates in Southern CA (warmest, southernmost)
            origin_region = self.config.disease_origin_region
            start, end = self.region_boundaries[origin_region]
            
            for i in range(start, end):
                if i not in self.refugia_sites:
                    new_prevalence[i] = 0.90 + self.rng.uniform(0, 0.08)
        else:
            # Spread through disease connectivity matrix
            for i in range(self.n_sites):
                if i in self.refugia_sites:
                    continue  # Refugia stay disease-free
                
                if self.disease_prevalence[i] < 0.1:
                    # Compute transmission pressure from infected sites
                    transmission_pressure = 0.0
                    for j in range(self.n_sites):
                        if self.disease_prevalence[j] > 0.1:
                            # Temperature of source affects spread
                            temp_mod = get_temperature_spread_modifier(
                                self.temperatures[j], self.config
                            )
                            transmission_pressure += (
                                self.disease_connectivity[j, i] *
                                self.disease_prevalence[j] *
                                self.config.disease_transmission_rate *
                                temp_mod
                            )
                    
                    # Within-region rapid spread on CONNECTED coast
                    # Fjords are isolated — disease doesn't spread fast between fjord sites
                    site_region = self.sites[i].region_id
                    region_cfg = self.config.regions[site_region]
                    r_start, r_end = self.region_boundaries[site_region]
                    region_infected = any(
                        self.disease_prevalence[j] > 0.1 
                        for j in range(r_start, r_end) if j != i
                    )
                    
                    if region_infected and region_cfg.region_type != RegionType.FJORD:
                        # Open coast: disease spreads fast within region
                        infection_prob = 0.90
                    else:
                        # Fjords or between-region: use connectivity matrix
                        # Boost the base transmission for coastal adjacency
                        infection_prob = 1.0 - np.exp(-transmission_pressure * 15.0)
                    
                    if self.rng.random() < infection_prob:
                        new_prevalence[i] = 0.80 + self.rng.uniform(0, 0.15)
                else:
                    # Decay toward endemic level
                    # After acute phase, disease drops rapidly
                    years_since_onset = year - self.config.disease_onset_year
                    temp_mod = get_temperature_spread_modifier(
                        self.temperatures[i], self.config
                    )
                    endemic = self.config.disease_endemic_prevalence * temp_mod
                    current = self.disease_prevalence[i]
                    
                    if years_since_onset <= self.config.disease_acute_years:
                        # During acute phase: slow decay (disease raging)
                        decay_rate = 0.3
                    else:
                        # Post-acute: disease subsides rapidly
                        decay_rate = 0.7
                    
                    target = endemic + (current - endemic) * (1 - decay_rate)
                    new_prevalence[i] = np.clip(
                        target + self.rng.normal(0, 0.01),
                        0.0, 0.95
                    )
        
        self.disease_prevalence = new_prevalence
    
    def _apply_disease_mortality(self):
        """Apply temperature-dependent disease mortality."""
        site_resistance = self._compute_site_resistance()
        
        for i in range(self.n_sites):
            if self.disease_prevalence[i] > 0.01:
                # Temperature modifier for this site
                temp_mod = get_temperature_disease_modifier(
                    self.temperatures[i], self.config
                )
                
                # Effective mortality = base * temp_mod * prevalence * (1 - resistance)
                effective_mortality = (
                    self.config.disease_base_mortality *
                    temp_mod *
                    self.disease_prevalence[i] *
                    (1 - site_resistance[i])
                )
                
                survival = 1 - effective_mortality
                survival = np.clip(survival, 0.01, 1.0)
                
                self.adults[i] *= survival
                self.juveniles[i] *= survival
    
    def _apply_selection(self):
        """Selection via differential survival at each locus."""
        max_r = self.config.max_resistance_freq
        
        for i in range(self.n_sites):
            if self.disease_prevalence[i] > 0.1:
                prevalence = self.disease_prevalence[i]
                temp_mod = get_temperature_disease_modifier(
                    self.temperatures[i], self.config
                )
                mortality = self.config.disease_base_mortality * temp_mod
                
                for locus in range(self.config.n_loci):
                    p = self.resistance_freqs[i, locus]
                    effect = self.locus_effects[locus]
                    
                    # Fitness of susceptible vs resistant alleles
                    w_S = 1.0 - mortality * prevalence
                    w_R = 1.0 - mortality * prevalence * (1.0 - effect)
                    
                    w_bar = p * w_R + (1.0 - p) * w_S
                    
                    if w_bar > 0:
                        p_prime = p * w_R / w_bar
                        self.resistance_freqs[i, locus] = min(p_prime, max_r)
    
    def _reproduce(self) -> np.ndarray:
        """Reproduction with Allee effect."""
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
                n_adults *
                self.config.recruitment_ratio *
                fertilization *
                breeding_success / self.config.breeding_success_ratio
            )
        
        return recruits
    
    def _disperse_larvae(self, recruits: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Disperse larvae via connectivity matrix."""
        settlers = self.larval_connectivity.T @ recruits
        
        # Gene flow
        settler_genes = np.zeros((self.n_sites, self.config.n_loci))
        for locus in range(self.config.n_loci):
            gene_flow = self.larval_connectivity.T @ (recruits * self.resistance_freqs[:, locus])
            for j in range(self.n_sites):
                if settlers[j] > 0:
                    settler_genes[j, locus] = gene_flow[j] / settlers[j]
                else:
                    settler_genes[j, locus] = self.resistance_freqs[j, locus]
        
        # Stochasticity
        settlers *= self.rng.uniform(0.8, 1.2, self.n_sites)
        
        # Density dependence
        density_effect = 1 - (self.populations / self.carrying_capacity)
        density_effect = np.clip(density_effect, 0.1, 1.0)
        settlers *= density_effect
        
        return np.maximum(settlers, 0), settler_genes
    
    def _update_genetics_from_settlers(self, settlers: np.ndarray, settler_genes: np.ndarray):
        """Update genetics from gene flow."""
        for i in range(self.n_sites):
            if settlers[i] > 0 and self.populations[i] > 0:
                old_weight = self.populations[i]
                new_weight = settlers[i]
                total = old_weight + new_weight
                
                for locus in range(self.config.n_loci):
                    self.resistance_freqs[i, locus] = (
                        (old_weight * self.resistance_freqs[i, locus] +
                         new_weight * settler_genes[i, locus]) / total
                    )
    
    def _apply_genetic_drift(self):
        """Random allele frequency changes."""
        max_r = self.config.max_resistance_freq
        
        for i in range(self.n_sites):
            n = self.populations[i]
            if n > 0:
                for locus in range(self.config.n_loci):
                    p = self.resistance_freqs[i, locus]
                    drift_var = p * (1 - p) / (2 * max(n, 10))
                    drift = self.rng.normal(0, np.sqrt(drift_var))
                    self.resistance_freqs[i, locus] = np.clip(p + drift, 0.01, max_r)


def run_pacific_coast_simulation(
    config: PacificCoastConfig = None,
    seed: int = 42,
    outplanting_sites: List[int] = None,
    outplanting_n: int = 0,
    outplanting_start: int = 15,
    outplanting_resistance: float = None,
) -> PacificCoastResult:
    """
    Convenience function to run Pacific Coast simulation.
    
    Args:
        config: Configuration (uses defaults if None)
        seed: Random seed
        outplanting_sites: Site indices to receive outplants
        outplanting_n: Number of outplants per site per event
        outplanting_start: Year to start outplanting
        outplanting_resistance: Resistance freq of outplants (None = wild)
    
    Returns:
        PacificCoastResult with full trajectories
    """
    if config is None:
        config = PacificCoastConfig()
    
    sim = PacificCoastSimulation(config, seed=seed)
    
    # Store outplanting config for use in simulation
    # (Would need to extend simulation to support this)
    # For now, just run base simulation
    
    return sim.run()
