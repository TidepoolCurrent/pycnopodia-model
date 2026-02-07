"""
Network-based metapopulation model with larval connectivity.

1000 sites × 1000 individuals = 1 million total population.
Connectivity scenarios rather than site-specific parameters.

All outputs in RATIOS relative to baseline.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Callable
from enum import Enum


class ConnectivityType(Enum):
    """Types of connectivity scenarios."""
    UNIFORM = "uniform"           # All sites equally connected
    STEPPING_STONE = "stepping"   # Only neighbors connected
    DISTANCE_DECAY = "decay"      # Connectivity decreases with distance
    ASYMMETRIC_FLOW = "asymmetric"  # Directional bias (current-like)
    HUB_NETWORK = "hub"           # Few highly-connected sites
    MODULAR = "modular"           # Clusters with weak between-links


@dataclass
class NetworkConfig:
    """Configuration for network model."""
    
    # Network structure
    n_sites: int = 1000
    n_per_site: int = 1000  # Initial individuals per site
    
    # Connectivity parameters
    connectivity_type: ConnectivityType = ConnectivityType.STEPPING_STONE
    total_connectivity: float = 0.3  # Fraction of larvae that disperse (vs self-recruit)
    asymmetry: float = 0.0  # 0 = symmetric, 1 = fully one-directional
    self_recruitment: float = 0.5  # Fraction staying at natal site
    dispersal_scale: int = 10  # Sites for distance-decay (e-folding distance)
    n_hubs: int = 10  # Number of hub sites (for hub network)
    n_modules: int = 10  # Number of clusters (for modular)
    
    # Demographics (simplified from full model)
    # Balanced for stable population without disease
    survival_adult: float = 0.90
    survival_juvenile: float = 0.50
    maturation_years: int = 3
    fecundity: float = 500  # Larvae per adult
    larval_survival: float = 0.01  # Fraction surviving to settlement
    carrying_capacity_per_site: int = 1500  # Density dependence ceiling
    
    # Allee effect
    allee_threshold: int = 50  # Below this, fertilization fails
    allee_half_sat: int = 100  # Half-saturation for fertilization
    
    # Disease
    disease_onset: int = 10
    disease_mortality: float = 0.6
    disease_peak_prevalence: float = 0.9
    disease_endemic_prevalence: float = 0.2
    disease_decay: float = 0.1
    
    # Intervention
    outplanting_per_site: int = 0  # Juveniles added per site per year
    outplanting_sites: List[int] = field(default_factory=list)  # Which sites get outplanting
    outplanting_start: int = 15
    
    # Simulation
    n_years: int = 100
    
    @property
    def total_initial_population(self) -> int:
        return self.n_sites * self.n_per_site


def build_connectivity_matrix(config: NetworkConfig) -> np.ndarray:
    """
    Build connectivity matrix C where C[i,j] = probability larvae from i settle at j.
    
    Rows sum to 1 (all larvae go somewhere or die).
    """
    n = config.n_sites
    C = np.zeros((n, n))
    
    if config.connectivity_type == ConnectivityType.UNIFORM:
        # All sites equally connected
        dispersal = (1 - config.self_recruitment) / (n - 1)
        C = np.full((n, n), dispersal)
        np.fill_diagonal(C, config.self_recruitment)
    
    elif config.connectivity_type == ConnectivityType.STEPPING_STONE:
        # Only adjacent sites connected
        for i in range(n):
            C[i, i] = config.self_recruitment
            dispersal = (1 - config.self_recruitment) / 2
            if i > 0:
                C[i, i-1] = dispersal * (1 - config.asymmetry)
            if i < n - 1:
                C[i, i+1] = dispersal * (1 + config.asymmetry)
        # Edge sites
        C[0, 0] += (1 - config.self_recruitment) / 2 * (1 - config.asymmetry)
        C[n-1, n-1] += (1 - config.self_recruitment) / 2 * (1 + config.asymmetry)
        # Normalize rows
        C = C / C.sum(axis=1, keepdims=True)
    
    elif config.connectivity_type == ConnectivityType.DISTANCE_DECAY:
        # Connectivity decreases exponentially with distance
        for i in range(n):
            for j in range(n):
                dist = abs(i - j)
                if i == j:
                    C[i, j] = config.self_recruitment
                else:
                    # Exponential decay with asymmetry
                    direction = 1.0
                    if j > i:  # "downstream"
                        direction = 1 + config.asymmetry
                    else:  # "upstream"
                        direction = 1 - config.asymmetry
                    C[i, j] = np.exp(-dist / config.dispersal_scale) * direction
            # Normalize
            C[i, :] = C[i, :] / C[i, :].sum()
            # Apply self-recruitment constraint
            C[i, i] = config.self_recruitment
            C[i, :] = C[i, :] / C[i, :].sum()
    
    elif config.connectivity_type == ConnectivityType.ASYMMETRIC_FLOW:
        # Strong directional flow (like ocean current)
        # Larvae mostly go "downstream" (increasing index)
        for i in range(n):
            C[i, i] = config.self_recruitment
            remaining = 1 - config.self_recruitment
            # Distribute remaining to downstream sites with decay
            for j in range(i+1, min(i + config.dispersal_scale * 3, n)):
                dist = j - i
                C[i, j] = remaining * np.exp(-dist / config.dispersal_scale)
            # Small amount goes upstream
            for j in range(max(0, i - config.dispersal_scale), i):
                dist = i - j
                C[i, j] = remaining * (1 - config.asymmetry) * 0.1 * np.exp(-dist / config.dispersal_scale)
        # Normalize
        C = C / C.sum(axis=1, keepdims=True)
    
    elif config.connectivity_type == ConnectivityType.HUB_NETWORK:
        # Some sites are highly connected hubs
        hub_indices = np.linspace(0, n-1, config.n_hubs, dtype=int)
        
        for i in range(n):
            C[i, i] = config.self_recruitment
            remaining = 1 - config.self_recruitment
            
            # Connect to nearest hub
            nearest_hub = hub_indices[np.argmin(np.abs(hub_indices - i))]
            
            if i in hub_indices:
                # Hub connects to all other hubs
                for h in hub_indices:
                    if h != i:
                        C[i, h] = remaining / len(hub_indices)
            else:
                # Non-hub connects to nearest hub
                C[i, nearest_hub] = remaining
        
        # Normalize
        C = C / C.sum(axis=1, keepdims=True)
    
    elif config.connectivity_type == ConnectivityType.MODULAR:
        # Clusters with high within-connectivity, low between
        module_size = n // config.n_modules
        
        for i in range(n):
            my_module = i // module_size
            C[i, i] = config.self_recruitment
            remaining = 1 - config.self_recruitment
            
            for j in range(n):
                if i == j:
                    continue
                their_module = j // module_size
                if my_module == their_module:
                    # Within-module: high connectivity
                    C[i, j] = remaining * 0.9 / (module_size - 1)
                else:
                    # Between-module: low connectivity
                    C[i, j] = remaining * 0.1 / (n - module_size)
        
        # Normalize
        C = C / C.sum(axis=1, keepdims=True)
    
    return C


@dataclass
class NetworkState:
    """State of the network at one point in time."""
    year: int
    populations: np.ndarray  # Shape: (n_sites,) - count per site
    
    # Baselines for ratio calculations
    N0_per_site: int = 1000
    N0_total: int = 1_000_000
    
    @property
    def total_population(self) -> int:
        return int(self.populations.sum())
    
    @property
    def n_ratio(self) -> float:
        """Total population as ratio of baseline."""
        return self.total_population / self.N0_total
    
    @property
    def occupied_sites(self) -> int:
        """Number of sites with population > 0."""
        return int((self.populations > 0).sum())
    
    @property
    def occupied_ratio(self) -> float:
        """Fraction of sites occupied."""
        return self.occupied_sites / len(self.populations)
    
    @property
    def site_n_ratios(self) -> np.ndarray:
        """Per-site population ratios."""
        return self.populations / self.N0_per_site
    
    @property
    def mean_site_n_ratio(self) -> float:
        """Mean population ratio across occupied sites."""
        occupied = self.populations[self.populations > 0]
        if len(occupied) == 0:
            return 0.0
        return float(occupied.mean() / self.N0_per_site)
    
    @property
    def extinct(self) -> bool:
        return self.total_population == 0
    
    def get_summary(self) -> Dict:
        return {
            "year": self.year,
            "n_ratio": self.n_ratio,
            "occupied_ratio": self.occupied_ratio,
            "mean_site_n_ratio": self.mean_site_n_ratio,
            "min_site_n_ratio": float(self.site_n_ratios.min()),
            "max_site_n_ratio": float(self.site_n_ratios.max()),
            "_total_population": self.total_population,
            "_occupied_sites": self.occupied_sites,
        }


@dataclass
class NetworkResult:
    """Results from network simulation."""
    config: NetworkConfig
    states: List[NetworkState] = field(default_factory=list)
    connectivity_matrix: np.ndarray = None
    
    @property
    def years(self) -> np.ndarray:
        return np.array([s.year for s in self.states])
    
    @property
    def n_ratio(self) -> np.ndarray:
        """Population trajectory as ratio of baseline."""
        return np.array([s.n_ratio for s in self.states])
    
    @property
    def occupied_ratio(self) -> np.ndarray:
        """Site occupancy trajectory."""
        return np.array([s.occupied_ratio for s in self.states])
    
    @property
    def extinct(self) -> bool:
        return len(self.states) > 0 and self.states[-1].extinct
    
    @property
    def final_n_ratio(self) -> float:
        return self.states[-1].n_ratio if self.states else 0.0
    
    @property
    def min_n_ratio(self) -> float:
        """Bottleneck depth."""
        if not self.states:
            return 0.0
        return min(s.n_ratio for s in self.states)
    
    @property
    def bottleneck_year(self) -> int:
        """Year of minimum population."""
        if not self.states:
            return 0
        return min(self.states, key=lambda s: s.n_ratio).year
    
    def recovered(self, threshold: float = 0.3) -> bool:
        """Whether population reached recovery threshold."""
        return self.final_n_ratio >= threshold


class NetworkSimulation:
    """
    Network-based metapopulation simulation.
    
    Each site has a population count (simplified from individual-based).
    Larvae disperse according to connectivity matrix.
    """
    
    def __init__(self, config: NetworkConfig = None, seed: int = None):
        self.config = config or NetworkConfig()
        self.rng = np.random.default_rng(seed)
        
        # Build connectivity matrix
        self.C = build_connectivity_matrix(self.config)
        
        # Initialize populations
        self.populations = np.full(self.config.n_sites, self.config.n_per_site, dtype=float)
        
        # Track adults vs juveniles (simplified age structure)
        self.adults = self.populations * 0.6  # Start with 60% adults
        self.juveniles = self.populations * 0.4
    
    def run(self) -> NetworkResult:
        """Run full simulation."""
        result = NetworkResult(
            config=self.config,
            connectivity_matrix=self.C
        )
        
        for year in range(self.config.n_years):
            state = self._simulate_year(year)
            result.states.append(state)
            
            if state.extinct:
                # Pad with extinction
                for y in range(year + 1, self.config.n_years):
                    result.states.append(NetworkState(
                        year=y,
                        populations=np.zeros(self.config.n_sites),
                        N0_per_site=self.config.n_per_site,
                        N0_total=self.config.total_initial_population
                    ))
                break
        
        return result
    
    def _simulate_year(self, year: int) -> NetworkState:
        """Simulate one year of network dynamics."""
        
        # 1. Natural mortality
        self.adults *= self.config.survival_adult
        self.juveniles *= self.config.survival_juvenile
        
        # 2. Disease mortality (after onset)
        if year >= self.config.disease_onset:
            prevalence = self._disease_prevalence(year)
            disease_survival = 1 - (prevalence * self.config.disease_mortality)
            self.adults *= disease_survival
            self.juveniles *= disease_survival
        
        # 3. Maturation (juveniles become adults)
        # Simplified: fraction of juveniles mature each year
        maturing = self.juveniles / self.config.maturation_years
        self.adults += maturing
        self.juveniles -= maturing
        
        # 4. Reproduction with Allee effect
        larvae = self._produce_larvae()
        
        # 5. Larval dispersal via connectivity matrix
        settlers = self._disperse_larvae(larvae)
        
        # 6. Add new recruits as juveniles
        self.juveniles += settlers
        
        # 7. Outplanting
        if (year >= self.config.outplanting_start and 
            self.config.outplanting_per_site > 0):
            for site in self.config.outplanting_sites:
                self.juveniles[site] += self.config.outplanting_per_site
        
        # Update total populations
        self.populations = self.adults + self.juveniles
        
        # Add stochasticity (demographic noise)
        noise = self.rng.normal(1.0, 0.05, self.config.n_sites)
        self.populations *= np.clip(noise, 0.9, 1.1)
        self.populations = np.maximum(self.populations, 0)
        
        # Split back to adults/juveniles (maintain ratio)
        total = self.populations.sum()
        if total > 0:
            adult_ratio = self.adults.sum() / (self.adults.sum() + self.juveniles.sum() + 1e-10)
            self.adults = self.populations * adult_ratio
            self.juveniles = self.populations * (1 - adult_ratio)
        
        return NetworkState(
            year=year,
            populations=self.populations.copy(),
            N0_per_site=self.config.n_per_site,
            N0_total=self.config.total_initial_population
        )
    
    def _disease_prevalence(self, year: int) -> float:
        """Calculate disease prevalence at given year."""
        if year < self.config.disease_onset:
            return 0.0
        years_since = year - self.config.disease_onset
        return (self.config.disease_endemic_prevalence + 
                (self.config.disease_peak_prevalence - self.config.disease_endemic_prevalence) * 
                np.exp(-self.config.disease_decay * years_since))
    
    def _produce_larvae(self) -> np.ndarray:
        """Produce larvae at each site with Allee effect."""
        larvae = np.zeros(self.config.n_sites)
        
        for i in range(self.config.n_sites):
            n_adults = self.adults[i]
            
            if n_adults < self.config.allee_threshold:
                # Below threshold: near-zero reproduction
                fertilization = 0.01
            else:
                # Saturating fertilization
                h = self.config.allee_half_sat
                fertilization = n_adults**2 / (n_adults**2 + h**2)
            
            # Larvae produced
            larvae[i] = (n_adults * self.config.fecundity * 
                        self.config.larval_survival * fertilization)
        
        return larvae
    
    def _disperse_larvae(self, larvae: np.ndarray) -> np.ndarray:
        """Disperse larvae according to connectivity matrix."""
        # settlers[j] = sum over i of larvae[i] * C[i,j]
        settlers = self.C.T @ larvae
        
        # Add some stochasticity
        settlers *= self.rng.uniform(0.8, 1.2, self.config.n_sites)
        
        # Density-dependent settlement (fewer survive at high density)
        K = self.config.carrying_capacity_per_site
        density_effect = 1 - (self.populations / K)
        density_effect = np.clip(density_effect, 0.1, 1.0)
        settlers *= density_effect
        
        return np.maximum(settlers, 0)


def run_scenario(
    connectivity_type: ConnectivityType,
    asymmetry: float = 0.0,
    self_recruitment: float = 0.5,
    n_replicates: int = 10,
    seed: int = 42,
    **kwargs
) -> Dict:
    """
    Run a connectivity scenario with multiple replicates.
    
    Returns summary statistics (all as ratios).
    """
    config = NetworkConfig(
        connectivity_type=connectivity_type,
        asymmetry=asymmetry,
        self_recruitment=self_recruitment,
        **kwargs
    )
    
    results = []
    for i in range(n_replicates):
        sim = NetworkSimulation(config, seed=seed + i)
        result = sim.run()
        results.append(result)
    
    # Summarize
    extinction_count = sum(1 for r in results if r.extinct)
    final_ratios = [r.final_n_ratio for r in results]
    min_ratios = [r.min_n_ratio for r in results]
    
    return {
        "connectivity_type": connectivity_type.value,
        "asymmetry": asymmetry,
        "self_recruitment": self_recruitment,
        "n_replicates": n_replicates,
        "extinction_probability": extinction_count / n_replicates,
        "mean_final_n_ratio": np.mean(final_ratios),
        "std_final_n_ratio": np.std(final_ratios),
        "mean_min_n_ratio": np.mean(min_ratios),
        "mean_bottleneck_year": np.mean([r.bottleneck_year for r in results]),
    }


def compare_scenarios(n_replicates: int = 10, seed: int = 42) -> List[Dict]:
    """Compare all connectivity scenarios."""
    scenarios = []
    
    for conn_type in ConnectivityType:
        for asymmetry in [0.0, 0.5, 0.9]:
            for self_rec in [0.3, 0.5, 0.7]:
                result = run_scenario(
                    connectivity_type=conn_type,
                    asymmetry=asymmetry,
                    self_recruitment=self_rec,
                    n_replicates=n_replicates,
                    seed=seed
                )
                scenarios.append(result)
                print(f"  {conn_type.value}, asym={asymmetry}, self={self_rec}: "
                      f"extinct={result['extinction_probability']:.0%}, "
                      f"final={result['mean_final_n_ratio']:.1%}")
    
    return scenarios
