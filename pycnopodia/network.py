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
    
    # Demographics - HIGH SURVIVAL without disease (disease is the killer)
    survival_adult: float = 0.95  # Very high when healthy (>3yr)
    survival_juvenile: float = 0.60  # Reasonable juvenile survival
    maturation_years: int = 3
    carrying_capacity_per_site: int = 1500  # Density dependence ceiling
    
    # Reproduction - RATIO-BASED (we don't know absolute offspring)
    # Sweepstakes: what fraction of breeders succeed?
    breeding_success_ratio: float = 0.08  # 8% of adults breed successfully
    # Successful breeders produce enough to replace population at K
    # This is implicit - we just track the ratio of successful reproduction
    recruitment_ratio: float = 0.10  # Recruits as fraction of adult population
    
    # Allee effect
    allee_threshold: int = 50  # Below this, fertilization fails
    allee_half_sat: int = 100  # Half-saturation for fertilization
    
    # Disease - SPREADS THROUGH NETWORK
    disease_onset_site: int = 500  # Initial outbreak site (middle of network)
    disease_onset_year: int = 10
    disease_spread_rate: float = 0.3  # Probability of spreading to connected site per year
    disease_mortality: float = 0.7  # Mortality rate when infected
    disease_recovery_rate: float = 0.05  # Rate of becoming disease-free
    
    # Genetics - TRACK ALLELES PER SITE
    n_loci: int = 10  # Resistance loci to track
    initial_resistance_freq: float = 0.1  # Starting resistance allele frequency
    resistance_effect: float = 0.5  # Mortality reduction for homozygous resistant
    
    # Intervention - OUTPLANTING
    outplanting_n: int = 100  # Number of stars per outplanting event
    outplanting_sites: List[int] = field(default_factory=list)  # Which sites receive outplants
    outplanting_start: int = 15  # Year outplanting begins
    outplanting_interval: int = 1  # Years between outplanting events
    
    # Outplant genetics
    # "wild" = same resistance freq as local wild population
    # "enhanced" = unusually high resistance (e.g., from selective breeding)
    outplant_resistance_mode: str = "wild"  # "wild" or "enhanced"
    outplant_enhanced_resistance: float = 0.5  # Resistance freq if enhanced mode
    
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
    """
    State of the network at one point in time.
    
    All metrics as RATIOS for interpretability.
    """
    year: int
    populations: np.ndarray  # Shape: (n_sites,) - count per site
    
    # Genetics per site - resistance allele frequencies
    resistance_freqs: np.ndarray = None  # Shape: (n_sites,) - mean resistance freq
    heterozygosity: np.ndarray = None  # Shape: (n_sites,) - He per site
    
    # Disease per site
    disease_prevalence: np.ndarray = None  # Shape: (n_sites,) - fraction infected
    
    # Baselines for ratio calculations
    N0_per_site: int = 1000
    N0_total: int = 1_000_000
    H0: float = 0.18  # Initial He (at p=0.1: 2*0.1*0.9 = 0.18)
    
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
    def mean_resistance_freq(self) -> float:
        """Network-wide mean resistance allele frequency."""
        if self.resistance_freqs is None:
            return 0.0
        # Weight by population size
        if self.total_population == 0:
            return 0.0
        return float(np.average(self.resistance_freqs, weights=self.populations + 1e-10))
    
    @property
    def mean_heterozygosity(self) -> float:
        """Network-wide mean heterozygosity."""
        if self.heterozygosity is None:
            return 0.0
        if self.total_population == 0:
            return 0.0
        return float(np.average(self.heterozygosity, weights=self.populations + 1e-10))
    
    @property
    def h_ratio(self) -> float:
        """Heterozygosity as ratio of baseline."""
        if self.H0 == 0:
            return 0.0
        return self.mean_heterozygosity / self.H0
    
    @property
    def infected_sites_ratio(self) -> float:
        """Fraction of sites with disease present."""
        if self.disease_prevalence is None:
            return 0.0
        return float((self.disease_prevalence > 0.01).sum() / len(self.disease_prevalence))
    
    @property
    def mean_disease_prevalence(self) -> float:
        """Mean disease prevalence across infected sites."""
        if self.disease_prevalence is None:
            return 0.0
        infected = self.disease_prevalence[self.disease_prevalence > 0.01]
        if len(infected) == 0:
            return 0.0
        return float(infected.mean())
    
    @property
    def extinct(self) -> bool:
        return self.total_population == 0
    
    def get_summary(self) -> Dict:
        return {
            "year": self.year,
            # Population ratios
            "n_ratio": self.n_ratio,
            "occupied_ratio": self.occupied_ratio,
            "mean_site_n_ratio": self.mean_site_n_ratio,
            # Genetic ratios
            "resistance_freq": self.mean_resistance_freq,
            "h_ratio": self.h_ratio,
            # Disease ratios
            "infected_sites_ratio": self.infected_sites_ratio,
            "mean_disease_prevalence": self.mean_disease_prevalence,
            # Raw counts (internal)
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
    
    Tracks per-site:
    - Population (adults + juveniles)
    - Genetics (resistance allele frequency, heterozygosity)
    - Disease (prevalence, spreading through network)
    
    All outputs as RATIOS.
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
        
        # Initialize genetics per site
        # Resistance allele frequency (starts uniform, will diverge)
        self.resistance_freqs = np.full(
            self.config.n_sites, 
            self.config.initial_resistance_freq
        )
        # Add small variation
        self.resistance_freqs += self.rng.normal(0, 0.02, self.config.n_sites)
        self.resistance_freqs = np.clip(self.resistance_freqs, 0.01, 0.99)
        
        # Heterozygosity per site (He = 2pq)
        self.heterozygosity = 2 * self.resistance_freqs * (1 - self.resistance_freqs)
        self.H0 = float(self.heterozygosity.mean())  # Baseline
        
        # Disease state per site (0 = disease-free, >0 = prevalence)
        self.disease_prevalence = np.zeros(self.config.n_sites)
    
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
                        resistance_freqs=np.zeros(self.config.n_sites),
                        heterozygosity=np.zeros(self.config.n_sites),
                        disease_prevalence=np.zeros(self.config.n_sites),
                        N0_per_site=self.config.n_per_site,
                        N0_total=self.config.total_initial_population,
                        H0=self.H0
                    ))
                break
        
        return result
    
    def _simulate_year(self, year: int) -> NetworkState:
        """Simulate one year of network dynamics."""
        
        # 1. Natural mortality (HIGH survival when healthy)
        self.adults *= self.config.survival_adult
        self.juveniles *= self.config.survival_juvenile
        
        # 2. Disease dynamics - SPREADS THROUGH NETWORK
        if year >= self.config.disease_onset_year:
            self._update_disease_spread(year)
            self._apply_disease_mortality()
        
        # 3. Selection on resistance (survivors have higher resistance)
        self._apply_selection()
        
        # 4. Maturation (juveniles become adults)
        maturing = self.juveniles / self.config.maturation_years
        self.adults += maturing
        self.juveniles -= maturing
        
        # 5. Reproduction with Allee effect and sweepstakes (RATIO-BASED)
        recruits = self._reproduce()
        
        # 6. Larval dispersal via connectivity matrix (carries genes!)
        settlers, settler_genes = self._disperse_larvae(recruits)
        
        # 7. Add new recruits as juveniles
        self.juveniles += settlers
        
        # 8. Update genetics from settlers
        self._update_genetics_from_settlers(settlers, settler_genes)
        
        # 9. Genetic drift (random changes in small populations)
        self._apply_genetic_drift()
        
        # 10. Outplanting
        if (year >= self.config.outplanting_start and 
            self.config.outplanting_n > 0 and
            len(self.config.outplanting_sites) > 0 and
            (year - self.config.outplanting_start) % self.config.outplanting_interval == 0):
            self._apply_outplanting()
        
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
        
        # Update heterozygosity
        self.heterozygosity = 2 * self.resistance_freqs * (1 - self.resistance_freqs)
        
        return NetworkState(
            year=year,
            populations=self.populations.copy(),
            resistance_freqs=self.resistance_freqs.copy(),
            heterozygosity=self.heterozygosity.copy(),
            disease_prevalence=self.disease_prevalence.copy(),
            N0_per_site=self.config.n_per_site,
            N0_total=self.config.total_initial_population,
            H0=self.H0
        )
    
    def _update_disease_spread(self, year: int):
        """
        Disease spreads through network from initial site.
        Uses connectivity matrix - disease follows larvae/currents.
        """
        # Initialize disease at onset
        if year == self.config.disease_onset_year:
            # Onset site, clipped to valid range
            onset_site = min(self.config.disease_onset_site, self.config.n_sites - 1)
            self.disease_prevalence[onset_site] = 0.9  # Initial outbreak
        
        # Spread to connected sites
        new_prevalence = self.disease_prevalence.copy()
        
        for i in range(self.config.n_sites):
            if self.disease_prevalence[i] < 0.01:
                # Site is disease-free, check if neighbors are infected
                # Probability of infection from connected sites
                infection_pressure = 0.0
                for j in range(self.config.n_sites):
                    if self.disease_prevalence[j] > 0.01:
                        # Infection spreads via connectivity
                        infection_pressure += (self.C[j, i] * 
                                              self.disease_prevalence[j] * 
                                              self.config.disease_spread_rate)
                
                # Stochastic infection
                if self.rng.random() < infection_pressure:
                    new_prevalence[i] = 0.5  # New outbreak starts at 50%
            else:
                # Site is infected - prevalence dynamics
                # Can increase (more transmission) or decrease (recovery/death)
                # Tends toward equilibrium based on resistance
                mean_resistance = self.resistance_freqs[i]
                equilibrium = 0.3 * (1 - mean_resistance)  # Lower if resistant
                
                # Move toward equilibrium
                new_prevalence[i] += 0.2 * (equilibrium - self.disease_prevalence[i])
                new_prevalence[i] += self.rng.normal(0, 0.05)
                new_prevalence[i] = np.clip(new_prevalence[i], 0, 0.95)
        
        self.disease_prevalence = new_prevalence
    
    def _apply_disease_mortality(self):
        """Apply disease mortality - modulated by resistance."""
        for i in range(self.config.n_sites):
            if self.disease_prevalence[i] > 0.01:
                # Mortality depends on prevalence and resistance
                resistance = self.resistance_freqs[i]
                # Resistant individuals have lower mortality
                effective_mortality = (self.config.disease_mortality * 
                                      self.disease_prevalence[i] * 
                                      (1 - resistance * self.config.resistance_effect))
                
                survival = 1 - effective_mortality
                self.adults[i] *= survival
                self.juveniles[i] *= survival
    
    def _apply_selection(self):
        """
        Selection increases resistance allele frequency.
        Survivors of disease are more resistant on average.
        """
        for i in range(self.config.n_sites):
            if self.disease_prevalence[i] > 0.1:
                # Strong selection when disease is present
                # Increase resistance frequency (simplified)
                selection_strength = self.disease_prevalence[i] * 0.02
                self.resistance_freqs[i] += selection_strength
                self.resistance_freqs[i] = min(self.resistance_freqs[i], 0.99)
    
    def _reproduce(self) -> np.ndarray:
        """
        Reproduction with Allee effect and sweepstakes.
        
        Returns recruits as RATIO of adult population.
        We don't know absolute offspring numbers - just ratios.
        """
        recruits = np.zeros(self.config.n_sites)
        
        for i in range(self.config.n_sites):
            n_adults = self.adults[i]
            
            if n_adults < self.config.allee_threshold:
                # Below threshold: near-zero reproduction
                fertilization_ratio = 0.01
            else:
                # Saturating fertilization
                h = self.config.allee_half_sat
                fertilization_ratio = n_adults**2 / (n_adults**2 + h**2)
            
            # Sweepstakes: only a fraction of adults breed successfully
            breeding_success = self.rng.beta(2, 20)  # Mean ~0.09, high variance
            
            # Recruits as ratio of adults × fertilization × sweepstakes
            recruits[i] = (n_adults * 
                          self.config.recruitment_ratio * 
                          fertilization_ratio * 
                          breeding_success / self.config.breeding_success_ratio)
        
        return recruits
    
    def _disperse_larvae(self, recruits: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Disperse larvae according to connectivity matrix.
        Also tracks gene flow - larvae carry parental alleles.
        
        Returns (settlers, settler_resistance_freqs)
        """
        # Settlers at each site
        settlers = self.C.T @ recruits
        
        # Gene flow: weighted average of source resistance frequencies
        # settler_genes[j] = sum_i (recruits[i] * C[i,j] * resistance[i]) / settlers[j]
        gene_flow = self.C.T @ (recruits * self.resistance_freqs)
        settler_genes = np.zeros(self.config.n_sites)
        for j in range(self.config.n_sites):
            if settlers[j] > 0:
                settler_genes[j] = gene_flow[j] / settlers[j]
            else:
                settler_genes[j] = self.resistance_freqs[j]  # No change
        
        # Add stochasticity
        settlers *= self.rng.uniform(0.8, 1.2, self.config.n_sites)
        
        # Density-dependent settlement
        K = self.config.carrying_capacity_per_site
        density_effect = 1 - (self.populations / K)
        density_effect = np.clip(density_effect, 0.1, 1.0)
        settlers *= density_effect
        
        return np.maximum(settlers, 0), settler_genes
    
    def _update_genetics_from_settlers(self, settlers: np.ndarray, settler_genes: np.ndarray):
        """Update site genetics based on incoming settlers (gene flow)."""
        for i in range(self.config.n_sites):
            if settlers[i] > 0 and self.populations[i] > 0:
                # Weight by relative population sizes
                old_weight = self.populations[i]
                new_weight = settlers[i]
                total = old_weight + new_weight
                
                # Weighted average of old and new allele frequencies
                self.resistance_freqs[i] = ((old_weight * self.resistance_freqs[i] + 
                                            new_weight * settler_genes[i]) / total)
    
    def _apply_genetic_drift(self):
        """
        Random changes in allele frequency due to finite population.
        Larger effect in small populations.
        """
        for i in range(self.config.n_sites):
            n = self.populations[i]
            if n > 0:
                # Drift variance inversely proportional to population size
                # Simplified: use normal approximation
                p = self.resistance_freqs[i]
                drift_var = p * (1 - p) / (2 * max(n, 10))
                drift = self.rng.normal(0, np.sqrt(drift_var))
                self.resistance_freqs[i] = np.clip(p + drift, 0.01, 0.99)
    
    def _apply_outplanting(self):
        """
        Add outplanted stars to designated sites.
        
        Outplants can have:
        - "wild" resistance: same as local wild population
        - "enhanced" resistance: higher resistance from selective breeding
        
        Updates both population counts and allele frequencies.
        """
        n_outplants = self.config.outplanting_n
        
        for site in self.config.outplanting_sites:
            if site >= self.config.n_sites:
                continue
            
            # Current population and genetics at site
            current_n = self.populations[site]
            current_p = self.resistance_freqs[site]
            
            # Determine outplant resistance frequency
            if self.config.outplant_resistance_mode == "enhanced":
                outplant_p = self.config.outplant_enhanced_resistance
            else:  # "wild" - match local population
                outplant_p = current_p
            
            # Add outplants as juveniles
            self.juveniles[site] += n_outplants
            
            # Update allele frequency (weighted average)
            new_n = current_n + n_outplants
            if new_n > 0:
                new_p = (current_n * current_p + n_outplants * outplant_p) / new_n
                self.resistance_freqs[site] = new_p
            
            # Update heterozygosity
            self.heterozygosity[site] = 2 * self.resistance_freqs[site] * (1 - self.resistance_freqs[site])


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
