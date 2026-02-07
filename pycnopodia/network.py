"""
Network-based metapopulation model with larval connectivity.

1000 sites × 1000 individuals = 1 million total population.
Connectivity scenarios rather than site-specific parameters.

All outputs in RATIOS relative to baseline.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Tuple
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
    # Real population: ~6 billion individuals (human input)
    # Using scaled representation for computational tractability
    n_sites: int = 1000
    n_per_site: int = 1200  # Start below K (carrying_capacity_per_site=1500)
    
    # LARVAL Connectivity parameters (oceanographic currents, planktonic dispersal)
    connectivity_type: ConnectivityType = ConnectivityType.STEPPING_STONE
    total_connectivity: float = 0.3  # Fraction of larvae that disperse (vs self-recruit)
    asymmetry: float = 0.0  # 0 = symmetric, 1 = fully one-directional
    self_recruitment: float = 0.5  # Fraction staying at natal site
    dispersal_scale: int = 10  # Sites for distance-decay (e-folding distance)
    n_hubs: int = 10  # Number of hub sites (for hub network)
    n_modules: int = 10  # Number of clusters (for modular)
    
    # DISEASE Connectivity parameters (separate from larval - different mechanism)
    # SSWD spreads ~3000km in ~2 years = ~500 sites/year with 1000 sites
    # Waterborne pathogen, potentially shorter range, faster than larval dispersal
    disease_connectivity_type: ConnectivityType = ConnectivityType.DISTANCE_DECAY
    disease_dispersal_scale: int = 50  # Wider range than larvae (faster pathogen spread)
    disease_asymmetry: float = 0.3  # Slight current bias
    disease_transmission_prob: float = 0.7  # Base transmission probability per connection
    
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
    
    # Disease - SPREADS THROUGH NETWORK (CALIBRATED TO SSWD)
    # Real SSWD: ~99% mortality, spread ~3000km in 2 years (Harvell 2019)
    # Human input: "mortality closer to 99%, maybe higher"
    disease_onset_site: int = 500  # Initial outbreak site (middle of network)
    disease_onset_year: int = 10
    disease_spread_rate: float = 0.90  # Rapid wave spread (near-simultaneous)
    disease_mortality: float = 0.99  # 99% mortality - matches expert input
    disease_recovery_rate: float = 0.01  # Disease persists - minimal recovery
    disease_endemic_prevalence: float = 0.25  # Long-term endemic level
    refugia_fraction: float = 0.05  # 5% of sites are refugia (disease-free)
    
    # Genetics - POLYGENIC RESISTANCE (multiple loci per site)
    # Matches individual-based model architecture
    n_loci: int = 10  # Number of resistance loci to track
    initial_resistance_freq: float = 0.03  # 3% - allows for some survivors
    resistance_effect: float = 0.70  # Total effect when all loci at 100%
    max_resistance_freq: float = 0.95  # Biological ceiling per locus
    
    # Per-locus effect sizes (if None, sampled from gamma distribution)
    # Effect sizes are normalized so sum = resistance_effect
    locus_effects: np.ndarray = None
    locus_effect_shape: float = 2.0  # Gamma shape for sampling effect sizes
    
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
    
    def get_locus_effects(self, rng: np.random.Generator = None) -> np.ndarray:
        """
        Get per-locus effect sizes.
        
        If locus_effects is set, returns it directly.
        Otherwise, samples from gamma distribution and normalizes.
        
        Returns shape (n_loci,) with sum = resistance_effect
        """
        if self.locus_effects is not None:
            return self.locus_effects
        
        if rng is None:
            rng = np.random.default_rng()
        
        # Sample from gamma distribution
        raw_effects = rng.gamma(self.locus_effect_shape, 1.0, self.n_loci)
        # Normalize so sum = resistance_effect
        normalized = raw_effects / raw_effects.sum() * self.resistance_effect
        return normalized


def build_connectivity_matrix(
    config: NetworkConfig,
    connectivity_type: ConnectivityType = None,
    dispersal_scale: int = None,
    asymmetry: float = None,
    self_recruitment: float = None,
    n_hubs: int = None,
    n_modules: int = None,
) -> np.ndarray:
    """
    Build connectivity matrix C where C[i,j] = probability from i to j.
    
    Can be used for both larval dispersal and disease transmission.
    Optional parameters override config values to allow separate disease matrix.
    
    Rows sum to 1 (all probability accounted for).
    """
    n = config.n_sites
    C = np.zeros((n, n))
    
    # Use provided values or fall back to config
    conn_type = connectivity_type if connectivity_type is not None else config.connectivity_type
    disp_scale = dispersal_scale if dispersal_scale is not None else config.dispersal_scale
    asym = asymmetry if asymmetry is not None else config.asymmetry
    self_rec = self_recruitment if self_recruitment is not None else config.self_recruitment
    hubs = n_hubs if n_hubs is not None else config.n_hubs
    modules = n_modules if n_modules is not None else config.n_modules
    
    if conn_type == ConnectivityType.UNIFORM:
        # All sites equally connected
        dispersal = (1 - self_rec) / (n - 1)
        C = np.full((n, n), dispersal)
        np.fill_diagonal(C, self_rec)
    
    elif conn_type == ConnectivityType.STEPPING_STONE:
        # Only adjacent sites connected
        for i in range(n):
            C[i, i] = self_rec
            dispersal = (1 - self_rec) / 2
            if i > 0:
                C[i, i-1] = dispersal * (1 - asym)
            if i < n - 1:
                C[i, i+1] = dispersal * (1 + asym)
        # Edge sites
        C[0, 0] += (1 - self_rec) / 2 * (1 - asym)
        C[n-1, n-1] += (1 - self_rec) / 2 * (1 + asym)
        # Normalize rows
        C = C / C.sum(axis=1, keepdims=True)
    
    elif conn_type == ConnectivityType.DISTANCE_DECAY:
        # Connectivity decreases exponentially with distance
        for i in range(n):
            for j in range(n):
                dist = abs(i - j)
                if i == j:
                    C[i, j] = self_rec
                else:
                    # Exponential decay with asymmetry
                    direction = 1.0
                    if j > i:  # "downstream"
                        direction = 1 + asym
                    else:  # "upstream"
                        direction = 1 - asym
                    C[i, j] = np.exp(-dist / disp_scale) * direction
            # Normalize
            C[i, :] = C[i, :] / C[i, :].sum()
            # Apply self-recruitment constraint
            C[i, i] = self_rec
            C[i, :] = C[i, :] / C[i, :].sum()
    
    elif conn_type == ConnectivityType.ASYMMETRIC_FLOW:
        # Strong directional flow (like ocean current)
        # Larvae mostly go "downstream" (increasing index)
        for i in range(n):
            C[i, i] = self_rec
            remaining = 1 - self_rec
            # Distribute remaining to downstream sites with decay
            for j in range(i+1, min(i + disp_scale * 3, n)):
                dist = j - i
                C[i, j] = remaining * np.exp(-dist / disp_scale)
            # Small amount goes upstream
            for j in range(max(0, i - disp_scale), i):
                dist = i - j
                C[i, j] = remaining * (1 - asym) * 0.1 * np.exp(-dist / disp_scale)
        # Normalize
        C = C / C.sum(axis=1, keepdims=True)
    
    elif conn_type == ConnectivityType.HUB_NETWORK:
        # Some sites are highly connected hubs
        hub_indices = np.linspace(0, n-1, hubs, dtype=int)
        
        for i in range(n):
            C[i, i] = self_rec
            remaining = 1 - self_rec
            
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
    
    elif conn_type == ConnectivityType.MODULAR:
        # Clusters with high within-connectivity, low between
        module_size = n // modules
        
        for i in range(n):
            my_module = i // module_size
            C[i, i] = self_rec
            remaining = 1 - self_rec
            
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


def build_disease_connectivity_matrix(config: NetworkConfig) -> np.ndarray:
    """
    Build disease connectivity matrix D where D[i,j] = transmission weight from i to j.
    
    SEPARATE from larval connectivity - disease spreads via waterborne pathogens,
    which may have different range and directionality than planktonic larvae.
    
    Uses the disease-specific parameters from config:
    - disease_connectivity_type
    - disease_dispersal_scale  
    - disease_asymmetry
    
    For disease, we don't have self-recruitment - disease spreads from infected
    to uninfected sites (no "staying" at same site makes sense for transmission).
    We use a small self-weight to represent local transmission.
    
    Returns matrix where rows sum to 1 (normalized transmission weights).
    """
    # Build using disease-specific parameters
    # Use a small self-recruitment (0.1) since disease "staying" at a site
    # is just local transmission continuing
    D = build_connectivity_matrix(
        config,
        connectivity_type=config.disease_connectivity_type,
        dispersal_scale=config.disease_dispersal_scale,
        asymmetry=config.disease_asymmetry,
        self_recruitment=0.1,  # Low self-weight for disease spread matrix
    )
    return D


@dataclass
class NetworkState:
    """
    State of the network at one point in time.
    
    All metrics as RATIOS for interpretability.
    
    Polygenic resistance:
    - resistance_freqs has shape (n_sites, n_loci) for per-locus tracking
    - mean_resistance_freq computes weighted average across loci using locus_effects
    """
    year: int
    populations: np.ndarray  # Shape: (n_sites,) - count per site
    
    # Genetics per site - polygenic resistance allele frequencies
    resistance_freqs: np.ndarray = None  # Shape: (n_sites, n_loci) - freq per locus per site
    heterozygosity: np.ndarray = None  # Shape: (n_sites,) - mean He per site
    locus_effects: np.ndarray = None  # Shape: (n_loci,) - effect size per locus
    
    # Disease per site
    disease_prevalence: np.ndarray = None  # Shape: (n_sites,) - fraction infected
    
    # Baselines for ratio calculations
    N0_per_site: int = 1000
    N0_total: int = 1_000_000
    H0: float = 0.01  # Initial He (at p=0.005: 2*0.005*0.995 ≈ 0.01)
    
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
    
    def _compute_site_resistance(self) -> np.ndarray:
        """
        Compute overall resistance per site from polygenic loci.
        
        For each site: resistance = sum(allele_freq[locus] * locus_effect[locus])
        This represents the expected resistance score in the population.
        
        Returns shape (n_sites,)
        """
        if self.resistance_freqs is None:
            return np.zeros(len(self.populations))
        
        # Handle both 1D (legacy) and 2D (polygenic) formats
        if self.resistance_freqs.ndim == 1:
            return self.resistance_freqs
        
        # Polygenic: sum(freq * effect) across loci
        # resistance_freqs: (n_sites, n_loci), locus_effects: (n_loci,)
        if self.locus_effects is None:
            # Fall back to simple mean if no effects specified
            return self.resistance_freqs.mean(axis=1)
        
        # Expected resistance = sum of (2 * freq * effect) for diploid
        # The factor of 2 accounts for diploid individuals (max 2 alleles per locus)
        # But since we're tracking frequencies, the expected genotype value is 2*p
        return (self.resistance_freqs * self.locus_effects).sum(axis=1)
    
    @property
    def mean_resistance_freq(self) -> float:
        """
        Network-wide mean resistance frequency.
        
        For polygenic model: weighted average of per-site resistance scores,
        then normalized by total possible resistance (sum of locus effects).
        This gives a value between 0 and 1 representing overall resistance level.
        """
        if self.resistance_freqs is None:
            return 0.0
        if self.total_population == 0:
            return 0.0
        
        # Handle both 1D (legacy) and 2D (polygenic) formats
        if self.resistance_freqs.ndim == 1:
            return float(np.average(self.resistance_freqs, weights=self.populations + 1e-10))
        
        # Polygenic: compute per-site resistance scores
        site_resistance = self._compute_site_resistance()
        
        # Weight by population size
        weighted_resistance = np.average(site_resistance, weights=self.populations + 1e-10)
        
        # Normalize by max possible resistance (sum of all locus effects)
        if self.locus_effects is not None and self.locus_effects.sum() > 0:
            max_resistance = self.locus_effects.sum()
            return float(weighted_resistance / max_resistance)
        
        return float(weighted_resistance)
    
    @property
    def per_locus_mean_freq(self) -> np.ndarray:
        """
        Mean allele frequency at each locus across the network.
        
        Returns shape (n_loci,) or empty array if not polygenic.
        """
        if self.resistance_freqs is None or self.resistance_freqs.ndim == 1:
            return np.array([])
        
        if self.total_population == 0:
            return np.zeros(self.resistance_freqs.shape[1])
        
        # Weight by population size
        return np.average(self.resistance_freqs, axis=0, weights=self.populations + 1e-10)
    
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
        summary = {
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
        
        # Add per-locus frequencies if polygenic
        per_locus = self.per_locus_mean_freq
        if len(per_locus) > 0:
            summary["per_locus_freqs"] = per_locus.tolist()
        
        return summary


@dataclass
class NetworkResult:
    """Results from network simulation."""
    config: NetworkConfig
    states: List[NetworkState] = field(default_factory=list)
    connectivity_matrix: np.ndarray = None  # Larval connectivity (legacy alias)
    larval_connectivity: np.ndarray = None  # Explicit larval connectivity
    disease_connectivity: np.ndarray = None  # Disease connectivity (separate)
    
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
    - Genetics (polygenic resistance allele frequencies, heterozygosity)
    - Disease (prevalence, spreading through network)
    
    Polygenic resistance:
    - resistance_freqs has shape (n_sites, n_loci)
    - Each locus evolves independently under selection, drift, and gene flow
    - Overall resistance = sum(freq * effect) across loci
    
    All outputs as RATIOS.
    """
    
    def __init__(self, config: NetworkConfig = None, seed: int = None):
        self.config = config or NetworkConfig()
        self.rng = np.random.default_rng(seed)
        
        # Build LARVAL connectivity matrix (for dispersal)
        self.C = build_connectivity_matrix(self.config)
        self.larval_connectivity = self.C  # Alias for clarity
        
        # Build DISEASE connectivity matrix (separate mechanism)
        self.D = build_disease_connectivity_matrix(self.config)
        self.disease_connectivity = self.D  # Alias for clarity
        
        # Initialize populations
        self.populations = np.full(self.config.n_sites, self.config.n_per_site, dtype=float)
        
        # Track adults vs juveniles (simplified age structure)
        self.adults = self.populations * 0.6  # Start with 60% adults
        self.juveniles = self.populations * 0.4
        
        # Initialize per-locus effect sizes
        self.locus_effects = self.config.get_locus_effects(self.rng)
        
        # Initialize polygenic resistance - shape (n_sites, n_loci)
        # Each locus starts at initial_resistance_freq with small variation
        self.resistance_freqs = np.full(
            (self.config.n_sites, self.config.n_loci),
            self.config.initial_resistance_freq
        )
        # Add small variation per site per locus
        noise = self.rng.normal(0, 0.02, (self.config.n_sites, self.config.n_loci))
        self.resistance_freqs += noise
        self.resistance_freqs = np.clip(self.resistance_freqs, 0.01, 0.99)
        
        # Heterozygosity per site (mean He across loci)
        # He = 2pq at each locus, then average
        self._update_heterozygosity()
        self.H0 = float(self.heterozygosity.mean())  # Baseline
        
        # Disease state per site (0 = disease-free, >0 = prevalence)
        self.disease_prevalence = np.zeros(self.config.n_sites)
        
        # Refugia: some sites never get infected (cooler waters, isolation, etc.)
        n_refugia = int(self.config.n_sites * self.config.refugia_fraction)
        self.refugia_sites = set(self.rng.choice(self.config.n_sites, n_refugia, replace=False))
    
    def _update_heterozygosity(self):
        """Update heterozygosity from polygenic resistance frequencies."""
        # He = 2pq at each locus, averaged across loci
        he_per_locus = 2 * self.resistance_freqs * (1 - self.resistance_freqs)
        self.heterozygosity = he_per_locus.mean(axis=1)
    
    def _compute_site_resistance(self) -> np.ndarray:
        """
        Compute overall resistance score per site from polygenic loci.
        
        Returns shape (n_sites,) with values in [0, 1] range.
        """
        # Sum of (allele_freq * locus_effect) across loci
        raw_resistance = (self.resistance_freqs * self.locus_effects).sum(axis=1)
        # Normalize by max possible
        max_resistance = self.locus_effects.sum()
        if max_resistance > 0:
            return raw_resistance / max_resistance
        return raw_resistance
    
    def run(self) -> NetworkResult:
        """Run full simulation."""
        result = NetworkResult(
            config=self.config,
            connectivity_matrix=self.C,  # Legacy alias
            larval_connectivity=self.larval_connectivity,
            disease_connectivity=self.disease_connectivity,
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
        
        # DISCRETE EXTINCTION: sites with <1 individual go extinct
        # (Can't have 0.3 individuals)
        self.populations[self.populations < 1] = 0
        self.adults[self.populations == 0] = 0
        self.juveniles[self.populations == 0] = 0
        
        # Split back to adults/juveniles (maintain ratio)
        total = self.populations.sum()
        if total > 0:
            adult_ratio = self.adults.sum() / (self.adults.sum() + self.juveniles.sum() + 1e-10)
            self.adults = self.populations * adult_ratio
            self.juveniles = self.populations * (1 - adult_ratio)
        
        # Update heterozygosity
        self._update_heterozygosity()
        
        return NetworkState(
            year=year,
            populations=self.populations.copy(),
            resistance_freqs=self.resistance_freqs.copy(),
            heterozygosity=self.heterozygosity.copy(),
            locus_effects=self.locus_effects.copy(),
            disease_prevalence=self.disease_prevalence.copy(),
            N0_per_site=self.config.n_per_site,
            N0_total=self.config.total_initial_population,
            H0=self.H0
        )
    
    def _update_disease_spread(self, year: int):
        """
        Disease spreads through the network via DISEASE CONNECTIVITY MATRIX.
        
        Uses a SEPARATE connectivity matrix from larval dispersal because:
        - Larval dispersal = planktonic, weeks-months in water column
        - Disease transmission = waterborne pathogen, different range/mechanism
        
        Spread pattern:
        - Transmission probability to site i = sum over infected sites j:
            D[j,i] * prevalence[j] * transmission_prob
        - Creates a spatial wave pattern matching SSWD spread
        
        SSWD calibration: ~3000km in ~2 years = ~500 sites/year with 1000 sites
        """
        years_since_onset = year - self.config.disease_onset_year
        
        if years_since_onset < 0:
            return
        
        new_prevalence = self.disease_prevalence.copy()
        
        if years_since_onset == 0:
            # Year 0: disease appears at onset site (single origin for wave)
            # Use disease_onset_site parameter (default: middle of network)
            # Clamp to valid range in case n_sites < default onset site
            onset_site = min(self.config.disease_onset_site, self.config.n_sites // 2)
            if onset_site not in self.refugia_sites:
                new_prevalence[onset_site] = 0.95
            else:
                # If onset site is refugia, pick nearest non-refugia
                for offset in range(1, self.config.n_sites):
                    for candidate in [onset_site + offset, onset_site - offset]:
                        if 0 <= candidate < self.config.n_sites and candidate not in self.refugia_sites:
                            new_prevalence[candidate] = 0.95
                            break
                    else:
                        continue
                    break
        else:
            # CONNECTIVITY-BASED SPREAD
            # For each uninfected site, compute transmission probability from all infected neighbors
            # Use disease connectivity matrix D (separate from larval connectivity C)
            
            for i in range(self.config.n_sites):
                # Skip refugia - they never get infected
                if i in self.refugia_sites:
                    continue
                
                if self.disease_prevalence[i] < 0.1:
                    # Uninfected site - compute transmission probability from neighbors
                    # Sum of: D[j,i] * prevalence[j] * base_transmission_prob for all j
                    # D[j,i] = connectivity weight from j to i
                    transmission_pressure = 0.0
                    for j in range(self.config.n_sites):
                        if self.disease_prevalence[j] > 0.1:
                            # Infected neighbor contributes to transmission
                            # D.T because D[j,i] means disease flows from j to i
                            # (rows of D are sources, we want incoming to i)
                            transmission_pressure += (
                                self.disease_connectivity[j, i] *
                                self.disease_prevalence[j] *
                                self.config.disease_transmission_prob
                            )
                    
                    # Probability of infection = 1 - exp(-pressure * scaling)
                    # Calibrated so SSWD spreads ~500 sites/year with 1000 sites
                    # Higher scaling = faster spread
                    infection_prob = 1.0 - np.exp(-transmission_pressure * 20.0)
                    
                    # Stochastic infection
                    if self.rng.random() < infection_prob:
                        new_prevalence[i] = 0.85 + self.rng.uniform(0, 0.12)
                else:
                    # Already infected - slow decay toward endemic
                    endemic = self.config.disease_endemic_prevalence
                    recovery = self.config.disease_recovery_rate
                    current = self.disease_prevalence[i]
                    target = endemic + (current - endemic) * (1 - recovery)
                    new_prevalence[i] = np.clip(target + self.rng.normal(0, 0.02), 
                                                 endemic * 0.5, 0.98)
        
        self.disease_prevalence = new_prevalence
    
    def _apply_disease_mortality(self):
        """Apply disease mortality - modulated by polygenic resistance."""
        # Pre-compute site-level resistance scores
        site_resistance = self._compute_site_resistance()
        
        for i in range(self.config.n_sites):
            if self.disease_prevalence[i] > 0.01:
                # Mortality depends on prevalence and overall resistance
                resistance = site_resistance[i]
                # Resistant individuals have lower mortality
                # resistance is already normalized to [0, 1]
                effective_mortality = (self.config.disease_mortality * 
                                      self.disease_prevalence[i] * 
                                      (1 - resistance))
                
                survival = 1 - effective_mortality
                self.adults[i] *= survival
                self.juveniles[i] *= survival
    
    def _apply_selection(self):
        """
        Selection via differential survival at each locus.
        
        For polygenic resistance, selection acts on each locus independently
        based on that locus's contribution to overall resistance.
        
        Uses standard single-locus selection model per locus:
        
            p' = p * w_R / w_bar
        
        where for locus j with effect e_j:
            p = allele frequency at locus j
            w_R = fitness of resistance allele (higher survival)
            w_S = fitness of susceptible allele (lower survival)
            w_bar = p * w_R + (1-p) * w_S  (mean fitness at this locus)
        
        The selection coefficient at each locus is proportional to its effect size.
        """
        max_r = self.config.max_resistance_freq
        mortality_rate = self.config.disease_mortality
        
        for i in range(self.config.n_sites):
            if self.disease_prevalence[i] > 0.1:
                prevalence = self.disease_prevalence[i]
                
                # Apply selection to each locus independently
                for locus in range(self.config.n_loci):
                    p = self.resistance_freqs[i, locus]
                    locus_effect = self.locus_effects[locus]
                    
                    # Fitness = survival probability during disease event
                    # The effect of this locus on mortality reduction
                    # Susceptible allele: base survival
                    # Resistant allele: mortality reduced by locus_effect
                    w_S = 1.0 - mortality_rate * prevalence
                    w_R = 1.0 - mortality_rate * prevalence * (1.0 - locus_effect)
                    
                    # Mean fitness at this locus
                    w_bar = p * w_R + (1.0 - p) * w_S
                    
                    if w_bar > 0:
                        # New allele frequency after selection
                        p_prime = p * w_R / w_bar
                        self.resistance_freqs[i, locus] = min(p_prime, max_r)
                    # If w_bar == 0, everyone dies - freq unchanged
    
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
        Also tracks gene flow - larvae carry parental alleles at all loci.
        
        Returns (settlers, settler_resistance_freqs)
        where settler_resistance_freqs has shape (n_sites, n_loci)
        """
        # Settlers at each site
        settlers = self.C.T @ recruits
        
        # Gene flow: weighted average of source resistance frequencies per locus
        # settler_genes[j, locus] = sum_i (recruits[i] * C[i,j] * resistance[i, locus]) / settlers[j]
        settler_genes = np.zeros((self.config.n_sites, self.config.n_loci))
        
        for locus in range(self.config.n_loci):
            gene_flow = self.C.T @ (recruits * self.resistance_freqs[:, locus])
            for j in range(self.config.n_sites):
                if settlers[j] > 0:
                    settler_genes[j, locus] = gene_flow[j] / settlers[j]
                else:
                    settler_genes[j, locus] = self.resistance_freqs[j, locus]  # No change
        
        # Add stochasticity
        settlers *= self.rng.uniform(0.8, 1.2, self.config.n_sites)
        
        # Density-dependent settlement
        K = self.config.carrying_capacity_per_site
        density_effect = 1 - (self.populations / K)
        density_effect = np.clip(density_effect, 0.1, 1.0)
        settlers *= density_effect
        
        return np.maximum(settlers, 0), settler_genes
    
    def _update_genetics_from_settlers(self, settlers: np.ndarray, settler_genes: np.ndarray):
        """Update site genetics based on incoming settlers (gene flow at each locus)."""
        for i in range(self.config.n_sites):
            if settlers[i] > 0 and self.populations[i] > 0:
                # Weight by relative population sizes
                old_weight = self.populations[i]
                new_weight = settlers[i]
                total = old_weight + new_weight
                
                # Weighted average of old and new allele frequencies at each locus
                for locus in range(self.config.n_loci):
                    self.resistance_freqs[i, locus] = (
                        (old_weight * self.resistance_freqs[i, locus] + 
                         new_weight * settler_genes[i, locus]) / total
                    )
    
    def _apply_genetic_drift(self):
        """
        Random changes in allele frequency due to finite population.
        Larger effect in small populations. Applied independently to each locus.
        """
        max_r = self.config.max_resistance_freq
        for i in range(self.config.n_sites):
            n = self.populations[i]
            if n > 0:
                # Drift variance inversely proportional to population size
                # Simplified: use normal approximation
                # Apply drift to each locus independently
                for locus in range(self.config.n_loci):
                    p = self.resistance_freqs[i, locus]
                    drift_var = p * (1 - p) / (2 * max(n, 10))
                    drift = self.rng.normal(0, np.sqrt(drift_var))
                    self.resistance_freqs[i, locus] = np.clip(p + drift, 0.01, max_r)
    
    def _apply_outplanting(self):
        """
        Add outplanted stars to designated sites.
        
        Outplants can have:
        - "wild" resistance: same as local wild population at each locus
        - "enhanced" resistance: higher resistance from selective breeding
        
        For enhanced mode, the outplant_enhanced_resistance is applied uniformly
        across all loci (representing consistent selective breeding).
        
        Updates both population counts and allele frequencies at each locus.
        """
        n_outplants = self.config.outplanting_n
        
        for site in self.config.outplanting_sites:
            if site >= self.config.n_sites:
                continue
            
            # Current population at site
            current_n = self.populations[site]
            
            # Add outplants as juveniles
            self.juveniles[site] += n_outplants
            
            # Update allele frequency at each locus (weighted average)
            new_n = current_n + n_outplants
            if new_n > 0:
                for locus in range(self.config.n_loci):
                    current_p = self.resistance_freqs[site, locus]
                    
                    # Determine outplant resistance frequency for this locus
                    if self.config.outplant_resistance_mode == "enhanced":
                        # Enhanced: all loci have high resistance
                        outplant_p = self.config.outplant_enhanced_resistance
                    else:  # "wild" - match local population at this locus
                        outplant_p = current_p
                    
                    # Weighted average
                    new_p = (current_n * current_p + n_outplants * outplant_p) / new_n
                    self.resistance_freqs[site, locus] = new_p
        
        # Update heterozygosity
        self._update_heterozygosity()


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
