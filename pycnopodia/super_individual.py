"""
Super-individual Pycnopodia population model.

Each super-individual represents `group_size` real stars and carries:
- Diploid genotype (n_loci × 2 alleles, 0 or 1)
- Age (in seasons)
- Disease status (susceptible/infected/recovered)
- Site assignment

This gives realistic genetics (heterozygosity, drift, inbreeding)
without tracking 6.1 billion individuals.

Target: ~60,000 super-individuals (6.1B / 100K per super-individual).
Runtime: minutes, not hours.
"""

import numpy as np
import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict
import sys
import os

from data.real_sites import ALL_SITES, RealSite, haversine_km

# Seasons
WINTER, SPRING, SUMMER, FALL = 0, 1, 2, 3


@dataclass
class SuperIndConfig:
    """Configuration for super-individual model."""
    # Time
    n_years: int = 80
    seasons_per_year: int = 4
    
    # Super-individual scaling
    group_size: int = 100_000  # Each super-individual = this many real stars
    
    # Demographics
    max_age_seasons: int = 80  # ~20 years max lifespan
    maturation_seasons: int = 12  # 3 years to maturity
    survival_adult_annual: float = 0.95
    survival_juvenile_annual: float = 0.60
    
    # Reproduction (winter only)
    fecundity: float = 0.35  # Fraction of adults producing recruits
    allee_threshold: int = 50  # Minimum adults for successful spawning (in real individuals)
    allee_half_sat: int = 150
    srs_breeding_fraction: float = 0.08  # Only 8% of adults breed (sweepstakes)
    
    # Disease
    disease_onset_year: int = 10
    disease_base_mortality: float = 0.85  # Annual (slightly lower than mean-field to compensate for stochastic variance)
    disease_acute_years: int = 3
    disease_transmission_rate: float = 1.50
    
    # Dispersal
    larval_dispersal_scale_km: float = 50.0
    disease_dispersal_scale_km: float = 200.0
    fjord_self_recruitment: float = 0.85
    
    # Sill/lens effects
    sill_disease_reduction: float = 0.5  # Stronger sill protection for stochastic model
    sill_depth_threshold: float = 50.0
    freshwater_lens_disease_reduction: float = 0.6  # Stronger lens protection 
    freshwater_lens_mortality_reduction: float = 0.4  # Gehman 2025: lens is THE key mechanism
    
    # Genetics
    n_loci: int = 50
    initial_resistance_freq: float = 0.02
    resistance_per_allele: float = 0.007  # Each resistance allele contributes this much
    # Max resistance = 50 loci × 2 alleles × 0.007 = 0.70
    
    # Temperature
    cold_temp_modifier_floor: float = 0.80
    warm_temp_modifier_ceiling: float = 1.20
    
    # Carrying capacity
    carrying_capacity_multiplier: float = 2.0


@dataclass
class SuperIndividual:
    """A super-individual representing group_size real stars."""
    site: int              # Site index
    age: int               # Age in seasons
    genotype: np.ndarray   # (n_loci, 2) diploid: 0=susceptible, 1=resistant
    infected: bool = False
    seasons_infected: int = 0
    
    @property
    def resistance(self) -> float:
        """Fraction of resistance alleles × per-allele effect."""
        return self.genotype.sum() * 0.007  # resistance_per_allele
    
    @property
    def is_mature(self) -> bool:
        return self.age >= 12  # maturation_seasons


class SuperIndSimulation:
    """Super-individual based simulation."""
    
    def __init__(self, config: SuperIndConfig = None, 
                 sites: List[RealSite] = None, seed: int = 0):
        self.config = config or SuperIndConfig()
        self.sites = sites or list(ALL_SITES)
        self.n_sites = len(self.sites)
        self.rng = np.random.default_rng(seed)
        
        # Build connectivity matrices (reuse geo_model logic)
        self._build_connectivity()
        
        # Initialize temperatures
        self.base_temperatures = np.array([s.base_temp_C for s in self.sites])
        self.temperatures = self.base_temperatures.copy()
        
        # Initialize super-individuals
        self._init_population()
        
        # Site carrying capacities (in real individuals)
        self._init_carrying_capacity()
        
        # Track disease prevalence per site
        self.site_disease_prevalence = np.zeros(self.n_sites)
    
    def _build_connectivity(self):
        """Build larval and disease connectivity matrices."""
        n = self.n_sites
        self.larval_connectivity = np.zeros((n, n))
        self.disease_connectivity = np.zeros((n, n))
        
        config = self.config
        
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                dist = haversine_km(
                    self.sites[i].lat, self.sites[i].lon,
                    self.sites[j].lat, self.sites[j].lon
                )
                
                # Larval connectivity
                larval_weight = math.exp(-dist / config.larval_dispersal_scale_km)
                if larval_weight < 0.001:
                    larval_weight = 0
                
                # Fjord sill effects on larval export
                src = self.sites[i]
                dst = self.sites[j]
                if src.site_type == "fjord" and src.sill_depth_m is not None:
                    if src.sill_depth_m < config.sill_depth_threshold:
                        larval_weight *= 0.5
                if dst.site_type == "fjord" and dst.sill_depth_m is not None:
                    if dst.sill_depth_m < config.sill_depth_threshold:
                        larval_weight *= 0.7
                if dst.has_freshwater_lens:
                    larval_weight *= (1 - config.freshwater_lens_disease_reduction * 0.3)
                
                self.larval_connectivity[i, j] = larval_weight
                
                # Disease connectivity (longer range)
                disease_weight = math.exp(-dist / config.disease_dispersal_scale_km)
                if disease_weight < 0.001:
                    disease_weight = 0
                if dst.site_type == "fjord" and dst.sill_depth_m is not None:
                    if dst.sill_depth_m < config.sill_depth_threshold:
                        disease_weight *= (1 - config.sill_disease_reduction)
                if dst.has_freshwater_lens:
                    disease_weight *= (1 - config.freshwater_lens_disease_reduction)
                
                self.disease_connectivity[i, j] = disease_weight
        
        # Self-recruitment
        for i in range(n):
            site = self.sites[i]
            if site.site_type == "fjord":
                self.larval_connectivity[i, i] = config.fjord_self_recruitment / (1 - config.fjord_self_recruitment) * self.larval_connectivity[i].sum()
            else:
                self.larval_connectivity[i, i] = self.larval_connectivity[i].sum() * 0.3
        
        # Row-normalize
        for i in range(n):
            row_sum = self.larval_connectivity[i].sum()
            if row_sum > 0:
                self.larval_connectivity[i] /= row_sum
            row_sum = self.disease_connectivity[i].sum()
            if row_sum > 0:
                self.disease_connectivity[i] /= row_sum
    
    def _init_population(self):
        """Create initial super-individuals at each site."""
        self.individuals: List[SuperIndividual] = []
        
        gs = self.config.group_size
        freq = self.config.initial_resistance_freq
        
        for i, site in enumerate(self.sites):
            # Determine number of super-individuals for this site
            base = 35000 * 1000  # Base real population
            if site.site_type == "fjord":
                base *= 0.7
            elif site.site_type == "island":
                base *= 0.8
            if site.base_temp_C < 8:
                base *= 1.2
            elif site.base_temp_C > 13:
                base *= 0.6
            
            n_super = max(1, int(base / gs))
            
            for _ in range(n_super):
                # Random age (0 to max, weighted toward younger)
                age = int(self.rng.exponential(20))  # Mean ~20 seasons = 5 years
                age = min(age, self.config.max_age_seasons - 1)
                
                # Diploid genotype: each allele independently has `freq` chance of being resistant
                genotype = (self.rng.random((self.config.n_loci, 2)) < freq).astype(np.int8)
                
                self.individuals.append(SuperIndividual(
                    site=i, age=age, genotype=genotype
                ))
        
        self.initial_pop_per_site = self._count_pop_per_site()
        print(f"  Initialized {len(self.individuals)} super-individuals "
              f"({len(self.individuals) * gs / 1e9:.2f}B real stars)", file=sys.stderr)
    
    def _init_carrying_capacity(self):
        """Set per-site carrying capacity in super-individuals."""
        gs = self.config.group_size
        self.carrying_capacity = np.zeros(self.n_sites)
        for i in range(self.n_sites):
            self.carrying_capacity[i] = self.initial_pop_per_site[i] * self.config.carrying_capacity_multiplier
    
    def _count_pop_per_site(self) -> np.ndarray:
        """Count super-individuals per site."""
        counts = np.zeros(self.n_sites)
        for ind in self.individuals:
            counts[ind.site] += 1
        return counts
    
    def _count_adults_per_site(self) -> np.ndarray:
        """Count mature super-individuals per site."""
        counts = np.zeros(self.n_sites)
        for ind in self.individuals:
            if ind.is_mature:
                counts[ind.site] += 1
        return counts
    
    def run(self):
        """Run full simulation."""
        total_steps = self.config.n_years * self.config.seasons_per_year
        
        states = []
        gs = self.config.group_size
        
        for step in range(total_steps):
            year = step // self.config.seasons_per_year
            season = step % self.config.seasons_per_year
            
            # 1. Update temperatures
            self._update_temperatures(step, year, season)
            
            # 2. Natural mortality
            self._natural_mortality(season)
            
            # 3. Aging
            for ind in self.individuals:
                ind.age += 1
            # Remove dead (age > max)
            self.individuals = [ind for ind in self.individuals 
                               if ind.age < self.config.max_age_seasons]
            
            # 4. Disease
            if year >= self.config.disease_onset_year:
                self._update_disease(year, season)
                self._disease_mortality(year, season)
            
            # 5. Reproduction (winter only)
            if season == WINTER:
                new_inds = self._reproduce()
                # Settlement in spring (add immediately for simplicity)
                self.individuals.extend(new_inds)
            
            # 6. Density regulation
            self._density_regulation()
            
            # Record state every season
            pop_per_site = self._count_pop_per_site()
            
            # Compute mean resistance per site
            resistance_per_site = np.zeros(self.n_sites)
            counts = np.zeros(self.n_sites)
            for ind in self.individuals:
                resistance_per_site[ind.site] += ind.resistance
                counts[ind.site] += 1
            for i in range(self.n_sites):
                if counts[i] > 0:
                    resistance_per_site[i] /= counts[i]
            
            states.append({
                'step': step,
                'year': year,
                'season': season,
                'populations': pop_per_site * gs,  # Scale to real numbers
                'disease_prevalence': self.site_disease_prevalence.copy(),
                'resistance': resistance_per_site,
                'n_individuals': len(self.individuals),
            })
            
            if step % 40 == 0:  # Every 10 years
                total_real = pop_per_site.sum() * gs
                print(f"  Year {year} ({2003+year}): {len(self.individuals)} super-inds, "
                      f"{total_real/1e9:.2f}B real stars", file=sys.stderr)
            
            if len(self.individuals) == 0:
                # Pad remaining
                for s in range(step + 1, total_steps):
                    states.append({
                        'step': s, 'year': s // 4, 'season': s % 4,
                        'populations': np.zeros(self.n_sites),
                        'disease_prevalence': np.zeros(self.n_sites),
                        'resistance': np.zeros(self.n_sites),
                        'n_individuals': 0,
                    })
                break
        
        return SuperIndResult(config=self.config, sites=self.sites, states=states)
    
    def _update_temperatures(self, step: int, year: int, season: int):
        """Seasonal temperatures + Blob anomaly."""
        warming_rate = 0.02
        climate_offset = warming_rate * year
        
        blob_anomaly = 0.0
        if 10 <= year <= 12:
            blob_anomaly = 2.5 if season in [SUMMER, FALL] else 1.5
        elif year == 13:
            blob_anomaly = 0.5 if season in [SUMMER, FALL] else 0.3
        
        seasonal_offsets = {WINTER: -2.0, SPRING: -0.5, SUMMER: 2.0, FALL: 0.5}
        
        for i in range(self.n_sites):
            self.temperatures[i] = (
                self.base_temperatures[i] + 
                seasonal_offsets[season] +
                climate_offset + blob_anomaly
            )
    
    def _natural_mortality(self, season: int):
        """Apply natural mortality (seasonal rate)."""
        adult_surv = self.config.survival_adult_annual ** (1/4)
        juv_surv = self.config.survival_juvenile_annual ** (1/4)
        
        survivors = []
        for ind in self.individuals:
            surv = adult_surv if ind.is_mature else juv_surv
            if self.rng.random() < surv:
                survivors.append(ind)
        self.individuals = survivors
    
    def _update_disease(self, year: int, season: int):
        """Update disease status of individuals."""
        years_since = year - self.config.disease_onset_year
        is_acute = years_since <= self.config.disease_acute_years
        epicenter_lat = 47.5
        
        # Count infected per site for prevalence
        infected_per_site = np.zeros(self.n_sites)
        total_per_site = np.zeros(self.n_sites)
        for ind in self.individuals:
            total_per_site[ind.site] += 1
            if ind.infected:
                infected_per_site[ind.site] += 1
        
        for i in range(self.n_sites):
            if total_per_site[i] > 0:
                self.site_disease_prevalence[i] = infected_per_site[i] / total_per_site[i]
            else:
                self.site_disease_prevalence[i] = 0
        
        # Geographic arrival delay
        for ind in self.individuals:
            if ind.infected:
                continue
            
            site = self.sites[ind.site]
            dist_from_epicenter = abs(site.lat - epicenter_lat)
            arrival_delay = dist_from_epicenter / 5.0
            if site.site_type == "fjord":
                arrival_delay += 2.0
                if site.has_freshwater_lens:
                    arrival_delay += 1.0
            
            seasons_since_onset = years_since * 4 + season - SUMMER
            if seasons_since_onset < arrival_delay:
                continue
            
            # Transmission pressure from neighbors
            pressure = 0.0
            for j in range(self.n_sites):
                if self.site_disease_prevalence[j] > 0.05:
                    pressure += (
                        self.disease_connectivity[j, ind.site] *
                        self.site_disease_prevalence[j] *
                        self.config.disease_transmission_rate
                    )
            
            # Initial seeding if just arrived
            if seasons_since_onset - arrival_delay < 2:
                temp = self.temperatures[ind.site]
                if temp >= 9:
                    base_prob = 0.7
                elif temp >= 7:
                    base_prob = 0.5
                else:
                    base_prob = 0.3
                if site.has_freshwater_lens:
                    base_prob *= 0.25  # Very strong lens protection (Gehman 2025)
                if site.site_type == "fjord" and site.sill_depth_m and site.sill_depth_m < self.config.sill_depth_threshold:
                    base_prob *= 0.3  # Shallow sill strongly blocks initial infection
                infection_prob = base_prob
            else:
                # Ongoing transmission
                seasonal_factor = {WINTER: 0.5, SPRING: 0.8, SUMMER: 1.3, FALL: 1.1}[season]
                infection_prob = 1.0 - math.exp(-pressure * 10.0 * seasonal_factor)
                
                if not is_acute and site.site_type == "fjord":
                    if site.sill_depth_m and site.sill_depth_m < self.config.sill_depth_threshold:
                        infection_prob *= 0.1
                if site.has_freshwater_lens:
                    infection_prob *= (1 - self.config.freshwater_lens_disease_reduction)
            
            # Resistance reduces infection probability
            infection_prob *= (1 - ind.resistance)
            
            if self.rng.random() < infection_prob:
                ind.infected = True
                ind.seasons_infected = 0
    
    def _disease_mortality(self, year: int, season: int):
        """Kill infected individuals based on disease severity."""
        years_since = year - self.config.disease_onset_year
        is_acute = years_since <= self.config.disease_acute_years
        
        # Pre-compute population per site
        total_per_site = self._count_pop_per_site()
        
        seasonal_factor = {WINTER: 0.7, SPRING: 0.9, SUMMER: 1.3, FALL: 1.1}[season]
        acute_mult = 1.3 if is_acute else 1.0
        
        base_seasonal_mort = 1.0 - (1.0 - self.config.disease_base_mortality) ** 0.25
        
        survivors = []
        for ind in self.individuals:
            if ind.infected:
                ind.seasons_infected += 1
                
                temp_mod = min(max(self.temperatures[ind.site] / 10.0, 
                                   self.config.cold_temp_modifier_floor),
                               self.config.warm_temp_modifier_ceiling)
                
                mortality = (
                    base_seasonal_mort * temp_mod * 
                    (1 - ind.resistance) * acute_mult * seasonal_factor
                )
                
                # Freshwater lens protection
                site = self.sites[ind.site]
                if site.has_freshwater_lens:
                    lens_red = self.config.freshwater_lens_mortality_reduction
                    if is_acute:
                        mortality *= (1 - lens_red * 0.7)  # Still protective during Blob
                    else:
                        mortality *= (1 - lens_red)
                
                # Disease clearance at low density
                pop = total_per_site[ind.site]
                init_pop = self.initial_pop_per_site[ind.site]
                if init_pop > 0 and pop / init_pop < 0.10 and not is_acute:
                    ind.infected = False
                    ind.seasons_infected = 0
                    survivors.append(ind)
                    continue
                
                if self.rng.random() < mortality:
                    continue  # Dead
            
            survivors.append(ind)
        self.individuals = survivors
    
    def _reproduce(self) -> List[SuperIndividual]:
        """Winter broadcast spawning with SRS and Allee effects."""
        new_individuals = []
        gs = self.config.group_size
        
        # Get adults per site
        adults_by_site: Dict[int, List[SuperIndividual]] = {}
        for ind in self.individuals:
            if ind.is_mature and not ind.infected:
                adults_by_site.setdefault(ind.site, []).append(ind)
        
        # Larval pool per source site
        larval_genotypes_by_source: Dict[int, List[np.ndarray]] = {}
        
        for site_idx, adults in adults_by_site.items():
            n_real_adults = len(adults) * gs
            
            # Allee effect
            if n_real_adults < self.config.allee_threshold:
                continue
            h = self.config.allee_half_sat
            fertilization = n_real_adults**2 / (n_real_adults**2 + h**2)
            
            # SRS: only a fraction breed
            n_breeders = max(2, int(len(adults) * self.config.srs_breeding_fraction))
            breeders = list(self.rng.choice(adults, size=min(n_breeders, len(adults)), replace=False))
            
            # Number of offspring (super-individuals)
            n_offspring = int(
                len(adults) * self.config.fecundity * fertilization *
                self.rng.beta(2, 20) / 0.091
            )
            
            offspring_genotypes = []
            for _ in range(n_offspring):
                # Pick two random breeders as parents
                if len(breeders) >= 2:
                    p1, p2 = self.rng.choice(breeders, size=2, replace=False)
                else:
                    p1 = p2 = breeders[0]
                
                # Mendelian inheritance: each parent contributes one allele per locus (vectorized)
                child_genotype = np.zeros((self.config.n_loci, 2), dtype=np.int8)
                allele_choices = self.rng.integers(2, size=(self.config.n_loci, 2))
                child_genotype[:, 0] = p1.genotype[np.arange(self.config.n_loci), allele_choices[:, 0]]
                child_genotype[:, 1] = p2.genotype[np.arange(self.config.n_loci), allele_choices[:, 1]]
                
                offspring_genotypes.append(child_genotype)
            
            if offspring_genotypes:
                larval_genotypes_by_source[site_idx] = offspring_genotypes
        
        # Larval dispersal using connectivity matrix
        for src_site, genotypes in larval_genotypes_by_source.items():
            for genotype in genotypes:
                # Determine settlement site
                probs = self.larval_connectivity[src_site]
                if probs.sum() == 0:
                    continue
                dest = self.rng.choice(self.n_sites, p=probs)
                
                new_individuals.append(SuperIndividual(
                    site=dest, age=0, genotype=genotype
                ))
        
        return new_individuals
    
    def _density_regulation(self):
        """Remove excess individuals at sites over carrying capacity."""
        pop_per_site = self._count_pop_per_site()
        
        for i in range(self.n_sites):
            if pop_per_site[i] > self.carrying_capacity[i]:
                # Random removal of excess
                site_inds = [idx for idx, ind in enumerate(self.individuals) if ind.site == i]
                n_remove = int(pop_per_site[i] - self.carrying_capacity[i])
                if n_remove > 0 and len(site_inds) > 0:
                    remove_idx = set(self.rng.choice(site_inds, size=min(n_remove, len(site_inds)), replace=False))
                    self.individuals = [ind for idx, ind in enumerate(self.individuals) if idx not in remove_idx]


@dataclass
class SuperIndResult:
    """Results from super-individual simulation."""
    config: SuperIndConfig
    sites: List[RealSite]
    states: List[dict]
    
    @property
    def n_years(self):
        return self.config.n_years
    
    def get_regional_trajectory(self, region: str) -> List[float]:
        """Get total population trajectory for a region."""
        site_indices = [i for i, s in enumerate(self.sites) if s.region == region]
        return [sum(s['populations'][i] for i in site_indices) for s in self.states]


def run_super_individual(config: SuperIndConfig = None, seed: int = 0) -> SuperIndResult:
    """Convenience function to run a single simulation."""
    sim = SuperIndSimulation(config=config, seed=seed)
    return sim.run()
