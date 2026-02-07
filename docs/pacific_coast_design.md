# Pacific Coast Geographic Configuration for Pycnopodia Network Model

**Author:** Research compiled for Friday Harbor Labs PhD project
**Date:** 2026-02-07
**Purpose:** Design document for spatially realistic network model based on real Pacific coast geography

---

## Executive Summary

This document provides the geographic foundation for building a Pycnopodia helianthoides metapopulation model that reflects real-world Pacific coast geography. Instead of abstract 1000-site grids, we define:

1. **8 distinct geographic regions** from Alaska to Southern California
2. **Asymmetric connectivity** driven by California Current system
3. **Temperature gradients** that affect disease severity
4. **Calibration targets** from observed SSWD decline patterns (2013-2017)

Key insight: **Fjord habitats in British Columbia act as refugia** due to cold temperatures and low connectivity—this should be captured in the model.

---

## Geographic Regions

### Outer Coast Network
Connected by California Current (southward) and Davidson Current (northward, winter):

| Region | Lat Range | Mean SST | Post-SSWD Survival |
|--------|-----------|----------|-------------------|
| SE Alaska | 55-60°N | 7°C | 40% |
| BC Outer Coast | 48.5-55°N | 9°C | 15% |
| WA Outer Coast | 46-48.5°N | 10°C | 1% |
| Oregon Coast | 42-46°N | 11°C | 0% |
| N. California | 37-42°N | 12°C | 0% |
| S. California | 32-37°N | 15°C | 0% |

### Inland Sea / Fjord Networks
Semi-enclosed with limited connectivity to outer coast:

| Region | Type | Mean SST | Post-SSWD Survival |
|--------|------|----------|-------------------|
| Salish Sea / Puget Sound | Inland sea | 11°C | 1% |
| BC Fjords | Fjord (refugia) | 7°C | 50% |

---

## Connectivity Matrix

Based on larval duration (2-10 weeks) and current systems:

```
                 From ↓  To →
                 SE_AK  BC_O  BC_F  SAL   WA_O  OR    N_CA  S_CA
SE Alaska        0.60   0.25  0.05  0.02  -     -     -     -
BC Outer         0.10   0.50  0.10  0.05  0.20  -     -     -
BC Fjords        -      0.15  0.80  -     -     -     -     -
Salish Sea       -      0.02  -     0.80  0.05  -     -     -
WA Outer         -      0.10  -     0.02  0.40  0.35  -     -
Oregon           -      -     -     -     0.10  0.40  0.35  -
N. California    -      -     -     -     -     0.10  0.50  0.30
S. California    -      -     -     -     -     -     0.15  0.70
```

Key patterns:
- **High retention** in isolated regions (fjords 80%, Salish 80%)
- **Asymmetric southward flow** on outer coast
- **Very low connectivity** between inland seas and outer coast

---

## Temperature-Disease Relationship

From Eisenlord et al. 2016, Kohl et al. 2016:

```python
def disease_mortality_rate(base_rate: float, temperature: float) -> float:
    """Higher temperature increases SSWD mortality."""
    temp_effect = 1.0 + 0.1 * max(0, temperature - 10.0)
    return base_rate * temp_effect
```

This explains:
- Alaska (coldest) had lowest mortality (~60% vs 99% elsewhere)
- Fjords (cold deep water) retain large adults
- Southern range (warmest) was hit hardest

---

## Calibration Targets

To validate model against observed SSWD outbreak:

| Metric | Target Value | Source |
|--------|--------------|--------|
| Global decline | 90.6% | IUCN 2021 |
| WA/Salish decline | 99.2% | Harvell 2019 |
| Alaska decline | ~60% | Multiple sources |
| CA/OR survival | ~0% | Functional extinction |
| Range contraction | 50% | Southern half lost |

Disease spread timeline to match:
- June 2013: First cases (WA)
- 2013-2014: Peak spread
- 2014: Reaches Alaska
- 2015: Peak mortality
- Oregon lagged by ~1 year

---

## Implementation: PacificCoastConfig

```python
from dataclasses import dataclass, field
from typing import Dict, List, Tuple
from enum import Enum

class RegionType(Enum):
    OUTER_COAST = "outer_coast"
    INLAND_SEA = "inland_sea"
    FJORD = "fjord"

@dataclass
class RegionConfig:
    """Configuration for a single geographic region."""
    name: str
    latitude_range: Tuple[float, float]
    region_type: RegionType
    base_temperature: float  # Mean annual SST (°C)
    temperature_variance: float
    warming_rate: float  # °C per decade
    historical_density: float  # Pre-SSWD relative density (0-1)
    post_sswd_survival: float  # Observed survival fraction
    habitat_area_km2: float
    
@dataclass  
class PacificCoastConfig:
    """Geographic configuration for Pycnopodia network model."""
    
    regions: Dict[str, RegionConfig] = field(default_factory=lambda: {
        "se_alaska": RegionConfig(
            name="Southeast Alaska",
            latitude_range=(55.0, 60.0),
            region_type=RegionType.OUTER_COAST,
            base_temperature=7.0,
            temperature_variance=3.0,
            warming_rate=0.2,
            historical_density=0.7,
            post_sswd_survival=0.40,
            habitat_area_km2=5000.0
        ),
        "bc_outer": RegionConfig(
            name="BC Outer Coast",
            latitude_range=(48.5, 55.0),
            region_type=RegionType.OUTER_COAST,
            base_temperature=9.0,
            temperature_variance=3.0,
            warming_rate=0.3,
            historical_density=0.8,
            post_sswd_survival=0.15,
            habitat_area_km2=3000.0
        ),
        "bc_fjords": RegionConfig(
            name="BC Fjords (Refugia)",
            latitude_range=(51.0, 54.0),
            region_type=RegionType.FJORD,
            base_temperature=7.0,
            temperature_variance=2.0,
            warming_rate=0.2,
            historical_density=0.5,
            post_sswd_survival=0.50,
            habitat_area_km2=500.0
        ),
        "salish_sea": RegionConfig(
            name="Salish Sea / Puget Sound",
            latitude_range=(47.0, 49.5),
            region_type=RegionType.INLAND_SEA,
            base_temperature=11.0,
            temperature_variance=4.0,
            warming_rate=0.4,
            historical_density=1.0,
            post_sswd_survival=0.01,
            habitat_area_km2=2000.0
        ),
        "wa_outer": RegionConfig(
            name="Washington Outer Coast",
            latitude_range=(46.0, 48.5),
            region_type=RegionType.OUTER_COAST,
            base_temperature=10.0,
            temperature_variance=3.0,
            warming_rate=0.3,
            historical_density=0.7,
            post_sswd_survival=0.01,
            habitat_area_km2=800.0
        ),
        "or_coast": RegionConfig(
            name="Oregon Coast",
            latitude_range=(42.0, 46.0),
            region_type=RegionType.OUTER_COAST,
            base_temperature=11.0,
            temperature_variance=3.0,
            warming_rate=0.3,
            historical_density=0.5,
            post_sswd_survival=0.0,
            habitat_area_km2=1200.0
        ),
        "n_california": RegionConfig(
            name="Northern California",
            latitude_range=(37.0, 42.0),
            region_type=RegionType.OUTER_COAST,
            base_temperature=12.0,
            temperature_variance=3.0,
            warming_rate=0.3,
            historical_density=0.4,
            post_sswd_survival=0.0,
            habitat_area_km2=1000.0
        ),
        "s_california": RegionConfig(
            name="Southern California",
            latitude_range=(32.0, 37.0),
            region_type=RegionType.OUTER_COAST,
            base_temperature=15.0,
            temperature_variance=3.0,
            warming_rate=0.4,
            historical_density=0.2,
            post_sswd_survival=0.0,
            habitat_area_km2=800.0
        ),
    })
    
    # Larval dispersal parameters
    larval_duration_days: Tuple[int, int] = (14, 70)  # 2-10 weeks
    mean_larval_duration: int = 45
    dispersal_distance_per_day_km: float = 10.0
    
    # Disease parameters
    disease_base_mortality_rate: float = 0.15
    disease_temperature_coefficient: float = 0.1
    disease_transmission_rate: float = 0.2
    
    def get_temperature_mortality_modifier(self, temperature: float) -> float:
        """Calculate mortality rate modifier based on temperature."""
        return 1.0 + self.disease_temperature_coefficient * max(0, temperature - 10.0)
```

---

## Usage in Network Model

To generate a spatially realistic network:

```python
def create_pacific_coast_network(config: PacificCoastConfig, sites_per_region: int = 50):
    """Create network with real geographic structure."""
    
    network = PycnopodiaNetwork()
    
    for region_id, region in config.regions.items():
        # Create sites within region
        for i in range(sites_per_region):
            site = Site(
                region=region_id,
                latitude=random.uniform(*region.latitude_range),
                temperature=region.base_temperature + random.gauss(0, region.temperature_variance/3),
                initial_population=region.historical_density * random.uniform(0.5, 1.5)
            )
            network.add_site(site)
    
    # Add connectivity based on region pairs
    for from_region, connections in config.connectivity.items():
        for to_region, probability in connections.items():
            # Create edges between sites in connected regions
            ...
    
    return network
```

---

## Key References

| Paper | Finding | DOI |
|-------|---------|-----|
| Harvell et al. 2019 | 80-100% biomass decline, heat wave correlation | 10.1126/sciadv.aau7042 |
| Gravem et al. 2021 | IUCN Critically Endangered, 90.6% decline | IUCN database |
| Montecino-Latorre et al. 2016 | Salish Sea subtidal impacts | 10.1371/journal.pone.0163190 |
| Eisenlord et al. 2016 | Temperature increases SSWD mortality | 10.1098/rstb.2015.0212 |
| Gehman et al. 2025 | BC fjords as refugia | 10.1098/rspb.2024.2770 |
| Prentice et al. 2025 | Vibrio pectenicida causes SSWD | 10.1038/s41559-025-02797-2 |

---

## Next Steps

1. [ ] Implement `RegionConfig` and `PacificCoastConfig` in `pycnopodia/config.py`
2. [ ] Create network generator that uses this structure
3. [ ] Add temperature dynamics to disease model
4. [ ] Calibrate connectivity parameters to match genetic data
5. [ ] Validate against observed SSWD decline timeline
6. [ ] Test reintroduction scenarios from refugia (Alaska, BC fjords)
