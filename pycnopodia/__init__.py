"""
Pycnopodia helianthoides population model.

Individual-based simulation integrating:
- Sweepstakes reproductive success (Hedgecock & Pudovkin 2011)
- Polygenic SSWD resistance (Schiebelhut et al. 2024)
- Broadcast spawner Allee effects
- Captive breeding and outplanting interventions
- Environmental stochasticity and climate projections
- Spatial structure and local interactions

Reference: Weertman 2026
"""

from .population import Population
from .simulation import Simulation
from .config import Config, DEFAULT_CONFIG
from .environment import EnvironmentConfig, EnvironmentState
from .spatial import SpatialConfig, SpatialPopulation
from .srs import HierarchicalSRSConfig, HierarchicalSRS
from .size_structured import SizeStructuredConfig, SizeStructuredPopulation
from .real_data import (
    HISTORICAL_SST, 
    get_historical_sst,
    get_sswd_parameters,
    generate_temperature_scenario,
    get_larval_connectivity,
)
from .pacific_coast import (
    PacificCoastConfig,
    PacificCoastSimulation,
    PacificCoastResult,
    PACIFIC_COAST_REGIONS,
    REGION_ORDER,
    RegionType,
    RegionConfig,
)

__version__ = "0.4.0"
__all__ = [
    "Population", 
    "Simulation", 
    "Config", 
    "DEFAULT_CONFIG",
    "EnvironmentConfig",
    "EnvironmentState", 
    "SpatialConfig",
    "SpatialPopulation",
    "HierarchicalSRSConfig",
    "HierarchicalSRS",
    "PacificCoastConfig",
    "PacificCoastSimulation",
    "PacificCoastResult",
    "PACIFIC_COAST_REGIONS",
    "REGION_ORDER",
    "RegionType",
]
