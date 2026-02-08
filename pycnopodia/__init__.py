"""
Pycnopodia helianthoides population model.

Site-based metapopulation simulation with:
- Geographic realism (56 named Pacific Coast sites)
- Fjord refugia with sill-depth and freshwater lens effects
- Polygenic SSWD resistance evolution
- Broadcast spawner Allee effects
- Distance-based larval and disease connectivity
- Marine heatwave (The Blob) disease trigger

Reference: Weertman 2026
"""

from .config import Config, DEFAULT_CONFIG
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
from .geo_model import (
    GeoConfig,
    GeoSimulation,
    GeoResult,
    GeoState,
)

__version__ = "0.5.0"
__all__ = [
    "Config", 
    "DEFAULT_CONFIG",
    "PacificCoastConfig",
    "PacificCoastSimulation",
    "PacificCoastResult",
    "PACIFIC_COAST_REGIONS",
    "REGION_ORDER",
    "RegionType",
    "GeoConfig",
    "GeoSimulation",
    "GeoResult",
    "GeoState",
]
