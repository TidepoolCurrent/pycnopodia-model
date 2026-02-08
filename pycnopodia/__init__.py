"""
Pycnopodia helianthoides population model.

Geography-based metapopulation simulation with:
- 142 real Pacific Coast sites (29°N-61°N)
- 75 fjords with sill-depth and freshwater lens effects
- Geographic disease spread from epicenter
- Seasonal dynamics (quarterly time steps, winter spawning)
- Polygenic SSWD resistance evolution (10 loci)
- Broadcast spawner Allee effects
- Distance-based larval and disease connectivity
- Marine heatwave (The Blob) disease trigger

Reference: Weertman 2026
"""

from .geo_model import (
    GeoConfig,
    GeoSimulation,
    GeoResult,
    GeoState,
    run_geo_ensemble,
)

__version__ = "0.6.0"
__all__ = [
    "GeoConfig",
    "GeoSimulation",
    "GeoResult",
    "GeoState",
    "run_geo_ensemble",
]
