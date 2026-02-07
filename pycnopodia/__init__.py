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

__version__ = "0.2.0"
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
]
