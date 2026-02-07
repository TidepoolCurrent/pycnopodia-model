"""
Pycnopodia helianthoides population model.

Individual-based simulation integrating:
- Sweepstakes reproductive success (Hedgecock & Pudovkin 2011)
- Polygenic SSWD resistance (Schiebelhut et al. 2024)
- Broadcast spawner Allee effects
- Captive breeding and outplanting interventions

Reference: Weertman 2026
"""

from .population import Population
from .simulation import Simulation
from .config import Config, DEFAULT_CONFIG

__version__ = "0.1.0"
__all__ = ["Population", "Simulation", "Config", "DEFAULT_CONFIG"]
