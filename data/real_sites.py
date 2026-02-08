"""
Real geographic sites for Pycnopodia model.

Sites derived from literature survey locations, known fjord refugia,
and monitoring station coordinates. Each site has:
- name: human-readable location
- lat, lon: decimal degrees
- region: model region ID
- site_type: 'fjord', 'outer_coast', 'inland_sea', 'island'
- sill_depth_m: approximate sill depth (fjords only, controls exchange)
- base_temp_C: approximate annual mean bottom temperature
- has_freshwater_lens: whether glacial/snow melt creates surface freshwater layer
- notes: relevant biology/survey info

Sources:
- Hamilton et al. 2021 (survey sites)
- Gehman et al. 2025 (BC fjord refugia)
- Konar et al. 2019 (Alaska intertidal)
- REEF/MARINe monitoring networks
- NOAA Status Review Appendix A (2022)
"""

from dataclasses import dataclass
from typing import Optional, List
import math


@dataclass
class RealSite:
    name: str
    lat: float
    lon: float
    region: str
    site_type: str  # fjord, outer_coast, inland_sea, island
    sill_depth_m: Optional[float] = None
    base_temp_C: float = 8.0
    has_freshwater_lens: bool = False
    notes: str = ""


# ============================================================
# SE ALASKA NORTH (fjord-dominated, >57°N)
# ============================================================
SE_ALASKA_NORTH = [
    # Fjord sites — confirmed Pycnopodia habitat
    RealSite("Glacier Bay (inner)", 58.75, -136.10, "se_alaska_north", "fjord",
             sill_depth_m=20, base_temp_C=5.5, has_freshwater_lens=True,
             notes="National Park, glacial fjord, very cold"),
    RealSite("Glacier Bay (outer)", 58.50, -136.00, "se_alaska_north", "fjord",
             sill_depth_m=60, base_temp_C=6.5, has_freshwater_lens=True),
    RealSite("Lynn Canal", 58.80, -135.30, "se_alaska_north", "fjord",
             sill_depth_m=300, base_temp_C=6.0, has_freshwater_lens=True,
             notes="Deepest fjord in NA, connects to Juneau"),
    RealSite("Icy Strait", 58.30, -135.80, "se_alaska_north", "outer_coast",
             base_temp_C=7.0, notes="Major current passage"),
    RealSite("Chatham Strait (north)", 57.80, -134.80, "se_alaska_north", "fjord",
             sill_depth_m=200, base_temp_C=6.5),
    RealSite("Tenakee Inlet", 57.78, -135.22, "se_alaska_north", "fjord",
             sill_depth_m=30, base_temp_C=6.5, has_freshwater_lens=True,
             notes="Chichagof Island, isolated fjord arm"),
    RealSite("Peril Strait", 57.45, -135.55, "se_alaska_north", "fjord",
             sill_depth_m=15, base_temp_C=7.0,
             notes="Narrow, tidal rapids limit exchange"),
    RealSite("Juneau area", 58.30, -134.42, "se_alaska_north", "fjord",
             sill_depth_m=100, base_temp_C=6.5, has_freshwater_lens=True),
]

# ============================================================
# SE ALASKA SOUTH (outer coast dominated, 55-57°N)
# ============================================================
SE_ALASKA_SOUTH = [
    RealSite("Sitka Sound", 57.05, -135.33, "se_alaska_south", "island",
             base_temp_C=7.5, notes="Major survey site, exposed"),
    RealSite("Baranof Island (south)", 56.30, -134.80, "se_alaska_south", "island",
             base_temp_C=7.0,
             notes="Bridges PoW west to Sitka Sound"),
    RealSite("Ketchikan area", 55.35, -131.64, "se_alaska_south", "island",
             base_temp_C=8.0, notes="Southernmost SE AK"),
    RealSite("Prince of Wales Island (west)", 55.90, -133.50, "se_alaska_south", "outer_coast",
             base_temp_C=7.5),
    RealSite("Clarence Strait", 55.60, -132.00, "se_alaska_south", "fjord",
             sill_depth_m=150, base_temp_C=7.5),
    RealSite("Behm Canal", 55.50, -131.30, "se_alaska_south", "fjord",
             sill_depth_m=50, base_temp_C=7.5, has_freshwater_lens=True),
    RealSite("Misty Fjords", 55.60, -130.60, "se_alaska_south", "fjord",
             sill_depth_m=30, base_temp_C=7.0, has_freshwater_lens=True,
             notes="National Monument, glacial"),
]

# ============================================================
# BC FJORDS (Central + North Coast, 51-54°N)
# ============================================================
BC_FJORDS = [
    # Confirmed refugia (Gehman et al. 2025)
    RealSite("Burke Channel", 52.15, -127.40, "bc_fjords", "fjord",
             sill_depth_m=40, base_temp_C=7.0, has_freshwater_lens=True,
             notes="CONFIRMED REFUGIUM - Gehman et al. 2025"),
    RealSite("Dean Channel", 52.35, -127.10, "bc_fjords", "fjord",
             sill_depth_m=50, base_temp_C=6.5, has_freshwater_lens=True,
             notes="170km total length, glacial input"),
    RealSite("Fisher Channel", 52.00, -127.80, "bc_fjords", "fjord",
             sill_depth_m=60, base_temp_C=7.0),
    RealSite("Bella Coola / North Bentinck Arm", 52.37, -126.75, "bc_fjords", "fjord",
             sill_depth_m=20, base_temp_C=7.0, has_freshwater_lens=True,
             notes="Head of Dean Channel system, glacial river input"),
    
    # Other Central Coast fjords
    RealSite("Rivers Inlet", 51.65, -127.25, "bc_fjords", "fjord",
             sill_depth_m=35, base_temp_C=7.5, has_freshwater_lens=True),
    RealSite("Smith Inlet", 51.40, -127.50, "bc_fjords", "fjord",
             sill_depth_m=25, base_temp_C=7.5, has_freshwater_lens=True),
    RealSite("Belize Inlet", 51.10, -127.10, "bc_fjords", "fjord",
             sill_depth_m=20, base_temp_C=8.0),
    
    # Southern BC / Discovery Islands transition zone
    # These sites bridge the Salish Sea and BC Central Coast
    RealSite("Discovery Islands", 50.25, -125.30, "bc_fjords", "island",
             base_temp_C=8.5,
             notes="Campbell River area, gateway between Salish Sea and open BC coast"),
    RealSite("Johnstone Strait", 50.50, -126.30, "bc_fjords", "inland_sea",
             base_temp_C=8.5,
             notes="Major tidal channel connecting Georgia Strait to Queen Charlotte Strait"),
    RealSite("Queen Charlotte Strait", 50.80, -127.00, "bc_fjords", "inland_sea",
             base_temp_C=8.5,
             notes="Northern end of Inside Passage, connects to open Pacific"),
    
    # North Coast fjords
    RealSite("Douglas Channel", 53.50, -129.20, "bc_fjords", "fjord",
             sill_depth_m=100, base_temp_C=7.0, has_freshwater_lens=True,
             notes="Long deep fjord, Kitimat terminal"),
    RealSite("Gardner Canal", 53.30, -128.60, "bc_fjords", "fjord",
             sill_depth_m=40, base_temp_C=6.5, has_freshwater_lens=True,
             notes="Remote glacial fjord"),
    RealSite("Portland Canal", 55.00, -130.00, "bc_fjords", "fjord",
             sill_depth_m=30, base_temp_C=7.0, has_freshwater_lens=True,
             notes="BC/AK border"),
    
    # Outer coast / islands
    RealSite("Calvert Island (Hakai)", 51.65, -128.13, "bc_fjords", "island",
             base_temp_C=8.5, notes="Hakai Institute field station"),
    RealSite("Central Coast outer islands", 52.00, -128.50, "bc_fjords", "outer_coast",
             base_temp_C=8.5),
]

# ============================================================
# BC OUTER COAST (52-55°N, exposed)
# ============================================================
BC_OUTER = [
    RealSite("Haida Gwaii (east)", 53.00, -131.80, "bc_outer", "island",
             base_temp_C=8.5, notes="Previously known as Queen Charlotte Is."),
    RealSite("Haida Gwaii (west)", 53.20, -132.80, "bc_outer", "outer_coast",
             base_temp_C=8.0, notes="Fully exposed Pacific"),
    RealSite("Hecate Strait", 53.00, -130.50, "bc_outer", "outer_coast",
             base_temp_C=8.5),
    RealSite("Princess Royal Island", 52.70, -128.70, "bc_outer", "island",
             base_temp_C=8.0),
    RealSite("Milbanke Sound", 52.15, -128.90, "bc_outer", "outer_coast",
             base_temp_C=8.5,
             notes="Bridges BC Fjords central coast to outer BC"),
    RealSite("Prince Rupert", 54.30, -130.30, "bc_outer", "outer_coast",
             base_temp_C=8.0,
             notes="Major port, bridges Douglas Channel to Portland Canal"),
]

# ============================================================
# SALISH SEA (48-49.5°N)
# ============================================================
SALISH_SEA = [
    RealSite("San Juan Islands", 48.55, -123.00, "salish_sea", "island",
             base_temp_C=9.5, notes="Friday Harbor Labs — human's institution!"),
    RealSite("Howe Sound", 49.40, -123.30, "salish_sea", "fjord",
             sill_depth_m=70, base_temp_C=9.0,
             notes="Fjord within Salish Sea, possible micro-refugium"),
    RealSite("Jervis Inlet", 49.80, -123.80, "salish_sea", "fjord",
             sill_depth_m=200, base_temp_C=8.5, has_freshwater_lens=True),
    RealSite("Desolation Sound", 50.10, -124.70, "salish_sea", "inland_sea",
             base_temp_C=9.0),
    RealSite("Georgia Strait (central)", 49.30, -123.80, "salish_sea", "inland_sea",
             base_temp_C=9.5),
    RealSite("Puget Sound (central)", 47.60, -122.50, "salish_sea", "inland_sea",
             base_temp_C=10.0),
    RealSite("Puget Sound (south)", 47.20, -122.60, "salish_sea", "inland_sea",
             base_temp_C=10.5),
    RealSite("Hood Canal", 47.70, -123.00, "salish_sea", "fjord",
             sill_depth_m=50, base_temp_C=10.0,
             notes="Fjord in Puget Sound, periodic hypoxia"),
]

# ============================================================
# WA/OR OUTER COAST (44-48°N)
# ============================================================
WA_OR = [
    RealSite("Olympic Coast", 47.80, -124.60, "wa_or_outer", "outer_coast",
             base_temp_C=9.5),
    RealSite("Grays Harbor area", 46.90, -124.20, "wa_or_outer", "outer_coast",
             base_temp_C=10.0),
    RealSite("Lincoln City", 44.95, -124.00, "wa_or_outer", "outer_coast",
             base_temp_C=10.0),
    RealSite("Tillamook", 45.50, -123.95, "wa_or_outer", "outer_coast",
             base_temp_C=10.0),
    RealSite("Astoria / Columbia River", 46.20, -124.00, "wa_or_outer", "outer_coast",
             base_temp_C=10.0,
             notes="Columbia River mouth, major freshwater influence"),
    RealSite("Oregon Central Coast", 44.60, -124.10, "wa_or_outer", "outer_coast",
             base_temp_C=10.0),
    RealSite("Coos Bay", 43.40, -124.30, "wa_or_outer", "outer_coast",
             base_temp_C=10.5),
    RealSite("Cape Blanco", 42.80, -124.55, "wa_or_outer", "outer_coast",
             base_temp_C=10.5,
             notes="Westernmost point in Oregon"),
    RealSite("Strait of Juan de Fuca", 48.20, -123.60, "wa_or_outer", "outer_coast",
             base_temp_C=9.5,
             notes="Connects Salish Sea to open Pacific, high tidal exchange"),
]

# ============================================================
# N. CALIFORNIA (39-42°N)
# ============================================================
N_CALIFORNIA = [
    RealSite("Brookings", 42.05, -124.30, "n_california", "outer_coast",
             base_temp_C=10.5,
             notes="OR/CA border, bridges WA/OR to N. California"),
    RealSite("Crescent City", 41.75, -124.20, "n_california", "outer_coast",
             base_temp_C=10.5),
    RealSite("Trinidad Head", 41.06, -124.15, "n_california", "outer_coast",
             base_temp_C=10.5),
    RealSite("Cape Mendocino", 40.40, -124.40, "n_california", "outer_coast",
             base_temp_C=10.5,
             notes="Bridges Trinidad Head to Mendocino"),
    RealSite("Mendocino", 39.30, -123.80, "n_california", "outer_coast",
             base_temp_C=11.0),
    RealSite("Point Arena", 38.95, -123.74, "n_california", "outer_coast",
             base_temp_C=11.0),
]

# ============================================================
# C. CALIFORNIA (35-39°N)
# ============================================================
C_CALIFORNIA = [
    RealSite("Point Reyes", 38.05, -123.02, "c_california", "outer_coast",
             base_temp_C=11.5),
    RealSite("Santa Cruz", 36.95, -122.00, "c_california", "outer_coast",
             base_temp_C=12.0,
             notes="Bridges Monterey Bay to SF area"),
    RealSite("Half Moon Bay", 37.46, -122.45, "c_california", "outer_coast",
             base_temp_C=12.0),
    RealSite("Farallon Islands", 37.70, -123.00, "c_california", "island",
             base_temp_C=11.5),
    RealSite("Monterey Bay", 36.80, -121.90, "c_california", "outer_coast",
             base_temp_C=11.0, notes="Upwelling zone — cold but exposed"),
    RealSite("Big Sur", 36.20, -121.80, "c_california", "outer_coast",
             base_temp_C=12.0),
    RealSite("San Luis Obispo", 35.18, -120.75, "c_california", "outer_coast",
             base_temp_C=12.5,
             notes="Bridges Santa Barbara Channel to Monterey"),
    RealSite("Morro Bay", 35.37, -120.86, "c_california", "outer_coast",
             base_temp_C=12.5),
]

# ============================================================
# S. CALIFORNIA (32-35°N)
# ============================================================
S_CALIFORNIA = [
    RealSite("Santa Barbara Channel", 34.40, -119.85, "s_california", "outer_coast",
             base_temp_C=13.5),
    RealSite("Channel Islands", 34.00, -119.50, "s_california", "island",
             base_temp_C=13.5),
    RealSite("Palos Verdes", 33.74, -118.40, "s_california", "outer_coast",
             base_temp_C=14.5),
    RealSite("La Jolla", 32.85, -117.27, "s_california", "outer_coast",
             base_temp_C=15.0),
]

# ============================================================
# All sites
# ============================================================
ALL_SITES = (
    SE_ALASKA_NORTH + SE_ALASKA_SOUTH + BC_FJORDS + BC_OUTER +
    SALISH_SEA + WA_OR + N_CALIFORNIA + C_CALIFORNIA + S_CALIFORNIA
)


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Great-circle distance between two points."""
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat/2)**2 +
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) *
         math.sin(dlon/2)**2)
    return R * 2 * math.asin(math.sqrt(a))


if __name__ == "__main__":
    print(f"Total sites: {len(ALL_SITES)}")
    for region in ["se_alaska_north", "se_alaska_south", "bc_fjords", "bc_outer",
                    "salish_sea", "wa_or_outer", "n_california", "c_california", "s_california"]:
        sites = [s for s in ALL_SITES if s.region == region]
        fjords = [s for s in sites if s.site_type == "fjord"]
        fw = [s for s in sites if s.has_freshwater_lens]
        print(f"  {region:20s}: {len(sites):2d} sites ({len(fjords)} fjord, {len(fw)} freshwater lens)")
