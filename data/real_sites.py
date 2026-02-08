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
    
    # === SE ALASKA NORTH ADDITIONAL INLETS ===
    RealSite("Tracy Arm", 57.91, -133.40, "se_alaska_north", "fjord",
             sill_depth_m=40, base_temp_C=5.5, has_freshwater_lens=True,
             notes="30mi long, tidewater glaciers"),
    RealSite("Endicott Arm", 57.75, -133.30, "se_alaska_north", "fjord",
             sill_depth_m=45, base_temp_C=5.5, has_freshwater_lens=True,
             notes="Parallel to Tracy Arm, Dawes Glacier"),
    RealSite("Holkham Bay", 57.78, -133.62, "se_alaska_north", "fjord",
             sill_depth_m=80, base_temp_C=6.0, has_freshwater_lens=True,
             notes="Entrance to Tracy/Endicott Arms"),
    RealSite("Taku Inlet", 58.40, -133.90, "se_alaska_north", "fjord",
             sill_depth_m=100, base_temp_C=6.0, has_freshwater_lens=True,
             notes="South of Juneau, Taku Glacier"),
    RealSite("Port Snettisham", 58.10, -133.70, "se_alaska_north", "fjord",
             sill_depth_m=50, base_temp_C=6.5, has_freshwater_lens=True),
    RealSite("Port Frederick", 58.10, -135.00, "se_alaska_north", "fjord",
             sill_depth_m=30, base_temp_C=7.0, has_freshwater_lens=True,
             notes="Chichagof Island, Hoonah area"),
    RealSite("Hoonah Sound", 58.10, -135.40, "se_alaska_north", "inland_sea",
             base_temp_C=7.0),
    RealSite("Kelp Bay", 57.25, -134.95, "se_alaska_north", "fjord",
             sill_depth_m=25, base_temp_C=7.0),
    RealSite("Freshwater Bay", 57.95, -135.70, "se_alaska_north", "fjord",
             sill_depth_m=20, base_temp_C=7.0, has_freshwater_lens=True),
    
    # Prince William Sound — northernmost significant Pycnopodia habitat
    RealSite("Prince William Sound (west)", 60.60, -148.20, "se_alaska_north", "fjord",
             sill_depth_m=150, base_temp_C=5.0, has_freshwater_lens=True,
             notes="Whittier/Passage Canal area, glacial fjords"),
    RealSite("Prince William Sound (central)", 60.70, -147.00, "se_alaska_north", "inland_sea",
             base_temp_C=5.5, has_freshwater_lens=True,
             notes="Largest embayment in Gulf of Alaska"),
    RealSite("Prince William Sound (east)", 60.50, -146.00, "se_alaska_north", "fjord",
             sill_depth_m=80, base_temp_C=5.5, has_freshwater_lens=True,
             notes="Cordova area, copper river delta influence"),
    RealSite("Resurrection Bay", 60.00, -149.40, "se_alaska_north", "fjord",
             sill_depth_m=100, base_temp_C=5.5, has_freshwater_lens=True,
             notes="Seward, Alaska SeaLife Center"),
    RealSite("Kachemak Bay", 59.60, -151.30, "se_alaska_north", "inland_sea",
             base_temp_C=6.0,
             notes="Homer area, KBNERR, strong tidal mixing"),
    RealSite("Yakutat Bay", 59.50, -139.80, "se_alaska_north", "fjord",
             sill_depth_m=60, base_temp_C=5.5, has_freshwater_lens=True,
             notes="Bridges SE AK to Gulf of Alaska"),
    
    # === Gulf of Alaska coast (bridging SE AK to PWS) ===
    RealSite("Lituya Bay", 58.65, -137.50, "se_alaska_north", "fjord",
             sill_depth_m=10, base_temp_C=5.5, has_freshwater_lens=True,
             notes="Extremely shallow sill, famous for megatsunamis"),
    RealSite("Icy Bay", 59.95, -141.40, "se_alaska_north", "fjord",
             sill_depth_m=50, base_temp_C=5.0, has_freshwater_lens=True,
             notes="Guyot/Tyndall glaciers, recently deglaciated"),
    RealSite("Controller Bay", 60.15, -144.20, "se_alaska_north", "outer_coast",
             base_temp_C=5.5,
             notes="Between Yakutat and Cordova, Bering Glacier outflow"),
    RealSite("Cordova / Orca Inlet", 60.55, -145.75, "se_alaska_north", "inland_sea",
             base_temp_C=5.5, has_freshwater_lens=True,
             notes="Copper River delta, links to PWS east"),
    
    # === Kenai Peninsula / Cook Inlet ===
    RealSite("Seward (outer coast)", 59.90, -149.50, "se_alaska_north", "outer_coast",
             base_temp_C=6.0,
             notes="Gulf of Alaska side of Kenai Peninsula"),
    RealSite("Aialik Bay", 59.75, -149.65, "se_alaska_north", "fjord",
             sill_depth_m=80, base_temp_C=5.5, has_freshwater_lens=True,
             notes="Kenai Fjords NP, tidewater glaciers"),
    RealSite("Kenai outer coast", 59.50, -150.50, "se_alaska_north", "outer_coast",
             base_temp_C=6.0),
    RealSite("Lower Cook Inlet", 59.20, -152.00, "se_alaska_north", "inland_sea",
             base_temp_C=6.5,
             notes="Between Kenai and Kodiak, strong tidal mixing"),
    
    # === Kodiak Island ===
    RealSite("Kodiak Island (east)", 57.80, -152.40, "se_alaska_north", "island",
             base_temp_C=6.5,
             notes="Large island, significant Pycnopodia habitat"),
    RealSite("Kodiak Island (west)", 57.50, -154.00, "se_alaska_north", "island",
             base_temp_C=6.0,
             notes="Shelikof Strait side"),
    
    # === Aleutian Islands (western range limit) ===
    RealSite("Unalaska / Dutch Harbor", 53.88, -166.53, "se_alaska_north", "island",
             base_temp_C=5.0,
             notes="Eastern Aleutians, confirmed Pycnopodia presence"),
    RealSite("Akutan Island", 54.13, -165.77, "se_alaska_north", "island",
             base_temp_C=5.0),
    RealSite("Sanak Islands", 54.42, -162.77, "se_alaska_north", "island",
             base_temp_C=5.5,
             notes="Between Kodiak and Unalaska"),
    RealSite("Sand Point / Shumagin Islands", 55.33, -160.50, "se_alaska_north", "island",
             base_temp_C=5.5,
             notes="Alaska Peninsula south side"),
    RealSite("Chignik Bay", 56.30, -158.40, "se_alaska_north", "inland_sea",
             base_temp_C=5.5,
             notes="Alaska Peninsula, sheltered bay"),
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
    RealSite("LeConte Bay", 56.80, -132.35, "se_alaska_south", "fjord",
             sill_depth_m=25, base_temp_C=6.5, has_freshwater_lens=True,
             notes="Southernmost tidewater glacier in N. America"),
    RealSite("Thomas Bay", 57.00, -132.60, "se_alaska_south", "fjord",
             sill_depth_m=35, base_temp_C=6.5, has_freshwater_lens=True,
             notes="North of Petersburg, Baird Glacier"),
    RealSite("Farragut Bay", 57.08, -132.75, "se_alaska_south", "fjord",
             sill_depth_m=30, base_temp_C=7.0),
    RealSite("Eastern Passage", 56.50, -132.10, "se_alaska_south", "fjord",
             sill_depth_m=80, base_temp_C=7.0, has_freshwater_lens=True,
             notes="Stikine River delta, major freshwater input"),
    RealSite("Blake Channel", 56.20, -132.70, "se_alaska_south", "fjord",
             sill_depth_m=60, base_temp_C=7.5),
    RealSite("Duncan Canal", 56.75, -133.00, "se_alaska_south", "fjord",
             sill_depth_m=40, base_temp_C=7.0),
    RealSite("Zimovia Strait", 56.35, -132.25, "se_alaska_south", "fjord",
             sill_depth_m=50, base_temp_C=7.5),
    RealSite("Sumner Strait", 56.30, -133.00, "se_alaska_south", "inland_sea",
             base_temp_C=7.5),
    RealSite("Ernest Sound", 55.85, -132.30, "se_alaska_south", "fjord",
             sill_depth_m=100, base_temp_C=7.5),
    RealSite("Bradfield Canal", 56.25, -131.50, "se_alaska_south", "fjord",
             sill_depth_m=45, base_temp_C=7.0, has_freshwater_lens=True,
             notes="Near BC border, glacial"),
    RealSite("Rudyerd Bay", 55.55, -130.85, "se_alaska_south", "fjord",
             sill_depth_m=30, base_temp_C=7.0, has_freshwater_lens=True,
             notes="Misty Fjords National Monument"),
    RealSite("Walker Cove", 55.45, -131.55, "se_alaska_south", "fjord",
             sill_depth_m=35, base_temp_C=7.5),
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
    
    # === BC CENTRAL COAST FJORDS (50-52°N) - major inlet systems ===
    RealSite("Knight Inlet", 50.80, -125.70, "bc_fjords", "fjord",
             sill_depth_m=70, base_temp_C=7.0, has_freshwater_lens=True,
             notes="125km long, one of longest BC fjords, major glacial input"),
    RealSite("Bute Inlet", 50.65, -124.89, "bc_fjords", "fjord",
             sill_depth_m=100, base_temp_C=7.0, has_freshwater_lens=True,
             notes="80km long, Homathko River glacial input, DFO monitoring"),
    RealSite("Toba Inlet", 50.25, -124.60, "bc_fjords", "fjord",
             sill_depth_m=120, base_temp_C=7.5, has_freshwater_lens=True,
             notes="35km long, >480m deep, jade-green glacial waters"),
    RealSite("Loughborough Inlet", 50.52, -125.53, "bc_fjords", "fjord",
             sill_depth_m=50, base_temp_C=7.5, has_freshwater_lens=True,
             notes="35km long, Discovery Islands region"),
    RealSite("Kingcome Inlet", 50.92, -126.50, "bc_fjords", "fjord",
             sill_depth_m=45, base_temp_C=7.5, has_freshwater_lens=True,
             notes="35km long, Kingcome River delta, First Nations territory"),
    RealSite("Seymour Inlet", 51.05, -127.25, "bc_fjords", "fjord",
             sill_depth_m=30, base_temp_C=7.5, has_freshwater_lens=True,
             notes="58km long, Nakwakto Rapids entrance"),
    RealSite("Wakeman Sound", 50.75, -126.45, "bc_fjords", "fjord",
             sill_depth_m=40, base_temp_C=7.5, has_freshwater_lens=True,
             notes="Tributary of Kingcome, grizzly habitat"),
    RealSite("Tribune Channel", 50.70, -126.60, "bc_fjords", "fjord",
             sill_depth_m=60, base_temp_C=7.5,
             notes="Between Gilford/Baker Islands"),
    RealSite("Frederick Sound (BC)", 51.04, -126.72, "bc_fjords", "fjord",
             sill_depth_m=35, base_temp_C=7.5, has_freshwater_lens=True),
    RealSite("Allison Sound", 51.05, -127.15, "bc_fjords", "fjord",
             sill_depth_m=25, base_temp_C=7.5),
    RealSite("Cascade Inlet", 52.50, -127.52, "bc_fjords", "fjord",
             sill_depth_m=35, base_temp_C=7.0, has_freshwater_lens=True),
    RealSite("Cousins Inlet", 52.32, -127.75, "bc_fjords", "fjord",
             sill_depth_m=30, base_temp_C=7.0),
    RealSite("Fitz Hugh Sound", 51.68, -127.92, "bc_fjords", "inland_sea",
             base_temp_C=8.0, notes="Mouth of Dean Channel system"),
    RealSite("Kwatna Inlet", 52.45, -127.65, "bc_fjords", "fjord",
             sill_depth_m=40, base_temp_C=7.0, has_freshwater_lens=True),
    RealSite("South Bentinck Arm", 52.30, -126.85, "bc_fjords", "fjord",
             sill_depth_m=20, base_temp_C=7.0, has_freshwater_lens=True),
    
    # === BC NORTH COAST FJORDS (53-55°N) ===
    RealSite("Grenville Channel", 53.40, -129.50, "bc_fjords", "fjord",
             sill_depth_m=200, base_temp_C=7.0,
             notes="83km long, 0.2nm wide at narrowest, Inside Passage"),
    RealSite("Devastation Channel", 53.67, -128.84, "bc_fjords", "fjord",
             sill_depth_m=80, base_temp_C=7.0),
    RealSite("Observatory Inlet", 55.33, -129.90, "bc_fjords", "fjord",
             sill_depth_m=60, base_temp_C=7.0, has_freshwater_lens=True),
    RealSite("Alice Arm", 55.45, -129.59, "bc_fjords", "fjord",
             sill_depth_m=40, base_temp_C=7.0, has_freshwater_lens=True,
             notes="East arm of Observatory Inlet"),
    RealSite("Hastings Arm", 55.50, -129.77, "bc_fjords", "fjord",
             sill_depth_m=50, base_temp_C=7.0, has_freshwater_lens=True),
    RealSite("Work Channel", 53.50, -128.70, "bc_fjords", "fjord",
             sill_depth_m=90, base_temp_C=7.0),
    RealSite("Ursula Channel", 53.55, -129.00, "bc_fjords", "fjord",
             sill_depth_m=70, base_temp_C=7.0),
    RealSite("Principe Channel", 53.80, -130.10, "bc_fjords", "fjord",
             sill_depth_m=150, base_temp_C=7.5),
    RealSite("Khutzeymateen Inlet", 54.50, -130.00, "bc_fjords", "fjord",
             sill_depth_m=35, base_temp_C=7.0, has_freshwater_lens=True,
             notes="Protected grizzly bear sanctuary"),
    RealSite("Europe Reach", 53.44, -128.42, "bc_fjords", "fjord",
             sill_depth_m=50, base_temp_C=7.0, notes="Gardner Canal complex"),
]

# ============================================================
# BC OUTER COAST (48-55°N, exposed Pacific including Vancouver Island)
# ============================================================
BC_OUTER = [
    # Vancouver Island outer coast — major exposed Pacific coastline
    RealSite("Barkley Sound", 48.85, -125.30, "bc_outer", "inland_sea",
             base_temp_C=9.0,
             notes="Bamfield Marine Sciences Centre, Pacific Rim NP area"),
    RealSite("Clayoquot Sound", 49.15, -125.90, "bc_outer", "inland_sea",
             base_temp_C=9.0,
             notes="Major inlet system on west Vancouver Island"),
    RealSite("Nootka Sound", 49.60, -126.60, "bc_outer", "inland_sea",
             base_temp_C=8.5,
             notes="Central west coast Vancouver Island"),
    RealSite("Kyuquot Sound", 50.00, -127.20, "bc_outer", "outer_coast",
             base_temp_C=8.5,
             notes="NW Vancouver Island"),
    RealSite("Quatsino Sound", 50.50, -128.05, "bc_outer", "inland_sea",
             base_temp_C=8.5,
             notes="Northern tip Vancouver Island, connects to Queen Charlotte Strait"),
    RealSite("Brooks Peninsula", 50.15, -127.80, "bc_outer", "outer_coast",
             base_temp_C=8.5,
             notes="Major biogeographic boundary, exposed headland"),
    
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
    RealSite("Saratoga Passage", 48.20, -122.55, "salish_sea", "inland_sea",
             base_temp_C=10.0,
             notes="Between Whidbey and Camano Islands, important Pycnopodia habitat"),
    RealSite("Admiralty Inlet", 48.15, -122.75, "salish_sea", "inland_sea",
             base_temp_C=9.5,
             notes="Entrance to Puget Sound from Strait of Juan de Fuca"),
    RealSite("Bellingham Bay", 48.75, -122.50, "salish_sea", "inland_sea",
             base_temp_C=9.5,
             notes="Northern Puget Sound, near BC border"),
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
    
    # === SALISH SEA FJORD ARMS ===
    RealSite("Indian Arm", 49.38, -122.88, "salish_sea", "fjord",
             sill_depth_m=60, base_temp_C=9.0, has_freshwater_lens=True,
             notes="20km long, extends north from Vancouver"),
    RealSite("Sechelt Inlet", 49.73, -123.73, "salish_sea", "fjord",
             sill_depth_m=20, base_temp_C=9.0, has_freshwater_lens=True,
             notes="40km long, Skookumchuck Rapids entrance"),
    RealSite("Salmon Inlet", 49.70, -123.55, "salish_sea", "fjord",
             sill_depth_m=15, base_temp_C=9.0),
    RealSite("Saanich Inlet", 48.62, -123.50, "salish_sea", "fjord",
             sill_depth_m=70, base_temp_C=9.0,
             notes="25km long, Vancouver Island, 60+ year DFO time series"),
    RealSite("Finlayson Arm", 48.57, -123.53, "salish_sea", "fjord",
             sill_depth_m=50, base_temp_C=9.0),
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
    RealSite("Coronado Islands", 32.42, -117.25, "s_california", "island",
             base_temp_C=15.5,
             notes="US/Mexico border area"),
    
    # Baja California — southern range limit
    RealSite("Ensenada", 31.85, -116.62, "s_california", "outer_coast",
             base_temp_C=16.0,
             notes="Northern Baja, historically had Pycnopodia"),
    RealSite("Isla Guadalupe", 29.05, -118.28, "s_california", "island",
             base_temp_C=16.0,
             notes="Offshore Baja island, southernmost confirmed Pycnopodia records"),
    RealSite("Punta Banda", 31.72, -116.72, "s_california", "outer_coast",
             base_temp_C=16.0,
             notes="South of Ensenada, kelp forests"),
    RealSite("San Quintín", 30.50, -115.95, "s_california", "outer_coast",
             base_temp_C=16.5,
             notes="Baja California, upwelling zone, historical range"),
    RealSite("Isla Cedros", 28.10, -115.20, "s_california", "island",
             base_temp_C=17.0,
             notes="Major upwelling zone, abalone/kelp habitat"),
    RealSite("Punta Eugenia", 27.85, -115.08, "s_california", "outer_coast",
             base_temp_C=17.0,
             notes="Biogeographic boundary, southern range limit area"),
    RealSite("Bahía Tortugas", 27.68, -114.90, "s_california", "outer_coast",
             base_temp_C=17.5,
             notes="Near southern historical Pycnopodia limit"),
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
