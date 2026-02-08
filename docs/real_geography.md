# Real Geography for Pycnopodia Model

## Fjord Counts
- **BC Coast: ~150 fjords** (locally called "inlets") — one of the major fjord coastlines in the world
  - Source: Geoscience World, 2010
  - BC coastline: 29,000 km (including islands)
  - Notable Pycnopodia refugia: **Burke Channel** (Gehman et al. 2025, Proc R Soc B)
  - Central Coast fjords specifically identified as sunflower star refugia by Hakai Institute
  - Freshwater lens from glacial/snow melt forces stars deeper into cold water (unique mechanism)

- **SE Alaska: thousands of narrow inlets**
  - Alaska total: 34,000 miles of indented tidal coastline, 15,000 sq mi of fjords
  - SE Alaska panhandle is heavily fjorded
  - Kenai Fjords had Pycnopodia density ~0.075/m² (Konar et al. 2019)
  - Prince William Sound: 0.233/m² (pre-SSWD)

## Key Paper: Gehman et al. 2025
"Threatened sea stars are finding refuge in BC fjords"
DOI: 10.1098/rspb.2024.2770

Key findings:
1. SSWD DID reach the fjords — they are NOT disease-free
2. Stars survived better due to:
   - Freshwater surface lens from glacial/snow melt (low salinity lid)
   - Forces animals deeper into cold water
   - Cold water reduces disease severity
3. Burke Channel specifically named as refuge
4. "Different way for the environment to be protecting an animal"
5. Warming/reduced snowfall threatens this refuge

## Model Implications
- Current abstract model has 40 sites per fjord region
- Real geography: ~150 fjords in BC, many more in SE AK
- Scaling up to real coastline would mean:
  - Each fjord = a discrete unit with its own connectivity
  - Sill depth determines exchange with outer coast
  - Freshwater lens depth modulates disease exposure
  - Real coordinates enable oceanographic connectivity models

## Possible Implementation
1. Use real fjord coordinates from GeoBC / USGS
2. Connectivity based on distance + sill characteristics
3. Each fjord gets temperature + salinity profile
4. Disease transmission modulated by both temperature and salinity
