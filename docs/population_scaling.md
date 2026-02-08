# Population Scaling Analysis for Pycnopodia Model

**Date:** 2026-02-08  
**Current Status:** 160 named sites, ~142 million total carrying capacity  
**Historical Estimate:** 5.75-6.1 billion sunflower stars (Gravem et al. 2021)

---

## 1. Gravem et al. 2021 Methodology

### Study Details
**Citation:** Gravem, S.A., W.N. Heady, V.R. Saccomanno, K.F. Alvstad, A.L.M. Gehman, T.N. Frierson, and S.L. Hamilton. 2021. *Pycnopodia helianthoides* IUCN Red List Assessment.

### Methodology Summary

**Data Sources:**
- 61,000+ population surveys from 31 datasets
- Geographic extent: Baja California, Mexico (29°N) to Gulf of Alaska (61°N)
- Temporal range: Historical baseline (1976-2013 pre-SSWD) vs. current (2017-2020)

**Approach:**
1. **Density Estimates:** Mean density per region calculated from survey data
2. **Habitat Availability:** Estimated available rocky reef/kelp forest habitat along coastline
3. **Population Calculation:** Density × Habitat Area = Population estimate

**Key Findings:**
- **Historical density (pre-SSWD):** Varied by region, shallow water (<25m) densities ranged from ~0.01-0.1 individuals/m²
- **Population decline:** 90.6% decline calculated
- **Mortality estimate:** 5.75 billion individuals died
- **Pre-SSWD population:** ~6.1 billion (implied from 5.75B deaths + survivors)

**Geographic Extent:**
- Pacific Coast from Baja (29°N) to Gulf of Alaska (61°N)
- Approximately 25,000 km of coastline
- Includes both exposed outer coast and protected fjord/inland sea habitats

---

## 2. Recommended K Scaling

### Current Model Parameters
- **Sites:** 160
- **Current total K:** ~142 million
- **Current K per site (average):** 887,500

### Scaling to Historical Baseline

**Target:** ~6.1 billion total population (pre-SSWD estimate)

**Calculation:**
- Scaling factor: 6,100,000,000 / 142,000,000 = **43×**
- Recommended K per site (average): **38,125,000**

### Regional Considerations

The uniform scaling above is a **first approximation**. For more realistic K values, consider:

1. **Site Type Weighting:**
   - **Fjords** (especially deep, glacial fjords): May have supported **lower** densities due to:
     - Cold temperatures
     - Freshwater lens stratification  
     - Restricted circulation (shallow sills)
     - Suggested multiplier: **0.3-0.7× average K**
   
   - **Outer Coast/Islands** (exposed, rocky reef): Likely supported **higher** densities:
     - Optimal temperatures
     - Strong mixing/nutrient delivery
     - Extensive rocky substrate
     - Suggested multiplier: **1.2-2.0× average K**
   
   - **Inland Seas** (Salish Sea, Strait of Georgia): **Moderate** densities:
     - Intermediate conditions
     - Some areas highly productive
     - Suggested multiplier: **0.8-1.2× average K**

2. **Latitudinal Gradient:**
   - **Alaska/BC fjords (north):** Colder, potentially lower baseline K
   - **California (south):** Warmer, but at southern range limit, variable K
   - **Optimal zone (~45-55°N):** Highest baseline K

3. **Coastline Representation:**
   - Each site represents ~156 km of coastline (25,000 km / 160 sites)
   - Some sites represent large regions (e.g., "Central Coast outer islands")
   - Others are specific inlets with limited extent

### Recommended K Scaling Strategy

**Option A: Uniform Scaling (Simple)**
- K = 38,125,000 per site
- Total: 6.1 billion across 160 sites

**Option B: Type-Weighted Scaling (Moderate Complexity)**
```
Fjord sites (60): K = 15,000,000 each → 900M
Outer coast/islands (60): K = 60,000,000 each → 3,600M  
Inland seas (40): K = 40,000,000 each → 1,600M
Total: 6.1 billion
```

**Option C: Detailed Scaling (High Complexity)**
- Calculate K based on:
  - Actual site area/coastline length
  - Historical survey densities (if available for that location)
  - Habitat quality metrics (temperature, substrate type, depth)
- Requires additional data compilation

**Recommended:** Start with **Option B** (type-weighted), then refine specific high-priority sites using local data.

---

## 3. Geographic Gap Analysis

### Current Site Distribution
| Region | Sites | Notes |
|--------|-------|-------|
| SE Alaska North | 38 | Includes PWS, Kodiak, Aleutians |
| BC Fjords | 40 | Central/North Coast inlets |
| SE Alaska South | 19 | Sitka to Ketchikan |
| Salish Sea | 16 | Puget Sound, Georgia Strait |
| BC Outer | 12 | Vancouver Island, Haida Gwaii |
| S. California | 12 | Santa Barbara to Baja |
| WA/OR Outer | 9 | Olympic Coast to Cape Blanco |
| C. California | 8 | Monterey to Point Reyes |
| N. California | 6 | Mendocino to OR border |

### Major Gaps Identified

**1. Aleutian Islands (CRITICAL GAP)**
- **Gap size:** 2,364 km (Unalaska to Principe Channel)
- **Current coverage:** 3 sites (Unalaska, Akutan, Sanak)
- **Missing:** Entire central/western Aleutian chain
- **Ecological importance:** Cold-water refugia, historically had populations
- **Recommended additions:** 15-20 sites

**2. Gulf of Alaska / Alaska Peninsula**
- **Gap size:** 1,229-1,449 km (Kodiak to PWS via coast)
- **Current coverage:** Limited (Kodiak, Chignik, Sand Point)
- **Missing:** Long stretches of Alaska Peninsula south coast
- **Recommended additions:** 10-12 sites

**3. British Columbia Outer Coast**
- **Current coverage:** 12 sites, but clustering around major islands
- **Missing:** Stretches between major features
- **Recommended additions:** 8-10 sites

**4. California Coast**
- **Current coverage:** 26 sites total across N/C/S regions
- **Gaps:** 
  - Scattered coverage in Central California
  - Some Southern California gaps
- **Recommended additions:** 10-15 sites

**5. Washington/Oregon Coast**
- **Current coverage:** 9 sites
- **Gaps:** Long stretches of outer coast
- **Recommended additions:** 5-8 sites

---

## 4. Recommended New Sites

### Priority 1: Aleutian Islands & Western Alaska (20 sites)

#### Central/Western Aleutians
1. **Adak Island**, 51.88°N, -176.65°W, se_alaska_north, island
2. **Atka Island**, 52.20°N, -174.20°W, se_alaska_north, island
3. **Amchitka Island**, 51.50°N, -179.00°W, se_alaska_north, island
4. **Kiska Island**, 51.97°N, 177.50°E, se_alaska_north, island
5. **Attu Island** (westernmost), 52.85°N, 173.00°E, se_alaska_north, island
6. **Tanaga Island**, 51.80°N, -178.00°W, se_alaska_north, island
7. **Kanaga Island**, 51.90°N, -177.10°W, se_alaska_north, island
8. **Great Sitkin Island**, 52.08°N, -176.10°W, se_alaska_north, island
9. **Seguam Island**, 52.32°N, -172.50°W, se_alaska_north, island
10. **Yunaska Island**, 52.65°N, -170.67°W, se_alaska_north, island

#### Near Islands / Rat Islands
11. **Buldir Island**, 52.35°N, 175.95°E, se_alaska_north, island
12. **Agattu Island**, 52.85°N, 173.60°E, se_alaska_north, island
13. **Semisopochnoi Island**, 51.95°N, 179.60°E, se_alaska_north, island

#### Eastern Aleutians / Alaska Peninsula (fill gaps)
14. **Unimak Island (east)**, 54.75°N, -163.50°W, se_alaska_north, island
15. **Unimak Island (west)**, 54.35°N, -165.00°W, se_alaska_north, island
16. **Cold Bay / Izembek**, 55.20°N, -162.70°W, se_alaska_north, outer_coast
17. **Port Moller**, 56.00°N, -160.50°W, se_alaska_north, inland_sea
18. **Pavlof Bay**, 55.50°N, -161.70°W, se_alaska_north, outer_coast
19. **Stepovak Bay**, 55.85°N, -159.80°W, se_alaska_north, inland_sea
20. **Kupreanof Peninsula**, 56.05°N, -159.00°W, se_alaska_north, outer_coast

### Priority 2: Gulf of Alaska Coast (12 sites)

21. **Wide Bay**, 57.08°N, -155.00°W, se_alaska_north, outer_coast
22. **Aniakchak Bay**, 56.88°N, -158.15°W, se_alaska_north, outer_coast
23. **Port Heiden**, 56.95°N, -158.65°W, se_alaska_north, inland_sea
24. **Ugashik Bay**, 57.50°N, -157.40°W, se_alaska_north, inland_sea
25. **Becharof Lake outlet**, 57.78°N, -156.50°W, se_alaska_north, outer_coast
26. **Katmai Coast**, 58.35°N, -155.00°W, se_alaska_north, outer_coast
27. **Kamishak Bay**, 59.50°N, -153.50°W, se_alaska_north, inland_sea
28. **Augustine Island**, 59.37°N, -153.42°W, se_alaska_north, island
29. **Iliamna Bay**, 59.72°N, -153.20°W, se_alaska_north, inland_sea
30. **Chinitna Bay**, 59.88°N, -152.50°W, se_alaska_north, inland_sea
31. **Tuxedni Bay**, 60.15°N, -152.45°W, se_alaska_north, inner_coast
32. **Port Graham**, 59.35°N, -151.83°W, se_alaska_north, inland_sea

### Priority 3: British Columbia Outer Coast (10 sites)

33. **Scott Islands**, 50.95°N, -129.00°W, bc_outer, island
34. **Brooks Bay**, 50.12°N, -127.65°W, bc_outer, outer_coast
35. **Kyuquot**, 50.05°N, -127.35°W, bc_outer, outer_coast
36. **Clayoquot outer coast**, 49.00°N, -125.75°W, bc_outer, outer_coast
37. **Ucluelet**, 48.93°N, -125.55°W, bc_outer, outer_coast
38. **Tofino outer coast**, 49.12°N, -125.90°W, bc_outer, outer_coast
39. **Vargas Island**, 49.22°N, -126.10°W, bc_outer, island
40. **Flores Island**, 49.25°N, -126.18°W, bc_outer, island
41. **Estevan Point**, 49.38°N, -126.54°W, bc_outer, outer_coast
42. **Hesquiat Peninsula**, 49.42°N, -126.45°W, bc_outer, outer_coast

### Priority 4: California Coast (15 sites)

#### Northern California
43. **Patrick's Point**, 41.13°N, -124.16°W, n_california, outer_coast
44. **Humboldt Bay**, 40.77°N, -124.22°W, n_california, inland_sea
45. **False Klamath Cove**, 41.55°N, -124.08°W, n_california, outer_coast

#### Central California  
46. **Fort Ross**, 38.52°N, -123.25°W, c_california, outer_coast
47. **Bodega Bay**, 38.33°N, -123.05°W, c_california, inland_sea
48. **Tomales Bay**, 38.20°N, -122.95°W, c_california, inland_sea
49. **Año Nuevo**, 37.11°N, -122.33°W, c_california, outer_coast
50. **Pigeon Point**, 37.18°N, -122.39°W, c_california, outer_coast
51. **Point Lobos**, 36.52°N, -121.95°W, c_california, outer_coast
52. **Point Sur**, 36.30°N, -121.90°W, c_california, outer_coast

#### Southern California
53. **Point Sal**, 34.90°N, -120.67°W, s_california, outer_coast
54. **Point Arguello**, 34.58°N, -120.65°W, s_california, outer_coast  
55. **Point Conception**, 34.45°N, -120.47°W, s_california, outer_coast
56. **San Miguel Island**, 34.03°N, -120.36°W, s_california, island
57. **Santa Cruz Island (west)**, 34.01°N, -119.95°W, s_california, island

### Priority 5: Washington/Oregon Coast (8 sites)

58. **Cape Flattery**, 48.38°N, -124.72°W, wa_or_outer, outer_coast
59. **Neah Bay**, 48.37°N, -124.62°W, wa_or_outer, outer_coast
60. **Cape Alava**, 48.17°N, -124.73°W, wa_or_outer, outer_coast
61. **La Push**, 47.91°N, -124.64°W, wa_or_outer, outer_coast
62. **Kalaloch**, 47.61°N, -124.37°W, wa_or_outer, outer_coast
63. **Willapa Bay**, 46.65°N, -124.05°W, wa_or_outer, inland_sea
64. **Newport, OR**, 44.64°N, -124.05°W, wa_or_outer, outer_coast
65. **Yaquina Head**, 44.68°N, -124.08°W, wa_or_outer, outer_coast

### Additional Fills: Alaska/BC Fjords (15 sites)

#### Southeast Alaska inlets
66. **Stephens Passage (south)**, 57.25°N, -133.75°W, se_alaska_south, fjord, 120m
67. **Taku River delta**, 58.32°N, -133.95°W, se_alaska_north, fjord, 80m
68. **Whiting River**, 58.10°N, -133.20°W, se_alaska_north, fjord, 60m
69. **Port Houghton**, 57.32°N, -133.50°W, se_alaska_south, fjord, 50m
70. **Windham Bay**, 57.55°N, -133.30°W, se_alaska_south, fjord, 40m

#### British Columbia north coast fjords
71. **Muir Inlet (Glacier Bay)**, 58.90°N, -136.00°W, se_alaska_north, fjord, 15m
72. **Chichagof Bay**, 58.05°N, -136.06°W, se_alaska_north, fjord, 40m
73. **Portland Inlet**, 55.15°N, -130.10°W, bc_fjords, fjord, 40m
74. **Pearse Canal**, 55.20°N, -130.00°W, bc_fjords, fjord, 35m
75. **Wales Passage**, 54.98°N, -130.22°W, bc_fjords, inland_sea
76. **Chatham Sound**, 54.40°N, -130.50°W, bc_fjords, inland_sea
77. **Banks Island area**, 53.35°N, -129.85°W, bc_fjords, island
78. **Aristazabal Island**, 52.50°N, -128.80°W, bc_fjords, island
79. **Milbanke Sound inner**, 52.25°N, -128.70°W, bc_fjords, inland_sea
80. **Hecate Strait (west)**, 53.50°N, -131.00°W, bc_outer, outer_coast

---

## 5. Implementation Recommendations

### Phase 1: Critical Gaps (20-25 sites)
Focus on Aleutian Islands and Gulf of Alaska coast to address the largest geographic gaps and ensure representation of coldwater, western populations.

### Phase 2: Outer Coast Enhancement (25-30 sites)
Add sites along exposed coastlines (BC, WA, OR, CA) to better represent high-productivity outer coast habitat.

### Phase 3: Fine-Scale Resolution (20-25 sites)
Fill remaining gaps and add detail to regions with complex bathymetry or high habitat heterogeneity.

### Total Recommended
**Current:** 160 sites  
**Additions:** 70-80 sites  
**Final:** 230-240 sites

This would give an average spacing of ~100-110 km per site, improving representation while remaining computationally tractable.

### K Assignment Strategy
For new sites, assign K based on:
1. Site type (fjord vs. outer coast vs. inland sea)
2. Regional productivity patterns
3. Historical survey data if available
4. Use type-weighted multipliers from Option B above

---

## 6. Data Needs

To refine K estimates further:
1. **Historical survey densities** by specific location (pre-SSWD)
2. **Habitat area estimates** (rocky reef extent, kelp forest coverage)
3. **Depth-stratified densities** (shallow <25m vs. deeper 25-100m)
4. **Substrate type** (bedrock, boulder, cobble suitability)
5. **Temperature/oceanographic** data for each site

---

## Summary

- **Gravem et al. 2021** used survey data + habitat availability to estimate ~6.1B historical population
- **Current model** has 160 sites, ~142M total K → needs ~43× scaling to match historical baseline
- **Recommended K per site:** 38M (uniform) or type-weighted (15-60M based on habitat)
- **Major geographic gaps:** Aleutian Islands, Gulf of Alaska, California outer coast
- **Recommended additions:** 70-80 new sites, focusing on Aleutians (20), Gulf of Alaska (12), BC/CA/WA coasts (38)

**Next Steps:**
1. Decide on K scaling approach (Option A, B, or C)
2. Prioritize site additions based on model goals (refugia identification? connectivity? full range representation?)
3. Compile site-specific K assignments
4. Update `real_sites.py` with new sites
5. Re-run population model with adjusted parameters
