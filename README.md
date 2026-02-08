# Pycnopodia helianthoides Population Model

Spatially explicit metapopulation model of sunflower sea star (*Pycnopodia helianthoides*) population dynamics under Sea Star Wasting Disease (SSWD).

## Overview

This model simulates the catastrophic decline (>94% range-wide) of sunflower stars following the 2013 SSWD outbreak, driven by the marine heatwave known as "The Blob." It explores recovery pathways through fjord refugia, resistance evolution, and larval connectivity.

### Key Features

- **142 real Pacific Coast sites** — named locations from Isla Guadalupe (29°N) to Prince William Sound (60.7°N)
- **75 fjord sites** with documented sill depths and freshwater lens effects
- **Seasonal dynamics** — 4 time steps per year (winter spawning, spring settlement, summer peak disease, fall persistence)
- **Distance-based connectivity** — separate larval dispersal and disease transmission networks (142×142 matrices)
- **Polygenic resistance evolution** — 10 loci with per-locus selection
- **Allee effects** on broadcast spawner reproduction
- **Marine heatwave** — The Blob (2013-2016) as disease trigger
- **Fjord refugia** — sill depth modulates disease transmission, freshwater lens reduces mortality

### Model Hierarchy

1. **`pycnopodia/geo_model.py`** — PRIMARY: 142 real sites, seasonal time steps, distance-based connectivity
2. **`pycnopodia/pacific_coast.py`** — ABSTRACT: 400 sites across 9 regions (annual steps)
3. **`pycnopodia/network.py`** — THEORETICAL: connectivity scenario exploration

## Quick Start

```bash
# Single run (142 sites, 80 years, 320 seasonal steps)
python run_geo.py

# Ensemble (10 seeds)
python run_geo.py --ensemble 10

# Interactive dashboard
cd dashboard && python serve.py
# Open http://localhost:8090
```

## Regions (9)

| Region | Sites | Lat Range | Key Feature |
|--------|-------|-----------|-------------|
| SE Alaska North | 23 | 57-61°N | Fjord refugia, PWS, glacial |
| SE Alaska South | 19 | 55-57°N | Inside passage, LeConte Bay |
| BC Fjords | 40 | 50-55°N | Dense fjord network, Gehman 2025 refugia |
| BC Outer Coast | 12 | 48-55°N | Vancouver Island, Haida Gwaii |
| Salish Sea | 16 | 47-50°N | Friday Harbor Labs, Saanich Inlet |
| WA/OR Coast | 9 | 43-48°N | Exposed outer coast |
| N. California | 6 | 39-42°N | Mendocino to Crescent City |
| C. California | 8 | 35-39°N | Monterey to Point Reyes |
| S. California | 9 | 29-35°N | Channel Islands to Baja |

## Seasonal Dynamics

| Season | Biology | Disease |
|--------|---------|---------|
| Winter | Broadcast spawning | Cold suppression (0.5×) |
| Spring | Larval settlement | Moderate (0.8×) |
| Summer | Growth, peak temps | **Peak transmission (1.3×)** |
| Fall | Pre-winter | Still active (1.1×) |

## Calibration

Validated against Hamilton et al. 2021 regional decline estimates (2013-2020):

| Region | Model (Y17) | Hamilton 2021 | Status |
|--------|-------------|---------------|--------|
| SE Alaska | 96.6-97.3% | 96.0% | ✓ |
| BC | 96.6-99.9% | 87.9% | ~ |
| Salish Sea | 99.6% | 92.4% | ~ |
| South of Salish | 99.9% | 99.2% | ✓ |

## Key References

- Hamilton et al. 2021 — Range-wide SSWD decline estimates
- Gehman et al. 2025 — BC fjord refugia confirmed (Proc R Soc B)
- Schiebelhut et al. 2024 — Polygenic resistance
- Hedgecock & Pudovkin 2011 — Sweepstakes reproductive success

## Interactive Dashboard

The dashboard (`dashboard/`) provides a Leaflet.js map visualization:
- Real map tiles (CartoDB dark) with zoom/pan
- 142 sites at actual coordinates
- Play/pause animation through 320 seasonal steps
- Regional statistics sidebar
- Disease (red), resistance (green ring), population (node size)

```bash
cd dashboard && python serve.py
# http://localhost:8090
```

## Repository Structure

```
pycnopodia-model/
├── pycnopodia/
│   ├── geo_model.py          # PRIMARY: 142 sites, seasonal
│   ├── pacific_coast.py      # ABSTRACT: 400 sites, annual
│   ├── network.py            # THEORETICAL: connectivity
│   ├── config.py             # Configuration
│   └── real_data.py          # Historical SST data
├── data/
│   ├── real_sites.py         # 142 named sites with coordinates
│   └── coastline.py          # Coastline generation
├── dashboard/                # Interactive Leaflet.js dashboard
├── figures/                  # Generated figures
│   ├── geo/                  # Geography model outputs
│   ├── pacific_coast/        # Abstract model outputs
│   └── investigation/        # Network analysis
├── tests/                    # Test suite (69 tests)
├── docs/                     # Design documents
├── run_geo.py                # Primary model runner
├── run_pacific_coast.py      # Abstract model runner
├── VALIDATION.md             # Validation report
└── INVESTIGATION.md          # Network model investigation
```

## License

Research code. Contact WLWeertman / Friday Harbor Labs for usage.
