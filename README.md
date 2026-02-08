# Pycnopodia Population Model

Individual-based metapopulation model for sunflower sea star (*Pycnopodia helianthoides*) population dynamics under Sea Star Wasting Disease (SSWD).

## What It Does

Simulates population crash (2013–2016) and potential recovery across the Pacific Coast, from Baja California to the Gulf of Alaska. Key features:

- **142 real geographic sites** spanning 29°N–61°N, including 75 named fjords
- **Geographic disease spread** from Washington epicenter (~47.5°N), not simultaneous onset
- **Seasonal dynamics** — quarterly time steps with winter spawning, summer disease peaks
- **Fjord refugia** — sill-depth modulated disease transmission, freshwater lens protection (Gehman et al. 2025)
- **Polygenic resistance** — 10 loci under frequency-dependent selection
- **Broadcast spawner Allee effects** — fertilization success depends on local adult density
- **Distance-based connectivity** — separate larval and disease dispersal networks

## Calibration

Validated against Hamilton et al. (2021) range-wide survey data:

| Region | Model Y17 Decline | Hamilton 2021 |
|--------|-------------------|---------------|
| SE Alaska | 92% | 96% |
| BC Fjords | 92% | 88% |
| Salish Sea | 97% | 92% |
| WA/OR Coast | 99% | >99% |
| California | 99% | >99% |

## Quick Start

```bash
# Run single simulation with figures
python3 run_geo.py --ensemble 1

# Run 10-seed ensemble
python3 run_geo.py --ensemble 10

# Output: figures/geo/ and dashboard/data.json
```

## Interactive Dashboard

Leaflet.js map with real site locations, disease spread animation, and population trajectories.

```bash
cd dashboard && python3 serve.py
# Open http://localhost:8090
```

Controls: play/pause, speed, year slider, connectivity toggle, region statistics.

## Project Structure

```
pycnopodia/
  __init__.py
  geo_model.py        # Core simulation (GeoConfig, GeoSimulation, GeoResult)
data/
  real_sites.py        # 142 named sites with lat/lon, sill depth, freshwater lens
  coastline.py         # Dense site network generation
dashboard/
  index.html           # Leaflet.js interactive map
  viz.js               # D3.js visualization
  serve.py             # Local HTTP server
run_geo.py             # CLI runner with ensemble + figures
docs/                  # Design documents and calibration notes
figures/               # Generated plots
```

## Key Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| Disease onset | Year 10 (= 2013) | Hewson et al. 2014 |
| Base mortality | 90% annual | Montecino-Latorre et al. 2016 |
| Acute phase | 3 years | Hamilton et al. 2021 |
| Allee threshold | 50 adults | Levitan 1991 (broadcast spawners) |
| Maturation | 3 years | Shivji et al. 1983 |
| Resistance loci | 10, max effect 70% | — |
| Freshwater lens reduction | 50% disease, 30% mortality | Gehman et al. 2025 |

## Key References

- Hamilton et al. 2021. *Biol. Conserv.* — Range-wide decline survey
- Gehman et al. 2025. *Proc. R. Soc. B.* — Fjord refugia, freshwater lens mechanism
- Hewson et al. 2014. *PNAS* — SSWD etiology
- Montecino-Latorre et al. 2016. *PLoS ONE* — SSWD epidemiology

## License

MIT
