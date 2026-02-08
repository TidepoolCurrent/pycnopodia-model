# Pycnopodia Geo Model Validation Report
**Date:** 2026-02-07
**Model:** geo_model.py (seasonal, 142 sites)

## Configuration
- 142 sites (75 fjords, 50 with freshwater lens)
- 9 regions: SE AK North/South, BC Fjords/Outer, Salish Sea, WA/OR, N/C/S California
- 320 time steps (80 years × 4 seasons)
- Seasonal dynamics: winter spawning, spring settlement, summer peak disease

## Hamilton 2021 Calibration (Year 17)

| Region | Model | Hamilton | Status |
|--------|-------|----------|--------|
| SE Alaska North | 96.6% | 96.0% | ✓ |
| SE Alaska South | 97.3% | 96.0% | ✓ |
| BC Fjords | 96.6% | 87.9% | ~ |
| BC Outer | 99.9% | 87.9% | ✗ |
| Salish Sea | 99.6% | 92.4% | ~ |
| WA/OR | 99.9% | 99.2% | ✓ |
| N. California | 99.9% | 99.2% | ✓ |
| C. California | 99.9% | 99.2% | ✓ |
| S. California | 99.9% | 99.2% | ✓ |

**Score: 8/9 within 10% of Hamilton targets**

## Recovery at Year 80
- SE AK North: 66.3% — strong fjord recovery
- SE AK South: 64.9% — good recovery with many inlets
- BC Fjords: 80.0% — strongest recovery (40 fjord sites)
- BC Outer: 12.9%
- Salish Sea: 42.5%
- WA/OR: 12.1%
- California: 7-17%

## Key Behaviors
1. Disease onset: Summer Y10 (The Blob), prevalence 90%
2. Population crash: 66.8M → 26.2M in one summer season
3. Acute phase: 4 years of sustained disease
4. Winter disease suppression in cold fjords
5. No sites go extinct (larval subsidy from dense fjord network)
6. Resistance evolution: 0.02 → 0.89 in SE AK North (44×)

## Known Issues
- BC Fjords and Salish decline too much vs Hamilton (need stronger cold/lens protection)
- BC Outer has no refugia mechanism (99.9% decline vs 87.9% target)
- Zero extinction with 142 sites may be too optimistic

