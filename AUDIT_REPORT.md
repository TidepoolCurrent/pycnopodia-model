# Pycnopodia Model Repository Audit Report

**Date:** 2026-02-07  
**Auditor:** Subagent repo-audit

## Executive Summary

The repository has evolved through multiple model iterations:
1. **Original individual-based model** (`population.py`, `simulation.py`) → Single population genetics
2. **Network metapopulation model** (`network.py`) → 1000 abstract sites with connectivity scenarios
3. **Pacific Coast abstract model** (`pacific_coast.py`) → 400 sites across 9 regions
4. **Geographic realism model** (`geo_model.py`) → 56 real named Pacific Coast sites ✅ **CURRENT PRIMARY**

Many files from earlier iterations remain and should be cleaned up.

---

## 🔴 HIGH PRIORITY: Delete These Files

### Obsolete Model Code

#### 1. **`pycnopodia/population.py`** (11,255 bytes)
- **Why:** Individual-based genetics model superseded by site-based models
- **Replaced by:** `geo_model.py` and `pacific_coast.py` handle population dynamics at site level
- **Still referenced by:** 
  - `simulation.py` (also obsolete)
  - `reproduction.py` (only used by simulation.py)
  - `disease.py` (only used by simulation.py)
  - `intervention.py` (only used by simulation.py)
- **Verdict:** DELETE (part of obsolete individual-based modeling stack)

#### 2. **`pycnopodia/simulation.py`** (13,415 bytes)
- **Why:** Orchestrates the obsolete individual-based model
- **Replaced by:** `GeoSimulation` in `geo_model.py`, `PacificCoastSimulation` in `pacific_coast.py`
- **Still used by:** `run_simulation.py` only
- **Verdict:** DELETE (only survives via obsolete runner script)

#### 3. **`pycnopodia/reproduction.py`** (9,938 bytes)
- **Why:** Individual-level Mendelian reproduction for population.py
- **Current models:** Use site-level allele frequencies, not individual genomes
- **Verdict:** DELETE (only used by obsolete simulation.py)

#### 4. **`pycnopodia/disease.py`** (4,964 bytes)
- **Why:** Individual-level SSWD mortality for population.py
- **Current models:** Site-level disease dynamics in geo_model.py and pacific_coast.py
- **Verdict:** DELETE (only used by obsolete simulation.py)

#### 5. **`pycnopodia/intervention.py`** (8,311 bytes)
- **Why:** Broodstock and outplanting for individual-based model
- **Current models:** Have their own intervention logic
- **Verdict:** DELETE (only used by obsolete simulation.py)

#### 6. **`pycnopodia/spatial.py`** (8,891 bytes)
- **Why:** Spatial connectivity module that appears unused
- **Check:** Not imported by any active model (network.py, geo_model.py, pacific_coast.py)
- **Grep results:** Only in `__init__.py` export and tests
- **Verdict:** DELETE (superseded by connectivity in current models)

#### 7. **`pycnopodia/environment.py`** (6,675 bytes)
- **Why:** Environmental variation module
- **Check:** Only used in tests, not in active models
- **Verdict:** DELETE (functionality absorbed into region-specific configs)

---

### Obsolete Runner Scripts

#### 8. **`run_simulation.py`** (6,615 bytes)
- **Why:** Runs the obsolete individual-based population.py model
- **Replaced by:** `run_geo.py` (primary) and `run_pacific_coast.py` (abstract)
- **Verdict:** DELETE

#### 9. **`run_network.py`** (5,215 bytes)
- **Why:** Runs the abstract 1000-site network model
- **Status:** network.py is still valid for theoretical connectivity studies
- **BUT:** Investigation already completed (see INVESTIGATION.md)
- **Verdict:** ⚠️ **KEEP** (network.py still useful for theory, but this runner is covered by investigate.py)
- **Alternative:** Could delete if `investigate.py` covers all use cases

#### 10. **`run_ensemble_check.py`** (3,117 bytes)
- **Why:** One-off ensemble test for pacific_coast.py parameter tuning
- **Status:** Test complete, results documented
- **Verdict:** DELETE (one-off analysis, no ongoing utility)

#### 11. **`investigate.py`** (28,789 bytes)
- **Why:** Comprehensive investigation script that generated INVESTIGATION.md
- **Status:** Investigation complete, all figures generated
- **Output:** 25 figures in `figures/investigation/`
- **Verdict:** ⚠️ **ARCHIVE or KEEP** (useful for re-running analyses, but completed its mission)
- **Recommendation:** Move to `scripts/investigate.py` if keeping

---

### Obsolete Test Files

#### 12. **`tests/test_model.py`** (6,118 bytes)
- **Why:** Tests the obsolete Population and simulation.py classes
- **Tests:** Config, Population, Simulation — all part of obsolete individual-based model
- **Verdict:** DELETE (tests obsolete code)

#### 13. **`tests/test_environment.py`** (3,646 bytes)
- **Why:** Tests the unused environment.py module
- **Verdict:** DELETE (tests obsolete code)

#### 14. **`tests/test_size_structured.py`** (5,351 bytes)
- **Why:** Tests size_structured.py which is not used in current models
- **Status:** size_structured.py only imported in tests and `__init__.py`
- **Verdict:** DELETE (along with pycnopodia/size_structured.py)

#### 15. **`tests/test_srs.py`** (4,433 bytes)
- **Why:** Tests srs.py (hierarchical sweepstakes reproductive success)
- **Status:** srs.py only imported in tests and `__init__.py`
- **Verdict:** DELETE (along with pycnopodia/srs.py)

---

### Unused Utility Modules

#### 16. **`pycnopodia/size_structured.py`** (11,269 bytes)
- **Why:** Size-structured population dynamics
- **Usage:** Only imported in tests and `__init__.py`, never used by active models
- **Verdict:** DELETE (future feature that was never integrated)

#### 17. **`pycnopodia/srs.py`** (6,931 bytes)
- **Why:** Hierarchical sweepstakes reproductive success
- **Usage:** Only imported in tests and `__init__.py`, never used by active models
- **Verdict:** DELETE (theoretical module never integrated)

---

## 🟡 MEDIUM PRIORITY: Stale Figures

### Top-Level Figures (Old Model Outputs)

#### Generated by obsolete `scripts/generate_figures.py`:

18. **`figures/allee_effect.png`** (104 KB, Feb 7 09:15)
- **From:** Individual-based model illustration
- **Verdict:** DELETE (use pacific_coast or geo figures instead)

19. **`figures/disease_dynamics.png`** (84 KB, Feb 7 09:15)
- **From:** Individual-based model
- **Verdict:** DELETE

20. **`figures/genetic_diversity.png`** (107 KB, Feb 7 09:15)
- **From:** Individual-based model
- **Verdict:** DELETE

21. **`figures/model_overview.png`** (79 KB, Feb 7 09:15)
- **From:** Network model overview (may still be useful)
- **Verdict:** ⚠️ **REVIEW** (if still accurate, keep; otherwise regenerate)

22. **`figures/population_trajectory.png`** (117 KB, Feb 7 09:15)
- **From:** Individual-based model
- **Verdict:** DELETE

23. **`figures/srs_effect.png`** (109 KB, Feb 7 09:15)
- **From:** Individual-based model
- **Verdict:** DELETE

24. **`figures/test.png`** (18 KB, Feb 7 00:16)
- **Why:** Test output file
- **Verdict:** DELETE

---

### Duplicate/Superseded Figures

25-32. **`figures/fig1_trajectory.png` through `fig8_stochasticity.png`** (Various dates)
- **Why:** Two sets of similarly-named figures from different model runs
- Newer versions (Feb 7 09:15-09:16): `fig1_trajectories.png`, `fig2_enhanced_vs_wild.png`, etc.
- Older versions (Feb 7 00:13-00:16): `fig1_trajectory.png`, `fig2_intensity.png`, etc.
- **Verdict:** DELETE older versions (00:13-00:16 timestamps)

**Specifically DELETE:**
- `figures/fig1_trajectory.png` (Feb 7 00:13) — superseded by `fig1_trajectories.png`
- `figures/fig2_intensity.png` (Feb 7 00:14) — different content than `fig2_enhanced_vs_wild.png`
- `figures/fig3_enhanced_vs_wild.png` (Feb 7 00:14) — duplicate naming
- `figures/fig4_refugia.png` (Feb 7 00:15) — superseded by `fig4_heatmap.png`
- `figures/fig5_mortality.png` (Feb 7 00:15) — duplicate
- `figures/fig6_connectivity.png` (Feb 7 00:16) — superseded by `fig6_mortality.png`
- `figures/fig7_timing.png` (Feb 7 00:16)
- `figures/fig8_stochasticity.png` (Feb 7 00:16) — different from `fig5_stochasticity.png`

---

### Investigation Figures (Completed Analysis)

**Status:** `figures/investigation/` contains 27 figures from the network model investigation (completed Feb 7 09:22-09:45)

**Decision:** ✅ **KEEP** — These document the theoretical network analysis and remain valid reference material

---

### Pacific Coast Figures

**Status:** `figures/pacific_coast/` contains 16 figures from the abstract 400-site model (Feb 7 13:09-15:54)

**Decision:** ✅ **KEEP** — Current abstract model outputs

---

### Geo Model Figures

**Status:** `figures/geo/` contains 5 figures from the real-geography model (Feb 7 17:11-17:24)

**Decision:** ✅ **KEEP** — Current primary model outputs

---

## 🟢 LOW PRIORITY: Documentation Updates

### Outdated Documentation

33. **`docs/NETWORK_MODEL_PROPOSAL.md`** (7,472 bytes)
- **Status:** Historical proposal that led to network.py implementation
- **Verdict:** ⚠️ **ARCHIVE** (move to `docs/archive/` or add "Historical" header)

### Completed Investigations

34. **`INVESTIGATION.md`** (12,568 bytes)
- **Status:** Completed analysis report (valuable reference)
- **Verdict:** ✅ **KEEP** (documents network model investigation)

---

### Documentation Needing Updates

35. **`README.md`**
- **Issue:** Still references individual-based model as co-equal with network model
- **Reality:** geo_model.py is now primary, pacific_coast.py is abstract reference
- **Action:** UPDATE to reflect current model hierarchy:
  1. **Primary:** `geo_model.py` (56 real sites)
  2. **Abstract:** `pacific_coast.py` (400 sites, 9 regions)
  3. **Theoretical:** `network.py` (connectivity scenarios)
  4. **Obsolete:** population.py/simulation.py (delete)

36. **`README.md` - Dashboard Section**
- **Status:** Dashboard section references `run_pacific_coast.py --export-json`
- **Reality:** Dashboard works but may need updating for geo_model.py
- **Action:** Verify dashboard works with geo_model.py or clarify it's pacific_coast only

---

## 🔵 KEEP: Active Code

### Core Models (Current)

✅ **`pycnopodia/geo_model.py`** (29,917 bytes) — PRIMARY MODEL (56 real sites)  
✅ **`pycnopodia/pacific_coast.py`** (51,572 bytes) — ABSTRACT MODEL (400 sites, 9 regions)  
✅ **`pycnopodia/network.py`** (45,693 bytes) — THEORETICAL MODEL (connectivity scenarios)

### Runner Scripts (Current)

✅ **`run_geo.py`** (10,159 bytes) — Primary model runner  
✅ **`run_pacific_coast.py`** (58,403 bytes) — Abstract model runner

### Supporting Modules

✅ **`pycnopodia/config.py`** (7,502 bytes) — Configuration classes  
✅ **`pycnopodia/__init__.py`** (1,533 bytes) — Package initialization  
✅ **`pycnopodia/real_data.py`** (10,153 bytes) — Real-world data integration

### Data

✅ **`data/real_sites.py`** (12,011 bytes) — Real Pacific Coast site coordinates  
✅ **`data/coastline.py`** (12,757 bytes) — Coastline and site network generation

### Tests (Keep These)

✅ **`tests/test_network.py`** (20,056 bytes) — Tests network.py  
✅ **`tests/test_pacific_coast.py`** (23,648 bytes) — Tests pacific_coast.py

### Scripts

✅ **`scripts/generate_figures.py`** (9,577 bytes) — Figure generation utility  
✅ **`scripts/quick_analysis.py`** (9,707 bytes) — Quick network analysis  
✅ **`scripts/exhaustive_analysis.py`** (24,387 bytes) — Comprehensive network analysis  
✅ **`scripts/visualize_disease_spread.py`** (7,803 bytes) — Disease visualization

### Documentation (Keep)

✅ **`docs/pacific_coast_design.md`** (10,279 bytes)  
✅ **`docs/real_geography.md`** (2,003 bytes)  
✅ **`docs/hamilton2021_numbers.md`** (2,511 bytes)  
✅ **`docs/fjord_subnetwork_design.md`** (2,737 bytes)  
✅ **`docs/allee_sensitivity.md`** (1,983 bytes)

### Dashboard

✅ **`dashboard/README.md`**  
✅ **`dashboard/TEST_RESULTS.md`**  
✅ **`dashboard/serve.py`**  
✅ **`dashboard/index.html`** (presumed)  
✅ **`dashboard/viz.js`** (presumed)

---

## Summary Statistics

### Files to DELETE: 34 files

**Python modules:** 11 files (~85 KB)
- population.py, simulation.py, reproduction.py, disease.py, intervention.py
- spatial.py, environment.py, size_structured.py, srs.py

**Runner scripts:** 3 files (~18 KB)
- run_simulation.py, run_ensemble_check.py

**Tests:** 4 files (~20 KB)
- test_model.py, test_environment.py, test_size_structured.py, test_srs.py

**Figures:** 15 files (~1.5 MB)
- 6 top-level figures from individual-based model
- 8 older fig*.png versions
- 1 test.png

**Documentation:** 1 file (7 KB)
- NETWORK_MODEL_PROPOSAL.md (or archive it)

### Files to UPDATE: 1 file

- README.md (clarify model hierarchy, update dashboard instructions)

### Files to KEEP: ~40 files

- 3 core models (geo, pacific_coast, network)
- 2 runner scripts (run_geo.py, run_pacific_coast.py)
- Supporting modules (config, real_data, data/)
- Valid tests (test_network.py, test_pacific_coast.py)
- All current documentation
- Dashboard
- Current figures (investigation/, pacific_coast/, geo/)

---

## Recommended Cleanup Sequence

### Phase 1: Safe Deletions (No Dependencies)

```bash
# Delete obsolete runner scripts
rm run_simulation.py run_ensemble_check.py

# Delete stale figures
rm figures/test.png
rm figures/allee_effect.png figures/disease_dynamics.png figures/genetic_diversity.png
rm figures/population_trajectory.png figures/srs_effect.png

# Delete old figure versions
rm figures/fig1_trajectory.png figures/fig2_intensity.png figures/fig3_enhanced_vs_wild.png
rm figures/fig4_refugia.png figures/fig5_mortality.png figures/fig6_connectivity.png
rm figures/fig7_timing.png figures/fig8_stochasticity.png

# Delete obsolete tests
rm tests/test_model.py tests/test_environment.py
rm tests/test_size_structured.py tests/test_srs.py

# Archive historical docs
mkdir -p docs/archive
mv docs/NETWORK_MODEL_PROPOSAL.md docs/archive/
```

### Phase 2: Delete Obsolete Model Stack

```bash
# Delete individual-based model modules
rm pycnopodia/population.py
rm pycnopodia/simulation.py
rm pycnopodia/reproduction.py
rm pycnopodia/disease.py
rm pycnopodia/intervention.py

# Delete unused utility modules
rm pycnopodia/spatial.py
rm pycnopodia/environment.py
rm pycnopodia/size_structured.py
rm pycnopodia/srs.py
```

### Phase 3: Update Package Exports

Edit `pycnopodia/__init__.py` to remove imports:
- Population, simulation classes
- SizeStructuredConfig, SizeStructuredPopulation
- HierarchicalSRSConfig, HierarchicalSRS
- SpatialConfig, SpatialPopulation

### Phase 4: Update Documentation

Update `README.md`:
- Remove references to individual-based model
- Clarify model hierarchy (geo > pacific_coast > network)
- Update dashboard instructions if needed

### Phase 5: Optional

Move or delete `investigate.py` (28 KB):
- Option A: Move to `scripts/investigate.py` (if keeping for re-runs)
- Option B: Delete (investigation complete, INVESTIGATION.md preserved)

---

## Final Repository Structure (After Cleanup)

```
pycnopodia-model/
├── pycnopodia/
│   ├── __init__.py           # Updated exports
│   ├── config.py             ✅
│   ├── geo_model.py          ✅ PRIMARY
│   ├── pacific_coast.py      ✅ ABSTRACT
│   ├── network.py            ✅ THEORETICAL
│   └── real_data.py          ✅
├── data/
│   ├── coastline.py          ✅
│   └── real_sites.py         ✅
├── tests/
│   ├── test_network.py       ✅
│   └── test_pacific_coast.py ✅
├── scripts/
│   ├── generate_figures.py   ✅
│   ├── quick_analysis.py     ✅
│   ├── exhaustive_analysis.py ✅
│   └── visualize_disease_spread.py ✅
├── figures/
│   ├── investigation/        ✅ (27 files)
│   ├── pacific_coast/        ✅ (16 files)
│   ├── geo/                  ✅ (5 files)
│   ├── fig1_trajectories.png ✅
│   ├── fig2_enhanced_vs_wild.png ✅
│   ├── fig3_refugia.png      ✅
│   ├── fig4_heatmap.png      ✅
│   ├── fig5_stochasticity.png ✅
│   └── fig6_mortality.png    ✅
├── dashboard/                ✅
├── docs/                     ✅ (+ archive/ subdir)
├── run_geo.py                ✅ PRIMARY
├── run_pacific_coast.py      ✅ ABSTRACT
├── INVESTIGATION.md          ✅
└── README.md                 📝 UPDATE

REMOVED:
- population.py, simulation.py, reproduction.py, disease.py, intervention.py
- spatial.py, environment.py, size_structured.py, srs.py
- run_simulation.py, run_ensemble_check.py, investigate.py
- test_model.py, test_environment.py, test_size_structured.py, test_srs.py
- 15 stale/duplicate figures
```

---

## Validation Checklist

After cleanup, verify:

- [ ] `python run_geo.py` still works
- [ ] `python run_pacific_coast.py` still works
- [ ] `python -m pytest tests/` passes (only test_network.py and test_pacific_coast.py)
- [ ] No import errors in `pycnopodia/__init__.py`
- [ ] Dashboard still works (if kept)
- [ ] README accurately reflects current models
- [ ] Git history preserved (use `git rm`, not `rm`)

---

**End of Audit Report**
