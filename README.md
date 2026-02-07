# Pycnopodia Population Model

**Ratio-based outputs for interpretable comparisons across scenarios.**

Individual-based simulation of *Pycnopodia helianthoides* (sunflower sea star) population dynamics, integrating:

- **Sweepstakes Reproductive Success (SRS)** — Hedgecock & Pudovkin 2011
- **Polygenic SSWD Resistance** — Based on Schiebelhut et al. 2024 (51 associated loci)
- **Broadcast Spawner Allee Effects** — Gamete dilution at low density
- **Captive Breeding & Outplanting** — Conservation intervention modeling

## Installation

```bash
pip install numpy
cd pycnopodia-model
pip install -e .
```

## Quick Start

```python
from pycnopodia import Simulation, Config

# Run single simulation with defaults
sim = Simulation()
result = sim.run()

# All outputs are RATIOS relative to baseline
print(f"Final N/N₀: {result.final_n_ratio:.1%}")    # Population vs baseline
print(f"Final H/H₀: {result.final_h_ratio:.1%}")    # Diversity retention
print(f"Min N/N₀: {result.min_n_ratio:.1%}")        # Bottleneck depth
print(f"Extinct: {result.extinct}")
print(f"Recovered to 30%: {result.recovered(0.30)}")
```

### Why Ratios?

Absolute numbers (e.g., "population = 847") are arbitrary—they depend on starting conditions. Ratios tell the real story:

- **N/N₀ = 0.05** means the population crashed to 5% of baseline
- **H/H₀ = 0.27** means 73% of genetic diversity was lost
- **Ne/N = 0.10** means effective size is 10× smaller than census

These are directly comparable across scenarios regardless of starting population size.

## Command Line

```bash
# Run baseline scenario (50 replicates)
python run_simulation.py

# No intervention (extinction risk)
python run_simulation.py --no-outplanting

# Metapopulation model
python run_simulation.py --preset metapopulation --n-replicates 100
```

## Configuration Presets

| Preset | Description |
|--------|-------------|
| `baseline` | Default parameters with outplanting |
| `no_disease` | Control without SSWD |
| `extreme_srs` | Only 2% of adults breed (extreme sweepstakes) |
| `wright_fisher` | All adults breed equally (no SRS) |
| `high_outplanting` | 500 juveniles/year, 50 broodstock |
| `limited_broodstock` | 6 parents only (genetic bottleneck) |
| `metapopulation` | 5 subpopulations (AK → CA) |
| `climate_warming` | Temperature-dependent SSWD dynamics |

## Model Components

### Demography

Age-structured population with stage-specific survival:

| Stage | Ages | Annual Survival |
|-------|------|-----------------|
| Settled juvenile | 0 | 15% |
| Juvenile | 1-4 | 70% |
| Adult | 5-18 | 88% |
| Senescent | 19+ | 26% |

### Sweepstakes Reproductive Success

Each year, only a fraction of adults successfully contribute offspring:

```
breeding_fraction ~ Beta(α, β)
mean = 8%, high variance
```

This creates effective population sizes (Ne) 10-100× smaller than census size.

### SSWD Disease Dynamics

Prevalence follows epidemic → endemic trajectory:

```
P(t) = P_endemic + (P_peak - P_endemic) × exp(-δ × (t - t_onset))
```

**Improvements over original:**
- Optional temperature-dependent mortality
- Optional density-dependent transmission
- Subpopulation-specific severity (thermal refugia)

### Polygenic Resistance

20 unlinked loci with additive effects:

```
mortality_reduction = Σ (genotype_i × effect_i)
```

Each resistance allele reduces SSWD mortality by 3.5% (default).

### Allee Effects

Two fertilization models:

1. **Saturating** (default): `f = N² / (N² + h²)`
2. **Levitan mechanistic**: Sperm-egg collision kinetics

### Outplanting Intervention

Captive-bred juveniles added annually:

- Configurable broodstock size and genetic diversity
- Optional inbreeding tracking
- Broodstock refreshment from wild population

## Example: Intervention Comparison

```bash
# With outplanting
python run_simulation.py --n-replicates 20

# Without outplanting  
python run_simulation.py --n-replicates 20 --no-outplanting
```

| Metric | With Outplanting | Without |
|--------|------------------|---------|
| Extinction | 0% | **75%** |
| Final N/N₀ | 141% | **5%** |
| Final H/H₀ | 179% | **27%** |
| Bottleneck (min N/N₀) | 23% | **1%** |
| Bottleneck year | 14 | 35 |

The ratio-based output makes the intervention effect immediately clear: without outplanting, populations crash to 1% of baseline and most go extinct.

## Key Findings

1. **SRS is the critical bottleneck** — Even modest SRS (8% breeding/year) reduces Ne to ~10% of census N, slowing adaptation.

2. **Outplanting is necessary** — Without intervention, extinction probability exceeds 80% under realistic SRS.

3. **Broodstock diversity matters** — ≥30 parents needed for evolutionary (not just demographic) rescue.

4. **Timing is critical** — Intervention within 5-10 years of epidemic onset is far more effective.

5. **Southern populations need direct intervention** — Adult migration (0.5%/year) cannot rescue collapsed California/Oregon populations from healthy Alaska.

## Improvements Over Original Model

| Feature | Original | This Version |
|---------|----------|--------------|
| Allee model | Saturating only | + Levitan mechanistic option |
| Disease | Fixed trajectory | + Temperature-dependent |
| Transmission | — | + Density-dependent option |
| Locus effects | Uniform | + Variable effect sizes |
| Broodstock | Static | + Inbreeding tracking |
| Thermal refugia | Implicit | + Explicit temperature gradient |

## References

1. Hedgecock D, Pudovkin AI (2011). Sweepstakes reproductive success in highly fecund marine fish and shellfish. *Bull Mar Sci* 87:971-1002.

2. Schiebelhut LM, Puritz JB, Dawson MN (2024). Reference genome of the sunflower sea star. *J Heredity* esae002.

3. Gravem SA et al. (2021). Pycnopodia helianthoides. *IUCN Red List* e.T178290276A197818455.

4. Hodin J et al. (2021). Culturing the sunflower sea star. *Aquaculture* 542:736873.

5. Levitan DR (1991). Influence of body size and population density on fertilization success. *Ecology* 72:1712-1728.

## License

MIT

## Authors

- Original model: Weertman 2026
- Code implementation: TidepoolCurrent

*Built for kelp forest restoration research at Friday Harbor Labs.*
