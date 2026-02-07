# Pycnopodia Population Model

Network-based metapopulation model for *Pycnopodia helianthoides* (sunflower sea star) recovery following Sea Star Wasting Disease (SSWD).

## Key Findings

### Without Intervention: Extinction

With parameters calibrated to observed SSWD dynamics:
- **99% mortality** (expert input: "closer to 99%, maybe higher")
- **0.5% initial resistance** (extremely rare before outbreak)
- **6 billion pre-outbreak population**

The model shows complete extinction within 4 years of disease onset:

| Year | Population | Notes |
|------|------------|-------|
| 10 | 100% | Disease arrives |
| 11 | 9% | Rapid collapse |
| 12 | 1.5% | Near-extinction |
| 14 | 0% | Extinct |

### Outplanting Prevents Extinction

Even minimal outplanting (100 individuals/year to 10 sites) prevents extinction, but population remains tiny:

| Scenario | Extinction | Final Pop | Resistance |
|----------|------------|-----------|------------|
| No intervention | 100% | 0% | — |
| 100/yr × 10 sites | 0% | 0.01% | evolving |
| 1k/yr × 50 sites | 0% | 1% | 99% |
| 5k/yr × 100 sites | 0% | 10% | 99% |
| 10k/yr × 100 sites | 0% | 20% | 98% |

### Natural Selection Drives Recovery

The most striking result: **resistance allele frequency evolves from 0.5% to 98%** under selection pressure. Survivors carry resistance alleles, and outplanting from wild survivors is more effective than artificial selection:

| Outplant Source | Final Pop | Final Resistance |
|-----------------|-----------|------------------|
| Wild survivors | 20% | 98% |
| Enhanced (50% resistance) | 17% | 52% |

**Why wild works better:** Surviving wild individuals have already been selected for resistance. Adding individuals with only 50% resistance *dilutes* the high-resistance gene pool.

### Policy Implications

1. **Outplanting is necessary** — without intervention, extinction is certain
2. **Scale matters** — 10,000/year to 100 sites needed for meaningful recovery (20% of baseline in 100 years)
3. **Source population matters** — captive breeding from survivors preferred over artificial selection
4. **Time horizon is long** — recovery to 20% takes ~100 years even with intensive intervention
5. **Resistance evolution is the mechanism** — the population that survives is genetically different from the pre-disease population

## Model Architecture

### Network Structure

- **1000 sites** (scaled representation of range)
- **10,000 individuals per site** (scaled from 6 billion total)
- **Connectivity matrix** controls larval dispersal between sites
- **Ratio-based outputs** — all results as fraction of baseline

### Connectivity Types

| Type | Description | Extinction Risk |
|------|-------------|-----------------|
| Uniform | All sites equally connected | Highest (dilution) |
| Stepping-stone | Only neighbors connected | Lowest |
| Distance-decay | Decreasing with distance | Moderate |
| Asymmetric | Directional (current-like) | Variable |
| Hub-network | Few highly-connected sites | Moderate |
| Modular | Clusters | Low |

**Key finding:** Stepping-stone connectivity is most protective. Uniform dispersal causes "dilution death spiral" where survivors spread too thin to rebuild.

### Disease Dynamics

Calibrated to observed SSWD:
- **Simultaneous multi-site outbreak** (10% of sites initially)
- **90% spread rate per year** (matches rapid coastal spread)
- **99% mortality** when infected
- **Endemic persistence** at 25% prevalence
- **Minimal recovery** (1% per year)

### Genetics

- **10 resistance loci** (simplified from ~51 in real genome)
- **Additive effects** — each resistance allele reduces mortality
- **Resistance effect:** 40% mortality reduction for homozygotes (99% → 59%)
- **Selection strength:** ~0.30 per locus (extremely strong)
- **Gene flow** via larval dispersal

## Installation

```bash
git clone https://github.com/TidepoolCurrent/pycnopodia-model.git
cd pycnopodia-model
pip install -e .
```

## Usage

### Basic Simulation

```python
from pycnopodia.network import NetworkSimulation, NetworkConfig

# Default parameters (99% mortality, 0.5% resistance)
config = NetworkConfig(n_sites=100, n_years=50)
sim = NetworkSimulation(config, seed=42)
result = sim.run()

print(f"Final population: {result.final_n_ratio:.1%} of baseline")
print(f"Bottleneck: {result.min_n_ratio:.1%}")
print(f"Extinct: {result.extinct}")
```

### Outplanting Scenario

```python
config = NetworkConfig(
    n_sites=1000,
    n_years=100,
    outplanting_n=5000,           # 5000 individuals per event
    outplanting_sites=list(range(0, 1000, 10)),  # 100 sites
    outplanting_start=15,          # Start 5 years after disease
    outplant_resistance_mode='wild',  # Use wild survivors
)

sim = NetworkSimulation(config, seed=42)
result = sim.run()
```

### Compare Connectivity Scenarios

```python
from pycnopodia.network import run_scenario, ConnectivityType

result = run_scenario(
    connectivity_type=ConnectivityType.STEPPING_STONE,
    self_recruitment=0.5,
    n_replicates=10,
    n_years=50,
)
print(f"Extinction probability: {result['extinction_probability']:.0%}")
```

## Parameters

### Demographics

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_sites` | 1000 | Number of sites in network |
| `n_per_site` | 10000 | Initial individuals per site |
| `survival_adult` | 0.95 | Annual adult survival (disease-free) |
| `survival_juvenile` | 0.60 | Annual juvenile survival |
| `maturation_years` | 3 | Years to reproductive maturity |

### Disease

| Parameter | Default | Description |
|-----------|---------|-------------|
| `disease_onset_year` | 10 | Year disease first appears |
| `disease_mortality` | 0.99 | Mortality rate when infected |
| `disease_spread_rate` | 0.90 | Fraction of sites infected per year |
| `disease_endemic_prevalence` | 0.25 | Long-term endemic level |

### Genetics

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_loci` | 10 | Number of resistance loci |
| `initial_resistance_freq` | 0.005 | Starting allele frequency |
| `resistance_effect` | 0.40 | Mortality reduction per allele |

### Intervention

| Parameter | Default | Description |
|-----------|---------|-------------|
| `outplanting_n` | 0 | Individuals per outplanting event |
| `outplanting_sites` | [] | Sites receiving outplants |
| `outplanting_start` | 15 | Year outplanting begins |
| `outplant_resistance_mode` | 'wild' | 'wild' or 'enhanced' |
| `outplant_enhanced_resistance` | 0.50 | Resistance freq if enhanced |

## Tests

```bash
pytest tests/ -v
```

## Citation

Model developed for sunflower sea star (*Pycnopodia helianthoides*) conservation research at Friday Harbor Labs, University of Washington.

## License

MIT

## References

- Harvell, C.D., et al. (2019). Disease epidemic and a marine heat wave are associated with the continental-scale collapse of a pivotal predator. *Science Advances*.
- Schiebelhut, L.M., et al. (2024). Genomic signatures of disease resistance in sea stars.
- IUCN Red List: *Pycnopodia helianthoides* — Critically Endangered (2021)
