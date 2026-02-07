# Pycnopodia Population Model

A computational model for *Pycnopodia helianthoides* (sunflower sea star) population dynamics, simulating recovery trajectories following Sea Star Wasting Disease (SSWD) and evaluating captive breeding intervention strategies.

## Overview

This model simulates Pycnopodia population dynamics with two complementary approaches:

1. **Network Metapopulation Model** — A 1000-site spatial model tracking population counts, genetics, and disease spread across connected habitat patches
2. **Individual-Based Model** — Detailed simulation tracking individual genomes, ages, and reproductive success for single-population analyses

Both models incorporate:
- **Sweepstakes Reproductive Success (SRS)** — Only ~8% of adults breed successfully each year
- **Polygenic SSWD Resistance** — Multiple loci contributing to disease survival
- **Broadcast Spawner Allee Effects** — Fertilization failure at low densities
- **Captive Breeding Interventions** — Outplanting scenarios with wild-caught or enhanced-resistance stock

## Model Architecture

### Network Model (`network.py`)

![Model Overview](figures/model_overview.png)

**Structure:**
- 1,000 sites representing habitat patches along the Pacific coast
- Connectivity via larval dispersal (stepping stone, distance decay, or other topologies)
- Per-site tracking of population size, resistance allele frequency, and disease prevalence

**Disease Dynamics:**
- Rapid multi-site outbreak (simulating 2013-2014 SSWD epidemic)
- 99% base mortality for susceptible individuals
- Resistance reduces mortality to ~30% for homozygous resistant genotypes
- Endemic persistence at 25% prevalence post-epidemic

**Genetics:**
- 3% initial resistance allele frequency (rare before outbreak)
- 70% resistance effect (RR mortality = 99% × 0.30 = 30%)
- 95% biological ceiling on resistance frequency
- Selection + drift with gene flow via larval dispersal

### Individual-Based Model (`population.py`, `simulation.py`)

**Structure:**
- Individual tracking with explicit genomes (20 loci)
- Age-structured population with stage-specific survival
- Mendelian inheritance during reproduction

**Key Components:**
- `Population` — Stores genomes, ages, sexes for all individuals
- `Simulation` — Orchestrates annual cycle (mortality → disease → reproduction → aging)
- `Broodstock` — Captive breeding population for outplanting

## Key Findings

### Population Trajectory

![Population Trajectory](figures/population_trajectory.png)

Without intervention, populations crash to ~1% of baseline but don't go extinct due to disease-free refugia (~5% of sites). With outplanting intervention, recovery to 10-15% of baseline is achievable over 100 years.

### Sweepstakes Reproductive Success

![SRS Effect](figures/srs_effect.png)

Only ~8% of adults successfully reproduce each year, dramatically reducing effective population size (Ne ≈ 10% of census N). This accelerates genetic drift and loss of rare alleles during bottlenecks.

### Disease Dynamics

![Disease Dynamics](figures/disease_dynamics.png)

SSWD follows an epidemic → endemic pattern:
- Peak prevalence: 90% at outbreak
- Decay rate: ~10% per year
- Endemic level: 20% long-term

### Allee Effects

![Allee Effect](figures/allee_effect.png)

Broadcast spawners require sufficient adult density for successful fertilization. Below ~30 adults, reproduction fails due to gamete dilution.

### Genetic Diversity

![Genetic Diversity](figures/genetic_diversity.png)

Outplanting helps maintain genetic diversity (H/H₀) by:
1. Providing demographic rescue during bottleneck
2. Introducing genetic variation from captive stock
3. Enabling continued selection for resistance

## Interactive Dashboard

**NEW!** 🌊 Explore the Pacific Coast model interactively with the visualization dashboard:

```bash
# Generate data
python3 run_pacific_coast.py --export-json

# Start dashboard server
cd dashboard
python3 -m http.server 8090

# Open http://localhost:8090 in your browser
```

**Features:**
- **Network graph** with sites positioned by lat/lon, colored by region
- **Time slider** to animate through 100 years (2010-2110)
- **Live statistics** showing population, disease, and resistance by region
- **Interactive nodes** showing population size, disease prevalence (red overlay), and resistance (green rings)
- **Playback controls** with adjustable speed (0.5× to 4×)

See [`dashboard/README.md`](dashboard/README.md) for full documentation.

## Outplanting Strategies

| Strategy | Final Population | Final Resistance |
|----------|------------------|------------------|
| No intervention | 0.3-1% | 5-10% |
| 5k/yr wild-caught | 9-12% | 90-95% |
| 5k/yr enhanced (50%) | 3-5% | 55-60% |
| 5k/yr enhanced (95%) | 11-15% | 94-95% |

**Key insight:** Wild-caught outplants (from natural survivors) often outperform moderately-enhanced stock because natural selection has already enriched resistance alleles. The 95% enhanced breeding scenario (near biological ceiling) provides the best outcomes, but wild-caught approaches similar performance due to natural selection pressure.

## Installation

```bash
git clone https://github.com/TidepoolCurrent/pycnopodia-model.git
cd pycnopodia-model
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

## Usage

### Network Model

```bash
# Single simulation with detailed output
python3 run_network.py --single --n-years 100

# Compare connectivity scenarios
python3 run_network.py --compare --n-replicates 10

# Specific scenario
python3 run_network.py --scenario stepping --self-recruitment 0.7
```

### Individual-Based Model

```bash
# Baseline scenario
python3 run_simulation.py --n-replicates 50 --n-years 80

# Without intervention
python3 run_simulation.py --no-outplanting

# Different preset
python3 run_simulation.py --preset high_outplanting
```

### Python API

```python
from pycnopodia.network import NetworkSimulation, NetworkConfig

# Configure simulation
config = NetworkConfig(
    n_sites=1000,
    n_years=100,
    disease_mortality=0.99,
    outplanting_n=5000,
    outplanting_sites=list(range(0, 1000, 10)),  # Every 10th site
    outplanting_start=15,
    outplant_resistance_mode='enhanced',
    outplant_enhanced_resistance=0.95,
)

# Run simulation
sim = NetworkSimulation(config, seed=42)
result = sim.run()

# Access results (ratio-based)
print(f"Final population: {result.final_n_ratio:.2%} of baseline")
print(f"Bottleneck: {result.min_n_ratio:.2%} at year {result.bottleneck_year}")
print(f"Final resistance: {result.states[-1].mean_resistance_freq:.1%}")
```

## Key Parameters

### Disease

| Parameter | Default | Description |
|-----------|---------|-------------|
| `disease_mortality` | 0.99 | Base mortality when infected |
| `disease_onset_year` | 10 | Year disease enters population |
| `disease_spread_rate` | 0.90 | Fraction of sites infected per year |
| `disease_endemic_prevalence` | 0.25 | Long-term prevalence |
| `refugia_fraction` | 0.05 | Fraction of sites that escape disease |

### Genetics (Polygenic Resistance)

The network model now supports **polygenic resistance** with multiple loci contributing to disease resistance, matching the individual-based model architecture.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `initial_resistance_freq` | 0.03 | Starting resistance allele frequency per locus |
| `resistance_effect` | 0.70 | Total mortality reduction when all loci at 100% |
| `max_resistance_freq` | 0.95 | Biological ceiling per locus |
| `n_loci` | 10 | Number of resistance loci to track |
| `locus_effects` | None | Per-locus effect sizes (if None, sampled from gamma) |
| `locus_effect_shape` | 2.0 | Gamma shape parameter for sampling effect sizes |

**Polygenic architecture:**
- Resistance allele frequencies tracked at each locus per site: shape `(n_sites, n_loci)`
- Selection acts independently on each locus based on its effect size
- Genetic drift affects each locus independently
- Gene flow transfers alleles at each locus during larval dispersal
- Overall resistance = sum of (frequency × effect) across loci, normalized to [0, 1]

### Intervention

| Parameter | Default | Description |
|-----------|---------|-------------|
| `outplanting_n` | 100 | Individuals outplanted per event |
| `outplanting_start` | 15 | Year outplanting begins |
| `outplanting_interval` | 1 | Years between outplanting events |
| `outplant_resistance_mode` | "wild" | "wild" or "enhanced" |
| `outplant_enhanced_resistance` | 0.50 | Resistance freq if enhanced mode |

### Connectivity (Network Model)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `connectivity_type` | STEPPING_STONE | Network topology |
| `self_recruitment` | 0.50 | Fraction of larvae staying at natal site |
| `asymmetry` | 0.0 | Directional bias (0=symmetric, 1=one-way) |

## Output Interpretation

All outputs are expressed as **ratios** relative to baseline for cross-scenario comparability:

- **N/N₀** — Population size relative to initial (e.g., 0.15 = 15% of baseline)
- **H/H₀** — Heterozygosity retention (genetic diversity)
- **Ne/N** — Effective vs. census population size (captures SRS effect)

## Project Structure

```
pycnopodia-model/
├── pycnopodia/
│   ├── __init__.py          # Package exports
│   ├── config.py             # Configuration and presets
│   ├── network.py            # Network metapopulation model
│   ├── population.py         # Individual-based population class
│   ├── simulation.py         # Simulation orchestration
│   ├── reproduction.py       # SRS and Mendelian inheritance
│   ├── disease.py            # SSWD dynamics
│   ├── intervention.py       # Captive breeding and outplanting
│   ├── environment.py        # Environmental stochasticity
│   ├── spatial.py            # Spatial structure (optional)
│   ├── srs.py                # Hierarchical SRS model
│   ├── size_structured.py    # Size-based demographics
│   └── real_data.py          # PNW environmental data
├── tests/                    # Test suite (53 tests)
├── scripts/
│   ├── generate_figures.py   # README figures
│   ├── quick_analysis.py     # Quick visualization
│   └── exhaustive_analysis.py # Full parameter sweeps
├── figures/                  # Generated figures
├── run_network.py            # Network model CLI
├── run_simulation.py         # Individual-based model CLI
└── README.md
```

## Testing

```bash
python3 -m pytest tests/ -v
```

All 53 tests verify model components including connectivity matrices, disease dynamics, genetic inheritance, and intervention mechanics.

## Scientific Background

### Pycnopodia helianthoides

The sunflower sea star is the largest sea star in the world (arm span up to 1m) and a keystone predator in Pacific kelp forest ecosystems. SSWD caused >90% population decline across its range from 2013-2017.

### SSWD (Sea Star Wasting Disease)

A wasting syndrome causing lesions, arm loss, and death within days. Associated with densovirus and environmental stressors. Mortality rates exceeded 90% in most populations.

### Conservation Context

This model supports recovery planning by:
1. Estimating natural recovery timelines
2. Evaluating captive breeding strategies
3. Identifying key parameters for intervention success
4. Quantifying genetic diversity retention under different scenarios

## References

- Harvell, C.D. et al. (2019). Disease epidemic and a marine heat wave are associated with the continental-scale collapse of a pivotal predator. *Science Advances*
- Schiebelhut, L.M. et al. (2024). Polygenic basis of SSWD resistance. *bioRxiv*
- Hedgecock, D. & Pudovkin, A.I. (2011). Sweepstakes reproductive success in highly fecund marine fish and shellfish. *Bull. Mar. Sci.*

## License

MIT

## Citation

If you use this model, please cite:

```
Pycnopodia Population Model (2026). 
https://github.com/TidepoolCurrent/pycnopodia-model
```
