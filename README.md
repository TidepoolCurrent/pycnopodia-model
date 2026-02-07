# Pycnopodia Population Model

Network-based metapopulation model for *Pycnopodia helianthoides* (sunflower sea star) recovery following Sea Star Wasting Disease (SSWD).

## Key Findings

### Calibrated to Reality

Parameters based on expert input (W. Weertman, Friday Harbor Labs):
- **99% mortality** when infected
- **~1-3% initial resistance** (rare before outbreak)  
- **95% maximum resistance** (biological ceiling)
- **5% refugia** (sites that escape disease)
- **~6 billion pre-outbreak population**

### Without Intervention: Severe Crash (Not Extinction)

With refugia, the population crashes but doesn't go fully extinct:

| Year | Population | Notes |
|------|------------|-------|
| 10 | 100% | Disease arrives |
| 11 | ~10% | Rapid collapse |
| 12 | ~1% | Bottleneck |
| 50 | 0.01% | Survivors in refugia |

**Key insight:** Complete extinction in 4 years was unrealistic. Real populations have refugia (cooler waters, isolation) that allow persistence.

### Outplanting Scenarios (100 years)

| Scenario | Final Pop | Resistance |
|----------|-----------|------------|
| No intervention | 0.01% | 10% |
| 1k/yr × 50 sites (wild) | 1.5% | 84% |
| 5k/yr × 100 sites (wild) | **15%** | 88% |
| 5k/yr × 100 (enhanced 50%) | 10% | 52% |
| 5k/yr × 100 (enhanced 80%) | 13% | 82% |
| 5k/yr × 100 (enhanced 95%) | **14%** | 95% |

### Key Insights

1. **Refugia prevent extinction** — Without intervention, ~0.01% survive in disease-free sites
2. **Wild-caught outplants work well** — Natural selection drives resistance to ~88%
3. **Enhanced breeding needs 80%+ resistance to match wild** — 50% enhanced is *worse* than wild
4. **95% enhanced ≈ wild performance** — At the biological ceiling, both strategies converge
5. **15% recovery is the practical ceiling** at 5k/yr scale over 100 years

### Policy Implications

1. **Outplanting accelerates recovery** but doesn't restore pre-outbreak population
2. **Source matters:** Wild survivors (high natural resistance) > low/moderate enhanced
3. **If breeding, push resistance close to 95%** — half-measures underperform wild
4. **Protect refugia** — They're the source of survivors
5. **100-year time horizon** for meaningful (10-15%) recovery

## Model Architecture

### Network Structure

- **1000 sites** (scaled representation of range)
- **10,000 individuals per site** (scaled from 6 billion total)
- **5% refugia** — sites that never get infected
- **95% resistance ceiling** — biological constraint

### Disease Dynamics

- **Simultaneous multi-site outbreak** (10% of non-refugia sites initially)
- **90% spread rate per year** to remaining sites
- **99% mortality** when infected (modulated by resistance)
- **Endemic persistence** at 25% prevalence

### Genetics

- **3% initial resistance** — rare but enough for some survivors
- **70% resistance effect** — RR mortality drops from 99% to 30%
- **95% ceiling** — resistance cannot exceed this
- **Selection + drift** — realistic evolutionary dynamics

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

config = NetworkConfig(n_sites=1000, n_years=100)
sim = NetworkSimulation(config, seed=42)
result = sim.run()

print(f"Final population: {result.final_n_ratio:.2%}")
print(f"Final resistance: {result.states[-1].mean_resistance_freq:.1%}")
```

### Outplanting with Enhanced Resistance

```python
config = NetworkConfig(
    n_sites=1000,
    n_years=100,
    outplanting_n=5000,
    outplanting_sites=list(range(0, 1000, 10)),  # 100 sites
    outplanting_start=15,
    outplant_resistance_mode='enhanced',
    outplant_enhanced_resistance=0.95,  # 95% resistance
)
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `disease_mortality` | 0.99 | Base mortality when infected |
| `initial_resistance_freq` | 0.03 | Starting resistance (3%) |
| `resistance_effect` | 0.70 | Mortality reduction for RR |
| `max_resistance_freq` | 0.95 | Biological ceiling |
| `refugia_fraction` | 0.05 | Fraction of sites that escape disease |
| `outplant_enhanced_resistance` | 0.95 | Resistance if enhanced mode |

## Tests

```bash
pytest tests/ -v
```

## Citation

Model developed for sunflower sea star conservation research at Friday Harbor Labs, University of Washington.

## License

MIT
