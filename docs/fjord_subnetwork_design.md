# Fjord Subnetwork Design

## Biology (from WLWeertman)

- Stars ARE recovering in AK, BC, parts of Salish Sea (post-2020)
- Recovery seeded from fjord refugia where disease bacteria couldn't penetrate
- Bacteria transmission: ~2 weeks exposure → symptoms → transmission
- Fjord geography physically blocks bacterial spread (not waterborne long-distance)
- Local larval retention in fjords = self-sustaining populations
- Hamilton et al. (through 2020) shows the CRASH phase; recovery is post-2020

## Model Architecture

### Current Problem
- Fjord regions are single well-mixed pools
- Either disease reaches all sites (crash too severe, no recovery)
- Or refugia protect too many sites (crash not severe enough)
- Can't get BOTH 96% decline AND subsequent recovery

### Proposed Fix: Subnetworks

Within each fjord region (BC Fjords, SE AK North), create **isolated subnetworks**:

```
Region: SE AK North (20 sites)
├── Fjord A (3 sites) — isolated pocket, high self-recruitment
├── Fjord B (2 sites) — isolated pocket  
├── Fjord C (3 sites) — isolated pocket
└── Outer sites (12 sites) — connected, disease-accessible
```

**Disease network:** Subnetworks nearly disconnected from each other and from outer sites.
Disease spreads freely among outer sites but barely enters fjord pockets.

**Larval network:** Fjord subnetworks have HIGH self-recruitment (0.8+) but also export 
some larvae to outer sites. Recovery happens when:
1. Outer sites crash from disease
2. Disease subsides (density-dependent decline)
3. Fjord populations export larvae to recolonize outer sites

### Key Parameters

- `n_subnetworks`: Number of isolated fjord pockets per region (3-5)
- `subnetwork_size`: Sites per pocket (2-4)
- `disease_isolation`: Factor reducing disease transmission into subnetwork (0.01-0.05)
- `larval_self_recruitment`: Within-subnetwork retention (0.7-0.9)
- `larval_export`: Fraction of larvae escaping subnetwork (0.1-0.3)

### Expected Behavior

Year 0-10: Pre-disease equilibrium
Year 10-17: Disease crashes outer sites (90%+ decline)
            Fjord pockets barely touched (isolated from disease)
            Regional average: ~96% decline (most sites are outer)
Year 17-30: Disease subsides (density-dependent)
            Fjord pockets begin exporting larvae
            Slow recovery in outer sites near fjords
Year 30-100: Gradual recovery, asymptoting at lower-than-original density
             (climate warming + periodic disease recurrence)

### Calibration Targets

Year 17 (2020): ~96% decline SE AK, ~92% Salish, ~88% BC
Year 25 (2028): Early signs of recovery in SE AK, BC fjord-adjacent sites  
Year 50 (2043): Partial recovery in northern regions, south still depressed
