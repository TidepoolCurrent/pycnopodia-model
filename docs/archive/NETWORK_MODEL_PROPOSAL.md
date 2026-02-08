# Network Model Proposal: Larval Connectivity

## The Problem

The current model treats populations as locally-recruiting. But echinoderm larvae don't work this way:

- **Planktonic duration:** Pycnopodia larvae spend 2-8 weeks in the water column
- **Passive dispersal:** Larvae go where currents take them
- **Non-local recruitment:** A site's recruits come from OTHER sites, not local spawning

This means a site's population depends on what's happening elsewhere in the network.

---

## Why This Matters for Recovery

### Current model assumption:
```
Site A spawns → Larvae → Recruits at Site A
```

### Reality:
```
Site A spawns → Larvae enter currents → Recruits at Sites B, C, D
Site A's recruits come from Sites X, Y, Z upstream
```

### Implications:

1. **A collapsed site can recover from external supply** - if upstream sites are healthy
2. **A healthy site can collapse from lack of supply** - if upstream sites fail
3. **Recovery depends on network structure** - not just local dynamics
4. **Outplanting location matters** - source sites vs sink sites

---

## The Connectivity Matrix

Define **C[i,j]** = probability that a larva produced at site i settles at site j

```
         To:
         AK    BC    WA    OR    CA
From:
  AK   [ 0.30  0.40  0.20  0.10  0.00 ]  # Alaska: mostly disperses south
  BC   [ 0.05  0.30  0.35  0.20  0.10 ]  # BC: strong southward flow
  WA   [ 0.00  0.05  0.30  0.40  0.25 ]  # WA: continues south
  OR   [ 0.00  0.00  0.10  0.40  0.50 ]  # OR: to California
  CA   [ 0.00  0.00  0.00  0.15  0.85 ]  # CA: retention (but isolated)
```

**Key features:**
- **Asymmetric:** Larvae flow south with California Current (summer)
- **Some self-recruitment:** Diagonal elements (eddies, bays)
- **Rows don't sum to 1:** Some larvae lost to unsuitable habitat
- **Northern sites are sources:** Alaska exports to the whole network
- **Southern sites are sinks:** California receives but exports less

---

## How Recruitment Changes

### Current model:
```python
# All recruits from local spawning
recruits = local_larvae * larval_survival
```

### Network model:
```python
# Recruits come from all connected sites
for site_j in sites:
    recruits[site_j] = sum(
        larvae_produced[site_i] * connectivity[site_i, site_j]
        for site_i in sites
    )
```

The number of larvae at site j depends on:
1. How many larvae each upstream site produced
2. The connectivity from each upstream site to j
3. Survival during transport (could add distance-dependent mortality)

---

## New Dynamics This Creates

### 1. Source-Sink Dynamics

Some sites are **net exporters** (sources):
- High local reproduction
- Currents carry larvae away
- Local recruitment < export

Some sites are **net importers** (sinks):
- May have low local reproduction
- Receive larvae from upstream
- Persist due to external supply

**For Pycnopodia:** Alaska/BC are likely sources. California is likely a sink.

### 2. Rescue Effects

A collapsed local population can recover if:
- Upstream populations are healthy
- Connectivity allows sufficient larval supply
- Local conditions support settlement and survival

**This is good news:** Don't need to restore every site directly.

### 3. Cascade Failures

If source populations collapse:
- Downstream sites lose larval supply
- Even healthy local adults can't compensate
- Collapse propagates through network

**This is bad news:** Losing Alaska could doom the whole coast.

### 4. Strategic Outplanting

Where you outplant matters:
- **Source sites:** Outplanted individuals supply the whole network
- **Sink sites:** Outplanted individuals only help locally

**Implication:** Focus outplanting on upstream source sites (Alaska, BC) for maximum network benefit.

---

## What We'd Need to Add

### Data Requirements:

1. **Oceanographic connectivity** 
   - Particle tracking models (ROMS, etc.)
   - Larval duration parameterization
   - Seasonal variation in currents

2. **Site locations and habitat**
   - Where are the actual populations?
   - What's the suitable habitat area at each?
   - Historical abundance estimates?

### Model Changes:

1. **Multi-site population tracking**
   - Already have metapopulation structure in spatial.py
   - Need to replace migration with larval connectivity

2. **Connectivity matrix**
   - Input from oceanographic model, or
   - Parameterized stepping-stone with asymmetry

3. **Settlement dynamics**
   - Larvae need suitable substrate
   - Density-dependent settlement? (avoid crowded sites)

4. **Site-specific parameters**
   - Temperature (already have for thermal refugia)
   - Habitat quality
   - Local SSWD severity

---

## Key Questions to Explore

1. **How sensitive is recovery to connectivity structure?**
   - Compare: stepping-stone vs. island model vs. realistic PNW currents
   - What if Alaska is lost first vs. California first?

2. **What's the minimum viable network?**
   - How many sites must persist for network recovery?
   - Are there keystone sites?

3. **Where should outplanting focus?**
   - Compare: uniform vs. source-focused vs. sink-focused
   - What's the optimal allocation?

4. **How does disease spread through the network?**
   - Currently assume independent outbreaks
   - Could model larval transmission
   - What if disease follows the same currents as larvae?

5. **Temporal matching:**
   - Spawning timing varies by latitude
   - Current patterns are seasonal
   - Does timing mismatch limit connectivity?

---

## Proposed Implementation

### Phase 1: Simple Network (days)
- Replace migration with connectivity matrix
- Use stylized stepping-stone pattern
- Test basic source-sink dynamics

### Phase 2: PNW-Specific (weeks)
- Get real oceanographic connectivity data
- 5 sites: Alaska, BC, Washington, Oregon, California
- Seasonal connectivity variation

### Phase 3: Full Realism (months)
- Higher resolution sites
- Couple to climate projections
- Larval behavior (vertical migration)

---

## Biological Grounding

### Pycnopodia larval biology:
- **Spawning:** Late winter to early summer (varies with latitude)
- **Larval stages:** Bipinnaria → Brachiolaria
- **Planktonic duration:** 2-8 weeks (temperature dependent)
- **Settlement cues:** Coralline algae, biofilms
- **Competency window:** Larvae must settle within time window or die

### PNW oceanography:
- **California Current:** Southward flow, especially summer
- **Davidson Current:** Northward, winter (weaker)
- **Coastal upwelling:** Seasonal, affects nearshore retention
- **Eddies:** Can create retention zones

### Empirical connectivity studies:
- No direct studies on Pycnopodia dispersal
- General echinoderm studies show 10-100km typical dispersal
- But long-distance events possible (100s of km)
- Genetic studies could inform connectivity

---

## Conclusion

The current model underestimates both the vulnerability (cascade failures from source loss) and the resilience (rescue from connected populations) of Pycnopodia populations.

A network model would:
1. Better capture the biology of broadcast spawning
2. Identify strategic sites for outplanting
3. Reveal network-level vulnerabilities
4. Guide monitoring priorities

The basic metapopulation structure exists in the code. The key change is replacing "local recruitment + adult migration" with "network-based larval recruitment."

---

*Prepared for critique. Key uncertainties: connectivity matrix values, larval survival during transport, temporal dynamics.*
