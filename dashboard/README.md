# Pycnopodia Population Model - Interactive Dashboard

An interactive web visualization of the Pacific Coast Pycnopodia helianthoides (sunflower sea star) population dynamics model, showing disease spread, population decline, and resistance evolution from 2010-2110.

## Features

- **Network Graph**: Sites positioned by latitude/longitude, colored by region
- **Population Dynamics**: Node size represents population at each site
- **Disease Spread**: Red overlay shows disease prevalence (darker = more disease)
- **Resistance Evolution**: Green rings show resistance allele frequency
- **Connectivity**: Edge opacity shows larval connectivity strength between sites
- **Time Animation**: Scrub through 100 years of simulation data
- **Playback Controls**: Play/pause, step forward/back, speed control (0.5× to 4×)
- **Statistics Panel**: Real-time regional and overall population statistics

## Quick Start

### 1. Generate Data

Run the simulation with the `--export-json` flag:

```bash
python3 run_pacific_coast.py --export-json
```

Or, to only export data without generating figures:

```bash
python3 run_pacific_coast.py --json-only
```

This creates `dashboard/data.json` with simulation results.

### 2. Start Local Server

**Option A: Python HTTP Server (Recommended)**

```bash
cd dashboard
python3 -m http.server 8090
```

**Option B: Using the provided server script**

```bash
cd dashboard
python3 serve.py
```

### 3. Open in Browser

Navigate to: **http://localhost:8090**

## Usage

### Controls

- **Play/Pause**: Start/stop animation (or press `Space`)
- **Step Forward/Back**: Advance/reverse by one year (or use arrow keys `←` `→`)
- **Year Slider**: Scrub to any year in the simulation
- **Speed Control**: Adjust playback speed (0.5×, 1×, 2×, 4×)

### Visualization

- **Node Size**: Larger nodes = higher population
- **Node Color**: Base color = region; red blend = disease prevalence
- **Green Rings**: Width and opacity = resistance frequency
- **Edges**: Opacity = connectivity strength (larval dispersal)
- **Hover**: Mouse over nodes for detailed statistics

### Key Patterns to Observe

1. **Disease Onset (Year 10)**: Disease begins in Southern California (warmest region)
2. **Northward Spread**: Disease propagates north along the coast following connectivity
3. **Regional Differences**: 
   - Southern regions (S. CA, N. CA) near-total collapse
   - BC Fjords maintain ~50% population (refugia)
   - SE Alaska North (coldest) highest survival
4. **Resistance Evolution**: Green rings grow stronger in surviving populations
5. **Temperature Gradient**: Cooler northern sites fare better

## Data Structure

The `data.json` file contains:

```json
{
  "sites": [
    {
      "idx": 0,
      "lat": 32.1,
      "lon": -125.3,
      "region": "s_california",
      "is_refugia": false,
      "temperature": 15.2
    },
    ...
  ],
  "connectivity": [[...], ...],  // 330×330 sparse matrix
  "years": [0, 1, 2, ..., 99],
  "timeseries": {
    "populations": [[...], ...],  // [site][year]
    "disease": [[...], ...],
    "resistance": [[...], ...]
  },
  "regions": {
    "s_california": {
      "name": "Southern CA / Channel Islands",
      "short_name": "S. CA",
      "color": "#d62728"
    },
    ...
  }
}
```

## Technology Stack

- **D3.js v7**: Network visualization and data binding
- **Vanilla JavaScript**: No build step required
- **CSS Grid**: Responsive layout
- **Dark Theme**: Optimized for scientific data visualization

## Browser Compatibility

Modern browsers with ES6 support:
- Chrome/Edge 90+
- Firefox 88+
- Safari 14+

## Performance

- **Data Size**: ~4 MB JSON (330 sites × 100 years)
- **Initial Load**: < 2 seconds on modern hardware
- **Animation**: 60 FPS with D3 transitions

## Customization

Edit `viz.js` to modify:
- Node size scaling (`sizeScale` domain/range)
- Color mapping (disease blending intensity)
- Connectivity edge threshold (`threshold` variable)
- Animation speed intervals

## Troubleshooting

**Issue**: Dashboard doesn't load / shows blank page
- Check browser console for errors
- Ensure `data.json` exists in the dashboard directory
- Verify you're accessing via HTTP server (not `file://`)

**Issue**: Animation is slow
- Reduce playback speed
- Close other browser tabs
- Try Chrome/Edge for best performance

**Issue**: Nodes too small/large
- Adjust `sizeScale` range in `viz.js` lines 140-142

## Citation

Model developed for Friday Harbor Labs Pycnopodia population dynamics research.
Dashboard visualization by OpenClaw AI (February 2026).

## License

Part of the pycnopodia-model repository. See parent directory for license information.
