// Pycnopodia Geo Model Dashboard - Leaflet Edition
// Interactive visualization with real map — 56 Pacific Coast sites

let data = null;
let currentYearIndex = 0;
let isPlaying = false;
let playSpeed = 1;
let playInterval = null;
let map = null;
let siteMarkers = [];
let resistanceRings = [];
let connectivityLines = [];
let showConnectivity = true;

const REGION_COLORS = {
    se_alaska_north: '#1a5276',
    se_alaska_south: '#2e86c1',
    bc_fjords: '#148f77',
    bc_outer: '#45b39d',
    salish_sea: '#7d3c98',
    wa_or_outer: '#e67e22',
    n_california: '#e67e22',
    c_california: '#d35400',
    s_california: '#a04000',
};

const REGION_NAMES = {
    se_alaska_north: 'SE Alaska N (fjords)',
    se_alaska_south: 'SE Alaska S',
    bc_fjords: 'BC Fjords',
    bc_outer: 'BC Outer Coast',
    salish_sea: 'Salish Sea',
    wa_or_outer: 'WA/OR Coast',
    n_california: 'N. California',
    c_california: 'C. California',
    s_california: 'S. California',
};

// Load data and initialize
d3.json('data.json').then(loadedData => {
    data = loadedData;
    console.log('Geo model loaded:', data.n_sites, 'sites,', (data.total_steps || data.n_years), 'years');
    initMap();
    initVisualization();
    setupControls();
    updateVisualization(0);
});

function initMap() {
    // Initialize Leaflet map centered on Pacific Coast
    map = L.map('map', {
        center: [48, -126],
        zoom: 5,
        zoomControl: true,
        minZoom: 4,
        maxZoom: 10,
    });

    // CartoDB dark matter tiles for dark theme
    L.tileLayer('https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png', {
        attribution: '© OpenStreetMap © CartoDB',
        subdomains: 'abcd',
        maxZoom: 19,
    }).addTo(map);

    // Alternative: OpenStreetMap standard tiles (lighter)
    // L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    //     attribution: '© OpenStreetMap contributors',
    //     maxZoom: 19,
    // }).addTo(map);
}

function initVisualization() {
    // Calculate initial populations for size scaling
    const initialPops = data.site_timeseries.map(s => s.population[0]);
    const maxInitPop = d3.max(initialPops);
    
    // Store size scale function for later use
    window.sizeScale = d3.scaleSqrt().domain([0, maxInitPop]).range([3, 20]);

    // Draw connectivity edges (behind nodes)
    const links = data.larval_connectivity.filter(l => l.from !== l.to && l.weight > 0.005);
    
    links.forEach(link => {
        const from = data.sites[link.from];
        const to = data.sites[link.to];
        
        const polyline = L.polyline(
            [[from.lat, from.lon], [to.lat, to.lon]],
            {
                color: '#2a3f5f',
                weight: 1,
                opacity: Math.min(link.weight * 3, 0.4),
                className: 'connectivity-edge'
            }
        ).addTo(map);
        
        connectivityLines.push(polyline);
    });

    // Create markers for each site
    data.sites.forEach((site, i) => {
        // Resistance ring (larger circle behind)
        const ring = L.circleMarker([site.lat, site.lon], {
            radius: 0,
            color: '#4caf50',
            fillColor: 'transparent',
            fillOpacity: 0,
            weight: 2,
            className: 'resistance-ring'
        }).addTo(map);
        resistanceRings.push(ring);

        // Main node marker
        const marker = L.circleMarker([site.lat, site.lon], {
            radius: 5,
            fillColor: REGION_COLORS[site.region] || '#888',
            color: '#0a0e27',
            weight: 1.5,
            fillOpacity: 1,
            className: 'pulsing-node'
        }).addTo(map);

        // Popup with site info
        marker.bindPopup(() => buildPopupContent(site, i), {
            closeButton: true,
            className: 'site-popup'
        });

        // Hover effect
        marker.on('mouseover', function() {
            this.setStyle({ weight: 2.5, color: '#4fc3f7' });
        });
        marker.on('mouseout', function() {
            this.setStyle({ weight: 1.5, color: '#0a0e27' });
        });

        siteMarkers.push(marker);
    });
}

function buildPopupContent(site, i) {
    const pop = data.site_timeseries[i].population[currentYearIndex];
    const disease = data.site_timeseries[i].prevalence[currentYearIndex];
    const resistance = data.site_timeseries[i].resistance[currentYearIndex];
    const initPop = data.site_timeseries[i].population[0];

    return `
        <div>
            <div><span class="tooltip-label">Site:</span> <span class="tooltip-value">${site.name}</span></div>
            <div><span class="tooltip-label">Region:</span> <span class="tooltip-value">${REGION_NAMES[site.region]}</span></div>
            <div><span class="tooltip-label">Type:</span> <span class="tooltip-value">${site.site_type}${site.has_freshwater_lens ? ' 🧊' : ''}${site.sill_depth ? ' (sill: '+site.sill_depth+'m)' : ''}</span></div>
            <div><span class="tooltip-label">Lat/Lon:</span> <span class="tooltip-value">${site.lat.toFixed(2)}°N, ${Math.abs(site.lon).toFixed(2)}°W</span></div>
            <div><span class="tooltip-label">Base Temp:</span> <span class="tooltip-value">${site.base_temp.toFixed(1)}°C</span></div>
            <hr style="border-color:#2a3f5f;margin:4px 0">
            <div><span class="tooltip-label">Population:</span> <span class="tooltip-value">${Math.round(pop).toLocaleString()} / ${Math.round(initPop).toLocaleString()}</span></div>
            <div><span class="tooltip-label">Disease:</span> <span class="tooltip-value" style="color:${disease>0.1?'#ff6b6b':'#8896ab'}">${(disease * 100).toFixed(1)}%</span></div>
            <div><span class="tooltip-label">Resistance:</span> <span class="tooltip-value" style="color:${resistance>0.05?'#4caf50':'#8896ab'}">${(resistance * 100).toFixed(1)}%</span></div>
        </div>
    `;
}

function updateVisualization(yearIndex) {
    currentYearIndex = Math.max(0, Math.min(yearIndex, (data.total_steps || data.n_years) - 1));

    document.getElementById('current-year').textContent = data.step_labels ? data.step_labels[currentYearIndex] : (2010 + currentYearIndex);
    document.getElementById('year-display').textContent = data.step_labels ? data.step_labels[currentYearIndex] : (2010 + currentYearIndex);
    document.getElementById('year-slider').value = currentYearIndex;

    const populations = data.site_timeseries.map(s => s.population[currentYearIndex]);
    const diseases = data.site_timeseries.map(s => s.prevalence[currentYearIndex]);
    const resistances = data.site_timeseries.map(s => s.resistance[currentYearIndex]);
    const initialPops = data.site_timeseries.map(s => s.population[0]);

    // Update each site marker
    data.sites.forEach((site, i) => {
        const marker = siteMarkers[i];
        const ring = resistanceRings[i];
        const pop = populations[i];
        const disease = diseases[i];
        const resistance = resistances[i];
        
        // Calculate node size based on population
        const radius = Math.max(3, window.sizeScale(pop));
        
        // Calculate color with disease blend
        let color = REGION_COLORS[site.region] || '#888';
        if (disease > 0.01) {
            const baseColor = d3.color(color);
            // Use bright magenta/pink for disease so it's visible on red (California) regions
            const diseaseColor = d3.rgb(255, 0, 200);
            color = d3.interpolateRgb(baseColor, diseaseColor)(disease * 0.8);
        }
        
        // Update marker appearance
        marker.setStyle({
            radius: radius,
            fillColor: color,
            fillOpacity: pop > 0 ? 1 : 0.15
        });

        // Update resistance ring
        if (pop > 0 && resistance > 0.05) {
            ring.setStyle({
                radius: radius + 4,
                weight: Math.max(1, resistance * 5),
                opacity: Math.min(resistance * 3, 1)
            });
        } else {
            ring.setStyle({
                radius: 0,
                opacity: 0
            });
        }

        // Update popup if open
        if (marker.isPopupOpen()) {
            marker.setPopupContent(buildPopupContent(site, i));
        }
    });

    updateStatistics(populations, diseases, resistances, initialPops);
}

function updateStatistics(populations, diseases, resistances, initialPops) {
    const totalPop = d3.sum(populations);
    const totalInit = d3.sum(initialPops);
    const decline = (1 - totalPop / totalInit) * 100;
    const liveSites = populations.filter(p => p > 0).length;
    const meanRes = d3.mean(resistances.filter((r, i) => populations[i] > 0)) || 0;
    const meanDis = d3.mean(diseases) || 0;

    document.getElementById('stat-total-pop').textContent = Math.round(totalPop).toLocaleString();
    document.getElementById('stat-decline').textContent = decline.toFixed(1) + '%';
    document.getElementById('stat-resistance').textContent = (meanRes * 100).toFixed(1) + '%';
    document.getElementById('stat-disease').textContent = (meanDis * 100).toFixed(1) + '%';

    // Regional breakdown
    const regStats = {};
    data.sites.forEach((site, i) => {
        const r = site.region;
        if (!regStats[r]) regStats[r] = { pop: 0, init: 0, dis: [], res: [] };
        regStats[r].pop += populations[i];
        regStats[r].init += initialPops[i];
        if (populations[i] > 0) {
            regStats[r].dis.push(diseases[i]);
            regStats[r].res.push(resistances[i]);
        }
    });

    const div = document.getElementById('regional-stats');
    div.innerHTML = '';
    const order = ['se_alaska_north','se_alaska_south','bc_fjords','bc_outer','salish_sea','wa_or_outer','n_california','c_california','s_california'];
    order.forEach(r => {
        const s = regStats[r];
        if (!s) return;
        const dec = ((1 - s.pop / s.init) * 100);
        const res = s.res.length ? d3.mean(s.res) : 0;
        const el = document.createElement('div');
        el.className = 'region-item';
        el.style.borderLeftColor = REGION_COLORS[r];
        el.innerHTML = `
            <div class="region-name">${REGION_NAMES[r]}</div>
            <div class="region-stats">
                <span>Pop: ${Math.round(s.pop).toLocaleString()}</span>
                <span>Dec: ${dec.toFixed(0)}%</span>
                <span>Res: ${(res * 100).toFixed(0)}%</span>
            </div>`;
        div.appendChild(el);
    });
}

function setupControls() {
    document.getElementById('year-slider').max = (data.total_steps || data.n_years) - 1;
    document.getElementById('play-btn').addEventListener('click', togglePlay);
    document.getElementById('step-back').addEventListener('click', () => stepYear(-1));
    document.getElementById('step-forward').addEventListener('click', () => stepYear(1));
    document.getElementById('year-slider').addEventListener('input', e => updateVisualization(parseInt(e.target.value)));

    // Connectivity toggle
    document.getElementById('toggle-connectivity').addEventListener('change', e => {
        showConnectivity = e.target.checked;
        connectivityLines.forEach(line => {
            if (showConnectivity) {
                line.addTo(map);
            } else {
                map.removeLayer(line);
            }
        });
    });

    // Speed controls
    document.querySelectorAll('.speed-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            document.querySelectorAll('.speed-btn').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            playSpeed = parseFloat(btn.dataset.speed);
            if (isPlaying) { stopPlay(); startPlay(); }
        });
    });

    // Keyboard shortcuts
    document.addEventListener('keydown', e => {
        if (e.code === 'Space') { 
            e.preventDefault(); 
            togglePlay(); 
        }
        else if (e.code === 'ArrowLeft') {
            e.preventDefault();
            stepYear(-1);
        }
        else if (e.code === 'ArrowRight') {
            e.preventDefault();
            stepYear(1);
        }
    });
}

function togglePlay() { 
    isPlaying ? stopPlay() : startPlay(); 
}

function startPlay() {
    isPlaying = true;
    document.getElementById('play-btn').textContent = '⏸ Pause';
    document.getElementById('app').classList.add('playing');
    playInterval = setInterval(() => {
        if (currentYearIndex >= (data.total_steps || data.n_years) - 1) { 
            stopPlay(); 
            return; 
        }
        updateVisualization(currentYearIndex + 1);
    }, 800 / playSpeed);
}

function stopPlay() {
    isPlaying = false;
    document.getElementById('play-btn').textContent = '▶ Play';
    document.getElementById('app').classList.remove('playing');
    if (playInterval) { 
        clearInterval(playInterval); 
        playInterval = null; 
    }
}

function stepYear(delta) {
    const n = currentYearIndex + delta;
    if (n >= 0 && n < (data.total_steps || data.n_years)) {
        updateVisualization(n);
    }
}

// Handle window resize
window.addEventListener('resize', () => {
    if (map) {
        map.invalidateSize();
    }
});
