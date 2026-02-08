// Pycnopodia Geo Model Dashboard
// Interactive D3.js visualization — 56 real Pacific Coast sites

let data = null;
let currentYearIndex = 0;
let isPlaying = false;
let playSpeed = 1;
let playInterval = null;

const REGION_COLORS = {
    se_alaska_north: '#1a5276',
    se_alaska_south: '#2e86c1',
    bc_fjords: '#148f77',
    bc_outer: '#45b39d',
    salish_sea: '#7d3c98',
    wa_or_outer: '#e67e22',
    n_california: '#e74c3c',
    c_california: '#c0392b',
    s_california: '#922b21',
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

const svg = d3.select('#network-svg');
const tooltip = d3.select('#tooltip');

d3.json('data.json').then(loadedData => {
    data = loadedData;
    console.log('Geo model loaded:', data.n_sites, 'sites,', data.n_years, 'years');
    initVisualization();
    setupControls();
    updateVisualization(0);
});

function initVisualization() {
    const container = document.getElementById('visualization');
    const width = container.clientWidth;
    const height = container.clientHeight;

    svg.attr('width', width).attr('height', height);
    svg.selectAll('g').remove();

    svg.append('g').attr('class', 'links');
    svg.append('g').attr('class', 'nodes');
    svg.append('g').attr('class', 'resistance-rings');
    svg.append('g').attr('class', 'labels');

    const latExtent = d3.extent(data.sites, d => d.lat);
    const lonExtent = d3.extent(data.sites, d => d.lon);

    const xScale = d3.scaleLinear().domain(lonExtent).range([width * 0.15, width * 0.85]);
    const yScale = d3.scaleLinear().domain(latExtent).range([height - 40, 60]);

    data.sites.forEach((site, i) => {
        site.x = xScale(site.lon);
        site.y = yScale(site.lat);
        site.idx = i;
    });

    // Draw larval connectivity edges
    const links = data.larval_connectivity.filter(l => l.from !== l.to && l.weight > 0.005);
    svg.select('.links')
        .selectAll('line')
        .data(links)
        .join('line')
        .attr('class', 'link')
        .attr('x1', d => data.sites[d.from].x)
        .attr('y1', d => data.sites[d.from].y)
        .attr('x2', d => data.sites[d.to].x)
        .attr('y2', d => data.sites[d.to].y)
        .attr('stroke-opacity', d => Math.min(d.weight * 3, 0.4));

    // Draw nodes
    svg.select('.nodes')
        .selectAll('circle')
        .data(data.sites)
        .join('circle')
        .attr('class', 'node')
        .attr('cx', d => d.x)
        .attr('cy', d => d.y)
        .attr('r', 5)
        .attr('fill', d => REGION_COLORS[d.region] || '#888')
        .on('mouseover', showTooltip)
        .on('mousemove', moveTooltip)
        .on('mouseout', hideTooltip);

    // Resistance rings
    svg.select('.resistance-rings')
        .selectAll('circle')
        .data(data.sites)
        .join('circle')
        .attr('class', 'resistance-ring')
        .attr('cx', d => d.x)
        .attr('cy', d => d.y)
        .attr('r', 0);

    // Site name labels for fjords
    svg.select('.labels')
        .selectAll('text')
        .data(data.sites.filter(s => s.site_type === 'fjord'))
        .join('text')
        .attr('x', d => d.x + 12)
        .attr('y', d => d.y + 3)
        .text(d => d.name)
        .attr('font-size', '8px')
        .attr('fill', '#556')
        .attr('opacity', 0.6);
}

function updateVisualization(yearIndex) {
    currentYearIndex = Math.max(0, Math.min(yearIndex, data.n_years - 1));

    document.getElementById('current-year').textContent = 2010 + currentYearIndex;
    document.getElementById('year-display').textContent = 2010 + currentYearIndex;
    document.getElementById('year-slider').value = currentYearIndex;

    const populations = data.site_timeseries.map(s => s.population[currentYearIndex]);
    const diseases = data.site_timeseries.map(s => s.prevalence[currentYearIndex]);
    const resistances = data.site_timeseries.map(s => s.resistance[currentYearIndex]);
    const initialPops = data.site_timeseries.map(s => s.population[0]);
    const maxInitPop = d3.max(initialPops);

    const sizeScale = d3.scaleSqrt().domain([0, maxInitPop]).range([2, 22]);

    svg.select('.nodes')
        .selectAll('circle')
        .data(data.sites)
        .transition().duration(250)
        .attr('r', (d, i) => Math.max(2, sizeScale(populations[i])))
        .attr('fill', (d, i) => {
            const base = d3.color(REGION_COLORS[d.region] || '#888');
            if (diseases[i] > 0.01) {
                return d3.interpolateRgb(base, d3.rgb(255, 30, 30))(diseases[i] * 0.8);
            }
            return base;
        })
        .attr('opacity', (d, i) => populations[i] > 0 ? 1 : 0.15);

    svg.select('.resistance-rings')
        .selectAll('circle')
        .data(data.sites)
        .transition().duration(250)
        .attr('r', (d, i) => {
            if (populations[i] > 0 && resistances[i] > 0.05) {
                return Math.max(2, sizeScale(populations[i])) + 3;
            }
            return 0;
        })
        .attr('stroke-width', (d, i) => Math.max(1, resistances[i] * 5))
        .attr('stroke-opacity', (d, i) => Math.min(resistances[i] * 3, 1));

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

function showTooltip(event, d) {
    const i = d.idx;
    const pop = data.site_timeseries[i].population[currentYearIndex];
    const disease = data.site_timeseries[i].prevalence[currentYearIndex];
    const resistance = data.site_timeseries[i].resistance[currentYearIndex];
    const initPop = data.site_timeseries[i].population[0];

    tooltip.html(`
        <div><span class="tooltip-label">Site:</span> <span class="tooltip-value">${d.name}</span></div>
        <div><span class="tooltip-label">Region:</span> <span class="tooltip-value">${REGION_NAMES[d.region]}</span></div>
        <div><span class="tooltip-label">Type:</span> <span class="tooltip-value">${d.site_type}${d.has_freshwater_lens ? ' 🧊' : ''}${d.sill_depth ? ' (sill: '+d.sill_depth+'m)' : ''}</span></div>
        <div><span class="tooltip-label">Lat/Lon:</span> <span class="tooltip-value">${d.lat.toFixed(2)}°N, ${d.lon.toFixed(2)}°W</span></div>
        <div><span class="tooltip-label">Base Temp:</span> <span class="tooltip-value">${d.base_temp.toFixed(1)}°C</span></div>
        <hr style="border-color:#2a3f5f;margin:4px 0">
        <div><span class="tooltip-label">Population:</span> <span class="tooltip-value">${Math.round(pop).toLocaleString()} / ${Math.round(initPop).toLocaleString()}</span></div>
        <div><span class="tooltip-label">Disease:</span> <span class="tooltip-value" style="color:${disease>0.1?'#ff6b6b':'#8896ab'}">${(disease * 100).toFixed(1)}%</span></div>
        <div><span class="tooltip-label">Resistance:</span> <span class="tooltip-value" style="color:${resistance>0.05?'#4caf50':'#8896ab'}">${(resistance * 100).toFixed(1)}%</span></div>
    `);
    tooltip.classed('visible', true);
}

function moveTooltip(event) {
    tooltip.style('left', (event.pageX + 15) + 'px').style('top', (event.pageY - 15) + 'px');
}

function hideTooltip() {
    tooltip.classed('visible', false);
}

function setupControls() {
    document.getElementById('year-slider').max = data.n_years - 1;
    document.getElementById('play-btn').addEventListener('click', togglePlay);
    document.getElementById('step-back').addEventListener('click', () => stepYear(-1));
    document.getElementById('step-forward').addEventListener('click', () => stepYear(1));
    document.getElementById('year-slider').addEventListener('input', e => updateVisualization(parseInt(e.target.value)));

    document.querySelectorAll('.speed-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            document.querySelectorAll('.speed-btn').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            playSpeed = parseFloat(btn.dataset.speed);
            if (isPlaying) { stopPlay(); startPlay(); }
        });
    });

    document.addEventListener('keydown', e => {
        if (e.code === 'Space') { e.preventDefault(); togglePlay(); }
        else if (e.code === 'ArrowLeft') stepYear(-1);
        else if (e.code === 'ArrowRight') stepYear(1);
    });
}

function togglePlay() { isPlaying ? stopPlay() : startPlay(); }

function startPlay() {
    isPlaying = true;
    document.getElementById('play-btn').textContent = '⏸ Pause';
    playInterval = setInterval(() => {
        if (currentYearIndex >= data.n_years - 1) { stopPlay(); return; }
        updateVisualization(currentYearIndex + 1);
    }, 800 / playSpeed);
}

function stopPlay() {
    isPlaying = false;
    document.getElementById('play-btn').textContent = '▶ Play';
    if (playInterval) { clearInterval(playInterval); playInterval = null; }
}

function stepYear(delta) {
    const n = currentYearIndex + delta;
    if (n >= 0 && n < data.n_years) updateVisualization(n);
}

window.addEventListener('resize', () => { if (data) { initVisualization(); updateVisualization(currentYearIndex); } });
