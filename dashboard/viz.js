// Pycnopodia Population Model Dashboard
// Interactive D3.js visualization

let data = null;
let currentYearIndex = 0;
let isPlaying = false;
let playSpeed = 1;
let playInterval = null;

const svg = d3.select('#network-svg');
const tooltip = d3.select('#tooltip');

// Load and initialize
d3.json('data.json').then(loadedData => {
    data = loadedData;
    console.log('Data loaded:', data.sites.length, 'sites,', data.years.length, 'years');
    
    initVisualization();
    setupControls();
    updateVisualization(0);
});

function initVisualization() {
    const container = document.getElementById('visualization');
    const width = container.clientWidth;
    const height = container.clientHeight;
    
    svg.attr('width', width).attr('height', height);
    
    // Create groups for layers
    svg.append('g').attr('class', 'links');
    svg.append('g').attr('class', 'nodes');
    svg.append('g').attr('class', 'resistance-rings');
    
    // Calculate positions based on lat/lon
    // Map latitude (32-60) to y-axis (height down to 50)
    // Map longitude variation to x-axis
    const latExtent = d3.extent(data.sites, d => d.lat);
    const lonExtent = d3.extent(data.sites, d => d.lon);
    
    const xScale = d3.scaleLinear()
        .domain(lonExtent)
        .range([width * 0.2, width * 0.8]);
    
    const yScale = d3.scaleLinear()
        .domain(latExtent)
        .range([height - 50, 100]);
    
    // Add positions to sites
    data.sites.forEach(site => {
        site.x = xScale(site.lon);
        site.y = yScale(site.lat);
    });
    
    // Draw connectivity edges (sparse)
    const links = [];
    const threshold = 0.01; // Only show strong connections
    
    for (let i = 0; i < data.connectivity.length; i++) {
        for (let j = 0; j < data.connectivity[i].length; j++) {
            const weight = data.connectivity[i][j];
            if (weight > threshold && i !== j) {
                links.push({
                    source: data.sites[i],
                    target: data.sites[j],
                    weight: weight
                });
            }
        }
    }
    
    console.log('Drawing', links.length, 'connectivity edges');
    
    // Draw links
    svg.select('.links')
        .selectAll('line')
        .data(links)
        .join('line')
        .attr('class', 'link')
        .attr('x1', d => d.source.x)
        .attr('y1', d => d.source.y)
        .attr('x2', d => d.target.x)
        .attr('y2', d => d.target.y)
        .attr('stroke-opacity', d => Math.min(d.weight * 2, 0.6));
    
    // Draw nodes (will update sizes/colors with updateVisualization)
    svg.select('.nodes')
        .selectAll('circle')
        .data(data.sites)
        .join('circle')
        .attr('class', 'node')
        .attr('cx', d => d.x)
        .attr('cy', d => d.y)
        .attr('r', 5)
        .attr('fill', d => data.regions[d.region].color)
        .on('mouseover', showTooltip)
        .on('mousemove', moveTooltip)
        .on('mouseout', hideTooltip);
    
    // Draw resistance rings
    svg.select('.resistance-rings')
        .selectAll('circle')
        .data(data.sites)
        .join('circle')
        .attr('class', 'resistance-ring')
        .attr('cx', d => d.x)
        .attr('cy', d => d.y)
        .attr('r', 0);
}

function updateVisualization(yearIndex) {
    currentYearIndex = Math.max(0, Math.min(yearIndex, data.years.length - 1));
    const year = data.years[currentYearIndex];
    
    // Update year display
    document.getElementById('current-year').textContent = year;
    document.getElementById('year-display').textContent = year;
    document.getElementById('year-slider').value = currentYearIndex;
    
    // Get current state
    const populations = data.timeseries.populations.map(siteSeries => siteSeries[currentYearIndex]);
    const diseases = data.timeseries.disease.map(siteSeries => siteSeries[currentYearIndex]);
    const resistances = data.timeseries.resistance.map(siteSeries => siteSeries[currentYearIndex]);
    
    // Calculate initial populations for scaling
    const initialPopulations = data.timeseries.populations.map(siteSeries => siteSeries[0]);
    const maxInitialPop = d3.max(initialPopulations);
    
    // Scale for node size (0-20px radius)
    const sizeScale = d3.scaleSqrt()
        .domain([0, maxInitialPop])
        .range([0, 20]);
    
    // Update nodes
    svg.select('.nodes')
        .selectAll('circle')
        .data(data.sites)
        .transition()
        .duration(300)
        .attr('r', (d, i) => Math.max(2, sizeScale(populations[i])))
        .attr('fill', (d, i) => {
            const baseColor = d3.color(data.regions[d.region].color);
            const disease = diseases[i];
            
            // Blend with red based on disease prevalence
            if (disease > 0.01) {
                const red = d3.rgb(255, 0, 0);
                return d3.interpolateRgb(baseColor, red)(disease * 0.8);
            }
            return baseColor;
        })
        .attr('opacity', (d, i) => populations[i] > 0 ? 1 : 0.2);
    
    // Update resistance rings
    svg.select('.resistance-rings')
        .selectAll('circle')
        .data(data.sites)
        .transition()
        .duration(300)
        .attr('r', (d, i) => {
            const pop = populations[i];
            const resist = resistances[i];
            if (pop > 0 && resist > 0.1) {
                return Math.max(2, sizeScale(pop)) + 3;
            }
            return 0;
        })
        .attr('stroke-width', (d, i) => {
            const resist = resistances[i];
            return Math.max(1, resist * 3);
        })
        .attr('stroke-opacity', (d, i) => resistances[i]);
    
    // Update statistics
    updateStatistics(populations, diseases, resistances, initialPopulations);
}

function updateStatistics(populations, diseases, resistances, initialPopulations) {
    // Overall stats
    const totalPop = d3.sum(populations);
    const initialTotalPop = d3.sum(initialPopulations);
    const decline = (1 - totalPop / initialTotalPop) * 100;
    const meanResistance = d3.mean(resistances.filter((r, i) => populations[i] > 0));
    const meanDisease = d3.mean(diseases);
    
    document.getElementById('stat-total-pop').textContent = Math.round(totalPop).toLocaleString();
    document.getElementById('stat-decline').textContent = decline.toFixed(1) + '%';
    document.getElementById('stat-resistance').textContent = (meanResistance * 100).toFixed(1) + '%';
    document.getElementById('stat-disease').textContent = (meanDisease * 100).toFixed(1) + '%';
    
    // Regional stats
    const regionalStats = {};
    
    data.sites.forEach((site, i) => {
        const region = site.region;
        if (!regionalStats[region]) {
            regionalStats[region] = {
                population: 0,
                initialPopulation: 0,
                disease: [],
                resistance: [],
                color: data.regions[region].color,
                name: data.regions[region].short_name
            };
        }
        regionalStats[region].population += populations[i];
        regionalStats[region].initialPopulation += initialPopulations[i];
        if (populations[i] > 0) {
            regionalStats[region].disease.push(diseases[i]);
            regionalStats[region].resistance.push(resistances[i]);
        }
    });
    
    // Render regional stats
    const regionalStatsDiv = document.getElementById('regional-stats');
    regionalStatsDiv.innerHTML = '';
    
    Object.entries(regionalStats).forEach(([regionId, stats]) => {
        const decline = (1 - stats.population / stats.initialPopulation) * 100;
        const meanDisease = stats.disease.length > 0 ? d3.mean(stats.disease) : 0;
        const meanResistance = stats.resistance.length > 0 ? d3.mean(stats.resistance) : 0;
        
        const div = document.createElement('div');
        div.className = 'region-item';
        div.style.borderLeftColor = stats.color;
        
        div.innerHTML = `
            <div class="region-name">${stats.name}</div>
            <div class="region-stats">
                <span>Pop: ${Math.round(stats.population)}</span>
                <span>Decline: ${decline.toFixed(0)}%</span>
                <span>Res: ${(meanResistance * 100).toFixed(0)}%</span>
            </div>
        `;
        
        regionalStatsDiv.appendChild(div);
    });
}

function showTooltip(event, d) {
    const i = d.idx;
    const year = currentYearIndex;
    
    const pop = data.timeseries.populations[i][year];
    const disease = data.timeseries.disease[i][year];
    const resistance = data.timeseries.resistance[i][year];
    
    tooltip.html(`
        <div><span class="tooltip-label">Site:</span> <span class="tooltip-value">#${i}</span></div>
        <div><span class="tooltip-label">Region:</span> <span class="tooltip-value">${data.regions[d.region].short_name}</span></div>
        <div><span class="tooltip-label">Latitude:</span> <span class="tooltip-value">${d.lat.toFixed(2)}°N</span></div>
        <div><span class="tooltip-label">Temperature:</span> <span class="tooltip-value">${d.temperature.toFixed(1)}°C</span></div>
        <div><span class="tooltip-label">Population:</span> <span class="tooltip-value">${Math.round(pop)}</span></div>
        <div><span class="tooltip-label">Disease:</span> <span class="tooltip-value">${(disease * 100).toFixed(1)}%</span></div>
        <div><span class="tooltip-label">Resistance:</span> <span class="tooltip-value">${(resistance * 100).toFixed(1)}%</span></div>
        ${d.is_refugia ? '<div style="color: #4caf50; margin-top: 5px;">🌊 Refugia Site</div>' : ''}
    `);
    
    tooltip.classed('visible', true);
}

function moveTooltip(event) {
    tooltip
        .style('left', (event.pageX + 15) + 'px')
        .style('top', (event.pageY - 15) + 'px');
}

function hideTooltip() {
    tooltip.classed('visible', false);
}

function setupControls() {
    const playBtn = document.getElementById('play-btn');
    const stepBackBtn = document.getElementById('step-back');
    const stepForwardBtn = document.getElementById('step-forward');
    const yearSlider = document.getElementById('year-slider');
    const speedButtons = document.querySelectorAll('.speed-btn');
    
    yearSlider.max = data.years.length - 1;
    
    playBtn.addEventListener('click', togglePlay);
    stepBackBtn.addEventListener('click', () => stepYear(-1));
    stepForwardBtn.addEventListener('click', () => stepYear(1));
    yearSlider.addEventListener('input', (e) => {
        updateVisualization(parseInt(e.target.value));
    });
    
    speedButtons.forEach(btn => {
        btn.addEventListener('click', () => {
            speedButtons.forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            playSpeed = parseFloat(btn.dataset.speed);
            
            // Restart play if currently playing
            if (isPlaying) {
                stopPlay();
                startPlay();
            }
        });
    });
    
    // Keyboard shortcuts
    document.addEventListener('keydown', (e) => {
        if (e.code === 'Space') {
            e.preventDefault();
            togglePlay();
        } else if (e.code === 'ArrowLeft') {
            stepYear(-1);
        } else if (e.code === 'ArrowRight') {
            stepYear(1);
        }
    });
}

function togglePlay() {
    if (isPlaying) {
        stopPlay();
    } else {
        startPlay();
    }
}

function startPlay() {
    isPlaying = true;
    document.getElementById('play-btn').textContent = '⏸ Pause';
    document.getElementById('app').classList.add('playing');
    
    const interval = 1000 / playSpeed; // Base: 1 year per second
    
    playInterval = setInterval(() => {
        if (currentYearIndex >= data.years.length - 1) {
            stopPlay();
            return;
        }
        updateVisualization(currentYearIndex + 1);
    }, interval);
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
    const newIndex = currentYearIndex + delta;
    if (newIndex >= 0 && newIndex < data.years.length) {
        updateVisualization(newIndex);
    }
}

// Handle window resize
window.addEventListener('resize', () => {
    if (data) {
        initVisualization();
        updateVisualization(currentYearIndex);
    }
});
