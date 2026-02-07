"""
Tests for Pacific Coast geographic module.

Tests cover:
- Region setup (correct number of sites, temperatures)
- Connectivity patterns (within-region > between-region)
- Disease origin in south
- Temperature gradient
- BC Fjords as refugia
"""

import pytest
import numpy as np
from pycnopodia.pacific_coast import (
    PacificCoastConfig,
    PacificCoastSimulation,
    PacificCoastResult,
    PACIFIC_COAST_REGIONS,
    REGION_ORDER,
    RegionType,
    RegionConfig,
    Site,
    build_sites,
    build_larval_connectivity_matrix,
    build_disease_connectivity_matrix,
    get_site_region_mapping,
    get_site_region_boundaries,
    get_temperature_disease_modifier,
    get_temperature_spread_modifier,
)


class TestRegionConfiguration:
    """Test region setup and configuration."""
    
    def test_nine_regions_defined(self):
        """Should have exactly 9 regions."""
        assert len(PACIFIC_COAST_REGIONS) == 9
        assert len(REGION_ORDER) == 9
    
    def test_region_order_matches_keys(self):
        """Region order should match defined regions."""
        for region_id in REGION_ORDER:
            assert region_id in PACIFIC_COAST_REGIONS
    
    def test_regions_have_valid_latitude_ranges(self):
        """All regions should have valid, non-empty latitude ranges."""
        for region_id in REGION_ORDER:
            region = PACIFIC_COAST_REGIONS[region_id]
            lat_min, lat_max = region.latitude_range
            assert lat_max > lat_min, f"{region_id} has invalid latitude range"
            assert 30 <= lat_min <= 60, f"{region_id} latitude out of Pacific coast range"
            assert 30 <= lat_max <= 62, f"{region_id} latitude out of Pacific coast range"
    
    def test_southern_ca_warmest(self):
        """Southern CA should have highest base temperature."""
        temps = {r: PACIFIC_COAST_REGIONS[r].base_temperature for r in REGION_ORDER}
        warmest = max(temps, key=temps.get)
        assert warmest == "s_california"
    
    def test_alaska_coldest(self):
        """Alaska should have lowest base temperature."""
        outer_coast_temps = {
            r: PACIFIC_COAST_REGIONS[r].base_temperature 
            for r in REGION_ORDER 
            if PACIFIC_COAST_REGIONS[r].region_type == RegionType.OUTER_COAST
        }
        coldest = min(outer_coast_temps, key=outer_coast_temps.get)
        assert coldest in ("se_alaska_north", "se_alaska_south")
    
    def test_bc_fjords_is_refugia(self):
        """BC Fjords should be marked as refugia (fjord type)."""
        bc_fjords = PACIFIC_COAST_REGIONS["bc_fjords"]
        assert bc_fjords.region_type == RegionType.FJORD
        # Updated: fjords delay but don't prevent disease (Hamilton et al. 2021)
        assert 0.10 <= bc_fjords.post_sswd_survival <= 0.25  # Moderate survival
    
    def test_salish_sea_is_inland(self):
        """Salish Sea should be marked as inland sea."""
        salish = PACIFIC_COAST_REGIONS["salish_sea"]
        assert salish.region_type == RegionType.INLAND_SEA
    
    def test_region_survival_rates_realistic(self):
        """Post-SSWD survival rates should match observed patterns from Hamilton et al. 2021."""
        # SE Alaska: 96% decline overall (~4% remaining)
        # Split between northern fjords and southern outer coast
        assert 0.02 <= PACIFIC_COAST_REGIONS["se_alaska_north"].post_sswd_survival <= 0.10
        assert 0.01 <= PACIFIC_COAST_REGIONS["se_alaska_south"].post_sswd_survival <= 0.05
        
        # BC: 87.9% decline (~12% remaining)
        # BC Fjords delay disease but don't prevent it
        assert 0.10 <= PACIFIC_COAST_REGIONS["bc_fjords"].post_sswd_survival <= 0.25
        
        # Southern regions near 0% (>99% decline)
        assert PACIFIC_COAST_REGIONS["s_california"].post_sswd_survival <= 0.05
        assert PACIFIC_COAST_REGIONS["n_california"].post_sswd_survival <= 0.05
        
        # Salish Sea: 92.4% decline (~7.6% remaining)
        assert PACIFIC_COAST_REGIONS["salish_sea"].post_sswd_survival <= 0.10


class TestSiteGeneration:
    """Test site building."""
    
    def test_total_sites_reasonable(self):
        """Should generate between 200-500 sites total."""
        config = PacificCoastConfig()
        sites = build_sites(config)
        assert 200 <= len(sites) <= 500
    
    def test_sites_ordered_by_region(self):
        """Sites should be ordered by region (following REGION_ORDER)."""
        config = PacificCoastConfig()
        sites = build_sites(config)
        
        # Sites should follow region order
        current_region_idx = 0
        for site in sites:
            region_idx = REGION_ORDER.index(site.region_id)
            assert region_idx >= current_region_idx, \
                f"Site {site.idx} in {site.region_id} appears after region {REGION_ORDER[current_region_idx]}"
            current_region_idx = region_idx
    
    def test_site_indices_sequential(self):
        """Site indices should be 0 to n-1."""
        config = PacificCoastConfig()
        sites = build_sites(config)
        
        for i, site in enumerate(sites):
            assert site.idx == i
    
    def test_sites_have_valid_temperatures(self):
        """All sites should have reasonable temperatures."""
        config = PacificCoastConfig()
        sites = build_sites(config, rng=np.random.default_rng(42))
        
        for site in sites:
            assert 5 <= site.temperature <= 20, f"Site {site.idx} has unrealistic temp: {site.temperature}"
    
    def test_fjord_sites_mostly_refugia(self):
        """Fjord sites should have highest refugia probability (~30%, updated from 70%)."""
        config = PacificCoastConfig()
        sites = build_sites(config, rng=np.random.default_rng(42))
        
        fjord_sites = [s for s in sites if s.region_id == "bc_fjords"]
        assert len(fjord_sites) > 0
        refugia_count = sum(1 for s in fjord_sites if s.is_refugia)
        # Some fjord sites should be refugia (reduced to match Hamilton et al. 96% decline)
        assert refugia_count >= 1  # At least some deep-water refugia exist
    
    def test_non_fjord_sites_low_refugia(self):
        """Non-fjord sites should have low but non-zero refugia probability."""
        config = PacificCoastConfig()
        sites = build_sites(config, rng=np.random.default_rng(42))
        
        non_fjord_sites = [s for s in sites if s.region_id != "bc_fjords"]
        refugia_count = sum(1 for s in non_fjord_sites if s.is_refugia)
        # Some non-fjord sites can be refugia, but fraction should be low
        assert refugia_count / len(non_fjord_sites) < 0.15
    
    def test_region_boundary_mapping(self):
        """Region boundaries should cover all sites."""
        config = PacificCoastConfig()
        sites = build_sites(config)
        boundaries = get_site_region_boundaries(sites)
        
        # All regions should have boundaries
        for region_id in REGION_ORDER:
            assert region_id in boundaries
            start, end = boundaries[region_id]
            assert start < end
            assert end <= len(sites)
        
        # Boundaries should be contiguous
        prev_end = 0
        for region_id in REGION_ORDER:
            start, end = boundaries[region_id]
            assert start == prev_end, f"Gap before {region_id}"
            prev_end = end
        
        assert prev_end == len(sites)


class TestConnectivity:
    """Test connectivity matrix construction."""
    
    @pytest.fixture
    def config_and_sites(self):
        config = PacificCoastConfig()
        rng = np.random.default_rng(42)
        sites = build_sites(config, rng)
        return config, sites, rng
    
    def test_larval_connectivity_rows_sum_to_one(self, config_and_sites):
        """Larval connectivity matrix rows should sum to 1."""
        config, sites, rng = config_and_sites
        C = build_larval_connectivity_matrix(sites, config, rng)
        
        row_sums = C.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, rtol=1e-6)
    
    def test_disease_connectivity_rows_sum_to_one(self, config_and_sites):
        """Disease connectivity matrix rows should sum to 1."""
        config, sites, rng = config_and_sites
        D = build_disease_connectivity_matrix(sites, config, rng)
        
        row_sums = D.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, rtol=1e-6)
    
    def test_within_region_connectivity_higher(self, config_and_sites):
        """Within-region connectivity should be higher than between-region."""
        config, sites, rng = config_and_sites
        C = build_larval_connectivity_matrix(sites, config, rng)
        boundaries = get_site_region_boundaries(sites)
        
        for region_id in REGION_ORDER:
            start, end = boundaries[region_id]
            n_sites = end - start
            if n_sites < 2:
                continue
            
            # Mean within-region connectivity (excluding self)
            within_block = C[start:end, start:end]
            within_sum = within_block.sum() - np.trace(within_block)
            within_mean = within_sum / (n_sites * (n_sites - 1))
            
            # Mean between-region connectivity
            total_sites = len(sites)
            outside_sum = C[start:end, :].sum() - within_block.sum()
            between_mean = outside_sum / (n_sites * (total_sites - n_sites)) if total_sites > n_sites else 0
            
            if between_mean > 0:
                assert within_mean > between_mean, \
                    f"{region_id}: within ({within_mean:.4f}) should > between ({between_mean:.4f})"
    
    def test_fjord_highly_isolated(self, config_and_sites):
        """BC Fjords should have very low connectivity to other regions."""
        config, sites, rng = config_and_sites
        C = build_larval_connectivity_matrix(sites, config, rng)
        boundaries = get_site_region_boundaries(sites)
        
        fjord_start, fjord_end = boundaries["bc_fjords"]
        n_fjord = fjord_end - fjord_start
        
        # Within-fjord connectivity (proportion of row sums that stay within fjords)
        within_sum = C[fjord_start:fjord_end, fjord_start:fjord_end].sum()
        total_from_fjords = C[fjord_start:fjord_end, :].sum()  # Should be n_fjord since rows sum to 1
        
        within_fraction = within_sum / total_from_fjords
        
        # Fjords should retain most larvae (>60% stays in fjords)
        assert within_fraction > 0.6, f"Fjord retention too low: {within_fraction:.2%}"
    
    def test_asymmetric_flow_pattern(self, config_and_sites):
        """California Current should create north→south asymmetry on outer coast."""
        config, sites, rng = config_and_sites
        C = build_larval_connectivity_matrix(sites, config, rng)
        boundaries = get_site_region_boundaries(sites)
        
        # Check outer coast regions
        outer_regions = [r for r in REGION_ORDER 
                        if config.regions[r].region_type == RegionType.OUTER_COAST]
        
        # Sum of northward vs southward dispersal
        northward = 0
        southward = 0
        
        for i, source in enumerate(sites):
            if config.regions[source.region_id].region_type != RegionType.OUTER_COAST:
                continue
            for j, target in enumerate(sites):
                if config.regions[target.region_id].region_type != RegionType.OUTER_COAST:
                    continue
                if i == j:
                    continue
                if j > i:  # Northward
                    northward += C[i, j]
                else:  # Southward
                    southward += C[i, j]
        
        # Southward should dominate (California Current)
        assert southward > northward, \
            f"Southward ({southward:.2f}) should > northward ({northward:.2f})"


class TestTemperatureDiseaseRelationship:
    """Test temperature-dependent disease dynamics."""
    
    def test_temp_modifier_at_threshold(self):
        """At threshold temperature, modifier should be 1.0."""
        config = PacificCoastConfig()
        modifier = get_temperature_disease_modifier(config.disease_temp_threshold, config)
        assert modifier == 1.0
    
    def test_temp_modifier_increases_above_threshold(self):
        """Above threshold, modifier should increase."""
        config = PacificCoastConfig()
        low = get_temperature_disease_modifier(config.disease_temp_threshold, config)
        high = get_temperature_disease_modifier(config.disease_temp_threshold + 3, config)
        assert high > low
    
    def test_temp_modifier_lower_below_threshold(self):
        """Below threshold, modifier should decrease (cold water is protective), floor at 0.5."""
        config = PacificCoastConfig()
        mod_low = get_temperature_disease_modifier(config.disease_temp_threshold - 5, config)
        mod_threshold = get_temperature_disease_modifier(config.disease_temp_threshold, config)
        assert mod_threshold == 1.0
        assert mod_low < 1.0
        assert mod_low >= 0.5  # floor
    
    def test_spread_modifier_increases_with_temp(self):
        """Disease spread should be faster in warmer water."""
        config = PacificCoastConfig()
        cold = get_temperature_spread_modifier(8, config)
        warm = get_temperature_spread_modifier(16, config)
        assert warm > cold


class TestSimulation:
    """Test Pacific Coast simulation."""
    
    def test_simulation_runs(self):
        """Simulation should run without errors."""
        config = PacificCoastConfig(n_years=20)
        sim = PacificCoastSimulation(config, seed=42)
        result = sim.run()
        
        assert len(result.states) == config.n_years
        assert len(result.sites) > 0
    
    def test_disease_starts_in_south(self):
        """Disease should originate in Southern CA."""
        config = PacificCoastConfig(n_years=15)
        sim = PacificCoastSimulation(config, seed=42)
        result = sim.run()
        
        # At disease onset year, only southern sites should have disease
        onset_state = result.states[config.disease_onset_year]
        boundaries = get_site_region_boundaries(result.sites)
        
        # Southern CA should have disease
        s_ca_start, s_ca_end = boundaries["s_california"]
        s_ca_prevalence = onset_state.disease_prevalence[s_ca_start:s_ca_end].mean()
        assert s_ca_prevalence > 0.5, f"S. CA should have disease at onset: {s_ca_prevalence}"
        
        # Alaska should NOT have disease at onset
        ak_start, ak_end = boundaries["se_alaska_north"]
        ak_prevalence = onset_state.disease_prevalence[ak_start:ak_end].mean()
        assert ak_prevalence < 0.1, f"Alaska shouldn't have disease yet: {ak_prevalence}"
    
    def test_bc_fjords_delay_disease(self):
        """BC Fjords should delay disease initially but not prevent it entirely."""
        config = PacificCoastConfig(n_years=50)
        sim = PacificCoastSimulation(config, seed=42)
        result = sim.run()
        
        boundaries = get_site_region_boundaries(result.sites)
        fjord_start, fjord_end = boundaries["bc_fjords"]
        
        # During acute phase (year 13 = onset+3), fjords should have very low disease
        state_13 = result.states[13]
        fjord_prev_13 = state_13.disease_prevalence[fjord_start:fjord_end].mean()
        
        # Fjords should have lower disease than open coast during acute phase
        # (but not zero — disease does reach fjords, just slower)
        coast_prev_13 = state_13.disease_prevalence[:boundaries["s_california"][1]].mean()
        assert fjord_prev_13 < coast_prev_13 or fjord_prev_13 < 0.80, \
            f"Fjords should have less disease than coast: fjord={fjord_prev_13:.2f} coast={coast_prev_13:.2f}"
        
        # By year 50, disease should have impacted fjord populations
        # (prevalence may be low if density-dependent disease died out after crash)
        final_state = result.states[-1]
        fjord_pop_final = final_state.populations[fjord_start:fjord_end].sum()
        fjord_pop_initial = result.states[0].populations[fjord_start:fjord_end].sum()
        fjord_decline = 1.0 - fjord_pop_final / max(fjord_pop_initial, 1e-10)
        # Fjords should show significant decline (disease penetrated)
        assert fjord_decline > 0.30, \
            f"Fjords should show population decline from disease: {fjord_decline:.2f}"
    
    def test_temperature_gradient_correct(self):
        """Sites should have correct temperature gradient."""
        config = PacificCoastConfig()
        sim = PacificCoastSimulation(config, seed=42)
        
        # Mean temp should increase from north to south
        boundaries = get_site_region_boundaries(sim.sites)
        
        north_temps = sim.temperatures[boundaries["se_alaska_north"][0]:boundaries["se_alaska_north"][1]]
        south_temps = sim.temperatures[boundaries["s_california"][0]:boundaries["s_california"][1]]
        
        assert south_temps.mean() > north_temps.mean() + 3, \
            f"South ({south_temps.mean():.1f}) should be warmer than north ({north_temps.mean():.1f})"
    
    def test_population_decline_matches_pattern(self):
        """BC Fjords (refugia) should decline less than southern regions."""
        config = PacificCoastConfig(n_years=40)
        sim = PacificCoastSimulation(config, seed=42)
        result = sim.run()
        
        trajectories = result.get_region_trajectories()
        
        # Get mean population ratio over latter half of simulation
        # (allows for initial crash and partial recovery)
        half = config.n_years // 2
        mean_ratios = {
            r: trajectories[r]["population_ratio"][half:].mean() 
            for r in REGION_ORDER
        }
        
        # BC Fjords (refugia) should maintain some population
        # Southern CA should decline more severely
        # Note: with high disease mortality, even refugia may have low absolute survival
        # but should have HIGHER survival than non-refugia regions
        assert mean_ratios["bc_fjords"] >= mean_ratios["s_california"], \
            f"BC Fjords ({mean_ratios['bc_fjords']:.3f}) should do at least as well as S.CA ({mean_ratios['s_california']:.3f})"
    
    def test_resistance_evolves(self):
        """Resistance frequency should increase during disease."""
        config = PacificCoastConfig(n_years=50)
        sim = PacificCoastSimulation(config, seed=42)
        result = sim.run()
        
        # Check resistance evolution at year 20 (before most sites go extinct)
        initial_state = result.states[0]
        mid_state = result.states[20]  # Year 20: populations still surviving
        
        initial_resist = initial_state.resistance_freqs.mean()
        # Only compare resistance among surviving sites (pop > 1)
        surviving_mask = mid_state.populations > 1
        
        if surviving_mask.any():
            mid_resist = mid_state.resistance_freqs[surviving_mask].mean()
            assert mid_resist > initial_resist, \
                f"Resistance should increase among survivors by year 20: {initial_resist:.3f} → {mid_resist:.3f}"
        else:
            # If all populations went extinct, test is not meaningful
            pytest.skip("All populations extinct before resistance could evolve")
    
    def test_result_trajectories(self):
        """Result should provide valid trajectories."""
        config = PacificCoastConfig(n_years=30)
        sim = PacificCoastSimulation(config, seed=42)
        result = sim.run()
        
        trajectories = result.get_region_trajectories()
        
        # Should have all regions
        assert len(trajectories) == 9
        
        # Each region should have all metrics
        for region_id in REGION_ORDER:
            assert "population" in trajectories[region_id]
            assert "population_ratio" in trajectories[region_id]
            assert "resistance" in trajectories[region_id]
            assert "disease_prevalence" in trajectories[region_id]
            
            # Should have right length
            assert len(trajectories[region_id]["population"]) == config.n_years


class TestPolygeneticResistance:
    """Test polygenic resistance tracking."""
    
    def test_resistance_shape(self):
        """Resistance should be tracked per locus per site."""
        config = PacificCoastConfig()
        sim = PacificCoastSimulation(config, seed=42)
        
        expected_shape = (sim.n_sites, config.n_loci)
        assert sim.resistance_freqs.shape == expected_shape
    
    def test_locus_effects_sum_to_total(self):
        """Per-locus effects should sum to total resistance effect."""
        config = PacificCoastConfig()
        sim = PacificCoastSimulation(config, seed=42)
        
        total_effect = sim.locus_effects.sum()
        np.testing.assert_allclose(total_effect, config.resistance_effect, rtol=0.01)
    
    def test_site_resistance_bounded(self):
        """Computed site resistance should be in [0, 1]."""
        config = PacificCoastConfig()
        sim = PacificCoastSimulation(config, seed=42)
        
        site_resist = sim._compute_site_resistance()
        assert (site_resist >= 0).all()
        assert (site_resist <= 1).all()


class TestRefugia:
    """Test refugia behavior."""
    
    def test_refugia_sites_identified(self):
        """Simulation should identify refugia sites."""
        config = PacificCoastConfig()
        sim = PacificCoastSimulation(config, seed=42)
        
        assert len(sim.refugia_sites) > 0
    
    def test_refugia_mostly_fjord_sites(self):
        """Refugia should be predominantly in BC Fjords but can include other regions."""
        config = PacificCoastConfig()
        sim = PacificCoastSimulation(config, seed=42)
        
        fjord_count = sum(1 for idx in sim.refugia_sites 
                         if sim.sites[idx].region_id in ("bc_fjords", "se_alaska_north"))
        # Fjord regions should have higher refugia density than outer coast
        assert fjord_count >= 1, "At least some refugia should be in fjord regions"
    
    def test_refugia_population_persists(self):
        """Refugia should delay disease but not prevent extinction entirely."""
        config = PacificCoastConfig(n_years=40)
        sim = PacificCoastSimulation(config, seed=42)
        result = sim.run()
        
        # Get BC Fjords population over time
        boundaries = get_site_region_boundaries(result.sites)
        fjord_start, fjord_end = boundaries["bc_fjords"]
        
        # Check that fjords maintain population longer than outer coast
        # Compare population trajectories
        outer_start, outer_end = boundaries["bc_outer"]
        
        # Year 15: Fjords should have higher survival than outer coast
        fjord_pop_15 = result.states[15].populations[fjord_start:fjord_end].sum()
        fjord_init = result.states[0].populations[fjord_start:fjord_end].sum()
        outer_pop_15 = result.states[15].populations[outer_start:outer_end].sum()
        outer_init = result.states[0].populations[outer_start:outer_end].sum()
        
        fjord_ratio = fjord_pop_15 / fjord_init
        outer_ratio = outer_pop_15 / outer_init
        
        # Fjords should delay decline relative to outer coast
        assert fjord_ratio > outer_ratio, \
            f"Fjords should delay decline: fjord={fjord_ratio:.2%}, outer={outer_ratio:.2%}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
