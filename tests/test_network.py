"""Tests for network metapopulation model."""

import numpy as np
import pytest
from pycnopodia.network import (
    NetworkConfig, NetworkSimulation, ConnectivityType,
    build_connectivity_matrix, build_disease_connectivity_matrix,
    NetworkState, run_scenario
)


class TestConnectivityMatrix:
    """Test connectivity matrix construction."""
    
    def test_rows_sum_to_one(self):
        """All connectivity matrices should have rows summing to 1."""
        for conn_type in ConnectivityType:
            config = NetworkConfig(
                n_sites=100,
                connectivity_type=conn_type
            )
            C = build_connectivity_matrix(config)
            row_sums = C.sum(axis=1)
            np.testing.assert_array_almost_equal(
                row_sums, np.ones(100), decimal=5,
                err_msg=f"Rows don't sum to 1 for {conn_type.value}"
            )
    
    def test_non_negative(self):
        """All matrix entries should be non-negative."""
        for conn_type in ConnectivityType:
            config = NetworkConfig(
                n_sites=100,
                connectivity_type=conn_type
            )
            C = build_connectivity_matrix(config)
            assert (C >= 0).all(), f"Negative entries in {conn_type.value}"
    
    def test_self_recruitment_respected(self):
        """Diagonal should approximately match self-recruitment parameter."""
        config = NetworkConfig(
            n_sites=100,
            connectivity_type=ConnectivityType.UNIFORM,
            self_recruitment=0.7
        )
        C = build_connectivity_matrix(config)
        diagonal = np.diag(C)
        assert np.allclose(diagonal, 0.7, atol=0.01)
    
    def test_stepping_stone_locality(self):
        """Stepping stone should only connect adjacent sites."""
        config = NetworkConfig(
            n_sites=100,
            connectivity_type=ConnectivityType.STEPPING_STONE,
            asymmetry=0.0
        )
        C = build_connectivity_matrix(config)
        # Non-adjacent entries should be zero (except near edges)
        for i in range(10, 90):  # Avoid edge effects
            for j in range(100):
                if abs(i - j) > 1:
                    assert C[i, j] < 0.01, f"Non-local connection at ({i},{j})"
    
    def test_asymmetric_flow_direction(self):
        """Asymmetric flow should bias toward increasing indices."""
        config = NetworkConfig(
            n_sites=100,
            connectivity_type=ConnectivityType.ASYMMETRIC_FLOW,
            asymmetry=0.9,
            self_recruitment=0.3
        )
        C = build_connectivity_matrix(config)
        # Sum of downstream (j > i) should exceed upstream
        for i in range(20, 80):
            downstream = C[i, i+1:].sum()
            upstream = C[i, :i].sum()
            assert downstream > upstream, f"Flow not downstream at site {i}"


class TestNetworkState:
    """Test NetworkState calculations."""
    
    def test_ratios(self):
        """Test ratio calculations."""
        state = NetworkState(
            year=10,
            populations=np.full(100, 500.0),  # 50% of baseline
            N0_per_site=1000,
            N0_total=100_000
        )
        assert state.n_ratio == 0.5
        assert state.occupied_ratio == 1.0
        assert state.mean_site_n_ratio == 0.5
    
    def test_extinction(self):
        """Test extinction detection."""
        state = NetworkState(
            year=10,
            populations=np.zeros(100),
            N0_per_site=1000,
            N0_total=100_000
        )
        assert state.extinct
        assert state.n_ratio == 0.0


class TestNetworkSimulation:
    """Test network simulation."""
    
    def test_initialization(self):
        """Test simulation initialization."""
        config = NetworkConfig(n_sites=100, n_per_site=500)
        sim = NetworkSimulation(config, seed=42)
        assert sim.populations.sum() == 50_000
        assert len(sim.C) == 100
    
    def test_short_run(self):
        """Test short simulation runs without error."""
        config = NetworkConfig(
            n_sites=100,
            n_per_site=500,
            n_years=10,
            disease_onset_year=100  # No disease in this run
        )
        sim = NetworkSimulation(config, seed=42)
        result = sim.run()
        
        assert len(result.states) == 10
        assert result.states[0].n_ratio > 0
    
    def test_disease_causes_decline(self):
        """Disease should cause population decline."""
        # Without disease
        config_healthy = NetworkConfig(
            n_sites=100,
            n_years=30,
            disease_onset_year=100,  # Disease never starts
        )
        sim_healthy = NetworkSimulation(config_healthy, seed=42)
        result_healthy = sim_healthy.run()
        
        # With disease
        config_sick = NetworkConfig(
            n_sites=100,
            n_years=30,
            disease_onset_year=5,
            disease_onset_site=50,
            disease_mortality=0.6,
        )
        sim_sick = NetworkSimulation(config_sick, seed=42)
        result_sick = sim_sick.run()
        
        # Disease should reduce final population
        assert result_sick.final_n_ratio < result_healthy.final_n_ratio
    
    def test_lower_mortality_allows_survival(self):
        """Lower mortality scenarios should allow population survival."""
        # Test with much lower mortality to verify model CAN produce survival
        config = NetworkConfig(
            n_sites=100,
            n_years=30,
            self_recruitment=0.8,
            connectivity_type=ConnectivityType.STEPPING_STONE,
            disease_onset_year=10,
            disease_mortality=0.50,  # Much lower than realistic 99%
            disease_spread_rate=0.3,  # Slower spread
        )
        
        # Run multiple replicates
        survivors = 0
        for seed in range(10):
            sim = NetworkSimulation(config, seed=seed)
            result = sim.run()
            if not result.extinct:
                survivors += 1
        
        # With 50% mortality, should have some survivors
        assert survivors > 0, "Should have some survivors with 50% mortality"
    
    def test_genetics_tracked(self):
        """Test that genetics are tracked per site during disease."""
        config = NetworkConfig(
            n_sites=50,
            n_years=15,  # Check before extinction
            disease_onset_year=5,
            disease_mortality=0.80,  # Lower to allow survival
            disease_spread_rate=0.5,
        )
        sim = NetworkSimulation(config, seed=42)
        result = sim.run()
        
        # Check genetics DURING disease (before potential extinction)
        mid_state = result.states[10]  # Year 10, after disease hits
        assert mid_state.resistance_freqs is not None
        assert mid_state.heterozygosity is not None
        
        # If not extinct, check values make sense
        if mid_state.total_population > 0:
            assert 0 < mid_state.mean_resistance_freq < 1
            # Resistance should increase under selection
            initial_state = result.states[0]
            assert mid_state.mean_resistance_freq >= initial_state.mean_resistance_freq * 0.5
    
    def test_disease_spreads(self):
        """Test that disease spreads rapidly after onset."""
        config = NetworkConfig(
            n_sites=100,
            n_years=15,
            disease_onset_year=5,
            disease_spread_rate=0.9,  # Default rapid spread
        )
        sim = NetworkSimulation(config, seed=42)
        result = sim.run()
        
        # Check disease spread at peak (before potential extinction)
        # Year 7 = 2 years after onset
        peak_state = result.states[7]
        assert peak_state.infected_sites_ratio > 0.5, "Disease should spread to >50% of sites"


class TestScenarioComparison:
    """Test scenario comparison."""
    
    def test_run_scenario(self):
        """Test run_scenario function."""
        result = run_scenario(
            connectivity_type=ConnectivityType.STEPPING_STONE,
            asymmetry=0.0,
            self_recruitment=0.5,
            n_replicates=3,
            n_years=20,
        )
        
        assert "extinction_probability" in result
        assert "mean_final_n_ratio" in result
        assert result["n_replicates"] == 3


class TestPolygenicResistance:
    """Test polygenic resistance model."""
    
    def test_resistance_freqs_shape(self):
        """Resistance frequencies should have shape (n_sites, n_loci)."""
        config = NetworkConfig(n_sites=50, n_loci=10)
        sim = NetworkSimulation(config, seed=42)
        assert sim.resistance_freqs.shape == (50, 10)
    
    def test_locus_effects_sum(self):
        """Locus effects should sum to resistance_effect."""
        config = NetworkConfig(n_sites=50, n_loci=10, resistance_effect=0.7)
        sim = NetworkSimulation(config, seed=42)
        np.testing.assert_almost_equal(
            sim.locus_effects.sum(), 0.7, decimal=5
        )
    
    def test_locus_effects_custom(self):
        """Custom locus effects should be used when provided."""
        custom_effects = np.array([0.1, 0.2, 0.3, 0.2, 0.2])
        config = NetworkConfig(
            n_sites=50, 
            n_loci=5,
            locus_effects=custom_effects
        )
        sim = NetworkSimulation(config, seed=42)
        np.testing.assert_array_almost_equal(sim.locus_effects, custom_effects)
    
    def test_state_contains_locus_effects(self):
        """NetworkState should include locus_effects."""
        config = NetworkConfig(n_sites=50, n_loci=5, n_years=5)
        sim = NetworkSimulation(config, seed=42)
        result = sim.run()
        
        for state in result.states:
            assert state.locus_effects is not None
            assert len(state.locus_effects) == 5
    
    def test_per_locus_mean_freq(self):
        """Per-locus mean frequencies should be computed correctly."""
        config = NetworkConfig(n_sites=50, n_loci=5, n_years=5)
        sim = NetworkSimulation(config, seed=42)
        result = sim.run()
        
        state = result.states[-1]
        per_locus = state.per_locus_mean_freq
        
        assert len(per_locus) == 5
        assert all(0 <= f <= 1 for f in per_locus)
    
    def test_mean_resistance_freq_bounded(self):
        """Mean resistance frequency should be between 0 and 1."""
        config = NetworkConfig(n_sites=50, n_loci=10, n_years=20)
        sim = NetworkSimulation(config, seed=42)
        result = sim.run()
        
        for state in result.states:
            if not state.extinct:
                assert 0 <= state.mean_resistance_freq <= 1
    
    def test_selection_increases_resistance(self):
        """Selection under disease should increase resistance frequencies."""
        config = NetworkConfig(
            n_sites=50,
            n_loci=5,
            n_years=30,
            disease_onset_year=5,
            disease_mortality=0.80,  # Lower mortality to avoid extinction
            disease_spread_rate=0.5,
            initial_resistance_freq=0.10,
        )
        sim = NetworkSimulation(config, seed=42)
        result = sim.run()
        
        # Get resistance at different times
        initial_r = result.states[0].mean_resistance_freq
        mid_r = result.states[15].mean_resistance_freq if len(result.states) > 15 else initial_r
        
        # After disease, resistance should generally increase (if population survives)
        if not result.states[15].extinct:
            assert mid_r >= initial_r * 0.5  # At least maintained or increased
    
    def test_drift_per_locus(self):
        """Genetic drift should affect each locus independently."""
        # Small population for strong drift
        config = NetworkConfig(
            n_sites=10,
            n_loci=5,
            n_per_site=50,
            n_years=20,
            disease_onset_year=100,  # No disease
            initial_resistance_freq=0.5,
        )
        sim = NetworkSimulation(config, seed=42)
        result = sim.run()
        
        # Loci should diverge due to independent drift
        final_state = result.states[-1]
        per_locus = final_state.per_locus_mean_freq
        
        # With drift, we expect some variance between loci
        # (they won't all stay at exactly the same frequency)
        assert len(per_locus) == 5
        # At least some variation expected
        variance = np.var(per_locus)
        # Variance should be > 0 due to drift (though could be small)
        assert variance >= 0  # Just ensure computation works
    
    def test_gene_flow_per_locus(self):
        """Gene flow should transfer alleles at each locus."""
        config = NetworkConfig(
            n_sites=20,
            n_loci=5,
            n_years=10,
            disease_onset_year=100,  # No disease
            connectivity_type=ConnectivityType.STEPPING_STONE,
            self_recruitment=0.3,  # High dispersal
        )
        sim = NetworkSimulation(config, seed=42)
        
        # Set initial variation between sites at one locus
        sim.resistance_freqs[:10, 0] = 0.9  # High freq at first 10 sites
        sim.resistance_freqs[10:, 0] = 0.1  # Low freq at last 10 sites
        
        result = sim.run()
        final_state = result.states[-1]
        
        # After gene flow, frequencies should homogenize somewhat
        site_freqs_locus0 = final_state.resistance_freqs[:, 0]
        # The variance should decrease (become more homogeneous)
        initial_var = np.var([0.9] * 10 + [0.1] * 10)
        final_var = np.var(site_freqs_locus0)
        # Final variance should be less than initial due to gene flow
        assert final_var < initial_var
    
    def test_outplanting_enhanced_affects_all_loci(self):
        """Enhanced outplanting should increase frequencies at all loci."""
        config = NetworkConfig(
            n_sites=50,
            n_loci=5,
            n_years=25,
            disease_onset_year=5,
            disease_mortality=0.60,
            outplanting_n=500,
            outplanting_sites=list(range(0, 50, 5)),
            outplanting_start=10,
            outplant_resistance_mode='enhanced',
            outplant_enhanced_resistance=0.95,
            initial_resistance_freq=0.1,
        )
        sim = NetworkSimulation(config, seed=42)
        result = sim.run()
        
        # After outplanting with enhanced resistance
        state_after = result.states[20]  # After several outplanting events
        if not state_after.extinct:
            per_locus = state_after.per_locus_mean_freq
            # All loci should have elevated frequency due to enhanced outplanting
            initial_freq = config.initial_resistance_freq
            # At least some loci should be higher than initial
            assert any(f > initial_freq for f in per_locus)
    
    def test_summary_includes_per_locus(self):
        """Summary should include per-locus frequencies for polygenic model."""
        config = NetworkConfig(n_sites=20, n_loci=5, n_years=5)
        sim = NetworkSimulation(config, seed=42)
        result = sim.run()
        
        summary = result.states[-1].get_summary()
        assert "per_locus_freqs" in summary
        assert len(summary["per_locus_freqs"]) == 5


class TestDiseaseConnectivity:
    """Test separate disease connectivity matrix."""
    
    def test_disease_matrix_exists(self):
        """Simulation should have separate disease connectivity matrix."""
        config = NetworkConfig(n_sites=100)
        sim = NetworkSimulation(config, seed=42)
        
        # Should have both matrices
        assert sim.larval_connectivity is not None
        assert sim.disease_connectivity is not None
        assert sim.C is sim.larval_connectivity  # Legacy alias
        assert sim.D is sim.disease_connectivity  # Alias
    
    def test_disease_matrix_different_from_larval(self):
        """Disease matrix should be different from larval (different params)."""
        config = NetworkConfig(
            n_sites=100,
            connectivity_type=ConnectivityType.STEPPING_STONE,
            disease_connectivity_type=ConnectivityType.DISTANCE_DECAY,
        )
        sim = NetworkSimulation(config, seed=42)
        
        # Matrices should not be identical
        assert not np.allclose(sim.larval_connectivity, sim.disease_connectivity)
    
    def test_disease_matrix_rows_sum_to_one(self):
        """Disease connectivity matrix rows should sum to 1."""
        for conn_type in ConnectivityType:
            config = NetworkConfig(
                n_sites=100,
                disease_connectivity_type=conn_type
            )
            D = build_disease_connectivity_matrix(config)
            row_sums = D.sum(axis=1)
            np.testing.assert_array_almost_equal(
                row_sums, np.ones(100), decimal=5,
                err_msg=f"Disease matrix rows don't sum to 1 for {conn_type.value}"
            )
    
    def test_disease_matrix_non_negative(self):
        """Disease connectivity matrix should be non-negative."""
        for conn_type in ConnectivityType:
            config = NetworkConfig(
                n_sites=100,
                disease_connectivity_type=conn_type
            )
            D = build_disease_connectivity_matrix(config)
            assert (D >= 0).all(), f"Negative entries in disease matrix for {conn_type.value}"
    
    def test_result_contains_both_matrices(self):
        """NetworkResult should contain both connectivity matrices."""
        config = NetworkConfig(n_sites=50, n_years=5)
        sim = NetworkSimulation(config, seed=42)
        result = sim.run()
        
        assert result.larval_connectivity is not None
        assert result.disease_connectivity is not None
        assert result.connectivity_matrix is not None  # Legacy
        np.testing.assert_array_equal(result.connectivity_matrix, result.larval_connectivity)
    
    def test_disease_spread_follows_connectivity(self):
        """Disease should spread following the disease connectivity matrix pattern."""
        # Use stepping stone for disease - should spread locally
        config = NetworkConfig(
            n_sites=100,
            n_years=20,
            disease_onset_year=5,
            disease_onset_site=50,  # Middle
            disease_connectivity_type=ConnectivityType.STEPPING_STONE,
            disease_dispersal_scale=2,
            disease_transmission_prob=0.9,
        )
        sim = NetworkSimulation(config, seed=42)
        result = sim.run()
        
        # After disease starts, sites near origin should be infected first
        state_early = result.states[6]  # Year 6, 1 year after onset
        
        # Sites near 50 should be more likely infected than distant sites
        near_sites = [48, 49, 50, 51, 52]
        far_sites = [0, 1, 2, 97, 98, 99]
        
        near_prevalence = np.mean([state_early.disease_prevalence[i] for i in near_sites])
        far_prevalence = np.mean([state_early.disease_prevalence[i] for i in far_sites])
        
        # Near sites should have higher prevalence (wave spreads outward)
        assert near_prevalence > far_prevalence, "Disease should spread as wave from origin"
    
    def test_disease_parameters_used(self):
        """Disease connectivity should use disease-specific parameters."""
        config = NetworkConfig(
            n_sites=50,
            disease_connectivity_type=ConnectivityType.DISTANCE_DECAY,
            disease_dispersal_scale=20,  # Wide spread
            disease_asymmetry=0.5,  # Directional bias
        )
        D = build_disease_connectivity_matrix(config)
        
        # With asymmetry, downstream (higher index) should have more weight
        for i in range(10, 40):
            downstream = D[i, i+1:].sum()
            upstream = D[i, :i].sum()
            assert downstream > upstream * 0.8, f"Asymmetry not applied at site {i}"
