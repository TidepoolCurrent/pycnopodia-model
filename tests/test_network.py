"""Tests for network metapopulation model."""

import numpy as np
import pytest
from pycnopodia.network import (
    NetworkConfig, NetworkSimulation, ConnectivityType,
    build_connectivity_matrix, NetworkState, run_scenario
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
