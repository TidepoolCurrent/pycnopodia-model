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
            disease_peak_prevalence=0.0  # No disease
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
            disease_peak_prevalence=0.0,
            disease_endemic_prevalence=0.0,
        )
        sim_healthy = NetworkSimulation(config_healthy, seed=42)
        result_healthy = sim_healthy.run()
        
        # With disease
        config_sick = NetworkConfig(
            n_sites=100,
            n_years=30,
            disease_onset=5,
            disease_peak_prevalence=0.8,
            disease_mortality=0.5,
        )
        sim_sick = NetworkSimulation(config_sick, seed=42)
        result_sick = sim_sick.run()
        
        # Disease should reduce final population
        assert result_sick.final_n_ratio < result_healthy.final_n_ratio
    
    def test_high_self_recruitment_prevents_extinction(self):
        """High self-recruitment should prevent extinction."""
        config = NetworkConfig(
            n_sites=100,
            n_years=50,
            self_recruitment=0.8,
            connectivity_type=ConnectivityType.STEPPING_STONE,
            disease_peak_prevalence=0.7,
            disease_mortality=0.5,
        )
        
        # Run multiple replicates
        extinctions = 0
        for seed in range(10):
            sim = NetworkSimulation(config, seed=seed)
            result = sim.run()
            if result.extinct:
                extinctions += 1
        
        # Should rarely go extinct with high self-recruitment
        assert extinctions < 3, "Too many extinctions with high self-recruitment"


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
