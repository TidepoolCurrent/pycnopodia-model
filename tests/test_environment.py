"""Tests for environmental stochasticity module."""
import pytest
import numpy as np


class TestEnvironmentState:
    """Test environmental time series generation."""
    
    def test_recruitment_multiplier_mean(self):
        """Mean recruitment multiplier should be ~1."""
        from pycnopodia.environment import EnvironmentConfig, EnvironmentState
        
        config = EnvironmentConfig(recruitment_cv=0.5)
        env = EnvironmentState(config, n_years=1000, seed=42)
        
        mean_mult = np.mean(env.recruitment_multipliers)
        assert 0.9 < mean_mult < 1.1
    
    def test_recruitment_variability(self):
        """Recruitment should have specified CV."""
        from pycnopodia.environment import EnvironmentConfig, EnvironmentState
        
        target_cv = 0.8
        config = EnvironmentConfig(recruitment_cv=target_cv)
        env = EnvironmentState(config, n_years=1000, seed=42)
        
        actual_cv = np.std(env.recruitment_multipliers) / np.mean(env.recruitment_multipliers)
        assert 0.5 * target_cv < actual_cv < 1.5 * target_cv
    
    def test_temperature_trend(self):
        """Temperature should show warming trend."""
        from pycnopodia.environment import EnvironmentConfig, EnvironmentState
        
        config = EnvironmentConfig(warming_rate=0.03)
        env = EnvironmentState(config, n_years=100, seed=42)
        
        # End should be warmer than start
        early_mean = np.mean(env.temperature_anomaly[:20])
        late_mean = np.mean(env.temperature_anomaly[-20:])
        assert late_mean > early_mean + 1.0  # At least 1C warming
    
    def test_heatwaves_occur(self):
        """Some years should be heatwaves."""
        from pycnopodia.environment import EnvironmentConfig, EnvironmentState
        
        config = EnvironmentConfig(heatwave_probability=0.1)
        env = EnvironmentState(config, n_years=100, seed=42)
        
        n_heatwave_years = len(env.heatwave_years)
        assert 3 < n_heatwave_years < 25  # Roughly 10%
    
    def test_sswd_multiplier_during_heatwave(self):
        """SSWD multiplier should be elevated during heatwaves."""
        from pycnopodia.environment import EnvironmentConfig, EnvironmentState
        
        config = EnvironmentConfig(
            heatwave_probability=1.0,  # Force heatwave
            heatwave_sswd_multiplier=2.0
        )
        env = EnvironmentState(config, n_years=10, seed=42)
        
        mult = env.get_sswd_multiplier(5)
        assert mult >= 2.0


class TestConnectivity:
    """Test connectivity matrix generation."""
    
    def test_stepping_stone_rows_sum_to_one(self):
        """Each row should sum to 1."""
        from pycnopodia.environment import generate_connectivity_matrix
        
        C = generate_connectivity_matrix(5, "stepping_stone", 0.1)
        row_sums = C.sum(axis=1)
        
        np.testing.assert_allclose(row_sums, 1.0, rtol=1e-10)
    
    def test_stepping_stone_adjacency(self):
        """Non-adjacent should have zero connectivity."""
        from pycnopodia.environment import generate_connectivity_matrix
        
        C = generate_connectivity_matrix(5, "stepping_stone", 0.1)
        
        # Position 0 shouldn't connect to position 4
        assert C[0, 4] == 0
        assert C[4, 0] == 0
    
    def test_pnw_connectivity(self):
        """PNW connectivity matrix should be valid."""
        from pycnopodia.environment import PNW_CONNECTIVITY
        
        assert PNW_CONNECTIVITY.shape == (5, 5)
        row_sums = PNW_CONNECTIVITY.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, rtol=1e-10)
