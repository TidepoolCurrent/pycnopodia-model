"""Tests for size-structured population model."""
import pytest
import numpy as np


class TestSizeStructuredPopulation:
    """Test size-structured population."""
    
    def test_initialization(self):
        """Should initialize with size distribution."""
        from pycnopodia import Config
        from pycnopodia.size_structured import SizeStructuredConfig, SizeStructuredPopulation
        
        size_cfg = SizeStructuredConfig()
        base_cfg = Config(initial_population=100)
        pop = SizeStructuredPopulation(size_cfg, base_cfg)
        
        assert pop.n == 100
        assert pop.sizes.min() >= size_cfg.min_size
        assert pop.sizes.max() <= size_cfg.max_size
    
    def test_size_fecundity_relationship(self):
        """Larger individuals should have higher fecundity."""
        from pycnopodia import Config
        from pycnopodia.size_structured import SizeStructuredConfig, SizeStructuredPopulation
        
        size_cfg = SizeStructuredConfig()
        base_cfg = Config(initial_population=50)
        pop = SizeStructuredPopulation(size_cfg, base_cfg)
        
        # Set all to female and adult
        pop.sexes[:] = True
        pop.sizes[:25] = 25  # Small adults
        pop.sizes[25:] = 50  # Large adults
        
        fec = pop.fecundity()
        
        # Large should have higher fecundity
        assert fec[25:].mean() > fec[:25].mean()
    
    def test_sswd_size_susceptibility(self):
        """Larger individuals should be more susceptible to SSWD."""
        from pycnopodia import Config
        from pycnopodia.size_structured import SizeStructuredConfig, SizeStructuredPopulation
        
        size_cfg = SizeStructuredConfig(sswd_size_effect=0.3)
        base_cfg = Config(initial_population=20)
        pop = SizeStructuredPopulation(size_cfg, base_cfg)
        
        # Set sizes
        pop.sizes[:10] = 20  # Small
        pop.sizes[10:] = 50  # Large
        
        suscept = pop.sswd_susceptibility()
        
        # Large should be more susceptible
        assert suscept[10:].mean() > suscept[:10].mean()
    
    def test_growth(self):
        """Individuals should grow over time."""
        from pycnopodia import Config
        from pycnopodia.size_structured import SizeStructuredConfig, SizeStructuredPopulation
        
        size_cfg = SizeStructuredConfig()
        base_cfg = Config(initial_population=50)
        pop = SizeStructuredPopulation(size_cfg, base_cfg)
        
        # Set all to small size
        pop.sizes[:] = 10
        initial_mean = pop.mean_size()
        
        rng = np.random.default_rng(42)
        pop.grow(rng)
        
        # Should have grown
        assert pop.mean_size() > initial_mean
    
    def test_survival_size_dependent(self):
        """Small juveniles should have lower survival than adults."""
        from pycnopodia import Config
        from pycnopodia.size_structured import SizeStructuredConfig, SizeStructuredPopulation
        
        size_cfg = SizeStructuredConfig(
            maturation_size=20,
            juvenile_survival_intercept=0.3,
            adult_survival_max=0.92
        )
        base_cfg = Config(initial_population=20)
        pop = SizeStructuredPopulation(size_cfg, base_cfg)
        
        pop.sizes[:10] = 5   # Juveniles
        pop.sizes[10:] = 30  # Adults
        
        survival = pop.survival_probability()
        
        # Adults should have higher survival
        assert survival[10:].mean() > survival[:10].mean()


class TestRealData:
    """Test real environmental data integration."""
    
    def test_historical_sst(self):
        """Should generate SST time series."""
        from pycnopodia.real_data import get_historical_sst
        
        sst = get_historical_sst("washington", 2000, 2020)
        
        assert len(sst) == 21
        assert 8 < sst.mean() < 14  # Reasonable range for WA
    
    def test_sswd_parameters(self):
        """Should return valid SSWD parameters."""
        from pycnopodia.real_data import get_sswd_parameters
        
        params = get_sswd_parameters("california")
        
        assert "peak_prevalence" in params
        assert 0.8 < params["peak_prevalence"] <= 1.0
        assert 0 < params["endemic_prevalence"] < 0.5
    
    def test_temperature_scenarios(self):
        """Different scenarios should have different trends."""
        from pycnopodia.real_data import generate_temperature_scenario
        
        # Use longer time series for more robust trend detection
        stable, _ = generate_temperature_scenario(100, "stable", region="washington")
        warming, _ = generate_temperature_scenario(100, "ssp585", region="washington")
        
        # Calculate trend (end - start)
        stable_trend = stable[-20:].mean() - stable[:20].mean()
        warming_trend = warming[-20:].mean() - warming[:20].mean()
        
        # SSP5-8.5 should have stronger warming trend
        assert warming_trend > stable_trend + 1
    
    def test_larval_connectivity(self):
        """Connectivity matrix should be valid."""
        from pycnopodia.real_data import get_larval_connectivity
        
        C = get_larval_connectivity()
        
        assert C.shape == (5, 5)
        # Rows should sum to 1
        np.testing.assert_allclose(C.sum(axis=1), 1.0, rtol=1e-10)
