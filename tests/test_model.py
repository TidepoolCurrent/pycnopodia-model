"""Tests for Pycnopodia population model."""
import pytest
import numpy as np


class TestConfig:
    """Test configuration."""
    
    def test_default_config_valid(self):
        from pycnopodia import Config
        config = Config()
        assert config.validate()
    
    def test_presets_valid(self):
        from pycnopodia.config import PRESETS
        for name, config in PRESETS.items():
            assert config.validate(), f"Preset {name} invalid"


class TestPopulation:
    """Test population class."""
    
    def test_initialization(self):
        from pycnopodia import Population, Config
        config = Config(initial_population=100)
        pop = Population(config)
        
        assert pop.n == 100
        assert pop.genomes.shape == (100, config.n_loci)
    
    def test_resistance_calculation(self):
        from pycnopodia import Population, Config
        config = Config(initial_population=50, n_loci=10, per_allele_effect=0.05)
        pop = Population(config)
        
        scores = pop.resistance_scores()
        assert len(scores) == 50
        assert scores.min() >= 0
        assert scores.max() <= 1
    
    def test_mortality(self):
        from pycnopodia import Population, Config
        config = Config(initial_population=100)
        pop = Population(config)
        
        initial_n = pop.n
        mortality = np.full(100, 0.5)  # 50% mortality
        pop.apply_mortality(mortality)
        
        # Should lose roughly half
        assert 30 < pop.n < 70


class TestReproduction:
    """Test reproduction module."""
    
    def test_breeding_fraction(self):
        from pycnopodia.reproduction import sample_breeding_fraction
        from pycnopodia import Config
        
        config = Config(srs_mean=0.1)
        rng = np.random.default_rng(42)
        
        fractions = [sample_breeding_fraction(config, rng) for _ in range(100)]
        mean_frac = np.mean(fractions)
        
        # Should be close to configured mean
        assert 0.05 < mean_frac < 0.20
    
    def test_allee_effect(self):
        from pycnopodia.reproduction import fertilization_success_saturating
        
        # Low density = low success
        low = fertilization_success_saturating(5, 30)
        high = fertilization_success_saturating(100, 30)
        
        assert low < high
        assert low < 0.5
        assert high > 0.9


class TestDisease:
    """Test disease module."""
    
    def test_prevalence_before_onset(self):
        from pycnopodia.disease import calculate_prevalence
        from pycnopodia import Config
        
        config = Config(sswd_onset_year=10)
        prev = calculate_prevalence(5, config)
        
        assert prev == 0.0
    
    def test_prevalence_at_peak(self):
        from pycnopodia.disease import calculate_prevalence
        from pycnopodia import Config
        
        config = Config(sswd_onset_year=10, sswd_peak_prevalence=0.9)
        prev = calculate_prevalence(10, config)
        
        assert prev == pytest.approx(0.9, rel=0.01)
    
    def test_prevalence_decays(self):
        from pycnopodia.disease import calculate_prevalence
        from pycnopodia import Config
        
        config = Config(sswd_onset_year=10)
        prev_early = calculate_prevalence(11, config)
        prev_late = calculate_prevalence(30, config)
        
        assert prev_early > prev_late


class TestSimulation:
    """Test simulation class."""
    
    def test_short_simulation(self):
        from pycnopodia import Simulation, Config
        
        config = Config(n_years=10, initial_population=200)
        sim = Simulation(config, seed=42)
        result = sim.run()
        
        assert len(result.records) == 10
        assert result.records[0].n_total > 0
    
    def test_no_disease_baseline(self):
        from pycnopodia import Simulation
        from pycnopodia.config import PRESETS
        
        config = PRESETS["no_disease"]
        config.n_years = 20
        
        sim = Simulation(config, seed=42)
        result = sim.run()
        
        # Should maintain population without disease
        assert not result.extinct
        assert result.final_population > 100
    
    def test_extinction_without_intervention(self):
        from pycnopodia import Simulation, Config
        
        # Harsh conditions, no outplanting
        config = Config(
            n_years=50,
            initial_population=500,
            sswd_onset_year=5,
            sswd_peak_prevalence=0.95,
            sswd_adult_mortality=0.80,
            outplanting_n_per_year=0,
        )
        
        # Run multiple replicates - at least some should go extinct
        extinctions = 0
        for seed in range(10):
            sim = Simulation(config, seed=seed)
            result = sim.run()
            if result.extinct:
                extinctions += 1
        
        # Most should go extinct under harsh conditions
        assert extinctions >= 5


class TestIntervention:
    """Test intervention module."""
    
    def test_broodstock_initialization(self):
        from pycnopodia.intervention import Broodstock
        from pycnopodia import Config
        
        config = Config(broodstock_n_parents=20)
        broodstock = Broodstock(config)
        
        assert len(broodstock.genomes) == 20
        assert broodstock.n_females + broodstock.n_males == 20
    
    def test_broodstock_breeding(self):
        from pycnopodia.intervention import Broodstock
        from pycnopodia import Config
        
        config = Config(broodstock_n_parents=20, n_loci=10)
        broodstock = Broodstock(config)
        rng = np.random.default_rng(42)
        
        genomes, ages, sexes = broodstock.breed_offspring(50, rng)
        
        assert len(genomes) == 50
        assert genomes.shape == (50, 10)
        assert all(ages == 0)
