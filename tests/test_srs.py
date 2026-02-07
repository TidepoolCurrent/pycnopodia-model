"""Tests for hierarchical SRS module."""
import pytest
import numpy as np


class TestHierarchicalSRS:
    """Test hierarchical sweepstakes reproductive success."""
    
    def test_zero_adults(self):
        """Should handle zero adults gracefully."""
        from pycnopodia.srs import HierarchicalSRSConfig, HierarchicalSRS
        
        config = HierarchicalSRSConfig()
        srs = HierarchicalSRS(config, seed=42)
        
        weights, diag = srs.sample_breeding_success(0)
        
        assert len(weights) == 0
        assert diag["n_breeders"] == 0
    
    def test_some_adults_breed(self):
        """At least some adults should breed in favorable conditions."""
        from pycnopodia.srs import HierarchicalSRSConfig, HierarchicalSRS
        
        config = HierarchicalSRSConfig(
            temporal_match_prob=1.0,  # Always favorable
            patch_success_prob=0.5
        )
        srs = HierarchicalSRS(config, seed=42)
        
        weights, diag = srs.sample_breeding_success(100)
        
        # Should have some breeders
        assert diag["n_breeders"] > 0
        assert diag["n_breeders"] < 100  # But not all
    
    def test_complete_failure_possible(self):
        """Should sometimes have complete reproductive failure."""
        from pycnopodia.srs import HierarchicalSRSConfig, HierarchicalSRS
        
        config = HierarchicalSRSConfig(
            temporal_match_prob=0.0  # Never favorable
        )
        srs = HierarchicalSRS(config, seed=42)
        
        weights, diag = srs.sample_breeding_success(100)
        
        assert diag["n_breeders"] == 0
        assert not diag["temporal_match"]
    
    def test_weights_sum_to_one(self):
        """Breeding weights should sum to 1 when there are breeders."""
        from pycnopodia.srs import HierarchicalSRSConfig, HierarchicalSRS
        
        config = HierarchicalSRSConfig(
            temporal_match_prob=1.0,
            patch_success_prob=1.0
        )
        srs = HierarchicalSRS(config, seed=42)
        
        weights, diag = srs.sample_breeding_success(50)
        
        if diag["n_breeders"] > 0:
            assert np.isclose(weights.sum(), 1.0)
    
    def test_effective_breeders_less_than_actual(self):
        """Ne should be less than or equal to actual breeders."""
        from pycnopodia.srs import HierarchicalSRSConfig, HierarchicalSRS
        
        config = HierarchicalSRSConfig(
            temporal_match_prob=1.0,
            patch_success_prob=0.8
        )
        srs = HierarchicalSRS(config, seed=42)
        
        weights, diag = srs.sample_breeding_success(100)
        ne = srs.calculate_effective_breeders(weights)
        
        assert ne <= diag["n_breeders"]
    
    def test_variance_across_years(self):
        """Should produce high variance in breeding success across years."""
        from pycnopodia.srs import HierarchicalSRSConfig, HierarchicalSRS
        
        config = HierarchicalSRSConfig()
        srs = HierarchicalSRS(config, seed=42)
        
        breeding_fractions = []
        for _ in range(100):
            weights, diag = srs.sample_breeding_success(200)
            frac = diag.get("breeding_fraction", 0)
            breeding_fractions.append(frac)
        
        # Should have high variance (CV > 0.3)
        mean_frac = np.mean(breeding_fractions)
        if mean_frac > 0:
            cv = np.std(breeding_fractions) / mean_frac
            assert cv > 0.3  # Moderate to high variance


class TestParentPairs:
    """Test parent pair sampling."""
    
    def test_weighted_sampling(self):
        """Higher weights should produce more offspring."""
        from pycnopodia.srs import sample_parent_pairs
        
        rng = np.random.default_rng(42)
        
        female_idx = np.array([0, 1, 2])
        male_idx = np.array([3, 4, 5])
        
        # Unequal weights
        female_weights = np.array([0.8, 0.1, 0.1, 0, 0, 0])
        male_weights = np.array([0, 0, 0, 0.1, 0.1, 0.8])
        
        mothers, fathers = sample_parent_pairs(
            female_idx, male_idx,
            female_weights, male_weights,
            n_pairs=1000,
            rng=rng
        )
        
        # Female 0 should be most common mother
        assert np.sum(mothers == 0) > np.sum(mothers == 1)
        # Male 5 should be most common father
        assert np.sum(fathers == 5) > np.sum(fathers == 3)
