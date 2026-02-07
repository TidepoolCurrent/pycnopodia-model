"""
Hierarchical Sweepstakes Reproductive Success model.

More realistic SRS that accounts for:
- Spatial patchiness (some sites succeed, others fail)
- Temporal matching (spawn timing × conditions × settlement)
- Variance at multiple levels
"""

import numpy as np
from typing import Tuple
from dataclasses import dataclass


@dataclass
class HierarchicalSRSConfig:
    """Parameters for hierarchical SRS model."""
    
    # Patch-level success
    n_patches: int = 10  # Number of spawning patches
    patch_success_prob: float = 0.3  # Probability a patch has successful recruitment
    
    # Within-patch breeding
    within_patch_breeding_mean: float = 0.15  # Mean breeding fraction in successful patches
    within_patch_breeding_cv: float = 0.5  # CV of breeding fraction
    
    # Temporal matching
    temporal_match_prob: float = 0.4  # Probability of favorable conditions
    
    # Settlement success
    settlement_cv: float = 1.0  # Variation in settlement success


class HierarchicalSRS:
    """
    Hierarchical model of sweepstakes reproductive success.
    
    Three levels of stochasticity:
    1. Patch-level: Which spawning aggregations succeed?
    2. Individual-level: Which adults in successful patches breed?
    3. Settlement-level: How many larvae from each pairing settle?
    
    This creates variance in reproductive success that matches
    empirical observations better than a single Beta draw.
    """
    
    def __init__(self, config: HierarchicalSRSConfig, seed: int = None):
        self.config = config
        self.rng = np.random.default_rng(seed)
    
    def sample_breeding_success(self, n_adults: int) -> Tuple[np.ndarray, dict]:
        """
        Determine which adults breed and their relative reproductive success.
        
        Parameters
        ----------
        n_adults : int
            Number of reproductive adults
            
        Returns
        -------
        breeding_weights : ndarray, shape (n_adults,)
            Relative reproductive contribution of each adult (0 = doesn't breed)
        diagnostics : dict
            Information about the SRS draw
        """
        config = self.config
        
        if n_adults == 0:
            return np.array([]), {"n_breeders": 0}
        
        # Assign adults to patches
        patch_assignments = self.rng.integers(0, config.n_patches, n_adults)
        
        # Determine which patches succeed (temporal matching)
        temporal_success = self.rng.random() < config.temporal_match_prob
        if not temporal_success:
            # Complete reproductive failure this year
            return np.zeros(n_adults), {
                "temporal_match": False,
                "n_breeders": 0,
                "n_successful_patches": 0
            }
        
        # Which patches have successful recruitment?
        patch_success = self.rng.random(config.n_patches) < config.patch_success_prob
        n_successful = patch_success.sum()
        
        if n_successful == 0:
            return np.zeros(n_adults), {
                "temporal_match": True,
                "n_breeders": 0,
                "n_successful_patches": 0
            }
        
        # Within successful patches, determine who breeds
        breeding_weights = np.zeros(n_adults)
        
        for patch_id in range(config.n_patches):
            if not patch_success[patch_id]:
                continue
            
            # Adults in this patch
            in_patch = patch_assignments == patch_id
            n_in_patch = in_patch.sum()
            
            if n_in_patch == 0:
                continue
            
            # Sample breeding fraction for this patch
            mean = config.within_patch_breeding_mean
            cv = config.within_patch_breeding_cv
            
            # Lognormal for breeding fraction
            sigma = np.sqrt(np.log(1 + cv**2))
            mu = np.log(mean) - sigma**2/2
            breeding_frac = min(1.0, np.exp(self.rng.normal(mu, sigma)))
            
            # Select breeders
            n_breed = max(1, int(n_in_patch * breeding_frac))
            patch_indices = np.where(in_patch)[0]
            breeders = self.rng.choice(patch_indices, size=min(n_breed, n_in_patch), replace=False)
            
            # Assign weights (can vary - some breeders more successful)
            if len(breeders) > 0:
                # Exponential variation in success among breeders
                weights = self.rng.exponential(1.0, len(breeders))
                weights /= weights.sum()
                breeding_weights[breeders] = weights
        
        # Normalize total
        if breeding_weights.sum() > 0:
            breeding_weights /= breeding_weights.sum()
        
        n_breeders = (breeding_weights > 0).sum()
        
        return breeding_weights, {
            "temporal_match": True,
            "n_breeders": n_breeders,
            "n_successful_patches": n_successful,
            "breeding_fraction": n_breeders / n_adults if n_adults > 0 else 0
        }
    
    def calculate_effective_breeders(self, weights: np.ndarray) -> float:
        """
        Calculate effective number of breeders from weight distribution.
        
        Uses the formula: Ne = 1 / sum(p_i^2)
        where p_i is the relative contribution of individual i.
        """
        weights = weights[weights > 0]
        if len(weights) == 0:
            return 0.0
        
        # Normalize
        p = weights / weights.sum()
        
        # Effective number
        ne = 1.0 / np.sum(p**2)
        return ne


def sample_parent_pairs(
    female_idx: np.ndarray,
    male_idx: np.ndarray,
    female_weights: np.ndarray,
    male_weights: np.ndarray,
    n_pairs: int,
    rng: np.random.Generator
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Sample parent pairs weighted by breeding success.
    
    Parameters
    ----------
    female_idx, male_idx : ndarray
        Indices of breeding females/males in population
    female_weights, male_weights : ndarray
        Relative breeding success weights
    n_pairs : int
        Number of matings to generate
    rng : Generator
        Random number generator
        
    Returns
    -------
    mother_idx, father_idx : ndarray
        Indices of parents for each pairing
    """
    if len(female_idx) == 0 or len(male_idx) == 0 or n_pairs == 0:
        return np.array([]), np.array([])
    
    # Normalize weights
    f_weights = female_weights[female_idx]
    m_weights = male_weights[male_idx]
    
    if f_weights.sum() == 0 or m_weights.sum() == 0:
        return np.array([]), np.array([])
    
    f_probs = f_weights / f_weights.sum()
    m_probs = m_weights / m_weights.sum()
    
    # Sample parents
    mothers = rng.choice(female_idx, size=n_pairs, p=f_probs)
    fathers = rng.choice(male_idx, size=n_pairs, p=m_probs)
    
    return mothers, fathers
