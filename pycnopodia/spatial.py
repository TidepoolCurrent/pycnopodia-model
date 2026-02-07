"""
Spatial structure module.

Implements grid-based individual tracking for:
- Local density effects on fertilization
- Spatially-explicit disease transmission
- Clustering of resistant genotypes
"""

import numpy as np
from typing import Tuple, Optional, List
from dataclasses import dataclass


@dataclass
class SpatialConfig:
    """Spatial model parameters."""
    
    grid_size: int = 50  # Grid cells per side
    cell_capacity: int = 10  # Max individuals per cell
    
    # Movement
    adult_movement_cells: float = 1.0  # Cells moved per year (mean)
    juvenile_dispersal_cells: float = 5.0  # Settlement dispersal from parents
    
    # Local interactions
    fertilization_radius: int = 2  # Cells for local density calculation
    transmission_radius: int = 1  # Cells for disease transmission
    
    # Density-dependent transmission
    local_transmission_weight: float = 0.7  # Weight of local vs global transmission


class SpatialPopulation:
    """
    Spatially-explicit population on a 2D grid.
    
    Extends base Population with spatial coordinates.
    """
    
    def __init__(
        self,
        config: SpatialConfig,
        n_individuals: int,
        rng: np.random.Generator
    ):
        self.config = config
        self.rng = rng
        
        # Grid dimensions
        self.grid_size = config.grid_size
        
        # Individual positions (x, y coordinates)
        self.x = rng.integers(0, self.grid_size, n_individuals)
        self.y = rng.integers(0, self.grid_size, n_individuals)
        
        # Density grid (updated periodically)
        self._density_grid = None
        self._density_stale = True
    
    @property
    def n(self) -> int:
        return len(self.x)
    
    def update_density_grid(self):
        """Recalculate density per cell."""
        self._density_grid = np.zeros((self.grid_size, self.grid_size), dtype=int)
        for i in range(self.n):
            self._density_grid[self.x[i], self.y[i]] += 1
        self._density_stale = False
    
    def get_density_grid(self) -> np.ndarray:
        """Get current density grid."""
        if self._density_stale:
            self.update_density_grid()
        return self._density_grid
    
    def local_density(self, x: int, y: int, radius: int = None) -> int:
        """
        Count individuals within radius of (x, y).
        """
        if radius is None:
            radius = self.config.fertilization_radius
        
        # Count individuals in neighborhood
        count = 0
        for dx in range(-radius, radius + 1):
            for dy in range(-radius, radius + 1):
                nx = (x + dx) % self.grid_size
                ny = (y + dy) % self.grid_size
                # Count individuals at this cell
                count += np.sum((self.x == nx) & (self.y == ny))
        return count
    
    def local_density_vectorized(self, radius: int = None) -> np.ndarray:
        """
        Get local density for each individual (vectorized).
        """
        if radius is None:
            radius = self.config.fertilization_radius
        
        density = self.get_density_grid()
        
        # Convolution with neighborhood kernel
        from scipy.ndimage import uniform_filter
        kernel_size = 2 * radius + 1
        smoothed = uniform_filter(density.astype(float), size=kernel_size, mode='wrap')
        
        # Look up density at each individual's position
        local_d = np.array([smoothed[self.x[i], self.y[i]] for i in range(self.n)])
        return local_d * kernel_size**2  # Convert from mean to sum
    
    def move_adults(self, alive_mask: np.ndarray):
        """
        Random adult movement.
        """
        n_alive = alive_mask.sum()
        if n_alive == 0:
            return
        
        # Movement distance (Poisson)
        movement = self.config.adult_movement_cells
        dx = self.rng.poisson(movement, n_alive) * self.rng.choice([-1, 1], n_alive)
        dy = self.rng.poisson(movement, n_alive) * self.rng.choice([-1, 1], n_alive)
        
        # Update positions (with wrapping)
        self.x[alive_mask] = (self.x[alive_mask] + dx) % self.grid_size
        self.y[alive_mask] = (self.y[alive_mask] + dy) % self.grid_size
        
        self._density_stale = True
    
    def disperse_settlers(
        self,
        parent_x: np.ndarray,
        parent_y: np.ndarray,
        n_settlers: int
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Disperse settlers around parent locations.
        
        Returns new (x, y) positions for settlers.
        """
        if n_settlers == 0:
            return np.array([]), np.array([])
        
        # Sample parent locations (with replacement)
        parent_idx = self.rng.choice(len(parent_x), n_settlers)
        
        # Dispersal distance (exponential)
        dispersal = self.config.juvenile_dispersal_cells
        dist = self.rng.exponential(dispersal, n_settlers)
        angle = self.rng.uniform(0, 2 * np.pi, n_settlers)
        
        dx = (dist * np.cos(angle)).astype(int)
        dy = (dist * np.sin(angle)).astype(int)
        
        new_x = (parent_x[parent_idx] + dx) % self.grid_size
        new_y = (parent_y[parent_idx] + dy) % self.grid_size
        
        return new_x, new_y
    
    def add_individuals(self, x: np.ndarray, y: np.ndarray):
        """Add new individuals at given positions."""
        self.x = np.concatenate([self.x, x])
        self.y = np.concatenate([self.y, y])
        self._density_stale = True
    
    def remove_individuals(self, mask: np.ndarray):
        """Remove individuals where mask is True."""
        keep = ~mask
        self.x = self.x[keep]
        self.y = self.y[keep]
        self._density_stale = True
    
    def local_transmission_probability(
        self,
        global_prevalence: float,
        infected_mask: np.ndarray
    ) -> np.ndarray:
        """
        Calculate spatially-explicit transmission probability.
        
        Combines global prevalence with local infected density.
        """
        config = self.config
        
        # Global component
        global_weight = 1 - config.local_transmission_weight
        p_global = global_prevalence * global_weight
        
        # Local component - count infected neighbors
        local_weight = config.local_transmission_weight
        radius = config.transmission_radius
        
        # Build infected density grid
        infected_grid = np.zeros((self.grid_size, self.grid_size))
        infected_x = self.x[infected_mask]
        infected_y = self.y[infected_mask]
        for i in range(len(infected_x)):
            infected_grid[infected_x[i], infected_y[i]] += 1
        
        # Smooth with transmission kernel
        from scipy.ndimage import uniform_filter
        kernel_size = 2 * radius + 1
        infected_smooth = uniform_filter(infected_grid, size=kernel_size, mode='wrap')
        
        # Normalize by local density
        density = self.get_density_grid().astype(float)
        density_smooth = uniform_filter(density, size=kernel_size, mode='wrap')
        
        # Local prevalence at each cell
        local_prev = np.zeros_like(infected_smooth)
        nonzero = density_smooth > 0
        local_prev[nonzero] = infected_smooth[nonzero] / density_smooth[nonzero]
        
        # Look up for each individual
        p_local = np.array([local_prev[self.x[i], self.y[i]] for i in range(self.n)])
        p_local *= local_weight
        
        # Combine
        return np.clip(p_global + p_local, 0, 1)


def calculate_spatial_fertilization(
    spatial_pop: SpatialPopulation,
    female_idx: np.ndarray,
    male_idx: np.ndarray,
    base_half_saturation: float
) -> np.ndarray:
    """
    Calculate per-female fertilization success based on local male density.
    
    Returns array of fertilization probabilities for each female.
    """
    if len(female_idx) == 0 or len(male_idx) == 0:
        return np.array([])
    
    female_x = spatial_pop.x[female_idx]
    female_y = spatial_pop.y[female_idx]
    male_x = spatial_pop.x[male_idx]
    male_y = spatial_pop.y[male_idx]
    
    radius = spatial_pop.config.fertilization_radius
    grid_size = spatial_pop.grid_size
    
    # For each female, count local males
    fert_probs = np.zeros(len(female_idx))
    
    for i, (fx, fy) in enumerate(zip(female_x, female_y)):
        # Count males within radius
        local_males = 0
        for mx, my in zip(male_x, male_y):
            # Wrapped distance
            dx = min(abs(fx - mx), grid_size - abs(fx - mx))
            dy = min(abs(fy - my), grid_size - abs(fy - my))
            if dx <= radius and dy <= radius:
                local_males += 1
        
        # Saturating fertilization
        fert_probs[i] = local_males**2 / (local_males**2 + base_half_saturation**2)
    
    return fert_probs
