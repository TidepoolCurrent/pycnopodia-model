"""
Generate dense site network along Pacific coast from coastline waypoints.

Approach:
1. Define coastline as waypoints (lat/lon) from SE Alaska to S. California
2. Interpolate sites every ~25 km along the coast
3. Flag known fjord locations from real_sites.py
4. Add fjord arm sites branching off the main coast

This gives ~300-500 sites with real geography.
"""

import math
import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional

# Import the curated fjord sites as anchors
from data.real_sites import ALL_SITES as CURATED_SITES, RealSite, haversine_km


# ============================================================
# Coastline waypoints (main coast, excluding fjord interiors)
# These trace the outer coast from SE Alaska to S. California
# ============================================================
COASTLINE_WAYPOINTS = [
    # SE Alaska North — Inside Passage / outer coast
    (58.90, -136.20, "se_alaska_north"),  # Glacier Bay entrance
    (58.50, -135.80, "se_alaska_north"),  # Icy Strait
    (58.30, -134.90, "se_alaska_north"),  # Near Juneau
    (57.80, -135.50, "se_alaska_north"),  # Chichagof Island
    (57.40, -136.00, "se_alaska_north"),  # Kruzof Island
    (57.05, -135.40, "se_alaska_north"),  # Sitka area (boundary)
    
    # SE Alaska South — outer coast / panhandle
    (56.80, -135.00, "se_alaska_south"),  # Baranof Island south
    (56.30, -134.50, "se_alaska_south"),  # Kuiu Island
    (55.90, -133.50, "se_alaska_south"),  # Prince of Wales west
    (55.50, -133.00, "se_alaska_south"),  # Prince of Wales south
    (55.30, -131.60, "se_alaska_south"),  # Ketchikan area
    (55.00, -130.50, "se_alaska_south"),  # Portland Canal approach
    (54.70, -130.40, "se_alaska_south"),  # AK/BC border
    
    # BC North Coast
    (54.30, -130.50, "bc_outer"),  # Prince Rupert area
    (53.80, -130.20, "bc_outer"),  # Porcher Island
    (53.30, -129.50, "bc_outer"),  # Pitt Island
    (53.00, -131.80, "bc_outer"),  # Haida Gwaii east
    (52.80, -131.50, "bc_outer"),  # Haida Gwaii south
    
    # BC Central Coast (fjord country)
    (52.50, -128.50, "bc_fjords"),  # Milbanke Sound
    (52.20, -128.20, "bc_fjords"),  # Near Burke Channel entrance
    (52.00, -128.10, "bc_fjords"),  # Calvert Island area
    (51.80, -128.00, "bc_fjords"),  # Cape Caution
    (51.60, -127.80, "bc_fjords"),  # Rivers Inlet approach
    (51.40, -127.50, "bc_fjords"),  # Smith Inlet area
    (51.10, -127.80, "bc_fjords"),  # Queen Charlotte Sound
    (50.80, -128.20, "bc_fjords"),  # Cape Scott (N Vancouver Is)
    
    # BC South / North Salish Sea approach
    (50.50, -128.00, "bc_outer"),   # NW Vancouver Island
    (50.20, -127.80, "bc_outer"),   # Quatsino Sound
    (49.80, -127.10, "bc_outer"),   # Nootka Sound
    (49.30, -126.20, "bc_outer"),   # Clayoquot Sound
    (48.80, -125.20, "bc_outer"),   # Barkley Sound
    (48.50, -124.80, "bc_outer"),   # West coast Vancouver Is south
    
    # Salish Sea (inside)
    (48.50, -124.60, "salish_sea"),  # Juan de Fuca entrance
    (48.40, -123.50, "salish_sea"),  # Victoria area
    (48.50, -123.00, "salish_sea"),  # San Juan Islands
    (48.70, -122.80, "salish_sea"),  # Bellingham area
    (49.00, -123.20, "salish_sea"),  # Vancouver area
    (49.30, -123.50, "salish_sea"),  # Howe Sound
    (49.60, -124.00, "salish_sea"),  # Texada Island
    (50.00, -124.80, "salish_sea"),  # Discovery Islands
    (50.30, -125.30, "salish_sea"),  # Johnstone Strait
    
    # WA Outer Coast
    (48.40, -124.70, "wa_or_outer"),  # Cape Flattery
    (48.10, -124.70, "wa_or_outer"),  # Olympic Coast N
    (47.70, -124.50, "wa_or_outer"),  # Quinault
    (47.30, -124.30, "wa_or_outer"),  # Grays Harbor
    (46.90, -124.10, "wa_or_outer"),  # Willapa Bay
    (46.50, -124.05, "wa_or_outer"),  # Long Beach
    (46.20, -124.00, "wa_or_outer"),  # Columbia River mouth
    
    # OR Coast
    (46.00, -123.95, "wa_or_outer"),  # Astoria
    (45.60, -123.95, "wa_or_outer"),  # Tillamook
    (45.00, -124.00, "wa_or_outer"),  # Lincoln City
    (44.60, -124.05, "wa_or_outer"),  # Newport
    (44.00, -124.10, "wa_or_outer"),  # Florence
    (43.40, -124.30, "wa_or_outer"),  # Coos Bay
    (42.80, -124.50, "wa_or_outer"),  # Gold Beach
    (42.05, -124.30, "wa_or_outer"),  # Brookings / OR-CA border
    
    # N California
    (41.75, -124.20, "n_california"),  # Crescent City
    (41.20, -124.15, "n_california"),  # Trinidad Head
    (40.80, -124.16, "n_california"),  # Eureka
    (40.40, -124.35, "n_california"),  # Cape Mendocino
    (39.80, -123.80, "n_california"),  # Fort Bragg
    (39.30, -123.80, "n_california"),  # Mendocino
    (38.95, -123.70, "n_california"),  # Point Arena
    
    # C California
    (38.35, -123.10, "c_california"),  # Bodega Bay
    (38.05, -123.00, "c_california"),  # Point Reyes
    (37.80, -122.50, "c_california"),  # San Francisco
    (37.50, -122.50, "c_california"),  # Half Moon Bay
    (37.00, -122.20, "c_california"),  # Santa Cruz
    (36.60, -121.90, "c_california"),  # Monterey
    (36.20, -121.80, "c_california"),  # Big Sur
    (35.70, -121.30, "c_california"),  # San Simeon
    (35.37, -120.85, "c_california"),  # Morro Bay
    (34.70, -120.60, "c_california"),  # Point Conception
    
    # S California
    (34.40, -119.85, "s_california"),  # Santa Barbara
    (34.05, -118.80, "s_california"),  # Ventura / Channel Is
    (33.95, -118.50, "s_california"),  # Malibu
    (33.74, -118.40, "s_california"),  # Palos Verdes
    (33.45, -117.80, "s_california"),  # Dana Point
    (33.20, -117.40, "s_california"),  # Oceanside
    (32.85, -117.27, "s_california"),  # La Jolla
    (32.55, -117.13, "s_california"),  # Border
]


# Known fjord locations (branches off main coast)
# Each is: (entrance_lat, entrance_lon, head_lat, head_lon, name, region, sill_depth, has_lens, n_sites)
FJORD_BRANCHES = [
    # SE Alaska North
    (58.90, -136.20, 59.00, -136.80, "Glacier Bay", "se_alaska_north", 20, True, 4),
    (58.30, -134.90, 58.30, -134.00, "Lynn Canal", "se_alaska_north", 300, True, 3),
    (57.78, -135.50, 57.78, -135.00, "Tenakee Inlet", "se_alaska_north", 30, True, 2),
    (57.40, -135.50, 57.50, -135.20, "Peril Strait", "se_alaska_north", 15, False, 2),
    (58.30, -134.50, 58.50, -134.20, "Chatham Strait N", "se_alaska_north", 200, False, 2),
    
    # SE Alaska South
    (55.50, -131.30, 55.60, -130.60, "Behm Canal/Misty Fjords", "se_alaska_south", 30, True, 3),
    (55.60, -132.00, 55.40, -131.50, "Clarence Strait", "se_alaska_south", 150, False, 2),
    
    # BC Fjords (the big ones)
    (52.20, -128.20, 52.15, -127.40, "Burke Channel", "bc_fjords", 40, True, 4),
    (52.35, -128.00, 52.37, -126.75, "Dean Channel/Bella Coola", "bc_fjords", 50, True, 5),
    (52.00, -128.10, 52.00, -127.50, "Fisher Channel", "bc_fjords", 60, False, 2),
    (51.65, -127.80, 51.65, -127.00, "Rivers Inlet", "bc_fjords", 35, True, 3),
    (51.40, -127.50, 51.40, -126.80, "Smith Inlet", "bc_fjords", 25, True, 2),
    (51.10, -127.80, 51.10, -127.10, "Belize Inlet", "bc_fjords", 20, False, 2),
    (53.50, -129.80, 53.50, -129.20, "Douglas Channel", "bc_fjords", 100, True, 3),
    (53.30, -129.00, 53.30, -128.20, "Gardner Canal", "bc_fjords", 40, True, 3),
    (55.00, -130.40, 55.00, -130.00, "Portland Canal", "bc_fjords", 30, True, 2),
    
    # Salish Sea fjords
    (49.30, -123.50, 49.60, -123.20, "Howe Sound", "salish_sea", 70, False, 2),
    (49.80, -124.00, 50.10, -123.70, "Jervis Inlet", "salish_sea", 200, True, 3),
    (47.70, -123.20, 47.40, -123.00, "Hood Canal", "salish_sea", 50, False, 3),
]


def interpolate_coast_sites(waypoints: list, spacing_km: float = 25.0) -> List[RealSite]:
    """
    Generate sites along coastline at regular intervals.
    
    Args:
        waypoints: List of (lat, lon, region) tuples
        spacing_km: Target distance between sites in km
    
    Returns:
        List of RealSite objects along the coast
    """
    sites = []
    site_id = 0
    
    for i in range(len(waypoints) - 1):
        lat1, lon1, region1 = waypoints[i]
        lat2, lon2, region2 = waypoints[i + 1]
        
        # Skip if crossing regions (different region = gap)
        dist = haversine_km(lat1, lon1, lat2, lon2)
        n_segments = max(1, int(dist / spacing_km))
        
        for s in range(n_segments):
            frac = s / n_segments
            lat = lat1 + frac * (lat2 - lat1)
            lon = lon1 + frac * (lon2 - lon1)
            region = region1  # Use source region
            
            # Determine base temperature from latitude
            base_temp = estimate_temperature(lat)
            
            sites.append(RealSite(
                name=f"Coast_{site_id:03d}",
                lat=lat,
                lon=lon,
                region=region,
                site_type="outer_coast",
                base_temp_C=base_temp,
            ))
            site_id += 1
    
    # Add final waypoint
    lat, lon, region = waypoints[-1]
    sites.append(RealSite(
        name=f"Coast_{site_id:03d}",
        lat=lat,
        lon=lon,
        region=region,
        site_type="outer_coast",
        base_temp_C=estimate_temperature(lat),
    ))
    
    return sites


def generate_fjord_sites(fjords: list) -> List[RealSite]:
    """Generate sites along fjord branches."""
    sites = []
    
    for entrance_lat, entrance_lon, head_lat, head_lon, name, region, sill_depth, has_lens, n_sites in fjords:
        for s in range(n_sites):
            frac = (s + 1) / (n_sites + 1)  # Evenly spaced, not at entrance
            lat = entrance_lat + frac * (head_lat - entrance_lat)
            lon = entrance_lon + frac * (head_lon - entrance_lon)
            
            # Inner sites have shallower effective sill (more isolated)
            effective_sill = sill_depth * (1 - 0.3 * frac)  # Gets shallower toward head
            
            base_temp = estimate_temperature(lat)
            if has_lens:
                base_temp -= 1.0  # Freshwater lens pushes to colder water
            
            site_name = f"{name} ({s+1})" if n_sites > 1 else name
            
            sites.append(RealSite(
                name=site_name,
                lat=lat,
                lon=lon,
                region=region,
                site_type="fjord",
                sill_depth_m=effective_sill,
                base_temp_C=base_temp,
                has_freshwater_lens=has_lens,
                notes=f"Part of {name} fjord system",
            ))
    
    return sites


def estimate_temperature(lat: float) -> float:
    """Estimate base bottom temperature from latitude."""
    # Simple linear model: colder at higher latitudes
    # SE Alaska (~58°N) ≈ 6°C, S. California (~33°N) ≈ 15°C
    return 6.0 + (58.0 - lat) * 0.36


def generate_full_site_network(coast_spacing_km: float = 25.0) -> List[RealSite]:
    """Generate complete site network: coast + fjords."""
    coast_sites = interpolate_coast_sites(COASTLINE_WAYPOINTS, coast_spacing_km)
    fjord_sites = generate_fjord_sites(FJORD_BRANCHES)
    
    # Mark coast sites near fjord entrances as "inland_sea" for Salish Sea
    for cs in coast_sites:
        if cs.region == "salish_sea":
            cs.site_type = "inland_sea"
    
    all_sites = coast_sites + fjord_sites
    return all_sites


if __name__ == "__main__":
    sites = generate_full_site_network(spacing_km=25.0)
    
    print(f"Total sites: {len(sites)}")
    
    # Count by region and type
    regions = {}
    for s in sites:
        key = s.region
        if key not in regions:
            regions[key] = {"total": 0, "fjord": 0, "coast": 0, "inland": 0, "lens": 0}
        regions[key]["total"] += 1
        if s.site_type == "fjord":
            regions[key]["fjord"] += 1
        elif s.site_type == "inland_sea":
            regions[key]["inland"] += 1
        else:
            regions[key]["coast"] += 1
        if s.has_freshwater_lens:
            regions[key]["lens"] += 1
    
    region_order = ["se_alaska_north", "se_alaska_south", "bc_fjords", "bc_outer",
                    "salish_sea", "wa_or_outer", "n_california", "c_california", "s_california"]
    
    print(f"\n{'Region':<22} {'Total':>5} {'Fjord':>5} {'Coast':>5} {'Inland':>6} {'FW Lens':>7}")
    print("-" * 55)
    for r in region_order:
        if r in regions:
            d = regions[r]
            print(f"{r:<22} {d['total']:>5} {d['fjord']:>5} {d['coast']:>5} {d['inland']:>6} {d['lens']:>7}")
    
    print(f"\nLatitude range: {min(s.lat for s in sites):.1f}°N to {max(s.lat for s in sites):.1f}°N")
    print(f"Temperature range: {min(s.base_temp_C for s in sites):.1f}°C to {max(s.base_temp_C for s in sites):.1f}°C")
