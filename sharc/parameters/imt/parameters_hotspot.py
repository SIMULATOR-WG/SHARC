from dataclasses import dataclass

from sharc.parameters.parameters_base import ParametersBase


@dataclass
class ParametersHotspot(ParametersBase):
    intersite_distance: int = None
    # Enable wrap around
    wrap_around: bool = False
    # Number of clusters topology
    num_clusters: int = 1
    # Number of hotspots per macro cell (sector)
    num_hotspots_per_cell: float = 1.0
    # Maximum 2D distance between hotspot and UE [m]
    # This is the hotspot radius
    max_dist_hotspot_ue: float = 100.0
    # Minimum 2D distance between macro cell base station and hotspot [m]
    min_dist_bs_hotspot: float = 0.0

    def validate(self, ctx):
        if not isinstance(self.intersite_distance, int) and not isinstance(self.intersite_distance, float):
            raise ValueError(f"{ctx}.intersite_distance should be a number")

        if self.num_clusters not in [1, 7]:
            raise ValueError(f"{ctx}.num_clusters should be either 1 or 7")

        if not isinstance(self.num_hotspots_per_cell, int) or self.num_hotspots_per_cell < 0:
            raise ValueError("num_hotspots_per_cell must be non-negative")

        if (not isinstance(self.max_dist_hotspot_ue, float) and not isinstance(self.max_dist_hotspot_ue, int))\
                or self.max_dist_hotspot_ue < 0:
            raise ValueError("max_dist_hotspot_ue must be non-negative")

        if (not isinstance(self.min_dist_bs_hotspot, float) and not isinstance(self.min_dist_bs_hotspot, int))\
                or self.min_dist_bs_hotspot < 0:
            raise ValueError("min_dist_bs_hotspot must be non-negative")
