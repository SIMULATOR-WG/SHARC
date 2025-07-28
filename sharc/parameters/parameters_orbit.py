from dataclasses import dataclass, field
from sharc.parameters.parameters_base import ParametersBase


@dataclass
class ParametersOrbit(ParametersBase):
    n_planes: int = 8
    sats_per_plane: int = 6
    phasing_deg: float = 7.5
    long_asc_deg: float = 0.0
    omega_deg: float = 0.0
    inclination_deg: float = 52.0
    perigee_alt_km: float = 1414.0
    apogee_alt_km: float = 1414.0
    initial_mean_anomaly: float = 0.0

    def load_parameters_from_file(self, config_file: str):
        """Load parameters from file and validate."""
        super().load_parameters_from_file(config_file)

        self.validate("ParametersOrbit")

    def validate(self, ctx: str):
        if not (0 <= self.inclination_deg <= 180):
            raise ValueError(f"Invalid {ctx}.inclination_deg = {self.inclination_deg}. \
                             Must be in the range [0, 180] degrees.")
        if self.perigee_alt_km < 0:
            raise ValueError(f"Invalid {ctx}.perigee_alt_km = {self.perigee_alt_km}. \
                             Altitude must be non-negative.")
        if self.apogee_alt_km < self.perigee_alt_km:
            raise ValueError(f"Invalid {ctx}.apogee_alt_km = {self.apogee_alt_km}. \
                             Must be greater than or equal to perigee_alt_km.")
        if not (0 <= self.phasing_deg <= 360):
            raise ValueError(f"Invalid {ctx}.phasing_deg = {self.phasing_deg}. \
                             Must be in the range [0, 360] degrees.")
