from dataclasses import dataclass, field
from typing import List
import numpy as np
from sharc.parameters.parameters_base import ParametersBase
from sharc.parameters.parameters_orbit import ParametersOrbit
from sharc.parameters.antenna.parameters_antenna_s1528 import ParametersAntennaS1528
from sharc.parameters.parameters_p619 import ParametersP619


@dataclass
class ParametersNgsoConstellation(ParametersBase):
    """
    Combines parameters for NGSO constellations and Mobile Satellite Systems (MSS).
    """

    # General configuration for NGSO constellation
    section_name: str = "ngso_mss"  # Section name in configuration file
    name: str = "Default"  # Name of the NGSO constellation

    # Orbital configuration
    orbits: List[ParametersOrbit] = field(default_factory=list)  # List of orbit parameters

    # Antenna configuration
    antenna: str = None  # Antenna type for the constellation
    max_transmit_power_dB: float = 46.0  # Maximum transmit power in dB
    max_transmit_gain_dBi: float = 30.0  # Maximum transmit gain in dBi
    max_receive_gain_dBi: float = 30.0  # Maximum receive gain in dBi
    max_num_of_beams: int = 19  # Maximum number of beams per satellite
    cell_radius_m: float = 19000.0  # Cell radius in meters
    antenna_s1528: ParametersAntennaS1528 = field(default_factory=ParametersAntennaS1528)  # Antenna parameters

    # MSS-specific parameters
    is_space_to_earth: bool = True  # Direction of communication (space-to-earth or earth-to-space)
    x: float = 0.0  # X-coordinate of the satellite bore-sight
    y: float = 0.0  # Y-coordinate of the satellite bore-sight
    frequency: float = 2110.0  # System center frequency in MHz
    bandwidth: float = 20.0  # System bandwidth in MHz
    spectral_mask: str = "3GPP E-UTRA"  # Transmitter spectral mask type
    spurious_emissions: float = -13  # Out-of-band spurious emissions in dB/MHz
    adjacent_ch_leak_ratio: float = 45.0  # Adjacent channel leakage ratio in dB
    altitude: float = 1200000.0  # Satellite altitude above sea level in meters
    intersite_distance: float = field(init=False)  # Inter-site distance, calculated from cell radius
    tx_power_density: float = 40.0  # Satellite transmit power density in dBW/MHz
    azimuth: float = 45.0  # Azimuth angle in degrees
    elevation: float = 90.0  # Elevation angle in degrees
    num_sectors: int = 7  # Number of sectors served by the satellite
    antenna_pattern: str = "ITU-R-S.1528-LEO"  # Antenna pattern type
    antenna_diameter: float = 1.0  # Antenna diameter in meters
    antenna_l_s: float = -6.75  # Near-in-side-lobe level in dB relative to peak gain
    antenna_3_dB_bw: float = 4.4127  # 3 dB beamwidth angle in degrees

    # Channel model configuration
    param_p619: ParametersP619 = field(default_factory=ParametersP619)  # Parameters for the P.619 channel model
    space_station_alt_m: float = 35780000.0  # Altitude of the space station in meters
    earth_station_alt_m: float = 0.0  # Altitude of the earth station in meters
    earth_station_lat_deg: float = 0.0  # Latitude of the earth station in degrees
    earth_station_long_diff_deg: float = 0.0  # Longitude difference for the earth station in degrees
    season: str = "SUMMER"  # Season for atmospheric conditions
    channel_model: str = "P619"  # Channel model type (e.g., FSPL, SatelliteSimple, P619)

    def __post_init__(self):
        """
        Perform initialization calculations.
        """
        # Calculate inter-site distance from cell radius
        self.intersite_distance = np.sqrt(3) * self.cell_radius_m

    def load_parameters_from_file(self, config_file: str):
        """
        Load parameters from a configuration file and validate.

        Parameters
        ----------
        config_file : str
            Path to the configuration file.

        Raises
        ------
        ValueError
            If a parameter fails validation.
        """
        super().load_parameters_from_file(config_file)

        # Validate NGSO constellation parameters
        if not self.name or not isinstance(self.name, str):
            raise ValueError(f"ParametersNgsoMss: Invalid name = {self.name}. Must be a non-empty string.")

        if not self.orbits:
            raise ValueError("ParametersNgsoMss: No orbits defined. At least one orbit must be specified.")

        for orbit in self.orbits:
            if not isinstance(orbit, ParametersOrbit):
                raise ValueError("ParametersNgsoMss: Invalid orbit configuration. Must be instances of ParametersOrbit.")

        # Validate MSS-specific parameters
        if self.num_sectors not in [1, 7, 19]:
            raise ValueError(f"ParametersNgsoMss: Invalid number of sectors: {self.num_sectors}")

        if self.cell_radius_m <= 0:
            raise ValueError(f"ParametersNgsoMss: cell_radius must be greater than 0, but is {self.cell_radius_m}")
        else:
            self.intersite_distance = np.sqrt(3) * self.cell_radius_m

        if not (0 <= self.azimuth <= 360):
            raise ValueError("ParametersNgsoMss: azimuth must be between 0 and 360 degrees")

        if not (0 <= self.elevation <= 90):
            raise ValueError("ParametersNgsoMss: elevation must be between 0 and 90 degrees")

        if self.channel_model.upper() not in ["FSPL", "P619", "SATELLITESIMPLE"]:
            raise ValueError(f"ParametersNgsoMss: Invalid channel model name {self.channel_model}")

        if self.channel_model == "P619":
            self.param_p619.set_external_parameters(
                space_station_alt_m=self.altitude,
                earth_station_alt_m=self.earth_station_alt_m,
                earth_station_lat_deg=self.earth_station_lat_deg,
                earth_station_long_diff_deg=self.earth_station_long_diff_deg,
                season=self.season
            )


if __name__ == "__main__":

    try:
        # Create orbital parameters for the first orbit
        orbit_1 = ParametersOrbit(
            n_planes=20,                  # Number of orbital planes
            sats_per_plane=32,            # Satellites per plane
            phasing_deg=3.9,              # Phasing angle in degrees
            long_asc_deg=18.0,            # Longitude of ascending node
            inclination_deg=54.5,         # Orbital inclination in degrees
            perigee_alt_km=525.0,         # Perigee altitude in kilometers
            apogee_alt_km=525.0           # Apogee altitude in kilometers
        )

        # Create orbital parameters for the second orbit
        orbit_2 = ParametersOrbit(
            n_planes=12,                  # Number of orbital planes
            sats_per_plane=20,            # Satellites per plane
            phasing_deg=2.0,              # Phasing angle in degrees
            long_asc_deg=30.0,            # Longitude of ascending node
            inclination_deg=26.0,         # Orbital inclination in degrees
            perigee_alt_km=580.0,         # Perigee altitude in kilometers
            apogee_alt_km=580.0           # Apogee altitude in kilometers
        )

        # Instantiate the ParametersNgsoConstellation class with manually defined attributes
        constellation = ParametersNgsoConstellation(
            name="Acme-Star-1",           # Name of the constellation
            antenna="Taylor1.4",          # Antenna type
            max_transmit_power_dB=46.0,   # Maximum transmit power in dB
            max_transmit_gain_dBi=30.0,   # Maximum transmit gain in dBi
            max_receive_gain_dBi=30.0,    # Maximum receive gain in dBi
            max_num_of_beams=19,          # Maximum number of beams
            cell_radius_m=19000.0,        # Cell radius in meters
            orbits=[orbit_1, orbit_2]     # List of orbital parameters
        )

        # Display the manually set parameters
        print(f"\n#### Constellation Details ####")
        print(f"Name: {constellation.name}")
        print(f"Antenna Type: {constellation.antenna}")
        print(f"Maximum Transmit Power (dB): {constellation.max_transmit_power_dB}")
        print(f"Maximum Transmit Gain (dBi): {constellation.max_transmit_gain_dBi}")
        print(f"Maximum Receive Gain (dBi): {constellation.max_receive_gain_dBi}")
        print(f"Maximum Number of Beams: {constellation.max_num_of_beams}")
        print(f"Cell Radius (m): {constellation.cell_radius_m}")

        # Iterate through the orbits and display their parameters
        print("\n#### Orbital Parameters ####")
        for idx, orbit in enumerate(constellation.orbits, start=1):
            print(f"Orbit {idx}:")
            print(f"  Number of Orbital Planes: {orbit.n_planes}")
            print(f"  Satellites per Plane: {orbit.sats_per_plane}")
            print(f"  Phasing Angle (deg): {orbit.phasing_deg}")
            print(f"  Longitude of Ascending Node (deg): {orbit.long_asc_deg}")
            print(f"  Orbital Inclination (deg): {orbit.inclination_deg}")
            print(f"  Perigee Altitude (km): {orbit.perigee_alt_km}")
            print(f"  Apogee Altitude (km): {orbit.apogee_alt_km}")
    except ValueError as error:
        print(f"Validation Error: {error}")
    except Exception as unexpected_error:
        print(f"An unexpected error occurred: {unexpected_error}")
