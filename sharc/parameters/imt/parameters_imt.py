# -*- coding: utf-8 -*-
"""Parameters definitions for IMT systems
"""
from dataclasses import dataclass, field
import typing

from sharc.parameters.parameters_base import ParametersBase
from sharc.parameters.parameters_p619 import ParametersP619
from sharc.parameters.imt.parameters_antenna_imt import ParametersAntennaImt
from sharc.parameters.parameters_antenna import ParametersAntenna
from sharc.parameters.imt.parameters_imt_topology import ParametersImtTopology


@dataclass
class ParametersImt(ParametersBase):
    """Dataclass containing the IMT system parameters
    """
    section_name: str = "imt"
    # whether to enable recursive parameters setting on .yaml file
    nested_parameters_enabled: bool = True

    minimum_separation_distance_bs_ue: float = 0.0
    interfered_with: bool = False
    frequency: float = 24350.0
    bandwidth: float = 200.0
    rb_bandwidth: float = 0.180
    spectral_mask: str = "IMT-2020"
    spurious_emissions: float = -13.0
    guard_band_ratio: float = 0.1
    # Adjacent Interference filter reception used when IMT is victim. Possible values is ACS and OFF
    adjacent_ch_reception: str = "OFF"

    @dataclass
    class ParametersBS(ParametersBase):
        """Dataclass containing the IMT Base Station (BS) parameters."""

        load_probability = 0.2
        conducted_power = 10.0
        height: float = 6.0
        noise_figure: float = 10.0
        ohmic_loss: float = 3.0
        antenna: ParametersAntenna = field(default_factory=lambda: ParametersAntenna(
            pattern="ARRAY", array=ParametersAntennaImt(downtilt=0.0)
        ))
    bs: ParametersBS = field(default_factory=ParametersBS)

    topology: ParametersImtTopology = field(default_factory=ParametersImtTopology)

    @dataclass
    class ParametersUL(ParametersBase):
        """Dataclass containing the IMT Uplink (UL) parameters."""

        attenuation_factor: float = 0.4
        sinr_min: float = -10.0
        sinr_max: float = 22.0
    uplink: ParametersUL = field(default_factory=ParametersUL)

    # Antenna model for adjacent band studies.
    adjacent_antenna_model: typing.Literal["SINGLE_ELEMENT", "BEAMFORMING"] = "SINGLE_ELEMENT"

    @dataclass
    class ParametersUE(ParametersBase):
        """Dataclass containing the IMT User Equipment (UE) parameters."""

        k: int = 3
        k_m: int = 1
        indoor_percent: int = 5.0
        distribution_type: str = "ANGLE_AND_DISTANCE"
        distribution_distance: str = "RAYLEIGH"
        distribution_azimuth: str = "NORMAL"
        azimuth_range: tuple = (-60, 60)
        tx_power_control: bool = True
        p_o_pusch: float = -95.0
        alpha: float = 1.0
        p_cmax: float = 22.0
        power_dynamic_range: float = 63.0
        height: float = 1.5
        noise_figure: float = 10.0
        ohmic_loss: float = 3.0
        body_loss: float = 4.0
        adjacent_ch_selectivity: float = 33  # Adjacent Channel Selectivity in dB
        antenna: ParametersAntenna = field(default_factory=lambda: ParametersAntenna(
            pattern="ARRAY"
        ))

        def validate(self, ctx: str):
            """Validate the UE antenna beamsteering range parameters."""
            if self.antenna.array.horizontal_beamsteering_range != (-180., 180.)\
                    or self.antenna.array.vertical_beamsteering_range != (0., 180.):
                raise NotImplementedError(
                    "UE antenna beamsteering limit has not been implemented. Default values of\n"
                    "horizontal = (-180., 180.), vertical = (0., 180.) should not be changed"
                )

    ue: ParametersUE = field(default_factory=ParametersUE)

    @dataclass
    class ParamatersDL(ParametersBase):
        attenuation_factor: float = 0.6
        sinr_min: float = -10.0
        sinr_max: float = 30.0

    downlink: ParamatersDL = field(default_factory=ParamatersDL)

    noise_temperature: float = 290.0
    # Channel parameters
    # channel model, possible values are "FSPL" (free-space path loss),
    #                                    "CI" (close-in FS reference distance)
    #                                    "UMa" (Urban Macro - 3GPP)
    #                                    "UMi" (Urban Micro - 3GPP)
    #                                    "TVRO-URBAN"
    #                                    "TVRO-SUBURBAN"
    #                                    "ABG" (Alpha-Beta-Gamma)
    # TODO: check if we wanna separate the channel model definition in its own nested attributes
    channel_model: str = "UMi"
    # Parameters for the P.619 propagation model
    # For IMT NTN the model is used for calculating the coupling loss between
    # the BS space station and the UEs on Earth's surface.
    # For now, the NTN footprint is centered over the BS nadir point, therefore
    # the paramters imt_lag_deg and imt_long_diff_deg SHALL be zero.
    #    space_station_alt_m - altitude of IMT space station (meters)
    #    earth_station_alt_m - altitude of IMT earth stations (UEs) (in meters)
    #    earth_station_lat_deg - latitude of IMT earth stations (UEs) (in degrees)
    #    earth_station_long_diff_deg - difference between longitudes of IMT space and earth stations
    #      (positive if space-station is to the East of earth-station)
    #    season - season of the year.
    param_p619 = ParametersP619()
    season: str = "SUMMER"

    # TODO: create parameters for where this is needed
    los_adjustment_factor: float = 18.0
    shadowing: bool = True

    def load_parameters_from_file(self, config_file: str):
        """Load the parameters from file an run a sanity check

        Parameters
        ----------
        file_name : str
            the path to the configuration file

        Raises
        ------
        ValueError
            if a parameter is not valid
        """
        super().load_parameters_from_file(config_file)

        if self.spectral_mask not in ["IMT-2020", "3GPP E-UTRA"]:
            raise ValueError(
                f"""ParametersImt: Inavlid Spectral Mask Name {self.spectral_mask}""",
            )
        if self.adjacent_ch_reception not in ["ACS", "OFF"]:
            raise ValueError(
                f"""ParametersImt: Invalid Adjacent Channel Reception model {self.adjacent_ch_reception}""",
            )

        if self.channel_model not in ["FSPL", "CI", "UMa", "UMi", "TVRO-URBAN", "TVRO-SUBURBAN", "ABG", "P619"]:
            raise ValueError(f"ParamtersImt: \
                             Invalid value for parameter channel_model - {self.channel_model}. \
                             Possible values are \"FSPL\",\"CI\", \"UMa\", \"UMi\", \"TVRO-URBAN\", \"TVRO-SUBURBAN\", \
                             \"ABG\", \"P619\".")

        if self.topology.type == "NTN" and self.channel_model not in ["FSPL", "P619"]:
            raise ValueError(
                f"ParametersImt: Invalid channel model {self.channel_model} for topology NTN",
            )

        if self.season not in ["SUMMER", "WINTER"]:
            raise ValueError(f"ParamtersImt: \
                             Invalid value for parameter season - {self.season}. \
                             Possible values are \"SUMMER\", \"WINTER\".")

        if self.topology.type == "NTN":
            self.is_space_to_earth = True
        #     self.param_p619.load_from_paramters(self)

        self.frequency = float(self.frequency)

        self.topology.ntn.set_external_parameters(
            bs_height=self.bs.height
        )

        self.bs.antenna.set_external_parameters(
            adjacent_antenna_model=self.adjacent_antenna_model,
            frequency=self.frequency,
            bandwidth=self.bandwidth,
        )

        self.ue.antenna.set_external_parameters(
            adjacent_antenna_model=self.adjacent_antenna_model,
            frequency=self.frequency,
            bandwidth=self.bandwidth,
        )

        self.validate("imt")
