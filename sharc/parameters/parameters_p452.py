# -*- coding: utf-8 -*-
from dataclasses import dataclass

from sharc.parameters.parameters_base import ParametersBase


@dataclass
class ParametersP452(ParametersBase):
    """Dataclass containing the P.452 propagation model parameters
    """
    # Total air pressure in hPa
    atmospheric_pressure: float = 935.0
    # Temperature in Kelvin
    air_temperature: float = 300.0
    # Sea-level surface refractivity (use the map)
    N0: float = 352.58
    # Average radio-refractive (use the map)
    delta_N: float = 43.127
    # percentage p. Float (0 to 100) or RANDOM
    percentage_p: float = 0.2
    # Distance over land from the transmit and receive antennas to the coast (km)
    Dct: float = 70.0
    # Distance over land from the transmit and receive antennas to the coast (km)
    Dcr: float = 70.0
    # Effective height of interfering antenna (m)
    Hte: float = 20.0
    # Effective height of interfered-with antenna (m)
    Hre: float = 3.0
    # Latitude of transmitter
    tx_lat: float = -23.55028
    # Latitude of receiver
    rx_lat: float = -23.17889
    # Antenna polarization
    polarization: str = "horizontal"
    # determine whether clutter loss following ITU-R P.2108 is added (TRUE/FALSE)
    clutter_loss: bool = True

    def load_from_paramters(self, param: ParametersBase):
        """Used to load parameters of P.452 from IMT or system parameters

        Parameters
        ----------
        param : ParametersBase
            IMT or system parameters
        """
        self.atmospheric_pressure = param.atmospheric_pressure
        self.air_temperature = param.air_temperature
        self.N0 = param.N0
        self.delta_N = param.delta_N
        self.percentage_p = param.percentage_p
        self.Dct = param.Dct
        self.Dcr = param.Dcr
        self.Hte = param.Hte
        self.Hre = param.Hre
        self.tx_lat = param.tx_lat
        self.rx_lat = param.rx_lat
        self.polarization = param.polarization
        self.clutter_loss = param.clutter_loss
