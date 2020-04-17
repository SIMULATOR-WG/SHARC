from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_fs import ParametersFs

import numpy as np


class AntennaTacticalRadioRelays(Antenna):
    """
    This class implements the antenna pattern described in the report "THE POTENTIAL FOR ACCOMMODATING THIRD GENERATION
    MOBILE SYSTEMS IN THE 1710–1850 MHZ BAND: Federal Operations, Relocation Costs, and Operational Impacts", by the US
    National Telecommunications and Information Administration, Table 4-4 (AN/GRC-226 (MSE) Parameters).
    """
    def __init__(self, param: ParametersFs):
        super().__init__()

        self.gain_levels = np.array([20, 11, 2])
        self.angle_breaks = np.array([20, 90])

    def calculate_gain(self, *args, **kwargs) -> np.array:
        phi = np.absolute(kwargs["off_axis_angle_vec"])

        gain = np.zeros(phi.shape)

        idx_0 = np.where(phi < self.angle_breaks[0])[0]
        gain[idx_0] = self.gain_levels[0]

        idx_1 = np.where((self.angle_breaks[0] <= phi) & (phi < self.angle_breaks[1]))[0]
        gain[idx_1] = self.gain_levels[1]

        idx_2 = np.where((self.angle_breaks[1] <= phi) & (phi <= 180))[0]
        gain[idx_2] = self.gain_levels[2]

        return gain



if __name__ == '__main__':
    import matplotlib.pyplot as plt

    phi = np.linspace(0.1, 180, num=100000)
    antenna = AntennaTacticalRadioRelays(None)

    gain = antenna.calculate_gain(off_axis_angle_vec=phi)

    plt.plot(phi, gain)
    plt.title("AN/GRC-226 (MSE) Tactical Radio Relay antenna pattern")
    plt.xlabel("Off-axis angle $\phi$ [deg]")
    plt.ylabel("Gain [dBi]")
    plt.xlim((phi[0], phi[-1]))
    plt.ylim((0, 22))
    plt.grid()
    plt.show()
