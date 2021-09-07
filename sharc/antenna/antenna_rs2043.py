# -*- coding: utf-8 -*-
"""
@Created: Luciano Camilo on Tue Aug 10 19:42:00 2021

"""

import numpy as np
from sharc.parameters.parameters_fs import ParametersFs


class AntennaRS2043(object):
    """
     Implements the reference antenna pattern described in Table 9 from Recommendation ITU-R RS.2043.

    Attributes
    ----------
        antenna_gain (float): maximum gain of FS omni antenna
    """

    def __init__(self, param: ParametersFs):
        """
        Constructs an AntennaElement object.

        Parameters
        ---------
            param (ParametersFS): antenna EESS parameters
        """
        super().__init__()

    def calculate_gain(self, **kwargs) -> np.array:
        phi = np.asarray(kwargs["off_axis_angle_vec"])
        theta = np.asarray(kwargs["theta_vec"])
        gain = self.get_gain_az(phi) + self.get_gain_el(theta)

        return gain

    def get_gain_az(self, phi: np.array) -> np.array:
        """
        Returns the antenna gain in the azimuth plane for the given direction.
        """
        gh = np.zeros(phi.shape)

        for i in range(len(phi)):
            if -0.542<= phi[i] <= 0.542:
                gh[i] = 0 - 45.53 * ((phi[i]) ** 2)

            if 0.542 < phi[i] <= 5.053:
                gh[i] = -11.210 - 4.0220 * phi[i]

            if -0.542 > phi[i] >= -5.053:
                gh[i] = -11.210 + 4.0220 * phi[i]

            if 5.053 < phi[i] <= 14.708:
                gh[i] = -26.720 - 0.9530 * phi[i]

            if -5.053 > phi[i] >= -14.708:
                gh[i] = -26.720 + 0.9530 * phi[i]

            if 14.708 < phi[i] <= 30:
                gh[i] = -35.031 - 0.3880 * phi[i]

            if -14.708 > phi[i] >= -30:
                gh[i] = -35.031 + 0.3880 * phi[i]

            if 30 < phi[i] <= 59.915:
                gh[i] = -41.836 - 0.1580 * phi[i]

            if -30 > phi[i] >= -59.915:
                gh[i] = -41.836 + 0.1580 * phi[i]

            if phi[i] > 59.915:
                gh[i] = -51.387

            if phi[i] < -59.915:
                gh[i] = -51.387

        return gh

    def get_gain_el(self, theta: np.array) -> np.array:
        """
        Returns the antenna gain in the elevation plane for the given direction.
        """
        gv = np.zeros(len(theta))

        for i in range(len(theta)):
            if -1.149 < theta[i] < 1.149:
                gv[i] = 47 - 9.91 * ((theta[i]) ** 2)

            if 1.149 <= theta[i] <= 9.587:
                gv[i] = 35.189 - 1.9440 * theta[i]

            if -1.149 >= theta[i] >= -9.587:
                gv[i] = 35.189 + 1.9440 * theta[i]

            if 9.587 <= theta[i] <= 29.976:
                gv[i] = 21.043 - 0.4680 * theta[i]

            if -9.587 >= theta[i] >= -29.976:
                gv[i] = 21.043 + 0.4680 * theta[i]

            if 29.976 <= theta[i] <= 50:
                gv[i] = 12.562 - 0.1850 * theta[i]

            if -29.976 >= theta[i] >= -50:
                gv[i] = 12.562 + 0.1850 * theta[i]

            if theta[i] > 50:
                gv[i] = 3.291

            if theta[i] < -50:
                gv[i] = 3.291

        return gv

if __name__ == '__main__':

    from sharc.parameters.parameters_eess_passive import ParametersEessPassive
    from matplotlib import pyplot as plt
    """
    Test routine - Comparison between SHARC code and Table 9 of ITU R - RS.2043-0
    """

    param = ParametersEessPassive()
    antenna = AntennaRS2043(param)
    # phi = np.arange(-90, 90, step=0.01)
    # theta = np.arange(-40, 40, step=0.01)
    phi = np.arange(0, 180, step=0.01)
    theta = np.arange(0, 90, step=0.01)
    csfont = {'fontname': 'Times New Roman'}
    hfont = {'fontname': 'Times New Roman'}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))

    ax1.semilogx(theta, antenna.get_gain_el(theta=theta), '-',color='blue', label='ITU-R RS.2043-0 (Table 9)')
    ax2.semilogx(phi, antenna.get_gain_az(phi=phi), '-', color='darkorange', label='ITU-R RS.2043-0 (Table 9)')

    # ax1.plot(theta, antenna.get_gain_el(theta=theta), '-', color='blue', label='ITU-R RS.2043-0 (Table 9)')
    # ax2.plot(phi, antenna.get_gain_az(phi=phi)+47, '-', color='darkorange', label='ITU-R RS.2043-0 (Table 9)')
    ax1.grid()
    ax2.grid()
    ax1.legend()
    ax2.legend()

    ax1.set_title('ITU-R RS.2043-0 - Vertical Antenna Pattern')
    ax2.set_title('ITU-R RS.2043-0 - Horizontal Antenna Pattern')
    ax1.set_xlabel('Elevation angle (degrees)', fontsize=12, color='black')
    ax1.set_ylabel('Gain (dBi)', fontsize=12, color='black')
    ax2.set_xlabel('Azimuth (degrees)', fontsize=12, color='black')
    ax2.set_ylabel('Gain (dBi)', fontsize=12, color='black')

    plt.show()
