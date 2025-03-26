# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 14:49:01 2017

@author: edgar
"""
import sys
from sharc.antenna.antenna import Antenna
from sharc.parameters.antenna.parameters_antenna_s1528 import ParametersAntennaS1528
from sharc.parameters.constants import SPEED_OF_LIGHT
import math
import numpy as np
from scipy.special import jn, jn_zeros


class AntennaS1528Taylor(Antenna):
    """
    Implements Recommendation ITU-R S.1528-0 Section 1.4: Satellite antenna reference pattern given by an analytical
    function which models the side lobes of the non-GSO satellite operating in the fixed-satellite service below 30 GHz.
    It is proposed to use a circular Taylor illumination function which gives the maximum flexibility to
    adapt the theoretical pattern to the real one. It takes into account the side-lobes effect of an antenna
    diagram.
    """

    def __init__(self, param: ParametersAntennaS1528):
        # Gmax
        self.peak_gain = param.antenna_gain
        self.frequency_mhz = param.frequency
        self.bandwidth_mhz = param.bandwidth
        # Wavelength of the lowest frequency of the band of interest (in meters).
        self.lamb = (SPEED_OF_LIGHT / 1e6) / (self.frequency_mhz - self.bandwidth_mhz / 2)
        # SLR is the side-lobe ratio of the pattern (dB), the difference in gain between the maximum
        # gain and the gain at the peak of the first side lobe.
        self.slr = param.slr
        # Number of secondary lobes considered in the diagram (coincide with the roots of the Bessel function)
        self.n_side_lobes = param.n_side_lobes
        # Radial and transverse sizes of the effective radiating area of the satellite transmit antenna (m).
        self.l_r = param.l_r
        self.l_t = param.l_t

        self.roll_off = param.roll_off

        # Beam roll-off (difference between the maximum gain and the gain at the edge of the illuminated beam)
        # Possible values are 3, 5, and 7
        if int(self.roll_off) not in [3, 5, 7]:
            raise ValueError(
                f"AntennaS1528Taylor: Invalid value for roll_off factor {self.roll_off}")
        self.roll_off = int(self.roll_off)

    def calculate_gain(self, *args, **kwargs) -> np.array:
        phi = np.abs(np.radians(kwargs.get('phi', 0)))
        theta = np.abs(np.radians(kwargs.get('theta', 0)))

        # Intermediary variables
        A = (1 / np.pi) * np.arccosh(10 ** (self.slr / 20))
        j1_roots = jn_zeros(1, self.n_side_lobes) / np.pi
        sigma = j1_roots[-1] / np.sqrt(A ** 2 + (self.n_side_lobes - 1 / 2) ** 2)
        u = (np.pi / self.lamb) * np.sqrt((self.l_r * np.sin(theta) * np.cos(phi)) ** 2 +
                                          (self.l_t * np.sin(theta) * np.sin(phi)) ** 2)

        mu = jn_zeros(1, self.n_side_lobes - 1) / np.pi
        v = np.ones(u.shape + (self.n_side_lobes - 1,))

        for i, ui in enumerate(mu):
            v[..., i] = (1 - u ** 2 / (np.pi ** 2 * sigma ** 2 *
                                       (A ** 2 + (i + 1 - 0.5) ** 2))) / (1 - (u / (np.pi * ui)) ** 2)

        # Take care of divide-by-zero
        with np.errstate(divide='ignore', invalid='ignore'):
            gain = self.peak_gain + 20 * \
                np.log10(np.abs((2 * jn(1, u) / u) * np.prod(v, axis=-1)))

        # Replace undefined values with -inf (or other desired value)
        gain = np.nan_to_num(gain, nan=-np.inf)

        return gain


class AntennaS1528Leo(Antenna):
    """
    Implements Recommendation ITU-R S.1528-0 Section 1.3 - LEO: Satellite antenna radiation
    patterns for LEO orbit satellite antennas operating in the
    fixed-satellite service below 30 GHz.
    """

    def __init__(self, param: ParametersAntennaS1528):
        super().__init__()
        self.peak_gain = param.antenna_gain
        self.psi_b = param.antenna_3_dB / 2
        # near-in-side-lobe level (dB) relative to the peak gain required by
        # the system design
        self.l_s = -6.75
        # for elliptical antennas, this is the ratio major axis/minor axis
        # we assume circular antennas, so z = 1
        # self.z = 1
        # far-out side-lobe level [dBi]
        self.l_f = 5
        self.y = 1.5 * self.psi_b
        self.z = self.y * np.power(10, 0.04 * (self.peak_gain + self.l_s - self.l_f))

    def calculate_gain(self, *args, **kwargs) -> np.array:
        """
        Calculates the gain in the given direction.

        Parameters
        ----------
            off_axis_angle_vec (np.array): azimuth angles (phi_vec) [degrees]

        Returns
        -------
            gain (np.array): gain corresponding to each of the given directions
        """
        psi = np.absolute(kwargs["off_axis_angle_vec"])
        gain = np.zeros(len(psi))

        idx_0 = np.where((psi <= self.y))[0]
        gain[idx_0] = self.peak_gain - 3 * \
            np.power(psi[idx_0] / self.psi_b, 2)

        idx_1 = np.where((self.y < psi) & (psi <= self.z))[0]
        gain[idx_1] = self.peak_gain + self.l_s - \
            25 * np.log10(psi[idx_1] / self.y)

        idx_3 = np.where((self.z < psi) & (psi <= 180))[0]
        gain[idx_3] = self.l_f

        return gain


class AntennaS1528(Antenna):
    """
    Implements Recommendation ITU-R S.1528-0: Satellite antenna radiation
    patterns for non-geostationary orbit satellite antennas operating in the
    fixed-satellite service below 30 GHz.
    This implementation refers to the pattern described in S.1528-0 Section 1.2
    """

    def __init__(self, param: ParametersAntennaS1528):
        super().__init__()
        self.peak_gain = param.antenna_gain

        # near-in-side-lobe level (dB) relative to the peak gain required by
        # the system design
        self.l_s = param.antenna_l_s

        # for elliptical antennas, this is the ratio major axis/minor axis
        # we assume circular antennas, so z = 1
        self.z = 1

        # far-out side-lobe level [dBi]
        self.l_f = 0

        # back-lobe level
        self.l_b = np.maximum(
            0, 15 + self.l_s + 0.25 *
            self.peak_gain + 5 * math.log10(self.z),
        )
        # one-half the 3 dB beamwidth in the plane of interest
        self.psi_b = param.antenna_3_dB / 2

        if self.l_s == -15:
            self.a = 2.58 * math.sqrt(1 - 1.4 * math.log10(self.z))
        elif self.l_s == -20:
            self.a = 2.58 * math.sqrt(1 - 1.0 * math.log10(self.z))
        elif self.l_s == -25:
            self.a = 2.58 * math.sqrt(1 - 0.6 * math.log10(self.z))
        elif self.l_s == -30:
            self.a = 2.58 * math.sqrt(1 - 0.4 * math.log10(self.z))
        else:
            sys.stderr.write(
                "ERROR\nInvalid AntennaS1528 L_s parameter: " + str(self.l_s),
            )
            sys.exit(1)

        self.b = 6.32
        self.alpha = 1.5

        self.x = self.peak_gain + self.l_s + \
            25 * math.log10(self.b * self.psi_b)
        self.y = self.b * self.psi_b * \
            math.pow(10, 0.04 * (self.peak_gain + self.l_s - self.l_f))

    def calculate_gain(self, *args, **kwargs) -> np.array:
        psi = np.absolute(kwargs["off_axis_angle_vec"])

        gain = np.zeros(len(psi))

        idx_0 = np.where(psi < self.a * self.psi_b)[0]
        gain[idx_0] = self.peak_gain - 3 * \
            np.power(psi[idx_0] / self.psi_b, self.alpha)

        idx_1 = np.where((self.a * self.psi_b < psi) &
                         (psi <= 0.5 * self.b * self.psi_b))[0]
        gain[idx_1] = self.peak_gain + self.l_s + 20 * math.log10(self.z)

        idx_2 = np.where((0.5 * self.b * self.psi_b < psi) &
                         (psi <= self.b * self.psi_b))[0]
        gain[idx_2] = self.peak_gain + self.l_s

        idx_3 = np.where((self.b * self.psi_b < psi) & (psi <= self.y))[0]
        gain[idx_3] = self.x - 25 * np.log10(psi[idx_3])

        idx_4 = np.where((self.y < psi) & (psi <= 90))[0]
        gain[idx_4] = self.l_f

        idx_5 = np.where((90 < psi) & (psi <= 180))[0]
        gain[idx_5] = self.l_b

        return gain


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    ## Plot gains for ITU-R-S.1528-SECTION1.2
    # initialize antenna parameters
    param = ParametersAntennaS1528()
    param.antenna_gain = 30
    param.antenna_pattern = "ITU-R-S.1528-SECTION1.2"
    param.antenna_3_dB = 4.4127

    psi = np.linspace(0, 30, num=1000)

    param.antenna_l_s = -15
    antenna = AntennaS1528(param)
    gain15 = antenna.calculate_gain(off_axis_angle_vec=psi)

    param.antenna_l_s = -20
    antenna = AntennaS1528(param)
    gain20 = antenna.calculate_gain(off_axis_angle_vec=psi)

    param.antenna_l_s = -25
    antenna = AntennaS1528(param)
    gain25 = antenna.calculate_gain(off_axis_angle_vec=psi)

    param.antenna_l_s = -30
    antenna = AntennaS1528(param)
    gain30 = antenna.calculate_gain(off_axis_angle_vec=psi)

    ## Plot gains for ITU-R-S.1528-LEO
    # initialize antenna parameters
    param = ParametersAntennaS1528()
    param.antenna_gain = 30
    param.antenna_pattern = "ITU-R-S.1528-LEO"
    param.antenna_3_dB = 1.6
    psi = np.linspace(0, 20, num=1000)

    param.antenna_l_s = -6.75
    antenna = AntennaS1528Leo(param)
    gain_leo = antenna.calculate_gain(off_axis_angle_vec=psi)

    fig = plt.figure(figsize=(8, 7), facecolor='w',
                     edgecolor='k')  # create a figure object

    psi_norm = psi / (param.antenna_3_dB / 2)
    plt.plot(psi_norm, gain15 - param.antenna_gain, "-b", label="$L_S = -15$ dB")
    plt.plot(psi_norm, gain20 - param.antenna_gain, "-r", label="$L_S = -20$ dB")
    plt.plot(psi_norm, gain25 - param.antenna_gain, "-g", label="$L_S = -25$ dB")
    plt.plot(psi_norm, gain30 - param.antenna_gain, "-k", label="$L_S = -30$ dB")
    plt.plot(psi_norm, gain_leo - param.antenna_gain, "-c", label="$L_S = -6.75$ dB (R1.3 - LEO)")

    plt.ylim((-40, 10))
    plt.xlim((0, np.max(psi_norm)))
    plt.xticks(np.arange(np.floor(np.max(psi_norm))))
    plt.title("ITU-R S.1528-0 antenna radiation pattern")
    plt.xlabel(r"Relative off-axis angle, $\psi/\psi_{3dB}$")
    plt.ylabel(r"Gain relative to $G_{max}$ [dB]")
    plt.legend(loc="upper right")
    plt.grid()

    ## Plot gains for ITU-R-S.1528-LEO
    # initialize antenna parameters
    param = ParametersAntennaS1528()
    param.antenna_gain = 35
    param.antenna_pattern = "ITU-R-S.1528-LEO"
    param.antenna_3_dB = 1.6
    psi = np.linspace(0, 20, num=1000)

    param.antenna_l_s = -6.75
    antenna = AntennaS1528Leo(param)
    gain_leo = antenna.calculate_gain(off_axis_angle_vec=psi)

    fig = plt.figure(figsize=(8, 7), facecolor='w', edgecolor='k')  # create a figure object
    psi_norm = psi / (param.antenna_3_dB / 2)
    plt.plot(psi_norm, gain_leo, "-b", label="$L_S = -6.75$ dB")

    # plt.ylim((-40, 10))
    plt.xlim((0, np.max(psi_norm)))
    plt.xticks(np.arange(np.floor(np.max(psi_norm))))
    plt.title("ITU-R S.1528-0 LEO antenna radiation pattern")
    plt.xlabel(r"Relative off-axis angle, $\psi/\psi_{3dB}$")
    plt.ylabel(r"Gain relative to $G_{max}$ [dB]")
    plt.legend(loc="upper right")
    plt.grid()

    # Section 1.4 (Taylor)
    params = ParametersAntennaS1528(
        antenna_gain=0,
        frequency=6000,
        bandwidth=500,
        slr=20,
        n_side_lobes=4,
        l_r=0.5,
        l_t=0.5,
        roll_off=3
    )

    # Create an instance of AntennaS1528Taylor
    antenna = AntennaS1528Taylor(params)
    print(f"Taylor antenna.lamb = {antenna.lamb}")

    # Define phi angles from 0 to 60 degrees for plotting
    theta_angles = np.linspace(0, 60, 600)

    # Calculate gains for each phi angle at a fixed theta angle (e.g., theta=0)
    gain = antenna.calculate_gain(theta=theta_angles, phi=np.zeros_like(theta_angles))

    # Plot the antenna gain as a function of phi angle
    plt.figure(figsize=(10, 6))
    plt.plot(theta_angles, gain)
    plt.xlabel('Theta (degrees)')
    plt.ylabel('Gain (dB)')
    plt.title('Normalized Antenna - Section 1.4')
    plt.grid(True)

    plt.show()
