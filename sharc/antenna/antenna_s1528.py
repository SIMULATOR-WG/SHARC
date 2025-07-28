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

NUM_OF_BESSEL_ROOTS = 3


class AntennaS1528Taylor(Antenna):
    """
    Implements Recommendation ITU-R S.1528-0 Section 1.4: Satellite antenna reference pattern given by an analytical
    function which models the side lobes of the non-GSO satellite operating in the fixed-satellite service below 30 GHz.
    It is proposed to use a circular Taylor illumination function which gives the maximum flexibility to
    adapt the theoretical pattern to the real one. It takes into account the side-lobes effect of an antenna
    diagram.
    """

    def __init__(self, param: ParametersAntennaS1528):
        super().__init__()
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

        # Radial (l_r) and transverse (l_t) sizes of the effective radiating area of the satellite transmitt antenna (m)
        self.l_r = param.l_r
        self.l_t = param.l_t

        # Intermediary variables
        self.A = (1 / np.pi) * np.arccosh(10 ** (self.slr / 20))
        self.j1_roots = jn_zeros(1, self.n_side_lobes) / np.pi
        self.sigma = self.j1_roots[-1] / np.sqrt(self.A ** 2 + (self.n_side_lobes - 1 / 2) ** 2)
        self.mu = jn_zeros(1, NUM_OF_BESSEL_ROOTS) / np.pi

    def calculate_gain(self, *args, **kwargs) -> np.array:
        # The reference angles for the simulator and the antenna realisation are switched.
        # Local theta is simulator off_axis_angle and local phi is simulator theta_vec
        if 'off_axis_angle_vec' not in kwargs:
            raise ValueError("off_axis_angle_vec vector must be given")
        if 'theta_vec' not in kwargs:
            raise ValueError("theta vector must be given")
        theta = np.abs(np.radians(kwargs.get('off_axis_angle_vec', 0)))
        phi = np.abs(np.radians(kwargs.get('theta_vec', 0)))

        u = (np.pi / self.lamb) * np.sqrt((self.l_r * np.sin(theta) * np.cos(phi)) ** 2 +
                                          (self.l_t * np.sin(theta) * np.sin(phi)) ** 2)

        v = np.ones(u.shape + (NUM_OF_BESSEL_ROOTS,))

        for i, ui in enumerate(self.mu):
            v[..., i] = (1 - u ** 2 / (np.pi ** 2 * self.sigma ** 2 *
                                       (self.A ** 2 + (i + 1 - 0.5) ** 2))) / (1 - (u / (np.pi * ui)) ** 2)

        # Take care of divide-by-zero
        with np.errstate(divide='ignore', invalid='ignore'):
            gain = self.peak_gain + 20 * \
                np.log10(np.abs((2 * jn(1, u) / u) * np.prod(v, axis=-1)))

        # Replace undefined values with -inf (or other desired value)
        gain = np.nan_to_num(gain, nan=-np.inf)
        # Replace values that were substituded specifically
        # because u == 0 with Lim (gain)_(u -> 0) = peak_gain
        gain[u == 0] = self.peak_gain

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
        self.psi_b = param.antenna_3_dB_bw / 2
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
        self.psi_b = param.antenna_3_dB_bw / 2

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
    param.antenna_gain = 37
    param.antenna_pattern = "ITU-R-S.1528-SECTION1.2"
    param.antenna_3_dB_bw = 4.4127

    max_psi = 30 * param.antenna_3_dB_bw / 2
    psi = np.linspace(0, max_psi, num=100)

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
    param.antenna_pattern = "ITU-R-S.1528-LEO"

    param.antenna_l_s = -6.75
    antenna = AntennaS1528Leo(param)
    gain_leo = antenna.calculate_gain(off_axis_angle_vec=psi)

    fig = plt.figure(figsize=(8, 7), facecolor='w',
                     edgecolor='k')  # create a figure object

    psi_norm = psi / (param.antenna_3_dB_bw / 2)
    plt.plot(psi_norm, gain15 - param.antenna_gain, "-b", label="$L_S = -15$ dB")
    plt.plot(psi_norm, gain20 - param.antenna_gain, "-r", label="$L_S = -20$ dB")
    plt.plot(psi_norm, gain25 - param.antenna_gain, "-g", label="$L_S = -25$ dB")
    plt.plot(psi_norm, gain30 - param.antenna_gain, "-k", label="$L_S = -30$ dB")
    plt.plot(psi_norm, gain_leo - param.antenna_gain, "-c", label="$L_S = -6.75$ dB (R1.3 - LEO)")

    plt.ylim((-40, 10))
    plt.yticks(np.arange(-40, 11, 5))
    plt.xlim((0, np.max(psi_norm)))
    plt.xticks(np.arange(0, 35, 5))
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
    param.antenna_3_dB_bw = 1.6
    psi = np.linspace(0, 20, num=1000)

    param.antenna_l_s = -6.75
    antenna = AntennaS1528Leo(param)
    gain_leo = antenna.calculate_gain(off_axis_angle_vec=psi)

    fig = plt.figure(figsize=(8, 7), facecolor='w', edgecolor='k')  # create a figure object
    psi_norm = psi / (param.antenna_3_dB_bw / 2)
    plt.plot(psi_norm, gain_leo, "-b", label="$L_S = -6.75$ dB")

    # plt.ylim((-40, 10))
    plt.xlim((0, np.max(psi_norm)))
    plt.xticks(np.arange(np.floor(np.max(psi_norm))))
    plt.title("ITU-R S.1528-0 LEO antenna radiation pattern")
    plt.xlabel(r"Relative off-axis angle, $\psi/\psi_{3dB}$")
    plt.ylabel(r"Gain relative to $G_{max}$ [dB]")
    plt.legend(loc="upper right")
    plt.grid()

    # Section 1.4 (Taylor) - Compare to Fig 6
    frequency = 12000  # MHz
    bandwidth = 10  # MHz
    antenna_gain = 0  # dBi
    slr = 20  # dB
    n_side_lobes = 4
    lamb = (SPEED_OF_LIGHT / 1e6) / (frequency - bandwidth / 2)
    beam_radius = 350  # km
    sat_altitude = 1446  # km
    a = np.arctan(beam_radius / (sat_altitude))  # radians
    l_r = 0.74 * lamb / np.sin(a)
    l_t = l_r
    params_rolloff_7 = ParametersAntennaS1528(
        antenna_gain=0,
        frequency=12000,
        bandwidth=10,
        slr=20,
        n_side_lobes=4,
        l_r=l_r,
        l_t=l_t,
    )

    # Create an instance of AntennaS1528Taylor
    antenna_rolloff_7 = AntennaS1528Taylor(params_rolloff_7)

    # Define phi angles from 0 to 60 degrees for plotting
    theta_angles = np.arange(0, 60.1, 0.1)

    # Calculate gains for each phi angle at a fixed theta angle (e.g., theta=0)
    gain_rolloff_7 = antenna_rolloff_7.calculate_gain(off_axis_angle_vec=theta_angles,
                                                      theta_vec=np.zeros_like(theta_angles))

    l_r = 0.64 * lamb / np.sin(a)
    l_t = l_r
    params_rolloff_5 = ParametersAntennaS1528(
        antenna_gain=0,
        frequency=12000,
        bandwidth=10,
        slr=20,
        n_side_lobes=4,
        l_r=l_r,
        l_t=l_t,
    )

    # Create an instance of AntennaS1528Taylor
    antenna_rolloff_5 = AntennaS1528Taylor(params_rolloff_5)

    gain_rolloff_5 = antenna_rolloff_5.calculate_gain(off_axis_angle_vec=theta_angles,
                                                      theta_vec=np.zeros_like(theta_angles))

    l_r = 0.51 * lamb / np.sin(a)
    l_t = l_r
    params_rolloff_3 = ParametersAntennaS1528(
        antenna_gain=0,
        frequency=12000,
        bandwidth=10,
        slr=20,
        n_side_lobes=4,
        l_r=l_r,
        l_t=l_t,
    )

    # Create an instance of AntennaS1528Taylor
    antenna_rolloff_3 = AntennaS1528Taylor(params_rolloff_3)

    gain_rolloff_3 = antenna_rolloff_3.calculate_gain(off_axis_angle_vec=theta_angles,
                                                      theta_vec=np.zeros_like(theta_angles))

    # Plot the antenna gain as a function of phi angle
    plt.figure(figsize=(10, 6))
    plt.plot(theta_angles, gain_rolloff_3, label='roll_off=3')
    plt.plot(theta_angles, gain_rolloff_5, label='roll_off=5')
    plt.plot(theta_angles, gain_rolloff_7, label='roll_off=7')
    plt.xlabel('Theta (degrees)')
    plt.ylabel('Gain (dB)')
    plt.title('Normalized Antenna - Section 1.4')
    plt.legend()
    plt.xticks(np.arange(0, 60, 10))
    plt.minorticks_on()
    plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(2))
    plt.grid(True, which='both')

    plt.show()
