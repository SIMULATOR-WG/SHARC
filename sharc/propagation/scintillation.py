# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 14:59:27 2017

@author: Andre Noll Barreto
"""

from sharc.propagation.propagation import Propagation
from sharc.propagation.atmosphere import ReferenceAtmosphere

import numpy as np

class Scintillation():
    """
    Implements the scintillation attenuation according to ITU-R P619
    (currently only tropospheric scintillation is implemented, ionospheric scintillation not available)
    """

    def __init__(self, random_number_gen: np.random.RandomState):
        self.random_number_gen = random_number_gen

        self.atmosphere = ReferenceAtmosphere()


    def get_tropospheric_attenuation (self, *args, **kwargs) -> np.array:
        """
        Calculates tropospheric scintillation based on ITU-R P.619, Appendix D

        Parameters
        ----------
            elevation (np.array) : free space elevation (degrees)
            frequency_MHz (float) : carrier frequency (MHz)
            antenna_gain_dB (np.array) : antenna gains at earth station (dBi)
            time_ratio (np.array / string) : percentage time that gain is exceeded
                                             if "random", then random values are chosen for each link
                                             default = "random"
            wet_refractivity (float) : wet term of the radio refractivity - optional
                                       if not given, then it is calculated from sat_params
            sat_params (ParametersFss) : satellite channel parameters - optional
                                         needed only if wet_refractivity is not given


        Returns
        -------
            attenuation (np.array): attenuation (dB) with dimensions equal to "elevation"
        """

        f_GHz = kwargs["frequency_MHz"] / 1000.
        elevation_rad = kwargs["elevation"] / 180. * np.pi
        antenna_gain_dB = kwargs["antenna_gain_dB"]
        if not np.isscalar(antenna_gain_dB):
            antenna_gain_dB = antenna_gain_dB.flatten()
        time_ratio = kwargs.pop("time_ratio", "random")
        wet_refractivity = kwargs.pop("wet_refractivity", False)

        if not wet_refractivity:
            sat_params = kwargs["sat_params"]

            temperature, \
            pressure, \
            water_vapour_density = self.atmosphere.get_reference_atmosphere_p835(sat_params.imt_lat_deg,
                                                                                 sat_params.imt_altitude,
                                                                                 sat_params.season)

            # calculate saturation water vapour pressure according to ITU-R P.453-12
            # water coefficients (ice disregarded)
            coef_a = 6.1121
            coef_b = 18.678
            coef_c = 257.14
            coef_d = 234.5
            EF_water = 1 + 1e-4 * (7.2 + pressure * (.032 + 5.9e-6 * temperature ** 2))
            vapour_pressure = water_vapour_density * temperature / 216.7  # eq 10 in P453
            wet_refractivity = (72 * vapour_pressure / temperature
                                + 3.75e5 * vapour_pressure / temperature ** 2)

        sigma_ref = 3.6e-3 + 1e-4 * wet_refractivity
        h_l = 1000
        path_length = 2 * h_l / (np.sqrt(np.sin(elevation_rad) ** 2 + 2.35e-4) + np.sin(elevation_rad))
        eff_antenna_diameter = .3 * 10 ** (.05 * antenna_gain_dB) / (np.pi * f_GHz)

        x = 1.22 * eff_antenna_diameter ** 2 * (f_GHz / path_length)
        antenna_averaging_factor = np.sqrt(3.86 * (x ** 2 + 1) ** (11 / 12) *
                                           np.sin(11 / 6 * np.arctan(1 / x)) - 7.08 * x ** 5 / 6)
        scintillation_intensity = (sigma_ref * f_GHz ** (7 / 12) * antenna_averaging_factor
                                   / np.sin(elevation_rad) ** 1.2)

        if isinstance(time_ratio, str) and time_ratio.lower() == "random":
            time_ratio = self.random_number_gen.rand(len(elevation_rad))

        # tropospheric scintillation attenuation not exceeded for time_percentage percent time
        time_percentage = time_ratio * 100.

        num_el = 1
        if np.isscalar(scintillation_intensity):
            if not np.isscalar(time_percentage):
                num_el = time_percentage.size
                scintillation_intensity = np.ones(num_el) * scintillation_intensity
        else:
            num_el = scintillation_intensity.size
            if np.isscalar(time_percentage):
                time_percentage = np.ones(num_el) * time_percentage

        attenuation = np.empty(num_el)

        a_ste = (2.672 - 1.258 * np.log10(time_percentage) -
                 .0835 * np.log10(time_percentage) ** 2 -
                 .0597 * np.log10(time_percentage) ** 3)

        a_stf = (3. - 1.71 * np.log10(100 - time_percentage) +
                 .072 * np.log10(100 - time_percentage) ** 2 -
                 .061 * np.log10(100 - time_percentage) ** 3)

        gain_indices = np.where(time_percentage <= 50.)[0]
        if gain_indices.size:
            attenuation[gain_indices] = - scintillation_intensity[gain_indices] * a_ste[gain_indices]
        fade_indices = np.where(time_percentage > 50.)[0]
        if fade_indices.size:
            attenuation[fade_indices] = scintillation_intensity[fade_indices] * a_stf[fade_indices]

        return attenuation

if __name__ == '__main__':
    from matplotlib import pyplot as plt

    propagation = Scintillation(np.random.RandomState(101))

    #################################
    # Plot troposcatter scintillation attenuation
    # compare with benchmark from ITU-R P-619 Fig. 8
    percentage_fading_exceeded = 10 ** np.arange(-2, 1.1, .1)
    antenna_gain = 33.5
    frequency_MHz = 2700.
    wet_refractivity = 42.5

    elevation_vec = np.array([5., 10., 20., 35., 90.])

    print("Plotting troposcatter scintillation attenuation:")

    plt.figure()
    for elevation in elevation_vec:
        attenuation = propagation.get_tropospheric_attenuation(elevation=elevation,
                                                               frequency_MHz=frequency_MHz,
                                                               antenna_gain_dB=antenna_gain,
                                                               time_ratio=1 - (percentage_fading_exceeded / 100),
                                                               wet_refractivity=wet_refractivity)
        plt.semilogx(percentage_fading_exceeded, attenuation, label="{} deg".format(elevation))

    percentage_gain_exceeded = 10 ** np.arange(-2, 1.1, .1)
    for elevation in elevation_vec:
        attenuation = propagation.get_tropospheric_attenuation(elevation=elevation,
                                                               frequency_MHz=frequency_MHz,
                                                               antenna_gain_dB=antenna_gain,
                                                               time_ratio=percentage_gain_exceeded / 100,
                                                               wet_refractivity=wet_refractivity)
        plt.loglog(percentage_gain_exceeded, np.abs(attenuation), ':', label="{} deg".format(elevation))

    plt.legend(title='elevation')
    plt.grid(True)

    plt.title("Troposcatter Scintillation Attenuation")
    plt.xlabel("Percentage time fades and enhancements exceeded")
    plt.ylabel("Enhancement/Fade (dB)")

    plt.show()
