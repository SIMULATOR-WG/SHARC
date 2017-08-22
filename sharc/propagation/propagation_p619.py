# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 15:35:00 2017

@author: andre barreto
"""

from sharc.propagation.propagation import Propagation

import math
import numpy as np
from sharc.propagation.propagation_free_space import PropagationFreeSpace
from sharc.propagation.propagation_clutter_loss import PropagationClutterLoss
from sharc.support.enumerations import StationType
from scipy.stats import norm
import matplotlib.pyplot as plt
import os

class PropagationP619(Propagation):
    """
    Implements the earth-to-space channel model from ITU-R P.619

    Public methods:
        get_loss: Calculates path loss for earth-space link
    """

    def __init__(self):
        super().__init__()
        self.clutter = PropagationClutterLoss()
        self.free_space = PropagationFreeSpace()
        self.depolarization_loss = 1.5
        self.polarization_mismatch_loss = 3.
        self.building_loss = 20

    @staticmethod
    def _get_reference_atmosphere_p835 (latitude, altitude = 1000, season = "summer"):
        """
        Returns reference atmosphere parameters based on ITU-R P835-5

        Parameters
         ----------
             latitude (float): latitude (degrees)
             altitude (float): altitude (m)
             season (string): season of the year, "summer"/"winter"

         Returns
         -------
             temperature (float): temperature (K)
             pressure_hPa (float): air pressure (hPa)
             water_vapour_density (float): density of water vapour (g/m**3)
        """

        h_km = altitude / 1000

        if latitude <= 22:
            # low latitude
            if h_km < 17.:
                temperature = 300.4222 - 6.3533 * h_km + .005886 * h_km **2
            elif h_km < 47:
                temperature = 194 + (h_km - 17) * 2.533
            elif h_km < 52:
                temperature = 270.
            elif h_km < 80:
                temperature = 270. - (h_km - 52) * 3.0714
            elif h_km <= 100:
                temperature = 184.
            else:
                error_message = "altitude > 100km not supported"
                raise ValueError(error_message)

            if h_km < 10.:
                pressure_hPa = 1012.0306 - 109.0338 * h_km + 3.6316 * h_km ** 2
            elif h_km <= 72.:
                p10 = 1012.0306 - 109.0338 * 10 + 3.6316 * 10 ** 2
                pressure_hPa = p10 * np.exp(-.147 * (h_km - 10.))
            elif h_km <= 100:
                p72 = ( 1012.0306 - 109.0338 * 10 + 3.6316 * 10 ** 2 ) * np.exp(-.147 * (72 - 10))
                pressure_hPa = p72 * np.exp(-.165 * (h - 72))
            else:
                error_message = "altitude > 100km not supported"
                raise ValueError(error_message)

            if h_km <= 15.:
                water_vapour_density = 19.6542 * np.exp( -.2313 * h_km - .1122 * h_km ** 2
                                                          + .01351 * h_km **3 - .0005923 * h_km ** 4)
            else:
                water_vapour_density = 0

        elif latitude <= 45.:
            # mid-latitude
            if season == "summer":
                if h_km < 13.:
                    temperature = 294.9838 - 5.2159 * h_km - .07109 * h_km ** 2
                elif h_km < 17.:
                    temperature = 215.15
                elif h_km < 47:
                    temperature = 215.15 * np.exp((h_km - 17.) * .008128)
                elif h_km < 53.:
                    temperature = 275.
                elif h_km < 80.:
                    temperature = 275 + (1 - np.exp((h - 53.)*.06)) * 20.
                elif h_km <= 100:
                    temperature = 175.
                else:
                    error_message = "altitude > 100km not supported"
                    raise ValueError(error_message)

                if h_km <= 10:
                    pressure_hPa = 1012.8186 - 111.5569 * h_km + 3.8646 * h_km ** 2
                elif h_km <=72.:
                    p10 = 1012.8186 - 111.5569 * 10 + 3.8646 * 10 ** 2
                    pressure_hPa = p10 * np.exp(-.147 * (h_km-10))
                elif h_km <= 100.:
                    p72 = (1012.8186 - 111.5569 * 10 + 3.8646 * 10 ** 2) * np.exp(-.147 * (72-10))
                    pressure_hPa = p72 * np.exp(-.165 * (h_km - 72))
                else:
                    error_message = "altitude > 100km not supported"
                    raise ValueError(error_message)

                if h_km <= 15.:
                    water_vapour_density = 14.3542 * np.exp(-.4174 * h_km - .02290 * h_km ** 2
                                                             + .001007 * h_km ** 3)
                else:
                    water_vapour_density = 0

            elif season == "winter":
                if h_km < 13.:
                    temperature = 272.7241 - 3.6217 * h_km - .1759 * h_km ** 2
                elif h_km < 33.:
                    temperature = 218.
                elif h_km < 47:
                    temperature = 218. * (h_km - 33.) * 3.3571
                elif h_km < 53.:
                    temperature = 265.
                elif h_km < 80.:
                    temperature = 265 + (h_km - 53.) * 2.037
                elif h_km <= 100:
                    temperature = 210.
                else:
                    error_message = "altitude > 100km not supported"
                    raise ValueError(error_message)

                if h_km <= 10:
                    pressure_hPa = 1018.8627 - 124.2954 * h_km + 4.8307 * h_km ** 2
                elif h_km <=72.:
                    p10 = 1018.8627 - 124.2954 * 10 + 4.8307 * 10 ** 2
                    pressure_hPa = p10 * np.exp(-.147 * (h_km-10))
                elif h_km <= 100.:
                    p72 = (1018.8627 - 124.2954 * 10 + 4.8307 * 10 ** 2) * np.exp(-.147 * (72-10))
                    pressure_hPa = p72 * np.exp(-.155 * (h_km - 72))
                else:
                    error_message = "altitude > 100km not supported"
                    raise ValueError(error_message)

                if h_km <= 15.:
                    water_vapour_density = 3.4742 * np.exp(-.2697 * h_km - .03604 * h_km ** 2
                                                             + .0004489 * h_km ** 3)
                else:
                    water_vapour_density = 0
            else:
                error_message = ("season {} not supported".format(season))
                raise ValueError(error_message)
        else:
            # high latitude (>45 deg)
            if season == "summer":
                if h_km < 13.:
                    temperature = 286.8374 - 4.7805 * h_km - 0.1402 * h_km ** 2
                elif h_km < 23.:
                    temperature = 225.
                elif h_km < 48:
                    temperature = 225. * np.exp((h_km - 23.) * .008317)
                elif h_km < 53.:
                    temperature = 277.
                elif h_km < 79.:
                    temperature = 277 - (h_km - 53.) * 4.0769
                elif h_km <= 100:
                    temperature = 171.
                else:
                    error_message = "altitude > 100km not supported"
                    raise ValueError(error_message)

                if h_km <= 10:
                    pressure_hPa = 1008.0278 - 113.2494 * h_km + 3.9408 * h_km ** 2
                elif h_km <= 72.:
                    p10 = 1008.0278 - 113.2494 * 10 + 3.9408 * (10 ** 2)
                    pressure_hPa = p10 * np.exp(-.140 * (h_km - 10))
                elif h_km <= 100.:
                    p72 = (1008.0278 - 113.2494 * 10 + 3.9408 * (10 ** 2)) * np.exp(-.140 * (72 - 10))
                    pressure_hPa = p72 * np.exp(-.165 * (h_km - 72))
                else:
                    error_message = "altitude > 100km not supported"
                    raise ValueError(error_message)

                if h_km <= 15.:
                    water_vapour_density = 8.988 * np.exp(-.3614 * h_km - .005402 * h_km ** 2
                                                             + .001955 * h_km ** 3)
                else:
                    water_vapour_density = 0

            elif season == "winter":
                if h_km < 8.5:
                    temperature = (257.4345 + 2.3474 * h_km - .15479 * h_km ** 2
                                   + .08473 * h_km ** 3)
                elif h_km < 30.:
                    temperature = 217.5
                elif h_km < 50:
                    temperature = 217.5 * (h_km - 30.) * 2.125
                elif h_km < 54.:
                    temperature = 260.
                elif h_km <= 100:
                    temperature = 260. - (h_km - 54.) * 1.667
                else:
                    error_message = "altitude > 100km not supported"
                    raise ValueError(error_message)

                if h_km <= 10:
                    pressure_hPa = 1010.8828 - 122.2411 * h_km + 4.554 * h_km ** 2
                elif h_km <= 72.:
                    p10 = 1010.8828 - 122.2411 * 10 + 4.554 * 10 ** 2
                    pressure_hPa = p10 * np.exp(-.147 * (h_km - 10))
                elif h_km <= 100.:
                    p72 = (1010.8828 - 122.2411 * 10 + 4.554 * 10 ** 2) * np.exp(-.147 * (72 - 10))
                    pressure_hPa = p72 * np.exp(-.150 * (h_km - 72))
                else:
                    error_message = "altitude > 100km not supported"
                    raise ValueError(error_message)

                if h_km <= 15.:
                    water_vapour_density = 1.2319 * np.exp(.07481 * h_km - .0981 * h_km ** 2
                                                            + .00281 * h_km ** 3)
                else:
                    water_vapour_density = 0
            else:
                error_message = ("season {} not supported".format(season))
                raise ValueError(error_message)

        return temperature, pressure_hPa, water_vapour_density

    @staticmethod
    def _get_specific_attenuation(pressure, water_vapour_pressure, temperature, frequency_MHz):
        """
         Calculates specific attenuation (dB/km) of an atmosphere layer according to ITU-T P.676-11

         Parameters
         ----------
             pressure (float): dry-air partial pressure (hPa)
             water_vapour_pressure (float): water-vapour partial pressure (hPa)
             temperature (float): temperature(K)
             frequency_MHz (float) : carrier frequency (MHz)

         Returns
         -------
             specific_attenuation (float): specific gaseous attenuation (dB/km)
         """

        oxygen_f0 = np.array(
                    [50.474214, 50.987745, 51.503360, 52.021429, 52.542418, 53.066934, 53.595775,
                     54.130025, 54.671180, 55.221384, 55.783815, 56.264774, 56.363399, 56.968211,
                     57.612486, 58.323877, 58.446588, 59.164204, 59.590983, 60.306056, 60.434778,
                     61.150562, 61.800158, 62.411220, 62.486253, 62.997984, 63.568526, 64.127775,
                     64.678910, 65.224078, 65.764779, 66.302096, 66.836834, 67.369601, 67.900868,
                     68.431006, 68.960312, 118.750334, 368.498246, 424.763020, 487.249273, 715.392902,
                     773.839490, 834.145546])
        a = np.array(
            [[0.975, 9.651, 6.690, 0.0, 2.566, 6.850], [2.529, 8.653, 7.170, 0.0, 2.246, 6.800],
             [6.193, 7.709, 7.640, 0.0, 1.947, 6.729], [14.320, 6.819, 8.110, 0.0, 1.667, 6.640],
             [31.240, 5.983, 8.580, 0.0, 1.388, 6.526], [64.290, 5.201, 9.060, 0.0, 1.349, 6.206],
             [124.600, 4.474, 9.550, 0.0, 2.227, 5.085], [227.300, 3.800, 9.960, 0.0, 3.170, 3.750],
             [389.700, 3.182, 10.370, 0.0, 3.558, 2.654], [627.100, 2.618, 10.890, 0.0, 2.560, 2.952],
             [945.300, 2.109, 11.340, 0.0, -1.172, 6.135], [543.400, 0.014, 17.030, 0.0, 3.525, -0.978],
             [1331.800, 1.654, 11.890, 0.0, -2.378, 6.547], [1746.600, 1.255, 12.230, 0.0, -3.545, 6.451],
             [2120.100, 0.910, 12.620, 0.0, -5.416, 6.056], [2363.700, 0.621, 12.950, 0.0, -1.932, 0.436],
             [1442.100, 0.083, 14.910, 0.0, 6.768, -1.273], [2379.900, 0.387, 13.530, 0.0, -6.561, 2.309],
             [2090.700, 0.207, 14.080, 0.0, 6.957, -0.776], [2103.400, 0.207, 14.150, 0.0, -6.395, 0.699],
             [2438.000, 0.386, 13.390, 0.0, 6.342, -2.825], [2479.500, 0.621, 12.920, 0.0, 1.014, -0.584],
             [2275.900, 0.910, 12.630, 0.0, 5.014, -6.619], [1915.400, 1.255, 12.170, 0.0, 3.029, -6.759],
             [1503.000, 0.083, 15.130, 0.0, -4.499, 0.844], [1490.200, 1.654, 11.740, 0.0, 1.856, -6.675],
             [1078.000, 2.108, 11.340, 0.0, 0.658, -6.139], [728.700, 2.617, 10.880, 0.0, -3.036, -2.895],
             [461.300, 3.181, 10.380, 0.0, -3.968, -2.590], [274.000, 3.800, 9.960, 0.0, -3.528, -3.680],
             [153.000, 4.473, 9.550, 0.0, -2.548, -5.002], [80.400, 5.200, 9.060, 0.0, -1.660, -6.091],
             [39.800, 5.982, 8.580, 0.0, -1.680, -6.393], [18.560, 6.818, 8.110, 0.0, -1.956, -6.475],
             [8.172, 7.708, 7.640, 0.0, -2.216, -6.545], [3.397, 8.652, 7.170, 0.0, -2.492, -6.600],
             [1.334, 9.650, 6.690, 0.0, -2.773, -6.650], [940.300, 0.010, 16.640, 0.0, -0.439, 0.079],
             [67.400, 0.048, 16.400, 0.0, 0.000, 0.000], [637.700, 0.044, 16.400, 0.0, 0.000, 0.000],
             [237.400, 0.049, 16.000, 0.0, 0.000, 0.000], [98.100, 0.145, 16.000, 0.0, 0.000, 0.000],
             [572.300, 0.141, 16.200, 0.0, 0.000, 0.000], [183.100, 0.145, 14.700, 0.0, 0.000, 0.000]])

        vapour_f0 = np.array([22.235080, 67.803960, 119.995940, 183.310087, 321.225630, 325.152888,
                              336.227764, 380.197353, 390.134508, 437.346667, 439.150807, 443.018343,
                              448.001085, 470.888999, 474.689092, 488.490108, 503.568532, 504.482692,
                              547.676440, 552.020960, 556.935985, 620.700807, 645.766085, 658.005280,
                              752.033113, 841.051732, 859.965698, 899.303175, 902.611085, 906.205957,
                              916.171582, 923.112692, 970.315022, 987.926764, 1780.000000])
        vapour_indices = np.array([0, 3, 4, 5, 7, 12, 20, 24, 34])

        b = np.array(
            [[.1079, 2.144, 26.38, .76, 5.087, 1.00], [.0011, 8.732, 28.58, .69, 4.930, .82],
             [.0007, 8.353, 29.48, .70, 4.780, .79], [2.273, .668, 29.06, .77, 5.022, .85],
             [.0470, 6.179, 24.04, .67, 4.398, .54], [1.514, 1.541, 28.23, .64, 4.893, .74],
             [.0010, 9.825, 26.93, .69, 4.740, .61], [11.67, 1.048, 28.11, .54, 5.063, .89],
             [.0045, 7.347, 21.52, .63, 4.810, .55], [.0632, 5.048, 18.45, .60, 4.230, .48],
             [.9098, 3.595, 20.07, .63, 4.483, .52], [.1920, 5.048, 15.55, .60, 5.083, .50],
             [10.41, 1.405, 25.64, .66, 5.028, .67], [.3254, 3.597, 21.34, .66, 4.506, .65],
             [1.260, 2.379, 23.20, .65, 4.804, .64], [.2529, 2.852, 25.86, .69, 5.201, .72],
             [.0372, 6.731, 16.12, .61, 3.980, .43], [.0124, 6.731, 16.12, .61, 4.010, .45],
             [.9785, .158, 26.00, .70, 4.500, 1.00], [.1840, .158, 26.00, .70, 4.500, 1.00],
             [497.0, .159, 30.86, .69, 4.552, 1.00], [5.015, 2.391, 24.38, .71, 4.856, .68],
             [.0067, 8.633, 18.00, .60, 4.000, .50], [.2732, 7.816, 32.10, .69, 4.140, 1.00],
             [243.4, .396, 30.86, .68, 4.352, .84], [.0134, 8.177, 15.90, .33, 5.760, .45],
             [.1325, 8.055, 30.60, .68, 4.090, .84], [.0547, 7.914, 29.85, .68, 4.530, .90],
             [.0386, 8.429, 28.65, .70, 5.100, .95], [.1836, 5.110, 24.08, .70, 4.700, .53],
             [8.400, 1.441, 26.73, .70, 5.150, .78], [.0079, 10.293, 29.00, .70, 5.000, .80],
             [9.009, 1.919, 25.50, .64, 4.940, .67], [134.6, .257, 29.85, .68, 4.550, .90],
             [17506., .952, 196.3, 2.00, 24.15, 5.00]])

        theta = 300 / temperature
        f_GHz = frequency_MHz / 1000

        # line strength
        s_oxygen = a[:, 0] * 1e-7 * pressure * theta ** 3 * np.exp(a[:, 1] * (1 - theta))
        s_vapour = b[:, 0] * 1e-1 * water_vapour_pressure * theta ** 3.5 \
                   * np.exp(b[:, 1] * (1 - theta))
        # line width
        df_oxygen = a[:, 2] * 1e-4 * (pressure * theta ** (.8 - a[:, 3])
                                      + 1.1 * water_vapour_pressure * theta)
        df_vapour = b[:, 2] * 1e-4 * (pressure * theta ** (b[:, 3])
                                      + b[:, 4] * water_vapour_pressure * theta ** b[:, 5])
        # correction factor
        delta_oxygen = (a[:, 4] + a[:, 5] * theta) * 1e-4 * \
                       (pressure + water_vapour_pressure) * theta ** .8

        # line shape factor
        sf_oxygen = f_GHz / oxygen_f0 * \
                    ((df_oxygen - delta_oxygen * (oxygen_f0 - f_GHz)) / ((oxygen_f0 - f_GHz) ** 2 + df_oxygen ** 2) +
                     (df_oxygen - delta_oxygen * (oxygen_f0 + f_GHz)) / ((oxygen_f0 + f_GHz) ** 2 + df_oxygen ** 2))

        sf_vapour = f_GHz / vapour_f0 * \
                    (df_vapour / ((vapour_f0 - f_GHz) ** 2 + df_vapour ** 2) +
                     df_vapour / ((vapour_f0 + f_GHz) ** 2 + df_vapour ** 2))

        # Debye spectrum width
        dw = 5.6e-4 * (pressure + water_vapour_pressure) * theta ** .8
        # dry continuum due to pressure-induced nitrogen absorption and the Debye spectrum
        nd = f_GHz * pressure * theta ** 2 * \
             (
             6.14e-5 / (dw * (1 + (f_GHz / dw) ** 2)) + 1.4e-12 * pressure * theta ** 1.5 / (1 + 1.9e-5 * f_GHz ** 1.5))
        att_oxygen = 0.182 * f_GHz * (np.sum(s_oxygen * sf_oxygen) + nd)
        att_vapour = 0.182 * f_GHz * np.sum(s_vapour * sf_vapour)

        specific_attenuation = att_oxygen + att_vapour

        return specific_attenuation

    @staticmethod
    def _get_atmospheric_params(altitude_km, water_vapour_density_sea_level, f_MHz):
        """
         Calculates atmospheric parameters based on ITU-R P.619, Atttachment C.5

         Parameters
         ----------
             altitude_km (float) : altitude [km]
             water_vapour_density_sea_level (float) : water vapour density at sea level (g/m**3)
             f_MHz (float) : carrier frequency (MHz)

         Returns
         -------
             temperature (float): temperature(K) at altitude_km in reference atmosphere
             pressure (float): dry-air partial pressure at altitude_km in reference atmosphere (hPa)
             water_vapour_pressure (float): water-vapour partial pressure at altitude_km in reference atmosphere (hPa)
             refractive_index (float) : index of refraction of atmospheric layer at altitude_km
             specific_attenuation (float): specific gaseous attenuation (dB/km)
         """

        # Table C.1 of ITU-R P619 - Constants for the reference dry atmosphere
        ref_atmosphere_altitude_km = [-math.inf, 11, 20, 32, 47, 51, 71]
        ref_atmosphere_temp_grad = [-6.5, 0, 1., 2.8, 0, -2.8, -2]
        ref_atmosphere_temperature = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65]
        ref_atmosphere_pressure = [1013.25, 226.323, 54.750, 8.68, 1.109, 0.669, 0.04]

        index = np.where(altitude_km >= np.array(ref_atmosphere_altitude_km))[0][-1]

        Ti = ref_atmosphere_temperature[index]
        Li = ref_atmosphere_temp_grad[index]
        Hi = ref_atmosphere_altitude_km[index]
        if Hi < 0:
            Hi = 0
        Pi = ref_atmosphere_pressure[index]

        temperature = Ti + Li * (altitude_km - Hi)
        if Li:
            pressure = Pi * (Ti / (Ti + Li * (altitude_km - Hi))) ** (34.163 / Li)
        else:
            pressure = Pi * np.exp(-34.163 * (altitude_km - Hi) / Ti)

        water_vapour_density = water_vapour_density_sea_level * np.exp(-altitude_km / 2)
        water_vapour_pressure = water_vapour_density * temperature / 216.7
        refractive_index = 1 + 1e-6 * (77.6 / temperature *
                                       (pressure + water_vapour_pressure +
                                        4810 * water_vapour_pressure / pressure))
        specific_attenuation = PropagationP619._get_specific_attenuation(pressure, water_vapour_pressure,
                                                                         temperature, f_MHz)

        return [temperature, pressure, water_vapour_pressure, refractive_index, specific_attenuation]

    @staticmethod
    def _get_atmospheric_gasses_loss(*args, **kwargs) -> float:
        """
        Calculates atmospheric gasses loss based on ITU-R P.619, Attachment C

        Parameters
        ----------
            frequency_MHz (float) : center frequencies [MHz]
            apparent_elevation (float) : apparent elevation angle (degrees)
            sat_params (ParametersFss) : parameters of satellite system

        Returns
        -------
            path_loss (float): scalar with atmospheric loss
        """

        frequency_MHz = kwargs["frequency_MHz"]
        apparent_elevation = kwargs["apparent_elevation"]
        sat_params = kwargs["sat_params"]

        surf_water_vapour_density = kwargs.pop("surf_water_vapour_density", False)

        earth_radius_km = sat_params.EARTH_RADIUS/1000
        a_acc = 0. # accumulated attenuation (in dB)
        h = sat_params.imt_altitude/1000 # ray altitude in km
        beta = (90-abs(apparent_elevation)) * np.pi / 180. # incidence angle

        if not surf_water_vapour_density:
            dummy, dummy, surf_water_vapour_density = \
                PropagationP619._get_reference_atmosphere_p835(sat_params.imt_lat_deg,
                                                               0, season="summer")

        rho_s = surf_water_vapour_density * np.exp(h/2) # water vapour density at h
        if apparent_elevation < 0:
            # get temperature (t), dry-air pressure (p), water-vapour pressure (e),
            #     refractive index (n) and specific attenuation (gamma)
            t, p, e, n, gamma = PropagationP619._get_atmospheric_params(h, rho_s, frequency_MHz)
            delta = .0001 + 0.01 * max(h, 0) # layer thickness
            r = earth_radius_km + h - delta # radius of lower edge
            while True:
                m = (r + delta) * np.sin(beta) - r
                if m >= 0:
                    dh = 2 * np.sqrt(2*r*(delta-m)+delta**2-m**2) # horizontal path
                    a_acc += dh * gamma
                    break
                ds = (r+delta)*np.cos(beta)-np.sqrt((r+delta)**2 * np.cos(beta)**2 -
                                                    (2*r*delta + delta**2)) # slope distance
                a_acc += ds*gamma
                alpha = np.arcsin((r+delta)/r * np.sin(beta)) # angle to vertical
                h -= delta
                r -= delta
                t, p, e, n_new, gamma = PropagationP619._get_atmospheric_params(h, rho_s, frequency_MHz)
                delta = 0.0001 + 0.01 * max(h, 0)
                beta = np.arcsin(n/n_new * np.sin(alpha))
                n = n_new

        t, p, e, n, gamma = PropagationP619._get_atmospheric_params(h, rho_s, frequency_MHz)
        delta = .0001 + .01 * max(h, 0)
        r = earth_radius_km + h

        while True:
            ds = np.sqrt(r**2 * np.cos(beta)**2 + 2*r*delta + delta**2) - r * np.cos(beta)
            a_acc += ds * gamma
            alpha = np.arcsin(r/(r+delta) * np.sin(beta))
            h += delta
            if h >= 100:
                break
            r += delta
            t, p, e, n_new, gamma = PropagationP619._get_atmospheric_params(h, rho_s,
                                                                            frequency_MHz)
            beta = np.arcsin(n/n_new * np.sin(alpha))
            n = n_new

        return a_acc

    @staticmethod
    def _get_beam_spreading_att(elevation, altitude, earth_to_space) -> np.array:
        """
        Calculates beam spreading attenuation based on ITU-R P.619, Section 2.4.2

        Parameters
        ----------
            elevation (np.array) : free space elevation (degrees)
            altitude (float) : altitude of earth station (m)
            sat_params (ParametersFss) : parameters of satellite system

        Returns
        -------
            attenuation (np.array): attenuation (dB) with dimensions equal to "elevation"
        """

        # calculate scintillation intensity, based on ITU-R P.618-12
        altitude_km = altitude / 1000

        numerator = .5411 + .07446 * elevation + altitude_km * (.06272 + .0276 * elevation) \
                    + altitude_km ** 2 * .008288
        denominator = (1.728 + .5411 * elevation + .03723 * elevation **2 +
                       altitude_km * (.1815 + .06272 * elevation + .0138 * elevation ** 2) +
                       (altitude_km ** 2) * (.01727 + .008288 * elevation))**2

        attenuation = 10 * np.log10(1 - numerator/denominator)

        if earth_to_space:
            attenuation = -attenuation

        return attenuation

    @staticmethod
    def _get_tropospheric_scintillation(*args, **kwargs)-> np.array:
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
            water_vapour_density = PropagationP619._get_reference_atmosphere_p835 (sat_params.imt_lat_deg,
                                                                                   sat_params.imt_altitude,
                                                                                   sat_params.season)

            # calculate saturation water vapour pressure according to ITU-R P.453-12
            # water coefficients (ice disregarded)
            coef_a = 6.1121
            coef_b = 18.678
            coef_c = 257.14
            coef_d = 234.5
            EF_water = 1 + 1e-4 * (7.2 + pressure * (.032 + 5.9e-6*temperature**2))
            vapour_pressure = water_vapour_density * temperature / 216.7 # eq 10 in P453
            wet_refractivity = (72 * vapour_pressure / temperature
                                + 3.75e5 * vapour_pressure / temperature ** 2)

        sigma_ref = 3.6e-3 + 1e-4 * wet_refractivity
        h_l = 1000
        path_length = 2 * h_l / (np.sqrt(np.sin(elevation_rad)**2 + 2.35e-4) + np.sin(elevation_rad))
        eff_antenna_diameter = .3 * 10 ** (.05 * antenna_gain_dB) / ( np.pi * f_GHz )

        x = 1.22 * eff_antenna_diameter ** 2 * (f_GHz / path_length)
        antenna_averaging_factor = np.sqrt(3.86*(x**2  + 1)**(11/12) *
                                           np.sin(11/6 * np.arctan(1/x)) - 7.08 * x**5/6)
        scintillation_intensity = (sigma_ref * f_GHz**(7/12) * antenna_averaging_factor
                                   / np.sin(elevation_rad)**1.2)

        if time_ratio.lower() == "random":
            time_ratio = np.random.rand(len(elevation_rad))

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

    @staticmethod
    def _get_building_entry_loss(frequency_MHz, elevation, prob = "random",
                                 building_class = "traditional") -> np.array:
        """
        Calculates building loss according to ITU-R P.2109

        Parameters
        ----------
            frequency_MHz (np.array) : carrier frequencies (MHz)
            elevation (np.array) : apparent elevation angles
            prob (np.array / string) : the probability with which the loss is not exceeded;
                                    if "random", then different values are chosen for each user
            building_class (string) : type of construction material, "traditional" or "thermally-efficient"

        Returns
        -------
            array with building loss values with dimensions of elevation

        """

        f_GHz = frequency_MHz/1000

        if prob.lower() == "random":
            prob = np.random.random(elevation.shape)

        if building_class == "traditional":
            r = 12.64
            s = 3.72
            t = .96
            u = 9.6
            v = 2.
            w = 9.1
            x = -3.
            y = 4.5
            z = -2.
        elif building_class == "thermally-efficient":
            r = 28.19
            s = -3.
            t = 8.48
            u = 13.5
            v = 3.8
            w = 27.8
            x = -2.9
            y = 9.4
            z = -2.1
        else:
            error_message = "building_class not supported"
            raise ValueError(error_message)
        c_dB = -3.
        hor_loss = r + s * np.log10(f_GHz) + t * np.log10(f_GHz) ** 2
        angle_correction = .212 * np.abs(elevation)
        mu_1 = hor_loss + angle_correction
        mu_2 = w + x * np.log10(f_GHz)
        sigma_1 = u + v * np.log10(f_GHz)
        sigma_2 = y + z * np.log10(f_GHz)

        a_dB = norm.ppf(prob) * sigma_1 + mu_1
        b_dB = norm.ppf(prob) * sigma_2 + mu_2

        a_lin = 10**(a_dB/10)
        b_lin = 10**(b_dB/10)
        c_lin = 10**(c_dB/10)

        loss = 10*np.log10(a_lin + b_lin + c_lin)

        return loss


    def get_loss(self, *args, **kwargs) -> np.array:
        """
        Calculates path loss for earth-space link

        Parameters
        ----------
            distance_3D (np.array) : distances between stations [m]
            frequency (np.array) : center frequencies [MHz]
            indoor_stations (np.array) : array indicating stations that are indoor
            elevation (dict) : elevation["apparent"](array): apparent elevation angles (degrees)
                               elevation["free_space"](array): free-space elevation angles (degrees)
            sat_params (ParametersFss) : parameters of satellite system
            single_entry (bool): True for single-entry interference, False for multiple-entry (default = False)
            number_of_sectors (int): number of sectors in a node (default = 1)

        Returns
        -------
            array with path loss values with dimensions of distance_3D

        """
        d = kwargs["distance_3D"]
        f = kwargs["frequency"]
        indoor_stations = kwargs["indoor_stations"]
        elevation = kwargs["elevation"]
        sat_params = kwargs["sat_params"]
        earth_to_space = kwargs["earth_to_space"]
        earth_station_antenna_gain = kwargs["earth_station_antenna_gain"]
        single_entry = kwargs.pop("single_entry", False)
        number_of_sectors = kwargs["number_of_sectors"]

        free_space_loss = self.free_space.get_loss(distance_3D=d, frequency=f)

        freq_set = np.unique(f)
        if len(freq_set) > 1:
            error_message = "different frequencies not supported in P619"
            raise ValueError(error_message)

        atmospheric_gasses_loss = self._get_atmospheric_gasses_loss(frequency_MHz=freq_set,
                                                                    apparent_elevation=np.mean(elevation["apparent"]),
                                                                    sat_params=sat_params)
        beam_spreading_attenuation = self._get_beam_spreading_att(elevation["free_space"],
                                                                  sat_params.imt_altitude,
                                                                  earth_to_space)
        diffraction_loss = 0

        if single_entry:
            elevation_sectors = np.repeat(elevation["free_space"], number_of_sectors)
            tropo_scintillation_loss = \
                self._get_tropospheric_scintillation(elevation = elevation_sectors,
                                                     antenna_gain_dB = earth_station_antenna_gain,
                                                     frequency_MHz = freq_set,
                                                     sat_params = sat_params)

            loss = (free_space_loss + self.depolarization_loss +
                    atmospheric_gasses_loss + beam_spreading_attenuation + diffraction_loss)
            loss = np.repeat(loss, number_of_sectors, 1) + tropo_scintillation_loss
        else:
            clutter_loss = np.maximum(0, self.clutter.get_loss(frequency=f,
                                                               elevation=elevation["free_space"],
                                                               station_type=StationType.FSS_SS))
            building_loss = self._get_building_entry_loss(f, elevation["apparent"])* indoor_stations

            loss = (free_space_loss + clutter_loss + building_loss + self.polarization_mismatch_loss +
                    atmospheric_gasses_loss + beam_spreading_attenuation + diffraction_loss)
            loss = np.repeat(loss, number_of_sectors, 1)

        return loss

if __name__ == '__main__':
    from sharc.parameters.parameters import Parameters

    params = Parameters()

    propagation_path = os.getcwd()
    sharc_path = os.path.dirname(propagation_path)
    param_file = os.path.join(sharc_path, "parameters", "parameters.ini")

    params.set_file_name(param_file)
    params.read_params()

    sat_params = params.fss_ss

    propagation = PropagationP619()

    ##########################
    # Plot specific attenuation
    # compare with benchmark from ITU-R P-676-11 Figs. 1 and 2
    temperature = 15 + 273.15  # K
    vapour_density = 7.5  # g/m**3
    pressure_hPa = 1013.25
    vapour_pressure_hPa = vapour_density * temperature / 216.7

    # generate plot
    f_GHz_vec = range(1, 1000)
    specific_att = np.zeros(len(f_GHz_vec))

    print("Plotting specific attenuation:")

    for index in range(len(f_GHz_vec)):
        specific_att[index] = propagation._get_specific_attenuation(pressure_hPa,
                                                                    vapour_pressure_hPa,
                                                                    temperature,
                                                                    float(f_GHz_vec[index]) * 1000)
    plt.figure()
    plt.semilogy(f_GHz_vec, specific_att)
    plt.xlabel('frequency(GHz)')
    plt.xlabel('Specific attenuation (dB/km)')
    plt.title("Atmospheric Specific Attenuation")

    ##########################
    # Plot atmospheric loss
    # compare with benchmark from ITU-R P-619 Fig. 3
    frequency_MHz = 30000.
    sat_params.imt_altitude = 1000

    #apparent_elevation = range(-1, 90, 2)
    apparent_elevation = range(10, 90, 2)

    loss_2_5 = np.zeros(len(apparent_elevation))
    loss_12_5 = np.zeros(len(apparent_elevation))

    print("Plotting atmospheric loss:")
    for index in range(len(apparent_elevation)):
        print("\tApparent Elevation: {} degrees".format(apparent_elevation[index]))

        surf_water_vapour_density = 2.5
        loss_2_5[index] = propagation._get_atmospheric_gasses_loss(frequency_MHz=frequency_MHz,
                                                                   apparent_elevation=apparent_elevation[index],
                                                                   sat_params=sat_params,
                                                                   surf_water_vapour_density=surf_water_vapour_density)
        surf_water_vapour_density = 12.5
        loss_12_5[index] = propagation._get_atmospheric_gasses_loss(frequency_MHz=frequency_MHz,
                                                                    apparent_elevation=apparent_elevation[index],
                                                                    sat_params=sat_params,
                                                                    surf_water_vapour_density=surf_water_vapour_density)

    plt.figure()
    plt.semilogy(apparent_elevation, loss_2_5, label='2.5 g/m^3')
    plt.semilogy(apparent_elevation, loss_12_5, label='12.5 g/m^3')

    plt.grid(True)


    plt.xlabel("apparent elevation (deg)")
    plt.ylabel("Loss (dB)")
    plt.title("Atmospheric Gasses Attenuation")
    plt.legend()

    altitude_vec = np.arange(0, 6.1, .5) * 1000
    elevation_vec = np.array([0, .5, 1, 2, 3, 5])
    attenuation = np.empty([len(altitude_vec), len(elevation_vec)])

    #################################
    # Plot beam spreading attenuation
    # compare with benchmark from ITU-R P-619 Fig. 7

    earth_to_space = False
    print("Plotting beam spreading attenuation:")

    plt.figure()
    for index in range(len(altitude_vec)):
        attenuation[index, :] = propagation._get_beam_spreading_att(elevation_vec,
                                                                    altitude_vec[index],
                                                                    earth_to_space)

    handles = plt.plot(altitude_vec / 1000, np.abs(attenuation))
    plt.xlabel("altitude (km)")
    plt.ylabel("Attenuation (dB)")
    plt.title("Beam Spreading Attenuation")

    for line_handle, elevation in zip(handles, elevation_vec):
        line_handle.set_label("{}deg".format(elevation))

    plt.legend(title="Elevation")

    plt.grid(True)

    #################################
    # Plot troposcatter scintillation attenuation
    # compare with benchmark from ITU-R P-619 Fig. 8

    percentage_fading_exceeded = 10 ** np.arange(-2, 1.1, .1)
    antenna_gain = 0.
    frequency_MHz = 30000.
    wet_refractivity = 42.5

    elevation_vec = np.array([5., 10., 20., 35., 90.])

    print("Plotting troposcatter scintillation attenuation:")

    plt.figure()
    for elevation in elevation_vec:
        attenuation = propagation._get_tropospheric_scintillation(elevation=elevation,
                                                                  frequency_MHz=frequency_MHz,
                                                                  antenna_gain_dB=antenna_gain,
                                                                  time_ratio=1 - (percentage_fading_exceeded / 100),
                                                                  wet_refractivity=wet_refractivity)
        plt.semilogx(percentage_fading_exceeded, attenuation, label="{} deg".format(elevation))

    percentage_gain_exceeded = 10 ** np.arange(-2, 1.1, .1)
    for elevation in elevation_vec:
        attenuation = propagation._get_tropospheric_scintillation(elevation=elevation,
                                                                  frequency_MHz=frequency_MHz,
                                                                  antenna_gain_dB=antenna_gain,
                                                                  time_ratio=percentage_gain_exceeded / 100,
                                                                  wet_refractivity=wet_refractivity)
        plt.loglog(percentage_gain_exceeded, np.abs(attenuation), ':', label="{} deg".format(elevation))

    plt.legend(title='elevation')
    plt.grid(True)

    plt.show()
