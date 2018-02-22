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
from sharc.propagation.propagation_building_entry_loss import PropagationBuildingEntryLoss
from sharc.propagation.atmosphere import ReferenceAtmosphere
from sharc.support.enumerations import StationType
from sharc.propagation.scintillation import Scintillation


class PropagationP619(Propagation):
    """
    Implements the earth-to-space channel model from ITU-R P.619

    Public methods:
        get_loss: Calculates path loss for earth-space link
    """

    def __init__(self, random_number_gen: np.random.RandomState):
        super().__init__(random_number_gen)

        self.clutter = PropagationClutterLoss(self.random_number_gen)
        self.free_space = PropagationFreeSpace(self.random_number_gen)
        self.building_entry = PropagationBuildingEntryLoss(self.random_number_gen)
        self.scintillation = Scintillation(self.random_number_gen)
        self.atmosphere = ReferenceAtmosphere()

        self.depolarization_loss = 1.5
        self.polarization_mismatch_loss = 3.
        self.elevation_has_atmospheric_loss = []
        self.freq_has_atmospheric_loss = []
        self.surf_water_dens_has_atmospheric_loss = []
        self.atmospheric_loss = []
        self.elevation_delta = .01


    def _get_atmospheric_gasses_loss(self, *args, **kwargs) -> float:
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
                self.atmosphere.get_reference_atmosphere_p835(sat_params.imt_lat_deg,
                                                               0, season="summer")

        # first, check if atmospheric loss was already calculated
        if len(self.elevation_has_atmospheric_loss):
            elevation_diff = np.abs(apparent_elevation - np.array(self.elevation_has_atmospheric_loss))
            indices = np.where(elevation_diff <= self.elevation_delta)
            if indices[0].size:
                index = np.argmin(elevation_diff)
                if self.freq_has_atmospheric_loss[index] == frequency_MHz and \
                   self.surf_water_dens_has_atmospheric_loss[index] == surf_water_vapour_density:
                    loss = self.atmospheric_loss[index]
                    return loss

        rho_s = surf_water_vapour_density * np.exp(h/2) # water vapour density at h
        if apparent_elevation < 0:
            # get temperature (t), dry-air pressure (p), water-vapour pressure (e),
            #     refractive index (n) and specific attenuation (gamma)
            t, p, e, n, gamma = self.atmosphere.get_atmospheric_params(h, rho_s, frequency_MHz)
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
                t, p, e, n_new, gamma = self.atmosphere.get_atmospheric_params(h, rho_s, frequency_MHz)
                delta = 0.0001 + 0.01 * max(h, 0)
                beta = np.arcsin(n/n_new * np.sin(alpha))
                n = n_new

        t, p, e, n, gamma = self.atmosphere.get_atmospheric_params(h, rho_s, frequency_MHz)
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
            t, p, e, n_new, gamma = self.atmosphere.get_atmospheric_params(h, rho_s,
                                                                           frequency_MHz)
            beta = np.arcsin(n/n_new * np.sin(alpha))
            n = n_new

        self.atmospheric_loss.append(a_acc)
        self.elevation_has_atmospheric_loss.append(apparent_elevation)
        self.freq_has_atmospheric_loss.append(frequency_MHz)
        self.surf_water_dens_has_atmospheric_loss.append(surf_water_vapour_density)

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
                self.scintillation.get_tropospheric_attenuation(elevation = elevation_sectors,
                                                                antenna_gain_dB = earth_station_antenna_gain,
                                                                frequency_MHz = freq_set,
                                                                sat_params = sat_params)

            loss = (free_space_loss + self.depolarization_loss +
                    atmospheric_gasses_loss + beam_spreading_attenuation + diffraction_loss)
            loss = np.repeat(loss, number_of_sectors, 1) + tropo_scintillation_loss
        else:
            clutter_loss = self.clutter.get_loss(frequency=f, distance=d,
                                                 elevation=elevation["free_space"],
                                                 station_type=StationType.FSS_SS)
            building_loss = self.building_entry.get_loss(f, elevation["apparent"]) * indoor_stations

            loss = (free_space_loss + clutter_loss + building_loss + self.polarization_mismatch_loss +
                    atmospheric_gasses_loss + beam_spreading_attenuation + diffraction_loss)
            loss = np.repeat(loss, number_of_sectors, 1)

        return loss

if __name__ == '__main__':
    from sharc.parameters.parameters import Parameters
    import matplotlib.pyplot as plt
    import os

    params = Parameters()

    propagation_path = os.getcwd()
    sharc_path = os.path.dirname(propagation_path)
    param_file = os.path.join(sharc_path, "parameters", "parameters.ini")

    params.set_file_name(param_file)
    params.read_params()

    sat_params = params.fss_ss

    propagation = PropagationP619()

    ##########################
    # Plot atmospheric loss
    # compare with benchmark from ITU-R P-619 Fig. 3
    frequency_MHz = 30000.
    sat_params.imt_altitude = 1000

    apparent_elevation = range(-1, 90, 2)

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


    plt.show()
