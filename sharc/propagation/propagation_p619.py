# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 15:35:00 2017

@author: andre barreto
"""


import os
import csv
from scipy.interpolate import interp1d
import numpy as np
from multipledispatch import dispatch
from sharc.station_manager import StationManager
from sharc.parameters.parameters import Parameters
from sharc.propagation.propagation import Propagation
from sharc.propagation.propagation_free_space import PropagationFreeSpace
from sharc.propagation.propagation_clutter_loss import PropagationClutterLoss
from sharc.propagation.propagation_building_entry_loss import PropagationBuildingEntryLoss
from sharc.propagation.atmosphere import ReferenceAtmosphere
from sharc.support.enumerations import StationType
from sharc.propagation.scintillation import Scintillation
from sharc.parameters.constants import EARTH_RADIUS


class PropagationP619(Propagation):
    """
    Implements the earth-to-space channel model from ITU-R P.619

    Public methods:
        get_loss: Calculates path loss for earth-space link
    """

    def __init__(
        self,
        random_number_gen: np.random.RandomState,
        space_station_alt_m: float,
        earth_station_alt_m: float,
        earth_station_lat_deg: float,
        earth_station_long_diff_deg: float,
        season: str,
        mean_clutter_height: str,
        below_rooftop: float
    ):
        """Implements the earth-to-space channel model from ITU-R P.619

        Parameters
        ----------
        random_number_gen : np.random.RandomState
            randon number generator
        space_station_alt_m : float
            The space station altitude in meters
        earth_station_alt_m : float
            The Earth station altitude in meters
        earth_station_lat_deg : float
            The Earth station latitude in degrees
        earth_station_long_diff_deg : float
             The difference between longitudes of earth-station and and space-sation.
             (positive if space-station is to the East of earth-station)
        season: str
            Possible values are "SUMMER" or "WINTER".
        """

        super().__init__(random_number_gen)

        self.is_earth_space_model = True
        self.clutter = PropagationClutterLoss(self.random_number_gen)
        self.free_space = PropagationFreeSpace(self.random_number_gen)
        self.building_entry = PropagationBuildingEntryLoss(
            self.random_number_gen,
        )
        self.scintillation = Scintillation(self.random_number_gen)
        self.atmosphere = ReferenceAtmosphere()

        self.depolarization_loss = 0
        self.polarization_mismatch_loss = 0
        self.elevation_has_atmospheric_loss = []
        self.freq_has_atmospheric_loss = []
        self.surf_water_dens_has_atmospheric_loss = []
        self.atmospheric_loss = []
        self.elevation_delta = .01

        self.space_station_alt_m = space_station_alt_m
        self.earth_station_alt_m = earth_station_alt_m
        self.earth_station_lat_deg = earth_station_lat_deg
        self.earth_station_long_diff_deg = earth_station_long_diff_deg
        self.mean_clutter_height = mean_clutter_height
        self.below_rooftop = below_rooftop
        if season.upper() not in ["SUMMER", "WINTER"]:
            raise ValueError(
                f"PropagationP619: Invalid value for parameter season - {season}. Possible values are \"SUMMER\", \"WINTER\".",
            )
        self.season = season
        # Load city name based on latitude
        self.lookup_table = False
        self.city_name = self._get_city_name_by_latitude()
        if self.city_name != "Unknown":
            self.lookup_table = True

    def _get_city_name_by_latitude(self):
        localidades_file = os.path.join(
            os.path.dirname(__file__), 'Dataset/localidades.csv',
        )
        with open(localidades_file, mode='r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                if abs(float(row['Latitude']) - self.earth_station_lat_deg) < 1e-5:
                    return row['Cidade']
        return 'Unknown'

    def _get_atmospheric_gasses_loss(self, *args, **kwargs) -> float:
        """
        Calculates atmospheric gasses loss based on ITU-R P.619, Attachment C

        Parameters
        ----------
            frequency_MHz (float) : center frequencies [MHz]
            apparent_elevation (float) : apparent elevation angle at the Earth-based station (degrees)
        Returns
        -------
            path_loss (float): scalar with atmospheric loss
        """
        frequency_MHz = kwargs["frequency_MHz"]
        apparent_elevation = kwargs["apparent_elevation"]
        lookupTable = kwargs.pop("lookupTable", True)
        surf_water_vapour_density = kwargs.pop(
            "surf_water_vapour_density", False,
        )

        if lookupTable and self.city_name != 'Unknown':
            # Define the path to the CSV file
            output_dir = os.path.join(os.path.dirname(__file__), 'Dataset')
            csv_file = os.path.join(
                output_dir, f'{self.city_name}_{int(frequency_MHz)}_{int(self.earth_station_alt_m)}m.csv',
            )
            if os.path.exists(csv_file):
                elevations = []
                losses = []
                with open(csv_file, mode='r') as file:
                    reader = csv.reader(file)
                    next(reader)  # Skip header
                    for row in reader:
                        elevations.append(float(row[0]))
                        losses.append(float(row[1]))
                interpolation_function = interp1d(
                    elevations, losses, kind='linear', fill_value='extrapolate',
                )
                return interpolation_function(apparent_elevation)

        earth_radius_km = EARTH_RADIUS / 1000
        a_acc = 0.  # accumulated attenuation (in dB)
        h = self.earth_station_alt_m / 1000  # ray altitude in km
        beta = (90 - abs(apparent_elevation)) * np.pi / 180.  # incidence angle

        if not surf_water_vapour_density:
            _, _, surf_water_vapour_density = self.atmosphere.get_reference_atmosphere_p835(
                self.earth_station_lat_deg, 0, season=self.season,
            )

        if len(self.elevation_has_atmospheric_loss):
            elevation_diff = np.abs(
                apparent_elevation - np.array(self.elevation_has_atmospheric_loss),
            )
            elevation_diff = np.abs(
                apparent_elevation - np.array(self.elevation_has_atmospheric_loss),
            )
            indices = np.where(elevation_diff <= self.elevation_delta)
            if indices[0].size:
                index = np.argmin(elevation_diff)
                if self.freq_has_atmospheric_loss[index] == frequency_MHz and \
                   self.surf_water_dens_has_atmospheric_loss[index] == surf_water_vapour_density:
                    return self.atmospheric_loss[index]

        rho_s = surf_water_vapour_density * \
            np.exp(h / 2)  # water vapour density at h
        if apparent_elevation < 0:
            t, p, e, n, gamma = self.atmosphere.get_atmospheric_params(
                h, rho_s, frequency_MHz,
            )
            delta = .0001 + 0.01 * max(h, 0)  # layer thickness
            r = earth_radius_km + h - delta  # radius of lower edge
            while True:
                m = (r + delta) * np.sin(beta) - r
                if m >= 0:
                    dh = 2 * np.sqrt(
                        2 * r * (delta - m) +
                        delta ** 2 - m ** 2,
                    )  # horizontal path
                    a_acc += dh * gamma
                    break
                ds = (r + delta) * np.cos(beta) - np.sqrt((r + delta) **
                                                          2 * np.cos(beta) ** 2 - (2 * r * delta + delta ** 2))
                a_acc += ds * gamma
                # angle to vertical
                alpha = np.arcsin((r + delta) / r * np.sin(beta))
                h -= delta
                r -= delta
                t, p, e, n_new, gamma = self.atmosphere.get_atmospheric_params(
                    h, rho_s, frequency_MHz,
                )
                delta = 0.0001 + 0.01 * max(h, 0)
                beta = np.arcsin(n / n_new * np.sin(alpha))
                n = n_new

        t, p, e, n, gamma = self.atmosphere.get_atmospheric_params(
            h, rho_s, frequency_MHz,
        )
        delta = .0001 + .01 * max(h, 0)
        r = earth_radius_km + h

        while True:
            ds = np.sqrt(
                r ** 2 * np.cos(beta) ** 2 + 2 * r *
                delta + delta ** 2,
            ) - r * np.cos(beta)
            a_acc += ds * gamma
            alpha = np.arcsin(r / (r + delta) * np.sin(beta))
            h += delta
            if h >= 100:
                break
            r += delta
            t, p, e, n_new, gamma = self.atmosphere.get_atmospheric_params(
                h, rho_s, frequency_MHz,
            )
            beta = np.arcsin(n / n_new * np.sin(alpha))
            n = n_new

        self.atmospheric_loss.append(a_acc)
        self.elevation_has_atmospheric_loss.append(apparent_elevation)
        self.freq_has_atmospheric_loss.append(frequency_MHz)
        self.surf_water_dens_has_atmospheric_loss.append(
            surf_water_vapour_density,
        )

        return a_acc

    @staticmethod
    def _get_beam_spreading_att(elevation, altitude, earth_to_space) -> np.array:
        """
        Calculates beam spreading attenuation based on ITU-R P.619, Section 2.4.2

        Parameters
        ----------
            elevation (np.array) : free space elevation (degrees)
            altitude (float) : altitude of earth station (m)

        Returns
        -------
            attenuation (np.array): attenuation (dB) with dimensions equal to "elevation"
        """

        # calculate scintillation intensity, based on ITU-R P.618-12
        altitude_km = altitude / 1000

        numerator = .5411 + .07446 * elevation + altitude_km * \
            (.06272 + .0276 * elevation) + altitude_km ** 2 * .008288
        denominator = (
            1.728 + .5411 * elevation + .03723 * elevation ** 2 + altitude_km * (
                .1815 + .06272 *
                elevation + .0138 * elevation ** 2
            ) + (altitude_km ** 2) * (.01727 + .008288 * elevation)
        ) ** 2

        attenuation = 10 * np.log10(1 - numerator / denominator)

        if earth_to_space:
            attenuation = -attenuation

        return attenuation

    @classmethod
    def apparent_elevation_angle(cls, elevation_deg: np.array, space_station_alt_m: float) -> np.array:
        """Calculate apparent elevation angle according to ITU-R P619, Attachment B

        Parameters
        ----------
        elevation_deg : np.array
            free-space elevation angle
        space_station_alt_m : float
            space-station altitude

        Returns
        -------
        np.array
            apparent elevation angle
        """
        elev_angles_rad = np.deg2rad(elevation_deg)
        tau_fs1 = 1.728 + 0.5411 * elev_angles_rad + 0.03723 * elev_angles_rad**2
        tau_fs2 = 0.1815 + 0.06272 * elev_angles_rad + 0.01380 * elev_angles_rad**2
        tau_fs3 = 0.01727 + 0.008288 * elev_angles_rad

        # change in elevation angle due to refraction
        tau_fs_deg = 1 / (
            tau_fs1 + space_station_alt_m * tau_fs2 +
            space_station_alt_m**2 * tau_fs3
        )
        tau_fs = tau_fs_deg / 180. * np.pi

        return np.degrees(elev_angles_rad + tau_fs)

    @dispatch(Parameters, float, StationManager, StationManager, np.ndarray, np.ndarray)
    def get_loss(
        self,
        params: Parameters,
        frequency: float,
        station_a: StationManager,
        station_b: StationManager,
        station_a_gains=None,
        station_b_gains=None,
    ) -> np.array:
        """Wrapper for the get_loss function that satisfies the ABS class interface

        Calculates P.619 path loss between station_a and station_b stations.

        Parameters
        ----------
        station_a : StationManager
            StationManager container representing station_a
        station_b : StationManager
            StationManager container representing station_b
        params : Parameters
            Simulation parameters needed for the propagation class.

        Returns
        -------
        np.array
            Return an array station_a.num_stations x station_b.num_stations with the path loss
            between each station
        """
        distance = station_a.get_3d_distance_to(station_b)
        frequency = frequency * np.ones(distance.shape)
        indoor_stations = np.tile(
            station_b.indoor, (station_a.num_stations, 1),
        )

        # Elevation angles seen from the station on Earth.
        elevation_angles = {}
        if station_a.is_space_station:
            elevation_angles["free_space"] = station_b.get_elevation(station_a)
            earth_station_antenna_gain = station_b_gains
            # if (station_b_gains.shape != distance.shape):
            #     raise ValueError(f"Invalid shape for station_b_gains = {station_b_gains.shape}")
            elevation_angles["apparent"] = self.apparent_elevation_angle(
                elevation_angles["free_space"],
                station_a.height,
            )
            # Transpose it to fit the expected path loss shape
            elevation_angles["free_space"] = np.transpose(elevation_angles["free_space"])
            elevation_angles["apparent"] = np.transpose(elevation_angles["apparent"])
        elif station_b.is_space_station:
            elevation_angles["free_space"] = station_a.get_elevation(station_b)
            earth_station_antenna_gain = station_a_gains
            elevation_angles["apparent"] = self.apparent_elevation_angle(
                elevation_angles["free_space"],
                station_b.height,
            )
        else:
            raise ValueError(
                "PropagationP619: At least one station must be an space station",
            )

        # Determine the Earth-space path direction and if the interferer is single or multiple-entry.
        # The get_loss interface won't tell us who is the interferer, so we need to get it from the parameters.
        is_intra_imt = True if station_a.is_imt_station(
        ) and station_b.is_imt_station() else False
        if is_intra_imt:
            # Intra-IMT
            if params.general.imt_link == "UPLINK":
                is_earth_to_space_link = True
            else:
                is_earth_to_space_link = False
            # Intra-IMT is always multiple-entry interference
            is_single_entry_interf = False
        else:
            # System-IMT
            imt_station, sys_station = (
                station_a, station_b,
            ) if station_a.is_imt_station() else (station_b, station_a)
            if params.imt.interfered_with:
                is_earth_to_space_link = True if imt_station.is_space_station else False
                if sys_station.num_stations > 1:
                    is_single_entry_interf = False
                else:
                    is_single_entry_interf = True
            else:
                is_earth_to_space_link = False if imt_station.is_space_station else True
                is_single_entry_interf = False
        earth_station_height = station_b.height
        mean_clutter_height = self.mean_clutter_height
        below_rooftop = self.below_rooftop
        loss = self.get_loss(
            distance,
            frequency,
            indoor_stations,
            elevation_angles,
            is_earth_to_space_link,
            earth_station_antenna_gain,
            is_single_entry_interf,
            earth_station_height,
            mean_clutter_height,
            below_rooftop,
        )

        return loss

    @dispatch(np.ndarray, np.ndarray, np.ndarray, dict, bool, np.ndarray, bool, np.ndarray, str, int)
    def get_loss(
        self,
        distance: np.array,
        frequency: np.array,
        indoor_stations: np.array,
        elevation: dict,
        earth_to_space: bool,
        earth_station_antenna_gain: np.array,
        single_entry: bool,
        earth_station_height: np.array,
        mean_clutter_height: np.str_,
        below_rooftop: float
    ) -> np.array:
        """
        Calculates path loss for earth-space link

        Parameters
        ----------
            distance (np.array) : 3D distances between stations [m]
            frequency (np.array) : center frequencies [MHz]
            indoor_stations (np.array) : array indicating stations that are indoor
            elevation (dict) : dict containing free-space elevation angles (degrees), and aparent elevation angles
            earth_to_space (bool): whether the ray direction is earth-to-space. Affects beam spreadding Attenuation.
            single_entry (bool): True for single-entry interference, False for
            multiple-entry (default = False)
            earth_station_antenna_gain (np.array): Earth station atenna gain array.

        Returns
        -------
            array with path loss values with dimensions of distance_3D

        """
        free_space_loss = self.free_space.get_free_space_loss(
            frequency=frequency, distance=distance,
        )

        freq_set = np.unique(frequency)
        if len(freq_set) > 1:
            error_message = "different frequencies not supported in P619"
            raise ValueError(error_message)

        atmospheric_gasses_loss = self._get_atmospheric_gasses_loss(
            frequency_MHz=freq_set,
            apparent_elevation=np.mean(elevation["apparent"]),
        )
        beam_spreading_attenuation = self._get_beam_spreading_att(
            elevation["free_space"],
            self.earth_station_alt_m,
            earth_to_space,
        )
        diffraction_loss = 0

        if single_entry:
            tropo_scintillation_loss = self.scintillation.get_tropospheric_attenuation(
                elevation=elevation["free_space"],
                antenna_gain_dB=earth_station_antenna_gain,
                frequency_MHz=freq_set,
                earth_station_alt_m=self.earth_station_alt_m,
                earth_station_lat_deg=self.earth_station_lat_deg,
                season=self.season,
            )

            loss = \
                free_space_loss + self.depolarization_loss + atmospheric_gasses_loss + beam_spreading_attenuation + \
                diffraction_loss + tropo_scintillation_loss
        else:
            clutter_loss = \
                self.clutter.get_loss(
                    frequency=frequency,
                    distance=distance,
                    elevation=elevation["free_space"],
                    station_type=StationType.FSS_SS,
                    earth_station_height=earth_station_height,
                    mean_clutter_height=mean_clutter_height,
                    below_rooftop = below_rooftop,
                )
            building_loss = self.building_entry.get_loss(
                frequency, elevation["apparent"],
            ) * indoor_stations

            loss = (
                free_space_loss + clutter_loss + building_loss + self.polarization_mismatch_loss +
                atmospheric_gasses_loss + beam_spreading_attenuation + diffraction_loss
            )

        return loss


if __name__ == '__main__':

    import matplotlib
    matplotlib.use('Agg')

    import matplotlib.pyplot as plt

    # params = Parameters()

    # propagation_path = os.getcwd()
    # sharc_path = os.path.dirname(propagation_path)
    # param_file = os.path.join(sharc_path, "parameters", "parameters.ini") deprecated

    # params.set_file_name(param_file)
    # params.read_params()

    # sat_params = params.fss_ss

    space_station_alt_m = 20000.0
    earth_station_alt_m = 1000.0
    earth_station_lat_deg = -15.7801
    earth_station_long_diff_deg = 0.0
    season = "SUMMER"

    random_number_gen = np.random.RandomState(101)
    propagation = PropagationP619(
        random_number_gen=random_number_gen,
        space_station_alt_m=space_station_alt_m,
        earth_station_alt_m=earth_station_alt_m,
        earth_station_lat_deg=earth_station_lat_deg,
        earth_station_long_diff_deg=earth_station_long_diff_deg,
        season=season,
    )

    ##########################
    # Plot atmospheric loss
    # compare with benchmark from ITU-R P-619 Fig. 3
    frequency_MHz = 2680.
    earth_station_alt_m = 1000

    apparent_elevation = np.linspace(0, 90, 91)

    loss_false = np.zeros(len(apparent_elevation))
    loss_true = np.zeros(len(apparent_elevation))

    print("Plotting atmospheric loss:")
    for index in range(len(apparent_elevation)):
        print(
            "\tApparent Elevation: {} degrees".format(
                apparent_elevation[index],
            ),
        )

        surf_water_vapour_density = 14.121645412892681
        table = False
        loss_false[index] = propagation._get_atmospheric_gasses_loss(
            frequency_MHz=frequency_MHz,
            apparent_elevation=apparent_elevation[index],
            surf_water_vapour_density=surf_water_vapour_density,
            lookupTable=table,
        )
        table = True
        loss_true[index] = propagation._get_atmospheric_gasses_loss(
            frequency_MHz=frequency_MHz,
            apparent_elevation=apparent_elevation[index],
            surf_water_vapour_density=surf_water_vapour_density,
            lookupTable=table,
        )

    plt.figure()
    plt.semilogy(apparent_elevation, loss_false, label='No Table')
    plt.semilogy(apparent_elevation, loss_true, label='With Table')

    plt.grid(True)
    plt.xlabel("apparent elevation (deg)")
    plt.ylabel("Loss (dB)")
    plt.title("Atmospheric Gasses Attenuation")
    plt.legend()

    # Save the figures instead of showing them
    plt.savefig("atmospheric_gasses_attenuation.png")
