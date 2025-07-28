from sharc.antenna.antenna import Antenna

import numpy as np
import math


class AntennaMultipleTransceiver(Antenna):
    def __init__(
        self,
        *,
        num_beams: int,
        transceiver_radiation_pattern: Antenna,
        azimuths: np.array,
        elevations: np.array
    ):
        if num_beams != len(azimuths) or num_beams != len(elevations):
            raise ValueError(
                "When using AntennaMultipleTransceiver you need to pass each transceiver's elevation and azimuth\n"
                f"You have passed {num_beams} number of beams, but {len(azimuths)} azimuths and {len(elevations)} elevations"
            )

        super().__init__()
        # WARNING: duplicated state here. May lead to problems
        # in case of coordinate transformations. Needs to always be updated
        # along with related StationManager
        # this problem already exists with antenna imt beamforming
        # TODO?: always set elevation and azimuth as antenna attribute?
        # that way the antenna could always just receive the pointing vector

        # transforming 0 at horizon to 0 at zenith
        self.phi = azimuths
        self.theta = 90 - elevations
        self.transceiver_radiation_pattern = transceiver_radiation_pattern

    def calculate_gain(self, *args, **kwargs) -> np.array:
        """
        Calculates the antena gain.
        """
        phi_vec = kwargs["phi_vec"]
        theta_vec = kwargs["theta_vec"]

        off_axis_angles = self.get_off_axis_angle(self.theta, self.phi, theta_vec, phi_vec)

        gains = self.transceiver_radiation_pattern.calculate_gain(
            off_axis_angle_vec=off_axis_angles,
            theta_vec=theta_vec
        )

        return 10 * np.log10(
            np.sum(
                10**(gains / 10),
                axis=0  # only sum each antenna results from the same angle used
            )
        )

    @staticmethod
    def get_off_axis_angle(
        antenna_theta: np.array,
        antenna_phi: np.array,
        obj_theta: np.array,
        obj_phi: np.array
    ):
        """
        Calculates the off-axis angle when comparing the antenna pointing
        and an object relative elevation and azimuth.

        Parameters
        ----------
            antenna_theta:
                1d np.array with theta angles of where each transceiver points
            antenna_phi:
                1d np.array with phi angles of where each transceiver points
            obj_theta:
                2d np.array with phi angles of the line that passes through
                the antenna and the object
            obj_phi:
                2d np.array with theta angles of the line that passes through
                the antenna and the object

        theta angles are relative from the x axis, counter clockwise
            (same as in the rest of simulator)
        phi angles are relative from the z axis, with 0 being at the zenith

        Returns
        -------
            The angle between the lines defined by (theta, phi)
            for each pair of (antenna_theta[i],antenna_phi[i]), (obj_theta[j], antenna_phi[j])
        """
        relative_phi = antenna_phi[:, np.newaxis] - obj_phi
        antenna_theta = antenna_theta[:, np.newaxis]

        rel_cos = (np.cos(np.radians(antenna_theta)) * np.cos(np.radians(obj_theta)) +
            np.sin(np.radians(antenna_theta)) * np.sin(np.radians(obj_theta)) * np.cos(np.radians(relative_phi)))

        rel = np.arccos(
            # sometimes floating point error accumulates enough
            # for phi_cos to fall just out of range of arccos
            np.clip(rel_cos, -1.0, 1.0)
        )

        return np.degrees(rel)

    def rotate_antenna(self):
        """
        TODO: rotate phi angles when transforming base station
        """
        pass
