# -*- coding: utf-8 -*-
import numpy as np
import typing

from sharc.antenna.antenna_element_imt_m2101 import AntennaElementImtM2101


# Maybe we should consider using a single array to also make the subarrays
# Equations are basically the same
# To make it possible, we should probably update the Beamforming so that, when using subarray:
# - impl doesn't transform coordinates to local twice
# - impl allows fixed eletrical downtilt (no beamforming in this case?)
# - parameters are more clearly defined and passed from Parameters class
class AntennaSubarrayIMT(object):
    """
    Implements a Subarray "Element" for IMT antenna.
    Subarray implementation defined in R23-WP5D-C-0413, Annex 4.2, Table 8.
    Most equations are already implemented in SHARC beamforming antenna, and M2101 element

    Attributes
    ----------
        element: Antenna element to form subarray
        eletrical: Antenna element to form subarray
        n_rows: How many rows subarray has
    """

    def __init__(
        self,
        *,
        element: AntennaElementImtM2101,
        eletrical_downtilt: float,
        n_rows: int,
        element_vert_spacing: float
    ):
        """
        Constructs an AntennaElementImt object.

        Parameters
        ---------
        """
        self.element = element
        self.eletrical_downtilt = eletrical_downtilt
        self.n_rows = n_rows
        # self.n_columns = 1
        self.dv_sub = element_vert_spacing

    def _super_position_vector(
        self,
        theta: float,
        n_rows: int,
        dv_sub: float
    ) -> np.array:
        """
        Calculates super position vector.
        Angles are in the local coordinate system.

        Parameters
        ----------
            theta (float): elevation angle [degrees]
            phi (float): azimuth angle [degrees]

        Returns
        -------
            v_vec (np.array): superposition vector
        """
        r_theta = np.deg2rad(theta)

        m = np.arange(n_rows) + 1

        v_n = np.exp(
            1.0j * 2 * np.pi * (m[: np.newaxis] - 1) * dv_sub * np.cos(r_theta)
        )

        return v_n

    def _weight_vector(self, eletrical_downtilt: float, n_rows: int, dv_sub: float) -> np.array:
        """
        Calculates weight/sub-array excitation.
        Angles are in the local coordinate system.

        Parameters
        ----------
            phi_tilt (float): eletrical horizontal steering [degrees]
            theta_tilt (float): eletrical down-tilt steering [degrees]

        Side Effect
        -------
            w_vec (np.array): weighting vector
        Returns
        -------
            None
        """
        m = np.arange(self.n_rows) + 1

        return 1 / np.sqrt(n_rows) * np.exp(
            1.0j * 2 * np.pi * (m - 1) * dv_sub * np.sin(np.deg2rad(eletrical_downtilt))
        )

    def _calculate_single_dir_gain(self, phi: float, theta: float):
        """
        There is no beamforming, this is just for better code
        """
        elem_g = self.element.element_pattern(
            phi, theta
        )

        v_vec = self._super_position_vector(theta, self.n_rows, self.dv_sub)

        w_vec = self._weight_vector(self.eletrical_downtilt, self.n_rows, self.dv_sub)

        array_g = 10 * np.log10(abs(np.sum(np.multiply(v_vec, w_vec)))**2)

        return array_g + elem_g

    def calculate_gain(
        self,
        phi_arr: typing.Union[np.array, float],
        theta_arr: typing.Union[np.array, float]
    ) -> typing.Union[np.array, float]:
        """
        Calculates the subarray radiation pattern gain.
        Assumes the angles are already received as local coordinates

        Parameters
        ----------
            theta (np.array): elevation angle [degrees]
            phi (np.array): azimuth angle [degrees]

        Returns
        -------
            gain (np.array): element radiation pattern gain value [dBi]
        """
        if isinstance(phi_arr, float) and isinstance(theta_arr, float):
            return self._calculate_single_dir_gain(phi_arr, theta_arr)

        n_direct = len(phi_arr)

        gains = np.zeros(n_direct)

        for i in range(n_direct):
            gains[i] = self._calculate_single_dir_gain(phi_arr[i], theta_arr[i])

        return gains


if __name__ == '__main__':
    from sharc.parameters.imt.parameters_antenna_imt import ParametersAntennaImt
    from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt, PlotAntennaPattern

    figs_dir = "figs/"

    bs_param = ParametersAntennaImt()
    bs2_param = ParametersAntennaImt()
    bs_param.adjacent_antenna_model = "SINGLE_ELEMENT"
    bs2_param.adjacent_antenna_model = "SINGLE_ELEMENT"

    bs_param.normalization = False
    bs2_param.normalization = False
    bs_param.normalization_file = 'beamforming_normalization\\bs_indoor_norm.npz'
    bs2_param.normalization_file = 'beamforming_normalization\\bs2_norm.npz'
    bs_param.minimum_array_gain = -200
    bs2_param.minimum_array_gain = -200

    bs_param.element_pattern = "M2101"
    bs_param.element_max_g = 6.4
    bs_param.element_phi_3db = 90
    bs_param.element_theta_3db = 65
    bs_param.element_am = 30
    bs_param.element_sla_v = 30
    bs_param.n_rows = 16
    bs_param.n_columns = 8
    bs_param.subarray.is_enabled = True
    bs_param.subarray.element_vert_spacing = 0.7
    bs_param.subarray.eletrical_downtilt = 3
    bs_param.subarray.n_rows = 3
    # bs_param.n_rows = 8
    # bs_param.n_columns = 16
    bs_param.element_horiz_spacing = 0.5
    bs_param.element_vert_spacing = 2.1
    bs_param.multiplication_factor = 12
    bs_param.downtilt = 0

    bs2_param.element_pattern = "M2101"
    bs2_param.element_max_g = 6.4
    bs2_param.element_phi_3db = 90
    bs2_param.element_theta_3db = 65
    bs2_param.element_am = 30
    bs2_param.element_sla_v = 30
    bs2_param.n_rows = 3
    bs2_param.n_columns = 1
    bs2_param.element_horiz_spacing = 0.5
    bs2_param.element_vert_spacing = 0.7
    bs2_param.multiplication_factor = 12

    plot = PlotAntennaPattern(figs_dir)

    # Plot BS TX radiation patterns
    par = bs_param.get_antenna_parameters()
    bs_array = AntennaBeamformingImt(par, 0, 0, bs_param.subarray)
    f = plot.plot_element_pattern(bs_array, "BS", "ELEMENT")
    # f.savefig(figs_dir + "BS_element.pdf", bbox_inches='tight')
    f = plot.plot_element_pattern(bs_array, "BS", "SUBARRAY")
    f = plot.plot_element_pattern(bs_array, "BS", "ARRAY")
    # f.savefig(figs_dir + "BS_array.pdf", bbox_inches='tight')

    # Plot UE TX radiation patterns
    par = bs2_param.get_antenna_parameters()
    bs2_array = AntennaBeamformingImt(par, 0, 0, bs2_param.subarray)
    # plot.plot_element_pattern(bs2_array, "BS 2", "ELEMENT")
    # plot.plot_element_pattern(bs2_array, "BS 2", "ARRAY")

    print('END')
