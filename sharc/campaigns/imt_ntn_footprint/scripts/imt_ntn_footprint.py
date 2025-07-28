import numpy as np
import matplotlib.pyplot as plt

from sharc.topology.topology_ntn import TopologyNTN
from sharc.station_factory import StationFactory
from sharc.parameters.imt.parameters_imt import ParametersImt
from sharc.parameters.imt.parameters_antenna_imt import ParametersAntennaImt
from sharc.parameters.parameters_mss_ss import ParametersMssSs


if __name__ == "__main__":

    # Input parameters for MSS_SS
    param_mss = ParametersMssSs()
    param_mss.frequency = 2100.0  # MHz
    param_mss.bandwidth = 20.0  # MHz
    param_mss.altitude = 500e3  # meters
    param_mss.azimuth = 0
    param_mss.elevation = 90  # degrees
    param_mss.cell_radius = 19e3  # meters
    param_mss.intersite_distance = param_mss.cell_radius * np.sqrt(3)
    param_mss.num_sectors = 19
    param_mss.antenna_gain = 30  # dBi
    param_mss.antenna_3_dB_bw = 4.4127
    param_mss.antenna_l_s = 20  # in dB
    # Parameters used for the S.1528 antenna
    param_mss.antenna_pattern = "ITU-R-S.1528-Taylor"
    # param_mss.antenna_pattern = "ITU-R-S.1528-LEO"
    wavelength = 3e8 / (param_mss.frequency * 1e6)  # in meters
    # assuming 7dB roll-off factor and circular antenna
    l_r = 0.74 * wavelength / np.sin(param_mss.antenna_3_dB_bw / 2)  # in meters
    l_t = l_r
    param_mss.antenna_s1528.set_external_parameters(frequency=param_mss.frequency,
                                                    bandwidth=param_mss.bandwidth,
                                                    antenna_gain=param_mss.antenna_gain,
                                                    antenna_l_s=param_mss.antenna_l_s,
                                                    antenna_3_dB_bw=param_mss.antenna_3_dB_bw,
                                                    l_r=l_r,
                                                    l_t=l_t)
    beam_idx = 15  # beam index used for gain analysis

    seed = 100
    rng = np.random.RandomState(seed)

    # Parameters used for IMT-NTN and UE distribution
    param_imt = ParametersImt()
    param_imt.topology.type = "NTN"
    param_imt.ue.azimuth_range = (-180, 180)
    param_imt.bandwidth = 10  # MHz
    param_imt.frequency = 2100  # MHz
    param_imt.spurious_emissions = -13  # dB
    param_imt.ue.distribution_azimuth = "UNIFORM"
    param_imt.ue.k = 1000

    ntn_topology = TopologyNTN(param_mss.intersite_distance,
                               param_mss.cell_radius,
                               param_mss.altitude,
                               param_mss.azimuth,
                               param_mss.elevation,
                               param_mss.num_sectors)

    ntn_topology.calculate_coordinates()
    param_ue_ant = ParametersAntennaImt()
    ntn_ue = StationFactory.generate_imt_ue_outdoor(
        param_imt, param_ue_ant, rng, ntn_topology)

    ntn_ue.active = np.ones(ntn_ue.num_stations, dtype=bool)
    ntn_bs = StationFactory.generate_mss_ss(param_mss)
    phi, theta = ntn_bs.get_pointing_vector_to(ntn_ue)
    station_1_active = np.where(ntn_bs.active)[0]
    station_2_active = np.where(ntn_ue.active)[0]
    beams_idx = np.zeros(len(station_2_active), dtype=int)
    off_axis_angle = ntn_bs.get_off_axis_angle(ntn_ue)
    gains = np.zeros(phi.shape)
    for k in station_1_active:
        gains[k, station_2_active] = \
            ntn_bs.antenna[k].calculate_gain(
                off_axis_angle_vec=off_axis_angle[k, station_2_active], theta_vec=theta[k, station_2_active])
                # phi=off_axis_angle[k, station_2_active], theta=theta[k, station_2_active])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # ax.set_xlim([-200, 200])
    # ax.set_ylim([-200, 200])
    ntn_topology.plot_3d(ax, False)  # Plot the 3D topology
    im = ax.scatter(xs=ntn_ue.x / 1000, ys=ntn_ue.y / 1000,
                    c=gains[beam_idx] - np.max(param_mss.antenna_gain), vmin=-50, cmap='jet')
    ax.view_init(azim=0, elev=90)
    fig.colorbar(im, label='Normalized antenna gain (dBi)')

    plt.show()
    exit()
