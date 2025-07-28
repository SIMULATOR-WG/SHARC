# Helper script that calculates the coordintaes of the NTN footprint based on the distance between borders with IMT-TN
from sharc.parameters.parameters_mss_ss import ParametersMssSs
from sharc.topology.topology_ntn import TopologyNTN
from sharc.topology.topology_macrocell import TopologyMacrocell
import matplotlib.pyplot as plt
import math
import numpy as np

if __name__ == "__main__":

    # distance from topology boarders in meters
    border_distances_array = np.array(
        [0, 20e3, 50e3, 100e3, 200e3, 300e3, 400e3, 500e3, 600e3, 700e3, 1000e3])

    # Index in border_distances_array used for plotting
    dist_idx = 5

    # Input parameters for MSS_SS
    param_mss = ParametersMssSs()
    param_mss.frequency = 2160  # MHz
    # param_mss.altitude = 500e3   # meters
    param_mss.altitude = 1200e3   # meters
    param_mss.azimuth = 0
    param_mss.elevation = 90  # degrees
    # param_mss.cell_radius = 19e3  # meters
    param_mss.cell_radius = 45e3  # meters
    param_mss.intersite_distance = param_mss.cell_radius * np.sqrt(3)
    param_mss.num_sectors = 19

    # Input paramters for IMT-TN macrocell topology
    macro_cell_radius = 500  # meters

    macrocell_num_clusters = 1

    #######################################################

    # calculate the position of the NTN footprint center
    ntn_footprint_left_edge = - 4 * param_mss.cell_radius
    ntn_footprint_radius = 5 * param_mss.cell_radius * np.sin(np.pi / 3)
    macro_topology_radius = 4 * macro_cell_radius
    ntn_footprint_x_offset = macro_topology_radius + ntn_footprint_radius + border_distances_array
    param_mss.x = ntn_footprint_x_offset[dist_idx]

    ntn_topology = TopologyNTN(param_mss.intersite_distance,
                               param_mss.cell_radius,
                               param_mss.altitude,
                               param_mss.azimuth,
                               param_mss.elevation,
                               param_mss.num_sectors)
    ntn_topology.calculate_coordinates()
    ntn_topology.x = ntn_topology.x + param_mss.x

    seed = 100
    rng = np.random.RandomState(seed)

    intersite_distance = macro_cell_radius * 3 / 2
    macro_topology = TopologyMacrocell(intersite_distance, macrocell_num_clusters)
    macro_topology.calculate_coordinates()

    ntn_circle = plt.Circle((param_mss.x, 0), radius=ntn_footprint_radius, fill=False, color='green')
    macro_cicle = plt.Circle((0, 0), radius=macro_topology_radius, fill=False, color='red')

    ## Print simulation information
    print("\n######## Simulation scenario parameters ########")
    for i, d in enumerate(border_distances_array):
        print(f"Satellite altitude = {param_mss.altitude / 1000} km")
        print(f"Beam footprint radius = {param_mss.cell_radius / 1000} km")
        print(f"Satellite elevation w.r.t to boresight = {param_mss.elevation} deg")
        print(f"Satellite azimuth w.r.t to boresight = {param_mss.azimuth} deg")
        print(f"Border distance = {d / 1000} km")
        print(f"NTN nadir offset w.r.t. IMT-TN cluster center = {ntn_footprint_x_offset[i] / 1000} km")
        print(f"Satellite elevation w.r.t. IMT-TN cluster center = {np.round(
            np.degrees(np.arctan(param_mss.altitude / ntn_footprint_x_offset[i])))} deg")
        slant_path_len = np.sqrt(param_mss.altitude**2 + ntn_footprint_x_offset[i]**2)
        print(f"Slant path lenght w.r.t. IMT-TN cluster center = {slant_path_len / 1000} km")
        fspl = np.round(20 * np.log10(slant_path_len) + 20 * np.log10(param_mss.frequency * 1e6) - 147.55, 2)
        print(f"Free space pathloss w.r.t. IMT-TN cluster center = {fspl}\n")

    ## Plot the coverage areas for NTN and TN
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ntn_topology.plot(ax, scale=1)
    macro_topology.plot(ax)

    ax.add_patch(macro_cicle)
    ax.add_patch(ntn_circle)
    ax.arrow(macro_topology_radius, 0, param_mss.x - ntn_footprint_radius - macro_topology_radius, 0,
             width=0.1, shape='full', color='red')
    ax.annotate(f'Border distance\n{border_distances_array[dist_idx] / 1000} km',
                ((macro_topology_radius + border_distances_array[dist_idx] / 2), 1000))

    plt.title("IMT-NTN vs IMT-TN Footprints")
    plt.xlabel("x-coordinate [m]")
    plt.ylabel("y-coordinate [m]")
    plt.tight_layout()
    plt.grid()

    plt.show()
    exit()
