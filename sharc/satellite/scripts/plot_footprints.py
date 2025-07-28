# Generates a 3D plot of the Earth with the satellites positions
# https://geopandas.org/en/stable/docs/user_guide/io.html
import os
import geopandas as gpd
import numpy as np
import plotly.graph_objects as go

from sharc.station_factory import StationFactory
from sharc.parameters.parameters_mss_d2d import ParametersOrbit, ParametersMssD2d
from sharc.support.sharc_geom import GeometryConverter
from sharc.satellite.utils.sat_utils import ecef2lla
from sharc.station_manager import StationManager
from sharc.parameters.antenna.parameters_antenna_s1528 import ParametersAntennaS1528


geoconv = GeometryConverter()

sys_lat = -14.5
sys_long = -45
sys_alt = 1200

geoconv.set_reference(
    sys_lat, sys_long, sys_alt
)


def plot_back(fig):
    """back half of sphere"""
    clor = 'rgb(220, 220, 220)'
    # Create a mesh grid for latitude and longitude.
    # For the "front" half, we can use longitudes from -180 to 0 degrees.
    # 50 latitude points from -90 to 90 degrees.
    lat_vals = np.linspace(-90, 90, 50)
    # 50 longitude points for the front half.
    lon_vals = np.linspace(0, 180, 50)
    # lon and lat will be 2D arrays.
    lon, lat = np.meshgrid(lon_vals, lat_vals)

    # Flatten the mesh to pass to the converter function.
    lat_flat = lat.flatten()
    lon_flat = lon.flatten()

    # Convert the lat/lon grid to transformed Cartesian coordinates.
    # Ensure your converter function can handle vectorized (numpy array) inputs.
    x_flat, y_flat, z_flat = geoconv.convert_lla_to_transformed_cartesian(
        lat_flat, lon_flat, 0)

    # Reshape the converted coordinates back to the 2D grid shape.
    x = x_flat.reshape(lat.shape)
    y = y_flat.reshape(lat.shape)
    z = z_flat.reshape(lat.shape)

    # Add the surface to the Plotly figure.
    fig.add_surface(
        x=x / 1e3, y=y / 1e3, z=z / 1e3,
        # Uniform color scale for a solid color.
        colorscale=[[0, clor], [1, clor]],
        opacity=1.0,
        showlegend=False,
        lighting=dict(diffuse=0.1),
        colorbar=dict(
            tickmode='array',
            tickvals=[0],
            ticktext=[""],
            title=""
        ),
    )


def plot_front(fig):
    """front half of sphere"""
    clor = 'rgb(220, 220, 220)'

    # Create a mesh grid for latitude and longitude.
    # For the "front" half, we can use longitudes from -180 to 0 degrees.
    # 50 latitude points from -90 to 90 degrees.
    lat_vals = np.linspace(-90, 90, 50)
    # 50 longitude points for the front half.
    lon_vals = np.linspace(-180, 0, 50)
    # lon and lat will be 2D arrays.
    lon, lat = np.meshgrid(lon_vals, lat_vals)

    # Flatten the mesh to pass to the converter function.
    lat_flat = lat.flatten()
    lon_flat = lon.flatten()

    # Convert the lat/lon grid to transformed Cartesian coordinates.
    # Ensure your converter function can handle vectorized (numpy array) inputs.
    x_flat, y_flat, z_flat = geoconv.convert_lla_to_transformed_cartesian(
        lat_flat, lon_flat, 0)

    # Reshape the converted coordinates back to the 2D grid shape.
    x = x_flat.reshape(lat.shape)
    y = y_flat.reshape(lat.shape)
    z = z_flat.reshape(lat.shape)

    # Add the surface to the Plotly figure.
    fig.add_surface(
        x=x / 1e3, y=y / 1e3, z=z / 1e3,
        # Uniform color scale for a solid color.
        colorscale=[[0, clor], [1, clor]],
        opacity=1.0,
        showlegend=False,
        lighting=dict(diffuse=0.1),
        colorbar=dict(
            tickmode='array',
            tickvals=[0],
            ticktext=[""],
            title=""
        ),
    )


def plot_polygon(poly, div=1, alt=0):
    xy_coords = poly.exterior.coords.xy
    lon = np.array(xy_coords[0])
    lat = np.array(xy_coords[1])

    # lon = lon * np.pi / 180
    # lat = lat * np.pi / 180

    # R = EARTH_RADIUS_KM
    x, y, z = geoconv.convert_lla_to_transformed_cartesian(lat, lon, alt)

    return x / div, y / div, z / div


def plot_mult_polygon(mult_poly, div=1e3):
    if mult_poly.geom_type == 'Polygon':
        return [plot_polygon(mult_poly, div=div, alt=1000)]
    elif mult_poly.geom_type == 'MultiPolygon':
        return [plot_polygon(poly, div=div, alt=1000) for poly in mult_poly.geoms]


def plot_globe_with_borders():
    # Read the shapefile.  Creates a DataFrame object
    countries_borders_shp_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                              "../../data/countries/ne_110m_admin_0_countries.shp")
    gdf = gpd.read_file(countries_borders_shp_file)
    fig = go.Figure()
    fig.update_layout(
        scene=dict(
            aspectmode="data",
            xaxis=dict(showbackground=False),
            yaxis=dict(showbackground=False),
            zaxis=dict(showbackground=False)
        ),
        margin=dict(l=0, r=0, b=0, t=0)
    )
    # return fig
    plot_front(fig)
    plot_back(fig)
    x_all, y_all, z_all = [], [], []

    for i in gdf.index:
        # print(gdf.loc[i].NAME)            # Call a specific attribute

        polys = gdf.loc[i].geometry         # Polygons or MultiPolygons

        if polys.geom_type == 'Polygon':
            x, y, z = plot_polygon(polys)
            x_all.extend(x / 1e3)
            x_all.extend([None])  # None separates different polygons
            y_all.extend(y / 1e3)
            y_all.extend([None])
            z_all.extend(z / 1e3)
            z_all.extend([None])

        elif polys.geom_type == 'MultiPolygon':

            for poly in polys.geoms:
                x, y, z = plot_polygon(poly)
                x_all.extend(x / 1e3)
                x_all.extend([None])  # None separates different polygons
                y_all.extend(y / 1e3)
                y_all.extend([None])
                z_all.extend(z / 1e3)
                z_all.extend([None])

    fig.add_trace(go.Scatter3d(x=x_all, y=y_all, z=z_all, mode='lines',
                               line=dict(color='rgb(0, 0, 0)'), showlegend=False))

    return fig


if __name__ == "__main__":
    colors = ['#fff7ec', '#fee8c8', '#fdd49e', '#fdbb84', '#fc8d59', '#7f0000']
    # steps from previous section
    # e.g. [3, 3, 14] means that 1st clr will be for -3 dB, 2nd for -6dB and 3rd for -20dB
    # from normalized gains
    step = [3, 17, 20, 20]  # dB
    # choose a seed to use for getting the satellites
    # if seed is different, different random numbers are taken
    SEED = 31
    # choose if colors appear continuously or in discrete steps
    DISCRETIZE = True
    # choose if a surface point received gain should be all sat gains summed
    # or if only max gain will be considered
    SUM_GAINS = False

    # Define the orbit parameters for two satellite constellations
    orbit_1 = ParametersOrbit(
        n_planes=28, sats_per_plane=120, phasing_deg=1.5, long_asc_deg=0,
        inclination_deg=53.0, perigee_alt_km=525.0, apogee_alt_km=525.0
    )

    # Antenna parameters
    g_max = 34.1  # dBi
    l_r = l_t = 1.6  # meters
    slr = 20  # dB
    n_side_lobes = 2  # number of side lobes
    freq = 2e3  # MHz

    antenna_params = ParametersAntennaS1528(
        antenna_gain=g_max,
        frequency=freq,  # in MHz
        bandwidth=5,  # in MHz
        slr=slr,
        n_side_lobes=n_side_lobes,
        l_r=l_r,
        l_t=l_t
    )

    spotbeam_radius = 39475  # meters

    # Configure the MSS D2D system parameters
    params = ParametersMssD2d(
        name="Example-MSS-D2D",
        antenna_pattern="ITU-R-S.1528-Taylor",
        num_sectors=19,
        antenna_s1528=antenna_params,
        intersite_distance=np.sqrt(3) * spotbeam_radius,
        orbits=[orbit_1]
    )
    params.sat_is_active_if.conditions = [
        # "MINIMUM_ELEVATION_FROM_ES",
        "LAT_LONG_INSIDE_COUNTRY",
    ]
    params.sat_is_active_if.minimum_elevation_from_es = 5.0
    params.sat_is_active_if.lat_long_inside_country.country_names = ["Brazil"]
    # params.beams_load_factor = 0.1
    # params.center_beam_positioning.type = "ANGLE_AND_DISTANCE_FROM_SUBSATELLITE"
    # params.center_beam_positioning.angle_from_subsatellite_phi.type = "~U(MIN,MAX)"
    # params.center_beam_positioning.angle_from_subsatellite_phi.distribution.min = -60.0
    # params.center_beam_positioning.angle_from_subsatellite_phi.distribution.max = 60.0
    # params.center_beam_positioning.distance_from_subsatellite.type = "~SQRT(U(0,1))*MAX"
    # params.center_beam_positioning.distance_from_subsatellite.distribution.max = 1111000.0

    print("instantiating stations")
    # Create a topology with a single base station
    from sharc.topology.topology_single_base_station import TopologySingleBaseStation
    imt_topology = TopologySingleBaseStation(
        cell_radius=500, num_clusters=1
    )

    # Create a random number generator
    rng = np.random.RandomState(seed=SEED)

    # Lists to store satellite positions (all and visible)
    # Plot the ground station (blue marker)
    # ground_sta_pos = lla2ecef(sys_lat, sys_long, sys_alt)
    ground_sta_pos = geoconv.convert_lla_to_transformed_cartesian(
        sys_lat, sys_long, sys_alt)

    center_of_earth = StationManager(1)
    # rotated and then translated center of earth
    center_of_earth.x = np.array([0.0])
    center_of_earth.y = np.array([0.0])
    center_of_earth.z = np.array([-geoconv.get_translation()])

    mss_d2d_manager = StationFactory.generate_mss_d2d(params, rng, geoconv)

    # Extract satellite positions
    x_vec = mss_d2d_manager.x / 1e3  # (Km)
    y_vec = mss_d2d_manager.y / 1e3  # (Km)
    z_vec = mss_d2d_manager.z / 1e3  # (Km)
    # Store all positions

    # Plot the globe with satellite positions
    fig = plot_globe_with_borders()

    # Set the camera position in Plotly
    show_range = 1e4
    fig.update_layout(
        scene=dict(
            zaxis=dict(
                range=(-show_range / 2, show_range / 2)
            ),
            yaxis=dict(
                range=(-show_range / 2, show_range / 2)
            ),
            xaxis=dict(
                range=(-show_range / 2, show_range / 2)
            ),
            camera=dict(
                center=dict(x=0, y=0, z=center_of_earth.z[0] / show_range / 1e3),
            )
        )
    )

    # # Satellite as fp center
    # center_fp_at_sat = 0
    # # get original sat xyz
    # orx, ory, orz = geoconv.revert_transformed_cartesian_to_cartesian(
    #     station_1.x[center_fp_at_sat],
    #     station_1.y[center_fp_at_sat],
    #     station_1.z[center_fp_at_sat],
    # )
    # sat_lat, sat_long, sat_alt = ecef2lla(orx, ory, orz)

    # lat_vals = np.linspace(sat_lat - 10.0, sat_lat + 10.0, 200)
    # lon_vals = np.linspace(sat_long - 10.0, sat_long + 10.0, 200)

    # # Ground station as fp center
    # lat_vals = np.linspace(geoconv.ref_lat - 10.0, geoconv.ref_lat + 10.0, 50)
    # lon_vals = np.linspace(geoconv.ref_long - 10.0, geoconv.ref_long + 10.0, 50)

    # Arbitrary range for fp calulation
    lat_vals = np.linspace(-33.69111, 2.81972, 200)
    lon_vals = np.linspace(-72.89583, -34.80861, 200)

    # lon and lat will be 2D arrays.
    lon, lat = np.meshgrid(lon_vals, lat_vals)

    # Flatten the mesh to pass to the converter function.
    lat_flat = lat.flatten()
    lon_flat = lon.flatten()

    # Convert the lat/lon grid to transformed Cartesian coordinates.
    # Ensure your converter function can handle vectorized (numpy array) inputs.
    x_flat, y_flat, z_flat = geoconv.convert_lla_to_transformed_cartesian(lat_flat, lon_flat, 0)

    # creates a StationManager to calculate the gains on
    surf_manager = StationManager(len(x_flat))
    surf_manager.x = x_flat
    surf_manager.y = y_flat
    surf_manager.z = z_flat
    surf_manager.height = z_flat

    station_1 = mss_d2d_manager
    mss_active = np.where(station_1.active)[0]
    station_2 = surf_manager
    station_2_active = np.where(station_2.active)[0]

    print("Calculating gains (memory intensive)")
    # Calculate vector and apointment off_axis
    gains = np.zeros((len(mss_active), len(station_2_active)))
    phi, theta = station_1.get_pointing_vector_to(station_2)
    off_axis_angle = station_1.get_off_axis_angle(station_2)
    phi, theta = station_1.get_pointing_vector_to(station_2)
    for k in mss_active:
        gains[k, station_2_active] = \
            station_1.antenna[k].calculate_gain(
                off_axis_angle_vec=off_axis_angle[k, station_2_active],
                theta_vec=theta[k, station_2_active],
                phi_vec=phi[k, station_2_active],
        )

    if SUM_GAINS:
        gains = 10 ** (gains / 10)
        gains = 10 * np.log10(np.sum(gains, axis=0))
    else:
        gains = np.max(gains, axis=0)

    world_surf_x = x_flat.reshape(lat.shape)
    world_surf_y = y_flat.reshape(lat.shape)
    world_surf_z = z_flat.reshape(lat.shape)
    reshaped_gain = gains.reshape(lat.shape)
    clor = 'rgb(220, 220, 220)'

    mx_gain = np.max(reshaped_gain)
    mn_gain = np.min(reshaped_gain)
    rnge = mx_gain - mn_gain
    n_steps = len(step) - 1
    colorscale = []
    bins = []
    at = 0
    offset = len(colors) - len(step)

    for i in range(0, n_steps + 1):
        ci = offset + i

        bins.append(mx_gain - at * rnge)
        at += step[i] / rnge

        if DISCRETIZE and i / n_steps - 1 / n_steps / 2 > 0:
            colorscale.append([i / n_steps - 1 / n_steps / 2, "rgb(100,100,100)"])
            colorscale.append([i / n_steps - 1 / n_steps / 2 + 0.001, colors[ci]])
        colorscale.append([i / n_steps, colors[ci]])
        if DISCRETIZE and i / n_steps + 1 / n_steps / 2 < 1.0:
            colorscale.append([i / n_steps + 1 / n_steps / 2 - 0.001, colors[ci]])
    # bins.append(mn_gain)
    bins.reverse()

    if DISCRETIZE:
        # cumul_steps = np.cumsum(steps)
        surfacecolor = np.digitize(reshaped_gain, bins, right=True)
        colorbar = dict(
            tickmode='array',
            tickvals=np.arange(0, len(bins)),
            ticktext=[f"< {bins[0]:.2f} dB"] + [f"{bins[i]:.2f} to {bins[i + 1]:.2f} dB" for i in range(len(bins) - 1)],
            title="Gain (dB)"
        )
    else:
        surfacecolor = reshaped_gain
        colorbar = None

    fig.add_surface(
        x=world_surf_x / 1e3, y=world_surf_y / 1e3, z=world_surf_z / 1e3,
        surfacecolor=surfacecolor,
        # Uniform color scale for a solid color.
        colorscale=colorscale,
        opacity=1.0,
        showlegend=False,
        # lighting=dict(diffuse=0)
        colorbar=colorbar,
    )

    # Plot all satellites (red markers)
    print("adding sats")
    fig.add_trace(go.Scatter3d(
        x=mss_d2d_manager.x / 1e3,
        y=mss_d2d_manager.y / 1e3,
        z=mss_d2d_manager.z / 1e3,
        mode='markers',
        marker=dict(size=2, color='red', opacity=0.5),
        showlegend=False
    ))

    # Plot visible satellites (green markers)
    # print(visible_positions['x'][visible_positions['x'] > 0])
    fig.add_trace(go.Scatter3d(
        x=mss_d2d_manager.x[mss_d2d_manager.active] / 1e3,
        y=mss_d2d_manager.y[mss_d2d_manager.active] / 1e3,
        z=mss_d2d_manager.z[mss_d2d_manager.active] / 1e3,
        mode='markers',
        marker=dict(size=3, color='green', opacity=0.8),
        showlegend=False
    ))

    fig.add_trace(go.Scatter3d(
        x=[ground_sta_pos[0] / 1e3],
        y=[ground_sta_pos[1] / 1e3],
        z=[ground_sta_pos[2] / 1e3],
        mode='markers',
        marker=dict(size=5, color='blue', opacity=1.0),
        showlegend=False
    ))

    polygons_lim = plot_mult_polygon(
        params.sat_is_active_if.lat_long_inside_country.filter_polygon
    )
    from functools import reduce

    lim_x, lim_y, lim_z = reduce(
        lambda acc, it: (list(it[0]) + [None] + acc[0], list(it[1]) + [None] + acc[1], list(it[2]) + [None] + acc[2]),
        polygons_lim,
        ([], [], [])
    )

    fig.add_trace(go.Scatter3d(
        x=lim_x,
        y=lim_y,
        z=lim_z,
        mode='lines',
        line=dict(color='rgb(0, 0, 255)'),
        showlegend=False
    ))

    # fig.add_trace(go.Scatter3d(
    #     x=center_of_earth.x / 1e3,
    #     y=center_of_earth.y / 1e3,
    #     z=center_of_earth.z / 1e3,
    #     mode='markers',
    #     marker=dict(size=5, color='black', opacity=1.0),
    #     showlegend=False
    # ))

    # Display the plot
    print(len(fig.data))
    fig.show()
