# Generates a 3D plot of the Earth with the satellites positions
# https://geopandas.org/en/stable/docs/user_guide/io.html
import os
import geopandas as gpd
import numpy as np
import plotly.graph_objects as go

from sharc.support.sharc_geom import GeometryConverter
from sharc.station_manager import StationManager
geoconv = GeometryConverter()

sys_lat = -14.5
sys_long = -45
sys_alt = 1200

geoconv.set_reference(
    sys_lat, sys_long, sys_alt
)

from sharc.satellite.ngso.orbit_model import OrbitModel
from sharc.satellite.utils.sat_utils import calc_elevation, lla2ecef
from sharc.satellite.ngso.constants import EARTH_RADIUS_KM
from sharc.parameters.parameters_mss_d2d import ParametersOrbit, ParametersMssD2d
from sharc.station_factory import StationFactory


def plot_back(fig):
    """back half of sphere"""
    clor = 'rgb(220, 220, 220)'
    # Create a mesh grid for latitude and longitude.
    # For the "front" half, we can use longitudes from -180 to 0 degrees.
    lat_vals = np.linspace(-90, 90, 50)       # 50 latitude points from -90 to 90 degrees.
    lon_vals = np.linspace(0, 180, 50)         # 50 longitude points for the front half.
    lon, lat = np.meshgrid(lon_vals, lat_vals)  # lon and lat will be 2D arrays.

    # Flatten the mesh to pass to the converter function.
    lat_flat = lat.flatten()
    lon_flat = lon.flatten()

    # Convert the lat/lon grid to transformed Cartesian coordinates.
    # Ensure your converter function can handle vectorized (numpy array) inputs.
    x_flat, y_flat, z_flat = geoconv.convert_lla_to_transformed_cartesian(lat_flat, lon_flat, 0)

    # Reshape the converted coordinates back to the 2D grid shape.
    x = x_flat.reshape(lat.shape)
    y = y_flat.reshape(lat.shape)
    z = z_flat.reshape(lat.shape)

    # Add the surface to the Plotly figure.
    fig.add_surface(
        x=(x / 1e3), y=(y / 1e3), z=(z / 1e3),
        colorscale=[[0, clor], [1, clor]],  # Uniform color scale for a solid color.
        opacity=1.0,
        showlegend=False,
        lighting=dict(diffuse=0.1)
    )


def plot_front(fig):
    """front half of sphere"""
    clor = 'rgb(220, 220, 220)'

    # Create a mesh grid for latitude and longitude.
    # For the "front" half, we can use longitudes from -180 to 0 degrees.
    lat_vals = np.linspace(-90, 90, 50)       # 50 latitude points from -90 to 90 degrees.
    lon_vals = np.linspace(-180, 0, 50)         # 50 longitude points for the front half.
    lon, lat = np.meshgrid(lon_vals, lat_vals)  # lon and lat will be 2D arrays.

    # Flatten the mesh to pass to the converter function.
    lat_flat = lat.flatten()
    lon_flat = lon.flatten()

    # Convert the lat/lon grid to transformed Cartesian coordinates.
    # Ensure your converter function can handle vectorized (numpy array) inputs.
    x_flat, y_flat, z_flat = geoconv.convert_lla_to_transformed_cartesian(lat_flat, lon_flat, 0)

    # Reshape the converted coordinates back to the 2D grid shape.
    x = x_flat.reshape(lat.shape)
    y = y_flat.reshape(lat.shape)
    z = z_flat.reshape(lat.shape)

    # Add the surface to the Plotly figure.
    fig.add_surface(
        x=(x / 1e3), y=(y / 1e3), z=(z / 1e3),
        colorscale=[[0, clor], [1, clor]],  # Uniform color scale for a solid color.
        opacity=1.0,
        showlegend=False,
        lighting=dict(diffuse=0.1)
    )


def plot_polygon(poly):

    xy_coords = poly.exterior.coords.xy
    lon = np.array(xy_coords[0])
    lat = np.array(xy_coords[1])

    # lon = lon * np.pi / 180
    # lat = lat * np.pi / 180

    # R = EARTH_RADIUS_KM
    x, y, z = geoconv.convert_lla_to_transformed_cartesian(lat, lon, 0)

    return x, y, z


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
    # Define the orbit parameters for two satellite constellations
    orbit_1 = ParametersOrbit(
        n_planes=6, sats_per_plane=8, phasing_deg=7.5, long_asc_deg=0,
        inclination_deg=52, perigee_alt_km=1414, apogee_alt_km=1414
    )

    orbit_2 = ParametersOrbit(
        n_planes=4, sats_per_plane=10, phasing_deg=5.0, long_asc_deg=90,
        inclination_deg=60, perigee_alt_km=1200, apogee_alt_km=1200
    )

    # Configure the MSS D2D system parameters
    params = ParametersMssD2d(
        name="Example-MSS-D2D",
        antenna_pattern="ITU-R-S.1528-Taylor",
        orbits=[orbit_1, orbit_2],
        frequency=2100,
        bandwidth=5
    )

    params.antenna_s1528.antenna_gain = 30.0

    params.sat_is_active_if.conditions = [
        "LAT_LONG_INSIDE_COUNTRY",
        "MINIMUM_ELEVATION_FROM_ES",
    ]
    params.sat_is_active_if.minimum_elevation_from_es = 5
    # params.sat_is_active_if.lat_long_inside_country.country_name = "Colombia"
    params.sat_is_active_if.lat_long_inside_country.margin_from_border = 0
    # params.sat_is_active_if.lat_long_inside_country.margin_from_border = 0
    params.sat_is_active_if.lat_long_inside_country.country_names = ["Brazil"]

    params.propagate_parameters()

    # Create a topology with a single base station
    from sharc.topology.topology_single_base_station import TopologySingleBaseStation
    imt_topology = TopologySingleBaseStation(
        cell_radius=500, num_clusters=1
    )

    # Create a random number generator
    rng = np.random.RandomState(seed=42)

    # Number of iterations (drops)
    NUM_DROPS = 100

    # Lists to store satellite positions (all and visible)
    all_positions = {'x': [], 'y': [], 'z': []}
    visible_positions = {'x': [], 'y': [], 'z': []}

    # Plot the ground station (blue marker)
    # ground_sta_pos = lla2ecef(sys_lat, sys_long, sys_alt)
    ground_sta_pos = geoconv.convert_lla_to_transformed_cartesian(sys_lat, sys_long, 1200.0)

    center_of_earth = StationManager(1)
    # rotated and then translated center of earth
    center_of_earth.x = np.array([0.0])
    center_of_earth.y = np.array([0.0])
    center_of_earth.z = np.array([-geoconv.get_translation()])

    vis_elevation = []
    for _ in range(NUM_DROPS):
        # Generate satellite positions using the StationFactory
        mss_d2d_manager = StationFactory.generate_mss_d2d(params, rng, geoconv)

        # Extract satellite positions
        x_vec = mss_d2d_manager.x / 1e3  # (Km)
        y_vec = mss_d2d_manager.y / 1e3  # (Km)
        z_vec = mss_d2d_manager.z / 1e3  # (Km)
        # Store all positions
        all_positions['x'].extend(x_vec)
        all_positions['y'].extend(y_vec)
        all_positions['z'].extend(z_vec)

        # Identify visible satellites
        vis_sat_idxs = np.where(mss_d2d_manager.active)[0]

        # should be pointing at nadir
        off_axis = mss_d2d_manager.get_off_axis_angle(center_of_earth)
        if len(np.where(off_axis > 0.01)[0]):
            print("AOPA, off axis parece estar errado")
            print("onde?: ", np.where(off_axis > 0.01))

        visible_positions['x'].extend(x_vec[vis_sat_idxs])
        visible_positions['y'].extend(y_vec[vis_sat_idxs])
        visible_positions['z'].extend(z_vec[vis_sat_idxs])
        vis_elevation.extend(mss_d2d_manager.elevation[vis_sat_idxs])

    # Flatten arrays
    all_positions['x'] = np.concatenate([all_positions['x']])
    all_positions['y'] = np.concatenate([all_positions['y']])
    all_positions['z'] = np.concatenate([all_positions['z']])

    visible_positions['x'] = np.concatenate([visible_positions['x']])
    visible_positions['y'] = np.concatenate([visible_positions['y']])
    visible_positions['z'] = np.concatenate([visible_positions['z']])

    # Plot the globe with satellite positions
    fig = plot_globe_with_borders()

    # Plot all satellites (red markers)
    print("adding sats")
    fig.add_trace(go.Scatter3d(
        x=all_positions['x'],
        y=all_positions['y'],
        z=all_positions['z'],
        mode='markers',
        marker=dict(size=2, color='red', opacity=0.5),
        showlegend=False
    ))

    # Plot visible satellites (green markers)
    # print(visible_positions['x'][visible_positions['x'] > 0])
    # print("vis_elevation", vis_elevation)
    fig.add_trace(go.Scatter3d(
        x=visible_positions['x'],
        y=visible_positions['y'],
        z=visible_positions['z'],
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
