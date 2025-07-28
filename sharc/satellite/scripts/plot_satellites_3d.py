# Generates a 3D plot of the Earth with the satellites positions
# https://geopandas.org/en/stable/docs/user_guide/io.html
import os
import geopandas as gpd
import numpy as np
import plotly.graph_objects as go

from sharc.satellite.ngso.orbit_model import OrbitModel
from sharc.satellite.utils.sat_utils import calc_elevation, lla2ecef
from sharc.satellite.ngso.constants import EARTH_RADIUS_KM


def plot_back(fig):
    """back half of sphere"""
    clor = 'rgb(220, 220, 220)'
    R = np.sqrt(EARTH_RADIUS_KM)
    u_angle = np.linspace(0, np.pi, 25)
    v_angle = np.linspace(0, np.pi, 25)
    x_dir = np.outer(R * np.cos(u_angle), R * np.sin(v_angle))
    y_dir = np.outer(R * np.sin(u_angle), R * np.sin(v_angle))
    z_dir = np.outer(R * np.ones(u_angle.shape[0]), R * np.cos(v_angle))
    fig.add_surface(z=z_dir, x=x_dir, y=y_dir, colorscale=[[0, clor], [1, clor]],
                    opacity=1.0, showlegend=False, lighting=dict(
                    # opacity=fig.sphere_alpha, colorscale=[[0, fig.sphere_color], [1, fig.sphere_color]])
                    diffuse=0.1))


def plot_front(fig):
    """front half of sphere"""
    clor = 'rgb(220, 220, 220)'
    R = np.sqrt(EARTH_RADIUS_KM)
    u_angle = np.linspace(-np.pi, 0, 25)
    v_angle = np.linspace(0, np.pi, 25)
    x_dir = np.outer(R * np.cos(u_angle), R * np.sin(v_angle))
    y_dir = np.outer(R * np.sin(u_angle), R * np.sin(v_angle))
    z_dir = np.outer(R * np.ones(u_angle.shape[0]), R * np.cos(v_angle))
    fig.add_surface(z=z_dir, x=x_dir, y=y_dir, colorscale=[[0, clor], [1, clor]], opacity=1.0, showlegend=False,
                    lighting=dict(
                        # opacity=fig.sphere_alpha, colorscale=[[0, fig.sphere_color], [1, fig.sphere_color]])
                        diffuse=0.1))


def plot_polygon(poly):

    xy_coords = poly.exterior.coords.xy
    lon = np.array(xy_coords[0])
    lat = np.array(xy_coords[1])

    lon = lon * np.pi / 180
    lat = lat * np.pi / 180

    R = EARTH_RADIUS_KM
    x = R * np.cos(lat) * np.cos(lon)
    y = R * np.cos(lat) * np.sin(lon)
    z = R * np.sin(lat)

    return x, y, z


def plot_globe_with_borders():
    # Read the shapefile.  Creates a DataFrame object
    countries_borders_shp_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                              "../../data/countries/ne_110m_admin_0_countries.shp")
    gdf = gpd.read_file(countries_borders_shp_file)
    fig = go.Figure()
    plot_front(fig)
    plot_back(fig)

    for i in gdf.index:
        # print(gdf.loc[i].NAME)            # Call a specific attribute

        polys = gdf.loc[i].geometry         # Polygons or MultiPolygons

        if polys.geom_type == 'Polygon':
            x, y, z = plot_polygon(polys)
            fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines',
                                       line=dict(color='rgb(0, 0,0)'), showlegend=False))

        elif polys.geom_type == 'MultiPolygon':

            for poly in polys.geoms:
                x, y, z = plot_polygon(poly)
                fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines',
                                           line=dict(color='rgb(0, 0,0)'), showlegend=False))
    return fig


if __name__ == "__main__":

    # Plot Global Star orbit using OrbitModel object
    orbit = OrbitModel(
        Nsp=6,  # number of sats per plane
        Np=8,  # number of planes
        phasing=7.5,  # phasing in degrees
        long_asc=0,  # longitude of the ascending node in degrees
        omega=0,  # argument of perigee in degrees
        delta=52,  # inclination in degrees
        hp=1414,  # perigee altitude in km
        ha=1414,  # apogee altitude in km
        Mo=0  # mean anomaly in degrees
    )

    # Set parameters for ground station
    GROUND_STA_LAT = -15.7801
    GROUND_STA_LON = -42.9292
    MIN_ELEV_ANGLE_DEG = 5.0

    # Plot satellite traces from a time interval
    fig = plot_globe_with_borders()
    pos_vec = orbit.get_satellite_positions_time_interval(initial_time_secs=0, interval_secs=5, n_periods=1)
    fig.add_trace(go.Scatter3d(x=pos_vec['sx'].flatten(),
                               y=pos_vec['sy'].flatten(),
                               z=pos_vec['sz'].flatten(),
                               mode='lines',
                               showlegend=False))

    fig.update_layout(
        title_text='Satellite Traces',
        scene=dict(
            xaxis_title='X [km]',
            yaxis_title='Y [km]',
            zaxis_title='Z [km]',
            aspectmode='data'
        ),
        margin=dict(l=0, r=0, b=0, t=0)
    )

    fig.show()

    # Plot satellites positions taken randomly
    fig = plot_globe_with_borders()
    NUM_DROPS = 100
    rng = np.random.RandomState(seed=6)
    pos_vec = orbit.get_orbit_positions_random(rng=rng, n_samples=NUM_DROPS)
    fig.add_trace(go.Scatter3d(x=pos_vec['sx'].flatten(),
                               y=pos_vec['sy'].flatten(),
                               z=pos_vec['sz'].flatten(),
                               mode='markers',
                               marker=dict(size=2,
                                           color='red',
                                           opacity=0.8),
                               showlegend=False))
    fig.show()

    # Show visible satellites from ground-station
    fig = plot_globe_with_borders()
    NUM_DROPS = 1000
    rng = np.random.RandomState(seed=6)
    pos_vec = orbit.get_orbit_positions_random(rng=rng, n_samples=NUM_DROPS)
    look_angles = calc_elevation(GROUND_STA_LAT, pos_vec['lat'], GROUND_STA_LON, pos_vec['lon'], orbit.apogee_alt_km)
    elevation_angles_per_drop = look_angles[np.where(np.array(look_angles) > 0)]
    num_of_visible_sats_per_drop = np.sum(look_angles > MIN_ELEV_ANGLE_DEG, axis=0)

    # plot all satellites in drops
    fig.add_trace(go.Scatter3d(x=pos_vec['sx'].flatten(),
                               y=pos_vec['sy'].flatten(),
                               z=pos_vec['sz'].flatten(),
                               mode='markers',
                               marker=dict(size=2,
                                           color='red',
                                           opacity=0.8),
                               showlegend=False))

    # plot visible satellites
    fig.add_trace(go.Scatter3d(
        x=pos_vec['sx'][np.where(look_angles > MIN_ELEV_ANGLE_DEG)].flatten(),
        y=pos_vec['sy'][np.where(look_angles > MIN_ELEV_ANGLE_DEG)].flatten(),
        z=pos_vec['sz'][np.where(look_angles > MIN_ELEV_ANGLE_DEG)].flatten(),
        mode='markers',
        marker=dict(size=2,
                    color='green',
                    opacity=0.8),
        showlegend=False))

    # plot ground station
    groud_sta_pos = lla2ecef(GROUND_STA_LAT, GROUND_STA_LON, 0.0)
    fig.add_trace(go.Scatter3d(
        x=np.array(groud_sta_pos[0] / 1e3),
        y=np.array(groud_sta_pos[1] / 1e3),
        z=np.array(groud_sta_pos[2] / 1e3),
        mode='markers',
        marker=dict(size=4,
                    color='blue',
                    opacity=1.0),
        showlegend=False))

    fig.show()
