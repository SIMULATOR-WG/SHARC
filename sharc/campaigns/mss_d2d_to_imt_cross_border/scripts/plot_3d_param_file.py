# Generates a 3D plot of the Earth with the satellites positions
# https://geopandas.org/en/stable/docs/user_guide/io.html
import geopandas as gpd
import numpy as np
import plotly.graph_objects as go
from pathlib import Path

from sharc.support.sharc_geom import GeometryConverter
from sharc.parameters.parameters import Parameters
from sharc.topology.topology_factory import TopologyFactory
from sharc.station_factory import StationFactory


def plot_back(fig, geoconv):
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
        x=x, y=y, z=z,
        colorscale=[[0, clor], [1, clor]],  # Uniform color scale for a solid color.
        opacity=1.0,
        showlegend=False,
        lighting=dict(diffuse=0.1)
    )


def plot_front(fig, geoconv):
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
        x=x, y=y, z=z,
        colorscale=[[0, clor], [1, clor]],  # Uniform color scale for a solid color.
        opacity=1.0,
        showlegend=False,
        lighting=dict(diffuse=0.1)
    )


def plot_polygon(poly, geoconv):

    xy_coords = poly.exterior.coords.xy
    lon = np.array(xy_coords[0])
    lat = np.array(xy_coords[1])

    # lon = lon * np.pi / 180
    # lat = lat * np.pi / 180

    # R = EARTH_RADIUS_KM
    x, y, z = geoconv.convert_lla_to_transformed_cartesian(lat, lon, 0)

    return x, y, z


def plot_mult_polygon(mult_poly, geoconv):
    if mult_poly.geom_type == 'Polygon':
        return [plot_polygon(mult_poly, geoconv)]
    elif mult_poly.geom_type == 'MultiPolygon':
        return [plot_polygon(poly, geoconv) for poly in mult_poly.geoms]


def plot_globe_with_borders(opaque_globe: bool, geoconv):
    # Read the shapefile.  Creates a DataFrame object
    project_root = Path(__file__).resolve().parents[4]
    countries_borders_shp_file = project_root / "sharc/data/countries/ne_110m_admin_0_countries.shp"
    gdf = gpd.read_file(countries_borders_shp_file)
    fig = go.Figure()
    # fig.update_layout(
    #     scene=dict(
    #         aspectmode="data",
    #         xaxis=dict(showbackground=False),
    #         yaxis=dict(showbackground=False),
    #         zaxis=dict(showbackground=False)
    #     ),
    #     margin=dict(l=0, r=0, b=0, t=0)
    # )
    if opaque_globe:
        plot_front(fig, geoconv)
        plot_back(fig, geoconv)
    x_all, y_all, z_all = [], [], []

    for i in gdf.index:
        # print(gdf.loc[i].NAME)            # Call a specific attribute

        polys = gdf.loc[i].geometry         # Polygons or MultiPolygons

        if polys.geom_type == 'Polygon':
            x, y, z = plot_polygon(polys, geoconv)
            x_all.extend(x)
            x_all.extend([None])  # None separates different polygons
            y_all.extend(y)
            y_all.extend([None])
            z_all.extend(z)
            z_all.extend([None])

        elif polys.geom_type == 'MultiPolygon':

            for poly in polys.geoms:
                x, y, z = plot_polygon(poly, geoconv)
                x_all.extend(x)
                x_all.extend([None])  # None separates different polygons
                y_all.extend(y)
                y_all.extend([None])
                z_all.extend(z)
                z_all.extend([None])

    fig.add_trace(go.Scatter3d(x=x_all, y=y_all, z=z_all, mode='lines',
                               line=dict(color='rgb(0, 0, 0)'), showlegend=False))

    return fig


if __name__ == "__main__":
    geoconv = GeometryConverter()
    SELECTED_SNAPSHOT_NUMBER = 0
    OPAQUE_GLOBE = True
    print(f"Plotting drop {SELECTED_SNAPSHOT_NUMBER}")
    # even when using the same drop number as in a simulation,
    # since the random generators are not in the same state, there isn't
    # a direct relationship between a drop in this plot and a drop in the
    # simulation loop
    # NOTE: if you want to plot the actual simulation scenarios for debugging,
    # you should do so inside the simulation loop
    print("  (not the same drop number as in simulation)")

    script_dir = Path(__file__).parent
    param_file = script_dir / "base_input.yaml"
    # param_file = script_dir / "../input/parameters_mss_d2d_to_imt_cross_border_0km_random_pointing_1beam_dl.yaml"
    param_file = param_file.resolve()
    print("File at:")
    print(f"  '{param_file}'")

    parameters = Parameters()
    parameters.set_file_name(param_file)
    parameters.read_params()

    geoconv.set_reference(
        parameters.imt.topology.central_latitude,
        parameters.imt.topology.central_longitude,
        parameters.imt.topology.central_altitude,
    )
    print(
        "imt at (lat, lon, alt) = ",
        (geoconv.ref_lat, geoconv.ref_long, geoconv.ref_alt),
    )

    import random
    random.seed(parameters.general.seed)

    secondary_seeds = [None] * parameters.general.num_snapshots

    max_seed = 2**32 - 1

    for index in range(parameters.general.num_snapshots):
        secondary_seeds[index] = random.randint(1, max_seed)

    seed = secondary_seeds[SELECTED_SNAPSHOT_NUMBER]

    topology = TopologyFactory.createTopology(parameters, geoconv)

    random_number_gen = np.random.RandomState(seed)

    # In case of hotspots, base stations coordinates have to be calculated
    # on every snapshot. Anyway, let topology decide whether to calculate
    # or not
    topology.calculate_coordinates(random_number_gen)

    # Create the base stations (remember that it takes into account the
    # network load factor)
    bs = StationFactory.generate_imt_base_stations(
        parameters.imt,
        parameters.imt.bs.antenna.array,
        topology, random_number_gen,
    )

    # Create the other system (FSS, HAPS, etc...)
    system = StationFactory.generate_system(
        parameters, topology, random_number_gen,
        geoconv
    )

    # Create IMT user equipments
    ue = StationFactory.generate_imt_ue(
        parameters.imt,
        parameters.imt.ue.antenna.array,
        topology, random_number_gen,
    )

    # Plot the globe with satellite positions
    fig = plot_globe_with_borders(OPAQUE_GLOBE, geoconv)

    polygons_lim = plot_mult_polygon(
        parameters.mss_d2d.sat_is_active_if.lat_long_inside_country.filter_polygon,
        geoconv
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

    # Plot all satellites (red markers)
    fig.add_trace(go.Scatter3d(
        x=system.x,
        y=system.y,
        z=system.z,
        mode='markers',
        marker=dict(size=2, color='red', opacity=0.5),
        showlegend=False
    ))

    # Plot visible satellites (green markers)
    # print(visible_positions['x'][visible_positions['x'] > 0])
    # print("vis_elevation", vis_elevation)
    fig.add_trace(go.Scatter3d(
        x=system.x[system.active],
        y=system.y[system.active],
        z=system.z[system.active],
        mode='markers',
        marker=dict(size=3, color='green', opacity=0.8),
        showlegend=False
    ))

    fig.add_trace(go.Scatter3d(
        x=ue.x,
        y=ue.y,
        z=ue.z,
        mode='markers',
        marker=dict(size=4, color='blue', opacity=1.0),
        showlegend=False
    ))

    fig.add_trace(go.Scatter3d(
        x=bs.x,
        y=bs.y,
        z=bs.z,
        mode='markers',
        marker=dict(size=4, color='black', opacity=1.0),
        showlegend=False
    ))

    # Display the plot
    range = 3e6
    fig.update_layout(
        scene=dict(
            zaxis=dict(
                range=(-range, range)
            ),
            yaxis=dict(
                range=(-range, range)
            ),
            xaxis=dict(
                range=(-range, range)
            ),
            camera=dict(
                center=dict(x=0, y=0, z=-geoconv.get_translation() / (2 * range)),  # Look at Earth's center
                # eye=eye,   # Camera position
                # center=dict(x=0, y=0, z=0),  # Look at Earth's center
                # up=dict(x=0, y=0, z=1)  # Ensure the up direction is correct
            )
        )
    )

    fig.show()
