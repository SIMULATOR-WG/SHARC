
"""
This module implements an IMT Mobile Satellite Service (MSS) for Direct Connectivity (D2D) topology.

It consists of a group of NGSO satellites that provide direct connectivity to user equipment (UE) on the ground.
The Space Stations positions are generated from the Keplerian elements of the orbits in the OrbitModel class.
Only a subset of Space Stations are used, which are the ones that are visible to the UE.
After satellite visibility is calculated, the ECEF coordinates are transformed to a new cartesian coordinate system
centered at the reference point defined in the GeometryConverter object.
The azimuth and elevation angles are also rotated to the new coordinate system.
The visible Space Stations are then used to generate the IMT Base Stations.
"""

import numpy as np
import geopandas as gpd
import shapely as shp
import pyproj

from sharc.topology.topology import Topology
from sharc.parameters.imt.parameters_imt_mss_dc import ParametersImtMssDc
from sharc.parameters.parameters_orbit import ParametersOrbit
from sharc.satellite.ngso.orbit_model import OrbitModel
from sharc.support.sharc_geom import GeometryConverter, rotate_angles_based_on_new_nadir
from sharc.topology.topology_ntn import TopologyNTN
from sharc.satellite.utils.sat_utils import calc_elevation
from sharc.support.sharc_geom import lla2ecef, cartesian_to_polar, polar_to_cartesian


class TopologyImtMssDc(Topology):
    def __init__(self, params: ParametersImtMssDc, geometry_converter: GeometryConverter):
        """Implements a IMT Mobile Satellite Service (MSS) for Direct Connectivity (D2D) topology.

        Parameters
        ----------
        params : ParametersImtMssDc
            Input parameters for the IMT MSS-DC topology
        geometry_converter : GeometryConverter
            GeometryConverter object that converts the ECEF coordintate system to one
            centered at GeometryConverter.reference.
        """
        # That means the we need to pass the groud reference points to the base stations generator
        self.is_space_station = True
        self.num_sectors = params.num_beams

        # Specific attributes
        self.geometry_converter = geometry_converter
        self.orbit_params = params
        self.space_station_x = None
        self.space_station_y = None
        self.space_station_z = None

        self.cell_radius = params.beam_radius
        # TODO: check this:
        self.intersite_distance = self.cell_radius * np.sqrt(3)

        self.lat = None
        self.lon = None
        self.orbits = []

        # Iterate through each orbit defined in the parameters
        for param in self.orbit_params.orbits:
            # Instantiate an OrbitModel for the current orbit
            orbit = OrbitModel(
                Nsp=param.sats_per_plane,  # Satellites per plane
                Np=param.n_planes,  # Number of orbital planes
                phasing=param.phasing_deg,  # Phasing angle in degrees
                long_asc=param.long_asc_deg,  # Longitude of ascending node in degrees
                omega=param.omega_deg,  # Argument of perigee in degrees
                delta=param.inclination_deg,  # Orbital inclination in degrees
                hp=param.perigee_alt_km,  # Perigee altitude in kilometers
                ha=param.apogee_alt_km,  # Apogee altitude in kilometers
                Mo=param.initial_mean_anomaly  # Initial mean anomaly in degrees
            )
            self.orbits.append(orbit)

    @staticmethod
    def get_coordinates(
        geometry_converter: GeometryConverter,
        orbit_params: ParametersImtMssDc,
        random_number_gen=np.random.RandomState(),
        only_active=True,
    ):
        """
        Computes the coordintates of the visible space stations
        """
        orbit_params.sat_is_active_if.validate("orbit_params.sat_is_active_if")
        # Calculate the total number of satellites across all orbits
        total_satellites = sum(orbit.n_planes * orbit.sats_per_plane for orbit in orbit_params.orbits)
        if any([
            not hasattr(orbit_params, attr)
                for attr in ["sat_is_active_if", "orbits", "beam_radius", "num_beams", "center_beam_positioning"]
            ]):
            raise ValueError(
                "Parameter passed to TopologyImtMssDc needs to contain all of the attributes:\n"
                '["sat_is_active_if", "orbits", "beam_radius", "num_beams", "center_beam_positioning"]'
            )

        idx_orbit = np.zeros(total_satellites, dtype=int)  # Add orbit index array

        # List to store indices of active satellites
        active_satellite_idxs = []

        MAX_ITER = 10000  # Maximum iterations to find at least one visible satellite
        i = 0  # Iteration counter for ensuring satellite visibility
        while len(active_satellite_idxs) == 0:
            # Initialize arrays to store satellite positions, angles and distance from center of earth
            all_positions = {"R": [], "lat": [], "lon": [], "sx": [], "sy": [], "sz": []}
            all_elevations = []  # Store satellite elevations
            all_azimuths = []  # Store satellite azimuths

            current_sat_idx = 0  # Index tracker for satellites across all orbits

            # Iterate through each orbit defined in the parameters
            for orbit_idx, param in enumerate(orbit_params.orbits):
                orbit = OrbitModel(
                    Nsp=param.sats_per_plane,  # Satellites per plane
                    Np=param.n_planes,  # Number of orbital planes
                    phasing=param.phasing_deg,  # Phasing angle in degrees
                    long_asc=param.long_asc_deg,  # Longitude of ascending node in degrees
                    omega=param.omega_deg,  # Argument of perigee in degrees
                    delta=param.inclination_deg,  # Orbital inclination in degrees
                    hp=param.perigee_alt_km,  # Perigee altitude in kilometers
                    ha=param.apogee_alt_km,  # Apogee altitude in kilometers
                    Mo=param.initial_mean_anomaly  # Initial mean anomaly in degrees
                )
                # Generate random positions for satellites in this orbit
                pos_vec = orbit.get_orbit_positions_random(rng=random_number_gen)

                # Determine the number of satellites in this orbit
                num_satellites = len(pos_vec["sx"])

                # Assign orbit index to satellites
                idx_orbit[current_sat_idx:current_sat_idx + num_satellites] = orbit_idx

                # Extract satellite positions and calculate distances
                sx, sy, sz = pos_vec['sx'], pos_vec['sy'], pos_vec['sz']
                r = np.sqrt(sx**2 + sy**2 + sz**2)  # Distance from Earth's center

                # When getting azimuth and elevation, we need to consider sx, sy and sz points
                # from the center of earth to the satellite, and we need to point the satellite
                # towards the center of earth
                elevations = np.degrees(np.arcsin(-sz / r))  # Calculate elevation angles
                azimuths = np.degrees(np.arctan2(-sy, -sx))  # Calculate azimuth angles

                # Append satellite positions and angles to global lists
                all_positions['lat'].extend(pos_vec['lat'])  # Latitudes
                all_positions['lon'].extend(pos_vec['lon'])  # Longitudes
                all_positions['sx'].extend(sx)  # X-coordinates
                all_positions['sy'].extend(sy)  # Y-coordinates
                all_positions['sz'].extend(sz)  # Z-coordinates
                all_positions["R"].extend(r)
                all_elevations.extend(elevations)  # Elevation angles
                all_azimuths.extend(azimuths)  # Azimuth angles

                active_sats_mask = np.ones(len(pos_vec['lat']), dtype=bool)

                if "MINIMUM_ELEVATION_FROM_ES" in orbit_params.sat_is_active_if.conditions:
                    # Calculate satellite visibility from base stations
                    elev_from_bs = calc_elevation(
                        geometry_converter.ref_lat,  # Latitude of base station
                        pos_vec['lat'],  # Latitude of satellites
                        geometry_converter.ref_long,  # Longitude of base station
                        pos_vec['lon'],  # Longitude of satellites
                        orbit.perigee_alt_km  # Perigee altitude in kilometers
                    )

                    # Determine visible satellites based on minimum elevation angle
                    active_sats_mask = active_sats_mask & (
                        elev_from_bs.flatten() >= orbit_params.sat_is_active_if.minimum_elevation_from_es
                    )

                if "MAXIMUM_ELEVATION_FROM_ES" in orbit_params.sat_is_active_if.conditions:
                    # no need to recalculate if already calculated above
                    if not "MINIMUM_ELEVATION_FROM_ES" in orbit_params.sat_is_active_if.conditions:
                        # Calculate satellite visibility from base stations
                        elev_from_bs = calc_elevation(
                            geometry_converter.ref_lat,  # Latitude of base station
                            pos_vec['lat'],  # Latitude of satellites
                            geometry_converter.ref_long,  # Longitude of base station
                            pos_vec['lon'],  # Longitude of satellites
                            orbit.perigee_alt_km  # Perigee altitude in kilometers
                        )

                    # Determine visible satellites based on minimum elevation angle
                    active_sats_mask = active_sats_mask & (
                        elev_from_bs.flatten() <= orbit_params.sat_is_active_if.maximum_elevation_from_es
                    )

                if "LAT_LONG_INSIDE_COUNTRY" in orbit_params.sat_is_active_if.conditions:
                    flat_active_lon = pos_vec["lon"].flatten()[active_sats_mask]
                    flat_active_lat = pos_vec["lat"].flatten()[active_sats_mask]

                    # create points(lon, lat) to compare to country
                    sats_points = gpd.points_from_xy(flat_active_lon, flat_active_lat, crs="EPSG:4326")

                    # Check if the satellite is inside the country polygon
                    polygon_mask = np.zeros_like(active_sats_mask)
                    polygon_mask[active_sats_mask] = sats_points.within(
                        orbit_params.sat_is_active_if.lat_long_inside_country.filter_polygon
                    )

                    active_sats_mask = active_sats_mask & polygon_mask

                visible_sat_idxs = np.arange(
                    current_sat_idx, current_sat_idx + len(pos_vec['lat']), dtype=int
                )[active_sats_mask]
                active_satellite_idxs.extend(visible_sat_idxs)

                # Update the index tracker for the next orbit
                current_sat_idx += len(sx)

            i += 1  # Increment iteration counter
            if i >= MAX_ITER:  # Check if maximum iterations reached
                raise RuntimeError(
                    "Maximum iterations reached, and no satellite was selected within the minimum elevation criteria."
                )
        # We have the list of visible satellites, now create a Topolgy of this subset and move the coordinate system
        # reference.
        if only_active:
            total_active_satellites = len(active_satellite_idxs)
            space_station_x = np.squeeze(np.array(all_positions['sx']))[active_satellite_idxs] * 1e3  # Convert X-coordinates to meters
            space_station_y = np.squeeze(np.array(all_positions['sy']))[active_satellite_idxs] * 1e3  # Convert Y-coordinates to meters
            space_station_z = np.squeeze(np.array(all_positions['sz']))[active_satellite_idxs] * 1e3  # Convert Z-coordinates to meters
            elevation = np.squeeze(np.array(all_elevations))[active_satellite_idxs]  # Elevation angles
            azimuth = np.squeeze(np.array(all_azimuths))[active_satellite_idxs]  # Azimuth angles
            # Store the latitude and longitude of the visible satellites for later use
            lat = np.squeeze(np.array(all_positions['lat']))[active_satellite_idxs]
            lon = np.squeeze(np.array(all_positions['lon']))[active_satellite_idxs]
        else:
            total_active_satellites = total_satellites
            space_station_x = np.squeeze(np.array(all_positions['sx'])) * 1e3  # Convert X-coordinates to meters
            space_station_y = np.squeeze(np.array(all_positions['sy'])) * 1e3  # Convert Y-coordinates to meters
            space_station_z = np.squeeze(np.array(all_positions['sz'])) * 1e3  # Convert Z-coordinates to meters
            elevation = np.squeeze(np.array(all_elevations))  # Elevation angles
            azimuth = np.squeeze(np.array(all_azimuths))  # Azimuth angles
            # Store the latitude and longitude of the visible satellites for later use
            lat = np.squeeze(np.array(all_positions['lat']))
            lon = np.squeeze(np.array(all_positions['lon']))

        rx, ry, rz = lla2ecef(
            np.squeeze(lat),
            np.squeeze(lon),
            0
        )
        earth_radius = np.sqrt(rx * rx + ry * ry + rz * rz)
        all_r = np.squeeze(np.array(all_positions['R'])) * 1e3
        if only_active:
            all_r = all_r[active_satellite_idxs]
        sat_altitude = np.array(all_r - earth_radius)

        # Convert the ECEF coordinates to the transformed cartesian coordinates and set the Space Station positions
        # used to generetate the IMT Base Stations
        space_station_x, space_station_y, space_station_z = \
            geometry_converter.convert_cartesian_to_transformed_cartesian(space_station_x, space_station_y, space_station_z)

        # Rotate the azimuth and elevation angles off the center beam the new transformed cartesian coordinates
        r = 1
        # transform pointing vectors, without considering geodesical earth coord system
        pointing_vec_x, pointing_vec_y, pointing_vec_z = polar_to_cartesian(r, azimuth, elevation)
        pointing_vec_x, pointing_vec_y, pointing_vec_z = \
            geometry_converter.convert_cartesian_to_transformed_cartesian(
                pointing_vec_x, pointing_vec_y, pointing_vec_z, translate=0)
        _, azimuth, elevation = cartesian_to_polar(pointing_vec_x, pointing_vec_y, pointing_vec_z)

        # We borrow the TopologyNTN method to calculate the sectors azimuth and elevation angles from their
        # respective x and y boresight coordinates
        sx, sy = TopologyNTN.get_sectors_xy(
            intersite_distance=orbit_params.beam_radius * np.sqrt(3),
            num_sectors=orbit_params.num_beams
        )

        assert (len(sx) == orbit_params.num_beams)
        assert (len(sy) == orbit_params.num_beams)

        # we give num_beams sectors to each satellite
        sx = np.resize(sx, orbit_params.num_beams * total_active_satellites)
        sy = np.resize(sy, orbit_params.num_beams * total_active_satellites)

        if orbit_params.center_beam_positioning.type == "ANGLE_AND_DISTANCE_FROM_SUBSATELLITE":
            match orbit_params.center_beam_positioning.angle_from_subsatellite_phi.type:
                case "FIXED":
                    azim_add = np.repeat(orbit_params.center_beam_positioning.angle_from_subsatellite_phi.fixed, total_active_satellites)
                case "~U(MIN,MAX)":
                    azim_add = random_number_gen.uniform(
                        orbit_params.center_beam_positioning.angle_from_subsatellite_phi.distribution.min,
                        orbit_params.center_beam_positioning.angle_from_subsatellite_phi.distribution.max,
                        total_active_satellites
                    )
                case "~SQRT(U(0,1))*MAX":
                    azim_add = random_number_gen.uniform(
                        0, 1,
                        total_active_satellites
                    ) * orbit_params.center_beam_positioning.angle_from_subsatellite_phi.distribution.max
                case _:
                    raise ValueError(
                        f"mss_d2d_params.center_beam_positioning.angle_from_subsatellite_phi.type = \n"
                        f"'{orbit_params.center_beam_positioning.angle_from_subsatellite_phi.type}' is not recognized!"
                    )

            match orbit_params.center_beam_positioning.distance_from_subsatellite.type:
                case "FIXED":
                    subsatellite_distance_add = np.repeat(orbit_params.center_beam_positioning.distance_from_subsatellite.fixed, total_active_satellites)
                case "~U(MIN,MAX)":
                    subsatellite_distance_add = random_number_gen.uniform(
                        orbit_params.center_beam_positioning.distance_from_subsatellite.distribution.min,
                        orbit_params.center_beam_positioning.distance_from_subsatellite.distribution.max,
                        total_active_satellites
                    )
                case "~SQRT(U(0,1))*MAX":
                    subsatellite_distance_add = random_number_gen.uniform(
                        0, 1,
                        total_active_satellites
                    ) * orbit_params.center_beam_positioning.distance_from_subsatellite.distribution.max
                case _:
                    raise ValueError(
                        f"mss_d2d_params.center_beam_positioning.distance_from_subsatellite.type = \n"
                        f"'{orbit_params.center_beam_positioning.angle_from_subsatellite_theta.type}' is not recognized!"
                    )

        elif orbit_params.center_beam_positioning.type == "ANGLE_FROM_SUBSATELLITE":
            match orbit_params.center_beam_positioning.angle_from_subsatellite_theta.type:
                case "FIXED":
                    off_nadir_add = np.repeat(orbit_params.center_beam_positioning.angle_from_subsatellite_theta.fixed, total_active_satellites)
                case "~U(MIN,MAX)":
                    off_nadir_add = random_number_gen.uniform(
                        orbit_params.center_beam_positioning.angle_from_subsatellite_theta.distribution.min,
                        orbit_params.center_beam_positioning.angle_from_subsatellite_theta.distribution.max,
                        total_active_satellites
                    )
                case "~SQRT(U(0,1))*MAX":
                    off_nadir_add = random_number_gen.uniform(
                        0, 1,
                        total_active_satellites
                    ) * orbit_params.center_beam_positioning.angle_from_subsatellite_theta.distribution.max
                case _:
                    raise ValueError(
                        f"mss_d2d_params.center_beam_positioning.angle_from_subsatellite_theta.type = \n"
                        f"'{orbit_params.center_beam_positioning.angle_from_subsatellite_theta.type}' is not recognized!"
                    )
            subsatellite_distance_add = sat_altitude * np.tan(off_nadir_add)

            match orbit_params.center_beam_positioning.angle_from_subsatellite_phi.type:
                case "FIXED":
                    azim_add = np.repeat(orbit_params.center_beam_positioning.angle_from_subsatellite_phi.fixed, total_active_satellites)
                case "~U(MIN,MAX)":
                    azim_add = random_number_gen.uniform(
                        orbit_params.center_beam_positioning.angle_from_subsatellite_phi.distribution.min,
                        orbit_params.center_beam_positioning.angle_from_subsatellite_phi.distribution.max,
                        total_active_satellites
                    )
                case "~SQRT(U(0,1))*MAX":
                    azim_add = random_number_gen.uniform(
                        0, 1,
                        total_active_satellites
                    ) * orbit_params.center_beam_positioning.angle_from_subsatellite_phi.distribution.max
                case _:
                    raise ValueError(
                        f"mss_d2d_params.center_beam_positioning.angle_from_subsatellite_phi.type = \n"
                        f"'{orbit_params.center_beam_positioning.angle_from_subsatellite_phi.type}' is not recognized!"
                    )

        subsatellite_distance_add = np.repeat(subsatellite_distance_add, orbit_params.num_beams)
        azim_add = np.repeat(azim_add, orbit_params.num_beams)

        sx += subsatellite_distance_add

        # Calculate the azimuth and elevation angles for each beam
        # as though their nadir is at (0,0)
        # before rotating them
        beams_azim = np.rad2deg(np.arctan2(sy, sx)) + azim_add
        beams_elev = np.rad2deg(np.arctan2(np.sqrt(sy * sy + sx * sx),
                                           np.repeat(sat_altitude, orbit_params.num_beams))
                            ) - 90

        beams_azim = beams_azim.reshape(
            (total_active_satellites, orbit_params.num_beams)
        )

        beams_elev = beams_elev.reshape(
            (total_active_satellites, orbit_params.num_beams)
        )

        # Rotate and set the each beam azimuth and elevation angles - only for the visible satellites
        for i in range(total_active_satellites):
            # Rotate the azimuth and elevation angles based on the new nadir point
            beams_elev[i], beams_azim[i] = rotate_angles_based_on_new_nadir(
                beams_elev[i],
                beams_azim[i],
                elevation[i],
                azimuth[i]
            )

        # In SHARC each sector is treated as a separate base station, so we need to repeat the satellite positions
        # for each sector.
        space_station_x = np.repeat(space_station_x, orbit_params.num_beams)
        space_station_y = np.repeat(space_station_y, orbit_params.num_beams)
        space_station_z = np.repeat(space_station_z, orbit_params.num_beams)

        num_base_stations = orbit_params.num_beams * total_active_satellites
        elevation = beams_elev.flatten()
        azimuth = beams_azim.flatten()
        lat = np.repeat(lat, orbit_params.num_beams)
        lon = np.repeat(lon, orbit_params.num_beams)

        altitudes = np.repeat(sat_altitude, orbit_params.num_beams)

        assert (space_station_x.shape == (num_base_stations,))
        assert (space_station_y.shape == (num_base_stations,))
        assert (space_station_z.shape == (num_base_stations,))
        assert (lat.shape == (num_base_stations,))
        assert (lon.shape == (num_base_stations,))
        assert (altitudes.shape == (num_base_stations,))
        assert (elevation.shape == (num_base_stations,))
        assert (azimuth.shape == (num_base_stations,))
        assert (sx.shape == (num_base_stations,))
        assert (sy.shape == (num_base_stations,))

        # update indices (multiply by num_beams)
        # and make all num_beams of satellite active
        active_satellite_idxs = np.ravel(np.array(active_satellite_idxs)[:, np.newaxis] * orbit_params.num_beams +
                                         np.arange(orbit_params.num_beams))

        return {
            "num_satellites": num_base_stations,
            "num_active_satellites": len(active_satellite_idxs),
            "active_satellites_idxs": active_satellite_idxs,
            "sat_x": space_station_x,
            "sat_y": space_station_y,
            "sat_z": space_station_z,
            "sat_lat": lat,
            "sat_lon": lon,
            "sat_alt": altitudes,
            "sat_antenna_elev": elevation,
            "sat_antenna_azim": azimuth,
            "sectors_x": sx,
            "sectors_y": sy,
            "sectors_z": np.zeros_like(sx)
        }

    def calculate_coordinates(self, random_number_gen=np.random.RandomState()):
        """
        Computes the coordintates of the visible space stations
        """
        self.geometry_converter.validate()

        sat_values = self.get_coordinates(self.geometry_converter, self.orbit_params, random_number_gen)

        self.num_base_stations = sat_values["num_satellites"]

        self.space_station_x = sat_values["sat_x"]
        self.space_station_y = sat_values["sat_y"]
        self.space_station_z = sat_values["sat_z"]
        self.height = sat_values["sat_alt"]
        self.lat = sat_values["sat_lat"]
        self.lon = sat_values["sat_lon"]

        self.elevation = sat_values["sat_antenna_elev"]
        self.azimuth = sat_values["sat_antenna_azim"]

        self.x = sat_values["sectors_x"]
        self.y = sat_values["sectors_y"]
        self.z = sat_values["sectors_z"]

        self.indoor = np.zeros(self.num_base_stations, dtype=bool)  # ofcourse, all are outdoor

        return

    # We can factor this out if another topology also ends up needing this
    def transform_ue_xyz(self, bs_i, x, y, z):
        x, y, z = super().transform_ue_xyz(bs_i, x, y, z)

        # translate by earth radius on the lat long passed
        # this way we mantain the center of topology on surface of earth
        # considering a geodesic earth
        # since we expect the area to be small, we can just consider
        # the center of the topology for this translation
        rx, ry, rz = lla2ecef(self.lat[bs_i], self.lon[bs_i], 0)
        earth_radius_at_sat_nadir = np.sqrt(rx * rx + ry * ry + rz * rz)
        z += earth_radius_at_sat_nadir

        # get angle around y axis
        around_y = np.arctan2(
                np.sqrt(
                    self.space_station_x[bs_i] ** 2 +
                    self.space_station_y[bs_i] ** 2
                ),
                self.space_station_z[bs_i] + self.geometry_converter.get_translation()
            )

        # get around z axis
        around_z = np.arctan2(self.space_station_y[bs_i], self.space_station_x[bs_i])

        # rotating around y
        nx = x * np.cos(around_y) + z * np.sin(around_y)
        nz = x * -np.sin(around_y) + z * np.cos(around_y)
        x = nx
        z = nz

        # now we have (x,y) = (A, 0)
        # rotating around z
        nx = x * np.cos(around_z) + y * -np.sin(around_z)
        ny = x * np.sin(around_z) + y * np.cos(around_z)
        x = nx
        y = ny

        # translate ue back so other system is in (0,0,0)
        z -= self.geometry_converter.get_translation()

        return (x, y, z)


# Example usage
if __name__ == '__main__':
    from sharc.parameters.imt.parameters_imt_mss_dc import ParametersImtMssDc
    from sharc.support.sharc_geom import GeometryConverter

    # Define the parameters for the IMT MSS-DC topology
    # SystemA Orbit parameters
    orbit = ParametersOrbit(
        n_planes=28,
        sats_per_plane=120,
        phasing_deg=1.5,
        long_asc_deg=0.0,
        inclination_deg=53.0,
        perigee_alt_km=525,
        apogee_alt_km=525
    )
    params = ParametersImtMssDc(
        beam_radius=39745.0,
        num_beams=19,
        orbits=[orbit]
    )
    params.sat_is_active_if.conditions = [
        "LAT_LONG_INSIDE_COUNTRY",
        "MINIMUM_ELEVATION_FROM_ES",
    ]
    params.sat_is_active_if.minimum_elevation_from_es = 5
    params.sat_is_active_if.lat_long_inside_country.country_names = ["Brazil"]

    # Define the geometry converter
    geometry_converter = GeometryConverter()
    geometry_converter.set_reference(-15.0, -42.0, 1200)

    # Instantiate the IMT MSS-DC topology
    imt_mss_dc_topology = TopologyImtMssDc(params, geometry_converter)

    # Calculate the coordinates of the space stations
    rng = np.random.RandomState(101)
    imt_mss_dc_topology.calculate_coordinates(random_number_gen=rng)

    # Plot the IMT MSS-DC space stations after selecting the visible ones and transforming the coordinate system
    import plotly.graph_objects as go

    # Create a 3D scatter plot using Plotly
    fig = go.Figure(data=[go.Scatter3d(
        x=imt_mss_dc_topology.space_station_x / 1e3,
        y=imt_mss_dc_topology.space_station_y / 1e3,
        z=imt_mss_dc_topology.space_station_z / 1e3,
        mode='markers',
        marker=dict(
            size=4,
            color='red',
            opacity=0.8
        )
    )])

    # Set plot title and axis labels
    fig.update_layout(
        title='IMT MSS-DC Topology',
        scene=dict(
            xaxis_title='X [m]',
            yaxis_title='Y [m]',
            zaxis_title='Z [m]'
        )
    )

    # Add a point at the origin
    fig.add_trace(go.Scatter3d(
        x=[0],
        y=[0],
        z=[0],
        mode='markers',
        marker=dict(
            size=6,
            color='blue',
            opacity=1.0
        ),
        name='Reference'
    ))
    # Calculate the elevation with respect to the x-y plane
    elevation_xy_plane = np.degrees(
        np.arctan2(
            imt_mss_dc_topology.space_station_z,
            np.sqrt(imt_mss_dc_topology.space_station_x**2 + imt_mss_dc_topology.space_station_y**2)
        )
    )

    # Add the elevation with respect to the x-y plane to the plot
    fig.add_trace(go.Scatter3d(
        x=imt_mss_dc_topology.space_station_x / 1e3,
        y=imt_mss_dc_topology.space_station_y / 1e3,
        z=imt_mss_dc_topology.space_station_z / 1e3,
        mode='markers',
        marker=dict(
            size=4,
            color=elevation_xy_plane,
            colorscale='Viridis',
            colorbar=dict(title='Elevation (degrees)', x=-0.1),
            opacity=0.8
        ),
        name='Space Stations Positions'
    ))

    # Add lines between the origin and the IMT space stations
    for x, y, z in zip(imt_mss_dc_topology.space_station_x / 1e3, imt_mss_dc_topology.space_station_y / 1e3,
                       imt_mss_dc_topology.space_station_z / 1e3):
        fig.add_trace(go.Scatter3d(
            x=[0, x],
            y=[0, y],
            z=[0, z],
            mode='lines',
            line=dict(color='green', width=2, dash='dash'),
            name='Elevation Line'
        ))
    # Suppress the legend for the elevation plot
    fig.update_traces(showlegend=False, selector=dict(name='Elevation Line'))

    # Plot beam boresight vectors
    boresight_length = 100  # Length of the boresight vectors for visualization
    boresight_x, boresight_y, boresight_z = polar_to_cartesian(
        boresight_length,
        imt_mss_dc_topology.azimuth,
        imt_mss_dc_topology.elevation
    )
    # Add arrow heads to the end of the boresight vectors
    for x, y, z, bx, by, bz in zip(imt_mss_dc_topology.space_station_x / 1e3,
                                   imt_mss_dc_topology.space_station_y / 1e3,
                                   imt_mss_dc_topology.space_station_z / 1e3,
                                   boresight_x,
                                   boresight_y,
                                   boresight_z):
        fig.add_trace(go.Cone(
            x=[x + bx],
            y=[y + by],
            z=[z + bz],
            u=[bx],
            v=[by],
            w=[bz],
            colorscale=[[0, 'orange'], [1, 'orange']],
            sizemode='absolute',
            sizeref=40,
            showscale=False
        ))
    for x, y, z, bx, by, bz in zip(imt_mss_dc_topology.space_station_x / 1e3,
                                   imt_mss_dc_topology.space_station_y / 1e3,
                                   imt_mss_dc_topology.space_station_z / 1e3,
                                   boresight_x,
                                   boresight_y,
                                   boresight_z):
        fig.add_trace(go.Scatter3d(
            x=[x, x + bx],
            y=[y, y + by],
            z=[z, z + bz],
            mode='lines',
            line=dict(color='orange', width=2),
            name='Boresight'
        ))
    # Suppress the legend for the boresight plot
    fig.update_traces(showlegend=False, selector=dict(name='Boresight'))

    # Maintain axis proportions
    fig.update_layout(scene_aspectmode='data')

    fig.show()

    # Plot the interception of the boresight vectors with the x-y plane
    fig_intercept = go.Figure()

    # Calculate the interception points of the boresight vectors with the x-y plane
    t_intercept = -imt_mss_dc_topology.space_station_z / boresight_z
    intercept_x = imt_mss_dc_topology.space_station_x + t_intercept * boresight_x
    intercept_y = imt_mss_dc_topology.space_station_y + t_intercept * boresight_y

    # Add the interception points to the plot
    fig_intercept.add_trace(go.Scatter(
        x=intercept_x / 1e3,
        y=intercept_y / 1e3,
        mode='markers',
        marker=dict(
            size=6,
            color='purple',
            opacity=0.8
        ),
        name='Interception Points'
    ))

    # Add the space station positions for reference
    fig_intercept.add_trace(go.Scatter(
        x=imt_mss_dc_topology.space_station_x / 1e3,
        y=imt_mss_dc_topology.space_station_y / 1e3,
        mode='markers',
        marker=dict(
            size=4,
            color='red',
            opacity=0.8
        ),
        name='Space Stations'
    ))

    # Set plot title and axis labels
    fig_intercept.update_layout(
        title='Interception of Boresight Vectors with x-y Plane',
        xaxis_title='X [km]',
        yaxis_title='Y [km]',
        showlegend=True
    )

    # Maintain axis proportions
    fig_intercept.update_yaxes(scaleanchor="x", scaleratio=1)

    # Show the plot
    fig_intercept.show()

    # Plot the IMT MSS-DC space stations in a 2D plane
    fig_2d = go.Figure()

    # Add circles centered at the (x, y) coordinates of the space stations
    for x, y in zip(imt_mss_dc_topology.x, imt_mss_dc_topology.y):
        circle = go.Scatter(
            x=[x / 1e3 + imt_mss_dc_topology.orbit_params.beam_radius /
               1e3 * np.cos(theta) for theta in np.linspace(0, 2 * np.pi, 100)],
            y=[y / 1e3 + imt_mss_dc_topology.orbit_params.beam_radius /
               1e3 * np.sin(theta) for theta in np.linspace(0, 2 * np.pi, 100)],
            mode='lines',
            line=dict(color='blue')
        )
        fig_2d.add_trace(circle)

    # Set plot title and axis labels
    fig_2d.update_layout(
        title='IMT MSS-DC Topology in x-y Plane',
        xaxis_title='X [m]',
        yaxis_title='Y [m]',
        showlegend=False
    )

    # Maintain axis proportions
    fig_2d.update_yaxes(scaleanchor="x", scaleratio=1)

    # Show the plot
    fig_2d.show()

    # Print the elevation angles
    print('Elevation angles:', imt_mss_dc_topology.elevation)
    # Print the azimuth angles
    print('Azimuth angles:', imt_mss_dc_topology.azimuth)
    # Print the elevation w.r.t. the x-y plane
    print('Elevation w.r.t. XY plane:', elevation_xy_plane)
    # Print the slant range
    idxs = np.arange(imt_mss_dc_topology.num_base_stations // imt_mss_dc_topology.num_sectors) * \
        imt_mss_dc_topology.num_sectors
    print('Slant range:', np.sqrt(imt_mss_dc_topology.space_station_x[idxs]**2 +
                                  imt_mss_dc_topology.space_station_y[idxs]**2 +
                                  imt_mss_dc_topology.space_station_z[idxs]**2) / 1e3)
