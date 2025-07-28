import numpy as np
import shapely as shp
import typing
import pyproj

from sharc.satellite.utils.sat_utils import lla2ecef, ecef2lla
from sharc.station_manager import StationManager


def cartesian_to_polar(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple:
    """
    Converts cartesian coordinates to polar coordinates.

    Parameters
    ----------
    x : np.ndarray
        x coordinate in meters
    y : np.ndarray
        y coordinate in meters
    z : np.ndarray
        z coordinate in meters

    Returns
    -------
    tuple
        range, azimuth and elevation in meters, degrees and degrees
    """
    # range calculation
    r = np.sqrt(x**2 + y**2 + z**2)

    # azimuth calculation
    azimuth = np.arctan2(y, x)

    # elevation calculation
    elevation = np.arcsin(z / r)

    return r, np.degrees(azimuth), np.degrees(elevation)


def polar_to_cartesian(r: np.ndarray, azimuth: np.ndarray, elevation: np.ndarray) -> tuple:
    """
    Converts polar coordinates to cartesian coordinates.

    Parameters
    ----------
    r : np.ndarray
        range in meters
    azimuth : np.ndarray
        azimuth in degrees
    elevation : np.ndarray
        elevation in degrees

    Returns
    -------
    tuple
        x, y and z coordinates in meters
    """
    azimuth = np.radians(azimuth)
    elevation = np.radians(elevation)
    x = r * np.cos(elevation) * np.cos(azimuth)
    y = r * np.cos(elevation) * np.sin(azimuth)
    z = r * np.sin(elevation)

    return x, y, z


def get_rotation_matrix(around_z, around_y):
    """
    Rotates with the right hand rule around the z axis (similar to simulator azimuth)
    and with the right hand rule around the y axis (similar to simulator elevation)
    """
    alpha = np.deg2rad(around_z)
    beta = np.deg2rad(around_y)

    ry = np.matrix([
        [np.cos(beta), 0.0, np.sin(beta)],
        [0.0, 1.0, 0.0],
        [-np.sin(beta), 0.0, np.cos(beta)],
    ])
    rz = np.matrix([
        [np.cos(alpha), -np.sin(alpha), 0.0],
        [np.sin(alpha), np.cos(alpha), 0.0],
        [0.0, 0.0, 1.0],
    ])

    return rz * ry


def rotate_angles_based_on_new_nadir(elev, azim, nadir_elev, nadir_azim):
    """
    Receives elevation and azimuth 2d, rotates around base
    so that base_elev <> 0deg and base_azim <> 0deg
    elevation being 0 at horizon (xy plane) and azimuth 0 at x axis

    Returns
    ------
        elevation, azimuth
        (xy plane elevation)
    """
    # translating to normal polar coordinate system, with theta being angle from z axis
    # and phi being angle from x axis in the xy plane
    nadir_theta = 90 - nadir_elev
    nadir_phi = nadir_azim

    nadir_point = np.matrix([
        np.sin(np.deg2rad(nadir_theta)) * np.cos(np.deg2rad(nadir_phi)),
        np.sin(np.deg2rad(nadir_theta)) * np.sin(np.deg2rad(nadir_phi)),
        np.cos(np.deg2rad(nadir_theta)),
    ])
    # first rotate around y axis nadir_theta-180 to reach new theta
    # since nadir_theta in (0,180), rotation will end up to azimuth=0
    # so we rotate it around z axis nadir_phi
    rotation_matrix = get_rotation_matrix(nadir_phi, nadir_theta - 180)

    theta = 90 - elev
    phi = azim

    phi_rad = np.ravel(np.array([np.deg2rad(phi)]))
    theta_rad = np.ravel(np.array([np.deg2rad(theta)]))

    points = np.matrix([
        np.sin(theta_rad) * np.cos(phi_rad),
        np.sin(theta_rad) * np.sin(phi_rad),
        np.cos(theta_rad),
    ])

    rotated_points = rotation_matrix @ points

    rotated_phi = np.ravel(
        np.asarray(
            np.rad2deg(
                np.arctan2(rotated_points[1], rotated_points[0]),
            ),
        ),
    )
    rotated_theta = np.ravel(
        np.asarray(
            np.rad2deg(np.arccos(rotated_points[2])),
        ),
    )

    # back to elevation = 0 deg at xy plane, 90 at zenith and -90 at nadir
    res_elev = 90 - rotated_theta

    return res_elev, rotated_phi


class GeometryConverter():
    """
    This is a Singleton. set_reference should be called once per simulation/snapshot.

    Alert:
        Every conversion between polar and geodesic should be intermediated by cartesian.
        Every transformation should be done either at cartesian or polar.
        Ignore at your own risk (and sadness)
    """
    def __init__(self):
        # geodesical
        self.ref_lat = None
        self.ref_long = None
        self.ref_alt = None

        # cartesian
        self.ref_x = None
        self.ref_y = None
        self.ref_z = None

        # polar
        self.ref_r = None
        self.ref_azim = None
        self.ref_elev = None

    def get_translation(self):
        return self.ref_r

    def validate(self):
        if None in [self.ref_elev, self.ref_azim, self.ref_r]:
            raise ValueError("You need to set a reference for coordinate transformation before using it")

    def set_reference(self, ref_lat: float, ref_long: float, ref_alt: float):
        self.ref_lat = ref_lat
        self.ref_long = ref_long
        self.ref_alt = ref_alt
        ref_x, ref_y, ref_z = lla2ecef(self.ref_lat, self.ref_long, self.ref_alt)
        self.ref_x = ref_x
        self.ref_y = ref_y
        self.ref_z = ref_z

        # polar coordinates
        # geodesic doesn't necessarily translate one to one here
        # so we use cartesian as intermediary
        self.ref_r, self.ref_azim, self.ref_elev = cartesian_to_polar(ref_x, ref_y, ref_z)

    def convert_cartesian_to_transformed_cartesian(
        self, x, y, z, *, translate=None
    ):
        """
        Transforms points by the same transformation required to bring reference to (0,0,0)
        You can only rotate by specifying translate=0
        """
        ref_elev = self.ref_elev
        ref_azim = self.ref_azim
        ref_r = self.ref_r

        self.validate()

        # calculate distances to the centre of the Earth
        dist_sat_centre_earth_km, azim, elev = cartesian_to_polar(x, y, z)

        dist_imt_centre_earth = translate
        if translate is None:
            dist_imt_centre_earth = ref_r

        # calculate Cartesian coordinates of , with origin at centre of the Earth,
        # but considering the x axis at same longitude as the ref_long
        # so that we can rotate around y to bring reference to top
        sat_lat_rad = np.deg2rad(elev)
        # consider coordinates rotating ref_long clockwise around z axis
        imt_long_diff_rad = np.deg2rad(azim - ref_azim)
        x1 = dist_sat_centre_earth_km * \
            np.cos(sat_lat_rad) * np.cos(imt_long_diff_rad)
        y1 = dist_sat_centre_earth_km * \
            np.cos(sat_lat_rad) * np.sin(imt_long_diff_rad)

        # didn't transform, shoud eq height
        z1 = dist_sat_centre_earth_km * np.sin(sat_lat_rad)

        # rotate axis and calculate coordinates with origin at IMT system
        imt_lat_rad = np.deg2rad(ref_elev)
        x2 = (
            x1 * np.sin(imt_lat_rad) - z1 * np.cos(imt_lat_rad)
        )

        y2 = y1

        z2 = (
            z1 * np.sin(imt_lat_rad) + x1 * np.cos(imt_lat_rad)
        ) - dist_imt_centre_earth

        return (x2, y2, z2)

    def revert_transformed_cartesian_to_cartesian(
        self, x2, y2, z2, *, translate=None
    ):
        """
        Reverses transformed points by the same transformation required to bring reference to (0,0,0)
        You can only rotate by specifying translate=0. You need to use the same 'translate' value used
        in transformation if you wish to reverse the transformation correctly
        """
        ref_elev = self.ref_elev
        ref_azim = self.ref_azim
        ref_r = self.ref_r

        self.validate()

        # rotate axis and calculate coordinates with origin at IMT system
        imt_lat_rad = np.deg2rad(ref_elev)

        dist_imt_centre_earth = translate
        if translate is None:
            dist_imt_centre_earth = ref_r

        # transposed transformation matrix
        z2 = z2 + dist_imt_centre_earth
        y1 = y2
        x1 = x2 * np.sin(imt_lat_rad) + z2 * np.cos(imt_lat_rad)
        z1 = z2 * np.sin(imt_lat_rad) - x2 * np.cos(imt_lat_rad)

        dist_sat_centre_earth_km = np.sqrt(x1 * x1 + z1 * z1 + y1 * y1)
        sat_lat_rad = np.arcsin(z1 / dist_sat_centre_earth_km)

        imt_long_diff_rad = np.arctan2(
            y1, x1
        )

        # calculate distances to the centre of the Earth
        x, y, z = polar_to_cartesian(
            dist_sat_centre_earth_km,
            np.rad2deg(imt_long_diff_rad) + ref_azim,
            np.rad2deg(sat_lat_rad)
        )

        return (x, y, z)

    def convert_lla_to_transformed_cartesian(
        self, lat: np.array, long: np.array, alt: np.array
    ):
        """
        You cannot transform this back to lla and expect something useful...
        This rotates and translates every point considering the reference that was set
        and a geodesical coordinate system
        """
        # get cartesian position by geodesical
        x, y, z = lla2ecef(lat, long, alt)

        return self.convert_cartesian_to_transformed_cartesian(x, y, z)

    def convert_station_3d_to_2d(
        self, station: StationManager, idx=None
    ) -> None:
        """
        In place rotate and translate all coordinates so that reference parameters end up in (0,0,0)
        and stations end up in same relative position according to each other,
        adapting their angles to the rotation.

        if idx is specified, only stations[idx] will be converted
        """
        # transform positions
        # print("(station.x, station.y, station.z)", (station.x[idx], station.y[idx], station.z[idx]))
        if idx is None:
            nx, ny, nz = self.convert_cartesian_to_transformed_cartesian(station.x, station.y, station.z)
        else:
            nx, ny, nz = self.convert_cartesian_to_transformed_cartesian(station.x[idx], station.y[idx], station.z[idx])

        if idx is None:
            azim = station.azimuth
            elev = station.elevation
        else:
            azim = station.azimuth[idx]
            elev = station.elevation[idx]

        r = 1
        # then get pointing vec
        pointing_vec_x, pointing_vec_y, pointing_vec_z = polar_to_cartesian(r, azim, elev)

        # transform pointing vectors, without considering geodesical earth coord system
        pointing_vec_x, pointing_vec_y, pointing_vec_z = self.convert_cartesian_to_transformed_cartesian(
            pointing_vec_x, pointing_vec_y, pointing_vec_z, translate=0
        )

        if idx is None:
            station.x = nx
            station.y = ny
            station.z = nz

            _, station.azimuth, station.elevation = cartesian_to_polar(pointing_vec_x, pointing_vec_y, pointing_vec_z)
        else:
            station.x[idx] = nx
            station.y[idx] = ny
            station.z[idx] = nz

            _, azimuth, elevation = cartesian_to_polar(pointing_vec_x, pointing_vec_y, pointing_vec_z)

            station.azimuth[idx] = azimuth
            station.elevation[idx] = elevation

    def revert_station_2d_to_3d(
        self, station: StationManager, idx=None
    ) -> None:
        """
        In place rotate and translate all coordinates so that reference parameters end up in (0,0,0)
        and stations end up in same relative position according to each other,
        adapting their angles to the rotation.

        if idx is specified, only stations[idx] will be converted
        """
        # transform positions
        # print("(station.x, station.y, station.z)", (station.x[idx], station.y[idx], station.z[idx]))
        if idx is None:
            nx, ny, nz = self.revert_transformed_cartesian_to_cartesian(station.x, station.y, station.z)
        else:
            nx, ny, nz = self.revert_transformed_cartesian_to_cartesian(station.x[idx], station.y[idx], station.z[idx])

        if idx is None:
            azim = station.azimuth
            elev = station.elevation
        else:
            azim = station.azimuth[idx]
            elev = station.elevation[idx]

        r = 1
        # then get pointing vec
        pointing_vec_x, pointing_vec_y, pointing_vec_z = polar_to_cartesian(r, azim, elev)

        # transform pointing vectors, without considering geodesical earth coord system
        pointing_vec_x, pointing_vec_y, pointing_vec_z = self.revert_transformed_cartesian_to_cartesian(
            pointing_vec_x, pointing_vec_y, pointing_vec_z, translate=0
        )

        if idx is None:
            station.x = nx
            station.y = ny
            station.z = nz

            _, station.azimuth, station.elevation = cartesian_to_polar(pointing_vec_x, pointing_vec_y, pointing_vec_z)
        else:
            station.x[idx] = nx
            station.y[idx] = ny
            station.z[idx] = nz

            _, azimuth, elevation = cartesian_to_polar(pointing_vec_x, pointing_vec_y, pointing_vec_z)

            station.azimuth[idx] = azimuth
            station.elevation[idx] = elevation


def get_lambert_equal_area_crs(polygon: shp.geometry.Polygon):
    centroid = polygon.centroid
    return pyproj.CRS.from_user_input(
        f"+proj=laea +lat_0={centroid.y} +lon_0={centroid.x} +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    )


def shrink_country_polygon_by_km(
    polygon: shp.geometry.Polygon, km: float
) -> shp.geometry.Polygon:
    """
    Projects a Polygon in "EPSG:4326" to Lambert Azimuthal Equal Area projection,
    shrinks the polygon by x km,
    projects the polygon back to EPSG:4326.

    Hint:
        Check for polygon validity after transformation:
        if poly.is_valid: raise Exception("bad polygon")
        if not poly.is_empty and poly.area > 0: continue # ignore
        ...
    """
    # Lambert is more precise, but could prob. get UTM projection
    # Didn't see any practical difference for current use cases
    proj_crs = get_lambert_equal_area_crs(polygon)

    # Create transformer objects
    # NOTE: important always_xy=True to not mix lat lon up order
    to_proj = pyproj.Transformer.from_crs("EPSG:4326", proj_crs, always_xy=True).transform
    from_proj = pyproj.Transformer.from_crs(proj_crs, "EPSG:4326", always_xy=True).transform

    # Transform to projection where unit is meters
    polygon_proj = shp.ops.transform(to_proj, polygon)

    # Shrink (negative buffer in meters)
    polygon_proj_shrunk = polygon_proj.buffer(-km * 1000)

    # Return to EPSG:4326
    return shp.ops.transform(from_proj, polygon_proj_shrunk)


def shrink_countries_by_km(
    countries: list[shp.geometry.MultiPolygon],
    km: float
) -> list[shp.geometry.MultiPolygon]:
    """
    Receives a MultiPolygon containing multiple countries
    and diminishes
    """
    polys = []

    for ext_poly in countries:
        if ext_poly.geom_type == 'Polygon':
            polys.append(shrink_country_polygon_by_km(ext_poly, km))
        elif ext_poly.geom_type == 'MultiPolygon':
            polys.append(shp.ops.unary_union([
                shrink_country_polygon_by_km(poly, km) for poly in ext_poly.geoms
            ]))

    for poly in polys:
        if not poly.is_valid:
            # may be ignorable..?
            # TODO: check if this error can be safely removed
            # If you need to look into this, plot the erroring scenario
            # drops and check to see if unexpected behavior occurs
            raise ValueError("SOME BAD POLYGON")

    return [
        poly for poly in polys if poly.is_valid and not poly.is_empty and poly.area > 0
    ]


if __name__ == "__main__":
    # baixo, direita, frente, esquerda, atrás, cima, cima
    # elev = np.array([-90., 0., 0., 0., 0., 90., 90.])
    # azim = np.array([0., 0., 90., 180., -90., 0., 90.])
    elev = np.array([-89.93761622])
    azim = np.array([180.])
    n = len(azim)

    # print("pointing nadir to the right (rotated around earth to the left)")
    # print(np.concatenate((elev, azim)).reshape((2, n)).transpose())
    # res_elev, res_az = rotate_angles_based_on_new_nadir(elev, azim, 0, 0)
    # res = np.concatenate((res_elev, res_az)).reshape((2, n)).transpose()
    # print(res)

    # print("pointing nadir to the left (rotated around earth to the right)")
    # print(np.concatenate((elev, azim)).reshape((2, n)).transpose())
    # res_elev, res_az = rotate_angles_based_on_new_nadir(elev, azim, 0, 180)
    # res = np.concatenate((res_elev, res_az)).reshape((2, n)).transpose()
    # print(res)

    # print(get_rotation_matrix_between_vecs(np.array([0,1,0]), np.array([0,0,1])))

    geoconv = GeometryConverter()

    sys_lat = 89
    sys_long = 0
    sys_alt = 1200

    # geoconv.set_reference(
    #     sys_lat, sys_long, sys_alt
    # )
    # stat = StationManager(1)
    # # stat.x = np.array([-2000.])
    # # stat.y = np.array([0.])
    # stat.x = np.array([0.])
    # stat.y = np.array([-2000.])
    # stat.z = np.array([0.])
    # stat.elevation = 0
    # stat.azimuth = -90

    # print("stat.x", stat.x)
    # print("stat.y", stat.y)
    # print("stat.z", stat.z)

    # print("stat.azimuth", stat.azimuth)
    # print("stat.elevation", stat.elevation)
    # print("#########")
    # geoconv.convert_station_3d_to_2d(stat)
    # print("#########")

    # print("stat.x", stat.x)
    # print("stat.y", stat.y)
    # print("stat.z", stat.z)

    # print("stat.azimuth", stat.azimuth)
    # print("stat.elevation", stat.elevation)
