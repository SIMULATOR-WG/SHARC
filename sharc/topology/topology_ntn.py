from sharc.topology.topology import Topology
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.axes
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon


class TopologyNTN(Topology):
    """
    Class to generate and manage the topology of Non-Terrestrial Network (NTN) sites
    based on a specified macrocell network topology.
    """

    ALLOWED_NUM_SECTORS = [1, 7, 19]

    def __init__(
        self, intersite_distance: float, cell_radius: int, bs_height: float, bs_azimuth: float,
        bs_elevation: float, num_sectors=7,
    ):
        """
        Initializes the NTN topology with specific network settings.

        Parameters:
        intersite_distance: Distance between adjacent sites in meters.
        cell_radius: Radius of the coverage area for each site in meters.
        bs_height: Height of the base station (satellite) from the x-y plane.
        bs_azimuth: Azimuth angle of the base station in degrees.
        bs_elevation: Elevation angle of the base station in degrees.
        num_sectors: Number of sectors for the topology (default is 7).
        """

        if num_sectors not in self.ALLOWED_NUM_SECTORS:
            raise ValueError(
                f"Invalid number of sectors: {num_sectors}. Allowed values are {self.ALLOWED_NUM_SECTORS}.",
            )

        # Call to the superclass constructor to set common properties
        super().__init__(intersite_distance, cell_radius)
        self.is_space_station = True
        self.space_station_x = None
        self.space_station_y = None
        self.space_station_z = None
        self.bs_azimuth = np.radians(bs_azimuth)
        self.bs_elevation = np.radians(bs_elevation)
        self.bs_radius = bs_height / np.sin(self.bs_elevation)
        self.num_sectors = num_sectors

        # Calculate the base station coordinates

        self.space_station_x = self.bs_radius * \
            np.cos(self.bs_elevation) * np.cos(self.bs_azimuth)
        self.space_station_y = self.bs_radius * \
            np.cos(self.bs_elevation) * np.sin(self.bs_azimuth)
        self.space_station_z = bs_height

        self.calculate_coordinates()

    def calculate_coordinates(self, random_number_gen=np.random.RandomState()):
        """
        Computes the coordinates of each site. This is where the actual layout calculation would be implemented.
        """

        d = self.intersite_distance

        self.x = [0]
        self.y = [0]

        # First ring (6 points)
        if self.num_sectors == 7 or self.num_sectors == 19:

            for k in range(6):
                angle = k * 60
                self.x.append(d * np.cos(np.radians(angle)))
                self.y.append(d * np.sin(np.radians(angle)))

        if self.num_sectors == 19:
            # Coordinates with 19 sectors
            # Second ring (12 points)
            for k in range(6):
                angle = k * 60
                self.x.append(2 * d * np.cos(np.radians(angle)))
                self.y.append(2 * d * np.sin(np.radians(angle)))
                self.x.append(
                    d * np.cos(np.radians(angle)) +
                    d * np.cos(np.radians(angle + 60)),
                )
                self.y.append(
                    d * np.sin(np.radians(angle)) +
                    d * np.sin(np.radians(angle + 60)),
                )

        self.x = np.array(self.x)
        self.y = np.array(self.y)

        # Assuming all points are at ground level
        self.z = np.zeros(len(self.x))

        # Rotate the anchor points by 30 degrees
        theta = np.radians(30)
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        self.x_rotated = self.x * cos_theta - self.y * sin_theta
        self.y_rotated = self.x * sin_theta + self.y * cos_theta

        # Calculate azimuth and elevation for each point
        self.azimuth = np.arctan2(
            self.y_rotated - self.space_station_y,
            self.x_rotated - self.space_station_x,
        ) * 180 / np.pi
        distance_xy = np.sqrt(
            (self.x_rotated - self.space_station_x) **
            2 + (self.y_rotated - self.space_station_y)**2,
        )
        self.elevation = np.arctan2(
            self.z - self.space_station_z, distance_xy,
        ) * 180 / np.pi

        # Update the number of base stations after setup
        self.num_base_stations = len(self.x)
        self.indoor = np.zeros(self.num_base_stations, dtype=bool)

        self.x = self.x_rotated
        self.y = self.y_rotated

    def plot(self, axis: matplotlib.axes.Axes):
        r = self.cell_radius / 1000  # Convert to kilometers

        # Plot each sector
        for x, y in zip(self.x / 1000, self.y / 1000):  # Convert to kilometers
            hexagon = []
            for a in range(6):
                angle_rad = math.radians(a * 60)
                hexagon.append([
                    x + r * math.cos(angle_rad),
                    y + r * math.sin(angle_rad),
                ])
            hexagon.append(hexagon[0])  # Close the hexagon
            hexagon = np.array(hexagon)

            sector = plt.Polygon(hexagon, fill=None, edgecolor='k')
            axis.add_patch(sector)

        # Plot base stations
        axis.scatter(
            self.x / 1000, self.y / 1000, s=200, marker='v', c='k', edgecolor='k', linewidth=1, alpha=1,
            label="Anchor Points",
        )

        # Add labels and title
        axis.set_xlabel("x-coordinate [km]")
        axis.set_ylabel("y-coordinate [km]")
        axis.set_title(f"NTN Topology - {self.num_sectors} Sectors")
        axis.legend()
        plt.tight_layout()

    def plot_3d(self, axis: matplotlib.axes.Axes, map=False):
        r = self.cell_radius / 1000  # Convert to kilometers

        if map:
            # Load the map of Brazil using GeoPandas
            brazil = gpd.read_file(
                "${workspaceFolder}\\sharc\\topology\\countries\\ne_110m_admin_0_countries.shp",
            )
            brazil = brazil[brazil['NAME'] == "Brazil"]

            # Coordinates of the Federal District (Brasília)
            federal_district_coords = (-47.9292, -15.7801)

            # Approximate conversion factors (1 degree latitude = 111 km, 1 degree longitude = 111 km)
            lat_to_km = 111
            lon_to_km = 111

            # Convert Federal District coordinates to kilometers
            federal_district_coords_km = (
                federal_district_coords[0] * lon_to_km, federal_district_coords[1] * lat_to_km,
            )

            # Calculate the shift required to move the Federal District to (0, 0)
            x_shift = federal_district_coords_km[0]
            y_shift = federal_district_coords_km[1]

            # Manually plot the map of Brazil on the xy-plane
            for geom in brazil.geometry:
                if isinstance(geom, Polygon):
                    lon, lat = geom.exterior.xy
                    x = np.array(lon) * lon_to_km - x_shift
                    y = np.array(lat) * lat_to_km - y_shift
                    axis.plot(x, y, zs=0, zdir='z', color='lightgray')
                elif isinstance(geom, MultiPolygon):
                    for poly in geom:
                        lon, lat = poly.exterior.xy
                        x = np.array(lon) * lon_to_km - x_shift
                        y = np.array(lat) * lat_to_km - y_shift
                        axis.plot(x, y, zs=0, zdir='z', color='lightgray')

            # Add the Federal District location to the plot
            axis.scatter(0, 0, 0, color='red', zorder=5)
            axis.text(0, 0, 0, 'Federal District', fontsize=12, ha='right')

        # Plot each sector
        for x, y in zip(self.x / 1000, self.y / 1000):  # Convert to kilometers
            hexagon = []
            for a in range(6):
                angle_rad = math.radians(a * 60)
                hexagon.append([
                    x + r * math.cos(angle_rad),
                    y + r * math.sin(angle_rad),
                ])
            hexagon.append(hexagon[0])  # Close the hexagon
            hexagon = np.array(hexagon)

            # 3D hexagon
            axis.plot(
                hexagon[:, 0], hexagon[:, 1],
                np.zeros_like(hexagon[:, 0]), 'k-',
            )

        # Plot base stations
        axis.scatter(
            self.x / 1000, self.y / 1000, np.zeros_like(self.x), s=75, marker='v', c='k', edgecolor='k',
            linewidth=1, alpha=1,
            label="Anchor Points",
        )

        # Plot the satellite
        axis.scatter(
            self.space_station_x / 1000, self.space_station_y / 1000, self.space_station_z / 1000, s=75, c='r',
            marker='^', edgecolor='k', linewidth=1, alpha=1,
            label=f"Satellite (φ={np.degrees(self.bs_azimuth):.1f}°, θ={np.degrees(self.bs_elevation):.1f}°)",
        )

        # Plot the height line
        axis.plot(
            [self.space_station_x / 1000, self.space_station_x / 1000],
            [self.space_station_y / 1000, self.space_station_y / 1000],
            [0, self.space_station_z / 1000], 'b-', label=f'Height = {self.space_station_z / 1000:.1f} km',
        )

        # Plot the slant range line
        axis.plot(
            [0, self.space_station_x / 1000],
            [0, self.space_station_y / 1000],
            [0, self.space_station_z / 1000], 'g--', label=f'Slant range = {self.bs_radius / 1000:.1f} km',
        )

        # Add labels and title
        axis.set_xlabel("x-coordinate [km]")
        axis.set_ylabel("y-coordinate [km]")
        axis.set_zlabel("z-coordinate [km]")
        axis.set_title(f"NTN Topology - {self.num_sectors} Sectors")
        axis.legend()
        plt.tight_layout()


# Example usage
if __name__ == '__main__':

    bs_height = 1000e3  # meters
    bs_azimuth = 45  # degrees
    bs_elevation = 45  # degrees
    beamwidth = 10
    cell_radius = bs_height * \
        math.tan(np.radians(beamwidth)) / math.cos(bs_elevation)
    intersite_distance = cell_radius * np.sqrt(3)  # meters
    map = False

    # Test for 1 sector
    ntn_topology_1 = TopologyNTN(
        intersite_distance, cell_radius, bs_height, bs_azimuth, bs_elevation, num_sectors=1,
    )
    ntn_topology_1.calculate_coordinates()  # Calculate the site coordinates

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ntn_topology_1.plot_3d(ax, map)  # Plot the 3D topology
    plt.show()

    # Test for 7 sectors
    ntn_topology_7 = TopologyNTN(
        intersite_distance, cell_radius, bs_height, bs_azimuth, bs_elevation, num_sectors=7,
    )
    ntn_topology_7.calculate_coordinates()  # Calculate the site coordinates

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ntn_topology_7.plot_3d(ax, map)  # Plot the 3D topology
    plt.show()

    # Test for 19 sectors
    ntn_topology_19 = TopologyNTN(
        intersite_distance, cell_radius, bs_height, bs_azimuth, bs_elevation, num_sectors=19,
    )
    ntn_topology_19.calculate_coordinates()  # Calculate the site coordinates

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ntn_topology_19.plot_3d(ax, map)  # Plot the 3D topology
    plt.show()
