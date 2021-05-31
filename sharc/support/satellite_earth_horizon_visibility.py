"""
Author: Luciano Camilo Alexandre (luciano-camilo@hotmail.com)

Date: Thu 29 March 09:00:00 2021

    Python script to calculate satellite/pseudo-satellite horizon visibility according Rec. ITU-R P.619 Annex I
"""

import numpy as np
from matplotlib import pyplot as plt

x = 23 # longitude difference to calculate distance [degrees]

EARTH_RADIUS = 6371000              # [meters]
satellite_altitude = 510000         # [meters]
system_altitude = 0                 # [meters]


satellite_lat_deg = -15.809422      # [degrees]
satellite_long_diff_deg = np.linspace(0, x, 181)  # [degrees]
system_lat_deg = -15.809422         # [degrees]

# Creating vectors
dist = np.linspace(0,180,181)
height = np.linspace(0,180,181)
system_long_diff_rad = np.linspace(0, 180, 181)
free_space_angle = np.linspace(0,5,181)
theta_0 = np.linspace(0,5,181)

# Calculate distances to the centre of the Earth
for i in range(len(satellite_long_diff_deg)):
    dist_hibs_centre_earth_km = (EARTH_RADIUS + satellite_altitude) / 1000
    dist_system_centre_earth_km = (EARTH_RADIUS + system_altitude) /1000

    sat_lat_rad = satellite_lat_deg * np.pi / 180.
    system_long_diff_rad[i] = satellite_long_diff_deg[i] * np.pi / 180.
    x1 = dist_hibs_centre_earth_km * np.cos(sat_lat_rad) * np.cos(system_long_diff_rad[i])
    y1 = dist_hibs_centre_earth_km * np.cos(sat_lat_rad) * np.sin(system_long_diff_rad[i])
    z1 = dist_hibs_centre_earth_km * np.sin(sat_lat_rad)

# Rotate axis and calculate coordinates with origin at System
    sys_lat_rad = system_lat_deg * np.pi / 180.
    x = (x1 * np.sin(sys_lat_rad) - z1 * np.cos(sys_lat_rad))*1000
    y = y1*1000
    height[i] = (z1 * np.sin(sys_lat_rad) + x1 * np.cos(sys_lat_rad)- dist_system_centre_earth_km)

    dts= np.sqrt(np.power(x, 2) + np.power(y, 2) + np.power(height, 2))
    gts = np.sqrt(np.power(x, 2) + np.power(y, 2))

    theta_0 = np.arctan2(gts,height)  # free-space elevation angle
    free_space_angle[i] = np.degrees(theta_0[i])

    phi = np.arctan2(x,y)
    phi_azimuth = np.degrees((phi))


    dist[i] = 6378.388 * np.arccos(np.sin(sat_lat_rad) * np.sin(sys_lat_rad ) + np.cos(sat_lat_rad) * np.cos(sys_lat_rad ) * np.cos(system_long_diff_rad[i]))

    #if (height[i] <= 0):
     #   break
    #print(x,y,height)
    #print(height)
    #print(dts)
    #print(gts)
#print(dts)
    #print(phi_azimuth)
    #print(free_s

plt.figure(1)
plt.legend(loc='upper left')
plt.plot(dist, height, 'b--')
plt.grid(which='minor', alpha=0.2)
plt.grid(which='major', alpha=0.5)
plt.grid(True, color='b', linestyle='--', linewidth=0.2)
plt.title('Horizon Visibility')
plt.xlabel('Distance from satellite [km]')
plt.ylabel('Altitude difference [km]')
plt.legend()
plt.show()

