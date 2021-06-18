import numpy as np
from matplotlib import pyplot as plt

"""
Author: Luciano Camilo Alexandre (luciano-camilo@hotmail.com)

Date: Thu 17 Jun 22:00:00 2021

    Python script to calculate satellite/pseudo-satellite horizon visibility according Rec. ITU-R P.619 Annex I
"""

#x = -4.68 # longitude difference to calculate distance [degrees]
x = 1.8 # longitude difference to calculate distance [degrees]
#x = -30 # longitude difference to calculate distance [degrees]

EARTH_RADIUS = 6371000              # [meters]
satellite_altitude = 20000         # [meters]
#satellite_altitude = 510000         # [meters]
system_altitude = 8                 # [meters]


satellite_lat_deg = -15.80      # [degrees]
satellite_long_diff_deg = np.linspace(0, x, 181)  # [degrees]
system_lat_deg = -15.80      # [degrees]

# Creating vectors
dist = np.linspace(0,180,181)
height = np.linspace(0,180,181)
height_earth_space_angle = np.linspace(0, 180, 181)
system_long_diff_rad = np.linspace(0, 180, 181)
free_space_angle = np.linspace(0,5,181)
free_space_angle_earth_space = np.linspace(0, 5, 181)
tau_fs = np.linspace(0, 180, 181)
tau_fs1 = np.linspace(0, 180, 181)
tau_fs2 = np.linspace(0, 180, 181)
tau_fs3 = np.linspace(0, 180, 181)
tau_fs_deg = np.linspace(0, 180, 181)
aparent_angle = np.linspace(0, 180, 181)
theta_0 = np.linspace(0,180,181)
theta_earth_space = np.linspace(0,180,181)
angle_trigonometry = np.linspace(0,180,181)
maximum_distance = 500000 # maximum distance evaluated in trigonometry


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
    height_earth_space_angle[i] = satellite_altitude - height[i]*1000
    #print(height_earth_space_angle[i])

    dts= np.sqrt(np.power(x, 2) + np.power(y, 2) + np.power(height, 2))
    gts = np.sqrt(np.power(x, 2) + np.power(y, 2))
    #print(gts)

    theta_0[i] = np.arctan2(height[i]*1000, gts)  # free-space elevation angle
    #print(np.degrees(theta_0[i]))

    phi = np.arctan2(x,y)
    phi_azimuth = np.degrees((phi))

    dist[i] = 6378.388 * np.arccos(np.sin(sat_lat_rad) * np.sin(sys_lat_rad ) + np.cos(sat_lat_rad) * np.cos(sys_lat_rad ) * np.cos(system_long_diff_rad[i]))

    tau_fs1[i] = 1.728 + 0.5411 * theta_0[i] + 0.03723 * theta_0[i] ** 2
    tau_fs2[i] = 0.1815 + 0.06272 * theta_0[i] + 0.01380 * theta_0[i] ** 2
    tau_fs3[i] = 0.01727 + 0.008288 * theta_0[i]

    # change in elevation angle due to refraction
    tau_fs_deg[i] = 1 / (tau_fs1[i] + system_altitude * tau_fs2[i] + system_altitude**2 * tau_fs3[i])
    tau_fs[i] = tau_fs_deg[i] / 180. * np.pi

    aparent_angle[i] = np.degrees(theta_0[i] + tau_fs[i])
    #print(aparent_angle[i])
    #print(np.degrees(angle_trigonometry[i]))

    #print(dts)
    #print(gts)
    #print(dts)
    #print(phi_azimuth)
    #print(free_s

dist2= np.linspace(0,maximum_distance,181)

for i in range(len(dist2)):
     #print(dist2[i])
     #dist2[i]=np.sqrt((20000**2)+(dist2[i])**2)
     angle_trigonometry[i]=np.arctan2(np.sqrt(dist2[i]**2),satellite_altitude)
     #print(dist2[i])
     #print(np.degrees(angle_trigonometry[i]))

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
plt.title('Earth-to-space path angle')
plt.xlabel('Distance from satellite [km]')
plt.ylabel('Angle [degrees] - Earth to space path')
plt.plot(dist, np.degrees(theta_0), 'b--', label='Free space elevation angle')
plt.plot(dist, aparent_angle, 'r-', label='Apparent elevation angle')
plt.plot(dist, 90-np.degrees(angle_trigonometry), 'g-', label='Trigonometry')
plt.legend()
plt.show()

