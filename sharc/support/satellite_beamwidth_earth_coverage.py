import numpy as np
import matplotlib.pyplot as plt

"""
Author: Luciano Camilo Alexandre (luciano-camilo@hotmail.com)

Date: Thu 31 May 14:00:00 2021

    Python script to calculate Low Earth Orbit satellite/pseudo-satellite beamwidth Earth coverage.
"""

re = 6378                                     # Earth radius [km]
satellite_altitude = 510                      # Satellite altitude
satellite_elevation = np.linspace(0,90,100)   # Satellite elevation vector
#satellite_elevation = 10


SEC = satellite_elevation + 90                # Satellite-Earth-Center angle
rs = re + satellite_altitude
a = (np.sin(np.radians(SEC))/rs)*re
a=np.degrees(np.arcsin(a))
b= 180 - SEC - a
distance = re * np.radians(b)
coverage= 2 * distance
#print(f'LEO Satellite coverage is {coverage} km')
earth_diameter = np.pi *6378*2
#print(earth_diameter)
diameter_coverage_LEO = (2*np.radians(b))/(2*np.pi)
diameter_coverage_LEO = diameter_coverage_LEO*earth_diameter
angle = ((2*np.radians(b))/(2*np.pi))*360
plt.title("LEO Coverage diameter arc")
plt.plot(satellite_elevation,diameter_coverage_LEO,label=f'LEO altitude: {satellite_altitude} km ')
plt.xlabel('Elevation angle [degrees]', fontsize=10, color='black')
plt.ylabel('Coverage diameter arc [km]', fontsize=10, color='black')
plt.legend()
plt.grid()
plt.show()
plt.title("LEO Full angle of the antenna beamwidth coverage")
plt.plot(2*a,diameter_coverage_LEO, label=f'LEO altitude: {satellite_altitude}')
plt.xlabel('Full angle of the antenna beamwidth at the satellite [degrees]', fontsize=10, color='black')
plt.ylabel('Coverage diameter arc [km]', fontsize=10, color='black')
plt.grid()
plt.legend()
plt.show()
plt.title("LEO Full angle of the antenna beamwidth coverage")
plt.plot(2*a,satellite_elevation,label=f'LEO altitude: {satellite_altitude}')
plt.xlabel('Full angle of the antenna beamwidth at the satellite [degrees]', fontsize=10, color='black')
plt.ylabel('Elevation angle [degrees]', fontsize=10, color='black')
plt.grid()
plt.legend()
plt.show()


