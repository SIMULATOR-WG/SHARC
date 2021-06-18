"""
@Created: Luciano Camilo on Tue Fr 04 18:27:25 2021

"""
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(6, 5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height])

"""
Calculate the vertical antenna pattern (elevation) of Cosecant squared radar antenna.

    Construction:
            g_max: maximum gain
            theta_3db : 3dB beamwidth
            maximum_csc = maximum cosecant squared angle
            highbeam_csc2 = initial cosecant squared angle for high beam in CSC^2 pattern
"""
g_max = 33.5
theta_3db = 5
maximum_csc = 40
highbeam_csc2 = 5

theta = np.linspace(-90, 90, 100)
const = np.pi * 50.8
g = np.zeros(len(theta))
mi = np.zeros(len(theta))
theta = theta-highbeam_csc2
maximum_csc = maximum_csc-highbeam_csc2

for i in range(-90, len(theta)):
    if -theta_3db / 0.88 <= theta[i] <= theta_3db:
        mi[i] = const * np.sin(np.radians(theta[i])) / np.radians(theta_3db)
        g[i] = 20 * np.log10(np.sin(np.radians(mi[i])) / mi[i]) + 35.1623
    if theta_3db <= theta[i] <= maximum_csc:
        g1 = ((np.sin(np.radians(const * np.sin(np.radians(theta_3db)) / np.radians(theta_3db))))
              / (const * np.sin(np.radians(theta_3db)) / np.radians(theta_3db))) + 0.12
        g[i] = 20 * np.log10(
            g1 * (((1 / np.sin(np.radians(theta[i]))) / (1 / np.sin(np.radians(theta_3db)))) ** 2))
    if maximum_csc <= theta[i] <= 90:
        g[i] = -55
    if theta[i] < -theta_3db / 0.88 or theta[i] >= 90:
        g[i] = -55

"""
Calculate the horizontal antenna pattern (azimuth) of cosine squared radar antenna.

    Construction:
            phi_3db : 3dB beamwidth
"""
phi = np.linspace(-90, 90, 90)
const = np.pi * 68.8
g1 = np.zeros(len(phi))
mi1 = np.zeros(len(phi))
phi_3db = 1.5

for i in range(len(phi)):
    if -phi_3db <= phi[i] < phi_3db:
        mi1[i] = const * np.sin(np.radians(phi[i])) / phi_3db
        g1[i] = 20 * np.log10((np.pi / 2) * ((np.cos(mi[i])) / ((np.pi / 2) ** 2 - (mi[i]) ** 2))) + 4.32 - 0.392
    if phi_3db <= phi[i] <= 180 or -phi_3db >= phi[i] >= -180:
        g1[i] = -17.51 * np.log(2.33 * (np.abs((phi[i])) / phi_3db)) - 4.32  # cos2
        if g1[i] < -50:
            g1[i] = -50

# Plot 3D Contour Surface
x, y = np.meshgrid(np.linspace(-90, 90, 100), np.linspace(-90, 90, 90))
X, Y = np.meshgrid(g, g1, sparse=False)
clevs1 = np.arange(-75, 34.1, 0.1)
cs1 = plt.contourf(y, x, X+Y+g_max, clevs1, cmap='viridis', alpha=1)
plt.colorbar(cs1)
ax.set_title('Surveillance Radar Antenna Pattern - CSC2')
ax.set_xlabel('azimuth [degrees]')
ax.set_ylabel('elevation [degrees]')
plt.show()
