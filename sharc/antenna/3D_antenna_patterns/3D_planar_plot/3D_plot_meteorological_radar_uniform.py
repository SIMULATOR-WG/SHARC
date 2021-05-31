"""
@Created: Luciano Camilo on Tue Fr 04 18:27:25 2021

"""
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(6, 5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height])

"""
Calculate the vertical antenna pattern (elevation) of uniform radar antenna.

    Construction:
            g_max: maximum gain
            theta_3db : 3dB beamwidth
"""
g_max = 0
theta_3db = 2
theta = np.linspace(-90, 90, 300)
const = np.pi * 50.8
g = np.zeros(len(theta))
mi = np.zeros(len(theta))

for i in range(len(theta)):
    if -theta_3db <= theta[i] <= theta_3db:
        mi[i] = const * np.sin(np.radians(theta[i])) / theta_3db
        g[i] = 20 * np.log10((np.pi / 2) * ((np.cos(mi[i])) / ((np.pi / 2) ** 2 - (mi[i]) ** 2))) + 4.32 - 0.392
    if theta_3db <= theta[i] <= 180 or -theta_3db >= theta[i] >= -180:
        g[i] = -8.584 * np.log(2.876 * (np.abs((theta[i])) / theta_3db)) - 3.72  # cos2
        if g[i] < -30:
            g[i] = -30

"""
Calculate the horizontal antenna pattern (azimuth) of uniform radar antenna.

    Construction:
            g_max: maximum gain
            theta_3db : 3dB beamwidth
"""
phi = np.linspace(-90, 90, 300)
phi_3db = 2
const = np.pi * 50.8
g1 = np.zeros(len(phi))
mi = np.zeros(len(phi))

for i in range(len(phi)):
    if -phi_3db <= phi[i] <= phi_3db:
        mi[i] = const * np.sin(np.radians(phi[i])) / phi_3db
        g1[i] = 20 * np.log10((np.pi / 2) * ((np.cos(mi[i])) / ((np.pi / 2) ** 2 - (mi[i]) ** 2))) + 4.32 - 0.392
    if phi_3db <= phi[i] <= 180 or -phi_3db >= phi[i] >= -180:
        g1[i] = -8.584 * np.log(2.876 * (np.abs((phi[i])) / phi_3db)) - 3.72  # cos2
        if g1[i] < -30:
            g1[i] = -30

# Plot 3D Contour Surface
x, y = np.meshgrid(np.linspace(-90, 90, 300), np.linspace(-90, 90, 300))
X, Y = np.meshgrid(g, g1, sparse=False)
clevs1 = np.arange(-60, 1, 1)
cs1 = plt.contourf(y, x, X+Y+g_max, clevs1, cmap='viridis', alpha=1)
plt.colorbar(cs1)
ax.set_title('Meteorological Radar Antenna Pattern - Uniform')
ax.set_xlabel('azimuth [degrees]')
ax.set_ylabel('elevation [degrees]')
plt.show()
