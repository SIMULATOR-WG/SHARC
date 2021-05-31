"""
@Created: Luciano Camilo on Tue Fr 04 18:27:25 2021

"""

import matplotlib.pyplot as plt
import numpy as np

# Antenna pattern - Vertical (Elevation)
g_max = 38
theta_3db = 2

theta = np.linspace(-90, 90, 200)
const = np.pi * 68.8
g = np.zeros(len(theta))
mi = np.zeros(len(theta))

for i in range(len(theta)):
    if theta_3db <= theta[i] <= theta_3db:
        mi[i] = const * np.sin(np.radians(theta[i])) / theta_3db
        g[i] = 20 * np.log10((np.pi / 2) * ((np.cos(mi[i])) / ((np.pi / 2) ** 2 - (mi[i]) ** 2))) + 4.32 - 0.392
    if theta_3db <= theta[i] <= 180 or theta_3db >= theta[i] >= -180:
        g[i] = -17.51 * np.log(2.33 * (np.abs((theta[i])) / theta_3db)) - 5.6  # cos2
        if g[i] < -50:
            g[i] = -50

# Antenna pattern - Horizontal (Azimuth)
phi = np.linspace(-90, 90, 200)
const = np.pi * 68.8
g1 = np.zeros(len(phi))
mi1 = np.zeros(len(phi))
phi_3db = 2

for i in range(len(phi)):

    if -phi_3db <= phi[i] < phi_3db:
        mi1[i] = const * np.sin(np.radians(phi[i])) / phi_3db
        g1[i] = 20 * np.log10((np.pi / 2) * ((np.cos(mi[i])) / ((np.pi / 2) ** 2 - (mi[i]) ** 2))) + 4.32 - 0.392
    if phi_3db <= phi[i] <= 180 or -phi_3db >= phi[i] >= -180:
        g1[i] = -17.51 * np.log(2.33 * (np.abs((phi[i])) / phi_3db)) - 5.6  # cos2
        if g1[i] < -50:
            g1[i] = -50

# Plot 3D Surface
x, y = np.meshgrid(np.linspace(-90, 90, 200), np.linspace(-90, 90, 200))
X, Y = np.meshgrid(g, g1, sparse=False)
fig = plt.figure(figsize=(10, 10))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, X+Y+g_max, cstride=1, rstride=1,  cmap='rainbow', antialiased=True, alpha=1)
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=20, ticks=np.linspace(-60, 38, 5))
ax.set_title('Meteorological Radar Antenna Pattern - Cosine Raised n=1')
ax.set_xlabel('azimuth [degrees]')
ax.set_ylabel('elevation [degrees]')
plt.show()
