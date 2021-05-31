"""
@Created: Luciano Camilo on Tue Fr 04 18:27:25 2021

"""

import matplotlib.pyplot as plt
import numpy as np

# Antenna pattern - Vertical (Elevation)
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

# Antenna pattern - Horizontal (Azimuth)
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

# Plot 3D Surface
x, y = np.meshgrid(np.linspace(-90, 90, 300), np.linspace(-90, 90, 300))
X, Y = np.meshgrid(g, g1, sparse=False)
fig = plt.figure(figsize=(10, 10))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, X+Y+g_max, cstride=1, rstride=1,  cmap='rainbow', antialiased=True, alpha=1)
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=20, ticks=np.linspace(-60, 0, 20))
ax.set_title('Meteorological Radar Antenna Pattern - Uniform')
ax.set_xlabel('azimuth [degrees]')
ax.set_ylabel('elevation [degrees]')
plt.show()
