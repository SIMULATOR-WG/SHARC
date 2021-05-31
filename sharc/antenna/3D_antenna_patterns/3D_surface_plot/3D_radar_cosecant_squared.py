"""
@Created: Luciano Camilo on Tue Fr 04 18:27:25 2021

"""
import matplotlib.pyplot as plt
import numpy as np

# Antenna pattern - Vertical (Elevation)
g_max = 33.5
theta_3db = 5
maximum_csc = 40
highbeam_csc2 = 0
maximum_csc = maximum_csc-highbeam_csc2

theta = np.linspace(-90, 90, 100)
const = np.pi * 50.8
g = np.zeros(len(theta))
mi = np.zeros(len(theta))
theta = theta-highbeam_csc2

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

g = g


# Antenna pattern - Horizontal (Azimuth)
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

# Plot 2D Contour Surface
x, y = np.meshgrid(np.linspace(-90, 90, 100), np.linspace(-90, 90, 90))
X, Y = np.meshgrid(g, g1, sparse=False)
fig = plt.figure(figsize=(10, 10))
ax = fig.gca(projection='3d')
clevs1 = np.arange(-72, 40, 1)
surf = ax.plot_surface(x, y, X+Y+g_max, cstride=1, rstride=1,  cmap='rainbow', antialiased=True, alpha=1)
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=20, ticks=np.linspace(-72, 33, 20))
ax.grid(True)
ax.set_title('Surveillance Radar Antenna Pattern - CSC2')
ax.set_xlabel('azimuth [degrees]')
ax.set_ylabel('elevation [degrees]')
plt.show()
