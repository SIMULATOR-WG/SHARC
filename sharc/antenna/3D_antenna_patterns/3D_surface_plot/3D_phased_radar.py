"""
@Created: Luciano Camilo on Tue Fr 04 18:27:25 2021

"""

import matplotlib.pyplot as plt
import numpy as np
# Antenna pattern - Vertical (Elevation)
g_max = 41
theta_3db = 1.1
phi_3db = 1.1
frequency = 2700
l_ambda = 299792458 / (frequency * 1e6)
element_space = 0.5 * l_ambda
beamsteeringangle_az = 0
beamsteeringangle_el = 0

if g_max == 41:
    normalization = 0.833
    number_elements = 38
if g_max == 46:
    normalization = 3.48681
    number_elements = 70

theta = np.linspace(-89, 89, 200)
const = np.pi * 50.8
g = np.zeros(len(theta))
mi = np.zeros(len(theta))
af = np.zeros(len(theta))
psi = np.zeros(len(theta))

for i in range(len(theta)):
    mi[i] = (const * np.sin(np.radians(theta[i]))) / theta_3db
    g[i] = 20 * np.log10(((np.sin(np.radians(mi[i]))) / (mi[i])))  # uniform
    psi[i] = (2 * np.pi * (element_space / l_ambda) * (
        np.sin((np.radians(theta[i]))) - np.sin(np.radians(beamsteeringangle_az))))
    af[i] = np.sin(((number_elements * psi[i]) / 2)) / np.sin((psi[i] / 2))
    g[i] = g[i] + 10 * np.log10((1 / number_elements)) + 10 * np.log10(abs(af[i]) ** 2) + 20.2 - normalization

# Antenna pattern - Horizontal (Azimuth)
phi = np.linspace(-89, 89, 2000)
const = np.pi * 50.8
g1 = np.zeros(len(theta))
mi = np.zeros(len(theta))
af = np.zeros(len(theta))
psi = np.zeros(len(theta))

for i in range(len(theta)):
    mi[i] = (const * np.sin(np.radians(theta[i]))) / theta_3db
    g1[i] = 20 * np.log10(((np.sin(np.radians(mi[i]))) / (mi[i])))  # uniform
    psi[i] = (2 * np.pi * (element_space / l_ambda) * (np.sin((np.radians(theta[i]))) -
                                                       np.sin(np.radians(beamsteeringangle_el))))
    af[i] = np.sin(((number_elements * psi[i]) / 2)) / np.sin((psi[i] / 2))
    g1[i] = g1[i] + 10 * np.log10((1 / number_elements)) + 10 * np.log10(abs(af[i]) ** 2) + 20.2 - normalization

# Plot 3D Surface
x, y = np.meshgrid(np.linspace(-89, 89, 200), np.linspace(-89, 89, 200))
X, Y = np.meshgrid(g, g1, sparse=False)
fig = plt.figure(figsize=(10, 10))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, X+Y+g_max, cstride=1, rstride=1,  cmap='rainbow', antialiased=True, alpha=1)
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=20, ticks=np.linspace(-450, 46, 20))
ax.set_title('Phased Radar Antenna')
ax.set_xlabel('azimuth [degrees]')
ax.set_ylabel('elevation [degrees]')
plt.show()
