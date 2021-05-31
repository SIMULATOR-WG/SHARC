"""
@Created: Luciano Camilo on Tue Fr 04 18:27:25 2021

"""
import numpy as np
from matplotlib import pyplot as plt

pontos_theta = 250
pontos_phi = 250

theta = np.linspace(0, 180, pontos_theta)
phi = np.linspace(0, 360, pontos_phi)
gv = np.zeros(pontos_theta)
gh = np.zeros(pontos_phi)
gv2 = np.zeros(200)

"""
Table 9 - Table 8 - ITU-R RS.2043-0
"""

for i in range(len(theta)):
    if theta[i] < 1.149:
        gv[i] = 47 - 9.91*((theta[i])**2)
    if 1.149 <= theta[i] <= 9.587:
        gv[i] = 35.189 - 1.9440*theta[i]
    if 9.587 <= theta[i] <= 29.976:
        gv[i] = 21.043 - 0.4680 * theta[i]
    if 29.976 <= theta[i] <= 50:
        gv[i] = 12.562 - 0.1850 * theta[i]
    if theta[i] > 50:
        gv[i] = 3.291


gv_flip = gv[::-1]
gv2 = np.append(gv_flip, gv)

for i in range(len(phi)):
    if phi[i] <= 0.542:
        gh[i] = 0 - 45.53*((phi[i])**2)
    if 0.542 < phi[i] <= 5.053:
        gh[i] = -11.210 - 4.0220 * phi[i]
    if 5.053 < phi[i] <= 14.708:
        gh[i] = -26.720 - 0.9530 * phi[i]
    if 14.708 < phi[i] <= 30:
        gh[i] = -35.031 - 0.3880 * phi[i]
    if 30 < phi[i] <= 59.915:
        gh[i] = -41.836 - 0.1580 * phi[i]
    if phi[i] > 59.915:
        gh[i] = -51.387

gh_flip = gh[::-1]
gh2 = np.append(gh_flip, gh)

csfont = {'fontname': 'Times New Roman'}
hfont = {'fontname': 'Times New Roman'}
plt.grid()
plt.title('ITU-R RS.2043-0 - Vertical Pattern')
plt.xlabel('Elevation angle (degrees)', fontsize=12, color='black')
plt.ylabel('Gain (dBi)', fontsize=12, color='black')
theta = np.linspace(-90, 90, 500)
plt.plot(theta, gv2, "-", color='darkorange', label="ITU-R RS.2043-0 (Table 9)")
plt.legend()
plt.show()
plt.title('ITU-R RS.2043-0 (Table 9) - Horizontal Antenna Pattern')
plt.xlabel('Azimuth (degrees)', fontsize=12, color='black')
plt.ylabel('Gain (dBi)', fontsize=12, color='black')
phi = np.linspace(-90, 90, 500)
plt.plot(phi, gh2, "--b", label="")
plt.grid()
x, y = np.meshgrid(np.linspace(-90, 90, 500), np.linspace(-90, 90, 500))
X, Y = np.meshgrid(gv2, gh2, sparse=False)
fig = plt.figure(figsize=(10, 10))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x, y, X+Y, cstride=1, rstride=1,  cmap='rainbow', antialiased=True, alpha=1)
fig.colorbar(surf, ax=ax, shrink=1, aspect=10, ticks=np.linspace(-50, 60, 20))
ax.set_title('3D Antenna Pattern ITU-R RS.2043 - Table 9')
ax.set_xlabel('azimuth [degrees]')
ax.set_ylabel('elevation [degrees]')
plt.show()
