"""
@Created: Luciano Camilo on Tue Fr 04 18:27:25 2021

"""
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(6, 5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height])

"""
Calculate the vertical antenna pattern (elevation) of Phased array radar antenna.

    Construction:
            frequency: frequency of operation (MHz)
            element_space : value * lambda
            beamsteeringangle_el: elevation beamsteering angle
            beamsteeringangle_az: azimuth beamsteering angle
            number_elements = numer of antenna elements
            theta_3db : 3dB theta beamwidth (fixed)
            phi_3db : 3dB phi beamwidth (fixed)
            g_max = maximum gain

"""
theta_3db = 1.1
phi_3db = 1.1
frequency = 2700
l_ambda = 299792458 / (frequency * 1e6)
element_space = l_ambda*0.5
beamsteeringangle_el = 0
beamsteeringangle_az = 0
g_max = 46

if g_max == 41:
    normalization = 0.833
    number_elements = 38
if g_max == 46:
    normalization = 3.48681
    number_elements = 70

theta = np.linspace(-90, 90, 200)
const = np.pi * 50.8
gv = np.zeros(len(theta))
mi = np.zeros(len(theta))
af = np.zeros(len(theta))
psi = np.zeros(len(theta))

for i in range(len(theta)):
    mi[i] = (const * np.sin(np.radians(theta[i]))) / theta_3db
    gv[i] = 20 * np.log10(((np.sin(np.radians(mi[i]))) / (mi[i])))  # uniform
    psi[i] = (2 * np.pi * (element_space / l_ambda) * (
        np.sin((np.radians(theta[i]))) - np.sin(np.radians(beamsteeringangle_az))))
    af[i] = np.sin(((number_elements * psi[i]) / 2)) / np.sin((psi[i] / 2))
    gv[i] = gv[i] + 10 * np.log10((1 / number_elements)) + 10 * np.log10(abs(af[i]) ** 2) + 20.2 - normalization

"""
Calculate the horizontal antenna pattern (azimuth) of Phased array radar antenna.

    Construction:
            frequency: frequency of operation (MHz)
            element_space : value * lambda
            beamsteeringangle_el: elevation beamsteering angle
            beamsteeringangle_az: azimuth beamsteering angle
            number_elements = numer of antenna elements
            theta_3db : 3dB theta beamwidth
            phi_3db : 3dB phi beamwidth
            g_max = maximum gain

"""
phi = np.linspace(-90, 90, 200)
const = np.pi * 50.8
gh = np.zeros(len(theta))
mi = np.zeros(len(theta))
af = np.zeros(len(theta))
psi = np.zeros(len(theta))

for i in range(len(phi)):
    mi[i] = (const * np.sin(np.radians(phi[i]))) / phi_3db
    gh[i] = 20 * np.log10(((np.sin(np.radians(mi[i]))) / (mi[i])))  # uniform
    psi[i] = (2 * np.pi * (element_space / l_ambda) * (
        np.sin((np.radians(phi[i]))) - np.sin(np.radians(beamsteeringangle_el))))
    af[i] = np.sin(((number_elements * psi[i]) / 2)) / np.sin((psi[i] / 2))
    gh[i] = gh[i] + 10 * np.log10((1 / number_elements)) + 10 * np.log10(abs(af[i]) ** 2) + 20.2 - normalization

# Plot 3D Contour Surface
x, y = np.meshgrid(np.linspace(-90, 90, 200), np.linspace(-90, 90, 200))
X, Y = np.meshgrid(gh, gv, sparse=False)
clevs1 = np.arange(-150, 46, 0.5)
cs1 = plt.contourf(y, x, X+Y+g_max, clevs1, cmap='viridis', alpha=1)
plt.colorbar(cs1)
ax.set_title('Radar Phased Array Antenna')
ax.set_xlabel('azimuth [degrees]')
ax.set_ylabel('elevation [degrees]')
plt.show()
