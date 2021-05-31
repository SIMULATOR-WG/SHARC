import numpy as np
from matplotlib import pyplot as plt

"""
Author: Agostinho Linhares de Souza Filho (linhares@anatel.gov.br)
        Luciano Camilo Alexandre (luciano-camilo@hotmail.com)

Date: Thu 29 April 09:29:25 2021

    Python script to calculate Total Radiated Efficiency of the Antenna Patterns.
    Antenna patterns used in this script: ITU-R RS.2043-0 Tables 8 and 9.
"""

# Antenna pattern points
theta_points = 18000
phi_points = 36000
# Scanning angles (theta/phi)
theta_angle = 18000
phi_angle = 36000
# Step size of antenna patterns
delta_theta = 180 / theta_points
delta_phi = 360 / phi_points
# Scanning angles vectors
theta = np.linspace(0, 180, theta_points)
phi = np.linspace(0, 360, phi_points)
# Antenna gain vectors (gv = vertical, gh = horizontal)
gv = np.zeros(theta_points)
gh = np.zeros(phi_points)

#################################################################################################################
"""
Antenna pattern according Recommendation ITU-R RS.2043-0 Table 8
"""
# Vertical antenna pattern (according theta angle vector)
for i in range(theta_angle):
    if theta[i] < 1.1:
        gv[i] = 47 - 9.91*((theta[i])**2)
    if 1.1 <= theta[i] <= 30:
        gv[i] = 35.9 - 0.83*theta[i]
    if theta[i] > 30:
        gv[i] = 11
# Horizontal antenna pattern (according phi angle vector)
for i in range(phi_angle):
    if phi[i] <= 0.5:
        gh[i] = 0 - 45.53*((phi[i])**2)
    if 0.5 < phi[i] <= 12:
        gh[i] = -10.97 - 2 * phi[i]
    if phi[i] > 12:
        gh[i] = -35

# Initialize gain vectors to transform gain from logaritmic to real scale
gv_ = np.zeros(theta_points)
gh_ = np.zeros(phi_points)

# Transform gain vectors from logaritmic to real scale

# Vertical antenna gain
for i in range(theta_points):
    gv_[i] = 10**(gv[i]/10)

# Horizontal antenna gain
for i in range(phi_points):
    gh_[i] = 10**(gh[i]/10)

# Calculate the Total Radiated Efficiency
eff_t = (np.sum(gv_*np.sin(theta*np.pi/180)+np.sin(theta*np.pi/180)*np.sum(gh_)
                * (delta_phi*np.pi/180))*(delta_theta*np.pi/180))/(np.pi*4)
print(f'Total efficiency ITU-R RS.2043-0 Table 8 = {eff_t*100:0.2f} %')

# Plot the antenna pattern according theta angle and phi angle defined in variables.

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))
ax1.fill_between(theta, gv, color='black', alpha=0.70)
csfont = {'fontname': 'Times New Roman'}
hfont = {'fontname': 'Times New Roman'}
ax1.grid()
ax1.set_title('ITU-R RS.2043-0 - Vertical Antenna Pattern')
ax1.set_xlabel('Elevation angle (degrees)', fontsize=12, color='black')
ax1.set_ylabel('Gain (dBi)', fontsize=12, color='black')
ax1.semilogx(theta, gv, "-", color='black', label="ITU-R RS.2043-0 (Table 8)")
ax1.legend()
ax1.text(0.01, 2.7, f'ITU-R RS.2043-0 (Table 8) Efficiency = ' '{0:3.2f}' '%' .format(eff_t*100, 2))

ax2.fill(phi, gh, color='black', alpha=0.70)
ax2.set_title('ITU-R RS.2043-0 - Horizontal Antenna Pattern')
ax2.set_xlabel('Azimuth (degrees)', fontsize=12, color='black')
ax2.set_ylabel('Gain (dBi)', fontsize=12, color='black')
ax2.semilogx(phi, gh, color='black', label="ITU-R RS.2043-0 (Table 8)")
ax2.legend(loc='upper right')
ax2.grid()

"""
Antenna pattern according Recommendation ITU-R RS.2043-0 Table 9

"""
# Vertical antenna pattern (according theta angle vector)
for i in range(theta_angle):
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

# Horizontal antenna pattern (according phi angle vector)
for i in range(phi_angle):
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

# Initialize gain vectors to transform gain from logaritmic to real scale
gv_ = np.zeros(theta_points)
gh_ = np.zeros(phi_points)

# Transform gain vectors from logaritmic to real scale

# Vertical antenna gain vector
for i in range(theta_points):
    gv_[i] = 10**(gv[i]/10)

# Horizontal antenna gain vector
for i in range(phi_points):
    gh_[i] = 10**(gh[i]/10)

# Inititialize Efficiency vector
eff_t = np.zeros(phi_points)

# Calculate the Total Radiated Efficiency
eff_t = (np.sum(gv_*np.sin(theta*np.pi/180)+np.sin(theta*np.pi/180)*np.sum(gh_)
                * delta_phi*np.pi/180)*(delta_theta*np.pi/180))/(np.pi*4)

print(f'Total efficiency ITU-R RS.2043-0 Table 9 = {eff_t*100:0.2f} %')

# Plot the antenna pattern according theta angle and phi angle defined in variables.
ax1.fill_between(theta, gv, color='darkorange', alpha=0.7)
csfont = {'fontname': 'Times New Roman'}
hfont = {'fontname': 'Times New Roman'}
ax1.grid()
ax1.semilogx(theta, gv, "-", color='darkorange', label='ITU-R RS.2043-0 (Table 9)')
ax1.text(0.01, 0.5, f'ITU-R RS.2043-0 (Table 9) Efficiency = ''{0:3.2f}' '%' .format(eff_t*100, 2))
ax1.legend()
ax2.fill(phi, gh, color='darkorange', alpha=0.7)
ax2.semilogx(phi, gh, "-",  color='darkorange', label="ITU-R RS.2043-0 (Table 9)")
ax2.legend(loc='upper right')
ax2.grid()
plt.show()
