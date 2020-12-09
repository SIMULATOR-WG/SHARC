"""
Module used to import output simulation files in format(.txt) to plot the results in one or multiples graphics
Author: Luciano Camilo Tue Nov 20 10:00:25 2020

"""

import pandas as pd
import numpy as np

data = pd.read_csv('[SYS] CDF of system INR.txt',skiprows=1, sep='	',header=None)
#data = pd.read_csv('[SYS] CDF of system INR2.txt',skiprows=1, sep='	',header=None)
#data1 = pd.read_csv('[IMT] CDF of IMT station antenna gain towards system.txt',skiprows=1, sep='	',header=None)
#data2 = pd.read_csv('[SYS] CDF of IMT to system path loss.txt',skiprows=1, sep='	',header=None)
#data3 = pd.read_csv('[IMT] CDF of IMT station antenna gain towards system.txt',skiprows=1, sep='	',header=None)
#data = pd.DataFrame(data)

import matplotlib.pyplot as plt
plt.figure(figsize=(8, 6))

x = data[0]
y = data[1]

#x1 = data1[0]
#y1 = data1[1]

#x2 = data2[0]
#y2 = data2[1]

#x3 = data3[0]
#y3 = data3[1]

x4 = np.linspace(-6,-6, 50)
y4 = np.linspace(0,1, 50)

#plt.plot(x, y, 'g^')
#plt.plot(x3, y3,'-.', linewidth=1, color = 'black', label= 'Limit')
#plt.ylabel('Probabilidade de INR < X', fontsize = 10, color ='black')
#plt.xlabel('INR [dB]', fontsize = 10, color ='black')
#plt.title('CDF of RAS INR')

#plt.subplot(2,2,1)
plt.plot(x, y,'r--', linewidth=1.0, color = 'blue', label= ' RAS (x=90000,y=0) / HIBS(x=0,y=0)')
plt.ylabel('Probability of INR < X', fontsize = 10, color ='black')
plt.xlabel('INR [dB]', fontsize = 10, color ='black')
plt.title('CDF of system INR')

#plt.subplot(2,2,2)
#plt.plot(x1, y1,'r--', linewidth=1.0, color = 'red', label = '90 degree elevation')
#plt.ylabel('Probabilidade de INR < X', fontsize = 10, color ='black')
#plt.xlabel('INR [dB]', fontsize = 10, color ='black')
#plt.title('[IMT] CDF of IMT station antenna gain towards system.')

#plt.subplot(2,2,3)
#plt.title('CDF-IMT To System Pathloss')
##plt.ylabel('Probabilidade de INR < X', fontsize = 10, color ='black')
#plt.xlabel('INR [dB]', fontsize = 10, color ='black')
#plt.plot(x2, y2,'r--', linewidth=1.0, color = 'orange', label= 'CDF - IMT to system Pathloss')

#plt.subplot(2,2,4)
#plt.figure(figsize=(16,22))

#plt.title('CDF of IMT Station antenna gain towards system')
#plt.ylabel('Probabilidade de do ganho da antena < X', fontsize = 10, color ='black')
#plt.xlabel('Gain [dBi]', fontsize = 10, color ='black')
#plt.plot(x3, y3,'r--', linewidth=1.0, color = 'orange', label= 'CDF - IMT Station Antenna gain towards system')

plt.legend(loc= 'upper left')
#plt.annotate('*(a)', xy=(-63.4, 0.79), xytext=(-63, 0.67),
#             arrowprops=dict(facecolor='black', shrink=0.005, linewidth=0.5),
#             )
plt.grid (which='minor', alpha=0.2)
plt.grid (which='major', alpha=0.5)
plt.grid(True, color='b', linestyle='--', linewidth=0.2)

plt.show()

