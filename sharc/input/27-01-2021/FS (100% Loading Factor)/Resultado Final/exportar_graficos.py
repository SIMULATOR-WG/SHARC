"""
Module used to import output simulation files in format(.txt) to plot the results in one or multiples graphics
Author: Luciano Camilo Tue Nov 20 10:00:25 2020

"""

import pandas as pd
import numpy as np
from matplotlib.ticker import FuncFormatter
import math

#data = pd.read_csv('[SYS] CDF of system INR2.txt',skiprows=1, sep='	',header=None)
data = pd.read_csv('[SYS] CDF of system INR_0.txt',skiprows=1, sep='	',header=None)
data1 = pd.read_csv('[SYS] CDF of system INR_50000.txt',skiprows=1, sep='	',header=None)
data2 = pd.read_csv('[SYS] CDF of system INR_100000.txt',skiprows=1, sep='	',header=None)
data3= pd.read_csv('[SYS] CDF of system INR_150000.txt',skiprows=1, sep='	',header=None)
data4= pd.read_csv('[SYS] CDF of system INR_300000.txt',skiprows=1, sep='	',header=None)
data5= pd.read_csv('[SYS] CDF of system INR_500000.txt',skiprows=1, sep='	',header=None)
#data6= pd.read_csv('[SYS] CDF of system INR_450000.txt',skiprows=1, sep='	',header=None)
#data7= pd.read_csv('[SYS] CDF of system interference power from IMT DL-7000000.txt',skiprows=1, sep='	',header=None)
#data8= pd.read_csv('[SYS] CDF of system interference power from IMT DL-10000000.txt',skiprows=1, sep='	',header=None)

import matplotlib.pyplot as plt
plt.figure(figsize=(8, 6))

x = data[0]
y = data[1]

x1 = data1[0]
y1 = data1[1]

x2 = data2[0]
y2 = data2[1]

x3 = data3[0]
y3 = data3[1]

x4 = data4[0]
y4 = data4[1]

x5 = data5[0]
y5 = data5[1]

#x6 = data6[0]
#y6 = data6[1]

#x7 = data7[0]-10*np.log10(1000)
#y7 = data7[1]

#x8= data8[0]-10*np.log10(1000)
#y8 = data8[1]

x_limite = np.linspace(-6,-6, 50)
y_limite = np.linspace(0,1, 50)


plt.subplot(2,2,1)
#plt.xlim(-250, -180)
#plt.ylim(0, 1)
#plt.legend(prop={'family': 'Times New Roman'})
#plt.xscale('linear')
#plt.yscale('log')
#plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
plt.plot(x, y,'r-', linewidth=1.5, color = 'brown', label= ' 0 km ')


###################################################################################################################


#plt.ylim(0, 1)
#plt.legend(prop={'family': 'Times New Roman'})
#plt.xscale('linear')
#plt.yscale('log')
#plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
plt.plot(x1, y1,'r-', linewidth=1.5, color = 'orange', label= ' 50 km ')


####################################################################################################################


#plt.ylim(0, 1)
#plt.legend(prop={'family': 'Times New Roman'})
#plt.xscale('linear')
#plt.yscale('log')
#plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
plt.plot(x2, y2,'r-', linewidth=1.5, color = 'green', label= ' 100 km ')


#######################################################################################################################

#plt.ylim(0, 1)
#plt.legend(prop={'family': 'Times New Roman'})
#plt.xscale('linear')
#plt.yscale('log')
#plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
plt.plot(x3, y3,'r-', linewidth=1.5, color = 'red', label= ' 150 km ')


#######################################################################################################################

#plt.ylim(0, 1)
#plt.legend(prop={'family': 'Times New Roman'})
#plt.xscale('linear')
#plt.yscale('log')
#plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
plt.plot(x4, y4,'r-', linewidth=1.5, color = 'blue', label= ' 300 km')



#######################################################################################################################


#plt.ylim(0, 1)
#plt.legend(prop={'family': 'Times New Roman'})
#plt.xscale('linear')
#plt.yscale('log')
#plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
plt.plot(x5, y5,'r-', linewidth=1.5, color = 'cyan', label= ' 500 km')


#######################################################################################################################

#plt.ylim(0, 1)
#plt.legend(prop={'family': 'Times New Roman'})
#plt.xscale('linear')
#plt.yscale('log')
#plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
#plt.plot(x6, y6,'r-', linewidth=1.5, color = 'yellow', label= ' 450 km ')


#######################################################################################################################

#plt.ylim(0, 1)
#plt.legend(prop={'family': 'Times New Roman'})
#plt.xscale('linear')
#plt.yscale('log')
#plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
#plt.plot(x7, y7,'r-', linewidth=1.5, color = 'red', label= ' 7000 km ')


#######################################################################################################################

#plt.ylim(0, 1)
#plt.legend(prop={'family': 'Times New Roman'})
#plt.xscale('linear')
#plt.yscale('log')
#plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
#plt.plot(x8, y8,'r-', linewidth=1.5, color = 'black', label= ' 10000 km ')


#######################################################################################################################
#plt.ylim(0, 1)
#plt.legend(prop={'family': 'Times New Roman'})
#plt.xscale('linear')
#plt.yscale('log')
#plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
plt.plot(x_limite, y_limite,'r--', linewidth=1.5, color = 'black', label= ' Protection Criteria ')
#plt.title('CDF of Interference Power')
plt.ylabel('Probability of INR < X', fontsize = 10, color ='black')
plt.xlabel('INR [dB]', fontsize = 10, color ='black')

plt.legend(loc= 'upper right', bbox_to_anchor=(1.33,1),
          fancybox=True, shadow=True, ncol=1)
plt.grid (which='minor', alpha=0.5)
plt.grid (which='major', alpha=0.5)
plt.grid(True, color='b', linestyle='--', linewidth=0.2)






















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

#plt.legend(loc= 'upper left')
#plt.annotate('*(a)', xy=(-63.4, 0.79), xytext=(-63, 0.67),
#             arrowprops=dict(facecolor='black', shrink=0.005, linewidth=0.5),
#             )
#plt.grid (which='minor', alpha=0.2)
#plt.grid (which='major', alpha=0.5)
#plt.grid(True, color='b', linestyle='--', linewidth=0.2)

plt.show()

