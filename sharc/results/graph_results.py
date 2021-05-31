"""
Module used to import output simulation files in format(.txt) to plot the results in one or multiples graphics
Author: Luciano Camilo Tue Nov 20 10:00:25 2020

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import math
"""
Copy the result file in output folder from SHARC to this folder
"[SYS] CDF of system INR.txt" = replace by the results that you want to plot

Any questions please contact Luciano Camilo, luciano-camilo@hotmail.com or lcamilo@tmgtelecom.com
"""
data = pd.read_csv('[SYS] CDF of system INR.txt', skiprows=1, sep='	', header=None)
# data = pd.read_csv('file path',skiprows=1, sep='	', header=None)
# data1 = pd.read_csv('file path',skiprows=1, sep='	', header=None)
# data2 = pd.read_csv('file path',skiprows=1, sep='	', header=None)
# data3= pd.read_csv('file path',skiprows=1, sep='	', header=None)
# data4= pd.read_csv('file path',skiprows=1, sep='	', header=None)
# data5= pd.read_csv('file path',skiprows=1, sep='	', header=None)
# data6= pd.read_csv('file path',skiprows=1, sep='	', header=None)
# data7= pd.read_csv('file path',skiprows=1, sep='	', header=None)
# data8= pd.read_csv('file path',skiprows=1, sep='	', header=None)

plt.figure(figsize=(8, 6))

x = data[0]
y = data[1]

# x1 = data1[0]
# y1 = data1[1]

# x2 = data2[0]
# y2 = data2[1]

# x3 = data3[0]
# y3 = data3[1]

# x4 = data4[0]
# y4 = data4[1]

# x5 = data5[0]
# y5 = data5[1]

# x6 = data6[0]
# y6 = data6[1]

# x7 = data7[0]
# y7 = data7[1]

# x8= data8[0]
# y8 = data8[1]

# Put in variable xlim, the value of system protection criteria
xlim = -10
x_limite = np.linspace(-xlim, -xlim, 50)
y_limite = np.linspace(0, 1, 50)


plt.subplot(2, 2, 1)
# plt.xlim(-250, -180)
plt.ylim(0.9951, 1)
# plt.legend(prop={'family': 'Times New Roman'})
# plt.xscale('linear')
plt.yscale('log')
# plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
plt.semilogy(x, y, 'r-', linewidth=1.5, color='brown', label=' 0 km ')


###################################################################################################################


# plt.ylim(0, 1)
# plt.legend(prop={'family': 'Times New Roman'})
# plt.xscale('linear')
# plt.yscale('log')
# plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
# plt.plot(x1, y1,'r-', linewidth=1.5, color = 'orange', label= ' 50 km ')


####################################################################################################################


# plt.ylim(0, 1)
# plt.legend(prop={'family': 'Times New Roman'})
# plt.xscale('linear')
# plt.yscale('log')
# plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
# plt.plot(x2, y2,'r-', linewidth=1.5, color = 'green', label= ' 100 km ')


#######################################################################################################################

# plt.ylim(0, 1)
# plt.legend(prop={'family': 'Times New Roman'})
# plt.xscale('linear')
# plt.yscale('log')
# plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
# plt.plot(x3, y3,'r-', linewidth=1.5, color = 'red', label= ' 120 km ')


#######################################################################################################################

# plt.ylim(0, 1)
# plt.legend(prop={'family': 'Times New Roman'})
# plt.xscale('linear')
# plt.yscale('log')
# plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
# plt.plot(x4, y4,'r-', linewidth=1.5, color = 'blue', label= ' 250 km ')


#######################################################################################################################


# plt.ylim(0, 1)
# plt.legend(prop={'family': 'Times New Roman'})
# plt.xscale('linear')
# plt.yscale('log')
# plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
# plt.plot(x5, y5,'r-', linewidth=1.5, color = 'cyan', label= ' 500 km ')


#######################################################################################################################

# plt.ylim(0, 1)
# plt.legend(prop={'family': 'Times New Roman'})
# plt.xscale('linear')
# plt.yscale('log')
# plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
# plt.plot(x6, y6,'r-', linewidth=1.5, color = 'yellow', label= ' 450 km ')


#######################################################################################################################

# plt.ylim(0, 1)
# plt.legend(prop={'family': 'Times New Roman'})
# plt.xscale('linear')
# plt.yscale('log')
# plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
# plt.plot(x7, y7,'r-', linewidth=1.5, color = 'red', label= ' 7000 km ')


#######################################################################################################################

# plt.ylim(0, 1)
# plt.legend(prop={'family': 'Times New Roman'})
#plt.xscale('linear')
#plt.yscale('log')
# plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
# plt.plot(x8, y8,'r-', linewidth=1.5, color = 'black', label= ' 10000 km ')


#######################################################################################################################
# plt.ylim(0, 1)
# plt.legend(prop={'family': 'Times New Roman'})
# plt.xscale('linear')
# plt.yscale('log')
# plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{}'.format(y)))
plt.plot(x_limite, y_limite, 'r--', linewidth=1.5, color='black', label=' Protection Criteria ')
# plt.title('CDF of Interference Power')
plt.ylabel('Probability of INR < X', fontsize=10, color='black')
plt.xlabel('INR [dB]', fontsize=10, color='black')

plt.legend(loc='upper right', bbox_to_anchor=(1.33, 1), fancybox=True, shadow=True, ncol=1)
plt.grid(which='minor', alpha=0.5)
plt.grid(which='major', alpha=0.5)
plt.grid(True, color='b', linestyle='--', linewidth=0.2)
plt.show()
