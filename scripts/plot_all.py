"""
Created on Wed Mar 14 19:03:36 2018

@author: Calil
"""

print("_____PLOT SCRIPT_____")
print("Setting up..")

from numpy import loadtxt
import matplotlib.pyplot as plt
from glob import glob
import os

cases_folder = os.path.join('..','cases')
subfolders = [f.path for f in os.scandir(cases_folder) if f.is_dir()] 

plt.figure()

print("Plotting all...")

for folder in subfolders:
    
    case = os.path.basename(os.path.normpath(folder))
    print("\n\nPloting case " + case)
    
    # Search all files with txt extension in folder
    files = glob(os.path.join(folder,'*.txt'))
    
    for file in files:
        case = os.path.basename(os.path.normpath(file))
        print("Ploting file \"" + case + "\"")
        # Collect and plot data
        data = loadtxt(file,skiprows=1)
        plt.plot(data[:,0],data[:,1])
        # Save figure
        plot_name = os.path.basename(os.path.normpath(file))
        plot_name = plot_name[:-4]
        plt.title(plot_name)
        
        # Save file
        file_name = (plot_name + ".png")
        file_name = file_name.replace(' ','_')
        file_name = file_name.replace('[','')
        file_name = file_name.replace(']','')
        file_path = os.path.join(folder,file_name)
        plt.savefig(file_path)
        plt.clf()

print('DONE PLOTTING')