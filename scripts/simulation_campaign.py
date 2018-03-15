# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:03:36 2018

@author: Calil
"""

print("_____SIMULATION SCRIPT_____")
print("Setting up..")

import os
import glob
from traceback import print_tb

from sharc.main_cli import main

# Setup paths
cases_folder = os.path.join('..','cases')
subfolders = [f.path for f in os.scandir(cases_folder) if f.is_dir()] 

print("Beginning simulation cases...")

for k, folder in enumerate(subfolders):
    case = os.path.basename(os.path.normpath(folder))
    print("\n\nCURRENT CASE: " + case)
    # Get ini file
    files = glob.glob(os.path.join(folder,'*.ini'))
    if len(files) == 0:
        print('\nWarning: no configuration file in case ' + case + 
              '. Going to next folder')
        continue
    elif len(files) > 1: 
        print('\nWarning: more than one configuration file in case ' + case + 
              '. Using file: ' + os.path.basename(os.path.normpath(files[0])))
    
    # Run simulation
    file = files[0]
    try:
        main(['-p',file,'-o',folder])
    except Exception as e:
        print(str(e) + "\nTraceback: ")
        print_tb(e.__traceback__)
        print("\n Moving on...")
        
    print("\n" + case + " case finished." + str(len(subfolders) - k - 1)\
          + " cases to go.\n")
    
print('DONE SIMULATING')
