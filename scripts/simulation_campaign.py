# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 19:03:36 2018

@author: Calil
"""

import os
import glob

from sharc.main_cli import main

# Setup paths
cases_folder = 'cases'
subfolders = [f.path for f in os.scandir(cases_folder) if f.is_dir()] 

for folder in subfolders:
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
    main(['-p',file,'-o',folder])
