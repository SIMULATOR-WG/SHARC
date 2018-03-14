#!/bin/bash

echo "INSTALLING PYTHON, PIP AND GIT"
sudo apt-get update
sudo apt-get install python3=3.5.3*
sudo apt-get install python3-pip
sudo python3 -m pip install --upgrade pip
sudo apt-get install python-support
sudo apt-get install python3-tk
sudo apt-get install git-core

echo "INSTALLING PYTHON PACKAGES"
sudo python3 -m pip install numpy
sudo python3 -m pip install shapely
sudo python3 -m pip install scipy
sudo python3 -m pip install matplotlib

echo "CLONING REPOSITORY"
git clone https://github.com/MWSL-UnB/FD-SHARC.git
git pull origin dev_simulation_scripts
