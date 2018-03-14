#!/bin/bash

echo "INSTALLING PYTHON, PIP AND GIT"
sudo apt-get update
sudo apt-get install python3=3.5.3*
sudo apt-get install python3-pip
python3 -m pip install --upgrade pip
sudo apt-get install git-core

echo "INSTALLING PYTHON PACKAGES"
python3 -m pip install numpy
python3 -m pip install shapely