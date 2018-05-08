#!/bin/bash

echo "Begin siulations..."
python3 simulation_campaign.py
echo "Plot results..."
python3 plot_all.py

echo "Save tar..."
cd ..
tar -czvf scripts/simulation_resutls.tar.gz cases
cd scripts

echo "Notify you..."
python3 notify_me.py