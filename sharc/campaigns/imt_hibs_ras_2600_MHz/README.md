
# Spectrum Sharing Study: IMT HIBS vs. Radio Astronomy (2.6 GHz)
This directory holds the code and data for a simulation study investigating spectrum sharing between:

IMT HIBS as a (NTN) system and Radio Astronomy Service (RAS) operating in the 2.6 GHz band.

Each campaing puts the RAS station farther from the IMT BS nadir point over the Earth's surface.

Main campaign parameters:
- IMT topolgy: NTN with IMT BS at 20km of altitude
- IMT @2680MHz/20MHz BW
- RAS @2695/10MHz BW
- Channel model: P.619 for both IMT and IMT-RAS links.

# Folder Structure
inputs: This folder contains parameter files used to configure the simulation campaigns. Each file defines specific scenarios for the NTN and RAS systems.
scripts: This folder holds post-processing and plotting scripts used to analyze the simulation data. These scripts generate performance metrics and visualizations based on the simulation outputs.

# Dependencies
This project may require additional software or libraries to run the simulations and post-processing scripts. Please refer to the individual script files for specific dependencies.

# Running the Simulations
`python3 main_cli.py -p campaigns/imt_hibs_ras_2600_MHz/input/<campaign-parameter-file>`
or on root
`python3 sharc/main_cli.py -p sharc/campaigns/imt_hibs_ras_2600_MHz/input/<campaign-parameter-file>`

# Running the scripts

## For plotting 

`python3 sharc/campaigns/imt_hibs_ras_2600_MHz/scripts/plot_results.py`

## For starting simulation multi thread

`sharc/campaigns/imt_hibs_ras_2600_MHz/scripts/start_simulations_multi_thread.py`

## For starting simulation single thread

`sharc/campaigns/imt_hibs_ras_2600_MHz/scripts/start_simulations_single_thread.py`