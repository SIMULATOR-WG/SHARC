MSS-SS to IMT-DL Simulation – 2300 MHz Band

📄 Overview

This folder contains the simulation setup for evaluating interference from a Mobile Satellite Service - Single Station (MSS-SS) to an IMT downlink (IMT-DL) system operating in the 2300 MHz frequency band.

The scenario models a single MSS-SS satellite footprint located near the IMT coverage area. The simulation varies the distance between the MSS-SS footprint border and the IMT coverage edge, and computes the resulting Interference-to-Noise Ratio (INR).

⸻

🚀 Running the Simulation
	1.	Generate simulation parameters
Run the parameter generation script:

./scripts/parameter_gen.py


	2.	Start the simulation
	•	For parallel execution (multi-threaded):

./scripts/start_simulations_multi_thread.py


	•	For serial execution (single-threaded):

./scripts/start_simulations_single_thread.py



⸻

📊 Generating Results

After the simulations are complete, generate the result plots by running:

./scripts/plot_resutls.py

This will produce interactive Plotly HTML graphs summarizing the simulation outcomes.

