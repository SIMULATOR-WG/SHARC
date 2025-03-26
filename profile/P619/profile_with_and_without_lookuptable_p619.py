# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 10:00:00 2024

Profile script for PropagationP619
"""

import os
import time
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from sharc.propagation.propagation_p619 import PropagationP619  # noqa: E402

# Constants
frequency_MHz = 2680.0
space_station_alt_m = 35786.0 * 1000  # Geostationary orbit altitude in meters
earth_station_alt_m = 1000.0
earth_station_lat_deg = -15.7801
earth_station_long_diff_deg = 0.0  # Not used in the loss calculation
season = "SUMMER"
apparent_elevation = np.arange(0, 90)  # Elevation angles from 0 to 90 degrees
surf_water_vapour_density = 2.5

# Initialize the propagation model
random_number_gen = np.random.RandomState(101)
propagation = PropagationP619(random_number_gen=random_number_gen,
                              space_station_alt_m=space_station_alt_m,
                              earth_station_alt_m=earth_station_alt_m,
                              earth_station_lat_deg=earth_station_lat_deg,
                              earth_station_long_diff_deg=earth_station_long_diff_deg,
                              season=season)

# Profile with lookup table
start_time = time.time()
losses_with_table = []
for elevation in apparent_elevation:
    loss = propagation._get_atmospheric_gasses_loss(frequency_MHz=frequency_MHz,
                                                    apparent_elevation=elevation,
                                                    surf_water_vapour_density=surf_water_vapour_density)
    losses_with_table.append(loss)
time_with_table = time.time() - start_time

# Profile without lookup table
start_time = time.time()
losses_without_table = []
for elevation in apparent_elevation:
    loss = propagation._get_atmospheric_gasses_loss(frequency_MHz=frequency_MHz,
                                                    apparent_elevation=elevation,
                                                    surf_water_vapour_density=surf_water_vapour_density,
                                                    lookupTable=False)
    losses_without_table.append(loss)
time_without_table = time.time() - start_time

# Save profiling results to CSV
output_dir = os.path.join(os.path.dirname(__file__), 'profile_results')
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, 'profiling_results_table.csv')

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Method', 'Time (s)'])
    writer.writerow(['With Lookup Table', time_with_table])
    writer.writerow(['Without Lookup Table', time_without_table])

print(f"Profiling results saved to {output_file}")

# Save atmospheric loss plots
plt.figure()
plt.plot(apparent_elevation, losses_with_table, label='With Lookup Table')
plt.plot(apparent_elevation, losses_without_table, label='Without Lookup Table')
plt.xlabel("Apparent Elevation (deg)")
plt.ylabel("Loss (dB)")
plt.title("Atmospheric Gasses Attenuation")
plt.legend()
plt.grid(True)
plot_file = os.path.join(output_dir, 'atmospheric_loss_comparison.png')
plt.savefig(plot_file)
print(f"Atmospheric loss plot saved to {plot_file}")


import matplotlib  # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

# Load profiling results from CSV
input_file = os.path.join(os.path.dirname(__file__), 'profile_results', 'profiling_results_table.csv')
methods = []
times = []

with open(input_file, mode='r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip header
    for row in reader:
        methods.append(row[0])
        times.append(float(row[1]))

# Plot profiling results
plt.figure()
plt.bar(methods, times, color=['blue', 'orange'])
plt.xlabel("Method")
plt.ylabel("Time (s)")
plt.title("Profiling Results for PropagationP619")
plt.grid(True)
plot_file = os.path.join(os.path.dirname(__file__), 'profile_results', 'profiling_results_bar_table.png')
plt.savefig(plot_file)
print(f"Profiling results bar plot saved to {plot_file}")
