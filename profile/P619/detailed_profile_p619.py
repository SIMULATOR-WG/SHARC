# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 10:00:00 2024

Profile script for PropagationP619
"""

import os
import csv
import numpy as np
import matplotlib
from sharc.propagation.propagation_p619 import PropagationP619
import cProfile
import pstats
import io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

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


def profile_get_atmospheric_gasses_loss():
    losses = []
    for elevation in apparent_elevation:
        loss = propagation._get_atmospheric_gasses_loss(frequency_MHz=frequency_MHz,
                                                        apparent_elevation=elevation,
                                                        surf_water_vapour_density=surf_water_vapour_density,
                                                        lookupTable=False)
        losses.append(loss)
    return losses


# Profile the function
pr = cProfile.Profile()
pr.enable()
profile_get_atmospheric_gasses_loss()
pr.disable()

# Create profile results directory
output_dir = os.path.join(os.path.dirname(__file__), 'profile_results')
os.makedirs(output_dir, exist_ok=True)

# Save profiling results to CSV
profile_data = io.StringIO()
ps = pstats.Stats(pr, stream=profile_data).sort_stats(pstats.SortKey.CUMULATIVE)
ps.print_stats()

profile_data.seek(0)
lines = profile_data.readlines()

output_file = os.path.join(output_dir, 'profiling_results.csv')
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Function', 'Calls', 'Total Time', 'Per Call', 'Cumulative Time', 'Per Call (Cum)'])
    for line in lines[1:]:
        if line.strip() == "" or "percall" in line:
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        func_name = " ".join(parts[5:])
        writer.writerow([func_name, parts[0], parts[2], parts[3], parts[4], parts[1]])

print(f"Profiling results saved to {output_file}")

# Load profiling results from CSV
functions = []
cumulative_times = []

with open(output_file, mode='r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip header
    for row in reader:
        try:
            cumulative_times.append(float(row[4]))
            functions.append(row[0])
        except ValueError:
            continue

# Sort the data
sorted_indices = sorted(range(len(cumulative_times)), key=lambda k: cumulative_times[k], reverse=True)
sorted_functions = [functions[i] for i in sorted_indices][:5]
sorted_cumulative_times = [cumulative_times[i] for i in sorted_indices][:5]

# Plot profiling results
plt.figure(figsize=(14, 8))
plt.barh(sorted_functions, sorted_cumulative_times, color='blue')
plt.xlabel("Cumulative Time (s)")
plt.ylabel("Function")
plt.title("Top 5 Functions by Cumulative Time")
plt.grid(True)
plt.tight_layout()  # Adjust layout to fit labels
plot_file = os.path.join(output_dir, 'profiling_results_bar.png')
plt.savefig(plot_file)
print(f"Profiling results bar plot saved to {plot_file}")
