# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 12:59:12 2017

"""

import os
import csv
import numpy as np
from sharc.propagation.propagation_p619 import PropagationP619

# Constants
frequency_MHz = 2680.0
space_station_alt_m = 35786.0 * 1000  # Geostationary orbit altitude in meters
earth_station_alt_m = 1000.0
earth_station_lat_deg = -15.7801
earth_station_long_diff_deg = 0.0  # Not used in the loss calculation
season = "SUMMER"
apparent_elevation = np.arange(0, 90)  # Elevation angles from 0 to 90 degrees
city_name = "BRASILIA"

# Initialize the propagation model
random_number_gen = np.random.RandomState(101)
propagation = PropagationP619(
    random_number_gen=random_number_gen,
    space_station_alt_m=space_station_alt_m,
    earth_station_alt_m=earth_station_alt_m,
    earth_station_lat_deg=earth_station_lat_deg,
    earth_station_long_diff_deg=earth_station_long_diff_deg,
    season=season,
)

# Calculate the loss for each elevation angle
losses = []
for elevation in apparent_elevation:
    loss = propagation._get_atmospheric_gasses_loss(
        frequency_MHz=frequency_MHz,
        apparent_elevation=elevation,
    )
    losses.append(loss)

# Save results to CSV file
output_dir = os.path.join(os.path.dirname(__file__), 'BRASILIA')
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(
    output_dir, f'{city_name}_{int(frequency_MHz)}_{int(earth_station_alt_m)}m.csv',
)

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['apparent_elevation', 'loss'])
    for elevation, loss in zip(apparent_elevation, losses):
        writer.writerow([elevation, loss])

print(f"Results saved to {output_file}")
