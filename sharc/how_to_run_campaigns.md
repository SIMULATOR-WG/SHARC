

## Setup the campaigns folder

1. Create a folder for your campaign:
   - Example: `imt_hibs_ras_2600_MHz`

2. Inside the campaign folder, create the following subfolders:
   - `input` -> This folder is for defining the simulation parameters. For example, look in `campaigns/imt_hibs_ras_2600_MHz/input`.
   - `output` -> SHARC will use this folder to save the results and plots.
   - `scripts` -> This folder is used to implement scripts for post-processing.

## Configure your simulation

1. Create a parameter file:
   - Example: `campaigns/imt_hibs_ras_2600_MHz/input/parameters_hibs_ras_2600_MHz_0km.ini`.

2. Set the configuration for your study in the parameter file.

3. Set the output folder in the parameter file:
    ```ini
    ###########################################################################
    # output destination folder - this is relative to the SHARC/sharc directory
    output_dir = campaigns/imt_hibs_ras_2600_MHz/output/
    ###########################################################################
    # output folder prefix
    output_dir_prefix = output_imt_hibs_ras_2600_MHz_0km
    ```
4. You can create multiple simulation parameters:
   - Check the folder `sharc/campaigns/imt_hibs_ras_2600_MHz/input` for examples.

## Run simulations and view the results

1. **Run simulations:**
   - In the `scripts` folder, you can create a file to start the simulation. Check the example: `campaigns/imt_hibs_ras_2600_MHz/scripts/start_simulations_multi_thread.py`.

   - If you want to run a single-threaded simulation, check the file: `campaigns/imt_hibs_ras_2600_MHz/scripts/start_simulations_single_thread.py`.

2. **Generate plots:**
   - You can create a file to read the data and generate the plots. SHARC has a function called `plot_cdf` to make plotting easy. Check the example: `campaigns/imt_hibs_ras_2600_MHz/scripts/plot_results.py`.
"""

