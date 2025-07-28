from sharc.run_multiple_campaigns_mut_thread import run_campaign_re
import argparse
import subprocess

# Set the campaign name
# The name of the campaign to run. This should match the name of the campaign directory.
name_campaign = "mss_d2d_to_imt"

# Run the campaign in single-thread mode
# This function will execute the campaign with the given name in a single-threaded manner.
# It will look for the campaign directory under the specified name and start the necessary processes.
# run_campaign_re(name_campaign, r'^parameters_mss_d2d_to_imt_co_channel_system_A.yaml')
# run_campaign_re(name_campaign, r'^parameters_mss_d2d_to_imt_(dl,ul)_co_channel_system_A.yaml')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run SHARC MSS-D2D simulations.")
    parser.add_argument(
        "--scenario",
        type=int,
        choices=range(0, 2),
        required=True,
        help="""Specify the campaign scenario as a number (0 or 1).
                0 - MSS-D2D to IMT-UL/DL
                1 - MSS-D2D to IMT-UL/DL with varying latitude."""
    )

    # Parse the arguments
    args = parser.parse_args()

    # Update the campaign regex with the provided scenario
    scenario = args.scenario
    print(f"Running scenario {scenario}...")
    if scenario == 0:
        regex_pattern = r'^parameters_mss_d2d_to_imt_(dl|ul)_co_channel_system_A.yaml'
    elif scenario == 1:
        print("Generating parameters for varying latitude campaign...")
        subprocess.run(["python", "parameter_gen_lat_variation.py"],
                       check=True)
        regex_pattern = r'^parameters_mss_d2d_to_imt_lat_.*_deg.yaml'

    # Run the campaign with the updated regex pattern
    print("Executing campaign with regex pattern:", regex_pattern)
    run_campaign_re(name_campaign, regex_pattern)
