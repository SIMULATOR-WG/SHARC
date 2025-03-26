import subprocess
import os
import sys


def run_campaign(campaign_name):
    # Get the current working directory
    workfolder = os.path.dirname(os.path.abspath(__file__))

    # Path to the main_cli.py script
    main_cli_path = os.path.join(workfolder, "main_cli.py")

    # Campaign directory
    campaign_folder = os.path.join(
        workfolder, "campaigns", campaign_name, "input",
    )

    # List of parameter files
    parameter_files = [
        os.path.join(campaign_folder, f) for f in os.listdir(
            campaign_folder,
        ) if f.endswith('.yaml')
    ]

    # Run the command for each parameter file
    for param_file in parameter_files:
        command = [sys.executable, main_cli_path, "-p", param_file]
        subprocess.run(command)


if __name__ == "__main__":
    # Example usage
    run_campaign("imt_hibs_ras_2600_MHz")
