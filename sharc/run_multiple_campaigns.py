import subprocess
import os
import sys
import re


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


def run_campaign_re(campaign_name, param_name_regex):
    """
    Run a campaign for parameter files matching a given regular expression.

    Args:
        campaign_name (str): Name of the campaign.
        param_name_regex (str): Regular expression to filter parameter files.
    """
    # Get the current working directory
    workfolder = os.path.dirname(os.path.abspath(__file__))

    # Path to the main_cli.py script
    main_cli_path = os.path.join(workfolder, "main_cli.py")

    # Campaign directory
    campaign_folder = os.path.join(
        workfolder, "campaigns", campaign_name, "input",
    )

    # List of parameter files
    pat = re.compile(param_name_regex)
    parameter_files = [
        os.path.join(campaign_folder, f) for f in os.listdir(
            campaign_folder,
        ) if pat.match(f)
    ]
    if len(parameter_files) == 0:
        raise ValueError("No parameter files where found with given filter")

    # Run the command for each parameter file
    for param_file in parameter_files:
        command = [sys.executable, main_cli_path, "-p", param_file]
        subprocess.run(command)


if __name__ == "__main__":
    # Example usage
    run_campaign("imt_hibs_ras_2600_MHz")
