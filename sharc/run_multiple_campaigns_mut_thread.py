import subprocess
import os
import sys
import re
from concurrent.futures import ThreadPoolExecutor


def run_command(param_file, main_cli_path):
    command = [sys.executable, main_cli_path, "-p", param_file]
    subprocess.run(command)


def run_campaign(campaign_name):
    # Path to the working directory
    workfolder = os.path.dirname(os.path.abspath(__file__))
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

    if len(parameter_files) == 0:
        raise ValueError(
            f"No parameter files were found in {campaign_folder}"
        )

    # Number of threads (adjust as needed)
    num_threads = min(len(parameter_files), os.cpu_count())

    # Run the commands in parallel
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        executor.map(
            run_command, parameter_files, [
                main_cli_path,
            ] * len(parameter_files),
        )


def run_campaign_re(campaign_name, param_name_regex):
    """
    Runs a campaign by executing main_cli.py for each parameter file in the specified campaign's input directory
    whose filename matches the given regular expression.

    Args:
        campaign_name (str): The name of the campaign.
        param_name_regex (str): Regular expression to filter parameter file names.
    """
    # Path to the working directory
    workfolder = os.path.dirname(os.path.abspath(__file__))
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

    # Number of threads (adjust as needed)
    num_threads = min(len(parameter_files), os.cpu_count())

    # Run the commands in parallel
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        executor.map(
            run_command, parameter_files, [
                main_cli_path,
            ] * len(parameter_files),
        )


if __name__ == "__main__":
    # Example usage
    run_campaign("imt_hibs_ras_2600_MHz")
