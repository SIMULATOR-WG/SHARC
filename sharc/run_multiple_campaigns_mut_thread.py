import subprocess
import os
import sys
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
