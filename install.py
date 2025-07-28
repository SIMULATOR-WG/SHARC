import os
import sys
import subprocess
from pathlib import Path
import platform


def run_command(command, shell=False):
    try:
        subprocess.check_call(command, shell=shell)
    except subprocess.CalledProcessError as e:
        print(f"Error: Command '{' '.join(command)}' failed with exit code {e.returncode}")
        sys.exit(e.returncode)


def main():
    if int("".join(sys.version.split(" ")[0].split("."))) > 3120:
        print("- Compatible Python version successfully recognized, starting installation process...")
        start_install()
    else:
        print(f"- Error: Try updating Python to version 3.12 or above! \n - Current Python version: {sys.version}")


def start_install():
    print("- SHARC Installer Starting...")

    # Set up path
    project_root = Path(__file__).resolve().parent
    venv_path = project_root / ".venv"

    # Create venv
    if not venv_path.exists():
        print("- Creating virtual environment...")
        run_command([sys.executable, "-m", "venv", str(venv_path)])
    else:
        print("- Virtual environment already exists.")

    # Activate venv
    if platform.system() == "Windows":
        pip_executable = venv_path / "Scripts" / "pip.exe"
    else:
        pip_executable = venv_path / "bin" / "pip"

    if not pip_executable.exists():
        print("- Could not locate pip inside virtual environment.")
        sys.exit(1)

    requirements_file = project_root / "requirements.txt"
    if not requirements_file.exists():
        print("- requirements.txt not found.")
        sys.exit(1)

    print("- Installing dependencies from requirements.txt...")
    run_command([str(pip_executable), "install", "-r", str(requirements_file)])

    print("- Installing SHARC in editable mode...")
    run_command([str(pip_executable), "install", "-e", str(project_root)])

    print("\n- SHARC installation complete!")
    print("- To activate your virtual environment, run:")
    if platform.system() == "Windows":
        print(r".venv\Scripts\activate")
    else:
        print("source .venv/bin/activate")


if __name__ == "__main__":
    main()
