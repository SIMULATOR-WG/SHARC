"""
_summary_
    Script to convert all .ini files found in the current directory and its subdirectories to .yaml,
    saving them in the same directory.
    Usage: `python3 script.py`
    Caveats: currently, it does not parse .ini lists.
"""

import re
import os

num_spaces_ident = 4


def convert_ini_to_yaml(input_file, output_file):

    print("Reading from file: ", input_file)
    print("Writing to file: ", output_file)

    input_file_content = open(input_file, "r").readlines()

    with open(output_file, 'w') as output_file:
        current_ident: int = 0  # Identation to be used in yaml file
        current_section: str = ""  # Section of ini file
        current_attr_name: str = ""  # Attribute name of ini file
        current_attr_comments: list[str] = []
        for line in input_file_content:
            if (line.isspace() or not line):
                continue
            line.strip()
            if (line.startswith("#")):
                current_attr_comments.append(line)
                continue
            elif (line.startswith("[")):
                # Line is section
                # go back 1 identation
                current_ident = max(current_ident - num_spaces_ident, 0)
                section_name = re.findall(r'\[(.+)\]', line)[0]
                section_name = section_name.lower()
                current_section = section_name
                for comment in current_attr_comments:
                    output_file.write(' ' * current_ident + comment)
                current_attr_comments = []
                output_file.write(
                    ' ' * current_ident +
                    f"{current_section}:\n",
                )
                current_ident += num_spaces_ident
            else:
                # line is attribute
                try:
                    current_attr_name = re.findall(r"(.+) *=(.+)", line)[0][0]
                    current_attr_value = re.findall(r"(.+) *=(.+)", line)[0][1]
                    if (current_attr_value == 'TRUE' or current_attr_value == 'True'):
                        current_attr_value = 'true'
                    if (current_attr_value == 'FALSE' or current_attr_value == 'False'):
                        current_attr_value = 'false'
                    for comment in current_attr_comments:
                        output_file.write(' ' * current_ident + comment)
                    current_attr_comments = []
                    output_file.write(
                        ' ' * current_ident + f"{current_attr_name} :{current_attr_value}\n",
                    )
                except IndexError:
                    print(input_file, 'did not match the expected format. Skipping...')

    print(f"Conversion complete: {output_file.name}")


if __name__ == "__main__":

    root_dir = os.getcwd()  # Use the current working directory as the root
    # Add your environment folder names here
    exclude_dirs = {'env', 'venv', '.venv'}

    ini_files = []
    for root, dirs, files in os.walk(root_dir):

        # Modify dirs in-place to exclude specific directories
        dirs[:] = [d for d in dirs if d not in exclude_dirs]

        for file in files:
            if file.endswith('.ini'):
                ini_files.append(os.path.join(root, file))

    if not ini_files:
        print("No .ini files found in the current directory or its subdirectories.")
        exit(1)

    print("Detected .ini files:")
    for ini_file in ini_files:
        print(ini_file)

    for ini_file in ini_files:
        output_file_name = os.path.splitext(ini_file)[0] + '.yaml'
        convert_ini_to_yaml(ini_file, output_file_name)
