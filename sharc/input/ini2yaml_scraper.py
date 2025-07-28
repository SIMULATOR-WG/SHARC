"""
_summary_
    Script to be run when you need to convert between old .ini parameters and new .yaml parameters.
    Usage: `python3 -i input-file.ini -o output-file.yaml`
    Caveats: currently, it does not parse .ini lists.
"""

if __name__ == "__main__":

    import re
    from argparse import ArgumentParser

    num_spaces_ident = 4

    parser = ArgumentParser()
    parser.add_argument(
        "-i", "--input", dest="input_file",
        help="Input file. Should be .ini. Absolute Path", required=True,
    )
    parser.add_argument(
        "-o", "--output", dest="output_file",
        help="Output file. Should be .yaml. Absolute Path", required=True,
    )

    args = parser.parse_args()

    print("Reading from file: ", args.input_file)
    print("Writing to file: ", args.output_file)

    input_file_content = open(args.input_file, "r").readlines()

    # if os.path.exists(args.output_file):
    #     response = input(f"File {args.output_file} already exists, are you sure that you want to override it? Y/n")
    # if(response == "n" or response == "N"):
    #     print("**** Operation Cancelled ****")
    #     exit(1)
    with open(args.output_file, 'w') as output_file:
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
                    ' ' * current_ident +
                    f"{current_attr_name} :{current_attr_value}\n",
                )
