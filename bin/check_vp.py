import argparse


# Parse command line arguments
def parse_commandline():
    """Parse the arguments given on the command-line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check_vp", help="Path to the check_vp file")
    parser.add_argument(
        "--simulated_phenos", help="path to the simulated phenotypes file"
    )
    args = parser.parse_args()

    return args


def get_vp(check_vp_file):
    """
    Read in the check_vp file a tab delim file with the following columns:
    1. Source
    2. Variance
    3. Standard error
    And pull out variance for the Vp source.
    Raises ValueError if Vp is not found.
    """
    vp_value = None
    with open(check_vp_file, "r") as f:
        for line_content in f:
            parts = line_content.strip().split("\t")
            if len(parts) > 1 and parts[0] == "Vp":
                vp_value = parts[1]
                break  # Found Vp, no need to read further

    if vp_value is None:
        raise ValueError(f"Vp not found in {check_vp_file}")
    return vp_value


def increase_pheno_var(phenos_file):
    """
    Read in the simulated phenotypes file (space delimited).
    Multiply the phenotype value (3rd column) by 1000.
    Returns a list of strings, where each string is a processed line for the output file.
    """
    processed_lines = []
    with open(phenos_file, "r") as f_in:
        for line_content in f_in:
            parts = line_content.strip().split(" ")
            if len(parts) >= 3:
                try:
                    pheno = float(parts[2])
                    new_pheno = pheno * 1000
                    processed_lines.append(f"{parts[0]} {parts[1]} {str(new_pheno)}\n")
                except ValueError:
                    print(
                        f"Warning: Could not convert phenotype in line: {line_content.strip()} - Skipping line."
                    )
            else:
                print(
                    f"Warning: Malformed line in {phenos_file}: {line_content.strip()} - Skipping line."
                )
    return processed_lines


def read_original_phenos(phenos_file):
    """
    Reads the simulated phenotypes file and returns its lines as a list of strings.
    Each string includes its original newline character.
    """
    original_lines = []
    with open(phenos_file, "r") as f_in:
        for line_content in f_in:
            original_lines.append(line_content)
    return original_lines


def write_phenos_to_temp(lines_to_write, output_filename="new_phenos.temp"):
    """
    Writes the given list of lines to the specified output file,
    overwriting the file if it exists.
    """
    with open(output_filename, "w") as f_out:
        for line in lines_to_write:
            f_out.write(line)


if __name__ == "__main__":
    args = parse_commandline()
    check_vp_file = args.check_vp
    phenos_file = args.simulated_phenos

    try:
        vp = get_vp(check_vp_file)
        print(f"Vp: {vp}")

        output_lines_for_temp_file = []

        if float(vp) < 0.000001:
            print(
                "Vp is less than 0.000001. Increasing phenotype variance and writing to new_phenos.temp."
            )
            output_lines_for_temp_file = increase_pheno_var(phenos_file)
        else:
            print(
                "Vp is greater than or equal to 0.000001. Writing original phenotypes to new_phenos.temp."
            )
            output_lines_for_temp_file = read_original_phenos(phenos_file)

        write_phenos_to_temp(output_lines_for_temp_file)
        print(f"Output written to new_phenos.temp")

    except FileNotFoundError as e:
        print(f"Error: File not found. {e}")
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
