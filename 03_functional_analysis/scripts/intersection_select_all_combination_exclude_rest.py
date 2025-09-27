### The script is to find intersections of GO terms from specified number of TXT files, excluding others
### The script is created by Jialin Wei at University of Bristol in 2023, if use, please cite.

import os
import argparse
import itertools

def read_file_to_set(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        return set(line.strip() for line in file)

def extract_identifier(file_name):
    return file_name.split('_')[0]

def find_intersections(input_directory, num_files, output_directory):
    all_files = os.listdir(input_directory)
    all_files = [f for f in all_files if f.endswith('.txt')] 

    file_combinations = itertools.combinations(all_files, num_files)

    for combination in file_combinations:
        selected_sets = [read_file_to_set(os.path.join(input_directory, file_name)) for file_name in combination]
        intersection = set.intersection(*selected_sets)

        # Files not in the current combination
        non_selected_files = set(all_files) - set(combination)
        non_selected_sets = [read_file_to_set(os.path.join(input_directory, file_name)) for file_name in non_selected_files]
        union_non_selected = set.union(*non_selected_sets) if non_selected_sets else set()

        # Exclude elements present in non-selected files
        exclusive_intersection = intersection - union_non_selected

        identifiers = [extract_identifier(file_name) for file_name in combination]
        output_file_name = f"intersection_exclusive_{'_'.join(identifiers)}.txt"
        output_file_path = os.path.join(output_directory, output_file_name)

        with open(output_file_path, 'w') as out_file:
            for element in exclusive_intersection:
                out_file.write(f"{element}\n")

        print(f"Exclusive intersection written to {output_file_path}")

def main():
    parser = argparse.ArgumentParser(description='Find intersections of GO terms from specified number of files, excluding others.')
    parser.add_argument('input_directory', help='Path to the input directory containing text files')
    parser.add_argument('output_directory', help='Path to the directory where output files will be saved')
    parser.add_argument('num_files', type=int, help='Number of files to select for each combination')
    
    args = parser.parse_args()

    find_intersections(args.input_directory, args.num_files, args.output_directory)

if __name__ == "__main__":
    main()
