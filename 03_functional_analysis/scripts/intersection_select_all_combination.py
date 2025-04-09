### The script is to find intersections of GO terms from specified number of TXT files, the intersected ones might also present in other TXT files, that is, not excluding
### The script is created by Jialin Wei at University of Bristol in 2023, if use, please cite.

import os
import argparse
import itertools

def read_file_to_set(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        return set(line.strip() for line in file)

# extract the unique identifier from a filename
def extract_identifier(file_name):
    
    return file_name.split('_')[0]

def find_intersections(input_directory, num_files, output_directory):
    all_files = os.listdir(input_directory)
    all_files = [f for f in all_files if f.endswith('.txt')] 

    # all combinations of the specified number of files
    file_combinations = itertools.combinations(all_files, num_files)

    # process each combination
    for combination in file_combinations:
        # read the files and find the intersection
        sets = [read_file_to_set(os.path.join(input_directory, file_name)) for file_name in combination]
        intersection = set.intersection(*sets)

        # create output file name based on the unique identifiers of the combination
        identifiers = [extract_identifier(file_name) for file_name in combination]
        output_file_name = f"intersection_{'_'.join(identifiers)}.txt"
        output_file_path = os.path.join(output_directory, output_file_name)

        # write the intersection to the output file
        with open(output_file_path, 'w') as out_file:
            for element in intersection:
                out_file.write(f"{element}\n")

        print(f"Intersection written to {output_file_path}")

def main():
    parser = argparse.ArgumentParser(description='Find intersections of GO terms from all combinations of a specified number of files.')
    parser.add_argument('input_directory', help='Path to the input directory containing text files')
    parser.add_argument('output_directory', help='Path to the directory where output files will be saved')
    parser.add_argument('num_files', type=int, help='Number of files to select for each combination')
    
    args = parser.parse_args()

    find_intersections(args.input_directory, args.num_files, args.output_directory)

if __name__ == "__main__":
    main()
