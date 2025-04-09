### The script is to find the bottom GO terms (most specific GO terms) from a GO terms list, using goatools package and database.
### The script is created by Jialin Wei at University of Bristol in 2023, if use, please cite.

import argparse
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag

parser = argparse.ArgumentParser(description="Find bottom GO terms from a list.")
parser.add_argument("input_file", help="Path to the input file containing GO terms.")
parser.add_argument("output_file", help="Path to the output file to write bottom GO terms.")
args = parser.parse_args()

# Load all relationships
godag = get_godag("goatools/go-basic.obo", optional_attrs={'relationship'})

# Load GO terms list
with open(args.input_file, 'r') as file:
    go_terms_list = [line.strip() for line in file.readlines()]

# Initialize a list to hold bottom terms and descriptions
bottom_terms_with_desc = []

# For each GO term in the list, check if there are children within the list
for go_term in go_terms_list:
    if go_term not in godag:
        print(f"Skipping {go_term}: Not found in GO DAG.")
        continue
    children = [child.id for child in godag[go_term].children if child.id in go_terms_list]
    if not children:  # If the GO term has no children in the list, then it is a bottom GO
        # Append the GO term and its description
        bottom_terms_with_desc.append((go_term, godag[go_term].name))

# Write to the output file
with open(args.output_file, 'w') as file:
    for term, desc in bottom_terms_with_desc:
        file.write(f"{term}\t{desc}\n")

print(f"Bottom GO terms and their descriptions written to {args.output_file}")

