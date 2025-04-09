### The script is to catogorise a GO terms list, using goatools package and database.
### The script is created by Jialin Wei at University of Bristol in 2023, if use, please cite.

from goatools.obo_parser import GODag
obo_fname = "/goatools/go-basic.obo"
# Load the GO DAG
godag = GODag(obo_fname)

# Load the input GO terms list
go_terms_path = "GO_list.txt"
with open(go_terms_path, 'r') as f:
    go_terms = [line.strip().split('\t')[0] for line in f.readlines() if line.strip()]

# Prepare categorized results with descriptions
results = {
    'biological_process': [],
    'cellular_component': [],
    'molecular_function': []
}

for go_id in go_terms:
    go_term = godag.get(go_id, None)
    if go_term:
        category = go_term.namespace.replace(" ", "_") 
        results[category].append(f"{go_id}\t{go_term.name}")
    else:
        print(f"GO term {go_id} not found in the ontology.")

# Output file path
output_file_path = "GOterms_categorized.txt"
# Write categorized lists with descriptions to the output file
with open(output_file_path, 'w') as out_file:
    for category, terms in results.items():
        out_file.write(f"## {category.replace('_', ' ').title()}\n")
        for term in terms:
            out_file.write(term + "\n")
        out_file.write("\n") 

print(f"Categorized GO terms and their descriptions written to {output_file_path}")
