import sys

def read_partial_names(filename):
    """
    Read the list of header prefixes from the file.
    Only include lines that start with '>' and remove the '>' character.
    """
    partial_names = set()
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                partial_names.add(line[1:])
    print("Partial Names Extracted:", partial_names) 
    return partial_names

def is_match(header, partial):
    """
    Returns True if the FASTA header matches the partial header.
    It requires that either the header is exactly equal to the partial,
    or that it starts with the partial and the next character is not a digit.
    """
    if header == partial:
        return True
    if header.startswith(partial):
        if len(header) > len(partial):
            next_char = header[len(partial)]
            if not next_char.isdigit():
                return True
    return False

def find_and_export_genes(partial_names, fasta_filename, output_filename):
    """
    Extract sequences from the FASTA file if the header matches one of the partial names,
    using the is_match() function for precise matching.
    """
    with open(fasta_filename, 'r') as fasta_file, open(output_filename, 'w') as out:
        write_sequence = False  
        for line in fasta_file:
            if line.startswith('>'):  # This is a header line
                header = line.strip()[1:]  # Remove the '>' for matching purposes
                matched = False
                for partial in partial_names:
                    if is_match(header, partial):
                        matched = True
                        print(f"Match Found: {header} (matched with: {partial})") 
                        break
                if matched:
                    write_sequence = True
                    out.write(line)  # Write the header line
                else:
                    write_sequence = False
            elif write_sequence:
                out.write(line)  # Write sequence lines if the header matched

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py headers_file fasta_file output_file")
        sys.exit(1)

    headers_file = sys.argv[1]   # For example, geneids.txt
    fasta_file = sys.argv[2]     # Input FASTA file
    output_file = sys.argv[3]    # The output file for extracted sequences

    partial_names = read_partial_names(headers_file)
    find_and_export_genes(partial_names, fasta_file, output_file)
