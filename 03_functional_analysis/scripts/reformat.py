import re
import sys

if len(sys.argv) != 3:
    print("Usage: python script.py input_file output_file")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# read the input file
with open(input_file, "r") as file:
    content = file.read()

# replace spaces and tab and "," with newlines, and delete all ":"
content = content.replace(" ", "\n").replace("\t", "\n").replace(":", "").replace(",", "\n")

# add ">" at the beginning of each line, except for the OG/HG lines
lines = content.split('\n')
result = []
for line in lines:
    if line.startswith("HG_") or line.startswith("OG"):
        result.append(line)
    else:
        result.append(">" + line)

# Write the result to the output file
with open(output_file, "w") as file:
    file.write('\n'.join(result))
