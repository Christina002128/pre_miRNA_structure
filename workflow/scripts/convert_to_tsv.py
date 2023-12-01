### Script called by Snakemake in rule 3 convert_to_tsv
### Last modified on 1 Dec 2023, Christina Zhang

import gzip
# File paths from Snakemake variables (only under "script" directive)
# input_file = snakemake.input
# output_file = snakemake.output[0]
# input variables with "shell" directive
import sys
#print(sys.argv)
input_file = sys.argv[1]
output_file = sys.argv[2]

# convert format function
def convert_to_tsv(input_file, output_file):
    # read input 
    with gzip.open(input_file, 'rt') as file:
        lines = file.readlines()
    # write output
    with gzip.open(output_file, 'wt') as out:
        # table header
        out.write("ID\tSequence\tStructure\tMFE\n")
        # each miRNA has 3 lines
        for i in range(0, len(lines), 3):
            id_line = lines[i].strip()
            sequence_line = lines[i + 1].strip()
            structure_line, mfe = lines[i + 2].strip().rsplit(" ", 1)
            id = id_line[1:]  # Remove the '>' character
            out.write(f"{id}\t{sequence_line}\t{structure_line}\t{mfe}\n")

# run the function
convert_to_tsv(input_file, output_file)
print("Conversion finished.")