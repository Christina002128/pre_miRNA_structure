import gzip
# File paths from Snakemake variables
org = snakemake.wildcards.org
input_file = snakemake.input
output_file = snakemake.output[0]

# convert format function
def convert_to_tsv(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

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


convert_to_tsv(input_file, output_file)
