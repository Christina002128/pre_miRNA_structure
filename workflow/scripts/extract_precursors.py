# import gzip

# with gzip.open(snakemake.output[0], "wt") as ofh:
#     ofh.write(f"Here we extract precursors for {snakemake.wildcards.org}\n")


import gzip
from Bio import SeqIO
from BCBio import GFF

def extract_miRNA_precursors(ann_file, seq_file):
    # Read the genome sequence
    genome_seq = SeqIO.to_dict(SeqIO.parse(seq_file, "fasta"))

    # Open the GFF file and parse it
    with open(ann_file) as in_handle:
        for rec in GFF.parse(in_handle, base_dict=genome_seq):
            # Extract miRNA precursors
            for feature in rec.features:
                if feature.type == "miRNA":
                    # Extract the precursor sequence
                    start = feature.location.start.position
                    end = feature.location.end.position
                    precursor_seq = rec.seq[start:end]
                    yield feature.id, precursor_seq

# File paths from Snakemake variables
org = snakemake.wildcards.org
seq_file = snakemake.input.seq
ann_file = snakemake.input.ann
output_file = snakemake.output[0]

# Extract and write the precursors
with gzip.open(output_file, "wt") as ofh:
    for miRNA_id, seq in extract_miRNA_precursors(ann_file, seq_file):
        ofh.write(f">{miRNA_id}\n{seq}\n")
