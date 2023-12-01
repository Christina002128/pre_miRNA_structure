from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF
import gzip
import snakemake
# File paths from Snakemake variables
org = snakemake.wildcards.org
seq_file = snakemake.input.seq
ann_file = snakemake.input.ann
nseqs = snakemake.wildcards.nseqs
output_file = snakemake.output[0]

#Concatenate fasta file into a whole genome string
with open(seq_file) as f:
    # Parse the FASTA file
    sequences = [str(record.seq) for record in SeqIO.parse(f, "fasta")]
    # Concatenate all sequences into a single string
    genome = ''.join(sequences).upper


mirna_dict={} # store miRNA id and sequence

# GFF file parsing
with open(ann_file) as in_handle:
    for rec in GFF.parse(in_handle):
        i=0
        # Extract DNA sequences that encodes miRNA precursors 
        for feature in rec.features:
            # only precursors, not matured miRNA
            if feature.type == "miRNA_primary_transcript":
                # Extract the sequence position
                start = feature.location.start.position - 1
                end = feature.location.end.position
                # forward or reverse strand
                strand = feature.strand
                if strand==int(1):
                    precursor = Seq(genome[start:end]).transcribe()
                elif strand==int(-1):
                    precursor = Seq(genome[start:end]).reverse_complement().transcribe()
                else:
                    print(f"Ambiguous strand ({feature.id})! Treat as + in default")
                    precursor = Seq(genome[start:end]).transcribe()
                #print(feature.id,len(precursor),precursor)
                mirna_dict[feature.id] = precursor
                i+=1
                if i>=10:
                    break

# print(str(mirna_dict['MI0013204']))
# output the miRNA sequence as fasta.gz format
with gzip.open(output_file, "wt") as f:
    for miRNA_id, seq in mirna_dict.items():
        f.write(f">{miRNA_id}\n{str(seq)}\n")
