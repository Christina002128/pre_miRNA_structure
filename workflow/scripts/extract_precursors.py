### Script called by Snakemake in rule 1 extract_precursors
### Last odified on 1 Dec 2023, Christina Zhang

from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF
import gzip
from BCBio.GFF import GFFExaminer

# #input Snakemake variables (need to use "script" directive)
# import snakemake 
# org = snakemake.wildcards.org
# seq_file = snakemake.input.seq
# ann_file = snakemake.input.ann
# nseqs = snakemake.wildcards.nseqs
# output_file = snakemake.output[0]

# input variables with "shell" directive
import sys
# print(sys.argv)
seq_file = sys.argv[1]
ann_file = sys.argv[2]
nseqs = sys.argv[3]
output_file = sys.argv[4]

# examin types in the GFF file
with open(ann_file) as in_handle:
    examiner = GFFExaminer()
    type_dict = examiner.available_limits(in_handle)["gff_type"]
# if precusor type is not found or is not in 'miRNA_primary_transcript' format, stop the program
if ('miRNA_primary_transcript',) not in type_dict:
    print(f"Type 'miRNA_primary_transcript' is not found in the GFF file! \n Fix it then run it again")
    print(f"The {ann_file} file include types : {type_dict}")
    exit
    

#Concatenate fasta file into a whole genome string
with open(seq_file) as f:
    # Parse the FASTA file
    sequences = [str(record.seq) for record in SeqIO.parse(f, "fasta")]
    # Concatenate all sequences into a single string
    genome = ''.join(sequences).upper()


# only precursors, not matured miRNA 
limit_info = dict(gff_type=[("miRNA_primary_transcript",)])
mirna_dict={} # store miRNA id and sequence

# GFF file parsing
with open(ann_file) as in_handle:
    i=0 
    print("Getting miRNA sequences...")
    for rec in GFF.parse(in_handle, limit_info=limit_info):
        # only extract required number of precursors
        if i>=int(nseqs):
            break
        # Extract DNA sequences that encodes miRNA precursors 
        for feature in rec.features:
            i+=1
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
            # print(i,feature.id,precursor,sep="\t")
            if i>=int(nseqs):
                break

# print(str(mirna_dict['MI0013204']))
# output the miRNA sequence as fasta.gz format
with gzip.open(output_file, "wt") as f:
    for miRNA_id, seq in mirna_dict.items():
        f.write(f">{miRNA_id}\n{str(seq)}\n")

print("Extraction Finished.")




