import gzip

with gzip.open(snakemake.output[0], "wt") as ofh:
    ofh.write(f"Here we extract precursors for {snakemake.wildcards.org}\n")
