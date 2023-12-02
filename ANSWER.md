## Answers for Theoretical questions

#### Q1. One of Tropic's main activity revolves sequencing small RNA for analysis of miRNA expression in different tissues. Please outline what key challenges can be encountered when handling sRNAseq data.

miRNA expression can be regulated at many steps including transcription, Drosha and Dicer processing. Since only the matured miRNA has actual function, I would focus on the matured miRNA expression analysis.

1. Multimaping. Because miRNA are very short sequences (18-22nt), when mapping reads to reference genome, the ambiguous mapping can be very high (up to 40%). The solutions would be either throwing away these reads, or using different strategies to assign them to their multimaping areas.

2. IsomiRs Detection. IsomiRs are the miRNA variants, which contain mutations in their sequences. The mutation might be one or two nucleotides mutation, insertion or deletion, whick make it hard to detect or causes false mapping. 

3. The degradation of long RNA sequences can be similar to small RNA sequences, which is hard to distinguish, but using quality control tools like "mirtrance"" can tell the overall degradation level of the sample to minimize this effect. 

4. For non-model organisms, miRNA annotations might be limited, which means identifying novel miRNAs would be a necessary process. 

5. miRNA expression results can depend on many factors, including the nature of samples, experimental condition and choice of technical methods for miRNA detection. However, under good experimental design and well-controlled experiments, along with normalization methods for data, the problem can be overcome.

In summary, although the challenges listed above makes it hard to get consistent conclusion from miRNA expression analysis, there are piplines or tools for miRNA analysis are developped to help reduce the effect and make the best use of the data, such as miRDeep2(https://www.mdc-berlin.de/content/mirdeep2-documentation?mdcbl%5B0%5D=/n-rajewsky%23t-data%2Csoftware%26resources&mdctl=0&mdcou=20738&mdcot=6&mdcbv=71nDTh7VzOJOW6SFGuFySs4mus4wnovu-t2LZzV2dL8), miRge3(https://mirge3.readthedocs.io/en/master/installation.html) and so on. 

#### Q2. Following identification of a novel non-transgenic, gene edited banana plant, we aim to submit this candidate to a regulatory body for exemption from regulation. To this end, we are required to demonstrate the complete absence of transgenic sequences derived from our vector. To this extent, whole genome sequencing was performed on the plant of interest. Please describe the process you would follow from obtaining the raw sequencing data to proving absence of transgenic material in the genome.

In order to prove the genome sequence doesn’t have transgenic seuquences,  reference genome of the banana species is required. If the whole genome for specific cultivar is not applicable, same 
species such as Musa acuminata can be also used for referencing. After quality control and necessary trimming, align the data to reference genome and save the unmapped sequences. 
If there are many unmapped reads, try assemble them into long sequences. Then use alignment tool such as blast(https://blast.ncbi.nlm.nih.gov/Blast.cgi) to check if the sequences exist in other organisms. 

If there is no gene from other species aligned to them, there’s no transgenic sequences from other organisms exist in our sample. If there is, can try to find out the reason, calculate how likely it is caused by gene editing or happened naturally, in order to prove there is no transgenic sequences.



#### Additionall Content

*Apart from the above two questions I chosed for this assignment, the following two questions are also interesting. Therefore, I did some research and would like to share some of my thoughts for these questions.*


#### Q3. Whilst scouring the literature, you found a gene of interest providing drought resistance in Coconut (_Cocos nucifera_). You want to know if a similar gene is present in the coffee genome (_Coffea arabica_). Please detail the steps you would undertake.

Get the gene (or protein, mRNA) name or identification code from the literature, download the gene sequence from database like NCBI (https://www.ncbi.nlm.nih.gov), then blast(https://blast.ncbi.nlm.nih.gov/Blast.cgi) the sequence to protein sequences in Coffea arabica. 

After getting the homologous proteins with high similarities, check the detailed information of these sequences stored in NCBI database. If there is one sequences that have the same function, then we can say the coffee also have the gene of interest found in coconut, and potentially also provide drought resistance.


#### Q4. GEiGS Biocompute (https://www.geigs.com/) predicts solutions allowing to redirect endogenous miRNAs against novel endogenous or pathogenic targets. Please explain potential challenges you might encounter during implementation of GEiGS solutions _in planta_.
1. Challenge of protoplast regeneration. Some species are challenging to do protoplast regeneration technique. Tobacco, Arabidopsis and other model plants has developed this technique, while other plants like wheat, coconut hasn’t.

2. Gene knockout often is easier than mutations or insertions while using gene editing tools like CRISPR Cas9. If the goal is to mutate the sequence into a desired version, the success rate might be low, which requires more samples and experiments to achieve.

3. Gene accessibility problem. Some genes are more difficult to target due to their location in the genome, sequence similarity with other genes, or because they are part of gene families.

