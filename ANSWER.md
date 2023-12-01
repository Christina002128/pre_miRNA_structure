1challenge for handling sRNAseq data
1. multimapping. short length, ambiguous mapping, small RNAs come from multiple precursors over the genome so multimapping is high (up to 40%)
2. isomiRs detection, which are the variants of a miRNA sequence, have mutated nucleotides and hard to detect
3. short length, is it degradation or actual small RNA
4. non-model organism, not enough annotations exist
5. identify novel: hard to tell novel small RNA or degredation
6. predicting targets: multiple candidates
7. quality control?? normalization?

short sequence, ambiguous mapping—low specificity, different stage of small RNA (pri,pre,matured)
PCR bias, sequence specific bias
low expression, 

ligation bias in libary construction, randomize adaptor sequencing (in wet-lab)

2prove the genome sequence has no transgenic material in genome, describe process
in order to prove the genome sequence doesn’t have transgenic seuquences,  reference genome of the banana species is required. If the whole genome for the specific cultivar is not applicable, same 
specises such as Musa acuminata can be also used for referencing. 
After quality control and necessary trimming, align the data to reference genome. check the unmapped sequences. 
If there’s many unmapped reads, try assemble them into long sequences, then use blast to align them against all other organisms. If there is no gene from other species aligned to them, there’s no 
transgenic sequences from other organisms exist in our sample. 
If there is, can try to find out the reason, calculate how likely it is caused by gene editing or happened naturally, in order to prove there is no transgenic sequences.

3found gene of drought resistance in Coconut (Cocos nucifera), if a similar gene is present in the coffee genome (Coffea arabica). Please detail the steps you would undertake
get sequence of gene, blast to coffee, find high significant genes, protein database, structure and function prediction using alphafold

4redirect endogenous miRNA against novel endogenous or pathogenic targets, what are the challenges to implement this in plant
direct mutation is hard. 
if use new templet sequence, protoplast regernation is hard, cause some species are challenging to do protoplast regernation technique. like tobacco, Arabidopsis and other model plants has 
developped technique. but wheat, coconut hasn’t.
