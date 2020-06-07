# TPM_Scripts
A series of python scripts (and one perl script) that will format .gtf and genbank files for the purposes of resulting TPM calculations. Also includes scripts that will average TPM values across biological replicates as well as a script to approximate the TPM contribution of cis-genes for those sRNAs that are cis-anti.

The genbank_locus_text file for input into the TPM_Calculator is in the format Start\tEnd\tLocus_Name where the start and end nucleotides have running 0's before nucleotide position to reach 7 digits (so, 0000001 would be base 1 of the genome). 

The genbank_Locus_Text file was made with the corresponding perl script uploaded here. The genbank file was first amended with the ranges of all sRNAs manually

A .gtf file was also made including all sRNAs for use in featureCounts and DESeq2. This was made with the Sorting_new_gtf.py script included here.

It was necessary to make sure the genbank file (used for TPM calculations) and the .gtf file (used for DESeq2 analysis) contained the same genes. I used the Comparing_genbank_gtf.py script included here to compare the two files and see what's missing where.

After TPM calculations, biological replicate TPMs were averaged, sRNA TPMs were extracted, and the contribution of Cis-genes for Cis-antisense sRNA TPMs was calculated with the No_Filter_BB_TPMs.py script included here.
