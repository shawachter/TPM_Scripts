#!/usr/bin/env python3

#This is me being lazy and not wanting to manually input all sRNA genes into the .gtf file for input into featureCounts.
#Inserting sRNAs into a .gtf file, when sRNAs are in the format: Start	End	Name

import re
from collections import OrderedDict

FASTA = open('/media/shaun/My Passport/Bart_RNA_Seq/sRNA_TPMs/Bartonella_genbank_locus_text_with_sRNAs', 'r')
FASTA2 = open('/home/shaun/R/R-3.4.4/Bart_KC583.gtf','r')
FHOUT = open('Bart_KC583_sRNAs.gtf', 'w')

Dict_BB_Start_End = {}
Dict_BB_Start_Name = {}
Dict_BB_Start_Strand = {}
Dict_Genes_Start_End = {}
Dict_Genes_Start_Name = {}
Dict_Genes_Start_Strand = {}
Sorted_Dict_Genes_Start_End = {}
Sorted_Dict_Genes_Start_Name = {}

#First, get dicts of sRNAs from the genbank file used in the TPM calculator
for line in FASTA:
	line = line.strip('\n')
	fields = re.split('\t', line)
	Name = '"'+str(fields[2])+'";'
	Start = int(fields[0])
	End = int(fields[1])
	Strand = "?"
	if "BB" in Name:
		Dict_BB_Start_End[Start] = End
		Dict_BB_Start_Name[Start] = Name
		Dict_BB_Start_Strand[Start] = Strand

#Now, get dicts of all genes from existing .gtf file for use in featureCounts
for line in FASTA2:
	line = line.strip('\n')
	fields2 = re.split('\t', line)
	Start2 = int(fields2[3])
	End2 = int(fields2[4])
	Strand2 = str(fields2[6])
	Name_Field = fields2[8]
	fields3 = re.split('\s', Name_Field)
	Name2 = str(fields3[5])
	Dict_Genes_Start_End[Start2] = End2
	Dict_Genes_Start_Name[Start2] = Name2
	Dict_Genes_Start_Strand[Start2] = Strand2

#Stupid definition, but it works
def Merge(dict1, dict2): 
    return(dict2.update(dict1)) 

#Merge the sRNA and Gene dicts to produce single dicts containing new sRNA genes
Merge(Dict_BB_Start_End, Dict_Genes_Start_End)
Merge(Dict_BB_Start_Name, Dict_Genes_Start_Name)
Merge(Dict_BB_Start_Strand, Dict_Genes_Start_Strand)

#Order the dicts by starting nucleotide (key)
#The only thing I could get to actually work is the OrderedDict command...
Sorted_Dict_Genes_Start_End = OrderedDict(sorted(Dict_Genes_Start_End.items()))
Sorted_Dict_Genes_Start_Name = OrderedDict(sorted(Dict_Genes_Start_Name.items()))
Sorted_Dict_Genes_Start_Strand = OrderedDict(sorted(Dict_Genes_Start_Strand.items()))


N = 0
for key in Sorted_Dict_Genes_Start_End:
	Value_End = Sorted_Dict_Genes_Start_End[key]
	for key2 in Sorted_Dict_Genes_Start_Name:
		if key == key2:
			Value_Name = Sorted_Dict_Genes_Start_Name[key]
			for key3 in Sorted_Dict_Genes_Start_Strand:
				if key2 == key3:
					Value_Strand = Sorted_Dict_Genes_Start_Strand[key]
					if key == key2 == key3:
						print("NC_008783.1","RefSeq","CDS",key,Value_End,'.',Value_Strand,'0','transcript_id "gene'+str(N)+'"; gene_id "gene'+str(N)+'"; gene_name '+Value_Name, sep = '\t', file = FHOUT)
						N = N+1
