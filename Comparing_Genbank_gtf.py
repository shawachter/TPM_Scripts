#!/usr/bin/env python3

#There are some genes in the genbank file not in the gtf file (tRNAs for one, which is fine)
#Likewise, there are some genes in the gtf file I've removed from the genbank file (pseudogenes, rRNA genes...)
#This is to see what's different between the two, so I can manually fix the gtf file and make it identical to the genbank file I used for sRNA TPMs (UpSet plot, heatmap, etc)

import re
from collections import OrderedDict

FASTA = open('/media/shaun/My Passport/Bart_RNA_Seq/sRNA_TPMs/Bartonella_genbank_locus_text_with_sRNAs', 'r')
FASTA2 = open('/media/shaun/My Passport/Bart_RNA_Seq/sRNA_TPMs/Bart_KC583_sRNAs.gtf','r')
FHOUT = open('Whats_Missing.txt', 'w')

Dict_BB_Start_Name = {}
Dict_Genes_Start_Name = {}
Sorted_Dict_BB_Start_Name = {}
Sorted_Dict_Genes_Start_Name = {}

for line in FASTA:
	line = line.strip('\n')
	fields = re.split('\t', line)
	Name = fields[2]
	Start = int(fields[0])
	End = int(fields[1])
	Dict_BB_Start_Name[Start] = Name

for line in FASTA2:
	line = line.strip('\n')
	fields2 = re.split('\t', line)
	Start2 = int(fields2[3])
	End2 = int(fields2[4])
	Strand2 = str(fields2[6])
	Name_Field = fields2[8]
	fields3 = re.split('\s', Name_Field)
	Name2 = str(fields3[5])
	Dict_Genes_Start_Name[Start2] = Name2

Sorted_Dict_BB_Start_Name = OrderedDict(sorted(Dict_BB_Start_Name.items()))
Sorted_Dict_Genes_Start_Name = OrderedDict(sorted(Dict_Genes_Start_Name.items()))

for key in Sorted_Dict_BB_Start_Name:
	Name_Value = Sorted_Dict_BB_Start_Name[key]
	if key in Sorted_Dict_Genes_Start_Name:
		continue
	else:
		print("Genbank not in gtf", key, Name_Value, file = FHOUT)

for key2 in Sorted_Dict_Genes_Start_Name:
	Name_Value2 = Sorted_Dict_Genes_Start_Name[key2]
	if key2 in Sorted_Dict_BB_Start_Name:
		continue
	else:
		print("gtf not in genbank", key2, Name_Value2, file = FHOUT)
