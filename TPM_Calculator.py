#!/usr/bin/env python3

#TPM calculator in python3
#We will need all features (nucleotide ranges, lengths, and sum of reads mapping to them based on userplot)
#That's pretty much it...

import re

FASTA = open('/media/shaun/My Passport/Bart_RNA_Seq/sRNA_TPMs/Bartonella_genbank_locus_text_with_sRNAs', 'r')
FASTA2 = open('/media/shaun/My Passport/Bart_RNA_Seq/sRNA_TPMs/Converted_Pl302_Unambiguous.txt','r')
FHOUT = open('ALL_Pl302_TPMs_NEW', 'w') 

Dict_Range = {} #Gene Name and nucleotide range
Dict_Length = {} #Gene name and length
Dict_Start = {} #Gene name and starting nucleotide
Dict_End = {} #Gene name and ending nucleotide
Bb_Gene_List = []
for line in FASTA:
	if re.match(r'^\d', line):
		line = line.strip('\n')
		fields = re.split('\t', line)
		Name = str(fields[2])
		Start = int(fields[0])
		End = int(fields[1])
		Length = int(End - Start)
		Range = range(Start, End)
	Dict_Range[Name] = Range
	Dict_Length[Name] = Length
	Dict_Start[Name] = Start
	Dict_End[Name] = End
	Bb_Gene_List.append(Name)

Dict_Cds_Reads = {}
Dict_Gene_T = {}
Dict_Nuc_Reads = {}

#iterate over lines and make a dictionary of nucleotide and read number...
for line in FASTA2:
	line = line.strip('\n')
	fields2 = re.split(',', line)
	Nuc = int(fields2[0])
	Reads = int(fields2[1])
	Dict_Nuc_Reads[Nuc] = Reads


for key in Dict_Range: #for every gene name
	Gene_Reads = 0
	value = Dict_Range[key] #we have the range of nucs
	for nuc_key in Dict_Nuc_Reads: #for every nucleotide
		Nuc_Reads = Dict_Nuc_Reads[nuc_key] #we have reads
		if nuc_key in value:
			Gene_Reads = Nuc_Reads + Gene_Reads
		else:
			continue
	Dict_Cds_Reads[key] = Gene_Reads #dictionary of gene names and sum of reads	
	for key2 in Dict_Length:
		value2 = Dict_Length[key2] #Lengths of genes
		if key == key2:
			Gene_Length = value2
		else:
			continue
	Gene_T = (int(Gene_Reads) * 150) / int(Gene_Length)
	Dict_Gene_T[key] = Gene_T # Dictionary of gene names and their T values	

Overall_T = 0
for key3 in Dict_Gene_T:
	value3 = Dict_Gene_T[key3]
	Overall_T = int(value3) + int(Overall_T) #calculate overall T value. Have to iterate here by itself to make this work

print("The overall T for this sample is", Overall_T, sep = ' ')

for key4 in Bb_Gene_List: #This is a list of gene names, because if you use a dictionary the genes will be random and all out of order.
	if key4 in Dict_Start:
		Final_Gene_Start = Dict_Start[key4] #Gene start
	if key4 in Dict_End:
		Final_Gene_End = Dict_End[key4] #Gene end
	if key4 in Dict_Length:
		Final_Gene_Length = Dict_Length[key4] #Gene length
	if key4 in Dict_Cds_Reads:
		Final_Gene_Reads = Dict_Cds_Reads[key4] #Reads mapping to gene.
	Numerator = int(Final_Gene_Reads) * 150 * 1000000
	Denominator = (int(Final_Gene_Length) * int(Overall_T))
	TPM = float(Numerator) / float(Denominator)
	print(Final_Gene_Start, Final_Gene_End, Final_Gene_Length, TPM, key4, sep = '\t', file = FHOUT)

