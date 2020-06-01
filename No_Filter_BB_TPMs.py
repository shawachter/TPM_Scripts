#!/usr/bin/env python3
#Here we are simply averaging 2 values across multiple files and combining the data into one file. 
#Then, we are opening that new file and filtering out the TPMs <300
#We are also modifying Cis-anti sRNA TPMs to take into account the contribution of the cis-gene.

import re

FASTA = open('BB_HUVE1_TPMs','r') #sRNA TPMs
FASTA2 = open('BB_HUVE2_TPMs', 'r') #sRNA TPMs
FHOUT = open('END_BB_HUVE_AVERAGE_TPMs', 'w') #sRNA average TPMs across biological replicates

Dict1 = {} #BbsR Name and TPM

for line in FASTA:
	if re.match(r'^\d', line):
		line = line.strip('\n')
		fields = re.split('\s', line)
		TPM = float(fields[2])
		Name = str(fields[3])
		Start = int(fields[0])
		End = int(fields[1])
		BB_Range = range(Start, End)
	Dict1[Name] = TPM

for line in FASTA2:
	line = line.strip('\n')
	fields2 = re.split('\s', line)
	Name2 = str(fields2[3])
	TPM2 = float(fields2[2])
	for key in Dict1:
		value = Dict1[key]
		if key == Name2:
			AVG_TPM = (value + TPM2)/2
			print(Name2, AVG_TPM, sep = '\t', file = FHOUT)
		else:
			continue
	continue
FHOUT.close()	

#Next, we need to average the 2-3 relevant whole_TPM files for use later

FASTA4 = open('ALL_HUVE1_TPMs','r') #Total RNA TPMs
FASTA5 = open('ALL_HUVE2_TPMs','r') #Total RNA TPMs
FHOUT3 = open('END_Total_HUVE_AVERAGE_TPMs', 'w') #Total RNA average TPMs across biological replicates

Dict2 = {} 
Dict4 = {} # gene name and average TPM

for line in FASTA4:
	if re.match(r'^\d', line):
		line = line.strip('\n')
		fields3 = re.split('\t', line)
		TPM3 = float(fields3[3])
		Name3 = str(fields3[4])
		Start = int(fields3[0])
		End = int(fields3[1])		
	else:
		continue
	Dict2[Name3] = TPM3

for line in FASTA5:
	if re.match(r'^\d', line):
		line = line.strip('\n')
		fields4 = re.split('\t', line)
		Name4 = str(fields4[4])
		TPM4 = float(fields4[3])
		for key2 in Dict2:
			value2 = Dict2[key2]
			if key2 == Name4:
				AVG_TPM2 = (value2 + TPM4)/2
				print(Name4, AVG_TPM2, sep = '\t', file = FHOUT3)
				Dict4[Name4] = AVG_TPM2 #Making our dictionary of genes and average TPMs for use later... 
			else:
				continue
		continue
FHOUT3.close()	

#Here we go calculating the overlap between the two ranges stored as values in two separate dictionaries
#Please note that TPM files do not display ranges properly; we have to use the genbank file to populate ranges.

Cis_Dict = {'BB001':'BARBAKC583_0002', 'BB002':'BARBAKC583_0004', 'BB004':'BARBAKC583_0027', 'BB008':'BARBAKC583_0064', 'BB009':'BARBAKC583_0076', 'BB010':'BARBAKC583_0089', 'BB012':'BARBAKC583_0094', 'BB013':'BARBAKC583_0099', 'BB017':'BARBAKC583_0154', 'BB028':'BARBAKC583_0271', 'BB032':'BARBAKC583_0296', 'BB033':'BARBAKC583_0298', 'BB034':'BARBAKC583_0306', 'BB037':'BARBAKC583_0329', 'BB040':'BARBAKC583_0367', 'BB042':'BARBAKC583_0387', 'BB043':'BARBAKC583_0391', 'BB044':'BARBAKC583_0414', 'BB047':'BARBAKC583_0432', 'BB049':'BARBAKC583_0452', 'BB050':'BARBAKC583_0483', 'BB052':'BARBAKC583_0492', 'BB053':'BARBAKC583_0495', 'BB054':'BARBAKC583_0506', 'BB056':'BARBAKC583_0548', 'BB058':'BARBAKC583_0555', 'BB059':'BARBAKC583_0559', 'BB061':'BARBAKC583_0579', 'BB062':'BARBAKC583_0589', 'BB063':'BARBAKC583_0600', 'BB064':'BARBAKC583_0601', 'BB065':'BARBAKC583_0622', 'BB067':'BARBAKC583_0636', 'BB068_1':'BARBAKC583_0718', 'BB068_2':'BARBAKC583_0718', 'BB069':'BARBAKC583_0735', 'BB070':'BARBAKC583_0774', 'BB072':'BARBAKC583_0780', 'BB073':'BARBAKC583_0781', 'BB074':'BARBAKC583_0792', 'BB076':'BARBAKC583_0814', 'BB077':'BARBAKC583_0816', 'BB078_1':'BARBAKC583_0825', 'BB078_2':'BARBAKC583_0825', 'BB079':'BARBAKC583_0838', 'BB081':'BARBAKC583_0871', 'BB082':'BARBAKC583_0872', 'BB084':'BARBAKC583_0890', 'BB085_1':'BARBAKC583_0894', 'BB085_2':'BARBAKC583_0895', 'BB087':'BARBAKC583_0916', 'BB090':'BARBAKC583_0998', 'BB091':'BARBAKC583_1009', 'BB094':'BARBAKC583_1021', 'BB096':'BARBAKC583_1043', 'BB097':'BARBAKC583_1048', 'BB098':'BARBAKC583_1077', 'BB099':'BARBAKC583_1089', 'BB100':'BARBAKC583_1091', 'BB101_1':'BARBAKC583_1093', 'BB106':'BARBAKC583_1108', 'BB108':'BARBAKC583_1122', 'BB110':'BARBAKC583_1165', 'BB111':'BARBAKC583_1166', 'BB114':'BARBAKC583_1190', 'BB116':'BARBAKC583_1192', 'BB117':'BARBAKC583_1208', 'BB118':'BARBAKC583_1213', 'BB120':'BARBAKC583_1217','BB123':'BARBAKC583_1306',  'BB125_1':'BARBAKC583_1330', 'BB126':'BARBAKC583_1332', 'BB127':'BARBAKC583_1347', 'BB128':'BARBAKC583_1355', 'BB129':'BARBAKC583_1357'}
#I need Cis_Dict here because I need to grab the ranges of the value genes in this Dict

FASTA6 = open('Bartonella_genbank_locus_text_with_sRNAs','r')
Dict3 = {} #ranges of everything including genes and BbsRs
Dict7 = {} #Name of BbsR and overlap length with relevant genes
Dict5 = {} #Lengths of every gene
for line in FASTA6:
	if re.match(r'^\d', line):
		line = line.strip('\n')
		fields5 = re.split('\t', line)
		Name5 = str(fields5[2])
		Start5 = int(fields5[0])
		End5 = int(fields5[1])
		Length = int(End5 - Start5)
		All_Range = range(Start5, End5)	
		Dict3[Name5] = All_Range
		Dict5[Name5] = Length
		
for key in Cis_Dict:
	xs = set(Dict3[key]) #set of BbsR ranges
	Overlap = len(xs.intersection(Dict3[Cis_Dict[key]]))
	Dict7[key] = Overlap #dictionary containing overlap lengths for each individual BbsR

#Now we have average sRNA and Total RNA TPM Files all in the same folder
#So here we are filtering the average BbsRs based on a TPM cut-off (300 based on manual parsing of Artemis, which itself was arbitrary)
#for every Cis-anti sRNA, let's see the TPM contribution of its gene and determine its influence on the sRNA TPM
#We want to pair a gene with each Cis-anti sRNA, print the gene name, look up the gene name in the whole TPM file, get out its TPM, and average that TPM across the biological replicates. We will print that information in the BB_XXXX_TPMs_Filtered file.

FASTA3 = open('END_BB_HUVE_AVERAGE_TPMs', 'r') #sRNA TPMs
FHOUT2 = open('END_BB_HUVE_TPMs_Not_Filtered','w') #sRNA 

#Our final filtered file should have the following header:
#Name	BbsR_TPM	Cis-Anti_Designation	Cis-Gene_Name	Cis-Gene_TPM	BbsR_Length	Gene_Length	Overlap_Length	True_sRNA_TPM

print("BbsR_Name", "BbsR_TPM", "Location", "Cis-Gene_Name", "Cis-Gene_TPM", "Overlap_Length", "True_sRNA_TPM", sep = '\t', file=FHOUT2)
for line in FASTA3:
	if re.match(r'^BB', line):
		line = line.strip('\n')
		fields6 = re.split('\t', line)
		TPM6 = float(fields6[1])
		Name6 = str(fields6[0])
		if Name6 in Cis_Dict:
			for key in Cis_Dict:
				value3 = Cis_Dict[key] #Dict of sRNA:cis-gene pair, so value 3 is the cis-gene
				value4 = Dict4[value3] #Dict of Cis-gene: TPM pair, so value 4 is the gene's TPM
				value5 = Dict5[key] #This is the length of each BbsR
				value6 = Dict5[value3] #This is the length of each relevent gene
				value7 = Dict7[key] #This is the length of the overlap between sRNA and relevent gene
				value8 = (value7 / value6) * value4
				True_TPM = TPM6 - value8
				if key == Name6:
					if True_TPM > 0:
						print(Name6, TPM6, 'Cis-Anti', value3, value4, value7, True_TPM, sep = '\t', file=FHOUT2)
					else:
						print(Name6, TPM6, 'Cis-Anti', value3, value4, value7, "LOW", sep = '\t', file = FHOUT2)
		else:
			if TPM6 > 0:
				print(Name6, TPM6, "IGR/Leader", "N/A", "N/A", "N/A", TPM6, sep = '\t', file=FHOUT2)
			else:
				print(Name6, TPM6, "IGR/Leader", "N/A", "N/A", "N/A", "LOW", sep = '\t', file = FHOUT2)