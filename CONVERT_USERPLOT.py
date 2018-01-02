#!/usr/bin/env python3

#Take the output of Nesoni ambiguous/unambiguous userplots and convert them to CSV format for TPM.pl script...

FASTA3 = open('Vero_SCV_2_NM2_Ambiguous.userplot','r')
FHOUT3 = open('Converted_Vero_SCV_2_NM2_Ambiguous.txt','w')

N=1

for line in FASTA3:
	line = line.strip('\n')
	print(N,line,sep=',', file = FHOUT3)
	N = N+1
	continue
