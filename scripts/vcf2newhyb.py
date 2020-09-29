#!/usr/bin/env python3

"""
Matt Gibson
Aug. 2019
Indiana University Departement of Biology
Moyle Lab

Converts vcf to NewHybrids format

Arg 1: input vcf file
Arg 2: input popmap file
Arg 3: output file name
"""


import sys

input_vcf = sys.argv[1]
popmap = sys.argv[2]
output = sys.argv[3]

vcf = open(input_vcf, "r")
popmap = open(popmap, "r")
output = open(output, 'w')

pops = {}
names = []
for i, line in enumerate(popmap):
	l = line.replace("\n", "").split()
	indv = l[0]
	idd = l[1]
	names.append(indv)
	if idd in pops.keys():
		pops[idd].append(indv)
	else:
		pops[idd] = [indv]
popmap.close()


print(pops)

indv_idxs = {}

genotypes = {}

loci = []
loci_markers = []

for i, line in enumerate(vcf):
	l = line.replace("\n", "").split()
	if line.startswith("#CHROM"):
		#at header
		indvs = l[9:]
		for x in names:
			idx = indvs.index(x)
			indv_idxs[x] = idx
	elif line.startswith("SL3.0"):
		#genotypes
		locus = l[2]
		marker = locus.split(":")[0]
		calls = l[9:]
		if marker not in loci_markers:
			loci.append(locus)
			loci_markers.append(marker)
			for key, value in indv_idxs.items():
				gt = calls[value].split(":")[0]
				if key in genotypes.keys():
					genotypes[key].append(gt)
				else:
					genotypes[key] = [gt]

		
vcf.close()


##OUTPUT

def convert_to_int_codes(gts):
	new = []
	for call in gts:
		if call == '0/0':
			new.append('11')
		elif call == '0/1' or call == '1/0':
			new.append('12')
		elif call == '1/1':
			new.append('22')
		elif call == './.':
			new.append('0')
	return(new)

def which_cat(indv, dict):
	for key, val in dict.items():
		if indv in val:
			return(key)

output.write("NumIndivs " + str(len(names)) + '\n')
output.write("NumLoci " + str(len(genotypes['MG120-10'])) + '\n')
output.write("Digits 1\n")
output.write("Format lumped\n\n")

output.write("LocusNames " + " ".join(loci) + '\n')

count = 1
for key, val in genotypes.items():
	val = convert_to_int_codes(val)
	cat = which_cat(key, pops)
	if cat == 'ad':
		output.write(str(count) + ' ' + " ".join(val) + '\n')
	else:
		output.write(str(count) + ' ' + cat + ' ' + " ".join(val) + '\n')
	count += 1

