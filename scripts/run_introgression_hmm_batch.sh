#!/bin/bash

#Matt Gibson
#Aug. 2019
#Indiana University Departement of Biology
#Moyle Lab

#Make batch calls to hmm

#Arg 1: Directory to input probabilities (generated with code from hmm_calculate_probabilities.R)
#Arg 2: Prefix to search for in directory (e.g., to run HMM on all chromosomes of MG114.5, arg 2 will be MG114.5)

#Usage:
#./run_introgression_hmm_batch.sh ./hmm_data/ MG114.5

pref=$2

cd $1

for f in *$pref*
do
	python ../introgression_hmm.py $f $f.tsv 1
done

