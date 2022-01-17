#!/bin/bash


grep "UPI0*" $1> uniparc-ID
grep "UPI0*" -v $1 > uniprot-ID
outputname=$2

for ID in $(cat uniparc-ID); 

do 
	curl https://www.uniprot.org/uniparc/$ID.fasta >> $outputname.out

done 

for ID in $(cat uniprot-ID);


do 
	
	curl https://www.uniprot.org/uniprot/$ID.fasta >> $outputname.out
done 

