#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 15:36:53 2019

@author: sofiastroustrup
"""
import numpy as np
from Bio import SeqIO
import csv
import pandas as  pd


#%%
#load list with IDs
speclist=pd.read_csv("/Users/sofiastroustrup/Desktop/MPI-CBG/phylogenetics/IDs/speclist.csv", header=None)
speclist.columns=["taxon_4letter", "Kingdom", "NCBI-tax"]

#%%
input_file="pfam_m1G-MT_proteins.fasta"
output_file="taxID_pfam_m1G-MT_proteins.fasta"

#Get taxon ID from fasta file  
with open(input_file) as original, open(output_file, 'w') as corrected:
    records = SeqIO.parse(input_file, 'fasta')
    for record in records:
        ID, _=record.id.split("/")
        uniprot,tax=ID.split("_")
        record.description=tax
        NCBI_tax_row=speclist[speclist['taxon_4letter']==tax]
        NCBI_tax=NCBI_tax_row['NCBI-tax']
        if sum(speclist['taxon_4letter']==tax)==1:
            record.id=str(int(NCBI_tax))
            print(int(NCBI_tax))
        if sum(speclist['taxon_4letter']==tax)==0:
            print("no ref")
            record.id="no ref"
        SeqIO.write(record, corrected, "fasta")
#%%
input_file="pfam_proteins_met10.fasta"
output_file="taxID_pfam_proteins_met10.fasta"

#Get taxon ID from fasta file  
with open(input_file) as original, open(output_file, 'w') as corrected:
    records = SeqIO.parse(input_file, 'fasta')
    for record in records:
        ID, _=record.id.split("/")
        uniprot,tax=ID.split("_")
        record.description=tax
        NCBI_tax_row=speclist[speclist['taxon_4letter']==tax]
        NCBI_tax=NCBI_tax_row['NCBI-tax']
        if sum(speclist['taxon_4letter']==tax)==1:
            record.id=str(int(NCBI_tax))
            print(int(NCBI_tax))
        if sum(speclist['taxon_4letter']==tax)==0:
            print("no ref")
            record.id="no ref"
        SeqIO.write(record, corrected, "fasta")
        
        #%%

# make csv file for mapping the IDs onto the tree

#load taxID
intersection=pd.read_csv("m1G-MT_mapping.tsv", sep=';', header=None)




#%%

input_file="pfam_m1G-MT_proteins.fasta"
output_file="kingdom_pfam_m1G-MT_proteins.fasta"

#Get taxon ID from fasta file  
with open(input_file) as original, open(output_file, 'w') as corrected:
    records = SeqIO.parse(input_file, 'fasta')
    for record in records:
        ID, _=record.id.split("/")
        uniprot,tax=ID.split("_")
        record.description=tax
        NCBI_tax_row=speclist[speclist['taxon_4letter']==tax]
        NCBI_tax=NCBI_tax_row['NCBI-tax']
        kingdom = NCBI_tax_row.loc[NCBI_tax_row['Kingdom']]
        print(str(kingdom))
        if sum(speclist['taxon_4letter']==tax)==1:
            record.id=str(int(NCBI_tax))
            record.name=kingdom
            #print(int(NCBI_tax))
        if sum(speclist['taxon_4letter']==tax)==0:
            #print("no ref")
            record.id="no ref"
        
        #SeqIO.write(record, corrected, "fasta")













