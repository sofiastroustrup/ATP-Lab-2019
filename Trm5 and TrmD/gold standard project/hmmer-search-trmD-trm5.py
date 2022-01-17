#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 15:26:07 2019

@author: sofiastroustrup
"""

import csv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#%%
#Loading the data and preparing the dataset
trmD_hmm_search=pd.read_csv('trmDprofile-against-trmD-swissprot.hmmsearch.csv', sep=',', header=None)
#trmD_hmm_search.drop(trmD_hmm_search.columns[[3, 4, 5, 6, 7]], axis=1, inplace=True) #remove the columns we are not interested in
trmD_hmm_search.columns = ["E-val full seq", "score full seq", "bias", "sequence", "Description"]

#%%
#Normalization factor was read from the hmm profile: clustalo-swissprot-trmD-19-07-30.cln.msa.hmmbuild
norm1=245
trmD_hmm_search['norm score full seq']=trmD_hmm_search['score full seq']/norm1
#%%
#plot histogram of trmD-profile alignment (from cleaned sequences) to trmD swissprot (uncleaned)
plt.hist(trmD_hmm_search['norm score full seq'])
plt.title('Gold standad trmD distribution of scores \n total number of sequences 641')
plt.ylabel("frequency")
plt.xlabel("normalized full sequence score")
plt.savefig("Goldstd_trmD_distribution.png")

#%%
#select column based on the value in another column
trmD_hmm_search.loc[trmD_hmm_search['norm score full seq']<1, 'sequence']
#drop rows that dont have trmD
trmD_hmm_search=trmD_hmm_search.drop(list(range(642, 649)), axis=0)

#%%
#plot the density distribution within goldstandard trm5
#Normalization factor was read from the hmm profile: clustalo-swissprot-trm5-19-07-30.cln.msa.hmmbuild
norm2=431
#load data 
trm5_hmm_search=pd.read_csv('trm5profile-against-trm5-swissprot.hmmsearch.csv', sep=',', header=None)
trm5_hmm_search.columns = ["E-val full seq", "score full seq", "bias", "sequence", "Description"]
#%%
# add column with normalized the scores 
trm5_hmm_search['norm score full seq']=trm5_hmm_search['score full seq']/norm2


#how do the lowest scores look?
trm5_hmm_search.loc[trm5_hmm_search['norm score full seq']<0.6, 'sequence']
#%%
#remove the last row because it is not a trm5
trm5_hmm_search=trm5_hmm_search.drop(61, axis=0)

#%%
#plot data normalized
plt.hist(trm5_hmm_search['norm score full seq'])
plt.title('Gold standad trm5 distribution of scores \n total number of sequences {}'.format(len(trm5_score)))
plt.xlabel("Normalized score")
plt.ylabel("Frequency")
plt.savefig("Goldstd_trm5_distribution.png")

#%%

#Plot the trm5profile against trmD jackhmmer hits 
#Import data

trm5prof_trmDdata=pd.read_csv("trm5profile-against_cdhit-0.7-jackhmmer-ecoli-E0.01-uniref100-bitscore.hmmsearch.csv2", sep=',', header=None)
#add column names 
trm5prof_trmDdata.columns = ["E-val full seq", "score full seq", "bias", "sequence", "Description"]
#Normalization factor was read from the hmm profile: clustalo-swissprot-trm5-19-07-30.cln.msa.hmmbuild
norm2=431
#add normalized column to dataframe 
trm5prof_trmDdata['norm score full seq']=trm5prof_trmDdata['score full seq']/norm2

#plot
plt.hist(trm5prof_trmDdata['norm score full seq'])
plt.title('Trm5 HMM profile against the trmD jackhmmer hits \n number of reported sequences 109 out of 4645')
plt.xlabel('normalized full sequence score')
plt.ylabel('frequency')
plt.savefig('Trm5profile_against_trmD_jackhmmer_hits.png')

#%%
#Plot trmD profile against trm5 jackhmmer hits 

#load data 

trmDprof_trm5data= pd.read_csv('trmDprofile-against_cdhit-0.7-jackhmmer-trm5-E0.0000001-uniref100-bitscore.hmmsearch.csv1', sep=',')
trmDprof_trm5data.columns = ["E-val full seq", "score full seq", "bias", "sequence", "Description"]

#Normalization factor was read from the hmm profile: clustalo-swissprot-trmD-19-07-30.cln.msa.hmmbuild
norm1=245
#Add normalized column to the dataframe 
trmDprof_trm5data['norm score full seq']=trmDprof_trm5data['score full seq']/norm1

#plot histogram
plt.hist(trmDprof_trm5data['norm score full seq'])
plt.title('TrmD profile against trm5 jackhmmer hits \n number of reported sequences 1199 out of 34774')
plt.xlabel('normalized score full sequence')
plt.ylabel('frequency')
plt.savefig('TrmDprofile_against_trm5_jackhmmer_hits.png')

#%%

trmDprof_trm5data.loc[trmDprof_trm5data['norm score full seq']>0.01, 'Description']