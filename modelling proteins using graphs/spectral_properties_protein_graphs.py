#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 10:51:00 2019

@author: sofiastroustrup
"""

import matplotlib.pyplot as plt
import random
import numpy as np
import pandas as pd
import networkx as nx #used to draw the graph
import statistics as stat
import scipy
import functions_graphs as fg
#%%
#import graphs 

trm5=fg.get_graph("5WT1", "A", 332)
nx.draw(trm5, node_color="skyblue", node_size=50, alpha=0.7, edge_color="grey")
plt.title("trm5")
plt.show()

trmD=fg.get_graph("4YVJ", "A", 246, -1)
nx.draw(trm5, node_color="skyblue", node_size=50, alpha=0.7, edge_color="grey")
plt.title("trmD")
plt.show()

random_trmD=nx.gnm_random_graph(239, 1156)
nx.draw(random_trmD, node_color="skyblue", node_size=50, alpha=0.7, edge_color="grey")
plt.title("random graph (trmD no of edges and nodes)")
plt.show()

random_trm5=nx.gnm_random_graph(332, 1707)
nx.draw(random_trm5, node_color="skyblue", node_size=50, alpha=0.7, edge_color="grey")
plt.title("random graph (trm5 no of edges and nodes)")
plt.show()

#%%
adj_trm5=nx.adjacency_matrix(trm5)
eigval_trm5=scipy.linalg.eigh(adj_trm5.todense(), eigvals_only=True)
plt.hist(eigval_trm5, density=0)
plt.title("pdf showing the distribution of eigenvalues \n for the graph representing trm5")
plt.xlabel("eigenvalues")
plt.ylabel("frequency")
plt.savefig("spectral/pdf_eigval_trm5.png")
#%%
adj_trmD=nx.adjacency_matrix(trmD)
eigval_trmD=scipy.linalg.eigh(adj_trmD.todense(), eigvals_only=True)
plt.hist(eigval_trmD)
plt.title("pdf showing the distribution of eigenvalues \n for the graph representing trmD")
plt.xlabel("eigenvalues")
plt.ylabel("frequency")
plt.savefig("spectral/pdf_eigval_trmD.png")
#%%
adj_random_trm5=nx.adjacency_matrix(random_trm5)
eigval_random_trm5=scipy.linalg.eigh(adj_random_trm5.todense(), eigvals_only=True)
plt.hist(eigval_random_trm5)
plt.title("pdf showing the distribution of eigenvalues \n for the graph representing random_trm5")
plt.xlabel("eigenvalues")
plt.ylabel("frequency")
plt.savefig("spectral/pdf_eigval_random_trm5.png")
#%%
adj_random_trmD=nx.adjacency_matrix(random_trmD)
eigval_random_trmD=scipy.linalg.eigh(adj_random_trmD.todense(), eigvals_only=True)
plt.hist(eigval_random_trmD)
plt.title("pdf showing the distribution of eigenvalues \n for the graph representing random_trmD")
plt.xlabel("eigenvalues")
plt.ylabel("frequency")
plt.savefig("spectral/pdf_eigval_random_trmD.png")
#%%


