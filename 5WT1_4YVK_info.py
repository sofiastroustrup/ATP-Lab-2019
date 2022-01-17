#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 11:22:16 2019

@author: sofiastroustrup
"""
#Example of how to convert specific pdb structures to graphs
#%%
#Metrics and robustness for the crystalstructures I will be working with and a random graph

trm5=fg.get_graph("5WT1", "A", 332)
nx.draw(trm5, node_color="skyblue", node_size=50, alpha=0.7, edge_color="grey")
plt.title("trm5")
get_metrics(trm5, "trm5", "images/", "5WT1")
randrem_edge(trm5, 100, "Trm5", "images/")

trmD=fg.get_graph("4YVJ", "A", 246, -1)
nx.draw(trm5, node_color="skyblue", node_size=50, alpha=0.7, edge_color="grey")
plt.title("trmD")
get_metrics(trmD, "trmD", "images/", "4YVJ")
randrem_edge(trmD, 100, "TrmD (4YVJ)", "images/")

#%%
#look at a random graph
random_trmD=nx.gnm_random_graph(239, 1707)
nx.draw(random_trmD, node_color="skyblue", node_size=50, alpha=0.7, edge_color="grey")
plt.show()
get_metrics(random_trmD, "random_trmD", "images/", "4YVJ")

random_trm5=nx.gnm_random_graph(332, 1707)
nx.draw(random_trm5, node_color="skyblue", node_size=50, alpha=0.7, edge_color="grey")
plt.show()
get_metrics(random_trm5, "random_trm5", "images/", "5WT1")
TarRem_edges(random_trm5, 50, "random_trm5", "images/")

#%%
nx.graph_edit_distance(trm5, trmD)

#%%

trm5_2=get_graph("3AY0", "A", 336, 2)

trm5_5WT3=get_graph("5WT3", "A", 332, 1)