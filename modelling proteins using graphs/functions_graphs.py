#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 13:38:48 2019

@author: sofiastroustrup
"""

import matplotlib.pyplot as plt
import random
import numpy as np
import pandas as pd
import networkx as nx #used to draw the graph
import statistics as stat


from evcouplings.compare import (
    PDB, DistanceMap, SIFTS, intra_dists,
    multimer_dists, coupling_scores_compared) #used to obtain distance matrices 

from evcouplings.visualize import plot_contact_map, plot_context

#%%

def get_graph(PDB_ID, chain, filter_end, filter_start=0, threshold=5):
    """Make a graph based on the pdb structure. Uses the evcouplings packages to calculate distance matrices. \n
    This function assign distance 100 to the diagonal in the adjacency matrix. filter_end: should be exactly the end \n
    position of your chain of interest in the pdb file."""
    #coords=PDB.from_id(PDB_ID).get_chain(chain).filter_atoms() can be used in more general code to only get the c-alfa
    coords=PDB.from_id(PDB_ID).get_chain(chain).filter_positions(range(filter_start, filter_end+1))
    distmap=DistanceMap.from_coords(coords)
    adjacency=pd.DataFrame.from_records(distmap.dist_matrix)
    adjacency[adjacency==0]=100
    adjacency[adjacency<threshold]=1
    adjacency[adjacency>threshold]=0
    graph = nx.from_pandas_adjacency(adjacency)
    return graph

#%%

def get_metrics(graph, name, path, PDBid):
    """uses graph to compute and plot standard metrics. Name is used for standard title\
    and path is used for saving the image to specific folder"""
    degrees=list(nx.degree(graph))
    pos,deg=zip(*degrees)
    min_deg=min(deg)
    edges=len(graph.edges)
    nodes=len(graph.nodes)
    print("Number of edges is {} ".format(edges))
    print("Number of nodes is {} ".format(nodes))
    print("Minimum degree is {}".format(min_deg))
    max_deg=max(deg)
    print("Maximum degree is {}".format(max_deg))
    median_deg=stat.median(sorted(deg))
    print("Median degree is {}".format(median_deg))
    #mode_deg=stat.mode(deg)
    #print("Mode degree is {}".format(mode_deg))
    
    fig, axes = plt.subplots(3, 2, figsize=(15,15))
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    plt.suptitle("Metrics for {} ({}) \n total no of edges is {} and total no of nodes is {}".format(name, PDBid, edges, nodes), size=16)
    ax0, ax1, ax2, ax3, ax4, ax5 = axes.flatten()
    
    ax0.hist(deg, density=0)
    ax0.set_title("Degree distribution for {}".format(name))
    ax0.set_xlabel("Degrees")
    ax0.set_ylabel("Frequency")
    #plt.savefig("{}{}_degree_distribution.png".format(path, name))
    #plt.show()
    
    ax1.hist(deg, density=0, cumulative=True, histtype='step')
    ax1.set_title("Cumulative degree distribution for {}".format(name))
    ax1.set_xlabel("Degrees")
    ax1.set_ylabel("Cumulative frequency")
    #plt.show()
    
    betweenness_edge=nx.edge_betweenness_centrality(graph)
    ax2.hist(betweenness_edge.values())
    ax2.set_title("Distribution of betweenness centrality for edges {}".format(name))
    ax2.set_xlabel("Betweenness centrality")
    ax2.set_ylabel("Frequency")
    #plt.show()
    
    betweenness_node=nx.betweenness_centrality(graph)
    ax3.hist(betweenness_node.values())
    ax3.set_title("Distribution of betweenness centrality for nodes {}".format(name))
    ax3.set_xlabel("Betweenness centrality")
    ax3.set_ylabel("Frequency")
    #plt.show()
    
    closeness=nx.closeness_centrality(graph)
    ax4.hist(closeness.values())
    ax4.set_title("Distribution of closeness centrality for nodes {}".format(name))
    ax4.set_xlabel("Closeness centrality")
    ax4.set_ylabel("Frequency")
    #plt.show()
    
    ax5.hist(closeness.values(), cumulative=True, histtype="step", density=True)
    ax5.set_title("Cumulative closeness centrality for nodes {}".format(name))
    ax5.set_ylabel("Relative frequency")
    ax5.set_xlabel("Cumulative closeness centrality")
    
    plt.savefig("{}metrics_{}".format(path, name))


#%%
    
def randrem_edge(graph, repeats, name, path):
    """function that plots the largest remaining component as a function of removed edges. repeats=how many different sequences of removed edges/n
    name=name for the protein for title, path = absolute or relative path for were to store the output file"""
    G=graph.copy()
    no_edges=len(graph.edges())
    giant_component=[[] for i in range(no_edges)]    
    for j in range(0, repeats):
        #if j>0: 
            #del rand_index
        G=graph.copy()
        no_edges=len(graph.edges)
        rand_index=random.sample(range(0, no_edges), no_edges)
        for i in range(0, no_edges):
            u, v=list(graph.edges)[rand_index[i]]
            G.remove_edge(u,v)
            giant_component[i]=len(max(nx.connected_components(G), key=len))
        x_coord=list(range(no_edges, 0, -1))
        plot=plt.plot(x_coord, giant_component, ".")
        plt.title("Random removal of edges for {}, no of graphs {}".format(name, repeats))
        plt.xlabel("no. of edges left \n total no. of edges {} ".format(no_edges))
        plt.ylabel("largest remaining component (nodes)")
        plt.savefig("{}random_removal_edges_{}".format(path, name))
    plt.show()
    
    return plot
#%%

def TarRem_edges(graph, percentage_edges_remove, name, path):
    """removes edges by betweenness, edges with the highest betweenness are removed first"""
    G = graph.copy()
    no_edges_remove = round((len(G.edges)/100)*percentage_edges_remove)
    no_edges=len(G.edges)
    incidence = no_edges-no_edges_remove
    connected_component=[[] for i in range(0, no_edges_remove)]
    
    for i in range(0, no_edges_remove):
        betweenness=nx.edge_betweenness_centrality(G)
        values_betweenness = list(betweenness.values())
        keys_betweenness = list(betweenness.keys())
        max_index=values_betweenness.index(max(values_betweenness))
        u, v=keys_betweenness[max_index]
        G.remove_edge(u, v)
        connected_component[i]=len(max(nx.connected_components(G), key=len))
        
    #connected_component
    x_coord = list(range(no_edges, incidence, -1))
    plt.plot(x_coord, connected_component, '.')
    plt.title("targeted removal of {} edges ({}%) by highest betweenness for {}".format(no_edges_remove, percentage_edges_remove, name))
    plt.xlabel("no of edges left (total no of edges {})".format(no_edges))
    plt.ylabel("largest remaining component (nodes)")
    plt.savefig("{}{}_targeted_removal_betweeness{}%.png".format(path, name, percentage_edges_remove))
    plt.show()

#%%
    
def TarRem_edges_percent(graph, percentage_edges_remove, name, path):
    """removes edges by betweenness, edges with the highest betweenness are removed first"""
    G = graph.copy()
    no_edges_remove = round((len(G.edges)/100)*percentage_edges_remove)
    no_edges=len(G.edges)
    incidence = no_edges-no_edges_remove
    connected_component=[[] for i in range(0, no_edges_remove)]
    
    for i in range(0, no_edges_remove):
        betweenness=nx.edge_betweenness_centrality(G)
        values_betweenness = list(betweenness.values())
        keys_betweenness = list(betweenness.keys())
        max_index=values_betweenness.index(max(values_betweenness))
        u, v=keys_betweenness[max_index]
        G.remove_edge(u, v)
        connected_component[i]=len(max(nx.connected_components(G), key=len))
        
    #connected_component
    coords=np.linspace(0, 100, no_edges)
    incidence=coords[-no_edges_remove:]
    incidence=incidence[::-1]
   #x_coord = round(list(range(100, incidence, -stepsize)))
    #print(x_coord)
    #x_coord = list(range(no_edges, incidence, -1))
    plt.plot(incidence, connected_component, '.')
    plt.title("targeted removal of {} edges ({}%) by highest betweenness for {}".format(no_edges_remove, percentage_edges_remove, name))
    plt.xlabel("no of edges left in % (total no of edges {})".format(no_edges))
    plt.ylabel("largest remaining component (nodes)")
    plt.savefig("{}{}_targeted_removal_x-axis-percent_betweeness{}%.png".format(path, name, percentage_edges_remove))
    plt.show()    
    
    
#%%

def randrem_edge_percent(graph, repeats, name, path):
    """function that plots the largest remaining component as a function of removed edges. repeats=how many different sequences of removed edges/n
    name=name for the protein for title, path = absolute or relative path for were to store the output file"""
    G=graph.copy()
    no_edges=len(graph.edges())
    giant_component=[[] for i in range(no_edges)]    
    for j in range(0, repeats):
        #if j>0: 
            #del rand_index
        G=graph.copy()
        no_edges=len(graph.edges)
        rand_index=random.sample(range(0, no_edges), no_edges)
        for i in range(0, no_edges):
            u, v=list(graph.edges)[rand_index[i]]
            G.remove_edge(u,v)
            giant_component[i]=len(max(nx.connected_components(G), key=len))
        coords=np.linspace(0, 100, no_edges)
        x_coords=coords[::-1]
        #x_coord=list(range(no_edges, 0, -1))
        plot=plt.plot(x_coords, giant_component, ".")
        plt.title("Random removal of edges for {}, no of graphs {}".format(name, repeats))
        plt.xlabel("percent of edges left \n total no. of edges {} ".format(no_edges))
        plt.ylabel("largest remaining component (nodes)")
        plt.savefig("{}random_removal_edges_x-axis-percent_{}_repeats_{}".format(path, name, repeats))
    plt.show()
    
    return plot