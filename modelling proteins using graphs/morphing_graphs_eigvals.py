#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 13:41:00 2019

@author: sofiastroustrup
"""

#morphing script 
import random
import networkx as nx
import matplotlib.pyplot as plt
import random
import numpy as np
import pandas as pd
import networkx as nx #used to draw the graph
import statistics as stat
import scipy
import itertools

#%%
Graph1=nx.dense_gnm_random_graph(13,20, seed=10)
Graph2=nx.dense_gnm_random_graph(20,35, seed=2)

#%%

Graph1=nx.dense_gnm_random_graph(20, 20, seed=10)
Graph2=nx.dense_gnm_random_graph(20,25, seed=11)

#%%

def eigval_graph_morphing(graph1, graph2, no_best_graph=3, changes=20):
    G1=graph1.copy()
    G2=graph2.copy()
    start_dist=dist_eig(G1, G2)
    start_graph=G1
    #print("start edit distance is: {}".format(nx.graph_edit_distance(G1, G2)))
    s=0
    t=0
    i=0
    while i==0 or new_dist<=start_dist:
        if i==0:
            print("i={}".format(i))
            print("start max eigenvalue distance = {}".format(start_dist))
            print("The number of best graphs to use = {}".format(no_best_graph))
            list_with_two_edit_graphs=pos_graphs2(pos_graphs1(G1, changes), changes)
            flat_two_edit_graphs = []
            for sublist in list_with_two_edit_graphs:
                for item in sublist:
                    flat_two_edit_graphs.append(item)
                
            total_dist=[[]for j in range(len(flat_two_edit_graphs))]
            for p in range(len(flat_two_edit_graphs)):
                total_dist[p]=dist_eig(flat_two_edit_graphs[p], G2)

            sorted_unique_total_dist=sorted(np.unique(total_dist))
        
            indices = [total_dist.index(sorted_unique_total_dist[p]) for p in range(len(sorted_unique_total_dist))] #get indices for the 3 graphs with the lowest different eigenvalue
            new_dist = sorted_unique_total_dist[0]
            print("new distance is {}".format(new_dist))
            #print("new edit distance is {}".format(nx.graph_edit_distance(flat_two_edit_graphs[indices[0]], G2)))
        
            pass
        
        if i > 0:
            print("i={}".format(i))
            print("the number of changes in pos_graphs1 {}".format(changes))
            print("The number of best graphs to use = {}".format(no_best_graph))
            start_dist = sorted_unique_total_dist[0]
            print("start max eigenvalue distance = {}".format(start_dist))
            start_graph = flat_two_edit_graphs[indices[0]]
            
            if no_best_graph>len(indices):
                print("Number of best graphs = {} is too high, max is {} and is used instead".format(no_best_graph, len(indices)))
                no_best_graph=len(indices)
            
            G1=[flat_two_edit_graphs[indices[i]] for i in range(no_best_graph)] #put the best 3 graphs into a list
            #print(nx.graph_edit_distance(G1[0], G2))
            
            list_with_two_edit_graphs=[[] for j in range(len(G1))]
            for g in range(len(G1)):
                list_with_two_edit_graphs[g]=pos_graphs2(pos_graphs1(G1[g], changes), changes)
            
            flat_two_edit_graphs = []
            for sublist in list_with_two_edit_graphs:
                for subsublist in sublist:
                    for item in subsublist:
                        flat_two_edit_graphs.append(item)
                
            total_dist=[[]for j in range(len(flat_two_edit_graphs))]
            for p in range(len(flat_two_edit_graphs)):
                total_dist[p]=dist_eig(flat_two_edit_graphs[p], G2)
                
            sorted_unique_total_dist=sorted(np.unique(total_dist)) #sorted_unique_total_dist=sorted(np.unique(total_dist))
        
            indices = [total_dist.index(sorted_unique_total_dist[p]) for p in range(len(sorted_unique_total_dist))] #get indices for the 3 graphs with the lowest different eigenvalue
            #print(sorted_unique_total_dist)
            new_dist = sorted_unique_total_dist[0]
            print("new eigenvalue distance is {}".format(new_dist))
            #print("new diameter is {}".format(nx.algorithms.distance_measures.diameter()))
    
            pass
        
        #if i!=0 and int(i/5)-(i/5)==0:
         #   print("5 rounds")
          #  return start_graph
        
        i=i+1
        
        stepsize=abs(new_dist-start_dist)
        
        if i==10:
            print("i==10")
            return start_graph
        
        if stepsize==0.1:
            no_best_graph=20
        
        if start_dist < new_dist:
            print("New distance is greater than old distance, optimization has been stopped")
            print("G1 has been turned into G2: {}".format(nx.is_isomorphic(start_graph, G2)))
            s=s+1
            changes=100
            no_best_graph=20
            print("s={}".format(s))
            if s==2:
                print("hey")
                return start_graph
                break 
            else:
                pass
        
        if start_dist == new_dist:
            print("Start distance equals new distance")
            print("G1 has been turned into G2: {}".format(nx.is_isomorphic(start_graph, G2)))
            changes=100
            no_best_graph=20
            t=t+1
            print("t={}".format(t))
            if t==2: 
                return start_graph
                break
            else:
                pass

    
    
    #%%
    
def pos_graphs1(graph, changes):
    numbers=range(1,5)
    g1=graph.copy()
    G1_nodes=list(g1.nodes())
    G1_edges=list(g1.edges())
    len_edges=range(len(G1_edges))
    len_nodes=range(len(G1_nodes)) 
    
    which_changes = np.random.choice(numbers, changes)
    if sum(which_changes==1)<len(g1.edges):
        no_remove_edge=sum(which_changes==1)
    else:
        no_remove_edge=len(g1.edges)
    if sum(which_changes==2)<len(g1.nodes):
        no_remove_node=sum(which_changes==2)
    else:
        no_remove_node=len(g1.nodes)
        
    no_add_edge=sum(which_changes==3)
      
    removed_edges_new_graph=[[] for i in range(no_remove_edge)]
    which_edges_remove=np.random.choice(len_edges, no_remove_edge, replace=False)
    #print(which_edges_remove)
    for j in range(no_remove_edge):
        G_cur = g1.copy()
        u,v=G1_edges[which_edges_remove[j]]
        G_cur.remove_edge(u,v)
        removed_edges_new_graph[j]=G_cur
    
    removed_nodes_new_graph=[[] for i in range(no_remove_node)]
    which_nodes_remove=np.random.choice(len_nodes, no_remove_node, replace=False)
    for j in range(no_remove_node):
        G_cur = g1.copy()
        G_cur.remove_node(G1_nodes[which_nodes_remove[j]])
        removed_nodes_new_graph[j]=G_cur

    #add nodes 
    no_add_node=5
    added_nodes_new_graph=[[] for i in range(no_add_node)]
    #g1=nx.Graph()
    for j in range(1, no_add_node+1):
        G_cur = g1.copy()
        if j==1: 
            which_node=len(G_cur.nodes)+1
            G_cur.add_node(which_node)
        if j>1:
            which_nodes=range(len(G_cur.nodes)+1,len(G_cur.nodes)+1+j)
            G_cur.add_nodes_from(which_nodes)
        added_nodes_new_graph[j-1]=G_cur
        
    #add edges 
    pos_com=list(itertools.combinations(G1_nodes, 2)) #code that samples all possible pair from the 10 nodes 
    for p in range(len(G1_edges)):
       edge=G1_edges[p]
       pos_com.remove(edge)
    edges_to_add=pos_com
       
    #print("no edges that can be added {}".format(edges_to_add))
    if no_add_edge<=len(edges_to_add):
        added_edges=[[] for p in range(no_add_edge)]
        which_edges_add=np.random.choice(len(edges_to_add), no_add_edge, replace=False)
        pass 

    if no_add_edge>len(edges_to_add):
        added_edges=[[] for p in range(len(edges_to_add))]
        which_edges_add=range(len(edges_to_add))
        no_add_edge=len(edges_to_add)
        pass    
    
    for p in range(no_add_edge):
       G_cur=g1.copy()
       u, v=edges_to_add[which_edges_add[p]]
       G_cur.add_edge(u,v)
       added_edges[p]=G_cur
    
    all_new_G = list()
    all_new_G.extend(removed_edges_new_graph)
    all_new_G.extend(removed_nodes_new_graph)
    all_new_G.extend(added_nodes_new_graph)
    all_new_G.extend(added_edges)
    
    #print("no of removed edges {}, no of removed nodes {}, no of added edges {}".format(len(removed_edges_new_graph), len(removed_nodes_new_graph), len(added_edges)))
    
    return all_new_G
    
#%%
    
def pos_graphs2(list_of_graphs, changes):
    new_list_of_graphs= [[] for j in range(len(list_of_graphs))]
    for i in range(len(list_of_graphs)):
        new_list_of_graphs[i]=pos_graphs1(list_of_graphs[i], changes)
    return new_list_of_graphs
    
    
    #%%
    
def dist_eig(g1, g2):
    eigval1=max(abs(scipy.linalg.eigh(nx.adjacency_matrix(g1).todense(), eigvals_only=True)))
    eigval2=max(abs(scipy.linalg.eigh(nx.adjacency_matrix(g2).todense(), eigvals_only=True)))
    dist=abs(eigval1 - eigval2)
    return dist

#%%      

def dist_eig(g1, g2):
   eigval1=sorted(abs(scipy.linalg.eigh(nx.adjacency_matrix(g1).todense(), eigvals_only=True)))
   eigval2=sorted(abs(scipy.linalg.eigh(nx.adjacency_matrix(g2).todense(), eigvals_only=True)))
   if len(eigval1)==len(eigval2):
       dist=scipy.spatial.distance.euclidean(eigval1, eigval2)
       pass
   else: 
       minimum=min([len(eigval1), len(eigval2)])
       dist=scipy.spatial.distance.euclidean(eigval1[:minimum], eigval2[:minimum])+abs(len(eigval1)-len(eigval2))
   return dist   

#%%
   
