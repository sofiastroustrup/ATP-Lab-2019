#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 13:42:22 2019

@author: sofiastroustrup
"""

#generate graphs at different threshold
AY_list=[[] for i in range(0, 10)]
for i in range(0, 10):
    no=i+1
    #print(i)
    #print(no)
    AY=get_graph("3AY0", "A", no)
    AY_list[i]=AY


#plotting the graphs depending on Å
for i in range(7):
    edges=len(AY_list[i].edges)
    nx.draw(AY_list[i], node_color="skyblue", node_size=50, alpha=0.7, edge_color="grey")
    plt.title("graph with threshold {} Å \n total no. of edges {} ".format(i+1, edges))
    plt.show()
    #edges=len(AY_list[i].edges)
    #print("no. of edges {} \n total no. of edges {}.format(edges))
    degrees=list(nx.degree(AY_list[i]))
    pos,deg=zip(*degrees)

    
    plt.hist(deg, density=0)
    plt.title("Degree distribution for graph with threshold {} Å \n total no. of edges {}".format(i+1, edges))
    plt.xlabel("Degree")
    plt.ylabel("Frequency")
    plt.show()


#grid = plt.GridSpec(2, 1, wspace=0.7, hspace=0.4)

degrees=list(nx.degree(AY_list[0]))
plt.subplot(grid[0,0])
plt.hist(deg, density=0)
plt.title("Degree distribution for graph with threshold {} Å \n total no. of edges {}".format(1, edges))
plt.xlabel("Degree")
plt.ylabel("Frequency")

degrees=list(nx.degree(AY_list[3]))
plt.subplot(grid[0,1:])
plt.hist(deg, density=0)
plt.title("Degree distribution for graph with threshold {} Å \n total no. of edges {}".format(1, edges))
plt.xlabel("Degree")
plt.ylabel("Frequency")