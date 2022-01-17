#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 13:23:44 2019

@author: sofiastroustrup
"""
import matplotlib.pyplot as plt
import random
from scipy.spatial import distance, distance_matrix

#%%
pop=(-1, 0, 1) #make the list from where to sample randomly

#delete all lists created when running the full script several times
del x_coord
del y_coord
del z_coord
#%%
#make lists
x_coord=list()
y_coord=list()
z_coord=list()

#Loop and sample values from "pop" by random 
#The length of the loop with determine how many coordinates you will end up having in your final walk
#With the following loops we create the datapoints for the x y, and z coordinate

for i in range(1350): 
    r=random.sample(pop, 1) 

    if i==0:
        if r==[-1]:
            x_coord.append(-1)        
        if r==[1]:
            x_coord.append(1)  
        if r==[0]:
            x_coord.append(0)

    else:
        if r==[-1]:
            new=x_coord[i-1]-1
            x_coord.append(new)
        if r==[1]:
            new=x_coord[i-1]+1
            x_coord.append(new)   
        if r==[0]:
            new=x_coord[i-1]+0
            x_coord.append(new) 
            

for i in range(1350): 
    r=random.sample(pop, 1)

    if i==0:
        if r==[-1]:
            y_coord.append(-1)        
        if r==[1]:
            y_coord.append(1) 
        if r==[0]:
            y_coord.append(0)

    else:
        if r==[-1]:
            new=y_coord[i-1]-1
            y_coord.append(new)
        if r==[1]:
            new=y_coord[i-1]+1
            y_coord.append(new)
        if r==[0]:
            new=y_coord[i-1]+0
            y_coord.append(new)


for i in range(1350): 
    r=random.sample(pop, 1)

    if i==0:
        if r==[-1]:
            z_coord.append(-1)        
        if r==[1]:
            z_coord.append(1)  
        if r==[0]:
            z_coord.append(0)

    else:
        if r==[-1]:
            new=z_coord[i-1]-1
            z_coord.append(new)
        if r==[1]:
            new=z_coord[i-1]+1
            z_coord.append(new)
        if r==[0]:
            new=z_coord[i-1]+0
            z_coord.append(new)
            
            
    #%%
#plotting the random walk in 2D

plt.plot(x_coord, y_coord)
plt.title("random walk")
plt.xlim(min(x_coord), max(x_coord))
plt.ylim(min(y_coord), max(y_coord))
plt.savefig("/Users/sofiastroustrup/Desktop/MPI-CBG/Graphs_and_proteins/random_walk_2D.png")

#%%

#plotting random walk 3D

from mpl_toolkits import mplot3d

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(x_coord, y_coord, z_coord, s=2)
ax.set_xlim(min(x_coord), max(x_coord))
ax.set_ylim(min(y_coord), max(y_coord))
ax.set_zlim(min(x_coord), max(z_coord))
plt.savefig("Random walk 3D.png")

#%%
#making adjacency matrix form the 3D random walk

coords_3D =list(zip(x_coord, y_coord, z_coord)) #create a list with coordinates in 3D
dive=round(len(coords_3D)/335) #335 is the  number of nodes we want in the final random walk graph
few_coords_3D=coords_3D[0::dive] #select every X coordinate in the coordinate list

#%%
#Create a dataframe where each row is one datapoint
df=pd.DataFrame(few_coords_3D, columns=['x', 'y', 'z'])

#%%
adjacency=pd.DataFrame(distance_matrix(df.values, df.values)) #make distance matrix
threshold=6 #set threshold for how close two points need to be for creating and edge
adjacency[adjacency==0]=100 #remove the zeroes from the diagonal in the matrix
adjacency[adjacency<=threshold]=1 #all entries equal to or smaller than the treshold is set to 1 = will become an edge
adjacency[adjacency>threshold]=0 #all entries bigger than the threshold, is set to 0
graph = nx.from_pandas_adjacency(adjacency) #generate the graph from the adjacency matrix
#%%
#plot the graph
nx.draw(graph, node_color="skyblue", node_size=50, alpha=0.7, edge_color="grey") 
plt.title("Network representing random walk in 3D") 
plt.savefig("random_walk_null_model.png" )

#%%
import functions_graphs as fg #use robustness function made in another script 
fg.get_metrics(graph,"random walk graph", "", " " )
fg.randrem_edge_percent(graph, 100, "random walk graph4", "")

#%%
fg.TarRem_edges(graph, 70, "random walk graph4", "")

#%%

