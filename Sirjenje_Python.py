#Import knjižnic
import matplotlib.pyplot as plt
import networkx as nx
from random import uniform, seed
import numpy as np
import pandas as pd
import time
from collections import Counter


#Create the graph
plt.figure(figsize=[12,10])
H = nx.DiGraph(); #define initial graph
E = [[0,1],[0,2],[0,3],[0,4],[0,5],[1,18],[2,6],[2,7],[3,4],[3,5],[3,7],[3,8],[3,9],[3,10],[3,11],
     [4,3],[4,5],[4,11],[4,12],[4,13],[4,15],[5,3],[5,4],[5,15],[5,16],[5,17],[5,18],[7,14],[7,2],
     [7,8],[8,7],[8,14],[14,7],[14,8],[18,5],[18,19],[18,20],[18,23],[19,20],[19,21],[19,22]]; #edges

w = np.random.uniform(0.3,1,len(E)); #create weights
w = [round(i, 2) for i in w]; #round to 2 decimals

for i in range(len(E)): #add the weights to the edges list to includei n the graph
    E[i].append(w[i]);

H.add_weighted_edges_from(E); #create the initial graph

G = nx.DiGraph() #create an identical graph with lists of nodes and edges that are sorted (important for work with color map)
G.add_nodes_from(sorted(H.nodes(data=True)))
G.add_edges_from(H.edges(data=True))

pos = {22:[0,1], 21:[1,0], 20:[2,0], 19:[1,1], 18:[2,1], 23:[3,0], 17:[4,1], 16:[5,1],
      15:[5,2], 5:[3,2], 1:[2,2], 0:[2,3], 2:[1,4], 6:[0,3], 3:[2,4], 4:[4,3], 13:[5,3], 
      12:[5,4], 11:[4,4], 7:[1,5], 14:[1,6], 8:[2,6], 9:[3,5], 10:[3,6]}; #position

nx.draw(G, pos); #draw the graph
nx.draw_networkx_labels(G, pos);
labels = nx.get_edge_attributes(G, 'weight');
nx.draw_networkx_edge_labels(G,pos,edge_labels=labels);


#Function for the search of new targets and their edges
def propagate_nx(g,new_active):
    
    targets = [] #reset the list of all new targets for this step
    targetEdges = [] #reset all the edges leading towards the new targets for this step
    nodeTargets = [] #reset all the targets of specific node for this step
    for node in new_active: #for every newly activated node
        targets += g.neighbors(node) #add its neighbours to the list of all targets
        nodeTargets = list((Counter(targets) - Counter(nodeTargets)).elements()) #get targets of the current node from the list of all targets with the help of the list of the previous node's targets
        nodeTargets = list(set(nodeTargets)) #remove duplicates from the current node's targets list
        for i in range(len(nodeTargets)): #create a list of the edges leading towards the new target nodes
            newEdge = (node, nodeTargets[i]); #candidate for a new edge to be included into the list
            if (newEdge in list(g.edges)): #if the connection between nodes actually exists
                targetEdges.append((node, nodeTargets[i])) #add it to the list of target edges
        
    return(targets, targetEdges)


#Function of influence propagation
def IC(graph_object,S,w):
    """
    Inputs: graph_object: networkx object
            S:  List of seed nodes
            w:  Influence propagation probability (weights)
    Output: Number of infected nodes, all infected nodes, steps
    """
    
    spread = []
    edges = list(graph_object.edges) #list of all connections in a graph
        
    # Simulate propagation process      
    new_active, A = S[:], S[:]
    steps = []
    while new_active:
                       
        # 1. Find out-neighbors for each newly active node
        targets, targetEdges = propagate_nx(graph_object,new_active) #get new target nodes and edges leading towards them
        
        index = [] #reset the indexes 
        for i in range(len(edges)): #for every edge in the graph
            for j in range(len(targetEdges)): #check which indexes the target edges have in the list of all edges
                if targetEdges[j] == edges[i]:
                    index.append(edges.index(edges[i])) #and save the indexes in a list
                
        p = [w[i] for i in index] #use the indexes to determine the correct weights
    
        # 2. Determine newly activated neighbors from w
        success = np.random.uniform(0,1,len(targets)) < p #check the success of the current step's propagation and save the logical values into a list
        new_ones = list(np.extract(success, sorted(targets))) #use the success list to determine which nodes got activated and save them to a list
            
        # 3. Find newly activated nodes and add to the set of activated nodes
        new_active = list(set(new_ones) - set(A))
        A += new_active
        steps.append(new_active)
            
    spread.append(len(A))
    
    steps.remove([]) #remove the last empty step
    return(np.mean(spread),A,steps) #return the end result: the number of all activated nodes, a list of all activated nodes and the propagation steps


#Results
S = 0
rezultati = IC(G, [S], w)
print('rezultati: ', rezultati)


#Exporting the results in a gexf format for visualization purposes
color_map = ['blue' for node in range(G.number_of_nodes())] #create a "blank" color map (all nodes blue - not active)
for node in G: #add the blue color attributes to the nodes in the graph (for gexf export)
    G.nodes[node]['viz'] = {'color': {'r': 0, 'g': 0, 'b': 255, 'a': 0}}
color_map[S] = 'red' #color the starting node in red
G.nodes[S]['viz'] = {'color': {'r': 255, 'g': 0, 'b': 0, 'a': 0}} #add red color attribute on starting node
nx.draw(G, pos, node_color=color_map, with_labels=True) #draw the graph of step 0
plt.show()
Fname = "korak 0" + '.gexf'
nx.write_gexf(G, Fname) #save step 0 in a gexf format
color_map[S] = 'brown' #after the save color the starting node brown (activated in the previous steps)
G.nodes[S]['viz'] = {'color': {'r': 150, 'g': 75, 'b': 0, 'a': 0}} #add brown color attribute on starting node
for i in range(len(rezultati[2])): #for every step
    for node in G: #check every node to determine
        if(node in rezultati[2][i]): #which nodes are in the current step
            color_map[node] = 'red'; #if a node is in the current step color it red
            G.nodes[node]['viz'] = {'color': {'r': 255, 'g': 0, 'b': 0, 'a': 0}} #add red color attribute to node
    nx.draw(G, pos, node_color=color_map, with_labels=True) #draw the graph
    plt.show()
    Fname = "korak" + str(i+1) + '.gexf'
    nx.write_gexf(G, Fname) #save the step in gexf format
    for j in range(len(color_map)): #for every node
        if(color_map[j] == 'red'): #check if it's red (if it got activated in this step)
            color_map[j] = 'brown'; #and color it brown
            G.nodes[j]['viz'] = {'color': {'r': 150, 'g': 75, 'b': 0, 'a': 0}} #add brown color attribute to node

nx.draw(G, pos, node_color=color_map, with_labels=True) #final state
plt.show() #show the final state of the graph
nx.write_gexf(G, "koncno_stanje.gexf") #save the final state in gexf format