import networkx as nx

######################### graphTools.py - helpful functions for the graph building algorithm ##########################

### def get_node(G, attribute, val)
##### graph G - the graph to search
##### string attribute - the name of the attribute you are searching
##### string/int/value val - the value of the atribute you are searching for
### returns:
##### the ID of the node with the matching attribute value
##### or None if the value is not found
def get_node(G, attribute, val):
    for node in G.nodes(data=True):
        if node[1][attribute]==val:
            return node[0]
    print("No such node found")
    return None

### def get_edge(G, node1, nod2)
##### graph G - the graph to search
##### string node1 - the ID of the first node
##### string node2- the Id of the second node
### returns:
##### The edge that connects node1 to node2
##### or None if the value is not found
def get_edge(G, node1, node2):
    for edge in G.edges(data=True):
        if edge[0]==node1 and edge[1]==node2:
            return edge
    print("No such edge found")
    return None