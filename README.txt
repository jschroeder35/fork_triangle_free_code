WRITTEN AND MAINTENED BY JOSHUA SCHROEDER
for questions/queries, email: joschroed@gmail.com
This  repository contains the code to verify parts of the proof that triangle-
free and fork-free graphs are 3-colorable. The code to set up the classes,
methods, and forbidden subgraph lists to perform the algorithms are contained 
in a module named graph_algs.py and can be imported with the statement

import graph_algs as *

if you are fine with importing names to the global namespace, else

import graph_algs

Notebooks such as "Lem2.4proof.ipynb" contain the code used to verify the proof 
of the corresponding lemma in the paper. Note that all class objects and 
methods have corresponding docstring documentation, so for more detailed 
information about any class or method, you can enter something of the form

help(<name of object/method>)
ex: help(graph.contains_subgraphs)

For now, a few important things to know getting started:

vert: the class of vertices of graphs each vertex has .edges, .non_edges, and 
    .unknown containing the indices of other vertices it has edges, known 
    non-edges, or not known whether there is an edge or not an edge 
    respectively. Vertices can also have .max_degree or .color properties, 
    though they are by default initialized to None. 

graph: the class of graphs, containing the methods to do the graph algorithms 
    and a vertex list comprising the set of vertices in the graph.

Most relevant methods for graph (for most methods, verbose=True/False will 
toggle printed output on vs off, and an output of 'True' indicates a 
contradiction was obtained):

.clean(): will update edge/non_edge/unknown lists so that if a vertex index is
    added to the edge list of one vertex, the corresponding lists of other 
    vertices will be updated appropriately. It is best to run this after adding 
    edges to the graph.

.contains_subgraphs(subgraph_list): searches the graph for any instances of any 
    subgraphs from the subgraph list within the graph. Will return True if a 
    subgraph from the list is present, and False otherwise. 

.update_graph(subgraph_list): will try to do a shallow search for any forced 
    structure implied by the subgraphs in subgraph list being forbidden in the 
    graph, adding any forced edges/non_edges along the way. Will return True if 
    in the process of updating the graph, a contradiction is forced, otherwise 
    will return False and the graph will be updated in place.

.update_graph2_2(subgraph_list, min_subgraph): Heavier graph update 
    function that attempts to search deeper in the space of adding possible 
    edges/non_edges in a smart way. If min_subgraph is specified, it will try 
    to add structure implied by the assumption that the first 
    len(min_subgraph.vertex_list) vertices in the graph are a choice of 
    min_subgraph in the graph such that the number of vertices exactly distance 
    2 is minimal.

.update_graph3(subgraph_list, min_subgraph): Even heavier graph update that 
    will check every unknown vertex pair by creating two temporary graphs, one 
    with the edge added, and one with the non-edge added between the vertex 
    pair, and running .update_graph2_2() to check for a contradiction. The 
    standard default for updating edges/non-edges without adding vertices.

.update_graph4(subgraph_list, min_subgraph): Same as update_graph3, but after 
    each iteration, will attempt to add new vertices to the graph whose 
    existence is implied by a vertex's neighborhood being entirely contained 
    in the neighborhood of another vertex and their unknown list being empty. 
    If min_subgraph is specified, will first attempt to add vertices whose 
    existence is implied by the assumption that the first 
    len(min_subgraph.vertex_list) vertices in the graph are a choice of 
    min_subgraph in the graph such that the number of vertices exactly 
    distance 2 is minimal.

.color_argument_iter(subgraph_list, start_vert, max_distance, min_subgraph
    carg_update_speed): For a partially colored graph starting from start_vert 
    (which is assumed to be unable to be colored), and for each color i, seeks 
    out all possible neighbors that could prevent start_vert from being colored 
    with color i. If this list is empty, then start_vert can be recolored with 
    i; returns True. If this list has length 1, containing only w, the relevant 
    coloring/edges are added to the graph (to prevent the graph from being able 
    to be colored) and (in appropriate cases) it recurses on w. 
    'carg_update_speed' determines which graph udpate function to run while 
    testing potential neighbors of each color (options are 'fast', 'medium', 
    'slow' and run .update_graph2_2(), .update_graph3(), and update_graph4() 
    respectively).

update_graph5(subgraph_list, min_subgraph, carg_update_speed, check_deg_three):
    Alternates between .update_graph4() and the following: sweeps through 
    triples of vertices v1, v2, and v3, and if N(v1) is contained in N(v2) 
    union N(v3), it will run .color_argument_iter() with start_vert as v1 and 
    with v2.color=0, v3.color=1. If check_deg_three=True, it will test if 
    vertices that are currently degree three are forced to be degree 3, and
    also if running .color_argument_iter() with start_vert at that vertex
    under the assumption it is degree 3 produces a contradiction (if it does,
    a new neighbor of that vertex is added to the graph). 
