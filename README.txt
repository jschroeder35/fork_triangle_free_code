This  repository contains the code to verify parts of the proof that triangle-free and fork-free graphs are 3-colorable. The code to set up the classes, methods, and forbidden subgraph lists to perform the algorithms are now contained in a module named graph_algs.py and can be imported with the statement

import graph_algs as *

graph_algs will always have the most up to date version of the code going forward, but an older version with some now-deprecated methods is saved in test_build6.0.ipynb for referencing purposes. That notebook's code, however, has not been kept up to date, contains many verified bugs and should not be used. Notebooks such as "Lem2.4proof.ipynb" contain the code used to verify the proof of the corresponding lemma in the paper. 

This README file may be expanded with further documentation and explanation of the methods at a later date. For now, a few quick things to know:

vert: the class of vertices of graphs each vertex has .edges, .non_edges, and .unknown containing the indices of other vertices it has edges, known non-edges, or not known whether there is an edge or not an edge respectively. Vertices can also have .max_degree or .color properties, though they are by default initialized to None. 

graph: the class of graphs, containing the methods to do the graph algorithms. contains a vertex list containing all the vertices of the graph, as well as a list of indices and a hidden property .depth for computational purposes. 

Most relevant methods for graph:

.clean(): will update edge/non_edge/unknown lists so that if a vertex index is added to the edge list of one vertex, the corresponding lists of other vertices will be updated appropriately. It is best to run this after adding edges to the graph.

.contains_subgraphs(subgraph_list): searches the graph for any instances of any subgraphs from the subgraph list within the graph. Will return True if a subgraph from the list is present, and False otherwise

.update_graph(subgraph_list): will try to do a shallow search for any forced structure implied by the subgraphs in subgraph list being forbidden in the graph, adding any forced edges/non_edges along the way. Will return True if in the process of updating the graph, a contradiction is forced, otherwise will return False and the graph will be updated in place.

.update_graph2(subgraph_list, max_depth): Heavier and slower version of the update_graph function that does a deeper search for forced structure, but can solve more cases. Return type is the same. max_depth indicates the maximum recursion to which a type of logic is allowed to run for this method, though 1 is generally.

.update_graph3(subgraph_list): Even heavier graph update that will check every unknown vertex pair by creating two temporary graphs, one with the edge added, and one with the non-edge added between the vertex pair, and running update_graph2 to check for a contradiction. Often can find more structure than update_graph2 but may take >30mins to run for graphs with many unknown vertex pairs.

.update_graph4(subgraph_list): Same as update_graph3, but after each iteration, will attempt to add vertices to the graph by searching for vertices whose neighborhood is entirely contained in the neighborhood of another vertex and whose unknown list is empty. If it finds such a vertex pair, it will add a new neighbor for that vertex not contained in the neighborhood of the other vertex, then run again.

.color_argument2_iter(subgraph_list, start_vert, max_distance): An expansion upon update_graph2(), this method also uses structure that is forced by partial coloring in the graph to try to build out the graph and seek a contradiction. subgraph_list is the set of forbidden subgraphs, start_vert is the single vertex that we assume is forced to be unable to be colored within the graph, and max_distance is the bound on how far the algorithm will recurse away from start_vert within the graph to try to obtain a contradiction. It is advised that runtime grows very very fast as a function of max_distance. 




[TODO: update this algorithm explanation to factor in changes/updates]

Here's a breakdown of the color_argument2 algorithm (return True means contradiction):

takes in: current vertex v_0, done_list, current_distance, max_distance, max_depth
Make color_list = (0, 1, 2) - color(v_0) (if v_0) is colored

If current_distance > 1 and you've looped back to the start, return False
Alternatively, if current_distance >= max_distance, return False (do not proceed). 

Supposing neither of these is the case, for each color c_0 in color_list:

    make 5 lists:

        list0: Check all the neighbors of v_0. If v_1 in N(v_0) is colored with c_0, run color_argument2(v_1, done_list + v_0, current_distance +1, max_distance, max_depth). If False (no contradiction) then append v_1 to list0. 

        list1: Check all the neighbors of v_0 again. If v_1 in N(v_0) is colored with 'None', check if it has any neighbors colored with c_0. If not, make a copy of the graph G', and within G' color v_1 with c_0. Run the methods that extend colorings appropriately through the graph, and check for forks or triangles. If these methods indicate no contradiction, run color_argument2(v_1, done_list + v_0, current_distance +1, max_distance, max_depth). If False (no contradiction) then append v_1 to list1. 

        list2: Check all vertices in the 'unknown' list of v_0. If v_1 in the 'unknown' list of v_0 has no neighbors colored c_0, make a copy of the graph G', and within G' color v_1 with c_0 and add an edge from v_0 to v_1. Run the methods that extend colorings appropriately through the graph, and check for forks or triangles. If these methods indicate no contradiction, run color_argument2(v_1, done_list + v_0, current_distance +1, max_distance, max_depth). If False (no contradiction) then append v_1 to list2.

        list3: If v_0 is not max degree, then make a copy G' of the graph. add a new vertex v_1 to G' and color it with c_0 and add an edge between v_1 and v_0. Run the methods that extend colorings appropriately through the graph, and check for forks or triangles. If these methods indicate no contradiction, run color_argument2(v_1, done_list + v_0, current_distance +1, max_distance, max_depth). If False (no contradiction) then append v_1 to list3.

        list4: Check all vertices in the 'unknown' list of v_0 again. If v_1 in the 'unknown' list of v_0 is colored with c_0, make a copy of the graph G', and within G' add an edge from v_0 to v_1. Run the methods that extend colorings appropriately through the graph, and check for forks or triangles. If these methods indicate no contradiction, run color_argument2(v_1, done_list + v_0, current_distance +1, max_distance, max_depth). If False (no contradiction) then append v_1 to list4.

    Check how many elements are in these lists:
        If all the lists are empty (for any color) return True
        If exactly one list is nonempty and has exactly 1 element, then add the according structure (NOTE: this is done on the current graph- in the case of recursion, this may be done on a copy, not the original):
            If list1 is nonempty, color the vertex of that index with c_0
            If list2 is nonempty, add an edge from v_0 to the vertex of that index and color that vertex with c_0
            If list3 is nonempty, add a new vertex to the graph adjacent to v_0 of color c_0
            If list4 is nonempty, add an edge from v_0 to the vertex of that index

return False (i.e., the default is to return False unless something prompts us to return 'True')