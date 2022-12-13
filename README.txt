This  repository contains the code to verify parts of the proof that triangle-free and fork-free graphs are 3-colorable. The most recent test build of the code can be found within 'testing6.1.ipynb'. This README file may be expanded with further documentation and explanation of the methods at a later date. For now, a few quick things to know:

vert: the class of vertices of graphs each vertex has .edges, .non_edges, and .unknown containing the indices of other vertices it has edges, known non-edges, or not known whether there is an edge or not an edge respectively. Vertices can also have .max_degree or .color properties, though they are by default initialized to None. 

graph: the class of graphs, containing the methods to do the graph algorithms. contains a vertex list containing all the vertices of the graph, as well as a list of indices and a hidden property .depth for computational purposes. 

Most relevant methods for graph:

.clean(): will update edge/non_edge/unknown lists so that if a vertex index is added to the edge list of one vertex, the corresponding lists of other vertices will be updated appropriately. It is best to run this after adding edges to the graph.

.contains_subgraphs(subgraph_list): searches the graph for any instances of any subgraphs from the subgraph list within the graph. Will return True if a subgraph from the list is present, and False otherwise

.update_graph(subgraph_list): will try to do a shallow search for any forced structure implied by the subgraphs in subgraph list being forbidden in the graph, adding any forced edges/non_edges along the way. Will return True if in the process of updating the graph, a contradiction is forced, otherwise will return False and the graph will be updated in place.

.update_graph2(subgraph_list, max_depth): Heavier and slower version of the update_graph function that does a deeper search for forced structure, but can solve more cases. Return type is the same. max_depth indicates the maximum recursion to which a type of logic is allowed to run for this method, though 1 is generally.

.color_argument2_iter(subgraph_list, start_vert, max_distance): An expansion upon update_graph2(), this method also uses structure that is forced by partial coloring in the graph to try to build out the graph and seek a contradiction. subgraph_list is the set of forbidden subgraphs, start_vert is the single vertex that we assume is forced to be unable to be colored within the graph, and max_distance is the bound on how far the algorithm will recurse away from start_vert within the graph to try to obtain a contradiction. It is advised that runtime grows very very fast as a function of max_distance. 