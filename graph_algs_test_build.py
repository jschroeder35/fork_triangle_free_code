import collections
import copy
from itertools import permutations
from itertools import combinations

class trienode:
    def __init__(self):
        # whether this can be the end of a subgraph
        self.is_end = False
        # a dictionary of children nodes
        self.children = {}

class trie:
    def __init__(self):
        ###start corresponds to empty subgraph
        self.root=trienode()

    def add_subgraph(self, graph):
        temp=self.root
        for i in range(len(graph.vertex_list)):
            vert=graph.vertex_list[i]
            ###Keep only edges to vertices that precede it in the graph, and convert to a key for the children dictionary
            temp_edges=[j for j in vert.edges if j<i]
            temp_non_edges=[j for j in vert.non_edges if j<i]
            key=(frozenset(temp_edges), frozenset(temp_non_edges))
            ###Add key to children dictionary if not already present
            if key not in temp.children:
                temp.children[key]=trienode()
            ###Move down the trie
            temp=temp.children[key]
            if i==len(graph.vertex_list)-1:
                temp.is_end=True

    def initialize(self, subgraph_list):
        ###Adds all the subgraphs to the list.
        ###I envision maybe one day adding code to search for the optimal ordering of vertices within the subgraphs for a maximally efficient trie, but that is not yet implemented.Checking for such an ordering by hand may be faster for this problem.
        for graph in subgraph_list:
            self.add_subgraph(graph)

    def show(self):
        ###This method displays the nodes of the trie.
        ###Nodes of equal depth are printed on the same line.
        queue=collections.deque([(self.root, 0)])
        while queue:
            temp, depth=queue.popleft()
            for key in temp.children:
                print(key, depth, end=', ')
                queue.append((temp.children[key], depth+1))
                if queue[0][1] > depth:
                    print('')

trielist=[]
etrielist=[]
graph_hash={}
color_hash={}

class etrie:
    ###This class is for tries containing keys with only edges (not non-edges) for the purpose of generating non-edge lists
    def __init__(self):
        ###start corresponds to empty subgraph
        self.root=trienode()
        self.non_edge_list=[]

    def add_subgraph(self, graph):
        temp=self.root
        for i in range(len(graph.vertex_list)):
            vert=graph.vertex_list[i]
            ###Keep only edges to vertices that precede it in the graph, and convert to a key for the children dictionary
            temp_edges=[j for j in vert.edges if j<i]
            temp_non_edges=[j for j in vert.non_edges if j<i]
            key=frozenset(temp_edges)
            ###Add key to children dictionary if not already present
            if key not in temp.children:
                temp.children[key]=trienode()
            ###Move down the trie
            temp=temp.children[key]
            temp.non_edge_list=temp_non_edges
            if i==len(graph.vertex_list)-1:
                temp.is_end=True

    def initialize(self, subgraph_list):
        ###Adds all the subgraphs to the list.
        ###I envision maybe one day adding code to search for the optimal ordering of vertices within the subgraphs for a maximally efficient trie, but that is not yet implemented.Checking for such an ordering by hand may be faster for this problem.
        for graph in subgraph_list:
            self.add_subgraph(graph)

    def show(self):
        ###This method displays the nodes of the trie.
        ###Nodes of equal depth are printed on the same line.
        queue=collections.deque([(self.root, 0)])
        while queue:
            temp, depth=queue.popleft()
            for key in temp.children:
                print(key, temp.children[key].non_edge_list, depth, end=', ')
                queue.append((temp.children[key], depth+1))
                if queue[0][1] > depth:
                    print('')


class vert: #####The class of vertices #####
    def __init__(self, edges, non_edges, index):
        self.edges = edges #####List of indices of vertices adjacent to this vertex #####
        self.non_edges = non_edges ##### List of indices of vertices with a non-edge to this vertex #####
        self.index= index #####label unique to each vertex #####
        self.max_degree=None #####max degree, can be assigned to vertices later#####
        self.color=None ##### Color assigned to this vertex, can be assigned later #####
        self.color_start=False ##### This is set to 'True' if a partial coloring of G originates around this vertex being unable to be colored #####
        self.color_list=None ###Used only within algorithms

class graph:
    def __init__(self, vertex_list):
        self.vertex_list=vertex_list ##### List of vertex objects #####
        self.index_list=[vertex.index for vertex in self.vertex_list] ##### List of indices of the vertices #####
        self.depth1=0 ##### Used only within algorithms
        self.depth2=0 ##### Used only within algorithms

    def show(self): ##### Method to display the graph #####
        for vertex in self.vertex_list:
            print('Vertex ' + str(vertex.index) + ':')
            print('    Edges:     ' + str(vertex.edges))
            print('    Non-edges: ' + str(vertex.non_edges))
            print('    Unknown:   ' + str(vertex.unknown))
            print('    Color:   ' + str(vertex.color))

    def contains_subgraphs(self, subgraph_list):
        ###This method wants to check if any graph from subgraph_list appears in self.
        ###If so, it returns True, otherwise it returns false.
        
        ###check if graph is in hash map first
        subgraphs_key=[]
        for subgraph in subgraph_list:
            subgraphs_key.append(tuple([(frozenset(vertex.edges), frozenset(vertex.non_edges), vertex.color) for vertex in subgraph.vertex_list]))
        graph_key = (tuple([(frozenset(vertex.edges), frozenset(vertex.non_edges), vertex.color) for vertex in self.vertex_list]), tuple(subgraphs_key))
        if graph_key in graph_hash:
            return graph_hash[graph_key]
        
        ### the method generates a trie from the subgraph list, then performs a dfs starting from each vertex using the trie
        ttrie=None
        for graph_list, temp in trielist:
            if subgraph_list==graph_list:
                ttrie=temp
        if not ttrie:
            ttrie=trie()
            ttrie.initialize(subgraph_list)
            trielist.append((subgraph_list, ttrie))
        #ttrie.show()
        ###Make a queue to store all remaining positions to search from while doing the dfs
        queue=collections.deque()
        for vertex in self.vertex_list:
            queue.append((vertex, ttrie.root, []))

        while queue:
            vertex, trienode, past_vertex_list = queue.pop()
            tpast_vertex_list=copy.copy(past_vertex_list)
            ###need to convert vertex indices to match the indices of the subgraph, and filter down to only those corresponding to vertices already added to our partial list
            temp_edges=[i for i, vert in enumerate(past_vertex_list) if vert in vertex.edges]
            temp_non_edges=[i for i, vert in enumerate(past_vertex_list) if vert in vertex.non_edges]
            key=(frozenset(temp_edges), frozenset(temp_non_edges))
            #print(vertex.index, past_vertex_list)
            ###now we can check if the vertex is a valid addition to build out one of our subgraphs.
            if key in trienode.children:
                trienode=trienode.children[key]
                tpast_vertex_list.append(vertex.index)
                ###check if a subgraph has been completed. If so, return True
                if trienode.is_end==True:
                    graph_hash[graph_key]=True
                    return True
                ###Continue the search on vertices not yet added to the subgraph (if it's not completed)
                for vertex2 in self.vertex_list:
                    if (vertex2.index not in tpast_vertex_list):
                        queue.append((vertex2, trienode, tpast_vertex_list))
        graph_hash[graph_key]=False
        return False

    def get_edgelists(self, subgraph_list):
        ###This method seeds to return a list of edgelists that are forced to prevent any of the subgraphs in the list from occurring.
        ###returns a list of frozensets of of edge pairs (that are frozensets)
        ### the method generates a trie from the subgraph list, then performs a dfs starting from each vertex using the trie, similar to the previous method, but the dfs also tracks added non-edges
        ttrie=None
        for graph_list, temp in etrielist:
            if subgraph_list==graph_list:
                ttrie=temp
        if not ttrie:
            ttrie=etrie()
            ttrie.initialize(subgraph_list)
            etrielist.append((subgraph_list, ttrie))
        #ttrie.show()
        out=[]
        ###queue stores positions left to search in the dfs
        queue=collections.deque()
        for vertex in self.vertex_list:
            queue.append((vertex, ttrie.root, [], []))

        while queue:
            vertex, trienode, past_vertex_list, edgelist = queue.pop()
            ###need to convert vertex indices to match the indices of the subgraph, and filter down to only those corresponding to vertices already added to our partial list
            temp_edges=[i for i, vert in enumerate(past_vertex_list) if vert in vertex.edges]
            temp_non_edges=[i for i, vert in enumerate(past_vertex_list) if vert in vertex.non_edges]
            key=frozenset(temp_edges)
            #print(vertex.index, past_vertex_list, key, trienode.children)
            ###now we can check if the vertex is a valid addition to build out one of our subgraphs.
            if key in trienode.children:
                trienode=trienode.children[key]
                ###Need copies of these lists to pass to the queue, otherwise things break because of the way python passes lists by reference
                tpast_vertex_list=copy.copy(past_vertex_list)
                tedgelist=copy.copy(edgelist)
                for j in set(trienode.non_edge_list)-set(temp_non_edges):
                    tedgelist.append(frozenset([past_vertex_list[j], vertex.index]))
                tpast_vertex_list.append(vertex.index)
                ###check if a subgraph has been completed. If so, append the temporary edgelist
                if trienode.is_end==True:
                    out.append(tedgelist)
                    #print(past_vertex_list)
                for vertex2 in self.vertex_list:
                    if vertex2.index not in tpast_vertex_list:
                        queue.append((vertex2, trienode, tpast_vertex_list, tedgelist))

        #####Now we want to condense the list as much as possible #####
        temp_list=[]
        for edge_list in out:
            edge_set = frozenset(edge_list)
            temp_list.append(edge_set)
        temp_list=list(set(temp_list))
        out=copy.copy(temp_list)
        for i in range(len(temp_list)):
            for j in range(i+1, len(temp_list)):
                if temp_list[i].issubset(temp_list[j]):
                    if temp_list[j] in out:
                        out.remove(temp_list[j])
        return out

    def get_subgraph_list(self, subgraph):
        ###this method generates a list of all distinct copies of subgraph in the graph, up to reordering
        
        ### the method generates a trie from the subgraph list, then performs a dfs starting from each vertex using the trie
        ttrie=None
        for graph_list, temp in trielist:
            if [subgraph]==graph_list:
                ttrie=temp
        if not ttrie:
            ttrie=trie()
            ttrie.initialize([subgraph])
            trielist.append(([subgraph], ttrie))
        #ttrie.show()
        ###Make a queue to store all remaining positions to search from while doing the dfs
        queue=collections.deque()
        out=set()
        for vertex in self.vertex_list:
            queue.append((vertex, ttrie.root, []))

        while queue:
            vertex, trienode, past_vertex_list = queue.pop()
            tpast_vertex_list=copy.copy(past_vertex_list)
            ###need to convert vertex indices to match the indices of the subgraph, and filter down to only those corresponding to vertices already added to our partial list
            temp_edges=[i for i, vert in enumerate(past_vertex_list) if vert in vertex.edges]
            temp_non_edges=[i for i, vert in enumerate(past_vertex_list) if vert in vertex.non_edges]
            key=(frozenset(temp_edges), frozenset(temp_non_edges))
            #print(vertex.index, past_vertex_list)
            ###now we can check if the vertex is a valid addition to build out one of our subgraphs.
            if key in trienode.children:
                trienode=trienode.children[key]
                tpast_vertex_list.append(vertex.index)
                ###check if a subgraph has been completed. If so, return True
                if trienode.is_end==True:
                    out.add(frozenset(tpast_vertex_list))
                ###Continue the search on vertices not yet added to the subgraph (if it's not completed)
                for vertex2 in self.vertex_list:
                    if (vertex2.index not in tpast_vertex_list):
                        queue.append((vertex2, trienode, tpast_vertex_list))
        return out
    
    def clean(self):
        #####This method fixes the edge/non-edge lists of the graph. #####
        #####If v1 has v2 in its edge list, but v2 doesn't have v1 in its edge list, then clean will append v1 to v2's edge list. #####
        #####After all lists have been updated, it will make an 'unknown' list for each vertex, of neither edges nor non-edges. #####
        #####It is recommended to call this method after making changes to a graph. ##### '
        self.index_list=[vertex.index for vertex in self.vertex_list]
        i=0
        while True:
            i+=1
            if i==1000:
                print('Clean got stuck in a loop. There was an invalid graph update. You need to interrupt the kernel and do some debugging.')
            updated=False
            for vertex in self.vertex_list:
                #####Making sure edge lists match up #####
                for vertex2_index in vertex.edges:
                    vertex2=self.vertex_list[vertex2_index]
                    if vertex.index not in vertex2.edges:
                        vertex2.edges.append(vertex.index)
                        updated=True
                #####Making sure non-edge lists match up #####
                for vertex2_index in vertex.non_edges:
                    vertex2=self.vertex_list[vertex2_index]
                    if vertex.index not in vertex2.non_edges:
                        vertex2.non_edges.append(vertex.index)
                        updated=True
                #####If any vertex becomes max degree, we prohibit any further edges for that vertex#####
                if len(vertex.edges)==vertex.max_degree:
                    old=copy.copy(vertex.non_edges)
                    vertex.non_edges = list(set(self.index_list) - set(vertex.edges))
                    if set(vertex.non_edges)!=set(old):
                        updated=True
            #####Now we check if the graph was updated#####
            if not updated:
                break
        #####Once all edge/non-edge lists are updated, we set the unknown lists for each vertex #####
        for vertex in self.vertex_list:
            vertex.unknown = list(set(self.index_list).difference(vertex.edges + vertex.non_edges + [vertex.index]))



    def check_for_subgraphs(self, subgraph_list, depth_cutoff=1):
        #####This method scans the graph for graphs in the sugraph list, as well as seeing if appending certain edges/nonedges is forced by forks or triangles. #####
        #####Will test all combinations of n edges/non-edges, where n is the depth_cutoff#####
        #####In practice, little is gained by a depth cutoff of more than one#####
        #####Method returns True if a fork/triangle is forced, False otherwise #####
        if self.contains_subgraphs(subgraph_list):
            return True

        if self.depth1 < depth_cutoff:
            for vertex in self.vertex_list:
                if len(vertex.unknown) > 0:
                    for edge in vertex.unknown:
                        if edge > vertex.index: #####Prevents running each case twice, to cut down the runtime #####
                            #####Make a temporary graph with the appended edge#####
                            temp_graph = copy.deepcopy(self)
                            temp_graph.vertex_list[vertex.index].edges.append(edge)
                            temp_graph.clean()
                            temp_graph.depth1=self.depth1+1 #####Need to increase the depth, so it doesn't recurse forever

                            #####Make another temporary graph with the appended non-edge #####
                            temp_graph2 = copy.deepcopy(self)
                            temp_graph2.vertex_list[vertex.index].non_edges.append(edge)
                            temp_graph2.clean()
                            temp_graph2.depth1=self.depth1+1 ####Need to increase the depth, so it doesn't recurse forever

                            ##### Check both cases. If there is a contradiction in either case, return True #####
                            bool1 = temp_graph.check_for_subgraphs(subgraph_list, depth_cutoff)
                            bool2 = temp_graph2.check_for_subgraphs(subgraph_list, depth_cutoff)
                            if (bool1 and bool2):
                                return True
                            ##### In this case, only and edge leads to a fork/triangle, so we append a non-edge #####
                            elif bool1:
                                vertex.non_edges.append(edge)
                                self.clean()
                            ##### In this case, an edge is added #####
                            elif bool2:
                                vertex.edges.append(edge)
                                self.clean()
        return False

    def check_for_subgraphs_iter(self, subgraph_list, depth_cutoff=1):
        #####This method runs the prior check_for_sugraphs method iteratively until the graph is no longer updated #####
        #####Like most methods, this returns True if there is a contradiction and False otherwise #####
        for i in sorted(list({1, depth_cutoff})): #####For runtime reasons, I set the script to run at depth 1 first, then at the specified depth. Now, though, I usually just run the code with depth 1, so not sure this should be used.
            while True:
                old_graph = copy.deepcopy(self)
                bool1=self.check_for_subgraphs(subgraph_list, i)
                if bool1:
                    return True
                updated=False
                for j, vertex in enumerate(old_graph.vertex_list):
                    if (vertex.edges != self.vertex_list[j].edges) or (vertex.non_edges != self.vertex_list[j].non_edges):
                        updated = True
                if not updated:
                    break
        return False

    def fill_in_graph(self, subgraph_list, depth_cutoff=1):
        #####This method runs check_for_subgraphs_iter, as well as appending new neighbors where there is neighborhood containment. #####
        #####I don't really use it currently, though, so I'm going to to go light on annotating it #####
        while True:
            old_graph = copy.deepcopy(self)
            ##### Run check_for_forks_triangles_iter and check for updates #####
            bool1 = self.check_for_subgraphs_iter(subgraph_list, depth_cutoff)
            if bool1:
                #print('Contradiction!')
                return True
            updated=False
            for j, vertex in enumerate(old_graph.vertex_list):
                if (vertex.edges != self.vertex_list[j].edges) or (vertex.non_edges != self.vertex_list[j].non_edges):
                    updated = True

            ##### Append vertices if there is neighborhood containment ######
            for vertex1 in self.vertex_list:
                for vertex2 in self.vertex_list:
                    if ((vertex1.index != vertex2.index) and (len(vertex1.unknown) ==0) and (set(vertex1.edges).issubset(set(vertex2.edges)))):
                        if (vertex1.max_degree == None) or ((vertex1.max_degree != None) and (len(vertex1.edges) < vertex1.max_degree)):
                            new_vertex=vert([vertex1.index], [vertex2.index], max(self.index_list) + 1)
                            self.vertex_list.append(new_vertex)
                            self.clean()
                            #print('New vertex added!')
                            updated=True
                        elif ((vertex1.max_degree != None) and (len(vertex1.edges) >= vertex1.max_degree)):
                            #print('Contradiction!')
                            return True
            if not updated:
                break
        return False

    def clean_colors(self): #####Adds non-edges according to the coloring of the graph #####
        for color in [0, 1, 2]:
            index_list=[]
            last_vert = None
            for vertex in self.vertex_list:
                if vertex.color==color:
                    vertex.non_edges = list(set(vertex.non_edges + index_list))
                    index_list.append(vertex.index)
        self.clean()

    def color_vert(self, start_vert):
        #####This method attempts to color a single vertex #####
        color_list = []
        uncolored_neighbors = []
        start_vertex=self.vertex_list[start_vert]

        #####Make the list of colors of neighbors, and indices of uncolored neighbors ########
        for neighbor in start_vertex.edges:
            vertex = self.vertex_list[neighbor]
            if vertex.color != None:
                color_list.append(vertex.color)
            else:
                uncolored_neighbors.append(vertex.index)
        color_list = list(set(color_list))

        ######Color start_vertex, if uncolored #########
        if (len(color_list) == 3) and (start_vertex.color_start==False):
            #print('Contradiction!')
            return True
        if (len(color_list) == 2) and (start_vertex.color_start==False):
            if start_vertex.color == None:
                start_vertex.color = list(set([0, 1, 2]) - set(color_list))[0]
            elif (start_vertex.color != list(set([0, 1, 2]) - set(color_list))[0]):
                return True

        #######If it is the first vertex being colored and it is degree 3, color all its neighbors with different colors #######
        if start_vertex.color_start==True:
            if (start_vertex.max_degree ==3):
                i=0
                for neighbor in start_vertex.edges:
                    self.vertex_list[neighbor].color=i
                    i=i+1

        ########Fix up the graph (non-edges) after adding new colors #########
        self.clean_colors()
        return False

    def fill_in_colors(self, start_vert=None):
        #####This method fills in a partial coloring of the graph, starting from start_vert, and assuming that start_vert cannot be colored.#####
        #####can be run with start_vert=None to just update colors for the graph #######
        #####Color neighborhood of starting vertex ##########

        if start_vert !=None:
            vertex=self.vertex_list[start_vert]
            vertex.color_start = True
            if self.color_vert(start_vert)==True:
                return True

        ######Iteratively add colors for the rest of the graph ########
        while True:
            old_graph = copy.deepcopy(self)
            #####Color the vertex while also checking for a contradiction #####
            for vertex in self.vertex_list:
                if self.color_vert(vertex.index)==True:
                    return True
            #####Check for updates
            updated =False
            for vertex1 in self.vertex_list:
                if vertex1.color != old_graph.vertex_list[vertex1.index].color:
                    updated=True
            if not updated:
                break
        return False

    def fill_in_colors2(self, subgraph_list, start_vert=None):
        #####Heavier, more powerful version of fill_in_colors imbued with the logic of forbidden-subgraph case-checking #####
        while True:
            #####First fill in colors normally, then get the edge lists implied by the subgraph_list
            if self.fill_in_colors(start_vert):
                return True
            old_graph = copy.deepcopy(self)
            edge_lists = self.get_edgelists(subgraph_list)
            for edge_list in edge_lists:
                #####reset the color_lists#####
                for vertex in self.vertex_list:
                    if vertex.color==None:
                        vertex.color_list=[]
                #####Iterate over the edges and check their possibilities
                for set1 in edge_list:
                    edge = list(set1)
                    temp_graph=copy.deepcopy(self)
                    temp_graph.vertex_list[edge[0]].edges.append(edge[1])
                    temp_graph.clean()
                    if temp_graph.update_graph(subgraph_list):
                        continue
                    #####Append any newly deduced colors to a list
                    for vertex in self.vertex_list:
                        if (vertex.color == None) and (temp_graph.vertex_list[vertex.index].color != None):
                            vertex.color_list.append(temp_graph.vertex_list[vertex.index].color)
                        elif (vertex.color == None) and (temp_graph.vertex_list[vertex.index].color == None):
                            vertex.color_list.append(-1)
                for vertex in self.vertex_list:
                    if (vertex.color==None) and (-1 not in vertex.color_list) and (len(set(vertex.color_list))==1):
                        vertex.color = vertex.color_list[0]
            ###Check for updates
            updated =False
            for vertex in self.vertex_list:
                if vertex.color != old_graph.vertex_list[vertex.index].color:
                    updated=True
            if not updated:
                break

    def update_graph(self, subgraph_list, start_vert=None, with_colors=True):
        ###Lightweight graph update function that fills in colors and checks for subgraphs without doing the heavier methods that incur longer runtime
        while True:
            old_graph = copy.deepcopy(self)
            if with_colors:
                if self.fill_in_colors(start_vert) or self.check_for_subgraphs_iter(subgraph_list):
                    return True
            else:
                if self.check_for_subgraphs_iter(subgraph_list):
                    return True
            updated =False
            for vertex in self.vertex_list:
                if ((vertex.color != old_graph.vertex_list[vertex.index].color) or
                    (set(vertex.edges) != set(old_graph.vertex_list[vertex.index].edges)) or
                    (set(vertex.non_edges) != set(old_graph.vertex_list[vertex.index].non_edges))
                   ):
                    updated=True
            if not updated:
                break
        return False

    def update_graph2(self, subgraph_list, max_depth, start_vert=None, current_depth=0, verbose=False, with_colors=True):
        #####This graph should append edges from the lists produced by forks, then get all structure forced by coloring/forks/triangles, then recurse up to max_depth and append all structure that holds in all cases#####
        while True:
            if self.update_graph(subgraph_list, start_vert, with_colors):
                return True
            old_graph = copy.deepcopy(self)
            updated=False
            edge_lists = self.get_edgelists(subgraph_list)
            if verbose:
                print('Edge lists to check:')
                for edge_list in edge_lists:
                    print(edge_list)
            for i, edge_list in enumerate(edge_lists):
                if verbose:
                    print('Completion={:8.2f}%'.format(100*i/len(edge_lists)))
                #print(edge_list)
                #####reset the graph_list#####
                graph_list=[]
                #####Iterate over the edges and check their possibilities
                for set1 in edge_list:
                    edge = list(set1)
                    temp_graph=copy.deepcopy(self)
                    temp_graph.vertex_list[edge[0]].edges.append(edge[1])
                    temp_graph.clean()
                    if temp_graph.update_graph(subgraph_list, start_vert, with_colors):
                        continue
                    #####Recurse if less than max_depth
                    if current_depth < max_depth -1:
                        if temp_graph.update_graph2(subgraph_list, max_depth, start_vert, current_depth+1):
                            continue
                    graph_list.append(temp_graph)
                #####This indicates there was a contradiction in all possible cases #####
                if len(graph_list)==0:
                    return True
                else:
                    for vertex in self.vertex_list:
                        #####Update colors #####
                        if vertex.color==None:
                            vertex.color_list=[]
                            for temp_graph in graph_list:
                                if temp_graph.vertex_list[vertex.index].color ==None:
                                    vertex.color_list.append(-1)
                                else:
                                    vertex.color_list.append(temp_graph.vertex_list[vertex.index].color)
                            if (-1 not in vertex.color_list) and (len(set(vertex.color_list))==1):
                                vertex.color=vertex.color_list[0]
                                updated=True
                        #####Update Edges #####
                        for index1 in vertex.unknown:
                            edge_list=[]
                            for temp_graph in graph_list:
                                if index1 in temp_graph.vertex_list[vertex.index].edges:
                                    edge_list.append(1)
                                elif index1 in temp_graph.vertex_list[vertex.index].non_edges:
                                    edge_list.append(-1)
                                else:
                                    edge_list.append(0)
                            if (len(set(edge_list))==1) and (0 not in edge_list):
                                if edge_list[0] ==1:
                                    vertex.edges.append(index1)
                                    updated=True
                                else:
                                    vertex.non_edges.append(index1)
                                    updated=True
                self.clean()
            if updated and verbose:
                print('Updated, running again...')
            if not updated:
                break
        return False

    def update_graph3(self, subgraph_list, start_vert=None, verbose=False, with_colors=True):
        ###Heaviest update graph function to date, checking every unknown edge to see if adding it as an edge/non-edge will obtain a contradiction with update_graph2()
        ###Under construction
        stored_index=0 ###Stores the index of the last visited vertex, so that as the graph is updated, vertice can be visited in a cyclical pattern.

        while True:
            if verbose:
                print('Setting up and and initial greedy updates...')
            ###Update graph with update_graph2 first, and again after adding any edge/non_edge
            if self.update_graph2(subgraph_list, 1, start_vert, 0, verbose=False, with_colors=with_colors):
                if verbose:
                    print('Contradiction obtained (if this is not the first line of output, there is an algorithmic error)')
                    self.show()
                return True
            if verbose:
                self.show()

            ###count number of unknown edge/non_edges in the graph, for progress printing purposes.
            if verbose:
                unknown_count=0
                for vertex in self.vertex_list:
                    for v2_index in vertex.unknown:
                        if v2_index>vertex.index:
                            unknown_count+=1
                print('There are {} unknown edges in the graph to check'.format(unknown_count))

            updated=False
            j=0
            for i in range(len(self.vertex_list)):
                vertex=self.vertex_list[(stored_index+i)%len(self.vertex_list)] ###picks up from where the last iteration left off, roughly
                if not updated and len(vertex.unknown) > 0:
                    for edge in vertex.unknown:
                        if not updated and edge > vertex.index: #####Prevents running each case twice, to cut down the runtime #####
                            if verbose:
                                print('Checking {}<->{}. Completion: {:8.2f}%'.format(vertex.index, edge, 100*j/unknown_count))
                            j+=1
                            #####Make a temporary graph with the appended edge#####
                            temp_graph = copy.deepcopy(self)
                            temp_graph.vertex_list[vertex.index].edges.append(edge)
                            temp_graph.clean()

                            #####Make another temporary graph with the appended non-edge #####
                            temp_graph2 = copy.deepcopy(self)
                            temp_graph2.vertex_list[vertex.index].non_edges.append(edge)
                            temp_graph2.clean()

                            ##### Check both cases. If there is a contradiction in either case, return True #####
                            bool1 = temp_graph.update_graph2(subgraph_list, 1, start_vert, 0, verbose=False, with_colors=with_colors)
                            bool2 = temp_graph2.update_graph2(subgraph_list, 1, start_vert, 0, verbose=False, with_colors=with_colors)
                            if (bool1 and bool2):
                                if verbose:
                                    print('Contradiction obtained when testing unknown edge/non-edge between {} and {}!'.format(vertex.index, edge))
                                return True
                            ##### In this case, only and edge leads to a fork/triangle, so we append a non-edge #####
                            elif bool1:
                                if verbose:
                                    print('Non-edge added from {} to {}; starting over.'.format(vertex.index, edge))
                                vertex.non_edges.append(edge)
                                self.clean()
                                stored_index=vertex.index
                                updated=True
                            ##### In this case, an edge is added #####
                            elif bool2:
                                if verbose:
                                    print('Edge added from {} to {}; starting over.'.format(vertex.index, edge))
                                vertex.edges.append(edge)
                                self.clean()
                                stored_index=vertex.index
                                updated=True
            if not updated:
                break
        return False

    def update_graph4(self, subgraph_list, start_vert=None, verbose=False, with_colors=True):
        ###Heaviest update graph function to date, checking every unknown edge to see if adding it as an edge/non-edge will obtain a contradiction with update_graph2()
        ###Under construction
        while True:
            updated=False
            if self.update_graph3(subgraph_list, start_vert, verbose, with_colors):
                return True
            for vertex1 in self.vertex_list:
                for vertex2 in self.vertex_list:
                    if ((vertex1.index != vertex2.index) and (len(vertex1.unknown) ==0) and (set(vertex1.edges).issubset(set(vertex2.edges)))):
                        if (vertex1.max_degree == None) or ((vertex1.max_degree != None) and (len(vertex1.edges) < vertex1.max_degree)):
                            new_vertex=vert([vertex1.index], [vertex2.index], max(self.index_list) + 1)
                            self.vertex_list.append(new_vertex)
                            self.clean()
                            if verbose:
                                print('Vertex {} added! '.format(new_vertex.index))
                            updated=True
                        elif ((vertex1.max_degree != None) and (len(vertex1.edges) >= vertex1.max_degree)):
                            return True
            if not updated:
                break
        return False
    
    def update_graph5(self, subgraph_list, add_vertices=True, quick_update=False, max_distance=10, verbose=False):
        ###update the graph and seek smart color argument updates in an automated fashion. 
        
        ###First, we need to clear away any partial coloring in the graph
        for vertex in self.vertex_list:
            vertex.color=None
            vertex.color_start=False
            
        while True:
            ###run either update_graph3 or update_graph4 based on whether adding vertices is allowed
            if add_vertices:
                if self.update_graph4(subgraph_list, start_vert=None, verbose=verbose, with_colors=False):
                    return True
            else:
                if self.update_graph3(subgraph_list, start_vert=None, verbose=verbose, with_colors=False):
                    return True
            
            ###now, search for pairs of vertices with neighborhood containment
            old_graph = copy.deepcopy(self)
            
            for v1 in self.vertex_list:
                for v2 in self.vertex_list:
                    for v3 in self.vertex_list:
                        ###check that vertices are distinct, and that N(v1) < N(v2) \cup N(v3)
                        if not (v1 !=v2 and v1 !=v3 and v2.index<v3.index and set(v1.edges + v1.unknown).issubset(set(v2.edges + v3.edges))):
                            continue
                        ###Add a new neighbor adjacent to v1 not adjacent to v2 or v3 and check if that yields a contradiction
                        temp_graph=copy.deepcopy(self)
                        new_vert=vert([v1.index], [v2.index, v3.index], len(temp_graph.vertex_list))
                        temp_graph.vertex_list.append(new_vert)
                        temp_graph.clean()
                        if verbose:
                            print('Checking if N({}) is contained in N({}) and N({})'.format(v1.index, v2.index, v3.index))
                        ###If the update returns a contradiction, proceed, otherwise break
                        if quick_update:
                            if not temp_graph.update_graph2(subgraph_list, 1, start_vert=None, verbose=False, with_colors=False):
                                continue
                        else:
                            if not temp_graph.update_graph3(subgraph_list, start_vert=None, verbose=False, with_colors=False):
                                continue

                        ###Now we run a coloring argument
                        if verbose:
                            print('N({}) is contained in N({}) and N({}). Attempting coloring...'.format(v1.index, v2.index, v3.index))
                        v2.color=0
                        v3.color=1
                        v1.color_start=True
                        if self.fill_in_colors():
                            if verbose:
                                print('Contradiction obtained filling in colors.')
                            return True
                        if self.color_argument_iter(subgraph_list, v1.index, max_distance, verbose=verbose, add_vertices=add_vertices):
                            return True
                        
                        ###reset the colors after
                        for vertex in self.vertex_list:
                            vertex.color=None
                            vertex.color_start=False
            
            ###check for updates
            updated=False
            if len(self.vertex_list) != len(old_graph.vertex_list):
                updated=True
            else:
                for vertex in self.vertex_list:
                    if ((set(vertex.edges) != set(old_graph.vertex_list[vertex.index].edges)) or
                        (set(vertex.non_edges) != set(old_graph.vertex_list[vertex.index].non_edges))
                       ):
                        updated=True
            if not updated:
                break
                
        return False
    
    def color_argument(self, subgraph_list, start_vert, done_list, current_distance, max_distance, max_depth, verbose=True, add_vertices=True):
        #####This method recursively builds out new structure to the graph #####
        ##### done_list is the set of vertices already visited by color_argument2 within the recursion ######
        ##### current_distance is the recursion depth/distance in the graph from the starting vertex
        ##### max_distance is the max distance allowed within the graph
        ##### max_depth is the max recursion depth to which the graph is allowed to be modified
        start_vertex=self.vertex_list[start_vert]

        #####color_list is colors to be forced #####
        color_list = list(set([0, 1, 2]) - set([start_vertex.color]))

        #####It leads to algorithmic errors if you don't check that you haven't looped back to the start######
        if current_distance > 1:
            for edge_index in start_vertex.edges:
                if (self.vertex_list[edge_index].color_start ==True):
                    return False

        #####Check for each color that it exists among neighbors#########
        done_list2=[pair[0] for pair in done_list]
        if (current_distance < max_distance) and (self.depth2 < max_depth):
            for temp_color in color_list:
                #####Check for already colored neighbors#####
                list0 = []
                for edge_index in start_vertex.edges:
                    vertex1=self.vertex_list[edge_index]
                    if (vertex1.color == temp_color) and (vertex1.index not in done_list2):
                        list0.append(vertex1.index)

                #####Check for uncolored neighbors that could be assigned that color#####
                list1=[]
                for edge_index in start_vertex.edges:
                    vertex1=self.vertex_list[edge_index]
                    if (vertex1.color == None) and (vertex1.index not in done_list2):
                        #####Check that vertex1 can be colored with temp_color #####
                        neighbor_colors = []
                        for edge_index2 in vertex1.edges:
                            neighbor_colors.append(self.vertex_list[edge_index2].color)
                        if temp_color not in neighbor_colors:
                            #####Copy graph and color vertex #####
                            temp_graph = copy.deepcopy(self)
                            temp_graph.depth2= self.depth2+1
                            temp_graph.vertex_list[vertex1.index].color=temp_color

                            #####Assuming no contradiction, recurse on vertex1#####
                            if not temp_graph.update_graph2(subgraph_list, 1):
                                list1.append(vertex1.index)

                #####Check for potential neighbors that could be that color#####
                list2=[]
                for edge_index in start_vertex.unknown:
                    vertex1=self.vertex_list[edge_index]
                    if (vertex1.color == None) and (vertex1.index not in done_list2):
                        #####Check that it can be assigned temp_color #####
                        neighbor_colors = []
                        for edge_index2 in vertex1.edges:
                            neighbor_colors.append(self.vertex_list[edge_index2].color)
                        if temp_color not in neighbor_colors:
                            #####Copy graph and add color and the new edge #####
                            temp_graph = copy.deepcopy(self)
                            temp_graph.depth2= self.depth2+1
                            vertex3=temp_graph.vertex_list[vertex1.index]
                            vertex3.color=temp_color
                            vertex3.edges.append(start_vert)
                            #####Clean the graph structure now that new structure has been added #####
                            temp_graph.clean()
                            #####Assuming no contradiction, recurse on vertex1 #####
                            if not temp_graph.update_graph2(subgraph_list, 1):
                                list2.append(vertex1.index)

                #####Check if a new vertex of that color could be added#####
                list3=[]
                if len(start_vertex.edges) != start_vertex.max_degree:
                    ##### Copy graph and add vertex #####
                    temp_graph = copy.deepcopy(self)
                    temp_graph.depth2= self.depth2+1
                    temp_vert = vert([start_vert], [], max(self.index_list) + 1)
                    temp_vert.color=temp_color
                    temp_graph.vertex_list.append(temp_vert)
                    #####Clean the graph #####
                    temp_graph.clean()
                    #####Assuming no contradiction, recurse on temp_vert #####
                    if not temp_graph.update_graph2(subgraph_list, 1):
                        list3.append(temp_vert.index)

                #####Check for potential neighbors that are already that color (poor order but I forgot and came back to this)#####
                list4=[]
                for edge_index in start_vertex.unknown:
                    vertex1=self.vertex_list[edge_index]
                    if (vertex1.color == temp_color) and (vertex1.index not in done_list2):
                        #####Copy graph and add the edge #####
                        temp_graph = copy.deepcopy(self)
                        temp_graph.depth2= self.depth2+1
                        temp_graph.vertex_list[vertex1.index].edges.append(start_vert)
                        #####Clean the graph #####
                        temp_graph.clean()
                        #####Assuming no contradiction, recurse on vertex1 #####
                        if not temp_graph.update_graph2(subgraph_list, 1):
                            list4.append(vertex1.index)
                                
                ###Check for vertices in the done_list that would be assigned that color in the recoloring strategy that are already neighbors.
                list5=[]
                for index, color in done_list:
                    if color==temp_color and index in start_vertex.edges:
                        list5.append(index)
                
                ###Check for vertices in the done_list that would be assigned that color in the recoloring strategy that could be neighbors.
                list6=[]
                for index, color in done_list:
                    if color==temp_color and index in start_vertex.unknown:
                        #####Copy graph and add the edge #####
                        temp_graph = copy.deepcopy(self)
                        temp_graph.depth2= self.depth2+1
                        temp_graph.vertex_list[index].edges.append(start_vert)
                        #####Clean the graph #####
                        temp_graph.clean()
                        #####Assuming no contradiction, recurse on vertex1 #####
                        if not temp_graph.update_graph2(subgraph_list, 1):
                            list6.append(index)
                
                #####Uncommenting this line will show much more detailed output of the algorithm, although it is admittedly somewhat difficult to read #####
                if verbose:
                    print('vertex: ' + str(start_vert) + ', done_list: ' + str(done_list) + ', colors: ' +  str(start_vertex.color) + '->' + str(temp_color) + ', ' + str([list0, list1, list2, list3, list4, list5, list6]) + ', num_vertices: ' + str(len(self.vertex_list)) )
                #####If there are no ways to force the color, return True #####
                if (len(list0) + len(list1) + len(list2) + len(list3) + len(list4)+ len(list5)+len(list6) ==0):
                    #self.show()
                    return True
                #####If there is only one way to force that color, append that structure #####
                elif (len(list0) + len(list1) + len(list2) + len(list3) + len(list4)+ len(list5)+len(list6) ==1) and self.depth2==0:
                    #####If a neighbor is already colored with temp_color #####
                    if len(list0) > 0:
                        if self.color_argument(subgraph_list, list0[0], done_list + [(start_vert, temp_color)], current_distance + 1, max_distance, max_depth, verbose, add_vertices):
                            return True
                    #####Coloring a vertex #####
                    if len(list1) > 0:
                        vertex1=self.vertex_list[list1[0]]
                        vertex1.color=temp_color
                        self.clean()
                        if self.update_graph2(subgraph_list, 1):
                            return True
                        if self.color_argument(subgraph_list, vertex1.index, done_list + [(start_vert, temp_color)], current_distance + 1, max_distance, max_depth, verbose, add_vertices):
                            return True
                    #####Adding an edge and coloring a vertex #####
                    if len(list2) > 0:
                        vertex1=self.vertex_list[list2[0]]
                        vertex1.color=temp_color
                        vertex1.edges.append(start_vert)
                        self.clean()
                        if self.update_graph2(subgraph_list, 1):
                            return True
                        if self.color_argument(subgraph_list, vertex1.index, done_list + [(start_vert, temp_color)], current_distance + 1, max_distance, max_depth, verbose, add_vertices):
                            return True

                    #####Appending a vertex to the graph, coloring it, and adding an edge #####
                    if len(list3) > 0 and add_vertices:
                        temp_vert = vert([start_vert], [], len(temp_graph.vertex_list))
                        temp_vert.color=temp_color
                        self.vertex_list.append(temp_vert)
                        self.clean()
                        if self.update_graph2(subgraph_list, 1):
                            return True
                        if self.color_argument(subgraph_list, temp_vert.index, done_list + [(start_vert, temp_color)], current_distance + 1, max_distance, max_depth, verbose, add_vertices):
                            return True
                    
                    ###still recurse, but on a copy of the graph if add_vertices is False
                    if len(list3) > 0 and not add_vertices:
                        temp_graph=copy.deepcopy(self)
                        temp_vert = vert([start_vert], [], len(temp_graph.vertex_list))
                        temp_vert.color=temp_color
                        temp_graph.vertex_list.append(temp_vert)
                        temp_graph.clean()
                        if temp_graph.update_graph2(subgraph_list, 1):
                            return True
                        if temp_graph.color_argument(subgraph_list, temp_vert.index, done_list + [(start_vert, temp_color)], current_distance + 1, max_distance, max_depth, verbose, add_vertices):
                            return True

                    #####Adding an edge #####
                    if len(list4) > 0:
                        vertex1=self.vertex_list[list4[0]]
                        vertex1.edges.append(start_vert)
                        self.clean()
                        if self.update_graph2(subgraph_list, 1):
                            return True
                        if self.color_argument(subgraph_list, vertex1.index, done_list + [(start_vert, temp_color)], current_distance + 1, max_distance, max_depth, verbose):
                            return True
                    
                    ###Add an edge ###
                    if len(list6) > 0:
                        vertex1=self.vertex_list[list4[0]]
                        vertex1.edges.append(start_vert)
                        self.clean()
                        if self.update_graph2(subgraph_list, 1):
                            return True
        return False
    
    def color_argument_iter(self, subgraph_list, start_vert, max_distance, verbose=True, add_vertices=True):
        #####This method iteratively applies the prior method, checks for updates, as well as checking for forks between applying coloring arguments. #####
        #####This method hides a lot of the inputs for color_argument2 to make it more 'user friendly' #####
        k=0
        while True:
            k=k+1
            if verbose:
                print('Iteration: ' + str(k))
                print('There are ' + str(len(self.vertex_list)) + ' vertices in the graph.')
            
            old_graph = copy.deepcopy(self)
            #####With the changes to the code, there is not much meaningful difference between the distance and depth. I can probably remove the second variable in the next update but I wanted to get this update out sooner. #####
            if self.color_argument(subgraph_list, start_vert, [], 0, max_distance, max_distance, verbose, add_vertices):
                return True
            if self.check_for_subgraphs_iter(subgraph_list, 1):
                return True
            #####Check for updates
            updated =False
            if len(self.vertex_list)!=len(old_graph.vertex_list):
                updated=True
            if updated==False:
                for vertex1 in self.vertex_list:
                    for vertex2 in old_graph.vertex_list:
                        if (vertex1.index==vertex2.index) and (vertex1.color !=vertex2.color):
                            updated=True
                        if (vertex1.index==vertex2.index) and (vertex1.edges !=vertex2.edges):
                            updated=True
                        if (vertex1.index==vertex2.index) and (vertex1.non_edges !=vertex2.non_edges):
                            updated=True
            if not updated:
                break
        return False
    
    
    def color_argument2(self, subgraph_list, start_vert, done_list, current_distance, max_distance, max_depth, verbose=True):
        #####This method recursively builds out new structure to the graph #####
        ##### done_list is the set of vertices already visited by color_argument2 within the recursion ######
        ##### current_distance is the recursion depth/distance in the graph from the starting vertex
        ##### max_distance is the max distance allowed within the graph
        ##### max_depth is the max recursion depth to which the graph is allowed to be modified
        start_vertex=self.vertex_list[start_vert]
        global color_hash

        #####color_list is colors to be forced #####
        color_list = list(set([0, 1, 2]) - set([start_vertex.color]))

        #####It leads to algorithmic errors if you don't check that you haven't looped back to the start######
        if current_distance > 1:
            for edge_index in start_vertex.edges:
                if (self.vertex_list[edge_index].color_start ==True):
                    return False

        #####Check for each color that it exists among neighbors#########
        done_list2=[pair[0] for pair in done_list]
        if (current_distance < max_distance) and (self.depth2 < max_depth):
            for temp_color in color_list:
                #####Check for already colored neighbors#####
                list0 = []
                list0_2 = []
                for edge_index in start_vertex.edges:
                    vertex1=self.vertex_list[edge_index]
                    if (vertex1.color == temp_color) and (vertex1.index not in done_list2):
                        list0_2.append(vertex1.index)
                        temp_done_list = done_list + [(start_vert, temp_color)]
                        temp_graph = copy.deepcopy(self)
                        temp_graph.depth2= self.depth2+1
                        #####Recurse on vertex1
                        bool1=False
                        #bool1=temp_graph.color_argument2(subgraph_list, vertex1.index, temp_done_list, current_distance + 1, max_distance, max_depth, verbose)
                        if bool1==False:
                            list0.append(vertex1.index)

                #####Check for uncolored neighbors that could be assigned that color#####
                list1=[]
                list1_2=[]
                for edge_index in start_vertex.edges:
                    vertex1=self.vertex_list[edge_index]
                    if (vertex1.color == None) and (vertex1.index not in done_list2):
                        #####Check that vertex1 can be colored with temp_color #####
                        neighbor_colors = []
                        for edge_index2 in vertex1.edges:
                            neighbor_colors.append(self.vertex_list[edge_index2].color)
                        if temp_color not in neighbor_colors:
                            #####Copy graph and color vertex #####
                            temp_graph = copy.deepcopy(self)
                            temp_graph.depth2= self.depth2+1
                            temp_graph.vertex_list[vertex1.index].color=temp_color

                            #####Assuming no contradiction, recurse on vertex1#####
                            if not temp_graph.update_graph2(subgraph_list, 1):
                                list1_2.append(vertex1.index)
                                temp_done_list = done_list + [(start_vert, temp_color)]
                                bool3=False
                                #bool3=temp_graph.color_argument2(subgraph_list, vertex1.index, temp_done_list, current_distance + 1, max_distance, max_depth, verbose)
                                if bool3==False:
                                    list1.append(vertex1.index)

                #####Check for potential neighbors that could be that color#####
                list2=[]
                list2_2=[]
                for edge_index in start_vertex.unknown:
                    vertex1=self.vertex_list[edge_index]
                    if (vertex1.color == None) and (vertex1.index not in done_list2):
                        #####Check that it can be assigned temp_color #####
                        neighbor_colors = []
                        for edge_index2 in vertex1.edges:
                            neighbor_colors.append(self.vertex_list[edge_index2].color)
                        if temp_color not in neighbor_colors:
                            #####Copy graph and add color and the new edge #####
                            temp_graph = copy.deepcopy(self)
                            temp_graph.depth2= self.depth2+1
                            vertex3=temp_graph.vertex_list[vertex1.index]
                            vertex3.color=temp_color
                            vertex3.edges.append(start_vert)
                            #####Clean the graph structure now that new structure has been added #####
                            temp_graph.clean()
                            #####Assuming no contradiction, recurse on vertex1 #####
                            if not temp_graph.update_graph2(subgraph_list, 1):
                                list2_2.append(vertex1.index)
                                temp_done_list = done_list + [(start_vert, temp_color)]
                                bool3=False
                                #bool3=temp_graph.color_argument2(subgraph_list, vertex1.index, temp_done_list, current_distance + 1, max_distance, max_depth, verbose)
                                if bool3==False:
                                    list2.append(vertex1.index)

                #####Check if a new vertex of that color could be added#####
                list3=[]
                list3_2=[]
                if len(start_vertex.edges) != start_vertex.max_degree:
                    ##### Copy graph and add vertex #####
                    temp_graph = copy.deepcopy(self)
                    temp_graph.depth2= self.depth2+1
                    temp_vert = vert([start_vert], [], max(self.index_list) + 1)
                    temp_vert.color=temp_color
                    temp_graph.vertex_list.append(temp_vert)
                    #####Clean the graph #####
                    temp_graph.clean()
                    #####Assuming no contradiction, recurse on temp_vert #####
                    if not temp_graph.update_graph2(subgraph_list, 1):
                        list3_2.append(temp_vert.index)
                        temp_done_list = done_list + [(start_vert, temp_color)]
                        bool3=False
                        #bool3=temp_graph.color_argument2(subgraph_list, temp_vert.index, temp_done_list, current_distance + 1, max_distance, max_depth, verbose)
                        if bool3==False:
                            list3.append(temp_vert.index)

                #####Check for potential neighbors that are already that color (poor order but I forgot and came back to this)#####
                list4=[]
                list4_2=[]
                for edge_index in start_vertex.unknown:
                    vertex1=self.vertex_list[edge_index]
                    if (vertex1.color == temp_color) and (vertex1.index not in done_list2):
                        #####Copy graph and add the edge #####
                        temp_graph = copy.deepcopy(self)
                        temp_graph.depth2= self.depth2+1
                        temp_graph.vertex_list[vertex1.index].edges.append(start_vert)
                        #####Clean the graph #####
                        temp_graph.clean()
                        #####Assuming no contradiction, recurse on vertex1 #####
                        if not temp_graph.update_graph2(subgraph_list, 1):
                            list4_2.append(vertex1.index)
                            temp_done_list = done_list + [(start_vert, temp_color)]
                            bool3=False
                            #bool3=temp_graph.color_argument2(subgraph_list, vertex1.index, temp_done_list, current_distance + 1, max_distance, max_depth, verbose)
                            if bool3==False:
                                list4.append(vertex1.index)
                                
                ###Check for vertices in the done_list that would be assigned that color in the recoloring strategy that are already neighbors.
                list5=[]
                for index, color in done_list:
                    if color==temp_color and index in start_vertex.edges:
                        list5.append(index)
                
                ###Check for vertices in the done_list that would be assigned that color in the recoloring strategy that could be neighbors.
                list6=[]
                for index, color in done_list:
                    if color==temp_color and index in start_vertex.unknown:
                        #####Copy graph and add the edge #####
                        temp_graph = copy.deepcopy(self)
                        temp_graph.depth2= self.depth2+1
                        temp_graph.vertex_list[index].edges.append(start_vert)
                        #####Clean the graph #####
                        temp_graph.clean()
                        #####Assuming no contradiction, recurse on vertex1 #####
                        if not temp_graph.update_graph2(subgraph_list, 1):
                            list6.append(index)
                
                ###Commit the information for this vertex to the hashmap
                color_hash[(start_vertex.index, temp_color, tuple(done_list))]= ([list0_2, list1_2, list2_2, list3_2, list4_2, list5, list6], [list0, list1, list2, list3, list4, list5, list6])
                
                #####Uncommenting this line will show much more detailed output of the algorithm, although it is admittedly pretty difficult to read #####
                if verbose:
                    print('vertex: ' + str(start_vert) + ', done_list: ' + str(done_list) + ', colors: ' +  str(start_vertex.color) + '->' + str(temp_color) + ', depth: ' +  str(self.depth2) + ', ' + str([list0, list1, list2, list3, list4, list5, list6]) + ', num_vertices: ' + str(len(self.vertex_list)) + ', color_hash: ' + str(color_hash[(start_vertex.index, temp_color, tuple(done_list))]) )
                #####If there are no ways to force the color, return True #####
                if (len(list0) + len(list1) + len(list2) + len(list3) + len(list4)+ len(list5)+len(list6) ==0):
                    #self.show()
                    return True
                #####If there is only one way to force that color, append that structure #####
                elif (len(list0) + len(list1) + len(list2) + len(list3) + len(list4)+ len(list5)+len(list6) ==1) and self.depth2==0:
                    #####If a neighbor is already colored with temp_color #####
                    if len(list0) > 0:
                        if self.color_argument2(subgraph_list, list0[0], done_list + [(start_vert, temp_color)], current_distance + 1, max_distance, max_depth, verbose):
                            return True
                    #####Coloring a vertex #####
                    if len(list1) > 0:
                        vertex1=self.vertex_list[list1[0]]
                        vertex1.color=temp_color
                        self.clean()
                        if self.update_graph2(subgraph_list, 1):
                            return True
                        if self.color_argument2(subgraph_list, vertex1.index, done_list + [(start_vert, temp_color)], current_distance + 1, max_distance, max_depth, verbose):
                            return True
                    #####Adding an edge and coloring a vertex #####
                    if len(list2) > 0:
                        vertex1=self.vertex_list[list2[0]]
                        vertex1.color=temp_color
                        vertex1.edges.append(start_vert)
                        self.clean()
                        if self.update_graph2(subgraph_list, 1):
                            return True
                        if self.color_argument2(subgraph_list, vertex1.index, done_list + [(start_vert, temp_color)], current_distance + 1, max_distance, max_depth, verbose):
                            return True

                    #####Appending a vertex to the graph, coloring it, and adding an edge #####
                    if len(list3) > 0:
                        temp_vert = vert([start_vert], [], max(self.index_list) + 1)
                        #print(temp_vert.index)
                        temp_vert.color=temp_color
                        self.vertex_list.append(temp_vert)
                        self.clean()
                        if self.update_graph2(subgraph_list, 1):
                            return True
                        if self.color_argument2(subgraph_list, temp_vert.index, done_list + [(start_vert, temp_color)], current_distance + 1, max_distance, max_depth, verbose):
                            return True

                    #####Adding an edge #####
                    if len(list4) > 0:
                        vertex1=self.vertex_list[list4[0]]
                        vertex1.edges.append(start_vert)
                        self.clean()
                        if self.update_graph2(subgraph_list, 1):
                            return True
                        if self.color_argument2(subgraph_list, vertex1.index, done_list + [(start_vert, temp_color)], current_distance + 1, max_distance, max_depth, verbose):
                            return True
                    
                    ###Add an edge ###
                    if len(list6) > 0:
                        vertex1=self.vertex_list[list4[0]]
                        vertex1.edges.append(start_vert)
                        self.clean()
                        if self.update_graph2(subgraph_list, 1):
                            return True
        return False

    def color_argument2_iter(self, subgraph_list, start_vert, max_distance, verbose=True):
        #####This method iteratively applies the prior method, checks for updates, as well as checking for forks between applying coloring arguments. #####
        #####This method hides a lot of the inputs for color_argument2 to make it more 'user friendly' #####
        k=0
        while True:
            k=k+1
            if verbose:
                print('Iteration: ' + str(k))
                print('There are ' + str(len(self.vertex_list)) + ' vertices in the graph.')
            global color_hash
            #color_hash={}
            old_graph = copy.deepcopy(self)
            #####With the changes to the code, there is not much meaningful difference between the distance and depth. I can probably remove the second variable in the next update but I wanted to get this update out sooner. #####
            if self.color_argument2(subgraph_list, start_vert, [], 0, max_distance, max_distance, verbose):
                return True
            if self.check_for_subgraphs_iter(subgraph_list, 1):
                return True
            #####Check for updates
            updated =False
            if len(self.vertex_list)!=len(old_graph.vertex_list):
                updated=True
            if updated==False:
                for vertex1 in self.vertex_list:
                    for vertex2 in old_graph.vertex_list:
                        if (vertex1.index==vertex2.index) and (vertex1.color !=vertex2.color):
                            updated=True
                        if (vertex1.index==vertex2.index) and (vertex1.edges !=vertex2.edges):
                            updated=True
                        if (vertex1.index==vertex2.index) and (vertex1.non_edges !=vertex2.non_edges):
                            updated=True
            if not updated:
                break
        return False


v_0 = vert([1, 4], [2, 3], 0)
v_1 = vert([0, 2], [3, 4], 1)
v_2 = vert([1, 3], [4, 0], 2)
v_3 = vert([2, 4], [0, 1], 3)
v_4 = vert([3, 0], [1, 2], 4)
C_5=graph([v_0, v_1, v_2, v_3, v_4])
    
v_0 = vert([1, 2], [], 0)
v_1 = vert([0, 2], [], 1)
v_2 = vert([0, 1], [], 2)
triangle = graph([v_0, v_1, v_2])
triangle.clean()

v_0 = vert([1, 2, 3, 4], [5, 6], 0)
v_1 = vert([0, 5], [2, 3, 4, 6], 1)
v_2 = vert([0, 6], [1, 3, 4, 5], 2)
v_3 = vert([0], [1, 2, 4, 5, 6], 3)
v_4 = vert([0], [1, 2, 3, 5, 6], 4)
v_5 = vert([1], [0, 2, 3, 4, 6], 5)
v_6 = vert([2], [0, 1, 3, 4, 5], 6)
fork = graph([v_0, v_1, v_2, v_3, v_4, v_5, v_6])

v_0 = vert([1, 2, 3, 4], [5, 6, 7, 8], 0)
v_1 = vert([0, 6], [2, 3, 4, 5, 7, 8], 1)
v_2 = vert([0, 5], [1, 3, 4, 6, 7, 8], 2)
v_3 = vert([0, 6, 7], [1, 2, 4, 5, 8], 3)
v_4 = vert([0, 5, 8], [1, 2, 3, 6, 7], 4)
v_5 = vert([2, 4, 6], [0, 1, 3, 7, 8], 5)
v_6 = vert([1, 3, 5], [0, 2, 4, 7, 8], 6)
v_7 = vert([3, 8], [0, 1, 2, 4, 5, 6], 7)
v_8 = vert([4, 7], [0, 1, 2, 3, 5, 6], 8)
lem2_3_graph= graph([v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8])

v_0 = vert([1, 2, 3, 4], [5, 6, 7, 8], 0)
v_1 = vert([0, 6], [2, 3, 4, 5, 7, 8], 1)
v_2 = vert([0, 5], [1, 3, 4, 6, 7, 8], 2)
v_3 = vert([0, 6, 7, 8], [1, 2, 4, 5], 3)
v_4 = vert([0, 5, 8], [1, 2, 3, 6, 7], 4)
v_5 = vert([2, 4, 6], [0, 1, 3, 7, 8], 5)
v_6 = vert([1, 3, 5], [0, 2, 4, 7, 8], 6)
v_7 = vert([3], [0, 1, 2, 4, 5, 6, 8], 7)
v_8 = vert([3, 4], [0, 1, 2, 5, 6, 7], 8)
lem2_4_graph= graph([v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8])

v_0 = vert([1, 2, 3, 4], [5, 6, 7, 8], 0)
v_1 = vert([0, 6], [2, 3, 4, 5, 7, 8], 1)
v_2 = vert([0, 5], [1, 3, 4, 6, 7, 8], 2)
v_3 = vert([0, 7, 8], [1, 2, 4, 5, 6], 3)
v_4 = vert([0, 7, 8], [1, 2, 3, 5, 6], 4)
v_5 = vert([2, 6], [0, 1, 3, 4, 7, 8], 5)
v_6 = vert([1, 5], [0, 2, 3, 4, 7, 8], 6)
v_7 = vert([3, 4], [0, 1, 2, 5, 6, 8], 7)
v_8 = vert([3, 4], [0, 1, 2, 5, 6, 7], 8)
lem2_5_graph1=graph([v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8])

v_0 = vert([1, 2, 3, 4], [5, 6, 7, 8], 0)
v_1 = vert([0, 6], [2, 3, 4, 5, 7, 8], 1)
v_2 = vert([0, 5], [1, 3, 4, 6, 7, 8], 2)
v_3 = vert([0, 6, 7, 8], [1, 2, 4, 5], 3)
v_4 = vert([0, 7, 8], [1, 2, 3, 5, 6], 4)
v_5 = vert([2, 6], [0, 1, 3, 4, 7, 8], 5)
v_6 = vert([1, 5, 3], [0, 2, 4, 7, 8], 6)
v_7 = vert([3, 4], [0, 1, 2, 5, 6, 8], 7)
v_8 = vert([3, 4], [0, 1, 2, 5, 6, 7], 8)
lem2_5_graph2=graph([v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8])

v_0 = vert([1, 2, 3, 4], [5, 6, 7, 8], 0)
v_1 = vert([0, 6], [2, 3, 4, 5, 7, 8], 1)
v_2 = vert([0, 5], [1, 3, 4, 6, 7, 8], 2)
v_3 = vert([0, 6, 7, 8], [1, 2, 4, 5], 3)
v_4 = vert([0, 5, 7, 8], [1, 2, 3, 6], 4)
v_5 = vert([2, 4, 6], [0, 1, 3, 7, 8], 5)
v_6 = vert([1, 3, 5], [0, 2, 4, 7, 8], 6)
v_7 = vert([3, 4], [0, 1, 2, 5, 6, 8], 7)
v_8 = vert([3, 4], [0, 1, 2, 5, 6, 7], 8)
lem2_5_graph3=graph([v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8])

v_0 = vert([1, 2, 3, 4], [5, 6, 7, 8], 0)
v_1 = vert([0, 6], [2, 3, 4, 5, 7, 8], 1)
v_2 = vert([0, 5], [1, 3, 4, 6, 7, 8], 2)
v_3 = vert([0, 6, 7, 8], [1, 2, 4, 5], 3)
v_4 = vert([0, 6, 7, 8], [1, 2, 3, 5], 4)
v_5 = vert([2, 6], [0, 1, 3, 4, 7, 8], 5)
v_6 = vert([1, 3, 4, 5], [0, 2, 7, 8], 6)
v_7 = vert([3, 4], [0, 1, 2, 5, 6, 8], 7)
v_8 = vert([3, 4], [0, 1, 2, 5, 6, 7], 8)
lem2_5_graph4=graph([v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8])

###Each list corresponds to the subgraphs you need to forbid to get the statement proven in the lemma and preceding lemmas
lem2_3_list=[fork, triangle, lem2_3_graph]
lem2_4_list=[fork, triangle, lem2_3_graph, lem2_4_graph]
lem2_4_2_list=[fork, triangle, lem2_3_graph, lem2_4_graph, lem2_5_graph2, lem2_5_graph3, lem2_5_graph4]
lem2_5_list=[fork, triangle, lem2_3_graph, lem2_4_graph, lem2_5_graph1, lem2_5_graph2, lem2_5_graph3, lem2_5_graph4] 
