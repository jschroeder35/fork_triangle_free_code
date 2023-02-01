import collections
import copy
from itertools import permutations

# A list of caches to store previous computations.
trielist={}
etrielist={}
graph_hash={}
edgelist_hash={}
subgraph_hash={}
max_hash_size=50000 # size for contains_subgraphs ie graph_hash
max_hash_size2=10000 # size for edgelist_hash and subgraph_hash

class trienode:
    """Node object within a prefix tree for subgraph searching. 
    
    Attributes:
        is_end: Indicates whether a node represent the last vertex of the 
            graph.
        children: Indicates valid edge/non-edge combinations to previously 
            added vertices in the subgraph that correspond to a 'next' 
            vertex towards completing one of the subgraphs in the trie.
    """
    def __init__(self):
        self.is_end = False
        self.children = {}

class trie:
    """Prefix tree object for subgraph searching. Tries of this type
    contain both edges and non-edges within their keys, as opposed to etries,
    which contain only edges while searching.
    
    Attributes:
        root: First node in the trie.
    """
    def __init__(self):
        # Start corresponds to empty subgraph
        self.root=trienode()

    def add_subgraph(self, graph):
        """Adds a subgraph to the trie.
        
        Args:
            subgraph (graph): Subgraph to be added to the prefix tree
        """
        temp=self.root
        for i in range(len(graph.vertex_list)):
            vertex=graph.vertex_list[i]
            # Keep only edges to vertices that precede it in the graph, and 
            # convert to a key for the children dictionary
            temp_edges=[j for j in vertex.edges if j<i]
            temp_non_edges=[j for j in vertex.non_edges if j<i]
            if vertex.fixed_index==False:
                key=(frozenset(temp_edges), frozenset(temp_non_edges))
            else:
                key=(frozenset(temp_edges), frozenset(temp_non_edges), i)
            # Add key to children dictionary if not already present
            if key not in temp.children:
                temp.children[key]=trienode()
            # Move down the trie
            temp=temp.children[key]
            if i==len(graph.vertex_list)-1:
                temp.is_end=True

    def initialize(self, subgraph_list):
        """Adds a subgraph from a subgraph list to a trie object.
        
        Args:
            subgraph_list (list[graph]): List of subgraphs to be added to the 
                prefix tree.
        """
        
        for graph in subgraph_list:
            self.add_subgraph(graph)

    def show(self):
        """Displays all the keys of all of the nodes in the trie.
        """
        
        queue=collections.deque([(self.root, 0)])
        while queue:
            popped=queue.popleft()
            temp=popped[0]
            depth=popped[1]
            for key in temp.children:
                print(key, depth, end=', ')
                queue.append((temp.children[key], depth+1))
                if queue[0][1] > depth:
                    print('')


class etrie:
    """Node object within a prefix tree with keys containing only edges, for 
    the purposes of producing edgelists. 
    
    Attributes:
        is_end: Indicates whether a node represent the last vertex of the 
            graph.
        children: Indicates valid lists of edges to previously added vertices 
            in the subgraph that correspond to a 'next' vertex towards 
            completing one of the subgraphs in the etrie.
    """
    def __init__(self):
        # start corresponds to empty subgraph
        self.root=trienode()
        self.non_edge_list=[]

    def add_subgraph(self, graph):
        """Adds a subgraph to the etrie.
        
        Args:
            subgraph (graph): Subgraph to be added to the prefix tree
        """
        temp=self.root
        for i in range(len(graph.vertex_list)):
            vertex=graph.vertex_list[i]
            # Keep only edges to vertices that precede it in the graph, and 
            # convert to a key for the children dictionary
            temp_edges=[j for j in vertex.edges if j<i]
            temp_non_edges=[j for j in vertex.non_edges if j<i]
            if vertex.fixed_index==False:
                key=frozenset(temp_edges)
            else:
                key=(frozenset(temp_edges), i)
            # Add key to children dictionary if not already present
            if key not in temp.children:
                temp.children[key]=trienode()
            # Move down the trie
            temp=temp.children[key]
            temp.non_edge_list=temp_non_edges
            if i==len(graph.vertex_list)-1:
                temp.is_end=True

    def initialize(self, subgraph_list):
        """Adds a subgraph from a subgraph List to a etrie object.
        
        Args:
            subgraph_list (list[graph]): List of subgraphs to be added to the 
                prefix tree.
        """
        for graph in subgraph_list:
            self.add_subgraph(graph)

    def show(self):
        """Displays all the keys of all of the nodes in the etrie.
        """
        queue=collections.deque([(self.root, 0)])
        while queue:
            popped=queue.popleft()
            temp=popped[0]
            depth=popped[1]
            for key in temp.children:
                print(key, temp.children[key].non_edge_list, depth, end=', ')
                queue.append((temp.children[key], depth+1))
                if queue[0][1] > depth:
                    print('')
                    
class vert: 
    """Object representing a vertex within a graph.
    
    Attributes:
        edges: list of integers corresponding to indices of vertices this 
            vertex has edges with.
        non_edges: list of integers corresponding to indices of vertices this 
            vertex is know to not have edges with.
        unknown: list of integers corresponding to indices of vertices with
            which it is not known whether this vertex has an edge (or not).
        index: integer corresponding to this vertex.
        max_degree: maximum number of edges this vertex may have, if any 
            (default is None).
        color: color assigned to this vertex, if any (default is None).
        color_start: boolean that, if True, denotes that this vertex is the 
            vertex in the graph that is assumed unable to be colored in a 
            partial coloring of the graph (default is False).
        fixed_index: boolean value for use in defining forbidden subgraphs 
            that, if True, indicates that the index of this vertex must be the
            same as the index of the vertex being matched to it in the graph
            being searched (default is False).
        clist: list of possible colors for this vertex, used only within
            algorithms (default is None).
    """
    def __init__(self, edges, non_edges, index):
        """Initializes a vertex object
        
        Args:
            edges (list[int]): List of indices of vertices that this vertex
                has an edge with
            non_edges (list[int]): List of indices of vertices that this vertex
                is known to not have an edge with
            index (int): Index corresponding to the vertex
        """
        self.edges = edges 
        self.non_edges = non_edges 
        self.index= index 
        self.max_degree=None 
        self.color=None 
        self.color_start=False 
        self.clist=None 
        self.fixed_index=False 

class graph:
    """Object representing a vertex within a graph.
    
    Attributes:
        vertex_list (list[vert]): List of vertices in the graph.
        index_list (list[vert]): List of indices of vertices in the graph.
        depth1 (int): Used to track recursion in algorithms (default is 0).
        depth2 (int): Used to track recursion in algorithms (default is 0).
    """
    
    def __init__(self, vertex_list):
        """Initializes a graph object. Note that algorithms assume that 
        vertices have indices from 0 to len(self.vertex_list) and are in order,
        so if this is not the case, most algorithms will fail in unexpected 
        ways. Make sure that vertex indices follow this convention.
        
        Args:
            vertex_list (list[vert]): List of vertices to compose the graph.
        """
        self.vertex_list=vertex_list 
        self.index_list=[vertex.index for vertex in self.vertex_list] 
        self.depth1=0 
        self.depth2=0 
        self.clean()

    def show(self): 
        """Displays all of the vertices in the graph.
        """
        for vertex in self.vertex_list:
            print('Vertex ' + str(vertex.index) + ':')
            print('    Edges:     ' + str(vertex.edges))
            print('    Non-edges: ' + str(vertex.non_edges))
            print('    Unknown:   ' + str(vertex.unknown))
            print('    Color:   ' + str(vertex.color))

    def contains_subgraphs(self, subgraph_list):
        """Checks if the graph contains any subgraphs in subgraph_list. 
        Returns True if at least one subgraph from the list is contained in 
        the graph, otherwise False.
        
        Args:
            subgraph_list (list[graph]): List of subgraphs to be searched for.
        
        Returns:
            bool: True if at least one subgraph from subgraph_list is 
                contained in the graph, else False.
        """

        # Check if graph is in the cache first
        subgraphs_key=[]
        for subgraph in subgraph_list:
            subgraphs_key.append(tuple([
                (frozenset(v.edges), frozenset(v.non_edges)) 
                if v.fixed_index==False 
                else (frozenset(v.edges), frozenset(v.non_edges), v.index) 
                for v in subgraph.vertex_list]))
        subgraphs_key=tuple(subgraphs_key)
        if subgraphs_key not in graph_hash:
            graph_hash[subgraphs_key] = collections.OrderedDict()
        graph_key = tuple([
            (frozenset(v.edges), frozenset(v.non_edges)) 
            for v in self.vertex_list])
        
        # Cache manipulations
        if graph_key in graph_hash[subgraphs_key]:
            graph_hash[subgraphs_key].move_to_end(graph_key)
            return graph_hash[subgraphs_key][graph_key]
        elif len(graph_hash[subgraphs_key])==max_hash_size:
            graph_hash[subgraphs_key].popitem(last=False)
            
        # The method generates a trie from the subgraph list, then 
        # performs a dfs starting from each vertex using the trie
        global trielist
        if subgraphs_key in trielist:
            ttrie=trielist[subgraphs_key]
        else:
            ttrie=trie()
            ttrie.initialize(subgraph_list)
            trielist[subgraphs_key]=ttrie

        # Make a queue to store all remaining positions to search from 
        # while doing the dfs
        queue=collections.deque()
        for vertex in self.vertex_list:
            queue.append((vertex, ttrie.root, []))

        while queue:
            v0, trienode, vlist = queue.pop()
            t_vlist=copy.copy(vlist)
            # Need to convert vertex indices to match the indices of the 
            # subgraph, and filter down to only those corresponding to 
            # vertices already added to our partial list
            t_edges=[i for i, v in enumerate(vlist) if v in v0.edges]
            t_non_edges=[i for i, v in enumerate(vlist) if v in v0.non_edges]
            key1=(frozenset(t_edges), frozenset(t_non_edges))
            key2=(frozenset(t_edges), frozenset(t_non_edges), v0.index)
            # print(v0.index, key1, vlist)
            # Now we can check if the vertex is a valid addition to build out 
            # one of our subgraphs.
            for key in [key1, key2]:
                if key in trienode.children:
                    trienode2=trienode.children[key]
                    t_vlist2=t_vlist+[v0.index]
                    # If a subgraph has been completed, return True
                    if trienode2.is_end==True:
                        graph_hash[subgraphs_key][graph_key]=True
                        graph_hash[subgraphs_key].move_to_end(graph_key)
                        return True
                    # Continue the search on vertices not yet added to the 
                    # subgraph (if it's not completed)
                    for v2 in self.vertex_list:
                        if (v2.index not in t_vlist2):
                            queue.append((v2, trienode2, t_vlist2))
        graph_hash[subgraphs_key][graph_key]=False
        graph_hash[subgraphs_key].move_to_end(graph_key)
        return False

    def get_edgelists(self, subgraph_list):
        """Searches the graph to return a list of lists of edges, such that at 
        least one edge within each list is forced to exist to prevent the 
        subgraphs in the list from occurring in the graph. 
        
        Args:
            subgraph_list (list[graph]): List of forbidden subgraphs.
        
        Returns:
            list[frozenset(frozenset(int))]: List of frozensets of frozensets. 
                Each frozenset within the list contains frozensets with pairs 
                of integers, corresponding to edges in the graph.
        """
        
        # Check if graph is in the cache first
        subgraphs_key=[]
        for subgraph in subgraph_list:
            subgraphs_key.append(tuple([
                (frozenset(v.edges), frozenset(v.non_edges)) 
                if v.fixed_index==False 
                else (frozenset(v.edges), frozenset(v.non_edges), v.index) 
                for v in subgraph.vertex_list]))
        subgraphs_key=tuple(subgraphs_key)
        if subgraphs_key not in edgelist_hash:
            edgelist_hash[subgraphs_key] = collections.OrderedDict()
        graph_key = tuple([
            (frozenset(vertex.edges), frozenset(vertex.non_edges)) 
            for vertex in self.vertex_list])
        
        # Cache manipulations
        if graph_key in edgelist_hash[subgraphs_key]:
            edgelist_hash[subgraphs_key].move_to_end(graph_key)
            return edgelist_hash[subgraphs_key][graph_key]
        elif len(edgelist_hash[subgraphs_key])==max_hash_size2:
            edgelist_hash[subgraphs_key].popitem(last=False)
        
        # Set up prefix tree
        global etrielist
        if subgraphs_key in etrielist:
            ttrie=etrielist[subgraphs_key]
        else:
            ttrie=etrie()
            ttrie.initialize(subgraph_list)
            etrielist[subgraphs_key]=ttrie

        out=[]
        # queue stores positions left to search in the dfs
        queue=collections.deque()
        for vertex in self.vertex_list:
            queue.append((vertex, ttrie.root, [], []))

        while queue:
            v0, trienode, vlist, edgelist = queue.pop()
            # Need to convert vertex indices to match the indices of the 
            # subgraph, and filter down to only those corresponding to 
            # vertices already added to our partial list.
            t_edges=[i for i, v in enumerate(vlist) if v in v0.edges]
            t_non_edges=[i for i, v in enumerate(vlist) if v in v0.non_edges]
            key1=frozenset(t_edges)
            key2=(frozenset(t_edges), v0.index)
            #print(vertex.index, past_vertex_list, key, trienode.children)
            # Now we can check if the vertex is a valid addition to build 
            # out one of our subgraphs.
            for key in [key1, key2]:
                if key in trienode.children:
                    #print(key)
                    trienode2=trienode.children[key]
                    # Need copies of these lists to pass to the queue.
                    t_vlist=copy.copy(vlist)
                    t_edgelist=copy.copy(edgelist)
                    for j in set(trienode2.non_edge_list)-set(t_non_edges):
                        t_edgelist.append(frozenset([vlist[j], v0.index]))
                    t_vlist.append(v0.index)
                    # Check if a subgraph has been completed. If so, 
                    # append the temporary edgelist
                    if trienode2.is_end==True:
                        out.append(t_edgelist)
                    for v2 in self.vertex_list:
                        if v2.index not in t_vlist:
                            queue.append((v2, trienode2, t_vlist, t_edgelist))

        # Now we want to condense the list as much as possible.
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
                        
        # Store output in the cache
        edgelist_hash[subgraphs_key][graph_key]=out
        edgelist_hash[subgraphs_key].move_to_end(graph_key)
        return out

    def get_subgraph_list(self, subgraph):
        """Generates a list of all distinct copies of subgraph in the graph, 
        up to reordering. Not designed to accept input subgraphs containing 
        vertices with fixed_indices.
        
        Args:
            subgraph (graph): Subgraph to be searched for.
        
        Returns:
            list[frozenset(int)]: List of frozensets of integers corresponding
                to copies of the subgraph within the graph.
        """
        
        # Check if graph is in hash map first
        subgraphs_key=tuple([
            (frozenset(v.edges), frozenset(v.non_edges)) 
            if v.fixed_index==False 
            else (frozenset(v.edges), frozenset(v.non_edges), v.index) 
            for v in subgraph.vertex_list])
        if subgraphs_key not in subgraph_hash:
            subgraph_hash[subgraphs_key] = collections.OrderedDict()
        graph_key = tuple([
            (frozenset(v.edges), frozenset(v.non_edges)) 
            for v in self.vertex_list])
        
        # Cache manipulations
        if graph_key in subgraph_hash[subgraphs_key]:
            subgraph_hash[subgraphs_key].move_to_end(graph_key)
            return subgraph_hash[subgraphs_key][graph_key]
        elif len(subgraph_hash[subgraphs_key])==max_hash_size2:
            subgraph_hash[subgraphs_key].popitem(last=False)
        
        # If it's not in the cache, make a trie and do a graph search
        # generate prefix tree
        global trielist
        if subgraphs_key in trielist:
            ttrie=trielist[subgraphs_key]
        else:
            ttrie=trie()
            ttrie.initialize([subgraph])
            trielist[subgraphs_key]=ttrie
        
        # Make a queue to store all remaining positions to search from while 
        # doing the dfs.
        queue=collections.deque()
        out=set()
        for vertex in self.vertex_list:
            queue.append((vertex, ttrie.root, []))

        while queue:
            v0, trienode, vlist = queue.pop()
            t_vlist=copy.copy(vlist)
            # Need to convert vertex indices to match the indices of the 
            # subgraph, and filter down to only those corresponding to 
            # vertices already added to our partial list.
            t_edges=[i for i, v in enumerate(vlist) if v in v0.edges]
            t_non_edges=[i for i, v in enumerate(vlist) if v in v0.non_edges]
            key=(frozenset(t_edges), frozenset(t_non_edges))
            # Now we can check if the vertex is a valid addition to build 
            # out one of our subgraphs.
            if key in trienode.children:
                trienode=trienode.children[key]
                t_vlist.append(v0.index)
                # If a subgraph has been completed, add to the list
                if trienode.is_end==True:
                    out.add(frozenset(t_vlist))
                # Continue the search on vertices not yet added to the 
                # subgraph (if it's not completed)
                for v2 in self.vertex_list:
                    if (v2.index not in t_vlist):
                        queue.append((v2, trienode, t_vlist))
                        
        # Store output in the cache
        subgraph_hash[subgraphs_key][graph_key]=out
        subgraph_hash[subgraphs_key].move_to_end(graph_key)
        return out

    def clean(self):
        """Fixes the edge/non-edge lists of the graph: if v1 has v2 in its edge 
        list, but v2 doesn't have v1 in its edge list, then clean will append 
        v1 to v2's edge list. After all lists have been updated, it will make 
        an 'unknown' list for each vertex, of neither edges nor non-edges. It 
        is recommended to call this method after making changes to a graph.
        """
        
        self.index_list=[vertex.index for vertex in self.vertex_list]
        i=0
        while True:
            i+=1
            if i==1000:
                print('Clean got stuck in a loop. There was an invalid graph' 
                      + ' update. You need to interrupt the kernel and do'
                      + ' some debugging.')
                self.show()
            updated=False
            for v in self.vertex_list:
                # Making sure edge lists match up 
                for v2_index in v.edges:
                    v2=self.vertex_list[v2_index]
                    if v.index not in v2.edges:
                        v2.edges.append(v.index)
                        updated=True
                # Making sure non-edge lists match up 
                for v2_index in v.non_edges:
                    v2=self.vertex_list[v2_index]
                    if v.index not in v2.non_edges:
                        v2.non_edges.append(v.index)
                        updated=True
                # If any vertex becomes max degree, we prohibit any further 
                # edges for that vertex
                if len(v.edges)==v.max_degree:
                    old=copy.copy(v.non_edges)
                    v.non_edges = list(set(self.index_list) - set(v.edges))
                    if set(v.non_edges)!=set(old):
                        updated=True
            # Now we check if the graph was updated
            if not updated:
                break
        # Once all edge/non-edge lists are updated, we set the unknown lists 
        # for each vertex 
        for v in self.vertex_list:
            v.unknown = list(set(self.index_list).difference(
                v.edges + v.non_edges + [v.index]))



    def check_for_subgraphs(self, subgraph_list, depth_cutoff=1):
        """Treating subgraph_list as a list of forbidden subgraphs, searches 
        the graph for any graph in subgraph_list with .contains_subgraphs(), 
        then for each unknown pair of vertices, it tries adding that pair as 
        an edge/non_edge in temporary graphs, then recurses up to depth_cutoff. 
        If the graph contains a forbidden subgraph in both cases, returns 
        True, otherwise, if a forbidden subgraph is found when an edge is 
        added between a pair of vertices, a non-edge is added to the graph, 
        or vice versa.
        
        Args:
            subgraph_list (list[graph]): List of forbidden subgraphs.
            depth_cutoff (int): max recursion depth (default is 1)
        
        Returns:
            bool: True if forbidden subgraph is found to exist in the graph, 
                otherwise False.
        """
        
        if self.contains_subgraphs(subgraph_list):
            return True

        if self.depth1 < depth_cutoff:
            for v in self.vertex_list:
                if len(v.unknown) == 0:
                    continue
                for edge in v.unknown:
                    # Prevents running each case twice
                    if edge > v.index: 
                        # Make a temporary graph with the appended edge
                        temp_graph = copy.deepcopy(self)
                        temp_graph.vertex_list[v.index].edges.append(edge)
                        temp_graph.clean()
                        # Depth tracks the amount of recursion
                        temp_graph.depth1=self.depth1+1 
                        # Make another temp graph with the appended non-edge 
                        temp_graph2 = copy.deepcopy(self)
                        temp_graph2.vertex_list[v.index].non_edges.append(edge)
                        temp_graph2.clean()
                        temp_graph2.depth1=self.depth1+1 

                        # Check both cases and store output
                        bool1 = temp_graph.check_for_subgraphs(
                            subgraph_list, depth_cutoff)
                        bool2 = temp_graph2.check_for_subgraphs(
                            subgraph_list, depth_cutoff)
                        # Contradiction in both cases should return True
                        if (bool1 and bool2):
                            return True
                        # If adding an edge creates a contradiction, add
                        # a non-edge
                        elif bool1:
                            v.non_edges.append(edge)
                            self.clean()
                        # Converse case
                        elif bool2:
                            v.edges.append(edge)
                            self.clean()
        return False

    def check_for_subgraphs_iter(self, subgraph_list, depth_cutoff=1):
        """Runs .check_for_sugraphs() iteratively until the graph is no longer 
        updated or .check_for_subgraphs() returns True.
        
        Args:
            subgraph_list (list[graph]): List of forbidden subgraphs.
            depth_cutoff (int): max recursion depth (default is 1)
        
        Returns:
            bool: True if forbidden subgraph is found to exist in the graph, 
                otherwise False.
        """
    
        # Performance is generally better if, for higher depth_cutoffs, this 
        # method is first run with a depth cutoff of 1
        for i in sorted(list({1, depth_cutoff})): 
            while True:
                old_graph = copy.deepcopy(self)
                bool1=self.check_for_subgraphs(subgraph_list, i)
                if bool1:
                    return True
                updated=False
                for j, v in enumerate(old_graph.vertex_list):
                    if ((v.edges != self.vertex_list[j].edges) or
                            (v.non_edges != self.vertex_list[j].non_edges)):
                        updated = True
                if not updated:
                    break
        return False

    def clean_colors(self):
        """Adds non-edges in accordance with the coloring of the graph.
        """
    
        for tcolor in [0, 1, 2]:
            index_list=[]
            for v in self.vertex_list:
                if v.color==tcolor:
                    v.non_edges = list(set(v.non_edges + index_list))
                    index_list.append(v.index)
        self.clean()

    def color_vert(self, start_vert):
        """Attempts to color a single vertex.
        
        Args:
            start_vert (int): index of the vertex to be colored
        
        Returns:
            bool: True if vertex with index start_vert is not color_start and 
                is unable to be colored, else False.
        """
        
        color_list = []
        uncolored_neighbors = []
        v0=self.vertex_list[start_vert]
        # Look at neighbors to make a list of colors and uncolored neighbors.
        for neighbor in v0.edges:
            vertex = self.vertex_list[neighbor]
            if vertex.color != None:
                color_list.append(vertex.color)
            else:
                uncolored_neighbors.append(vertex.index)
        color_list = list(set(color_list))

        # Color v0, if uncolored 
        if (len(color_list) == 3) and (v0.color_start==False):
            #print('Contradiction!')
            return True
        if (len(color_list) == 2) and (v0.color_start==False):
            if v0.color == None:
                v0.color = list(set([0, 1, 2]) - set(color_list))[0]
            elif (v0.color != list(set([0, 1, 2]) - set(color_list))[0]):
                return True

        # If v0 is the first vertex being colored and it is degree 3, color 
        # all its neighbors with different colors. 
        if v0.color_start==True:
            if ((v0.max_degree ==3) and len(v0.edges)==3 and 
                    (all([w.color==None for w in self.vertex_list]))):
                i=0
                for neighbor in v0.edges:
                    self.vertex_list[neighbor].color=i
                    i=i+1

        # Fix up the graph (non-edges) after adding new colors. 
        self.clean_colors()
        return False

    def fill_in_colors(self, start_vert=None):
        """Iteratively updates the colors of vertices in the graph as much as 
        possible, running .color_vert() on all of the vertices in the graph.
        
        Args:
            start_vert (int): if included, indicates the index of the vertex
                that cannot be colored (.color_start will be set to True). 
                Default is None.
        
        Returns:
            bool: True if .color_vert() returns True for some vertex, else 
                False.
        """

        # Set color_start to True for the start_vert (if not None)
        if start_vert !=None:
            vertex=self.vertex_list[start_vert]
            vertex.color_start = True
            if self.color_vert(start_vert)==True:
                return True

        # Iteratively add colors for the rest of the graph 
        while True:
            old_graph = copy.deepcopy(self)
            # Color the vertex while also checking for a contradiction 
            for v in self.vertex_list:
                if self.color_vert(v.index)==True:
                    return True
            # Check for updates
            updated =False
            for v1 in self.vertex_list:
                if v1.color != old_graph.vertex_list[v1.index].color:
                    updated=True
            if not updated:
                break
        return False

    def update_graph(self, subgraph_list, start_vert=None, with_colors=True):
        """Updates colors and structure (edges and non-edges) of the graph
        as much as possible using .fill_in_colors() and 
        .check_for_subgraphs_iter().
        
        Args:
            subgraph_list (list[graph]): list of forbidden subgraphs.
            start_vert (int): if included, start_vert index to be passed to 
                .fill_in_colors() (default is None)
            with_colors (bool): determines whether to update colors of 
                vertices in the graph as well as structure (Default is True).
        
        Returns:
            bool: True if .fill_in_colors() or .check_for_subgraphs_iter() 
                returns True at any point, else False.
        """

        while True:
            old_graph = copy.deepcopy(self)
            if with_colors:
                if (self.fill_in_colors(start_vert) or 
                        self.check_for_subgraphs_iter(subgraph_list)):
                    return True
            else:
                if self.check_for_subgraphs_iter(subgraph_list):
                    return True
            updated = False
            for v in self.vertex_list:
                v_old=old_graph.vertex_list[v.index]
                if ((v.color != v_old.color) or
                        (set(v.edges) != set(v_old.edges)) or
                        (set(v.non_edges) != set(v_old.non_edges))):
                    updated=True
            if not updated:
                break
        return False

    def update_graph2(self, subgraph_list, start_vert=None, verbose=False, 
            with_colors=True):
        """Updates colors and structure (edges and non-edges) of the graph
        as much as possible by taking the edgelists from .get_edgelists() and,
        for each edgelist, creating a temporary graph for each edge in the
        edgelist, adding the edge to the temporary graph, running 
        .update_graph(), and adding any edges/non_edges that are present in
        all cases that don't return True (a contradiction). If all cases 
        return True, True is returned overall. 
        
        Args:
            subgraph_list (list[graph]): list of forbidden subgraphs.
            start_vert (int): if included, start_vert index to be passed to 
                .fill_in_colors() (default is None)
            verbose (bool): determines whether to print detailed output about
                progress and updates to the graph (default is False).
            with_colors (bool): determines whether to update colors of 
                vertices in the graph as well as structure (Default is True).
        
        Returns:
            bool: True if .update_graph() returns a contradiction, or if a 
                contradiction occurs in all cases implied by an edgelist,
                else False.
        """
        
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
                # Reset the graph_list
                graph_list=[]
                # Iterate over the edges and check their possibilities
                for set1 in edge_list:
                    edge = list(set1)
                    tgraph=copy.deepcopy(self)
                    tgraph.vertex_list[edge[0]].edges.append(edge[1])
                    tgraph.clean()
                    if tgraph.update_graph(subgraph_list, start_vert, 
                            with_colors):
                        continue
                    graph_list.append(tgraph)
                # Indicates there was a contradiction in all possible cases 
                if len(graph_list)==0:
                    return True
                else:
                    for v in self.vertex_list:
                        # Update colors 
                        if with_colors and v.color==None:
                            v.clist=[]
                            for tgraph in graph_list:
                                v_temp=tgraph.vertex_list[v.index]
                                if v_temp.color ==None:
                                    v.clist.append(-1)
                                else:
                                    v.clist.append(v_temp.color)
                            if (-1 not in v.clist) and (len(set(v.clist))==1):
                                v.color=v.clist[0]
                                updated=True
                        # Update Edges 
                        for index1 in v.unknown:
                            edge_list=[]
                            for tgraph in graph_list:
                                v_temp=tgraph.vertex_list[v.index]
                                if index1 in v_temp.edges:
                                    edge_list.append(1)
                                elif index1 in v_temp.non_edges:
                                    edge_list.append(-1)
                                else:
                                    edge_list.append(0)
                            if ((len(set(edge_list))==1) and 
                                    (0 not in edge_list)):
                                if edge_list[0] ==1:
                                    v.edges.append(index1)
                                    updated=True
                                else:
                                    v.non_edges.append(index1)
                                    updated=True
                self.clean()
            if updated and verbose:
                print('Updated, running again...')
            if not updated:
                break
        return False

    def update_graph2_2(self, subgraph_list, max_depth=1, verbose=False, 
            with_colors=True, min_subgraph=None):
        """Extension of .update_graph2() that supports (optionally) obtaining 
        structure/contradictions through adding vertices that are forced 
        through the assumption that a choice of subgraph in the graph has a 
        minimal number of vertices distance 2 from this copy of the subgraph. 
        If a min_subgraph argument is passed, the code requires that the graph 
        is formatted such that the vertices with the smallest indices in the 
        graph correspond to the minimal subgraph. Can recurse (to add multiple 
        vertices forced in this way). If min_subgraph=None, will be the same 
        as update_graph2().
        
        Args:
            subgraph_list (list[graph]): list of forbidden subgraphs.
            max_depth (int): Max recursion depth; ie, number of vertices it 
                can try to add at a time to seek structure/contradictions 
                (default is 1). 
            start_vert (int): if included, start_vert index to be passed to 
                .fill_in_colors() (default is None)
            verbose (bool): determines whether to print detailed output about
                progress and updates to the graph (default is False).
            with_colors (bool): determines whether to update colors of 
                vertices in the graph as well as structure (Default is True).
            min_subgraph (graph): subgraph such that the first 
                len(min_subgraph.vertex_list) vertices in the graph are a 
                copy of this subgraph, and is assumed to be a choice of 
                this subgraph in the graph such that the number of vertices
                of distance exactly 2 is minimal (default is None).
        
        Returns:
            bool: True if .update_graph2() returns a contradiction, or if a 
                contradiction occurs in all cases implied by a min_subgraph,
                else False.
        """

        if min_subgraph != None:
            min_indices=[i for i in range(len(min_subgraph.vertex_list))]

        while True:
            # First run update_graph2
            if self.update_graph2(subgraph_list, start_vert=None, 
                    verbose=verbose, with_colors=with_colors):
                return True
            # At this point, there is no specificed min C_5, we are done.
            if min_subgraph==None or self.depth2>=max_depth:
                return False
            # Otherwise, we obtain a list of 5-cycles and compute a list of 
            # vertices distance 2 from min_indices
            subgraphs=self.get_subgraph_list(min_subgraph)
            edgelists=self.get_edgelists(subgraph_list)
            dist_two_verts=[]
            for v in self.vertex_list:
                # Check first that distance is at least 2 from min_indices
                if set(min_indices).issubset(set(v.non_edges)):
                    # Now, check that distance is at most 2 from min_indices
                    for w in v.edges:
                        w_edges=set(self.vertex_list[w].edges)
                        if len(w_edges.intersection(set(min_indices)))>0:
                            dist_two_verts.append(v.index)
                            break
            updated=False
            if verbose:
                print("Subgraphs: {}".format(subgraphs))
                print("dist_two_verts: {}".format(dist_two_verts))
            # Iterate over our subgraph list for vertices to add
            for tset in subgraphs:
                # Check if there are any vertices in dist_two_verts adjacent 
                # to a vertex in the subgraph
                if verbose:
                    print('Checking subgraph {}'.format(list(tset)))
                d2count=0
                for index in dist_two_verts:
                    t_edges=set(self.vertex_list[index].edges)
                    if len(t_edges.intersection(tset))>0:
                        d2count+=1
                # Now, need to subtract from d2count for each vertex adjacent 
                # (or potentially adjacent) to min_indices but not tset
                for v in self.vertex_list:
                    min_indices_only=set(min_indices)-tset
                    v_pos_edges=set(v.edges+v.unknown)
                    if (len(min_indices_only.intersection(v_pos_edges))>0 and 
                            len(tset.intersection(set(v.edges)))==0):
                        # Check first if an edge between v and tset is implied 
                        # by the edgelists of the graph
                        exists_edge=False
                        for tlist in edgelists:
                            # Check that v is in all edges, and that all edges 
                            # are to tset within the edge list. If that is 
                            # true, v must be adjacent to tset, even though no 
                            # particular edge is forced at the moment.
                            if (all([v.index in pair for pair in tlist]) and 
                                set([list(pair-{v.index})[0] for pair in tlist]
                                   ).issubset(tset)):
                                exists_edge=True
                                break
                        if not exists_edge:
                            d2count-=1
                # If d2count>0, then min_indices does not have the minimal 
                # number of vertices distance 2
                if verbose:
                    print('d2count:{}'.format(d2count))
                if d2count<=0:
                    continue
                # Now, iterate over possible ways to add a vertex to make 
                # min_indices minimal again
                graph_list=[]
                for idx in set(min_indices)-tset:
                    tgraph=copy.deepcopy(self)
                    tgraph.depth2+=1
                    # Add a vertex adjacent to index and non-adjacent to tset
                    new_vert=vert([idx], list(tset), len(tgraph.vertex_list))
                    tgraph.vertex_list.append(new_vert)
                    tgraph.clean()
                    if tgraph.update_graph2_2(subgraph_list, verbose=False, 
                            with_colors=with_colors, min_subgraph=min_subgraph):
                        if verbose:
                            print('{}: True'.format(idx))
                        continue
                    if verbose:
                        print('{}: False'.format(idx))
                    graph_list.append(tgraph)

                # Now we add any structure to the graph and check for updates
                if len(graph_list)==0:
                    return True

                for v in self.vertex_list:
                    # Update colors 
                    if with_colors and v.color==None:
                        v.clist=[]
                        for tgraph in graph_list:
                            v_temp=tgraph.vertex_list[v.index]
                            if v_temp.color ==None:
                                v.clist.append(-1)
                            else:
                                v.clist.append(v_temp.color)
                        if (-1 not in v.clist) and (len(set(v.clist))==1):
                            v.color=v.clist[0]
                            updated=True
                    # Update Edges 
                    for index1 in v.unknown:
                        edge_list=[]
                        for tgraph in graph_list:
                            v_temp=tgraph.vertex_list[v.index]
                            if index1 in v_temp.edges:
                                edge_list.append(1)
                            elif index1 in v_temp.non_edges:
                                edge_list.append(-1)
                            else:
                                edge_list.append(0)
                        if ((len(set(edge_list))==1) and 
                                (0 not in edge_list)):
                            if edge_list[0] ==1:
                                v.edges.append(index1)
                                self.clean()
                                updated=True
                            else:
                                v.non_edges.append(index1)
                                self.clean()
                                updated=True
                                
            if verbose and updated:
                self.show()
                print('Updated (from minimal subgraph); running again...')
            if not updated:
                break
        return False

    def update_graph3(self, subgraph_list, verbose=False, with_colors=True, 
            min_subgraph=None, max_depth=1):
        """Heavier graph update function, capturing more structure by checking 
        every unknown vertex pair to see if adding it as an edge/non-edge will 
        obtain a contradiction with .update_graph2_2().
        
        Args:
            subgraph_list (list[graph]): list of forbidden subgraphs. 
            verbose (bool): determines whether to print detailed output about
                progress and updates to the graph (default is False).
            with_colors (bool): determines whether to update colors of 
                vertices in the graph as well as structure (default is True).
            min_subgraph (graph): subgraph such that, if specified, the first 
                len(min_subgraph.vertex_list) vertices in the graph are a 
                copy of this subgraph, and is assumed to be a choice of 
                this subgraph in the graph such that the number of vertices
                of distance exactly 2 is minimal (default is None).
            max_depth (int): Max recursion depth argument to be passed to 
                .update_graph2_2() (default is 1)
        
        Returns:
            bool: True if a contradiction is forced in some case, else
                False.
        """
        # Stores the index of the last visited vertex, so that as the graph 
        # is updated, vertice can be visited in a cyclical pattern.
        stored_index=0 

        while True:
            if verbose:
                print('Setting up and and initial greedy updates...')
            # Update first with update_graph2_2
            if self.update_graph2_2(subgraph_list, max_depth, verbose=False, 
                    with_colors=with_colors, min_subgraph=min_subgraph):
                if verbose:
                    print('Contradiction obtained in initial updates.')
                    self.show()
                return True
            if verbose:
                self.show()

            # Count number of unknown vertex pairs for printing progress
            if verbose:
                unknown_count=0
                for vertex in self.vertex_list:
                    for v2_index in vertex.unknown:
                        if v2_index>vertex.index:
                            unknown_count+=1
                print('There are {} unknown edges'.format(unknown_count)
                    + ' in the graph to check')

            updated=False
            j=0
            for i in range(len(self.vertex_list)):
                # Picks up from where the last iteration left off
                v=self.vertex_list[(stored_index+i)%len(self.vertex_list)]
                # Want to exit and restart whenever an update occurs.
                if updated or len(v.unknown) == 0:
                    continue
                for edge in v.unknown:
                    if updated or edge < v.index:
                        continue
                    if verbose:
                        completion=100*j/unknown_count
                        print('Checking {}<->{}.'.format(v.index, edge)
                        + ' Completion: {:8.2f}%'.format(completion))
                    j+=1
                    # Make a temporary graph with the appended edge 
                    tgraph = copy.deepcopy(self)
                    tgraph.vertex_list[v.index].edges.append(edge)
                    tgraph.clean()

                    # Make another temporary graph with the appended non-edge 
                    tgraph2 = copy.deepcopy(self)
                    tgraph2.vertex_list[v.index].non_edges.append(edge)
                    tgraph2.clean()

                    # Check both cases for a contradiction
                    bool1 = tgraph.update_graph2_2(subgraph_list, max_depth, 
                            verbose=False, with_colors=with_colors, 
                            min_subgraph=min_subgraph)
                    bool2 = tgraph2.update_graph2_2(subgraph_list, max_depth, 
                            verbose=False, with_colors=with_colors, 
                            min_subgraph=min_subgraph)
                    if (bool1 and bool2):
                        if verbose:
                            print('Contradiction obtained when testing unknown' 
                                + ' edge/non-edge between' 
                                + ' {} and {}!'.format(v.index, edge))
                        return True
                    # In this case, we append a non-edge
                    elif bool1:
                        if verbose:
                            print('Non-edge added from' 
                                + ' {} to {};'.format(v.index, edge)
                                + ' starting over.')
                        v.non_edges.append(edge)
                        self.clean()
                        stored_index=v.index
                        updated=True
                    # In this case, an edge is added 
                    elif bool2:
                        if verbose:
                            print('Edge added from'
                                + ' {} to {};'.format(v.index, edge)
                                + ' starting over.')
                        v.edges.append(edge)
                        self.clean()
                        stored_index=v.index
                        updated=True
            if not updated:
                break
        return False

    def update_graph4(self, subgraph_list, verbose=False, with_colors=True, 
            min_subgraph=None, max_depth=1, max_new_vertices=None):
        """Extension of update_graph3() that will also attempt to add new
        vertices to the graph that are forced to exist by neighborhood 
        containment (an the assumption that the graph is a minimal 
        counterexample) as well as (optionally) new vertices forced to exist
        by the assumption that the first len(min_subgraph.vertex_list) vertices 
        in the graph are a choice of min_subgraph such that the number of
        vertices distance 2 is minimal as well as updating.
        
        Args:
            subgraph_list (list[graph]): list of forbidden subgraphs. 
            verbose (bool): determines whether to print detailed output about
                progress and updates to the graph (default is False).
            with_colors (bool): determines whether to update colors of 
                vertices in the graph as well as structure (Default is True).
            min_subgraph (graph): subgraph such that, if specified, the first 
                len(min_subgraph.vertex_list) vertices in the graph are a 
                copy of this subgraph, and is assumed to be a choice of 
                this subgraph in the graph such that the number of vertices
                of distance exactly 2 is minimal (default is None).
            max_depth (int): Max recursion depth argument to be passed to 
                .update_graph2_2() (default is 1).
            max_new_vertices (int): If specified, will put a maximum on the 
                number of new vertices that can be added to the graph 
                (default is None).
            
        
        Returns:
            bool: True if a contradiction is forced in some case, else
                False.
        """

        new_vert_count=0
        exit=False
        while True:
            updated=False
            if self.update_graph3(subgraph_list, verbose, with_colors, 
                    min_subgraph=min_subgraph, max_depth=max_depth):
                return True

            # If min_subgraph is set, try to add vertices forced by choice of 
            # minimal subgraph first.
            if min_subgraph != None and not exit:
                # Setup
                min_indices=[i for i in range(len(min_subgraph.vertex_list))]
                subgraphs=self.get_subgraph_list(min_subgraph)
                edgelists=self.get_edgelists(subgraph_list)
                dist_two_verts=[]
                for v in self.vertex_list:
                    if set(min_indices).issubset(set(v.non_edges)):
                        dist_two_verts.append(v.index)
                updated=False
                
                # Iterate over our subgraph list for vertices to add
                for tset in subgraphs:
                    if exit:
                        break
                    # Check if there are any vertices in dist_two_verts adjacent 
                    # to a vertex in the subgraph
                    d2count=0
                    for index in dist_two_verts:
                        t_edges=set(self.vertex_list[index].edges)
                        if len(t_edges.intersection(tset))>0:
                            d2count+=1
                            
                    # Now, need to subtract from d2count for each vertex adjacent 
                    # (or potentially adjacent) to min_indices but not tset
                    for v in self.vertex_list:
                        min_indices_only=set(min_indices)-tset
                        v_pos_edges=set(v.edges+v.unknown)
                        if (len(min_indices_only.intersection(v_pos_edges))>0 and 
                                len(tset.intersection(set(v.edges)))==0):
                            # Check first if an edge between v and tset is implied 
                            # by the edgelists of the graph
                            exists_edge=False
                            for tlist in edgelists:
                                # Check that v is in all edges, and that all edges 
                                # are to tset within the edge list. If that is 
                                # true, v must be adjacent to tset, even though no 
                                # particular edge is forced at the moment.
                                if (all([v.index in pair for pair in tlist]) and 
                                    set([list(pair-{v.index})[0] for pair in tlist]
                                       ).issubset(tset)):
                                    exists_edge=True
                                    break
                            if not exists_edge:
                                d2count-=1
                    
                    # If d2count>0, then min_indices does not have the minimal 
                    # number of vertices distance 2
                    if d2count<=0:
                        continue
                    # Now, iterate over possible ways to add a vertex to G to 
                    # make min_indices minimal again
                    index_list=[]
                    for idx in set(min_indices)-tset:
                        tgraph=copy.deepcopy(self)
                        tgraph.depth2+=1
                        # Add a vertex adjacent to idx and non-adjacent to tset
                        new_idx=len(tgraph.vertex_list)
                        new_vert=vert([idx], list(tset), new_idx)
                        tgraph.vertex_list.append(new_vert)
                        tgraph.clean()
                        if tgraph.update_graph2_2(subgraph_list, max_depth, 
                                verbose=False, with_colors=with_colors, 
                                min_subgraph=min_subgraph):
                            continue
                        index_list.append(idx)

                    # Now we return a contradiction or add a vertex if appropriate
                    if len(index_list)==0:
                        if verbose:
                            print('Min subgraph forced a contradiction.')
                        return True
                    if len(index_list)==1:
                        new_idx=len(self.vertex_list)
                        if verbose:
                            print('Vertex {} added adjacent'.format(new_idx)
                                + ' to {} that cannot be'.format(index_list[0])
                                + ' adjacent to {}.'.format(list(tset)))
                        new_vert=vert([index_list[0]], list(tset), new_idx)
                        self.vertex_list.append(new_vert)
                        self.clean()
                        updated=True
                        if self.update_graph2_2(subgraph_list, max_depth, 
                                verbose=False, with_colors=with_colors, 
                                min_subgraph=min_subgraph):
                            return True
                        new_vert_count+=1
                        if (max_new_vertices!=None and 
                                new_vert_count>=max_new_vertices):
                            exit=True
                            
            # Next, check for vertices implied by neighborhood containment.
            for v1 in self.vertex_list:
                if exit:
                    break
                for v2 in self.vertex_list:
                    if exit:
                        break 
                    if (v1.index != v2.index and len(v1.unknown) ==0 and 
                            set(v1.edges).issubset(set(v2.edges))):
                        # Checking that we don't exceed max_degree
                        if (v1.max_degree == None or (v1.max_degree != None and 
                                len(v1.edges) < v1.max_degree)):
                            new_idx=len(self.vertex_list)
                            new_vert=vert([v1.index], [v2.index], new_idx)
                            self.vertex_list.append(new_vert)
                            self.clean()
                            if verbose:
                                print('Vertex {} added'.format(new_vert.index)
                                    + ' adjacent to {} that'.format(v1.index)
                                    + ' cannot be in N({})!'.format(v2.index))
                            updated=True
                            new_vert_count+=1
                            if (max_new_vertices!=None and 
                                    new_vert_count>=max_new_vertices):
                                exit=True
                        # If we would be forced to exceed max_degree, 
                        # return True
                        elif (v1.max_degree != None and 
                                len(v1.edges) >= v1.max_degree):
                            return True
            if not updated:
                break
        return False
    
    def update_graph_flexible(self, subgraph_list, update_speed, verbose=False, 
            min_subgraph=None, with_colors=True, max_depth=1, max_new_vertices=None):
        """Takes an input variable 'update_speed' to determine which graph 
        update function to run, then runs that method.
        
        Args:
            subgraph_list (list[graph]): list of forbidden subgraphs. 
            update_speed (str): Determines which graph update function to run. 
                'fast' corresponds to .update_graph2_2(), 'medium' to 
                .update_graph3(), and 'slow' to .update_graph4().
            verbose (bool): determines whether to print detailed output about
                progress and updates to the graph (default is False).
            with_colors (bool): determines whether to update colors of 
                vertices in the graph as well as structure (Default is True).
            min_subgraph (graph): subgraph such that, if specified, the first 
                len(min_subgraph.vertex_list) vertices in the graph are a 
                copy of this subgraph, and is assumed to be a choice of 
                this subgraph in the graph such that the number of vertices
                of distance exactly 2 is minimal (default is None).
            max_depth (int): Max recursion depth argument to be passed to 
                .update_graph2_2() (default is 1).
            max_new_vertices (int): If specified, will put a maximum on the 
                number of new vertices that can be added to the graph 
                (default is None).
        
        Returns:
            bool: True if the graph update function return True, else False.
        """
        
        if update_speed=='fast':
            if self.update_graph2_2(subgraph_list, max_depth, verbose=verbose, 
                    with_colors=with_colors, min_subgraph=min_subgraph):
                return True
        elif update_speed=='medium':
            if self.update_graph3(subgraph_list, verbose=verbose, 
                    with_colors=with_colors, min_subgraph=min_subgraph, 
                    max_depth=max_depth):
                return True
        elif update_speed=='slow':
            if self.update_graph4(subgraph_list, verbose=verbose, 
                    with_colors=with_colors, min_subgraph=min_subgraph, 
                    max_depth=max_depth, max_new_vertices=max_new_vertices):
                return True
    
    def update_graph5(self, subgraph_list, add_vertices=True, 
            update_speed='medium', max_distance=10, verbose=False, 
            min_subgraph=None, check_deg_three=False, 
            carg_update_speed='fast'):
        """Updates the graph by searching for ways to obtain partial colorings
        implied by the neighborhood of one vertex being contained in the 
        neighborhood of two other vertices, then applying 
        color_argument_iter(). If check_deg_three is set to true, it will also
        check if vertices currently degree 3 must have new neighbors by
        applying color_argument_iter().
        
        Args:
            subgraph_list (list[graph]): list of forbidden subgraphs. 
            add_vertices (bool): Determines if vertices may be added to the 
                graph while the algorithm runs (default is True).
            update_speed (str): Determines which graph update function to run
                while checking for neighborhood containment. 'fast' corresponds 
                to .update_graph2_2(), 'medium' to .update_graph3(), and 'slow'
                to .update_graph4() (default is 'medium').
            carg_update_speed (str): Determines which graph update function to 
                run within color_argument_iter. 'fast', 'medium', and 'slow'
                are options, same as update_speed (default is 'fast').
            verbose (bool): determines whether to print detailed output about
                progress and updates to the graph (default is False).
            min_subgraph (graph): subgraph such that, if specified, the first 
                len(min_subgraph.vertex_list) vertices in the graph are a 
                copy of this subgraph, and is assumed to be a choice of 
                this subgraph in the graph such that the number of vertices
                of distance exactly 2 is minimal (default is None).
            check_deg_three (bool): if True, will check if new neighbors are
                forced for vertices in the graph that are currently degree 3
                (default is False).
            max_distance (int): max distance .color_argument() may recurse 
                (default is 10).
        
        Returns:
            bool: True if a contradiction is forced in all cases, else False.
        """

        # First, we need to clear away any partial coloring in the graph
        for v in self.vertex_list:
            v.color=None
            v.color_start=False

        while True:
            # Tun either update_graph3 or update_graph4 based on whether 
            # adding vertices is allowed
            if add_vertices:
                if self.update_graph4(subgraph_list, verbose=verbose, 
                        with_colors=False, min_subgraph=min_subgraph):
                    return True
            else:
                if self.update_graph3(subgraph_list, verbose=verbose, 
                        with_colors=False, min_subgraph=min_subgraph):
                    return True

            # Now, search for pairs of vertices with neighborhood containment
            old_graph = copy.deepcopy(self)

            for v1, v2, v3 in permutations(self.vertex_list, 3):
                if v2.index>v3.index:
                    continue
                # Check that that N(v1) < N(v2) \cup N(v3)
                pos_edges=set(v1.edges+v1.unknown)
                if not pos_edges.issubset(set(v2.edges+v3.edges)):
                    # For each v4 in pos_edges\{N(v2) + N(v3)}, if both edges 
                    # are unknown, try adding non-edges between v4 and {v2, v3} 
                    # and seeing if this produces a contradiction. If so, 
                    # v4 \sim {v2,v3}, even though neither is currently an edge 
                    # in the graph. 
                    escape=False
                    for vidx in pos_edges-set(v2.edges+v3.edges):
                        v4=self.vertex_list[vidx]
                        if not (v2.index in v4.unknown and 
                                v3.index in v4.unknown):
                            escape=True
                            break
                        tgraph=copy.deepcopy(self)
                        tgraph.vertex_list[vidx].non_edges.append(v2.index)
                        tgraph.vertex_list[vidx].non_edges.append(v3.index)
                        tgraph.clean()
                        if not tgraph.update_graph_flexible(subgraph_list, 
                                update_speed, verbose=False, with_colors=False, 
                                min_subgraph=min_subgraph):
                            escape=True
                            break
                    if escape:
                        continue
                        
                # Add a new neighbor adjacent to v1 not adjacent to v2 or v3 
                # and check if that yields a contradiction
                if v1.max_degree==None or len(v1.edges)<v1.max_degree:
                    tgraph=copy.deepcopy(self)
                    new_idx=len(tgraph.vertex_list)
                    new_vert=vert([v1.index], [v2.index, v3.index], new_idx)
                    tgraph.vertex_list.append(new_vert)
                    tgraph.clean()
                    if verbose:
                        print('Checking if N({})'.format(v1.index)
                            + ' is contained in N({})'.format(v2.index)
                            + ' and N({})'.format(v3.index))
                    # Break unless the update returns a contradiction
                    if not tgraph.update_graph_flexible(subgraph_list, 
                            update_speed, verbose=False, with_colors=False, 
                            min_subgraph=min_subgraph):
                        continue

                # Now we run a coloring argument
                if verbose:
                    print('N({}) is contained in'.format(v1.index)
                        + ' N({}) and N({}).'.format(v2.index, v3.index)
                        + ' Attempting coloring...')
                v2.color=0
                v3.color=1
                v1.color_start=True
                old_num_vertices=len(self.vertex_list)
                if self.fill_in_colors():
                    if verbose:
                        print('Contradiction obtained filling in colors.')
                    return True
                if self.color_argument_iter(subgraph_list, v1.index, 
                        max_distance, verbose=verbose, 
                        add_vertices=add_vertices, min_subgraph=min_subgraph, 
                        carg_update_speed=carg_update_speed):
                    return True
                if len(self.vertex_list)>old_num_vertices:
                    if verbose:
                        print('New vertex added in coloring-updating...')
                    if self.update_graph4(subgraph_list, verbose=False, 
                            with_colors=False, min_subgraph=min_subgraph):
                        return True

                # Reset the colors after
                for v in self.vertex_list:
                    v.color=None
                    v.color_start=False
                
            for v1 in self.vertex_list:
                # Try to color if the vertex is degree 3
                if check_deg_three and len(v1.edges)==3 and len(v1.unknown)==0:
                    # First check if it is forced to have degree only 3 by 
                    # adding a new neighbor and seeing if that creates a 
                    # contradiction.
                    if v1.max_degree != 3:
                        if verbose:
                            print('Checking if {}'.format(v1.index)
                                + ' is degree 3.')
                        tgraph=copy.deepcopy(self)
                        new_vert=vert([v1.index], [], len(tgraph.vertex_list))
                        tgraph.vertex_list.append(new_vert)
                        tgraph.clean()
                        deg_three=False
                        if tgraph.update_graph_flexible(subgraph_list, 
                                update_speed, verbose=False, with_colors=False, 
                                min_subgraph=min_subgraph):
                            deg_three=True
                    else:
                        deg_three=True

                    # In both cases, we try to color, but if the degree is 3, 
                    # then the coloring returns a contradiction if the degree is 
                    # not forced to be 3, then a coloring contradiction forces a 
                    # new neighbor.
                    if deg_three:
                        if verbose:
                            print('{} is degree 3.'.format(v1.index)
                                + ' Attempting coloring...')
                        v1.max_degree=3
                        if self.fill_in_colors(v1.index):
                            if verbose:
                                print('Contradiction obtained filling'
                                    + ' in colors.')
                            return True
                        if self.color_argument_iter(subgraph_list, v1.index, 
                                max_distance, verbose=verbose, 
                                add_vertices=add_vertices, 
                                min_subgraph=min_subgraph, 
                                carg_update_speed=carg_update_speed):
                            return True
                    # This is the case where the degree is not forced to be 3. 
                    elif add_vertices:
                        if verbose:
                            print('{} is not forced to be'.format(v1.index) 
                                + ' degree 3. Attempting coloring to try to'
                                + ' force a new neighbor'
                                + ' of {}.'.format(v1.index))
                        tgraph=copy.deepcopy(self)
                        tgraph.vertex_list[v1.index].max_degree=3
                        add_vert=False
                        if tgraph.fill_in_colors(v1.index):
                            add_vert=True
                        if (not add_vert) and tgraph.color_argument_iter(
                                subgraph_list, v1.index, max_distance, 
                                verbose=verbose, add_vertices=add_vertices, 
                                min_subgraph=min_subgraph, 
                                carg_update_speed=carg_update_speed):
                            add_vert=True
                        # Add a new vertex and update if appropriate
                        if add_vert:
                            if verbose:
                                print('{} Must have degree at'.format(v1.index)
                                    + ' least 4; adding a new neighbor and'
                                    + ' updating...')
                            new_idx=len(self.vertex_list)
                            new_vert=vert([v1.index], [], new_idx)
                            self.vertex_list.append(new_vert)
                            self.clean()
                            if self.update_graph4(subgraph_list, verbose=False, 
                                    with_colors=False, 
                                    min_subgraph=min_subgraph):
                                return True
                            self.show()

                    #Reset the colors after, for either case
                    for v in self.vertex_list:
                        v.color=None
                        v.color_start=False     
                    
            ###check for updates
            updated=False
            if len(self.vertex_list) != len(old_graph.vertex_list):
                updated=True
            else:
                for v in self.vertex_list:
                    v_old=old_graph.vertex_list[v.index]
                    if (set(v.edges) != set(v_old.edges) or
                            set(v.non_edges) != set(v_old.non_edges)):
                        updated=True
            if not updated:
                break

        return False

    def color_argument(self, subgraph_list, start_vert, done_list, 
            current_distance, max_distance, verbose=True, add_vertices=True, 
            min_subgraph=None, carg_update_speed='fast'):
        """Given a starting vertex of index start_vert, for each color i in 
        [0, 1, 2] - color(start_vert), will make 7 lists corresponding to 
        possible neighbors of start_vert that could be color i. the lists
        correspond to [0] current neighbors currently assigned i, [1] current 
        neighbors that could be assigned i, [2] possible neighbors that could 
        be assigned i, [3] new vertex added to the graph colored i, [4] 
        possible neighbors that could be assigned i, [5] current neighbors in 
        done_list that would be reassigned color i by a recoloring strategy, 
        and [6] possible neighbors in done_list that would be reassigned
        color i by a recoloring strategy. If all lists are empty, returns True.
        If there is only 1 element among lists and it is in list 2, 3, 4 or 6, 
        that structure is added to the graph, and if it is list 0, 1, 2, 3, or 
        4, then the algoritm recurses on that vertex, with (start_vert, i)
        appended to done_list. Else, if lists 0, 5, and 6 are empty, then 
        any structure common to all graphs obtained from adding the structure 
        required for lists 1, 2, 3, and 4 is added to the graph before 
        continuing.
        
        Args:
            subgraph_list (list[graph]): list of forbidden subgraphs. 
            start_vert (int): index of the vertex assumed to be unable to be 
                colored.
            done_list (list[tuple(int)]): list of vertex-color pairs, 
                corresponding to the recoloring strategy that would be applied
                to arrive at this vertex.
            current_distance (int): current number of recursion steps.
            max_distance (int): max distance .color_argument() may recurse.
            verbose (bool): determines whether to print detailed output about
                progress and updates to the graph (default is True).
            add_vertices (bool): Determines if vertices may be added to the 
                graph while the algorithm runs (default is True).
            min_subgraph (graph): subgraph such that, if specified, the first 
                len(min_subgraph.vertex_list) vertices in the graph are a 
                copy of this subgraph, and is assumed to be a choice of 
                this subgraph in the graph such that the number of vertices
                of distance exactly 2 is minimal (default is None).
            carg_update_speed (str): Determines which graph update function to 
                run within color_argument_iter. 'fast', 'medium', and 'slow'
                are options, same as update_speed (default is 'fast').
        
        Returns:
            bool: True if the set of all possible neighbors of a required color
                of start_vert is empty, else False.
        """
        
        start_vertex=self.vertex_list[start_vert]

        # color_list is colors to be forced 
        color_list = list(set([0, 1, 2]) - set([start_vertex.color]))

        # Check for each color that it exists among neighbors 
        done_list2=[pair[0] for pair in done_list]
        if (current_distance >= max_distance):
            return False
        
        for tcolor in color_list:
            graph_list=[]
            # Check for already colored neighbors
            list0 = []
            for e_idx in start_vertex.edges:
                v1=self.vertex_list[e_idx]
                if (v1.color == tcolor) and (v1.index not in done_list2):
                    list0.append(v1.index)

            # Check for uncolored neighbors that could be assigned that color
            list1=[]
            for e_idx in start_vertex.edges:
                v1=self.vertex_list[e_idx]
                if (v1.color == None) and (v1.index not in done_list2):
                    # Check that v1 can be colored with tcolor
                    neighbor_colors = []
                    for e_idx2 in v1.edges:
                        neighbor_colors.append(self.vertex_list[e_idx2].color)
                    if tcolor not in neighbor_colors:
                        # Copy graph and color vertex 
                        tgraph = copy.deepcopy(self)
                        tgraph.vertex_list[v1.index].color=tcolor
                        # If no contradiction upon update, append to the 
                        # appropriate list
                        if not tgraph.update_graph_flexible(subgraph_list, 
                                carg_update_speed, min_subgraph=min_subgraph):
                            list1.append(v1.index)
                            graph_list.append(tgraph)

            # Check for potential neighbors that could be that color
            list2=[]
            for e_idx in start_vertex.unknown:
                v1=self.vertex_list[e_idx]
                if (v1.color == None) and (v1.index not in done_list2):
                    #Check that it can be assigned tcolor 
                    neighbor_colors = []
                    for e_idx2 in v1.edges:
                        neighbor_colors.append(self.vertex_list[e_idx2].color)
                    if tcolor not in neighbor_colors:
                        # Copy graph and add color and the new edge 
                        tgraph = copy.deepcopy(self)
                        v3=tgraph.vertex_list[v1.index]
                        v3.color=tcolor
                        v3.edges.append(start_vert)
                        tgraph.clean()
                        # If no contradiction upon update, append to the 
                        # appropriate list
                        if not tgraph.update_graph_flexible(subgraph_list, 
                                carg_update_speed, min_subgraph=min_subgraph):
                            list2.append(v1.index)
                            graph_list.append(tgraph)

            # Check if a new vertex of that color could be added
            list3=[]
            if len(start_vertex.edges) != start_vertex.max_degree:
                # Copy graph and add vertex 
                tgraph = copy.deepcopy(self)
                tvert = vert([start_vert], [], max(self.index_list) + 1)
                tvert.color=tcolor
                tgraph.vertex_list.append(tvert)
                tgraph.clean()
                # If no contradiction upon update, append to the list
                if not tgraph.update_graph_flexible(subgraph_list, 
                        carg_update_speed, min_subgraph=min_subgraph):
                    list3.append(tvert.index)
                    graph_list.append(tgraph)

            # Check for potential neighbors that are already that color 
            list4=[]
            for e_idx in start_vertex.unknown:
                v1=self.vertex_list[e_idx]
                if (v1.color == tcolor) and (v1.index not in done_list2):
                    # Copy graph and add the edge 
                    tgraph = copy.deepcopy(self)
                    tgraph.vertex_list[v1.index].edges.append(start_vert)
                    tgraph.clean()
                    # If no contradiction upon update, append to the 
                    # appropriate list
                    if not tgraph.update_graph_flexible(subgraph_list, 
                            carg_update_speed, min_subgraph=min_subgraph):
                        list4.append(v1.index)
                        graph_list.append(tgraph)

            # Check for vertices in the done_list that would be assigned that 
            # color in the recoloring strategy that are already neighbors.
            list5=[]
            for idx, color in done_list:
                if color==tcolor and idx in start_vertex.edges:
                    list5.append(idx)

            # Check for vertices in the done_list that would be assigned that 
            # color in the recoloring strategy that could be neighbors.
            list6=[]
            for idx, color in done_list:
                if color==tcolor and idx in start_vertex.unknown:
                    # Copy graph and add the edge 
                    tgraph = copy.deepcopy(self)
                    tgraph.vertex_list[idx].edges.append(start_vert)
                    tgraph.clean()
                    # If no contradiction upon update, append to the 
                    # appropriate list
                    if not temp_graph.update_graph_flexible(subgraph_list, 
                            carg_update_speed, min_subgraph=min_subgraph):
                        list6.append(idx)

            if verbose:
                lists= [list0, list1, list2, list3, list4, list5, list6]
                print('Vertex {} is colored'.format(start_vert)
                    + ' {}. Possible neighbors'.format(start_vertex.color)
                    + ' colored {}: {}.'.format(tcolor, lists)
                    + ' done_list: {}'.format(done_list))
            # If there are no ways to force the color, return True 
            if (len(list0) + len(list1) + len(list2) + len(list3) + len(list4) 
                    + len(list5) + len(list6) == 0):
                return True
            # If there is only one way to force tcolor, append that structure 
            # and recurse on that vertex, if appropriate
            elif (len(list0) + len(list1) + len(list2) + len(list3) 
                    + len(list4) + len(list5) + len(list6) == 1):
                # If a neighbor is already colored with tcolor 
                if len(list0) > 0:
                    if self.color_argument(subgraph_list, list0[0], 
                            done_list + [(start_vert, tcolor)], 
                            current_distance + 1, max_distance, verbose, 
                            add_vertices, min_subgraph, carg_update_speed):
                        return True
                # Coloring a vertex 
                if len(list1) > 0:
                    v1=self.vertex_list[list1[0]]
                    v1.color=tcolor
                    self.clean()
                    if self.update_graph_flexible(subgraph_list, 
                            carg_update_speed, min_subgraph=min_subgraph):
                        return True
                    if self.color_argument(subgraph_list, v1.index, 
                            done_list + [(start_vert, tcolor)], 
                            current_distance + 1, max_distance, verbose, 
                            add_vertices, min_subgraph, carg_update_speed):
                        return True
                # Adding an edge and coloring a vertex 
                if len(list2) > 0:
                    v1=self.vertex_list[list2[0]]
                    v1.color=tcolor
                    v1.edges.append(start_vert)
                    self.clean()
                    if self.update_graph_flexible(subgraph_list, 
                            carg_update_speed, min_subgraph=min_subgraph):
                        return True
                    if self.color_argument(subgraph_list, v1.index, 
                            done_list + [(start_vert, tcolor)], 
                            current_distance + 1, max_distance, verbose, 
                            add_vertices, min_subgraph, carg_update_speed):
                        return True
                
                # Appending a vertex to the graph of the appropriate color
                # No recursion in this case
                if len(list3) > 0 and add_vertices:
                    tvert = vert([start_vert], [], len(self.vertex_list))
                    tvert.color=tcolor
                    self.vertex_list.append(tvert)
                    self.clean()
                    if self.update_graph_flexible(subgraph_list, 
                            carg_update_speed, min_subgraph=min_subgraph):
                        return True

                # Adding an edge 
                if len(list4) > 0:
                    v1=self.vertex_list[list4[0]]
                    v1.edges.append(start_vert)
                    self.clean()
                    if self.update_graph_flexible(subgraph_list, 
                            carg_update_speed, min_subgraph=min_subgraph):
                        return True
                    if self.color_argument(subgraph_list, v1.index, 
                            done_list + [(start_vert, tcolor)], 
                            current_distance + 1, max_distance, verbose, 
                            add_vertices, min_subgraph, carg_update_speed):
                        return True

                # Add an edge 
                # No recursion in this case
                if len(list6) > 0:
                    v1=self.vertex_list[list4[0]]
                    v1.edges.append(start_vert)
                    self.clean()
                    if self.update_graph_flexible(subgraph_list, 
                            carg_update_speed, min_subgraph=min_subgraph):
                        return True

            elif (len(list0) + len(list5) + len(list6) ==0):
                # In this case, there are multiple ways to obtain a neighbor of 
                # color tcolor, but all entail adding new structure to the 
                # graph, so we check if any structure is common among all cases 
                # where we add a neighbor of that color. If so, it is added to 
                # the original graph.
                for v in self.vertex_list:
                    # Update colors 
                    if v.color==None:
                        v.clist=[]
                        for tgraph in graph_list:
                            v_temp=tgraph.vertex_list[v.index]
                            if v_temp.color ==None:
                                v.clist.append(-1)
                            else:
                                v.clist.append(v_temp.color)
                        if (-1 not in v.clist) and (len(set(v.clist))==1):
                            v.color=v.clist[0]
                    # Update Edges 
                    for index1 in v.unknown:
                        edge_list=[]
                        for tgraph in graph_list:
                            v_temp=tgraph.vertex_list[v.index]
                            if index1 in v_temp.edges:
                                edge_list.append(1)
                            elif index1 in v_temp.non_edges:
                                edge_list.append(-1)
                            else:
                                edge_list.append(0)
                        if ((len(set(edge_list))==1) and 
                                (0 not in edge_list)):
                            if edge_list[0] ==1:
                                v.edges.append(index1)
                                self.clean()
                            else:
                                v.non_edges.append(index1)
                                self.clean()
                                
        return False

    def color_argument_iter(self, subgraph_list, start_vert, max_distance, 
            verbose=True, add_vertices=True, min_subgraph=None, 
            carg_update_speed='fast'):
        """Given a starting vertex of index start_vert and a max_distance, will
        iteratively run .color_argument() with those starting parameters until
        doing so produces no updates to the graph.
        
        Args:
            subgraph_list (list[graph]): list of forbidden subgraphs. 
            start_vert (int): index of the vertex assumed to be unable to be 
                colored.
            max_distance (int): max distance .color_argument() may recurse.
            verbose (bool): determines whether to print detailed output about
                progress and updates to the graph (default is True).
            add_vertices (bool): Determines if vertices may be added to the 
                graph while the algorithm runs (default is True).
            min_subgraph (graph): subgraph such that, if specified, the first 
                len(min_subgraph.vertex_list) vertices in the graph are a 
                copy of this subgraph, and is assumed to be a choice of 
                this subgraph in the graph such that the number of vertices
                of distance exactly 2 is minimal (default is None).
            carg_update_speed (str): Determines which graph update function to 
                run within color_argument_iter. 'fast', 'medium', and 'slow'
                are options, same as update_speed (default is 'fast').
        
        Returns:
            bool: True if .color_argument() ever returns True, else False.
        """

        k=0
        while True:
            k=k+1
            if verbose:
                num_verts=len(self.vertex_list)
                print('Iteration: {}'.format(k))
                print('There are {} vertices in the graph.'.format(num_verts))

            old_graph = copy.deepcopy(self)
            if self.color_argument(subgraph_list, start_vert, [], 0, 
                    max_distance, verbose, add_vertices, min_subgraph, 
                    carg_update_speed):
                return True
            if self.update_graph_flexible(subgraph_list, carg_update_speed, 
                    min_subgraph=min_subgraph):
                return True
            # Check for updates
            updated =False
            if len(self.vertex_list)!=len(old_graph.vertex_list):
                updated=True
            else:
                for v in self.vertex_list:
                    v_old=old_graph.vertex_list[v.index]
                    if (v.color != v_old.color or 
                            set(v.edges) != set(v_old.edges) or 
                            set(v.non_edges) != set(v_old.non_edges)):
                        updated=True
            if not updated:
                break
        return False


# Defining graphs for use in proofs later
                        
v_0 = vert([1, 4], [2, 3], 0)
v_1 = vert([0, 2], [3, 4], 1)
v_2 = vert([1, 3], [4, 0], 2)
v_3 = vert([2, 4], [0, 1], 3)
v_4 = vert([3, 0], [1, 2], 4)
C_5=graph([v_0, v_1, v_2, v_3, v_4])
C_5.clean()

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
fork.clean()

v_0 = vert([1, 2, 3, 4], [5, 6, 7, 8, 9], 0)
v_1 = vert([0, 6], [2, 3, 4, 5, 7, 8, 9], 1)
v_2 = vert([0, 5], [1, 3, 4, 6, 7, 8, 9], 2)
v_3 = vert([0, 6, 7], [1, 2, 4, 5, 8, 9], 3)
v_4 = vert([0, 5, 8, 9], [1, 2, 3, 6, 7], 4)
v_5 = vert([2, 4, 6], [0, 1, 3, 7, 8, 9], 5)
v_6 = vert([1, 3, 5], [0, 2, 4, 7, 8, 9], 6)
v_7 = vert([3, 8, 9], [0, 1, 2, 4, 5, 6], 7)
v_8 = vert([4, 7], [0, 1, 2, 3, 5, 6, 9], 8)
v_9 = vert([4, 7], [0, 1, 2, 3, 5, 6, 8], 9)
lem2_3_graph0= graph([v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8, v_9])
lem2_3_graph0.clean()

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
lem2_3_graph.clean()

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
lem2_4_graph.clean()

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
lem2_5_graph1.clean()

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
lem2_5_graph2.clean()

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
lem2_5_graph3.clean()

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
lem2_5_graph4.clean()

# Each list corresponds to the subgraphs you need to forbid to get the 
# statement proven in the lemma and preceding lemmas
lem2_3_list0=[fork, triangle, lem2_3_graph0]
lem2_3_list=[fork, triangle, lem2_3_graph]
lem2_4_list=[fork, triangle, lem2_3_graph, lem2_4_graph]
lem2_4_2_list=[fork, triangle, lem2_3_graph, lem2_4_graph, lem2_5_graph2, 
        lem2_5_graph3, lem2_5_graph4]
lem2_5_list=[fork, triangle, lem2_3_graph, lem2_4_graph, lem2_5_graph1, 
        lem2_5_graph2, lem2_5_graph3, lem2_5_graph4]

# For Lemma 3_5, the first 5 indices are fixed, since the choice of 5-cycle 
# must be minimal within this subgraph. 
lem3_4_list=[]
for j in range(5):
    v_0 = vert([1, 4], [], 0)
    v_1 = vert([2, 0], [], 1)
    v_2 = vert([3, 1], [], 2)
    v_3 = vert([4, 2], [], 3)
    v_4 = vert([0, 3], [], 4)
    v_5 = vert([j, (3+j)%5, 7], [], 5)
    v_6 = vert([j, (2+j)%5], [7], 6)
    v_7 = vert([5], [0, 1, 2, 3, 4], 7)
    lem3_4_graph=graph([v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7])
    lem3_4_graph.clean()
    lem3_4_graph.update_graph(lem2_5_list)
    for i in range(5):
        lem3_4_graph.vertex_list[i].fixed_index=True
    lem3_4_list.append(lem3_4_graph)
    
    v_0 = vert([1, 4], [], 0)
    v_1 = vert([2, 0], [], 1)
    v_2 = vert([3, 1], [], 2)
    v_3 = vert([4, 2], [], 3)
    v_4 = vert([0, 3], [], 4)
    v_5 = vert([j, (2+j)%5, 7], [], 5)
    v_6 = vert([j, (3+j)%5], [7], 6)
    v_7 = vert([5], [0, 1, 2, 3, 4], 7)
    lem3_4_graph=graph([v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7])
    lem3_4_graph.clean()
    lem3_4_graph.update_graph(lem2_5_list)
    for i in range(5):
        lem3_4_graph.vertex_list[i].fixed_index=True
    lem3_4_list.append(lem3_4_graph)

lem3_4_list0=lem3_4_list
lem3_4_list=lem3_4_list+lem2_5_list

lem3_5_list0=[]
for j in range(5):
    v_0 = vert([1, 4], [], 0)
    v_1 = vert([2, 0], [], 1)
    v_2 = vert([3, 1], [], 2)
    v_3 = vert([4, 2], [], 3)
    v_4 = vert([0, 3], [], 4)
    v_5 = vert([j, (3+j)%5, 7], [], 5)
    v_6 = vert([j, (2+j)%5, 7], [], 6)
    v_7 = vert([5, 6], [0, 1, 2, 3, 4], 7)
    lem3_5_graph=graph([v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7])
    lem3_5_graph.clean()
    lem3_5_graph.update_graph(lem2_5_list)
    for i in range(5):
        lem3_5_graph.vertex_list[i].fixed_index=True
    lem3_5_list0.append(lem3_5_graph)
lem3_5_list0=lem3_5_list0+lem3_4_list

lem3_5_list=[]
for j in range(5):
    v_0 = vert([1, 4], [], 0)
    v_1 = vert([2, 0], [], 1)
    v_2 = vert([3, 1], [], 2)
    v_3 = vert([4, 2], [], 3)
    v_4 = vert([0, 3], [], 4)
    v_5 = vert([j, (2+j)%5, 7], [(3+j)%5], 5)
    v_6 = vert([j, (2+j)%5, 7], [], 6)
    v_7 = vert([5, 6], [0, 1, 2, 3, 4], 7)
    lem3_5_graph=graph([v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7])
    lem3_5_graph.clean()
    lem3_5_graph.update_graph(lem2_5_list)
    for i in range(5):
        lem3_5_graph.vertex_list[i].fixed_index=True
    lem3_5_list.append(lem3_5_graph)
    
    v_0 = vert([1, 4], [], 0)
    v_1 = vert([2, 0], [], 1)
    v_2 = vert([3, 1], [], 2)
    v_3 = vert([4, 2], [], 3)
    v_4 = vert([0, 3], [], 4)
    v_5 = vert([j, 7], [(2+j)%5, (3+j)%5], 5)
    v_6 = vert([j, (2+j)%5, 7], [], 6)
    v_7 = vert([5, 6], [0, 1, 2, 3, 4], 7)
    lem3_5_graph=graph([v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7])
    lem3_5_graph.clean()
    lem3_5_graph.update_graph(lem2_5_list)
    for i in range(5):
        lem3_5_graph.vertex_list[i].fixed_index=True
    lem3_5_list.append(lem3_5_graph)
    
    v_0 = vert([1, 4], [], 0)
    v_1 = vert([2, 0], [], 1)
    v_2 = vert([3, 1], [], 2)
    v_3 = vert([4, 2], [], 3)
    v_4 = vert([0, 3], [], 4)
    v_5 = vert([j, (3+j)%5, 7], [], 5)
    v_6 = vert([j, (3+j)%5, 7], [(2+j)%5], 6)
    v_7 = vert([5, 6], [0, 1, 2, 3, 4], 7)
    lem3_5_graph=graph([v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7])
    lem3_5_graph.clean()
    lem3_5_graph.update_graph(lem2_5_list)
    for i in range(5):
        lem3_5_graph.vertex_list[i].fixed_index=True
    lem3_5_list.append(lem3_5_graph)
    
    v_0 = vert([1, 4], [], 0)
    v_1 = vert([2, 0], [], 1)
    v_2 = vert([3, 1], [], 2)
    v_3 = vert([4, 2], [], 3)
    v_4 = vert([0, 3], [], 4)
    v_5 = vert([j, (3+j)%5, 7], [], 5)
    v_6 = vert([j, 7], [(2+j)%5, (3+j)%5], 6)
    v_7 = vert([5, 6], [0, 1, 2, 3, 4], 7)
    lem3_5_graph=graph([v_0, v_1, v_2, v_3, v_4, v_5, v_6, v_7])
    lem3_5_graph.clean()
    lem3_5_graph.update_graph(lem2_5_list)
    for i in range(5):
        lem3_5_graph.vertex_list[i].fixed_index=True
    lem3_5_list.append(lem3_5_graph)

lem3_5_list=lem3_5_list+lem3_5_list0

lem3_6_list=[]
for j in range(5):
    v_0 = vert([1, 4], [], 0)
    v_1 = vert([2, 0], [], 1)
    v_2 = vert([3, 1], [], 2)
    v_3 = vert([4, 2], [], 3)
    v_4 = vert([0, 3], [], 4)
    v_5 = vert([j, (2+j)%5], [], 5)
    v_6 = vert([5], [0, 1, 2, 3, 4], 6)
    lem3_6_graph=graph([v_0, v_1, v_2, v_3, v_4, v_5, v_6])
    lem3_6_graph.clean()
    lem3_6_graph.update_graph(lem2_5_list)
    for i in range(5):
        lem3_6_graph.vertex_list[i].fixed_index=True
    lem3_6_list.append(lem3_6_graph)
lem3_6_list=lem3_6_list + lem3_5_list

lem3_7_list0=[]
for j in range(1, 5):
    v_0 = vert([1, 4], [], 0)
    v_1 = vert([2, 0], [], 1)
    v_2 = vert([3, 1], [], 2)
    v_3 = vert([4, 2], [], 3)
    v_4 = vert([0, 3], [], 4)
    v_5 = vert([j, (2+j)%5], [], 5)
    lem3_7_graph=graph([v_0, v_1, v_2, v_3, v_4, v_5])
    lem3_7_graph.clean()
    lem3_7_graph.update_graph(lem2_5_list)
    for i in range(5):
        lem3_7_graph.vertex_list[i].fixed_index=True
    lem3_7_list0.append(lem3_7_graph)
lem3_7_list0=lem3_7_list0 + lem3_6_list

lem3_7_list=[]
for j in range(5):
    v_0 = vert([1, 4], [], 0)
    v_1 = vert([2, 0], [], 1)
    v_2 = vert([3, 1], [], 2)
    v_3 = vert([4, 2], [], 3)
    v_4 = vert([0, 3], [], 4)
    v_5 = vert([j, (2+j)%5], [], 5)
    lem3_7_graph=graph([v_0, v_1, v_2, v_3, v_4, v_5])
    lem3_7_graph.clean()
    lem3_7_graph.update_graph(lem2_5_list)
    for i in range(5):
        lem3_7_graph.vertex_list[i].fixed_index=True
    lem3_7_list.append(lem3_7_graph)
lem3_7_list=lem3_7_list + lem2_5_list

