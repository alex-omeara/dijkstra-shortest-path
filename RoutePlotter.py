# -*- coding: utf-8 -*-
# 
# @author Alex O'Meara 119442536


class Vertex:
    """ A Vertex in a graph. """

    def __init__(self, element):
        """ Create a vertex, with a data element.

        Args:
            element - the data or label to be associated with the vertex
        """
        self._element = element

    def __str__(self):
        """ Return a string representation of the vertex. """
        return str(self._element)

    def __lt__(self, v):
        """ Return true if this element is less than v's element.

        Args:
            v - a vertex object
        """
        return self._element < v.element()

    def element(self):
        """ Return the data for the vertex. """
        return self._element


class Edge:
    """ An edge in a graph.

        Implemented with an order, so can be used for directed or undirected
        graphs. Methods are provided for both. It is the job of the Graph class
        to handle them as directed or undirected.
    """

    def __init__(self, v, w, element):
        """ Create an edge between vertices v and w, with a data element.

        Element can be an arbitrarily complex structure.

        Args:
            element - the data or label to be associated with the edge.
        """
        self._vertices = (v, w)
        self._element = element

    def __str__(self):
        """ Return a string representation of this edge. """
        return ('(' + str(self._vertices[0]) + '--'
                + str(self._vertices[1]) + ' : '
                + str(self._element) + ')')

    def vertices(self):
        """ Return an ordered pair of the vertices of this edge. """
        return self._vertices

    def start(self):
        """ Return the first vertex in the ordered pair. """
        return self._vertices[0]

    def end(self):
        """ Return the second vertex in the ordered pair. """
        return self._vertices[1]

    def opposite(self, v):
        """ Return the opposite vertex to v in this edge.

        Args:
            v - a vertex object
        """
        if self._vertices[0] == v:
            return self._vertices[1]
        elif self._vertices[1] == v:
            return self._vertices[0]
        else:
            return None

    def element(self):
        """ Return the data element for this edge. """
        return self._element


class Graph:
    """ Represent a simple graph.

    This version maintains only undirected graphs, and assumes no
    self loops.

    Implements the Adjacency Map style. Also maintains a top level
    dictionary of vertices.
    """

    #Implement as a Python dictionary
    #  - the keys are the vertices
    #  - the values are the sets of edges for the corresponding vertex.
    #    Each edge set is also maintained as a dictionary,
    #    with the opposite vertex as the key and the edge object as the value.

    def __init__(self):
        """ Create an initial empty graph. 
            open - APQ object of all vertices not yet visited 
            locs - dict where a key is a vertex object and a value is a reference to the object in open
            closed - dict where a key is a vertex and a value is a tuple of an element object's key, 
                and the previous element object used to reach it
            preds - dict where a key is a vertex object and a value is the vertex used to reach it
            graph - graph object of all vertices and edges to be operated on"""
            
        self._structure = dict()
        self._open = AdaptablePriorityQueue()
        self._locs = {}
        self._closed = {}
        self._preds = {}

    def __str__(self):
        """ Return a string representation of the graph. """
        hstr = ('|V| = ' + str(self.num_vertices())
                + '; |E| = ' + str(self.num_edges()))
        vstr = '\nVertices: '
        for v in self._structure:
            vstr += str(v) + '-'
        edges = self.edges()
        estr = '\nEdges: '
        for e in edges:
            estr += str(e) + ' '
        return hstr + vstr + estr

    #-----------------------------------------------------------------------#

    # ADT methods to query the graph

    def num_vertices(self):
        """ Return the number of vertices in the graph. """
        return len(self._structure)

    def num_edges(self):
        """ Return the number of edges in the graph. """
        num = 0
        for v in self._structure:
            num += len(self._structure[v])    # the dict of edges for v
        return num // 2     # divide by 2, since each edege appears in the
        # vertex list for both of its vertices

    def vertices(self):
        """ Return a list of all vertices in the graph. """
        return [key for key in self._structure]

    def get_vertex_by_label(self, element):
        """ Return the first vertex that matches element. """
        for v in self._structure:
            if v.element() == element:
                return v
        return None

    def edges(self):
        """ Return a list of all edges in the graph. """
        edgelist = []
        for v in self._structure:
            for w in self._structure[v]:
                #to avoid duplicates, only return if v is the first vertex
                if self._structure[v][w].start() == v:
                    edgelist.append(self._structure[v][w])
        return edgelist

    def get_edges(self, v):
        """ Return a list of all edges incident on v.

        Args:
            v - a vertex object
        """
        if v in self._structure:
            edgelist = []
            for w in self._structure[v]:
                edgelist.append(self._structure[v][w])
            return edgelist
        return None

    def get_edge(self, v, w):
        """ Return the edge between v and w, or None.

        Args:
            v - a vertex object
            w - a vertex object
        """
        if (self._structure is not None
            and v in self._structure
                and w in self._structure[v]):
            return self._structure[v][w]
        return None

    def degree(self, v):
        """ Return the degree of vertex v.

        Args:
            v - a vertex object
        """
        return len(self._structure[v])

    #----------------------------------------------------------------------#

    # ADT methods to modify the graph

    def add_vertex(self, element):
        """ Add a new vertex with data element.

        If there is already a vertex with the same data element,
        this will create another vertex instance.
        """
        v = Vertex(element)
        self._structure[v] = dict()
        return v

    def add_vertex_if_new(self, element):
        """ Add and return a vertex with element, if not already in graph.

        Checks for equality between the elements. If there is special
        meaning to parts of the element (e.g. element is a tuple, with an
        'id' in cell 0), then this method may create multiple vertices with
        the same 'id' if any other parts of element are different.

        To ensure vertices are unique for individual parts of element,
        separate methods need to be written.

        """
        for v in self._structure:
            if v.element() == element:
                return v
        return self.add_vertex(element)

    def add_edge(self, v, w, element):
        """ Add and return an edge between two vertices v and w, with  element.

        If either v or w are not vertices in the graph, does not add, and
        returns None.
            
        If an edge already exists between v and w, this will
        replace the previous edge.

        Args:
            v - a vertex object
            w - a vertex object
            element - a label
        """
        if v not in self._structure or w not in self._structure:
            return None
        e = Edge(v, w, element)
        self._structure[v][w] = e
        self._structure[w][v] = e
        return e

    def add_edge_pairs(self, elist):
        """ add all vertex pairs in elist as edges with empty elements.

        Args:
            elist - a list of pairs of vertex objects
        """
        for (v, w) in elist:
            self.add_edge(v, w, None)

    #---------------------------------------------------------------------#

    # Additional methods to explore the graph

    def highestdegreevertex(self):
        """ Return the vertex with highest degree. """
        hd = -1
        hdv = None
        for v in self._structure:
            if self.degree(v) > hd:
                hd = self.degree(v)
                hdv = v
        return hdv

    def dijkstra(self, vertexLabel):
        """Return the closed dict of vertices visited by dijkstra's algorithm 
        
        Args:
            vertexLabel - the entry of a vertex in graph object
        """
        srcVertex = self.get_vertex_by_label(vertexLabel)
        if srcVertex:
            vertexDict = {srcVertex: self._structure.get(srcVertex)}
            self._preds = {srcVertex: None}
            srcElement = self._open.add(0, vertexDict)
            self._locs[srcVertex] = srcElement

            while self._open._heap:
                minElem = self._open.removeMin()
                elemVertex = next(iter(minElem._value))
                self._locs.pop(elemVertex)
                pred = self._preds.pop(elemVertex)
                self._closed[elemVertex] = (minElem._key, pred)

                for edge in self.get_edges(elemVertex):
                    oppVertex = edge.opposite(elemVertex)

                    if oppVertex not in self._closed:
                        newCost = minElem._key + edge._element

                        if oppVertex not in self._locs:
                            self._preds[oppVertex] = elemVertex
                            vertexDict = {
                                oppVertex: self._structure.get(oppVertex)}
                            self._open.add(newCost, vertexDict)
                            self._locs[oppVertex] = minElem

                        elif newCost < self._open.getElement(oppVertex)._key:
                            self._preds[oppVertex] = elemVertex
                            predElement = self._open.getElement(oppVertex)
                            self._open.updateKey(predElement, newCost)

            return self._closed

    # End of class definition

#---------------------------------------------------------------------------#


class Element:
    """ An element in an adaptable priority queue """

    def __init__(self, key, value, index):
        """ Create an element, with a key, value, and index

        Args:
            key - the cost to reach the element
            value - the vertex sub dictionary
            index - the index of the element in the APQ
        """
        self._key = key  # cost to get to that vertex
        self._value = value
        self._index = index
    
    def __str__(self):
        """ Return a string representation of the element """
        return "(" + str(self._key) + ", " + str(self._value) + ", " + str(self._index) + ") "
    
    def __eq__(self, value):
        """ Return true if the element's key is equal to value's key 
        
        Args:
            value - an element object
        """
        return self._key == value._key
    
    def __lt__(self, value):
        """ Return true if this element's key is less than values's key.

        Args:
            value - an element object
        """
        return self._key < value._key

class AdaptablePriorityQueue:
    """ A min binary heap representation of an adaptable priority queue """

    def __init__(self):
        """ Create an adaptable priority queue (APQ) with a heap """
        self._heap = []
    
    def __str__(self):
        """ Return a string representation of the APQ """
        output = str(self._heap[0])
        i = 2
        for num in range(1, 7):
            if num == i-1:
                i *= 2
                output += "\n"   
            output += str(self._heap[num])
        return output
  
    def getParent(self, element):
        """ Return the parent of element or None if it has no parent 
        
        Args:
            element - element object in the APQ
        """
        if (element._index - 1) // 2 <= 0:
            return None
        else:
            parent = self._heap[(element._index-1)//2]
            return self._heap[(element._index-1)//2]
    
    def getLeftChild(self, element):
        """ Return the left child of element or None if it has no left child
        
        Args:
            element - element object in the APQ
        """
        if len(self._heap) <= (element._index * 2) + 1:
            return None
        return self._heap[(2 * element._index) + 1]
    
    def getRightChild(self, element):
        """ Return the right child of element or None if it has no right child
        
        Args:
            element - element object in the APQ
        """
        if len(self._heap) <= (element._index * 2) + 2:
            return None
        return self._heap[(2 * element._index) + 2]
    
    def minChild(self, element):
        """ Return the highest priority child of element 
        
        Args:
            element - element object in the APQ
        """
        leftChild = self.getLeftChild(element)
        rightChild = self.getRightChild(element)
        if leftChild and rightChild:
            if leftChild < rightChild:
                return leftChild
            else:
                return rightChild
        elif leftChild:
            return leftChild
        else:
            return rightChild
        return None

    def getKey(self, element):
        """ Return the key of element
        
        Args:
            element - element object in the APQ
         """
        return element._key
    
    def getElement(self, vertex):
        """ Return element of vertex in the AQP

        Args:
            vertex - vertex object
        """
        for element in self._heap:
            if next(iter(element._value)) == vertex:
                # if the key of the value dict == vertex
                return element
        return None

    def min(self):
        """ Return the element on top of the heap """
        return self._heap[0]

    def add(self, key, item):
        """ Return the element added to the APQ 

        Args:
            key - priority of the element to be added
            item - the vertex dictionary to added to the element
        """
        element = Element(key, item, len(self._heap))
        # create the element object
        self._heap.append(element)
        parent = self.getParent(element)
        if parent and element < parent:
            # if the element has a parent and has a higher priority than its parent
           return self.bubbleUp(element)
        else:
            return element
    
    def bubbleUp(self, element):
        """ Return the element after it has been recursively bubbled up the heap
        
        Args:
            element - element object in APQ to be bubbled up
        """
        parent = self.getParent(element)
        if parent and element < parent:
            # update element and parents indices
            element._index, parent._index = parent._index, element._index
            # swaps element and elements parent in heap
            self._heap[element._index], self._heap[parent._index] = element, parent
            self.bubbleUp(element)
        return element

    def removeMin(self):
        """ Return, and remove the element on top of the heap or None if the heap is empty """
        if self._heap:
            minElem = self.min()
            element = self._heap.pop()
            # get element at bottom of heap
            if len(self._heap) > 0:
                element._index = 0
                self._heap[0] = element
                # swap element at bottom of heap into top
                self.bubbleDown(element)
            return minElem
        else:
            return None
    
    def remove(self, element):
        """ Return, and remove the specified element from the APQ or none if it is not in the APQ 
        
        Args:
            element - element object in APQ to be removed
        """
        if element in self._heap:
            lastElem = self._heap.pop()
            # pop last elemnt in heap
            self._heap[element._index] = lastElem
            # put last element into where element is
            if lastElem < self.getParent(lastElem):
                self.bubbleUp(lastElem)
            else:
                self.bubbleDown(lastElem)
            return element
        else:
            return None
    
    def bubbleDown(self, element):
        """ Return the element after it has been recursively bubbled down the heap
        
        Args:
            element - element object in APQ to be bubbled down
        """
        minChild = self.minChild(element)
        if minChild and minChild < element:
            # if the element has lower priority than its child
            element._index, minChild._index = minChild._index, element._index
            # swap element's index, and child's index
            self._heap[element._index], self._heap[minChild._index] = element, minChild
            # swap element, and child in heap
            self.bubbleDown(element)
        return element
    
    def updateKey(self, element, newKey):
        """ Update the key of a given element

        Args:
            element - element object in APQ
            newKey - new key of element
        """
        element._key = newKey
        if self.getParent(element) and element < self.getParent(element):
            return self.bubbleUp(element)
        else:
            return self.bubbleDown(element)


class RouteMap(Graph):
    """ A simple undirected weighted subclass of graph """

    def __init__(self):
        """ Create an inherited empty graph where coords is a dict with a vertex object key, 
            and tuple of floating point numbers latitude and longitude.
            Vertex labels is a dict with an entry for vertex key and the vertex object as a value
        """
        super().__init__()
        self._coords = dict()
        self._vertexLabels = dict()

    def __str__(self):
        """ Return a string representation of the route graph"""
        if self.num_vertices() <= 100 and self.num_edges() <= 100:
            return super().__str__()
        else:
            return "|V| = " + str(self.num_vertices()) + "; |E| = " + str(self.num_edges())
            + "\nNumber of vertices and edges too great to represent"
        
    def get_vertex_by_label(self, element):
        """ Return the vertex that matches element 
        
        Args:
            element - element object
        """
        return self._vertexLabels.get(element)

    def add_vertex(self, element, latitude, longitude):
        """ Return the vertex object added to the route map 
        
        Args:
            element - element object to be added to route map
            latitude - floating point number of element's latitude on a map
            longitude - floating point number of element's longitude on a map
        """
        vertex = super().add_vertex(element)
        self._coords[vertex] = (latitude, longitude)
        self._vertexLabels[element] = vertex
        return vertex
    
    def sp(self, v, w):
        """ Returns the shortest path from v to w 
        
        Args:
            v - the label of the source vertex object
            w - the label of the desination vertex object
        """
        result = self.dijkstra(v)
        finalVert = self.get_vertex_by_label(v)
        curVertex = self.get_vertex_by_label(w)
        path = []
        while curVertex != finalVert:
            prevVertex = result.get(curVertex)
            path.append((curVertex, prevVertex[0]))
            curVertex = prevVertex[1]
        path.append((curVertex, 0))
        path.reverse()
        return path

    def printPath(self, path):
        """ Prints a representation of path
        
        Args:
            path - shortest path list of vertex objects from first vertex to last   
        """
        print("type,latitude,longitude,element,cost")
        for vertex in path:
            coords = self._coords.get(vertex[0])
            print("W", str(coords[0]), str(coords[1]), str(
                vertex[0]._element), str(vertex[1]), sep=",", end="\n")


# Test methods

def graphreader(filename):
    """ Read and return the route map in filename. """
    route = RouteMap()
    file = open(filename, 'r')
    entry = file.readline()  # either 'Node' or 'Edge'
    num = 0
    while entry == 'Node\n':
        num += 1
        nodeid = int(file.readline().split()[1])
        nodeGps = file.readline().split()
        vertex = route.add_vertex(nodeid, round(float(nodeGps[1]), 6), round(float(nodeGps[2]), 6))
        entry = file.readline()  # either 'Node' or 'Edge'
    print('Read', num, 'vertices and added into the graph')
    num = 0
    while entry == 'Edge\n':
        num += 1
        source = int(file.readline().split()[1])
        sv = route.get_vertex_by_label(source)
        target = int(file.readline().split()[1])
        tv = route.get_vertex_by_label(target)
        length = float(file.readline().split()[1])
        time = round(float(file.readline().split()[1]), 2)
        edge = route.add_edge(sv, tv, time)
        file.readline()  # read the one-way data
        entry = file.readline()  # either 'Node' or 'Edge'
    print('Read', num, 'edges and added into the graph')
    print(route)
    return route

def test_dijkstra(graph, src, dest):
    routeMap = graph.dijkstra(src)
    route = graph.sp(src, dest)
    graph.printPath(route)

if __name__ == "__main__":
    graph = graphreader("corkCityData.txt")

    print("Western Gateway to Neptune Stadium")
    test_dijkstra(graph, 1669466540, 1147697924)

    # print("The Old Oak to Cork University Hospital" )
    # test_dijkstra(graph, 358357, 860206013)

    # print("Cork City Gaol to Mahon Point")
    # test_dijkstra(graph, 3777201945, 330068634)
