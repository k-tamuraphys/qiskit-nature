from retworkx import PyGraph, NodeIndices, adjacency_matrix, WeightedEdgeList
from retworkx.visualization import mpl_draw
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Tuple, Union

class Lattice:
    # multigraph=False is assumed
    def __init__(self, graph:PyGraph):
        if not graph.multigraph:
            if graph.edges() == [None]*graph.num_edges(): # If weights are None, initialized as 1.0
                weighted_edges = [edge + (1.,) for edge in graph.edge_list()] 
                for start, end, weight in weighted_edges:
                    graph.update_edge(start, end, weight)
            self._graph = graph
        else:
            raise ValueError(f"Invalid `multigraph` {graph.multigraph} is given. `multigraph` must be `False`.")

    @property
    def graph(self) -> PyGraph:
        return self._graph 

    @property
    def num_nodes(self) -> int:
        return self._graph.num_nodes()

    @property
    def nodes(self) -> NodeIndices:
        return self._graph.node_indexes()
    
    @property
    def weighted_edge_list(self) -> WeightedEdgeList:
        return self._graph.weighted_edge_list()
    
    def copy(self) -> "Lattice":
        return Lattice(self._graph.copy())

    
    @classmethod
    def from_adjacency_matrix(cls, adjacency_matrix:np.ndarray) -> "Lattice":
        """returns an instance of Lattice from a given hopping_matrix
        Args:
            adjacency_matrix: adjacency_matrix with real or complex matrix elements

        Returns:
            Lattice made from a given adjacency_matrix
        """
        
        col_length, row_length = adjacency_matrix.shape
        graph = PyGraph(multigraph=False)
        graph.add_nodes_from(range(col_length))
        for source_index in range(col_length):
            for target_index in range(source_index, row_length):
                weight = adjacency_matrix[source_index, target_index]
                if not weight == 0.0:
                    graph.add_edge(source_index, target_index, weight)
                    
        return cls(graph)

    @classmethod
    def from_nodes_edges(cls, num_nodes:int, weighted_edges:List[Tuple[int, int, complex]]) -> "Lattice":
        """returns an instance of Lattice from the number of nodes and the list of edges
        """
        graph = PyGraph(multigraph=False)
        graph.add_nodes_from(range(num_nodes))
        graph.add_edges_from(weighted_edges)
        return cls(graph)

    def to_adjacency_matrix(self) -> np.ndarray:
        """returns the hopping matrix from weighted edges
        the weighted edge list is interpreted as the upper triangular matrix
        """
        real_part = adjacency_matrix(self._graph, weight_fn=lambda x :np.real(x))
        real_part = real_part - (1/2)*np.diag(real_part.diagonal())
        imag_part = adjacency_matrix(self._graph, weight_fn=lambda x :np.imag(x))
        imag_part = np.triu(imag_part) - np.triu(imag_part).T
        return real_part + 1.0j*imag_part

    def draw(self, self_loop:bool=False, pos=None, ax=None, arrows=True, with_labels=False, **kwargs):
        """draws a lattice
        Args:
            self_loop: draw self-loops in a lattice
        """
        if self_loop == True:
            mpl_draw(self._graph, pos, ax, arrows, with_labels, **kwargs)
            plt.draw()
        elif self_loop == False:
            graph_no_loop = self._graph.copy()
            self_loops = [(i, i) for i in range(self.num_nodes) if graph_no_loop.has_edge(i, i)]
            graph_no_loop.remove_edges_from(self_loops)
            mpl_draw(graph_no_loop, pos, ax, arrows, with_labels, **kwargs)
            plt.draw()
        else:
            raise ValueError({f"Invalid value for self_loop is given: {self_loop}"})

class LineLattice(Lattice):
    def __init__(
        self,
        num_nodes:int,
        edge_parameter:complex = 1.0,
        onsite_parameter:complex = 0.0,
        boundary_condition:str = "open"
    ) -> "Lattice":
        
        self.edge_parameter = edge_parameter
        self.onsite_parameter = onsite_parameter
        self.boundary_condition = boundary_condition
        graph = PyGraph(multigraph=False)
        if edge_parameter != 0.0:
            weighted_edge_list = [(i, i+1, edge_parameter) for i in range(num_nodes-1)]
        
        if boundary_condition == "open":
            pass
        elif boundary_condition == "periodic":
            if edge_parameter != 0.0:
                weighted_edge_list.append((num_nodes-1, 0, edge_parameter))
        else:
            raise ValueError(f"Invalid `boundary condition` {boundary_condition} is given. `boundary condition` must be `open` or `periodic`.")
        if onsite_parameter != 0.0:
            onsite_loops = [(i, i, onsite_parameter) for i in range(num_nodes)]
        
        weighted_edge_list = weighted_edge_list + onsite_loops
        graph.add_nodes_from(range(num_nodes))
        graph.add_edges_from(weighted_edge_list)
        super().__init__(graph)

    
    @classmethod
    def from_adjacency_matrix(cls):
        raise NotImplementedError()
        
class SquareLattice(Lattice):
    def __init__(
        self,
        rows:int, # length of x-direction
        cols:int, # length of y-direction
        edge_parameter:Union[complex, Tuple[complex, complex]] = 1.0, #0の時にエッジが無視されない！
        onsite_parameter:complex = 0.0, # 0の時にエッジが無視されない！
        boundary_condition:Union[str, Tuple[str, str]] = "open"
    ) -> "Lattice":

        self.rows = rows
        self.cols = cols

        if isinstance(boundary_condition, str):
            boundary_condition = (boundary_condition, boundary_condition)
        elif isinstance(boundary_condition, tuple):
            pass
        else:
            raise TypeError(f"Type of `boundary_condition` must be str, or tuple of length 2, not {type(edge_parameter)}.")

        self.boundary_conditions = boundary_condition
        
        if isinstance(edge_parameter, (int, float, complex)):
            edge_parameter = (edge_parameter, edge_parameter)
        elif isinstance(edge_parameter, tuple) and len(edge_parameter) == 2:
            pass
        else:
            raise TypeError(f"Type of `edge_parameter` must be int, float, complex, or tuple of length 2, not {type(edge_parameter)}.")

        self.edge_parameter = edge_parameter

        self.onsite_parameter = onsite_parameter

        num_nodes = np.product([rows, cols])
        graph = PyGraph(multigraph=False)
        graph.add_nodes_from(range(num_nodes))

        # x direction
        if edge_parameter[0] != 0.0:
            for x in range(rows-1):
                for y in range(cols):
                    node_a = y*rows + x
                    node_b = node_a + 1
                    graph.add_edge(node_a, node_b, edge_parameter[0])

        
        # y direction
        if edge_parameter[1] != 0.0:
            for x in range(rows):
                for y in range(cols-1):
                    node_a = y*rows + x
                    node_b = node_a + rows
                    graph.add_edge(node_a, node_b, edge_parameter[1])
            
        # self loop
        if onsite_parameter != 0.0:
            for x in range(rows):
                for y in range(cols):
                    node_a = y*rows + x
                    graph.add_edge(node_a, node_a, onsite_parameter)
                

        #boundary condition(x)
        if boundary_condition[0] == "open":
            pass
        elif boundary_condition[0] == "periodic":
            if edge_parameter[0] != 0.0:
                for y in range(cols):
                    node_a = (y+1)*rows - 1
                    node_b = node_a - (rows-1)
                    graph.add_edge(node_a, node_b, edge_parameter[0])
        else:
            raise ValueError(f"Invalid `boundary condition` {boundary_condition[0]} is given. `boundary condition` must be `open` or `periodic`.")

        #boundary condition(y)
        if boundary_condition[1] == "open":
            pass
        elif boundary_condition[1] == "periodic":
            if edge_parameter[1] != 0.0:
                for x in range(rows):
                    node_a = rows*(cols-1) + x
                    node_b = node_a % rows
                    graph.add_edge(node_a, node_b, edge_parameter[1])
        else:
            raise ValueError(f"Invalid `boundary condition` {boundary_condition[1]} is given. `boundary condition` must be `open` or `periodic`.")
                    
        super().__init__(graph)

    @classmethod
    def from_adjacency_matrix(cls):
        raise NotImplementedError()
class TriangularLattice(Lattice):
    def __init__(
        self,
        rows:int,
        cols:int,
        edge_parameter:Union[complex, Tuple[complex, complex, complex]],
        onsite_parameter:complex,
        boundary_condition:str
    ) -> "Lattice":
        self.rows = rows
        self.cols = cols
        self.boundary_condition = boundary_condition

        if isinstance(edge_parameter, (int, float, complex)):
            edge_parameter = (edge_parameter, edge_parameter, edge_parameter)
        elif isinstance(edge_parameter, tuple) and len(edge_parameter) == 3:
            pass
        else:
            raise TypeError(f"Type of `edge parameter` must be int, float, complex, or tuple of length 3, not {type(edge_parameter)}.")
        
        self.edge_parameter = edge_parameter
        self.onsite_parameter = onsite_parameter

        num_nodes = np.product([rows, cols])
        graph = PyGraph(multigraph=False)
        graph.add_nodes_from(range(num_nodes))

        # x direction
        if edge_parameter[0] != 0.0:
            for x in range(rows-1):
                for y in range(cols):
                    node_a = y*rows + x
                    node_b = node_a + 1
                    graph.add_edge(node_a, node_b, edge_parameter[0])

        # y direction
        if edge_parameter[1] != 0.0:
            for x in range(rows):
                for y in range(cols-1):
                    node_a = y*rows + x
                    node_b = node_a + rows
                    graph.add_edge(node_a, node_b, edge_parameter[1])

        # diagonal direction
        if edge_parameter[2] != 0.0:
            for x in range(rows-1):
                for y in range(cols-1):
                    node_a = y*rows + x
                    node_b = node_a + 1 + rows
                    graph.add_edge(node_a, node_b, edge_parameter[2])

        # self loop
        if onsite_parameter != 0.0:
            for x in range(rows):
                for y in range(cols):
                    node_a = y*rows + x
                    graph.add_edge(node_a, node_a, onsite_parameter)


        #boundary condition
        if boundary_condition == "open":
            pass
        elif boundary_condition == "periodic":
            # x direction
            if edge_parameter[0] != 0.0:
                for y in range(cols):
                    node_a = (y+1)*rows - 1
                    node_b = node_a - (rows-1)
                    graph.add_edge(node_a, node_b, edge_parameter[0])
            # y direction
            if edge_parameter[1] != 0.0:
                for x in range(rows):
                    node_a = rows*(cols-1) + x
                    node_b = node_a % rows
                    graph.add_edge(node_a, node_b, edge_parameter[1])
            # diagonal direction
            if edge_parameter[2] != 0.0:
                for y in range(cols-1):
                    node_a = (y+1)*rows - 1
                    node_b = node_a - (rows-1) + rows
                    graph.add_edge(node_a, node_b, edge_parameter[2])

                for x in range(rows-1):
                    node_a = rows*(cols-1) + x
                    node_b = node_a % rows + 1
                    graph.add_edge(node_a, node_b, edge_parameter[2])
            
                node_a = rows*cols - 1
                node_b = 0
                graph.add_edge(node_a, node_b, edge_parameter[2])
            
        else:
            raise ValueError(f"Invalid `boundary condition` {boundary_condition} is given. `boundary condition` must be `open` or `periodic`.")

        super().__init__(graph)

    @classmethod
    def from_adjacency_matrix(cls):
        raise NotImplementedError()
"""
class HexagonalLattice(Lattice):
    # boundary condition -> "open", "periodic", "irregular?"
    pass

class KagomeLattice(Lattice):
    pass

# hyperrectangle lattice in arbitrary d-dimension
class Hyperrectangle(Lattice):
    pass
"""
