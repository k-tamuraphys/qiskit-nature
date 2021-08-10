from retworkx import PyGraph, NodeIndices
from retworkx.visualization import mpl_draw
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Tuple, Union
import copy

class Lattice:
    # multigraph=False is assumed
    def __init__(self, graph:PyGraph):
        if not graph.multigraph:
            if graph.edges() == [None]*graph.num_edges(): # If weights are None, initialize as 1.0
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
    def weighted_edge_list(self) -> int:
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
        """
        adjacency_matrix = np.zeros((self.num_nodes, self.num_nodes), dtype=complex)
        for node_a, node_b, weight in self.weighted_edge_list:
            node_a, node_b = tuple(sorted([node_a, node_b]))
            adjacency_matrix[node_a, node_b] = weight

        adjacency_matrix = adjacency_matrix + np.conjugate(np.triu(adjacency_matrix, k=1).T)
        return adjacency_matrix

    def draw(self, pos=None, ax=None, arrows=True, with_labels=False, **kwargs):
        """draws a lattice
        """
        mpl_draw(self._graph, pos, ax, arrows, with_labels, **kwargs)
        plt.draw()

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
        weighted_edge_list = [(i, i+1, edge_parameter) for i in range(num_nodes-1)]
        if boundary_condition == "open":
            pass
        elif boundary_condition == "periodic":
            weighted_edge_list.append((num_nodes-1, 0, edge_parameter))
        else:
            raise ValueError(f"Invalid `boundary condition` {boundary_condition} is given. `boundary condition` must be `open` or `periodic`.")
        onsite_loops = [(i, i, onsite_parameter) for i in range(num_nodes)]
        weighted_edge_list = weighted_edge_list + onsite_loops
        graph.add_nodes_from(range(num_nodes))
        graph.add_edges_from(weighted_edge_list)
        super().__init__(graph)

    
    @classmethod
    def from_adjacency_matrix(cls, adjacency_matrix: np.ndarray) -> "Lattice":
        # 1次元系に適しているかのチェック？
        return Lattice.from_adjacency_matrix(adjacency_matrix)
        
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

        for x in range(rows-1):
            for y in range(cols-1):
                node_a = y*rows + x
                # x-direction
                node_b = node_a + 1
                graph.add_edge(node_a, node_b, edge_parameter[0])
                # y-direction
                node_b = node_a + rows
                graph.add_edge(node_a, node_b, edge_parameter[1])

                # on-site potential
                graph.add_edge(node_a, node_a, onsite_parameter)
        
        # boundary(top)
        for x in range(rows-1):
            node_a = rows*(cols-1) + x
            node_b = node_a + 1
            graph.add_edge(node_a, node_b, edge_parameter[0])
            graph.add_edge(node_a, node_a, onsite_parameter)

        # boundary(right)
        for y in range(cols-1):
            node_a = (y+1)*rows - 1
            node_b = node_a + rows
            graph.add_edge(node_a, node_b, edge_parameter[1])
            graph.add_edge(node_a, node_a, onsite_parameter)
        
        graph.add_edge(rows*cols-1, rows*cols-1, onsite_parameter)
        

        #boundary condition(x)
        if boundary_condition[0] == "open":
            pass
        elif boundary_condition[0] == "periodic":
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
            for x in range(rows):
                node_a = rows*(cols-1) + x
                node_b = node_a % rows
                graph.add_edge(node_a, node_b, edge_parameter[1])
        else:
            raise ValueError(f"Invalid `boundary condition` {boundary_condition[1]} is given. `boundary condition` must be `open` or `periodic`.")
                    
        super().__init__(graph)
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

        for x in range(rows-1):
            for y in range(cols-1):
                node_a = y*rows + x
                node_b = node_a + 1
                graph.add_edge(node_a, node_b, edge_parameter[0])
                node_b = node_a + rows
                graph.add_edge(node_a, node_b, edge_parameter[1])
                node_b = node_a + 1 + rows
                graph.add_edge(node_a, node_b, edge_parameter[2])
                graph.add_edge(node_a, node_a, onsite_parameter)

        # boundary(top)
        for x in range(rows-1):
            node_a = rows*(cols-1) + x
            node_b = node_a + 1
            graph.add_edge(node_a, node_b, edge_parameter[0])
            graph.add_edge(node_a, node_a, onsite_parameter)

        # boundary(right)
        for y in range(cols-1):
            node_a = (y+1)*rows - 1
            node_b = node_a + rows
            graph.add_edge(node_a, node_b, edge_parameter[1])
            graph.add_edge(node_a, node_a, onsite_parameter)

        graph.add_edge(rows*cols-1, rows*cols-1, onsite_parameter)

        #boundary condition
        if boundary_condition == "open":
            pass
        elif boundary_condition == "periodic":
            for y in range(cols-1):
                node_a = (y+1)*rows - 1
                node_b = node_a - (rows-1)
                graph.add_edge(node_a, node_b, edge_parameter[0])
                node_b = node_a - (rows-1) +rows
                graph.add_edge(node_a, node_b, edge_parameter[2])

            for x in range(rows-1):
                node_a = rows*(cols-1) + x
                node_b = node_a % rows
                graph.add_edge(node_a, node_b, edge_parameter[1])
                node_b = node_a % rows + 1
                graph.add_edge(node_a, node_b, edge_parameter[2])
            
            node_a = rows*cols - 1
            node_b = rows*(cols-1)
            graph.add_edge(node_a, node_b, edge_parameter[0])
            node_b = rows - 1
            graph.add_edge(node_a, node_b, edge_parameter[1])
            node_b = 0
            graph.add_edge(node_a, node_b, edge_parameter[2])
            
        else:
            raise ValueError(f"Invalid `boundary condition` {boundary_condition} is given. `boundary condition` must be `open` or `periodic`.")

        super().__init__(graph)

class HexagonalLattice(Lattice):
    # 境界条件?
    pass

class KagomeLattice(Lattice):
    pass


## 任意の並進対象な lattice も作れそう？
## hexiagonal の例などを考えると、 node と edge を remove できる仕組みがあると良い
class TranslationalSymLattice(Lattice):
    def __init__(
        self, 
        size:List[Union[float, complex]],
        num_nodes_in_unitcell:int,
        weighted_edge_in_unitcell:List[Tuple]
    ) -> "Lattice":
        self.dim = len(size)
    
