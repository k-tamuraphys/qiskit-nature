import retworkx
from retworkx.visualization import mpl_draw
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Tuple, Union

class Lattice:
    # edge に重複があってはいけないとする方がわかりやすい?
    def __init__(self, graph:retworkx.PyGraph):
        if graph.edges() == [None]*graph.num_edges(): # weight がない時は 1.0 に初期化
            weighted_edges = [edge + (1.,) for edge in graph.edge_list()] 
            for start, end, weight in weighted_edges:
                graph.update_edge(start, end, weight)
        self._graph = graph

    @property
    def graph(self) -> retworkx.PyGraph:
        return self._graph

    @property
    def num_nodes(self) -> int:
        return self._graph.num_nodes()

    @property
    def nodes(self) -> retworkx.NodeIndices:
        return self._graph.node_indexes()
    
    @property
    def weighted_edge_list(self) -> int:
        return self._graph.weighted_edge_list()

    @classmethod
    def from_adjacency_matrix(cls, adjacency_matrix:np.ndarray) -> "Lattice": 
        """returns an instance of Lattice from a given hopping_matrix
        """
        
        """
        # from_adjacency_matrix を使う場合（正の重みしか載せられない）
        graph  = retworkx.PyGraph.from_adjacency_matrix(adjacency_matrix)
        return cls(graph)
        """
        col_length, row_length = adjacency_matrix.shape
        graph = retworkx.PyGraph()
        graph.add_nodes_from(range(col_length))
        for source_index in range(col_length):
            for target_index in range(source_index, row_length):
                weight = adjacency_matrix[source_index, target_index]
                if not weight == 0.0:
                    graph.add_edge(source_index, target_index, weight)
                    
        return cls(graph)

    @classmethod
    def from_nodes_edges(cls, num_nodes, weighted_edges) -> "Lattice":
        """returns an instance of Lattice from the number of nodes and the list of edges
        """
        graph = retworkx.PyGraph()
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
        hopping_parameter:float,
        onsite_potential:float,
        boundary_condition:str
    ) -> "Lattice":
        
        self.hopping_parameter = hopping_parameter
        self.onsite_potential = onsite_potential
        self.boundary_condition = boundary_condition
        graph = retworkx.PyGraph()
        weighted_edge_list = [(i, i+1, hopping_parameter) for i in range(num_nodes-1)]
        if boundary_condition == "periodic":
            weighted_edge_list.append((num_nodes-1, 0, hopping_parameter))
        elif boundary_condition == "open":
            pass
        else:
            raise ValueError(f"Invalid `boundary condition` {boundary_condition} is given. `boundary condition` must be `open` or `periodic`.")
        onsite_loops = [(i, i, onsite_potential) for i in range(num_nodes)]
        weighted_edge_list = weighted_edge_list + onsite_loops
        graph.add_nodes_from(range(num_nodes))
        graph.add_edges_from(weighted_edge_list)
        super().__init__(graph)
class SquareLattice(Lattice):
    def __init__(
        self,
        size:List[int],
        hopping_parameter:List[float],
        onsite_potential:float,
        boundary_condition:List[str]
    ) -> "Lattice":

        self.size = size
        self.boundary_conditions = boundary_condition
        self.hopping_parameters = hopping_parameter
        self.onsite_potential = onsite_potential

        num_nodes = np.product(size)
        graph = retworkx.PyGraph()
        weighted_edge_list = []

        for x in range(size[0]):
            for y in range(size[1]):
                # x-direction
                node_a = y*size[0] + x
                if x == size[0] - 1:
                    if boundary_condition[0] == "periodic":
                        node_b = y*size[0]
                        weighted_edge_list.append((node_a, node_b, hopping_parameter[0]))
                    elif boundary_condition[0] == "open":
                        pass
                    else:
                        raise ValueError(f"Invalid `boundary condition` {boundary_condition[0]} is given. `boundary condition` must be `open` or `periodic`.")
                else:
                    node_b = node_a + 1
                    weighted_edge_list.append((node_a, node_b, hopping_parameter[0]))
                
                # y-direction
                if y == size[1] - 1:
                    if boundary_condition[1] == "periodic":
                        node_b = x
                        weighted_edge_list.append((node_a, node_b, hopping_parameter[1]))
                    elif boundary_condition[1] == "open":
                        pass
                    else:
                        raise ValueError(f"Invalid `boundary condition` {boundary_condition[1]} is given. `boundary condition` must be `open` or `periodic`.")
                else:
                    node_b = node_a + size[0]
                    weighted_edge_list.append((node_a, node_b, hopping_parameter[1]))

                # on-site potential
                weighted_edge_list.append((node_a, node_a, onsite_potential))

        graph.add_nodes_from(range(num_nodes))
        graph.add_edges_from(weighted_edge_list)
        super().__init__(graph)
class TriangularLattice(Lattice):
    pass

class HexagonalLattice(Lattice):
    pass

class KagomeLattice(Lattice):
    pass


## 任意の並進対象な lattice も作れそう？
class TranslationalSymLattice(Lattice):
    def __init__(self, dims:int, inner_nodes:int, sizes:list, unit_edge_list:list):
        pass
