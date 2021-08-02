import retworkx
from retworkx.visualization import mpl_draw
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Tuple, Union

class Lattice:
    # edge に重複があってはいけないとする方がわかりやすい
    def __init__(self, graph:retworkx.PyGraph):
        if graph.edges() == [None]*graph.num_edges(): # initialize the weights as 1.0 when graph does not have weights
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
    def from_hopping_matrix(cls, hopping_matrix:np.ndarray) -> "Lattice": # complex matrices are not allowed
        """returns an instance of Lattice from a given hopping_matrix
        """
        graph  = retworkx.PyGraph.from_adjacency_matrix(hopping_matrix)
        return cls(graph)

    @classmethod
    def from_nodes_edges(cls, num_nodes, weighted_edges) -> "Lattice":
        """returns an instance of Lattice from the number of nodes and the list of edges
        """
        graph = retworkx.PyGraph()
        graph.add_nodes_from(range(num_nodes))
        graph.add_edges_from(weighted_edges)
        return cls(graph)

    def to_hopping_matrix(self) -> np.ndarray:
        """returns the hopping matrix from weighted edges

        Returns:
            np.ndarray: hopping matrix
        """
        hopping_matrix = np.zeros((self.num_nodes, self.num_nodes))
        for node_a, node_b, weight in self.weighted_edge_list:
            node_a, node_b = tuple(sorted([node_a, node_b]))
            hopping_matrix[node_a, node_b] = weight

        hopping_matrix = hopping_matrix + np.triu(hopping_matrix, k=1).transpose()
        return hopping_matrix

    def draw(self, pos=None, ax=None, arrows=True, with_labels=False, **kwargs):
        """draws a lattice

        Args:
            pos (dict, optional): [description]. Defaults to None.
            ax (matplotlib.Axes, optional): [description]. Defaults to None.
            arrows (bool, optional): [description]. Defaults to True.
            with_labels (bool, optional): [description]. Defaults to False.
        """
        mpl_draw(self._graph, pos, ax, arrows, with_labels, **kwargs)
        plt.draw()

class Line(Lattice):
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
        self_loops = [(i, i, onsite_potential) for i in range(num_nodes)]
        weighted_edge_list = weighted_edge_list + self_loops
        graph.add_nodes_from(range(num_nodes))
        graph.add_edges_from(weighted_edge_list)
        super().__init__(graph)
class SquareLattice(Lattice):
    def __init__(
        self,
        sizes:List[int],
        hopping_parameter:List[float],
        onsite_potential:float,
        boundary_condition:List[str]
    ) -> "Lattice":

        self.sizes = sizes
        self.boundary_conditions = boundary_condition
        self.hopping_parameters = hopping_parameter,
        self.onsite_potential = onsite_potential,
        num_nodes = np.product(sizes)
        graph = retworkx.PyGraph()
        for x in range(sizes[0]):
            for y in range(sizes[1]):
                pass


class TranslationalSymLattice(Lattice):
    def __init__(self, dims:int, inner_nodes:int, sizes:list, unit_edge_list:list):
        pass
