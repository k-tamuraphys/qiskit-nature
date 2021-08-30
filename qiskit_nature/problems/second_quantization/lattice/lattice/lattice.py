from retworkx import PyGraph, NodeIndices, adjacency_matrix, WeightedEdgeList
from retworkx.visualization import mpl_draw
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Tuple

class Lattice:
    # multigraph=False is assumed
    def __init__(self, graph:PyGraph) -> None:
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
        return self._graph.copy()

    @property
    def num_nodes(self) -> int:
        return self.graph.num_nodes()

    @property
    def nodes(self) -> NodeIndices:
        return self.graph.node_indexes()
    
    @property
    def weighted_edge_list(self) -> WeightedEdgeList:
        return self.graph.weighted_edge_list()
    
    def copy(self) -> "Lattice":
        return Lattice(self.graph.copy())

    
    @classmethod
    def from_adjacency_matrix(cls, adjacency_matrix:np.ndarray) -> "Lattice":
        """returns an instance of Lattice from a given hopping_matrix
        Args:
            adjacency_matrix: adjacency_matrix with real or complex matrix elements

        Returns:
            Lattice generated from a given adjacency_matrix
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
        real_part = adjacency_matrix(self.graph, weight_fn=lambda x :np.real(x))
        real_part = real_part - (1/2)*np.diag(real_part.diagonal())
        imag_part = adjacency_matrix(self.graph, weight_fn=lambda x :np.imag(x))
        imag_part = np.triu(imag_part) - np.triu(imag_part).T
        return real_part + 1.0j*imag_part

    def draw(self, self_loop:bool=False, pos=None, ax=None, arrows=True, with_labels=False, **kwargs):
        """draws a lattice
        Args:
            self_loop: draw self-loops in a lattice
        """
        if self_loop == True:
            mpl_draw(self.graph, pos, ax, arrows, with_labels, **kwargs)
            plt.draw()
        elif self_loop == False:
            graph_no_loop = self.graph
            self_loops = [(i, i) for i in range(self.num_nodes) if graph_no_loop.has_edge(i, i)]
            graph_no_loop.remove_edges_from(self_loops)
            mpl_draw(graph_no_loop, pos, ax, arrows, with_labels, **kwargs)
            plt.draw()