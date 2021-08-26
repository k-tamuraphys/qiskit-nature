from retworkx import PyGraph
import numpy as np
from typing import Tuple, Union
from itertools import product
from .lattice import Lattice
from retworkx.visualization import mpl_draw
import matplotlib.pyplot as plt


class HyperCubic(Lattice):
    def __init__(
        self,
        size:Tuple[int, ...],
        edge_parameter: Union[complex, Tuple[complex, ...]] = 1.0,
        onsite_parameter:complex = 0.0,
        boundary_condition:Union[str, Tuple[str, ...]] = "open"
    ) -> "Lattice":

        self.dim = len(size)

        self.size = size

        # edge parameter
        if isinstance(edge_parameter, (int, float, complex)):
            edge_parameter = (edge_parameter, )*self.dim
        elif isinstance(edge_parameter, tuple) and len(edge_parameter) == self.dim:
            pass
        else:
            raise TypeError(f"Type of `edge_parameter` must be int, float, complex, or tuple of length {self.dim}, not {type(edge_parameter)}.")
        self.edge_parameter = edge_parameter
        # onsite parameter
        self.onsite_parameter = onsite_parameter

        # boundary condition
        if isinstance(boundary_condition, str):
            boundary_condition = (boundary_condition, )*self.dim
        elif isinstance(boundary_condition, tuple) and len(boundary_condition) == self.dim:
            pass
        else:
            raise TypeError(f"Type of `boundary_condition` must be str, or tuple of length {self.dim}, not {type(boundary_condition)}.")
        self.boundary_conditions = boundary_condition


        graph = PyGraph(multigraph=False)
        graph.add_nodes_from(range(np.prod(size)))

        coordinates = list(product(*map(range, size)))
        base = np.array([np.prod(size[:i]) for i in range(self.dim)], dtype=int)
        for d in range(self.dim):
            if edge_parameter[d] != 0.0:
                for coord in coordinates:
                    if coord[d] != size[d] - 1:
                        node_a = np.dot(coord, base)
                        node_b = node_a + base[d]
                        graph.add_edge(node_a, node_b, edge_parameter[d])


        # onsite parameter
        if onsite_parameter != 0.0:
            for node_a in range(np.prod(size)):
                graph.add_edge(node_a, node_a, onsite_parameter)

        self.boundary_edges = []
        # boundary condition
        for d in range(self.dim):
            if boundary_condition[d] == "open":
                pass
            elif boundary_condition[d] == "periodic":
                if size[d] > 2: # open and periodic boundary condition for a lattice with two sites are the same
                    size_list = list(size)
                    size_list[d] = 1
                    coordinates = list(product(*map(range, size_list)))
                    for coord in coordinates:
                        node_b = np.dot(coord, base)
                        node_a = node_b + base[d]*(size[d]-1)
                        graph.add_edge(node_a, node_b, edge_parameter[d])
                        self.boundary_edges.append((node_a, node_b))
            else:
                raise ValueError(f"Invalid `boundary condition` {boundary_condition[d]} is given. `boundary condition` must be `open` or `periodic`.")
        super().__init__(graph)

    @classmethod
    def from_adjacency_matrix(cls):
        raise NotImplementedError()

    def draw(self, boundary_edges:bool=False, self_loop:bool=False, pos=None, ax=None, arrows=True, with_labels=False, **kwargs):
        """draws a lattice
        Args:
            boundary_edges: draw edges from the boundaries
            self_loop: draw self-loops in a lattice
        """
        graph = self.graph

        if boundary_edges == True:
            pass
        elif boundary_edges == False:
            graph.remove_edges_from(self.boundary_edges)
        else:
            raise TypeError(f"Type of `boundary_edges` must be bool, not {type(boundary_edges)}.")

        if self_loop == True:
            mpl_draw(graph, pos, ax, arrows, with_labels, **kwargs)
            plt.draw()
        elif self_loop == False:
            self_loops = [(i, i) for i in range(self.num_nodes) if graph.has_edge(i, i)]
            graph.remove_edges_from(self_loops)
            mpl_draw(graph, pos, ax, arrows, with_labels, **kwargs)
            plt.draw()
        else:
            raise TypeError(f"Type of `self_loop` must be bool, not {type(self_loop)}.")