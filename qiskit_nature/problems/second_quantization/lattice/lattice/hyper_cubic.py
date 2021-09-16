"""The Hyper-Cubic Lattice"""
from typing import Tuple, Union
from itertools import product
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from retworkx import PyGraph
from retworkx.visualization import mpl_draw
from .lattice import Lattice


class HyperCubic(Lattice):
    """Hyper-cubic lattice in d dimension"""
    def __init__(
        self,
        size:Tuple[int, ...],
        edge_parameter: Union[complex, Tuple[complex, ...]] = 1.0,
        onsite_parameter:complex = 0.0,
        boundary_condition:Union[str, Tuple[str, ...]] = "open"
    ) -> None:
        """
        Args:
            size: lengths of each direction.
            edge_parameter: weights on the unit edges. Defaults to 1.0.
            onsite_parameter: weights on the self loop. Defaults to 0.0.
            boundary_condition: boundary condition for each direction. Defaults to "open".

        Raises:
            ValueError: given edge parameter or boundary condition are invalid values.
        """

        self.dim = len(size)

        self.size = size

        # edge parameter
        if isinstance(edge_parameter, (int, float, complex)):
            edge_parameter = (edge_parameter, )*self.dim
        elif isinstance(edge_parameter, tuple):
            if len(edge_parameter) == self.dim:
                pass
            else:
                raise ValueError(f"The length of `edge_parameter` must be the same as that of size, {self.dim}.")
        
        self.edge_parameter = edge_parameter
        # onsite parameter
        self.onsite_parameter = onsite_parameter

        # boundary condition
        if isinstance(boundary_condition, str):
            boundary_condition = (boundary_condition, )*self.dim
        elif isinstance(boundary_condition, tuple):
            if len(boundary_condition) == self.dim:
                pass
            else:
                raise ValueError(f"The length of `boundary_condition` must be the same as that of size, {self.dim}.")
        
        self.boundary_conditions = boundary_condition


        graph = PyGraph(multigraph=False)
        graph.add_nodes_from(range(np.prod(size)))

        coordinates = list(product(*map(range, size)))
        base = np.array([np.prod(size[:i]) for i in range(self.dim)], dtype=int)
        for dim in range(self.dim):
            if edge_parameter[dim] != 0.0:
                for coord in coordinates:
                    if coord[dim] != size[dim] - 1:
                        node_a = np.dot(coord, base)
                        node_b = node_a + base[dim]
                        graph.add_edge(node_a, node_b, edge_parameter[dim])


        # onsite parameter
        if onsite_parameter != 0.0:
            for node_a in range(np.prod(size)):
                graph.add_edge(node_a, node_a, onsite_parameter)

        self.boundary_edges = []
        # boundary condition
        for dim in range(self.dim):
            if boundary_condition[dim] == "open":
                pass
            elif boundary_condition[dim] == "periodic":
                if size[dim] > 2:
                    size_list = list(size)
                    size_list[dim] = 1
                    coordinates = list(product(*map(range, size_list)))
                    for coord in coordinates:
                        node_b = np.dot(coord, base)
                        node_a = node_b + base[dim]*(size[dim]-1)
                        graph.add_edge(node_a, node_b, edge_parameter[dim])
                        self.boundary_edges.append((node_a, node_b))
            else:
                raise ValueError(f"Invalid `boundary condition` {boundary_condition[dim]} is given. `boundary condition` must be `open` or `periodic`.")
        super().__init__(graph)

    @classmethod
    def from_adjacency_matrix(cls, input_adjacency_matrix):
        raise NotImplementedError()

    def draw(
        self, 
        boundary_edges:bool=False, 
        self_loop:bool=False, 
        pos:dict=None, 
        ax:Axes=None, 
        arrows:bool=True, 
        with_labels:bool=False, 
        **kwargs
    ):
        """draws a lattice
        Args:
            boundary_edges: draw edges from the boundaries
            self_loop: draw self-loops in a lattice
        """
        graph = self.graph

        if boundary_edges:
            pass
        elif not boundary_edges:
            graph.remove_edges_from(self.boundary_edges)

        if self_loop:
            mpl_draw(graph, pos, ax, arrows, with_labels, **kwargs)
            plt.draw()
        elif not self_loop:
            self_loops = [(i, i) for i in range(self.num_nodes) if graph.has_edge(i, i)]
            graph.remove_edges_from(self_loops)
            mpl_draw(graph, pos, ax, arrows, with_labels, **kwargs)
            plt.draw()
