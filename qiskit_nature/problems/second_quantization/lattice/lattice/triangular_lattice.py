"""Triangular lattice"""
from typing import Tuple, Union
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from retworkx import PyGraph
from retworkx.visualization import mpl_draw
from .lattice import Lattice
class TriangularLattice(Lattice):
    """Triangular lattice"""
    def __init__(
        self,
        rows:int,
        cols:int,
        edge_parameter:Union[complex, Tuple[complex, complex, complex]] = 1.0,
        onsite_parameter:complex = 0.0,
        boundary_condition:str = "open"
    ) -> None:
        """
        Args:
            rows: length of the x direction.
            cols: length of the y direction.
            edge_parameter: weights on the unit edges. Defaults to 1.0,
            onsite_parameter: weights on the self loop. Defaults to 0.0.
            boundary_condition: boundary condition for each direction. Defaults to "open".

        Raises:
            ValueError: given edge parameter or boundary condition are invalid values.
        """
        self.rows = rows
        self.cols = cols
        self.size = (rows, cols)
        self.dim = 2
        self.boundary_condition = boundary_condition

        if isinstance(edge_parameter, (int, float, complex)):
            edge_parameter = (edge_parameter, edge_parameter, edge_parameter)
        elif isinstance(edge_parameter, tuple):
            if len(edge_parameter) == 3:
                pass
            else:
                ValueError(f"The length of `edge_parameter` must be 3, not {len(edge_parameter)}.")
        
        self.edge_parameter = edge_parameter
        self.onsite_parameter = onsite_parameter

        graph = PyGraph(multigraph=False)
        graph.add_nodes_from(range(np.prod(self.size)))

        # edge in the x direction
        if edge_parameter[0] != 0.0:
            for x in range(rows-1):
                for y in range(cols):
                    node_a = y*rows + x
                    node_b = node_a + 1
                    graph.add_edge(node_a, node_b, edge_parameter[0])

        # edge in the y direction
        if edge_parameter[1] != 0.0:
            for x in range(rows):
                for y in range(cols-1):
                    node_a = y*rows + x
                    node_b = node_a + rows
                    graph.add_edge(node_a, node_b, edge_parameter[1])

        # edge in the diagonal direction
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

        self.boundary_edges = []
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
                    self.boundary_edges.append((node_a, node_b))
            # y direction
            if edge_parameter[1] != 0.0:
                for x in range(rows):
                    node_a = rows*(cols-1) + x
                    node_b = node_a % rows
                    graph.add_edge(node_a, node_b, edge_parameter[1])
                    self.boundary_edges.append((node_a, node_b))
            # diagonal direction
            if edge_parameter[2] != 0.0:
                for y in range(cols-1):
                    node_a = (y+1)*rows - 1
                    node_b = node_a - (rows-1) + rows
                    graph.add_edge(node_a, node_b, edge_parameter[2])
                    self.boundary_edges.append((node_a, node_b))

                for x in range(rows-1):
                    node_a = rows*(cols-1) + x
                    node_b = node_a % rows + 1
                    graph.add_edge(node_a, node_b, edge_parameter[2])
                    self.boundary_edges.append((node_a, node_b))
            
                node_a = rows*cols - 1
                node_b = 0
                graph.add_edge(node_a, node_b, edge_parameter[2])
                self.boundary_edges.append((node_a, node_b))
            
        else:
            raise ValueError(f"Invalid `boundary condition` {boundary_condition} is given. `boundary condition` must be `open` or `periodic`.")
        super().__init__(graph)

    @classmethod
    def from_adjacency_matrix(cls):
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
