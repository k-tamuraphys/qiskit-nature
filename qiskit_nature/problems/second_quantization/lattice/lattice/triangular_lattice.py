from retworkx import PyGraph
import numpy as np
from typing import Tuple, Union
from .lattice import Lattice
from retworkx.visualization import mpl_draw
import matplotlib.pyplot as plt
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
        self.size = (rows, cols)
        self.dim = 2
        self.boundary_condition = boundary_condition

        if isinstance(edge_parameter, (int, float, complex)):
            edge_parameter = (edge_parameter, edge_parameter, edge_parameter)
        elif isinstance(edge_parameter, tuple) and len(edge_parameter) == 3:
            pass
        else:
            raise TypeError(f"Type of `edge parameter` must be int, float, complex, or tuple of length 3, not {type(edge_parameter)}.")
        
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