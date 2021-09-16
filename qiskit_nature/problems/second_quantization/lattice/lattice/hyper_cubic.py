# This code is part of Qiskit.
#
# (C) Copyright IBM 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

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
    """Hyper-cubic lattice in d dimension."""

    def __init__(
        self,
        size: Tuple[int, ...],
        edge_parameter: Union[complex, Tuple[complex, ...]] = 1.0,
        onsite_parameter: complex = 0.0,
        boundary_condition: Union[str, Tuple[str, ...]] = "open",
    ) -> None:
        """
        Args:
            size: Lengths of each direction.
            edge_parameter: Weights on the unit edges. Defaults to 1.0.
            onsite_parameter: Weights on the self loop. Defaults to 0.0.
            boundary_condition: Boundary condition for each direction. Defaults to "open".

        Raises:
            ValueError: Given edge parameter or boundary condition are invalid values.
        """

        self.dim = len(size)

        self.size = size

        # edge parameter
        if isinstance(edge_parameter, (int, float, complex)):
            edge_parameter = (edge_parameter,) * self.dim
        elif isinstance(edge_parameter, tuple):
            if len(edge_parameter) == self.dim:
                pass
            else:
                raise ValueError(
                    f"The length of `edge_parameter` must be the same as that of size, {self.dim}."
                )

        self.edge_parameter = edge_parameter
        # onsite parameter
        self.onsite_parameter = onsite_parameter

        # boundary condition
        if isinstance(boundary_condition, str):
            boundary_condition = (boundary_condition,) * self.dim
        elif isinstance(boundary_condition, tuple):
            if len(boundary_condition) == self.dim:
                pass
            else:
                raise ValueError(
                    f"The length of `boundary_condition` must be the same as that of size, {self.dim}."
                )

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
                        node_a = node_b + base[dim] * (size[dim] - 1)
                        graph.add_edge(node_a, node_b, edge_parameter[dim])
                        self.boundary_edges.append((node_a, node_b))
            else:
                raise ValueError(
                    f"Invalid `boundary condition` {boundary_condition[dim]} is given. `boundary condition` must be `open` or `periodic`."
                )
        super().__init__(graph)

    @classmethod
    def from_adjacency_matrix(cls, input_adjacency_matrix: np.ndarray):
        """Not implemented.

        Args:
            input_adjacency_matrix: Adjacency matrix with real or complex matrix elements.

        Raises:
            NotImplementedError
        """
        raise NotImplementedError()

    def draw(
        self,
        pos: dict = None,
        ax: Axes = None,
        arrows: bool = True,
        with_labels: bool = False,
        self_loop: bool = False,
        boundary_edges: bool = False,
        **kwargs,
    ):
        """Draw the lattice.

        Args:
            pos: An optional dictionary (or
                a :class:`~retworkx.Pos2DMapping` object) with nodes as keys and
                positions as values. If not specified a spring layout positioning will
                be computed. See `layout_functions` for functions that compute
                node positions.
            ax: An optional Matplotlib Axes object to draw the
                graph in.
            arrows: For :class:`~retworkx.PyDiGraph` objects if ``True``
                draw arrowheads. (defaults to ``True``) Note, that the Arrows will
                be the same color as edges.
            arrowstyle (str): An optional string for directed graphs to choose
                the style of the arrowsheads. See
                :class:`matplotlib.patches.ArrowStyle` for more options. By default the
                value is set to ``'-\|>'``.
            arrow_size (int): For directed graphs, choose the size of the arrow
                head's length and width. See
                :class:`matplotlib.patches.FancyArrowPatch` attribute and constructor
                kwarg ``mutation_scale`` for more info. Defaults to 10.
            with_labels: Set to ``True`` to draw labels on the nodes. Edge
                labels will only be drawn if the ``edge_labels`` parameter is set to a
                function. Defaults to ``False``.

            node_list (list): An optional list of node indices in the graph to
                draw. If not specified all nodes will be drawn.
            edge_list (list): An option list of edges in the graph to draw. If not
                specified all edges will be drawn
            node_size (Union[int, list]): Optional size of nodes. If an array is
                specified it must be the same length as node_list. Defaults to 300
            node_color: Optional node color. Can be a single color or
                a sequence of colors with the same length as node_list. Color can be
                string or rgb (or rgba) tuple of floats from 0-1. If numeric values
                are specified they will be mapped to colors using the ``cmap`` and
                ``vmin``,``vmax`` parameters. See :func:`matplotlib.scatter` for more
                details. Defaults to ``'#1f78b4'``)
            node_shape (str): The optional shape node. The specification is the
                same as the :func:`matplotlib.pyplot.scatter` function's ``marker``
                kwarg, valid options are one of
                ``['s', 'o', '^', '>', 'v', '<', 'd', 'p', 'h', '8']``. Defaults to
                ``'o'``
            alpha (float): Optional value for node and edge transparency
            cmap (Colormap): An optional Matplotlib colormap
                object for mapping intensities of nodes
            vmin (float): Optional minimum value for node colormap scaling
            vmax (float): Optional minimum value for node colormap scaling
            linewidths (Union[float, sequence]): An optional line width for symbol
                borders. If a sequence is specified it must be the same length as
                node_list. Defaults to 1.0
            width (Union[float, sequence]): An optional width to use for edges. Can
                either be a float or sequence  of floats. If a sequence is specified
                it must be the same length as node_list. Defaults to 1.0
            edge_color (Union[str, sequence]): color or array of colors (default='k')
                Edge color. Can be a single color or a sequence of colors with the same
                length as edge_list. Color can be string or rgb (or rgba) tuple of
                floats from 0-1. If numeric values are specified they will be
                mapped to colors using the ``edge_cmap`` and ``edge_vmin``,
                ``edge_vmax`` parameters.
            edge_cmap (Colormap): An optional Matplotlib
                colormap for mapping intensities of edges.
            edge_vmin (float): Optional minimum value for edge colormap scaling
            edge_vmax (float): Optional maximum value for node colormap scaling
            style (str): An optional string to specify the edge line style.
                For example, ``'-'``, ``'--'``, ``'-.'``, ``':'`` or words like
                ``'solid'`` or ``'dashed'``. See the
                :class:`matplotlib.patches.FancyArrowPatch` attribute and kwarg
                ``linestyle`` for more details. Defaults to ``'solid'``.
            labels (func): An optional callback function that will be passed a
                node payload and return a string label for the node. For example::

                    labels=str

                could be used to just return a string cast of the node's data payload.
                Or something like::

                    labels=lambda node: node['label']

                could be used if the node payloads are dictionaries.
            edge_labels (func): An optional callback function that will be passed
                an edge payload and return a string label for the edge. For example::

                    edge_labels=str

                could be used to just return a string cast of the edge's data payload.
                Or something like::

                    edge_labels=lambda edge: edge['label']

                could be used if the edge payloads are dictionaries. If this is set
                edge labels will be drawn in the visualization.
            font_size (int): An optional fontsize to use for text labels, By
                default a value of 12 is used for nodes and 10 for edges.
            font_color (str): An optional font color for strings. By default
                ``'k'`` (ie black) is set.
            font_weight (str): An optional string used to specify the font weight.
                By default a value of ``'normal'`` is used.
            font_family (str): An optional font family to use for strings. By
                default ``'sans-serif'`` is used.
            label (str): An optional string label to use for the graph legend.
            connectionstyle (str): An optional value used to create a curved arc
                of rounding radius rad. For example,
                ``connectionstyle='arc3,rad=0.2'``. See
                :class:`matplotlib.patches.ConnectionStyle` and
                :class:`matplotlib.patches.FancyArrowPatch` for more info. By default
                this is set to ``"arc3"``.
            self_loop: Draw self-loops in a lattice. Defaults to False.
            boundary_edges: Draw edges from the boundaries. Defaults to False.

        Returns:
            A matplotlib figure for the visualization if not running with an
            interactive backend (like in jupyter) or if ``ax`` is not set.
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
