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
from math import pi
from typing import List, Tuple, Union, Optional, Sequence, Callable
from itertools import product
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import Colormap
from retworkx import PyGraph
from retworkx.visualization import mpl_draw
from .lattice import Lattice


class HyperCubicLattice(Lattice):
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
                        if node_a < node_b:
                            graph.add_edge(node_a, node_b, edge_parameter[dim])
                            self.boundary_edges.append((node_a, node_b))
                        elif node_a > node_b:
                            graph.add_edge(node_b, node_a, np.conjugate(edge_parameter[dim]))
                            self.boundary_edges.append((node_a, node_b))
            else:
                raise ValueError(
                    f"Invalid `boundary condition` {boundary_condition[dim]} is given."
                    "`boundary condition` must be `open` or `periodic`."
                )

        if self.dim == 1:
            if self.boundary_conditions == ("open",):
                self.position = {i: [i, 0] for i in range(self.size[0])}
            elif self.boundary_conditions == ("periodic",):
                theta = 2 * pi / self.size[0]
                self.position = {
                    i: [np.cos(i * theta), np.sin(i * theta)] for i in range(self.size[0])
                }
        elif self.dim == 2:
            position_dict = {}
            for index in range(np.prod(size)):
                x = index % self.size[0]
                y = index // self.size[0]
                if self.boundary_conditions[1] == "open":
                    return_x = x
                elif self.boundary_conditions[1] == "periodic":
                    return_x = x + 0.2 * np.sin(pi * y / (self.size[1] - 1))
                if self.boundary_conditions[0] == "open":
                    return_y = y
                elif self.boundary_conditions[0] == "periodic":
                    return_y = y + 0.2 * np.sin(pi * x / (self.size[0] - 1))
                position_dict[index] = [return_x, return_y]
            self.position = position_dict
        else:
            self.position = None

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

    # pylint: disable=arguments-differ
    def draw(
        self,
        pos: Optional[dict] = None,
        ax: Optional[Axes] = None,
        arrows: bool = True,
        arrowstyle: Optional[str] = None,
        arrow_size: int = 10,
        with_labels: bool = False,
        node_list: Optional[list] = None,
        edge_list: Optional[list] = None,
        node_size: Union[int, list] = 300,
        node_color: Union[
            str,
            Tuple[float, float, float],
            Tuple[float, float, float],
            List[Union[str, Tuple[float, float, float], Tuple[float, float, float, float]]],
        ] = "#1f78b4",
        node_shape: str = "o",
        alpha: Optional[float] = None,
        cmap: Optional[Colormap] = None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        linewidths: Union[float, Sequence] = 1.0,
        width: Union[float, Sequence] = 1.0,
        edge_color: Union[str, Sequence] = "k",
        edge_cmap: Optional[Colormap] = None,
        edge_vmin: Optional[float] = None,
        edge_vmax: Optional[float] = None,
        style: str = "solid",
        labels: Optional[Callable] = None,
        edge_labels: Optional[Callable] = None,
        font_size: int = 12,
        font_color: str = "k",
        font_weight: str = "normal",
        font_family: str = "sans-serif",
        label: Optional[str] = None,
        connectionstyle: str = "arc3",
        self_loop: bool = False,
        boundary_edges: bool = False,
        **kwargs,
    ):
        r"""Draw the lattice.

        Args:
            pos: An optional dictionary (or
                a :class:`~retworkx.Pos2DMapping` object) with nodes as keys and
                positions as values. If not specified a spring layout positioning will
                be computed. See `layout_functions` for functions that compute
                node positions.
            ax: An optional Axes object to draw the
                graph in.
            arrows: For :class:`~retworkx.PyDiGraph` objects if ``True``
                draw arrowheads. (defaults to ``True``) Note, that the Arrows will
                be the same color as edges.
            arrowstyle: An optional string for directed graphs to choose
                the style of the arrowheads. See
                :class:`matplotlib.patches.ArrowStyle` for more options. By default the
                value is set to ``'-\|>'``.
            arrow_size: For directed graphs, choose the size of the arrow
                head's length and width. See
                :class:`matplotlib.patches.FancyArrowPatch` attribute and constructor
                kwarg ``mutation_scale`` for more info. Defaults to 10.
            with_labels: Set to ``True`` to draw labels on the nodes. Edge
                labels will only be drawn if the ``edge_labels`` parameter is set to a
                function. Defaults to ``False``.

            node_list: An optional list of node indices in the graph to
                draw. If not specified all nodes will be drawn.
            edge_list: An option list of edges in the graph to draw. If not
                specified all edges will be drawn
            node_size: Optional size of nodes. If an array is
                specified it must be the same length as node_list. Defaults to 300
            node_color: Optional node color. Can be a single color or
                a sequence of colors with the same length as node_list. Color can be
                string or rgb (or rgba) tuple of floats from 0-1. If numeric values
                are specified they will be mapped to colors using the ``cmap`` and
                ``vmin``,``vmax`` parameters. See :func:`matplotlib.scatter` for more
                details. Defaults to ``'#1f78b4'``)
            node_shape: The optional shape node. The specification is the
                same as the :func:`matplotlib.pyplot.scatter` function's ``marker``
                kwarg, valid options are one of
                ``['s', 'o', '^', '>', 'v', '<', 'd', 'p', 'h', '8']``. Defaults to
                ``'o'``
            alpha: Optional value for node and edge transparency
            cmap: An optional Colormap
                object for mapping intensities of nodes
            vmin: Optional minimum value for node colormap scaling
            vmax: Optional minimum value for node colormap scaling
            linewidths: An optional line width for symbol
                borders. If a sequence is specified it must be the same length as
                node_list. Defaults to 1.0
            width: An optional width to use for edges. Can
                either be a float or sequence  of floats. If a sequence is specified
                it must be the same length as node_list. Defaults to 1.0
            edge_color: color or array of colors (default='k')
                Edge color. Can be a single color or a sequence of colors with the same
                length as edge_list. Color can be string or rgb (or rgba) tuple of
                floats from 0-1. If numeric values are specified they will be
                mapped to colors using the ``edge_cmap`` and ``edge_vmin``,
                ``edge_vmax`` parameters.
            edge_cmap: An optional Matplotlib
                colormap for mapping intensities of edges.
            edge_vmin: Optional minimum value for edge colormap scaling
            edge_vmax: Optional maximum value for node colormap scaling
            style: An optional string to specify the edge line style.
                For example, ``'-'``, ``'--'``, ``'-.'``, ``':'`` or words like
                ``'solid'`` or ``'dashed'``. See the
                :class:`matplotlib.patches.FancyArrowPatch` attribute and kwarg
                ``linestyle`` for more details. Defaults to ``'solid'``.
            labels: An optional callback function that will be passed a
                node payload and return a string label for the node. For example::

                    labels=str

                could be used to just return a string cast of the node's data payload.
                Or something like::

                    labels=lambda node: node['label']

                could be used if the node payloads are dictionaries.
            edge_labels: An optional callback function that will be passed
                an edge payload and return a string label for the edge. For example::

                    edge_labels=str

                could be used to just return a string cast of the edge's data payload.
                Or something like::

                    edge_labels=lambda edge: edge['label']

                could be used if the edge payloads are dictionaries. If this is set
                edge labels will be drawn in the visualization.
            font_size: An optional font size to use for text labels, By
                default a value of 12 is used for nodes and 10 for edges.
            font_color: An optional font color for strings. By default
                ``'k'`` (i.e. black) is set.
            font_weight: An optional string used to specify the font weight.
                By default a value of ``'normal'`` is used.
            font_family: An optional font family to use for strings. By
                default ``'sans-serif'`` is used.
            label: An optional string label to use for the graph legend.
            connectionstyle: An optional value used to create a curved arc
                of rounding radius rad. For example,
                ``connectionstyle='arc3,rad=0.2'``. See
                :class:`matplotlib.patches.ConnectionStyle` and
                :class:`matplotlib.patches.FancyArrowPatch` for more info. By default
                this is set to ``"arc3"``.
            self_loop: Draw self-loops in a lattice. Defaults to False.
            boundary_edges: Draw edges from the boundaries. Defaults to False.
            **kwargs: Kwargs for drawing the lattice.

        Returns:
            A matplotlib figure for the visualization if not running with an
            interactive backend (like in jupyter) or if ``ax`` is not set.
        """
        graph = self.graph

        if pos is None:
            pos = self.position

        if not boundary_edges:
            graph.remove_edges_from(self.boundary_edges)

        if not self_loop:
            self_loops = [(i, i) for i in range(self.num_nodes) if graph.has_edge(i, i)]
            graph.remove_edges_from(self_loops)

        mpl_draw(
            graph=graph,
            pos=pos,
            ax=ax,
            arrows=arrows,
            arrowstyle=arrowstyle,
            arrow_size=arrow_size,
            with_labels=with_labels,
            node_list=node_list,
            edge_list=edge_list,
            node_size=node_size,
            node_color=node_color,
            node_shape=node_shape,
            alpha=alpha,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            linewidths=linewidths,
            width=width,
            edge_color=edge_color,
            edge_cmap=edge_cmap,
            edge_vmin=edge_vmin,
            edge_vmax=edge_vmax,
            style=style,
            labels=labels,
            edge_labels=edge_labels,
            font_size=font_size,
            font_color=font_color,
            font_weight=font_weight,
            font_family=font_family,
            label=label,
            connectionstyle=connectionstyle,
            **kwargs,
        )
        plt.draw()
