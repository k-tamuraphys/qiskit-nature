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

"""Triangular lattice"""
from math import pi
from typing import Callable, List, Optional, Sequence, Tuple, Union
from itertools import product
import numpy as np
from matplotlib.axes import Axes
from matplotlib.colors import Colormap
from retworkx import PyGraph
from .lattice import Lattice


class TriangularLattice(Lattice):
    """Triangular lattice."""

    def __init__(
        self,
        rows: int,
        cols: int,
        edge_parameter: Union[complex, Tuple[complex, complex, complex]] = 1.0,
        onsite_parameter: complex = 0.0,
        boundary_condition: str = "open",
    ) -> None:
        """
        Args:
            rows: Length of the x direction.
            cols: Length of the y direction.
            edge_parameter: Weights on the unit edges.
                This is specified as a tuple of length 3 or a single number.Defaults to 1.0,
            onsite_parameter: Weight on the self-loops, which are edges connecting a node to itself.
                Defaults to 0.0.
            boundary_condition: Boundary condition for each direction.
                Boundary condition must be specified by "open" or "periodic".
                Defaults to "open".

        Raises:
            ValueError: Given size, edge parameter or boundary condition are invalid values.
        """
        self.rows = rows
        self.cols = cols
        self.size = (rows, cols)
        self.dim = 2
        self.boundary_condition = boundary_condition

        if rows < 2 or cols < 2 or (rows, cols) == (2, 2):
            # If it's True, triangular lattice can't be well defined.
            raise ValueError(
                "Both of `rows` and `cols` must not be (2, 2)"
                "and must be greater than or equal to 2."
            )

        if isinstance(edge_parameter, (int, float, complex)):
            edge_parameter = (edge_parameter, edge_parameter, edge_parameter)
        elif isinstance(edge_parameter, tuple):
            if len(edge_parameter) != 3:
                raise ValueError(
                    f"The length of `edge_parameter` must be 3, not {len(edge_parameter)}."
                )

        self.edge_parameter = edge_parameter
        self.onsite_parameter = onsite_parameter

        graph = PyGraph(multigraph=False)
        graph.add_nodes_from(range(np.prod(self.size)))

        # add edges excluding the boundary edges
        coordinates = list(product(*map(range, np.array([rows, cols]))))
        for x, y in coordinates:
            node_a = y * rows + x
            for i in range(3):
                # x direction
                if i == 0 and x != rows - 1:
                    node_b = node_a + 1
                # y direction
                elif i == 1 and y != cols - 1:
                    node_b = node_a + rows
                # diagonal direction
                elif i == 2 and x != rows - 1 and y != cols - 1:
                    node_b = node_a + 1 + rows
                else:
                    continue
                graph.add_edge(node_a, node_b, edge_parameter[i])

        # add self-loops
        for x, y in coordinates:
            node_a = y * rows + x
            graph.add_edge(node_a, node_a, onsite_parameter)

        # depend on boundary condition
        self.boundary_edges = []
        # add edges when the boundary condition is periodic.
        if boundary_condition == "periodic":
            # The periodic boundary condition in the x direction.
            # It makes sense only when rows is greater than 2.
            if rows > 2:
                for y in range(cols):
                    node_a = (y + 1) * rows - 1
                    node_b = node_a - (rows - 1)  # node_b < node_a
                    graph.add_edge(node_b, node_a, edge_parameter[0].conjugate())
                    self.boundary_edges.append((node_b, node_a))
            # The periodic boundary condition in the y direction.
            # It makes sense only when cols is greater than 2.
            if cols > 2:
                for x in range(rows):
                    node_a = rows * (cols - 1) + x
                    node_b = x  # node_b < node_a
                    graph.add_edge(node_b, node_a, edge_parameter[1].conjugate())
                    self.boundary_edges.append((node_b, node_a))
            # The periodic boundary condition in the diagonal direction.
            for y in range(cols - 1):
                node_a = (y + 1) * rows - 1
                node_b = node_a + 1  # node_b > node_a
                graph.add_edge(node_a, node_b, edge_parameter[2])
                self.boundary_edges.append((node_b, node_a))

            for x in range(rows - 1):
                node_a = rows * (cols - 1) + x
                node_b = x + 1  # node_b < node_a
                graph.add_edge(node_b, node_a, edge_parameter[2].conjugate())
                self.boundary_edges.append((node_b, node_a))

            node_a = rows * cols - 1
            node_b = 0  # node_b < node_a
            graph.add_edge(node_b, node_a, edge_parameter[2].conjugate())
            self.boundary_edges.append((node_a, node_b))

        elif boundary_condition != "open":
            raise ValueError(
                f"Invalid `boundary condition` {boundary_condition} is given."
                "`boundary condition` must be `open` or `periodic`."
            )

        # default position
        self.pos = {}
        for index in range(np.prod(self.size)):
            # maps an index to two-dimensional coordinate
            x = index % self.size[0]
            y = index // self.size[0]
            if self.boundary_condition == "open":
                return_x = x
                return_y = y
            elif self.boundary_condition == "periodic":
                # For the periodic boundary conditions,
                # the positionis are shifted so that the edges between boundaries can be seen.
                return_x = x + 0.2 * np.sin(pi * y / (self.size[1] - 1))
                return_y = y + 0.2 * np.sin(pi * x / (self.size[0] - 1))
            self.pos[index] = [return_x, return_y]
        super().__init__(graph)

    # pylint: disable=missing-param-doc
    def draw_without_boundary(
        self,
        self_loop: bool = False,
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
        **kwargs,
    ):
        r"""Draw the lattice.

        Args:
            self_loop: Draw self-loops in the lattice. Defaults to False.
            **kwargs : Kwargs for retworkx.visualization.mpl_draw.
                Please see
                https://qiskit.org/documentation/retworkx/stubs/retworkx.visualization.mpl_draw.html#retworkx.visualization.mpl_draw
                for details.

        Returns:
            A matplotlib figure for the visualization if not running with an
            interactive backend (like in jupyter) or if ``ax`` is not set.
        """
        graph = self.graph

        if pos is None:
            if self.dim == 1:
                pos = {i: [i, 0] for i in range(self.size[0])}
            elif self.dim == 2:
                pos = {i: [i % self.size[0], i // self.size[0]] for i in range(np.prod(self.size))}

        graph.remove_edges_from(self.boundary_edges)

        self._mpl(
            graph=graph,
            self_loop=self_loop,
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
