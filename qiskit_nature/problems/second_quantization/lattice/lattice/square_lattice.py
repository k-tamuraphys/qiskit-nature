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

"""Square lattice"""
from typing import Tuple, Union
import numpy as np
from .hyper_cubic import HyperCubic


class SquareLattice(HyperCubic):
    """Square lattice."""

    def __init__(
        self,
        rows: int,
        cols: int,
        edge_parameter: Union[complex, Tuple[complex, complex]] = 1.0,
        onsite_parameter: complex = 0.0,
        boundary_condition: Union[str, Tuple[str, str]] = "open",
    ) -> None:
        """
        Args:
            rows: Length of the x direction.
            cols: Length of the y direction.
            edge_parameter: Weights on the unit edges. Defaults to 1.0.
            onsite_parameter: Weights on the self loop. Defaults to 0.0.
            boundary_condition: Boundary condition for each direction. Defaults to "open".
        """
        self.rows = rows
        self.cols = cols
        super().__init__(
            size=(rows, cols),
            edge_parameter=edge_parameter,
            onsite_parameter=onsite_parameter,
            boundary_condition=boundary_condition,
        )

    @classmethod
    def from_adjacency_matrix(cls, input_adjacency_matrix: np.ndarray):
        """Not implemented.

        Args:
            input_adjacency_matrix: Adjacency matrix with real or complex matrix elements.

        Raises:
            NotImplementedError
        """
        raise NotImplementedError()
