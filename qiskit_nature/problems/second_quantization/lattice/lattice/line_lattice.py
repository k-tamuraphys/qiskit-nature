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

"""Line lattice"""
from .hyper_cubic_lattice import HyperCubicLattice


class LineLattice(HyperCubicLattice):
    """Line lattice."""

    def __init__(
        self,
        num_nodes: int,
        edge_parameter: complex = 1.0,
        onsite_parameter: complex = 0.0,
        boundary_condition: str = "open",
    ) -> None:
        """
        Args:
            num_nodes: The number of sites.
            edge_parameter: Weight on the edges. Defaults to 1.0.
            onsite_parameter: Weight on the self loops. Defaults to 0.0.
            boundary_condition: Boundary condition.
                Boundary condition must be specified by "open" or "periodic".
                Defaults to "open".
        """

        super().__init__(
            size=(num_nodes,),
            edge_parameter=edge_parameter,
            onsite_parameter=onsite_parameter,
            boundary_condition=boundary_condition,
        )
