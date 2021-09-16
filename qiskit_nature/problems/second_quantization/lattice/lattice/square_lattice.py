"""Square lattice"""
from typing import Tuple, Union
from .hyper_cubic import HyperCubic

class SquareLattice(HyperCubic):
    """Square lattice"""
    def __init__(
        self,
        rows:int,
        cols:int,
        edge_parameter:Union[complex, Tuple[complex, complex]] = 1.0, 
        onsite_parameter:complex = 0.0, 
        boundary_condition:Union[str, Tuple[str, str]] = "open"
    ) -> None:
        """
        Args:
            rows: length of the x direction.
            cols: length of the y direction.
            edge_parameter: weights on the unit edges. Defaults to 1.0.
            onsite_parameter: weights on the self loop. Defaults to 0.0.
            boundary_condition: boundary condition for each direction. Defaults to "open".
        """
        self.rows = rows
        self.cols = cols
        super().__init__(
            size=(rows, cols),
            edge_parameter=edge_parameter,
            onsite_parameter=onsite_parameter,
            boundary_condition=boundary_condition
        )
        