from .hyper_cubic import HyperCubic
from typing import Tuple, Union

class SquareLattice(HyperCubic):
    def __init__(
        self,
        rows:int, # length of x-direction
        cols:int, # length of y-direction
        edge_parameter:Union[complex, Tuple[complex, complex]] = 1.0, 
        onsite_parameter:complex = 0.0, 
        boundary_condition:Union[str, Tuple[str, str]] = "open"
    ) -> "HyperCubic":
        self.rows = rows
        self.cols = cols
        super().__init__(
            size=(rows, cols),
            edge_parameter=edge_parameter,
            onsite_parameter=onsite_parameter,
            boundary_condition=boundary_condition
        )