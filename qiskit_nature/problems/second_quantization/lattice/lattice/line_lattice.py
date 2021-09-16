"""Line lattice"""
from .hyper_cubic import HyperCubic

class LineLattice(HyperCubic):
    """Line lattice"""
    def __init__(
        self,
        num_nodes:int,
        edge_parameter:complex = 1.0,
        onsite_parameter:complex = 0.0,
        boundary_condition:str = "open"
    ) -> None:
        """
        Args:
            num_nodes: the number of sites.
            edge_parameter: weight on the edges. Defaults to 1.0.
            onsite_parameter: weights on the self loop. Defaults to 0.0.
            boundary_condition: boundary condition. Defaults to "open".
        """

        super().__init__(
            size=(num_nodes, ),
            edge_parameter=edge_parameter,
            onsite_parameter=onsite_parameter,
            boundary_condition=boundary_condition
        )
        