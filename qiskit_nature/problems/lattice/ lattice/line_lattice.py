from .hyper_cubic import HyperCubic

class LineLattice(HyperCubic):
    def __init__(
        self,
        num_nodes:int,
        edge_parameter:complex = 1.0,
        onsite_parameter:complex = 0.0,
        boundary_condition:str = "open"
    ) -> "HyperCubic":
        super().__init__(
            size=(num_nodes, ),
            edge_parameter=edge_parameter,
            onsite_parameter=onsite_parameter,
            boundary_condition=boundary_condition
        )