from typing import Optional
from qiskit_nature.problems.second_quantization.lattice.lattice.lattice import Lattice
from qiskit_nature.operators.second_quantization import FermionicOp
import numpy as np

class FermiHubbardModel:
    def __init__(self, lattice:Lattice, onsite_interaction:complex):
        """
        Args:
            lattice: lattice geometry on which the model is defined
            onsite_interaction: the strength of the on-site interaction
        """
        self._lattice = lattice 
        self.onsite_interaction = onsite_interaction
    
    @property
    def lattice(self) -> Lattice:
        return (self._lattice).copy() 

    def hopping_matrix(self) -> np.ndarray:
        """returns the hopping matrix

        Returns:
            hopping matrix
        """
        return self.lattice.to_adjacency_matrix()
    @classmethod
    def uniform_parameters(cls, lattice:Lattice, uniform_hopping:complex, uniform_onsite_potential:complex, onsite_interaction:complex) -> "FermiHubbardModel":
        """ set a uniform hopping parameter and on-site potential over a given lattice
        Args:
            lattice: lattice geometry on which the model is defined
            uniform_hopping: hopping parameter 
            uniform_onsite_potential: on-site potential
            onsite_interaction: the strength of the on-site interaction
        """
        graph = lattice.graph
        for node_a, node_b, _ in graph.weighted_edge_list():
            if node_a == node_b:
                pass
            else:
                graph.update_edge(node_a, node_b, uniform_hopping)

        for node_a in graph.node_indexes():
            if graph.has_edge(node_a, node_a):
                graph.update_edge(node_a, node_a, uniform_onsite_potential)
            else:
                graph.add_edge(node_a, node_a, uniform_onsite_potential)

        return cls(Lattice(graph), onsite_interaction)


    def second_q_ops(self, display_format:Optional[str]=None) -> FermionicOp:
        """returns the Hamiltonian of the Fermi-Hubbard model in terms of FermionicOp
        Args:
        display_format: If sparse, the label is represented sparsely during output. if dense, the label is represented densely during output. (default: dense)

        Returns:
            FermionicOp: Hamiltonian of the Fermi-Hubbard model
        """
        kinetic_ham = []
        interaction_ham = []
        weighted_edge_list = self._lattice.weighted_edge_list
        register_length = 2*self._lattice.num_nodes
        # kinetic terms
        for spin in range(2):
            for node_a, node_b, weight in weighted_edge_list: # no duplication in weighted_edge_list is assumed 
                if node_a == node_b:
                    index = 2*node_a + spin
                    kinetic_ham.append((f"N_{index}", weight)) # onsite potential

                else:
                    if node_a < node_b:
                        index_left = 2*node_a + spin
                        index_right = 2*node_b + spin
                        kinetic_ham.append((f"+_{index_left} -_{index_right}", weight))
                        kinetic_ham.append((f"-_{index_left} +_{index_right}", -np.conjugate(weight)))
                    elif node_a > node_b:
                        index_left = 2*node_b + spin
                        index_right = 2*node_a + spin
                        kinetic_ham.append((f"+_{index_left} -_{index_right}", np.conjugate(weight)))
                        kinetic_ham.append((f"-_{index_left} +_{index_right}", -weight))

        # on-site interaction terms
        for node in self._lattice.nodes:
            index_up = 2*node
            index_down = 2*node + 1
            interaction_ham.append((f"N_{index_up} N_{index_down}", self.onsite_interaction))

        ham = kinetic_ham + interaction_ham
        
        return FermionicOp(ham, register_length=register_length, display_format=display_format)

    @classmethod
    def from_parameters(cls, hopping_matrix:np.ndarray, onsite_interaction:float) -> "FermiHubbardModel":
        """returns the Hamiltonian of the Fermi-Hubbard model from the given hopping matrix and on-site interaction.

        Args:
            hopping_matrix (np.ndarray): a real or complex valued square matrix
            onsite_interaction (float): the strength of the on-site interaction

        Returns:
            FermiHubbardModel : Fermi-Hubbard model generated from the given hopping matrix and on-site interaction
        """
        shape = hopping_matrix.shape
        if len(shape) == 2 and shape[0] == shape[1]:
            lat = Lattice.from_adjacency_matrix(hopping_matrix)
            return cls(lat, onsite_interaction)
        else:
            raise ValueError(f"Invalid shape {shape} is given.")
