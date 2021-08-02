import lattice
from qiskit_nature.operators.second_quantization import FermionicOp

def Fermi_Hubbard(lattice:lattice.Lattice, onsite_int:float) -> FermionicOp:
    # ホッピング行列もオプショナルな引数として持たせていいかもしれない
    # もしホッピング行列が与えられないのならばlattice側が持ってる重みをホッピング行列として扱う
    # もしホッピング行列が与えられれば、それをホッピング行列として扱う（どういう形式のinputが適切？）例：[(0, 1, -1.0), (1, 1, 1.0), ...]
    """returns the Hamiltonian of the Fermi-Hubbard model in terms of FermionicOp

    Args:
        lattice (lattice.Lattice): lattice with weighted edges
        onsite_int (float): parameter of the on-site interaction

    Returns:
        FermionicOp: Hamiltonian of the Fermi-Hubbard model
    """
    kinetic_Ham = []
    interaction_Ham = []
    weighted_edge_list = lattice.weighted_edge_list
    register_length = 2*lattice.num_nodes
    # kinetic terms
    for spin in range(2):
        for node_a, node_b, weight in weighted_edge_list: # edge に重複があってはいけない
            if node_a == node_b:
                index = 2*node_a + spin
                kinetic_Ham.append((f"N_{index}", weight))

            else:
                if node_a < node_b:
                    index_left = 2*node_a + spin
                    index_right = 2*node_b + spin
                    kinetic_Ham.append((f"+_{index_left} -_{index_right}", weight))
                    kinetic_Ham.append((f"-_{index_left} +_{index_right}", -weight))
                    #when weight is a complex value
                    #kinetic_Ham.append((f"+_{index_left} -_{index_right}", weight))
                    #kinetic_Ham.append((f"-_{index_left} +_{index_right}", -np.conjugate(weight)))
                    


                elif node_a > node_b:
                    index_left = 2*node_b + spin
                    index_right = 2*node_a + spin
                    kinetic_Ham.append((f"+_{index_left} -_{index_right}", weight))
                    kinetic_Ham.append((f"-_{index_left} +_{index_right}", -weight))
                    #when weight is a complex value
                    # kinetic_Ham.append((f"+_{index_left} -_{index_right}", np.conjugate(weight)))
                    # kinetic_Ham.append((f"-_{index_left} +_{index_right}", -weight))

                
    # on-site interaction terms
    for node in lattice.nodes:
        index_up = 2*node
        index_down = 2*node + 1
        interaction_Ham.append((f"N_{index_up} N_{index_down}", onsite_int))

    Ham = kinetic_Ham + interaction_Ham
    return FermionicOp(Ham, register_length=register_length)