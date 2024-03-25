from __future__ import annotations

import numpy as np
from typing import Union
import pickle
from itertools import combinations


class Pauli:
    """
    Class to represent a Pauli operator.
    """

    def __init__(self, T: Union[list, np.ndarray]):
        """
        Initializes the Pauli operator.

        Args:
            T: The Pauli operator.
        """
        if isinstance(T, list):
            T = np.array(T)
        elif not isinstance(T, np.ndarray):
            raise ValueError("T should be either a list or a numpy array.")

        self.T = T
        self.n = T.shape[0] // 2

    def __str__(self) -> str:
        pauli_str = ""
        for i in range(self.n):
            if self.T[i] == 1 and self.T[i + self.n] == 1:
                pauli_str += "Y"
            elif self.T[i + self.n] == 1:
                pauli_str += "Z"
            elif self.T[i] == 1:
                pauli_str += "X"
            else:
                pauli_str += "I"

        return pauli_str
    
    def __add__(self, other) -> Pauli:
        T  = (self.T + other.T) % 2
        return Pauli(T)

    def __repr__(self) -> str:
        return self.__str__()

    def __hash__(self) -> int:
        return hash(tuple(self.T))

    def __eq__(self, other) -> bool:
        return np.array_equal(self.T, other.T)


def is_commute(pauli_1: Pauli, pauli_2: Pauli) -> bool:
    """
    Checks if the two given Pauli's commute.

    Args:
        pauli_1: The first Pauli operator.
        pauli_2: The second Pauli operator.

    Returns:
        bool: True if the two Pauli's commute, False otherwise.
    """
    if pauli_1.n != pauli_2.n:
        raise ValueError(
            f"The two Pauli's should be in the same dimension. Pauli 1: {pauli_1.n}, Pauli 2: {pauli_2.n}"
        )

    n = pauli_1.n

    T_a = pauli_1.T
    T_b = pauli_2.T

    T_ax = T_a[:n]
    T_az = T_a[n:]

    T_bx = T_b[:n]
    T_bz = T_b[n:]

    return (T_az.dot(T_bx) + T_ax.dot(T_bz)) % 2 == 0


def find_commutative_paulis(pauli: Pauli) -> set[Pauli]:
    """
    Finds all the Pauli's that commute with the given Pauli.

    Args:
        pauli: The given Pauli operator.

    Returns:
        set[Pauli]: A set of Pauli operators that commute with the given Pauli.
    """

    n = pauli.n * 2
    commutative_paulis = set()

    for i in range(2**n):
        T_b = np.array([int(x) for x in list(np.base_repr(i, base=2).zfill(n))])
        pauli_b = Pauli(T_b)
        if is_commute(pauli, pauli_b):
            commutative_paulis.add(pauli_b)

    return commutative_paulis


def find_all_commutative_paulis(n: int) -> dict[Pauli, set[Pauli]]:
    """
    Finds commutative Pauli's for all Pauli's and stores them in a dictionary.

    Args:
        n: The number of qubits.

    Returns:
        dict[Pauli, set[Pauli]]: A dictionary containing the Pauli's as keys and their commutative Pauli's as values.
    """

    # Pauli string for n qubit is of length 2n
    n = 2 * n

    commutative_paulis = {}
    for i in range(2**n):
        T_a = np.array([int(x) for x in list(np.base_repr(i, base=2).zfill(n))])
        pauli = Pauli(T_a)
        commutative_paulis[pauli] = find_commutative_paulis(pauli)

    return commutative_paulis


def save_commutative_paulis(n: int) -> None:
    """
    Calculates and saves the commutative Pauli's for all Pauli's in a file.

    Args:
        n: The number of qubits.
    """
    commutative_paulis = find_all_commutative_paulis(n)

    with open(f"commutative_paulis_{n}.pkl", "wb") as f:
        pickle.dump(commutative_paulis, f)


def load_commutative_paulis(n: int) -> dict[tuple, np.ndarray]:
    """
    Loads the commutative Pauli's from a file.

    Args:
        n: The number of qubits.

    Returns:
        dict[tuple, np.ndarray]: A dictionary containing the Pauli's as keys and their commutative Pauli's as values.
    """

    with open(f"commutative_paulis_{n}.pkl", "rb") as f:
        commutative_paulis = pickle.load(f)

    return commutative_paulis


def generate_canonical_isotropic_gens(n: int, m: int) -> list[Pauli]:
    """
    Generates the canonical isotropic generators in n qubits, that has n - m dimensions.

    Args:
        n: The number of qubits.
        m: Parameter that determines the dimension of the isotropic generators. Dimension is n - m.

    Returns:
        list[Pauli]: A list of Pauli operators that are the canonical isotropic generators.
    """

    if m > n:
        raise ValueError(
            "The dimension of the isotropic generators should be less than or equal to the number of qubits."
        )

    if m == n:
        return [Pauli(np.zeros(2 * n))]

    generators = []

    for i in range(n - m):
        x_pos = n - i - 1

        # Creates a Pauli where I * I * ... * X * I * ... * I where X is at the x_pos. Interpret * as tensor product.
        T = np.zeros(2 * n)
        T[x_pos] = 1

        pauli = Pauli(T)

        generators.append(pauli)

    return generators


def find_commutative_paulis_for_isotropic(
    isotropic_gens: list[Pauli], commutative_paulis: dict[Pauli, set[Pauli]] = None
) -> set[Pauli]:
    """
    Finds the commutative Pauli's for the isotropic generators.

    Args:
        isotropic_gens: The isotropic generators.
        commutative_paulis: The commutative Pauli's for all Pauli's.

    Returns:
        set[Pauli]: A set of Pauli operators that commute with all isotropic generators. Hence commute with the whole group.
    """

    if commutative_paulis is None:
        try:
            commutative_paulis = load_commutative_paulis(isotropic_gens[0].n)
        except FileNotFoundError:
            save_commutative_paulis(isotropic_gens[0].n)
            commutative_paulis = load_commutative_paulis(isotropic_gens[0].n)

    result = commutative_paulis[isotropic_gens[0]]

    for i in range(1, len(isotropic_gens)):
        result = result.intersection(commutative_paulis[isotropic_gens[i]])

    return result


def generate_isotropic_from_gens(isotropic_gens: list[Pauli]) -> set[Pauli]:
    """
    Generates the isotropic group from the isotropic generators.

    Args:
        isotropic_gens: The isotropic generators.

    Returns:
        set[Pauli]: The isotropic group generated by the isotropic generators.
    """
    iden = np.zeros(2 * isotropic_gens[0].n)
    isotropic_group = {Pauli(iden)}

    for r in range(1, len(isotropic_gens) + 1):
        for subset in combinations(isotropic_gens, r):
            if r == 1: 
                isotropic_group.add(subset[0])
            else:
                T = np.zeros(2 * isotropic_gens[0].n)
                for p in subset:
                    T = (T + p.T) % 2 
                isotropic_group.add(Pauli(T))

    return isotropic_group


def find_ksi_anticommuting_paulis(
    commutative_paulis_for_isotropic: set[Pauli], m: int
) -> set[Pauli]:
    """
    Finds the Pauli operators that commute with all isotropic generators and and anticommute with each other.

    Args:
        commutative_paulis_for_isotropic: The commutative Pauli's for the isotropic generators.
        m: Parameter that determines the number of anticommute Pauli's. Number of anticommute Pauli's is 2m + 1.

    Returns:
        set[Pauli]: A set of Pauli operators that anticommute with each other and commute with the isotropic generators.
    """

    ksi = 2 * m + 1

    commutative_paulis_for_isotropic = commutative_paulis_for_isotropic.copy() # Copy the set

    a = next(iter(commutative_paulis_for_isotropic)) # Get an element from the set
    
    anticommuting_paulis = {a}

    while len(anticommuting_paulis) < ksi:
        for p in anticommuting_paulis:
            # in each step remove the elements that commute with the anticommuting Paulis
            new_commutative_paulis_for_isotropic = set()
            for q in commutative_paulis_for_isotropic:
                if not is_commute(p, q):
                    new_commutative_paulis_for_isotropic.add(q)
            
            commutative_paulis_for_isotropic = new_commutative_paulis_for_isotropic

        anticommuting_paulis.add(next(iter(commutative_paulis_for_isotropic)))

    return anticommuting_paulis

class Cnc_Set:
    def __init__(self, isotropic_gens: list[Pauli], anticommuting_paulis: set[Pauli]) -> None:
        self.isotropic_gens = isotropic_gens.copy()
        self.anticommuting_paulis = anticommuting_paulis.copy()
        self.n = isotropic_gens[0].n
        self.m = (len(isotropic_gens) - 1) // 2

    def full_set(self) -> set[Pauli]:
        """
        Uses the isotropic generators and the anticommuting Paulis to generate the full closed noncontextual set.
        
        Returns:
            set[Pauli]: The full closed noncontextual set.
        """
        isotropic = generate_isotropic_from_gens(self.isotropic_gens)
        cnc = isotropic.copy()
        for p in self.anticommuting_paulis:
            for q in isotropic:
                cnc.add(p + q)
        
        return cnc

    def __str__(self) -> str:
        return f"Isotropic generators: {self.isotropic_gens}\nAnticommuting Paulis: {self.anticommuting_paulis}"
    
    def __repr__(self) -> str:
        return self.__str__()



def generate_cnc(n: int, m: int) -> Cnc_Set:
    """
    Generates the canonical noncontextual set for n qubits and with parameter m.

    Args:
        n: The number of qubits.
        m: The parameter that determines the dimension of the isotropic generators. Dimension is n - m.

    Returns:
        Cnc_Set: A closed noncontextual set for the given n and m.
    """
    isotropic_gens = generate_canonical_isotropic_gens(n, m)
    
    # Handle the edge case of m = 0
    if m == 0:
        iden = np.zeros(2 * n)
        return Cnc_Set(isotropic_gens, {Pauli(iden)})
    
    commutative_paulis_for_isotropic = find_commutative_paulis_for_isotropic(
        isotropic_gens,
    )

    isotropic = generate_isotropic_from_gens(isotropic_gens)  

    commutative_paulis_for_isotropic = commutative_paulis_for_isotropic.difference(isotropic)

    ksi_anticommuting_paulis = find_ksi_anticommuting_paulis(
        commutative_paulis_for_isotropic, m
    )    

    return Cnc_Set(isotropic_gens, ksi_anticommuting_paulis)

if __name__ == "__main__":
    cnc = generate_cnc(5, 2)
    print(cnc.full_set())
    print(cnc)