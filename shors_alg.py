# filename: shors_alg.py
# authors: Amanda Curtis and Arthur Diep-Nguyen

from qiskit.circuit import QuantumCircuit, QuantumRegister, AncillaRegister
from qiskit.quantum_info import Statevector
#import matplotlib.pyplot as plt

import numpy as np
import math

########################################################################
# Main oracle function
########################################################################

def oracle(a: int, N: int) -> QuantumCircuit:
    """Returns the oracle from Shor's algorithm.

    Parameters
    ----------
    a : int
        Positive integer
    N : int
        Positive integer strictly greater than `a`

    Returns
    -------
    QuantumCircuit
        Controlled circuit CU_a given by
            CU_a |x>_1 |y>_n
        that returns |x>_1 |a*y mod N>_n when x=1,
        where n = 2*ceil(log2(N))
    """
    #TODO
    pass

########################################################################
# Subcircuits
########################################################################

#TODO

########################################################################
# Other helper functions
########################################################################

#TODO
