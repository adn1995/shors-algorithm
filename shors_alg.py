# filename: shors_alg.py
# authors: Amanda Curtis and Arthur Diep-Nguyen

from qiskit.circuit import QuantumCircuit, QuantumRegister, AncillaRegister
from qiskit.quantum_info import Statevector
from qiskit.circuit.library import QFT
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
    #TODO see Figure 7
    #TODO don't forget to use to_gate() when composing subcircuits
    pass

########################################################################
# Subcircuits
########################################################################

def adder(a: int, N: int) -> QuantumCircuit:
    #TODO see Section 2.1 and Figure 3
    # Must first "solve" for n, set up registers needed
    # Based off of cited paper and QC Bootcamp 2's implementation of 
    # Draper's adder circuit 
    # and class lecture on 22 May 2025

    # We utilize the corollary that relates the QFT, A_k (Draper), 
    # and this (controlled) phase gate P_n(a) (phiADD(a)) 

    # "Solving" for n 
    n = int(np.ceil(np.log2(N)))

    # Setting up Quantum Register 
    quantum_register = QuantumRegister(size=n, name ='x')
    phi_add_a = QuantumCircuit(quantum_register, name="phi add a")

    # Building P_n(a) by making a phase gate p
    # for each qubit 
    for idx, q in enumerate(reversed(quantum_register)):
        phi_add_a.p(np.pi * a / (1 << idx), q)
    

    #####
    

    
    pass

def cc_adder_mod(a: int, N: int) -> QuantumCircuit:
    #TODO doubly controlled, see Section 2.2 and Figure 5
    pass

def c_mult_mod(a: int, N: int) -> QuantumCircuit:
    #TODO controlled, see Section 2.3 and Figure 6
    pass

#TODO we also need to write a controlled swap circuit, since we are
# probably not allowed to use the built in swap gate :(

########################################################################
# Other helper functions
########################################################################

#TODO anything that is not a circuit can go here
