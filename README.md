# shors-algorithm
A collaboration between Amanda Curtis and Arthur Diep-Nguyen

## Motivation

The oracle circuit we implement plays a key role is Shor's algorithm for prime factorization. When first developed, this landmark algorithm showcased the potential power of quantum computation, offering a superpolynomial speed up over the existing classical prime factorization algorithms. When presented in 1994, Shor's work had siginificant implications for cryptography - should a functioning quantum computer be able to run this algorithm, RSA could be broken. 

Many of the steps in this algorithm can be performed classically, but the quantum advantage appears in Step 4, where an order-finding sub-agorithm is needed.

## Contents

## Building Blocks and Outline 

In "Circuit for Shor's algorithm using $2n+3$ qubits", Beauregard offers a circuit construction for the order-finding sub-algorithm. We follow this construction, with some restrictions on the types of gates from which we are allowed to work.

We restrict our construction to:
1. 1-qubit gates
2. multi-controlled phase gates
3. X gates with any number of controls
4. Qiskit's built-in QFT (and its inverse)

The Oracle circuit requires a number of sub-circuits, and we first build those. The construction order starts from the most basic of component gates and their inverses - we begin with the standard adder gate, then create controled and doubly controlled versions of the adder and its inverse, which are used in construction of the modular adder gate (and its inverse). Finally, these are used to create the controlled multiplier gate and its inverse, which are the key componants of the circuit.

The best place to read more about our oracle and the process through which we built it is in the conversationally written jupyter notebook. The outline of the process is as follows:

1. Explanation of Process, Rationale
2. Building phiADD(a)
3. Building phiADD(a)mod(N)
4. Building CMult(a)mod(N)
5. The Complete Oracle

## Sources 

Our sources are also listed in section 6 of the jupyter notebook. 

## Points of Contact 

Amanda Curtis - dr.curtis.math (at) gmail (dot) com
Arthur Diep-Nguyen - 

