# Operator and state conventions

The current numerical workflow is documented in [numerics/README.md](../numerics/README.md).
Its definitions follow the resource-state and phase-space comparisons in
[arXiv:2410.23734](https://arxiv.org/abs/2410.23734).

The older files `data/locally_closed_L_3.h5` and
`data/locally_closed_L_prime_3.h5` describe, respectively, tensor-product
constructions from one-/two-qubit local Pauli vertices and the union of
deterministic and CNC operators. Their names do not denote the current L3,
which is C union D union M. They remain available for studying those sets.

The older state-generation scripts use E(G)T^U|+>^n and require `--legacy`.
The numerical comparison workflow uses T^n H^U E(G)|+>^n; U therefore labels
Hadamard decorations in that workflow.

Earlier max-weight generators incorrectly reduced the running sum modulo 3.
The corrected seed reduces each monomial modulo 3 and the total modulo 2.
Current max-weight files carry convention metadata; loaders reject unversioned
files. Some earlier CNC files had trace-minus-one columns; current C operators
are normalized to trace one. Whole-column sign changes preserve the signed l1
optimization, but normalization is necessary for the stated operator convention.
