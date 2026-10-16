# Numerical validation

The calculations accompanying [arXiv:2410.23734](https://arxiv.org/abs/2410.23734)
are checked with primal decompositions and dual optimality witnesses. The
numerical tolerance is 1e-8. Four small-system cases additionally have exact
certificates over Q(sqrt(2)).

Validation covers all 278 optima, dictionary membership, full-coordinate target
reconstruction, dual feasibility and primal-dual gaps. For four qubits the
reduced calculation is checked in all 256 physical Pauli coordinates. The
verifier also checks phase-space monotonicity and the 22 strict improvements
among the 24 independent-decoration four-qubit resources.

Independent enumeration checks the graph catalogues through n=4. Two independent
state constructions check all 417 stored targets through n=5. The regression
suite covers both modular-arithmetic errors, decoration matching, reconstruction
of missing dictionaries, rejection of corrupt dictionaries, gate order and
normalization. Julia tests additionally check LP certificates and their saved
representation. The stale-file fixture does not require historical data.

The release calculations were tested with Python 3.12.2, NumPy 1.26.4,
SciPy 1.17.1, SymPy 1.14.0, h5py 3.16.0 and Julia 1.11.4 on macOS.
The pinned dependencies are in `numerics/requirements-tested.txt`.
Standalone reproduction and intentional corruption checks passed: a changed
saved dual is rejected, and table export rejects results modified after
verification. `results/robustness/verification.json` contains the current
per-certificate diagnostics and hashes.

A subsequent public-path cleanup retained all numerical values and numerical/
exact certificates, and reran the Python regressions, certificate verifier and
exporter. It removed dependence on internal manuscript comparisons.

These checks do not constitute a proof audit of every theorem in the paper.
No n=5 robustness calculation or cross-platform test matrix is claimed.
Dependencies were installed for testing; fresh network installation is separate.
