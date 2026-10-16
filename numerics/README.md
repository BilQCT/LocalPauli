# Robustness calculations

These scripts reproduce the numerical comparisons in
[arXiv:2410.23734](https://arxiv.org/abs/2410.23734). No manuscript source or
internal revision files are required.

## Reproduce and verify

From the repository root, using Python 3.12 and the packages in
`numerics/requirements-tested.txt`:

```sh
python -B numerics/test_numerics.py
python -B numerics/reproduce.py --prepare-data
python -B numerics/reproduce.py --verify-only
python -B numerics/export_tables.py
```

Preparation enumerates the dictionaries, rebuilding the large ignored files,
and checks the saved targets and graph catalogues. Verification checks the
supplied certificates without solving the optimization problems again. To
recompute all results, run `python -B numerics/reproduce.py` instead.
Do not use Python's `-O` option: assertions are part of the verification.

## Conventions

- Resource states: `T^n H^U E(G)|+>^n`, where U labels Hadamard decorations.
- Coordinates: `Tr(P A)`, lexicographic I,X,Y,Z tensor order; identity is +1.
- D: all 8^n deterministic operators.
- S: all stabilizer projectors, including entangled states.
- C: all maximal CNC operators, normalized to trace one.
- M: the complete orbit of the quadratic max-weight seed under tensor products
  of single-qubit Clifford conjugations. X,Y,Z correspond to inputs 0,1,2.
  Each pairwise product is evaluated modulo 3, and the sum modulo 2.
  The orbit sizes for n=2,3,4 are 96, 13,824 and 165,888.
- L1=D union M; L2=S union D union M; L3=C union D union M.
- Objective: minimum sum of absolute real coefficients reconstructing the target
  in every Pauli coordinate, including the identity.

All 3+10+50 connected colored-graph representatives through n=4 are computed.
The paper's independent-decoration subset has 2+6+24 states. Catalogue record IDs
are zero-based; exports match graphs and decorations to the paper's ordering.

## Outputs

`results/robustness/` contains:

- `results.json` and `results.csv`: 278 optima, conventions and input provenance.
- `certificates.npz`: primal supports, coefficients, targets and full dual
  witnesses. Load with `allow_pickle=False`.
- `exact_small.json`: four exact certificates over Q(sqrt(2)) for P3/K3 in C/L3.
- `verification.json`: full-coordinate errors, duality gaps, and verified hashes.
- `preparation.json`: target and catalogue checks.
- `exports/`: generated table rows, four-qubit plot coordinates and a CSV.

`robustness_tables.tex` combines the generated small-system table and four-qubit
figure. The exporter rejects results whose hashes no longer match verification.
The four-qubit verifier explicitly averages primal operators under the target's
16-element Clifford symmetry group, reconstructs all 256 physical coordinates,
and checks every original dual inequality. See `four_core.py` for the reduction.

HDF5 byte hashes record serialization provenance. A regenerated file may have
identical mathematical columns but different bytes or ordering. Such differences
are reported as `serialization_changes`; complete column-set comparison against
independent enumeration remains mandatory. Numerical certificates are accepted
only when reconstruction, dual violations and gaps are below 1e-8.

Saved target states through n=5 are also checked, without claiming n=5 robustness
calculations. See [validation](../docs/validation.md) for test coverage.
