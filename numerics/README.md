# Corrected publication numerics

This is the current numerical pipeline for the v75 revision, using the mixed
modular arithmetic stated in v74. The reusable enumeration and exact-certificate
helpers were adapted from the independently verified `robustness-review-v75`
project. This standalone release includes the required helpers; the external review folder is not required.

## Reproduce

Use Python 3.11 or newer with the packages in `requirements.txt`:

```sh
python -m pip install -r numerics/requirements.txt
python -B numerics/reproduce.py
python -B numerics/reproduce.py --verify-only
```

Run from the LocalPauli repository root. Do not use `python -O`: assertions are
part of the verification. A full run writes corrected M files and trace-one C
files under `data/` and replaces only `results/manuscript_v75/` outputs. It checks
existing S/D dictionaries against independently enumerated pools. Missing
S/D/C dictionaries are rebuilt from the enumerators, including the large files
excluded from Git. The graph catalogues and saved target files under `data/`
are used for independent cross-checks. Original M/C files from this correction
are retained in the historical archive.

## Precisely defined calculation

- Target: `T^n H^U E(G)|+>^n`; U denotes Hadamard decorations.
- Coordinates: `Tr(P A)`, lexicographic tensor order I,X,Y,Z. Identity is +1.
- D: all 8^n products of the eight single-qubit deterministic operators.
- S: all stabilizer projectors, not just the 6^n product stabilizers.
- C: all maximal CNC operators, with trace +1.
- M: complete proper local-Clifford orbit of the symmetric p=2 seed. Each
  monomial is multiplied modulo 3; the sum is modulo 2. X,Y,Z label inputs 0,1,2.
- L1=D union M; L2=S union D union M; L3=C union D union M.
- M sizes at n=2,3,4: 96, 13,824, 165,888. Bell relabelings are not included.
- Objective: minimum sum of absolute real decomposition weights, with every
  physical Pauli-coordinate equality enforced, including the identity row.

All 3+10+50 catalogue states are recomputed. The manuscript subset consists of
2+6+24 independent-set decorations. Zero-based record IDs and vertex labels in
the JSON/CSV are not the one-based four-qubit indices in the manuscript; match
edges and decorations. The export helper handles this mapping explicitly.

## Results and evidence

`results.json` and `results.csv` contain 278 LP optima: six phase spaces per
small-system target and four per four-qubit target. `certificates.npz` stores
primal supports, coefficients, full dual witnesses, and targets. Load it with
`allow_pickle=False`. The four-qubit primal seeds must be averaged under the
target's 16-element Clifford symmetry group; the verifier explicitly performs
this average in all 256 physical coordinates. It also checks membership and
every original dual inequality. See the mathematical justification in
`four_core.py`.

`exact_small.json` contains exact Q(sqrt(2)) primal/dual proofs for P3/K3 in C
and L3. They are verified against freshly enumerated pools and symbolic targets.
Other values are numerical certificates with a 1e-8 acceptance threshold.

`verification.json` records all residuals, dual violations and gaps. Every result records input dictionary SHA-256 hashes, conventions and the original v74
source hash. Byte hashes are provenance; mathematical dictionary equality is checked independently (see below). The code checks phase-space monotonicity and the 22/24 improvement.

The graph catalogues are independently checked by exhaustive enumeration of
connected labeled graphs and color-preserving isomorphisms through n=4. Saved
target states through n=5 are checked against direct state-vector expectations;
this does not claim to recompute all five-qubit robustness values.

Julia's original GAP orbit generator and direct GLPK LP solver remain available
as independent cross-checks. The tests reject old M metadata and test both
modulo-arithmetic pitfalls. The old results are not silently reused.

## Standalone release workflow

Use `requirements-tested.txt` for the evaluated Python 3.12 environment;
`requirements.txt` gives broader, untested compatibility ranges. In a fresh
checkout, run `python -B numerics/reproduce.py --prepare-data` before
`--verify-only`. This creates the three intentionally ignored large files.
Full reproduction prepares them automatically.

No manuscript file is required. `reference-v74.json` holds only the original
printed values and source hash; it is comparison data, not a source of corrected
optima. `export_manuscript.py` writes inside
`results/manuscript_v75/manuscript/`.

HDF5 byte hashes record input serialization provenance. Rebuilt files can have
different metadata or column ordering with exactly the same operator sets.
The verifier reports such differences as `serialization_changes` and checks
complete column-set equality against independent enumeration. The actual
certificates are checked against those enumerated pools. A changed mathematical
dictionary is rejected. Historical unused files in an older hash manifest are
not certified by this verification; only the documented D/S/C/M pools are.
