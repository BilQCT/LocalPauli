# LocalPauli

Companion code for **Classical simulation of universal measurement-based quantum
computation using multipartite Bell scenarios**, Okay, Yücel and Ipek,
[arXiv:2410.23734](https://arxiv.org/abs/2410.23734).

## Corrected v75 publication numerics

The current pipeline and certified results are in `numerics/` and
`results/manuscript_v75/`. They implement the manuscript seed

```
gamma_p(u) = sum(product(u_i for i in I) mod 3 for |I|=p) mod 2
```

Each product is reduced modulo 3 and the sum modulo 2. Earlier generators
reduced the running sum modulo 3, producing different dictionaries.
See [the numerical audit](docs/numerical-audit.md),
[the change record](docs/v75-change-record.md), and
[the test evaluation](docs/v75-test-evaluation.md).

## Reproduce from a standalone checkout

Python 3.12 is the tested version. Use an isolated environment:

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r numerics/requirements-tested.txt
python -B numerics/test_numerics.py
python -B numerics/reproduce.py --prepare-data
python -B numerics/reproduce.py --verify-only
python -B numerics/export_manuscript.py
```

`--prepare-data` independently reconstructs the omitted large dictionaries and
checks the graph catalogues and saved targets. `--verify-only` checks the supplied
certificates without re-solving the optimization problems. To recompute all
278 optima and four exact certificates, run:

```sh
python -B numerics/reproduce.py
python -B numerics/export_manuscript.py
```

Full reproduction also prepares missing dictionaries. No manuscript source,
parent-directory files, historical archive, credentials, or network connection
is needed after dependencies are installed. The original displayed values and
manuscript hash are in `numerics/reference-v74.json`, solely for comparing the
revision. Copyable LaTeX and the comparison guide are exported to
`results/manuscript_v75/manuscript/`. Do not run with `python -O`:
assertions are part of the certificate checks.

## Conventions and scope

The resource state is `T^n H^U E(G)|+>^n`. Coordinates are `Tr(P rho)` in
lexicographic I,X,Y,Z order, with identity coefficient +1. The manuscript
uses independent Hadamard decorations: 2, 6 and 24 inequivalent states at
n=2,3,4. The computation covers all 3, 10 and 50 colored representatives.

D contains all deterministic operators; S all stabilizer projectors, including
entangled states; C all maximal CNC operators. M is the orbit of the quadratic
seed under tensor products of single-qubit Clifford conjugations, with sizes
96, 13,824 and 165,888. This local-Clifford orbit differs from the full Bell
relabeling orbit. Here "local Clifford" is distinct from "locally closed".
L1=D union M, L2=S union D union M, and L3=C union D union M.

Every result has a primal decomposition and dual certificate. The verifier
checks the full physical coordinates and every original dictionary inequality,
including for the four-qubit symmetry reduction. It accepts errors below 1e-8.
Four small-system results also have exact certificates over Q(sqrt(2)).
See [numerical documentation](numerics/README.md) for the data format.

## Julia cross-checks

The existing Julia project and manifest are retained (Julia 1.11.4).
Install its dependencies, prepare the large dictionaries using Python as above,
then run:

```sh
julia --project=MyProject -e 'using Pkg; Pkg.instantiate()'
julia --project=MyProject systematic_robustness/test/runtests.jl
julia --project=MyProject systematic_robustness/test/legacy_helpers.jl
```

`systematic_robustness/` supplies direct GLPK LPs with saved primal/dual witnesses
and dictionary hashes. `libs/MaxWeightSeed.jl` centralizes the seed and rejects
unversioned max-weight files. The GAP generator provides an independent
implementation of local-Clifford orbits.

## Historical data

`data/locally_closed_L_3.h5` and `locally_closed_L_prime_3.h5` are older,
different dictionaries; their filenames do not denote the current manuscript's
L3. Legacy `E(G)T^U|+>` runners require `--legacy`. The old README is retained
as [historical documentation](docs/historical-readme.md), not current usage.
Old numerical keys were removed from the publication tree and preserved locally
under ignored `archive/`; no current calculation requires that archive.
Large regenerable dictionaries are ignored; current results and certificates
are deliberately included. Keep Git metadata and local archives out of transfer
bundles. Integration instructions accompany the update bundle.

## Citation

Please cite arXiv:2410.23734 and identify the repository release or commit used.
The internal manuscript designation "v75" is not an arXiv version number.

# Patent Pending Notice

This repository contains material related to our arXiv preprint arXiv:2410.23734. Please note that specific aspects of the technology described herein is subject to a pending patent (application number: 18/925,447).
Licensed under MIT, which does not grant patent rights.

For inquiries: cihan.okay@bilkent.edu.tr.
