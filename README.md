# LocalPauli

This repository accompanies **Classical simulation of universal measurement-based
quantum computation using multipartite Bell scenarios**, by Okay, Yücel and Ipek
([arXiv:2410.23734](https://arxiv.org/abs/2410.23734)). It provides operator
constructions, robustness calculations, and numerical verification tools for
the phase-space simulation framework developed in the paper.

## Background

The paper develops a probabilistic representation of **measurement-based Pauli
computation (MBPC)**, a universal model of quantum computation using adaptive
single-qubit Pauli measurements on suitable resource states. This model connects
the measurement-based approach to quantum computation with magic-state resources;
see [Danos, Kashefi and Panangaden](https://arxiv.org/abs/0704.1263) and
[Raussendorf et al.](https://arxiv.org/abs/1511.08506).

Pauli measurements on a graph state alone are efficiently simulable by stabilizer
methods. To support universal computation, the resource states incorporate
non-stabilizer resources. Given a graph $G=(V,E)$ on $n$ vertices and a subset
$U\subseteq V$, the paper introduces *magic cluster states* of the form

```math
|\psi_{G,U}\rangle = E(G)T^{\otimes U}|+\rangle^{\otimes n},
\qquad E(G)=\prod_{(i,j)\in E}CZ_{ij},
```

where $T^{\otimes U}$ applies a $T$ gate to each qubit in $U$ and the identity
elsewhere. For the numerical comparisons, the repository uses the reduced
resource family

```math
|\phi_{G,U}\rangle=T^{\otimes n}H^{\otimes U}E(G)|+\rangle^{\otimes n}.
```

In this second expression, $U$ labels the Hadamard decorations. The two
conventions are distinguished in the calculation scripts and
[numerical documentation](numerics/README.md).

The central geometric object is the **local Pauli polytope** $\mathrm{LP}_n$,
a polytope of trace-one Hermitian operators that contains the quantum state
space. It is identified with the nonsignaling polytope of the $(n,3,2)$ Bell
scenario: $n$ parties, three measurement settings per party, and two outcomes
per setting. The three settings correspond to the single-qubit Pauli
measurements $X$, $Y$, and $Z$.

A useful structural property is that tensor products preserve vertices:
if $A_1$ and $A_2$ are vertices of $\mathrm{LP}_{n_1}$ and
$\mathrm{LP}_{n_2}$, respectively, then $A_1\otimes A_2$ is a vertex of
$\mathrm{LP}_{n_1+n_2}$.

The paper introduces **locally closed operators**. Two Pauli operators locally
commute when their factors commute on every qubit. A set of Pauli operators is
locally closed when it contains the product of every pair of locally commuting
operators in the set. Together with compatible value assignments, these sets
define phase-space operators adapted to local Pauli measurements.

For a single qubit, $\mathrm{LP}_1$ is a cube enclosing the Bloch sphere. Tensor
products of its eight vertices give the deterministic vertices of the
$(n,3,2)$ Bell polytope. All vertices of $\mathrm{LP}_2$ are also locally closed
operators, providing further constructions through tensor products.

Efficiently represented locally closed operators, together with their measurement
update rules, support classical simulation through phase-space decompositions.
This framework includes the closed noncontextual (CNC) representations of
[Raussendorf et al.](https://arxiv.org/abs/1905.05374), which in turn include
stabilizer representations.

For an operator set $`\mathcal{A}=\{A_\alpha\}`$, the associated robustness is

```math
\mathfrak{R}_{\mathcal{A}}(\rho)
=\min\left\{\sum_\alpha |q_\alpha|:
\rho=\sum_\alpha q_\alpha A_\alpha,\quad q_\alpha\in\mathbb{R}\right\}.
```

The coefficients may be negative. Their minimum total absolute weight enters
the sampling cost of the estimation algorithms studied in the paper.

## Numerical comparisons

The calculations compare the following operator sets:

- $D$: deterministic vertices.
- $S$: all stabilizer projectors, including entangled stabilizer states.
- $C$: CNC operators.
- $M$: the local-Clifford orbit of the quadratic max-weight operator defined
  in the paper.

The combined sets are $L_1=D\cup M$, $L_2=S\cup D\cup M$, and
$L_3=C\cup D\cup M$. Here **local Clifford** refers to conjugation by tensor
products of single-qubit Clifford unitaries; it is distinct from **locally
closed**, the property defined above.

The repository reproduces the paper's two-, three-, and four-qubit robustness
comparisons and supplies decompositions and optimality certificates for checking
the results. The older sets named $L$ and $L'$ are described in the
[operator conventions](docs/operator-conventions.md).

## Getting started

Run the following commands from the repository's root directory. Python 3.12 is
the tested version; Julia is optional for the Python workflow.

Create an environment and install the tested dependencies (macOS/Linux):

```sh
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r numerics/requirements-tested.txt
```

Run the regression tests, generate the required operator data, and verify the
supplied results:

```sh
python -B numerics/test_numerics.py
python -B numerics/reproduce.py --prepare-data
python -B numerics/reproduce.py --verify-only
```

Preparation rebuilds the large operator files omitted from the repository.
Verification checks the saved decompositions and optimality witnesses without
solving all optimization problems again. The supplied results contain 278
numerical certificates and four exact certificates.

To recompute the results from scratch:

```sh
python -B numerics/reproduce.py
```

Numerical results and optimality certificates are available in
[results/robustness/](results/robustness/). Export tables and plots with
`python -B numerics/export_tables.py`; outputs are written to
`results/robustness/exports/`.
The calculations run independently of the manuscript source. For the theoretical
framework and interpretation of the results, see the
[paper on arXiv](https://arxiv.org/abs/2410.23734).

Run these scripts without Python's `-O` flag, because assertions are used for
verification.

## Repository contents

| Directory | Contents |
|---|---|
| [numerics/](numerics/) | Python enumeration, robustness calculations, and certificate verification |
| [data/](data/) | Graph catalogues, resource states, and operator dictionaries |
| [results/](results/) | Numerical results and optimality certificates |
| [systematic_robustness/](systematic_robustness/) | Julia linear-programming implementation and tests |
| [libs/](libs/) | Julia tools for Pauli operators and local-Clifford actions |
| [local_robustness/](local_robustness/) | Julia generators and calculation scripts |
| [docs/](docs/) | Operator conventions and numerical validation |

See the [numerical documentation](numerics/README.md) for the resource-state
conventions, operator sets, file formats, and scope of the calculations. The
[validation documentation](docs/validation.md) describes the checks and their
limitations. Earlier datasets and state conventions are documented separately
in the [operator conventions](docs/operator-conventions.md).

## Julia implementation

The Julia tools provide an additional implementation of operator generation and
robustness optimization. The project was tested with Julia 1.11.4.

After preparing the operator data using the Python instructions above, install
the Julia dependencies and run the tests:

```sh
julia --project=MyProject -e 'using Pkg; Pkg.instantiate()'
julia --project=MyProject systematic_robustness/test/runtests.jl
julia --project=MyProject systematic_robustness/test/legacy_helpers.jl
```

See [systematic_robustness/README.md](systematic_robustness/README.md) for running
calculations with the Julia implementation.

## Citation

If you use this code or its numerical results, please cite
[arXiv:2410.23734](https://arxiv.org/abs/2410.23734) and identify the repository
release or commit used.

# Patent Pending Notice

This repository contains material related to our arXiv preprint arXiv:2410.23734. Please note that specific aspects of the technology described herein is subject to a pending patent (application number: 18/925,447).
Licensed under MIT, which does not grant patent rights.

For inquiries: cihan.okay@bilkent.edu.tr.
