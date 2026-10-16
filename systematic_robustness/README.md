> Current publication results and reproduction instructions are in
> [../numerics/README.md](../numerics/README.md). The former M dictionaries were
> replaced on 2026-10-12; old M files are rejected by the loader.

# Systematic bicolored-resource robustness

This folder provides a self-contained workflow for computing robustness of the
reduced resource states emphasized in the manuscript,

```math
|\widetilde\psi_{G,U}\rangle
  = T^{\otimes n}H^{\otimes U}E(G)|+\rangle^{\otimes n},
```

relative to a **fixed, explicitly supplied collection of input vectors**. It
does not use the older `E(G)T^U|+>^n` resource-state convention.

## Files

- `BicoloredRobustness.jl`: state generation, dictionary loading, LP solving,
  validation, batch evaluation, and result persistence.
- `run_batch.jl`: command-line batch runner.
- `test/runtests.jl`: checks state conventions against the existing saved
  states and reproduces a known two-qubit robustness value.

The Pauli convention is lexicographic `I,X,Y,Z`, and a state vector contains
`Tr(P rho)`. Thus the identity entry is one.

## Run the tests

From `LocalPauli`:

```sh
julia --project=MyProject systematic_robustness/test/runtests.jl
```

## Run a batch

The runner takes a qubit count, output file, dictionary label, and one or more
HDF5 matrices. Every input file must contain a `data` matrix whose columns are
the fixed input vectors.

```sh
julia --project=MyProject systematic_robustness/run_batch.jl \
  4 results/L2_n4.h5 L2 \
  data/stabilizer_states_4.h5 \
  data/deterministic_4.h5 \
  data/max_weight_operators_p_2_4.h5
```

The dictionary is loaded once and is not changed between states. We use the
same compact names as the manuscript:

- `L0 = D` (written $\mathcal L_0$ in the manuscript);
- `L1 = conv(D ∪ M2)` (written $\mathcal L_1$);
- `L2 = conv(S ∪ D ∪ M2)` (written $\mathcal L_2$);
- `L3 = conv(C ∪ D ∪ M2)` (written $\mathcal L_3$).

Here `D`, `S`, `C`, and `M2` are the deterministic, stabilizer, CNC, and
degree-two max-weight dictionaries. Combining a previous phase space with new
locally closed vectors makes the comparison nested; in particular, `L2` can
never have larger robustness than `S` on the same target.

The result file records:

- the dictionary label, sources, and column count;
- the state-family and Pauli conventions;
- graph ids, hashes, edges, and decorated vertices;
- objectives and numerical residuals;
- solver statuses, checked sparse decompositions, full dual witnesses, dual feasibility and duality gaps;
- SHA-256 hashes of every input dictionary.

An objective computed with a restricted dictionary is the exact LP optimum for
that dictionary, but only an upper bound on robustness relative to a larger
phase space.
