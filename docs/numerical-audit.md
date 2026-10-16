> Original scientific audit from the earlier LocalPauli workspace. Its paths identify that historical workspace. For this standalone release use the root README and v75-test-evaluation.md.

# Numerical audit for the v75 revision

Reference: `v74/main.tex`, supplied 12 October 2026. The reference manuscript is unchanged. Apply the accompanying **v75-manuscript-changes.txt** to prepare v75; it contains replacement text, equations, table rows, and plot coordinates. `manuscript-value-changes.csv` records all 104 displayed values before and after correction.

## Findings and corrections

1. **The seed arithmetic differed between code and manuscript.** The manuscript evaluates each monomial modulo 3 and adds the resulting terms modulo 2:
   `gamma_p(u) = [sum_{i1<...<ip} [u_i1 ... u_ip]_3]_2`.
   The old generator reduced the running sum modulo 3 before taking its parity. For p=2 and u=(1,1,1), these give 1 and 0, respectively. The corrected implementation uses the manuscript formula. A common implementation now supplies the Julia generators; independently implemented Python enumeration reproduces the corrected pools.

   Hoban et al., *Generalized Bell Inequality Experiments and Computation*, Section II, establishes the arithmetic conventions and general framework: https://arxiv.org/pdf/1108.4798. It does not uniquely select this particular elementary-polynomial seed. Accordingly, the correction is to make the numerical implementation agree with the manuscript's explicitly specified seed, rather than to attribute a unique seed formula to Hoban et al.

2. **The numerical symmetry group must be explicit.** The revised calculations use the complete orbit under products of single-qubit Clifford conjugations. The corrected M sizes are 96, 13,824, and 165,888 for n=2,3,4. The old n=4 pool had 55,296 columns. Equal orbit sizes at n=3 did not imply equal pools: the old and corrected three-qubit pools are disjoint. The Bell relabeling group is different; its two-qubit seed orbit has 192 elements, whereas the proper local-Clifford orbit has 96. The guide defines M without ambiguity. No claim that Bell-orbit robustness equals these LC-orbit values is made.

3. **The stored three-qubit CNC pool contained 36,864 trace-minus-one columns.** Multiplying each such entire column by -1 produces precisely the independently enumerated trace-one CNC pool. Whole-column sign changes leave the signed l1 optimization unchanged, explaining why C values survive. The C files are now normalized to the manuscript convention.

4. **Legacy state and coordinate helpers were inconsistent with the current convention.** The reusable helpers now use normalized state vectors and coordinates Tr(P rho), including identity coefficient +1. Undefined Pauli-matrix dependencies were removed. The current resource family is T^n H^U E(G)|+>^n. Historical E(G)T^U|+>^n runners require an explicit legacy flag, so they cannot silently supply the current manuscript targets. The combined-phase example now uses all stabilizer projectors for S, including entangled stabilizers.

5. **Historical results needed separation from current outputs.** Old sources, dictionaries, and numerical keys are preserved under `LocalPauli/archive/pre_termwise_seed_20261012/`. The archive is excluded from Git because of its size; preserve it separately if transferring the workspace. Current publication results live only under `LocalPauli/results/manuscript_v75/`. The current Julia loaders reject unversioned max-weight files. New saved Julia certificates include dictionary hashes, primal decompositions, dual witnesses, and residuals. The deterministic generator defaults to the required n=2:4 range rather than attempting n=7.

6. **A manuscript parity argument is false.** Total truth-table parity is not invariant under independent outcome relabelings. Flipping one setting negates three entries in the two-qubit table, changing its parity while staying in the same Bell orbit. `parity-counterexample.json` records the explicit example. The guide replaces the erroneous orbit-distinguishing argument. The orbit sizes and robustness computations do not depend on this argument.

## Numerical impact

Of 104 displayed robustness entries, **31 change**: two small-table L3 entries, fourteen four-qubit L1 entries, and fifteen four-qubit L2 entries. The displayed S, C, and D values remain unchanged at their printed precision.

For undecorated three-qubit states:

| State | Old displayed L3 | Correct L3 | Exact expression |
|---|---:|---:|---|
| P3 | 1.183 | 1.15533008589 | (5 + 3 sqrt(2))/8 |
| K3 | 1.190 | 1.13639610307 | 1/2 + 9 sqrt(2)/20 |

Examples of revised four-qubit L1=L2 values are 1.445431759768 for the undecorated star, 1.957106781187 for undecorated P4, and 1.260200459405 for undecorated K4. The complete replacement series is in the text guide and `v75-four-qubit-plots.tex`.

The strict L2 improvement over S for **22 of the 24** displayed four-qubit resources remains true. Cases 5 and 10 have S=L2=1.

## Verification completed

- Recomputed 278 LP optima: six dictionaries for all 3 two-qubit and 10 three-qubit colored-graph representatives, and four dictionaries for all 50 four-qubit representatives. The manuscript uses the independent-set subsets of sizes 2, 6, and 24.
- Every optimum has a saved primal decomposition and dual witness. Verification rebuilds the dictionaries, checks membership, reconstructs the physical target in every Pauli coordinate, checks every full-dictionary dual inequality, and checks the primal-dual gap. The four-qubit symmetry reduction is checked by explicitly averaging primal operators in the full 256-coordinate space.
- Four additional exact certificates over Q(sqrt(2)) verify P3 and K3 in C and L3. Other reported optima have numerical certificates; they are not claimed as symbolic exact results.
- The final independent verification passed all 278 numerical and four exact certificates. Maximum reconstruction residual: 1.068e-13; maximum dual-constraint violation: 1.248e-13; maximum primal-dual gap: 1.954e-14. The acceptance tolerance is 1e-8.
- All 417 saved target states through n=5 agree with an independent normalized state-vector construction, with maximum coordinate discrepancy 1.11e-15. This is a state check, not a recomputation of n=5 robustness.
- Exhaustive connected colored-graph enumeration through n=4 confirms both unrestricted and independent-set catalogue counts and the manuscript-to-record mapping.
- Corrected Julia and Python generators agree on the entire n=2 and n=3 max-weight column sets. The Julia tests additionally check the n=4 seed and expected pool dimensions. A second complete Julia n=4 orbit enumeration was not required or performed.
- Julia regression tests cover state construction, saved targets, robustness certificates, certificate serialization, mixed seed arithmetic, stale-file rejection, and trace-one CNC operators. Separate legacy-helper tests check normalized gates and Pauli coordinates.
- Exact two-qubit nonsignaling enumeration yields 1,408 vertices in orbits of sizes 64, 576, 288, 288, and 192. Exact enumeration of the two-qubit Pauli polytope from the 60 stabilizer facets yields 22,320 vertices, agreeing with the manuscript.

## Reproduction and evidence

From the LocalPauli directory, with the dependencies in `numerics/requirements.txt` installed:

```text
python -B numerics/reproduce.py
python -B numerics/reproduce.py --verify-only
python -B numerics/export_manuscript.py
```

Do not use Python's optimization flag: assertions form part of the verification. The first command regenerates corrected dictionaries and solves the LPs. Large ignored dictionary files can be rebuilt; existing S/D dictionaries are checked against independent enumerations. See `LocalPauli/numerics/README.md` for the conventions and saved-certificate format.

Principal evidence:

- `LocalPauli/results/manuscript_v75/results.json` and `results.csv`: corrected optima, input hashes, runtime provenance, and target metadata.
- `certificates.npz`, `exact_small.json`, and `verification.json` in that directory: numerical and exact certificates and verification diagnostics.
- `submission-audit/final-verification.log`, `julia-tests.log`, `legacy-helper-tests.log`, `julia-python-orbit-agreement.json`, `ns232-verification.log`, and `lambda2-verification.log`: independent check outputs.
- `submission-audit/provenance-sha256.json`: hashes of the revised numerical sources, publication outputs, guide, and reference manuscript.

## Remaining manuscript work and limits

Apply the guide to create v75, then compile and inspect the revised manuscript, especially the replaced table and plot. This audit did not edit or compile a revised full manuscript. The guide preserves the resource numbering and graph drawings.

The conclusion's unspecified claim about “initial numerical results” at n=3 needs separately identifiable supporting evidence; the robustness certificates do not by themselves establish that claim. Remove it if that evidence cannot be supplied. The audit checks the numerical pipeline and the specific mathematical issues identified above; it is not a proof audit of all simulation theorems. Robustness improvements support the stated additive-error estimation bounds and should not be advertised as an independently established total-runtime improvement or an exact sampler for signed decompositions.
