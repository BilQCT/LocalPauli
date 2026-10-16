# v75 update record

Date: 2026-10-12. Target: **LocalPauli-git**, not the earlier LocalPauli audit
workspace. Baseline commit: `6ada97e`. The supplied target already had local
changes to `.gitignore`, the combined-phase runner, and Finder metadata. The
source-code changes matched the preserved pre-audit snapshot; they were carried
forward. No commit, push, remote change, or manuscript edit was performed.

## Scientific corrections transferred

- Centralized termwise-modulo-3 / sum-modulo-2 seed arithmetic in
  `libs/MaxWeightSeed.jl`; updated the GAP orbit generator and its metadata.
- Corrected M dictionaries for n=2,3,4 and trace-one CNC dictionaries for n=2,3.
  Complete proper local-Clifford orbit sizes are 96, 13,824, 165,888.
- Corrected normalized state-vector and Tr(P rho) helpers and Pauli matrices.
  Current resources use T^n H^U E(G)|+>^n. Legacy runners require --legacy.
- Corrected the combined-phase example to use all stabilizer projectors.
- Added the direct Julia LP implementation and checked sparse primal/dual
  certificates, including serialization, dictionary hashes and residual checks.
- Added independent Python enumeration, 278 LP optima, four exact certificates,
  catalogue/state checks, and full-coordinate verification.
- Regenerated publication values and certificates in the standalone release
  copy and copied them into this folder. All values agree with the prior audit
  to nine decimals; the largest floating-point difference is 4.441e-16.
- Preserved the manuscript's 24-state order. There are 31 changed printed cells:
  two small-table L3 values, 14 four-qubit L1 values and 15 L2 values.

## Publication and portability work

- Removed runtime dependencies on ../v74 and ../submission-audit. A small
  `numerics/reference-v74.json` comparison fixture records the old printed
  values and reference hash; it is not used to calculate corrected optima.
- Added --prepare-data to rebuild the intentionally omitted large dictionaries
  without solving LPs. Full reproduction also rebuilds missing dictionaries.
- Made mathematical column-set equality authoritative for rebuilt HDF5 pools;
  byte-hash differences are explicitly reported as serialization provenance.
  Corrupt mathematical dictionaries still fail validation.
- Bound the manuscript exporter to hashes of the verified result, numerical
  certificate and exact certificate files, rejecting stale verification.
- Added fast regression tests covering both modular-arithmetic errors, graph
  matching, missing-file reconstruction, dictionary corruption, gate order,
  and normalization. Julia stale-file testing now uses a synthetic fixture
  instead of depending on a private archive.
- Rewrote README and updated numerical documentation, pinned the tested Python
  packages, retained the Julia Project and Manifest, and exposed the existing
  MyProject/src stub previously hidden by the broad ignore rule.
- Replaced the broad MyProject ignore rule with targeted exclusions; ignored
  environments, bytecode, OS files and local historical archives. Current
  results/certificates are deliberately not ignored.
- Moved old keys and tracked .DS_Store files into the local ignored archive;
  retained the old README as explicitly historical documentation. The two
  older locally_closed_L data files remain available and are clearly separated
  from current manuscript L3.
- Preserved LICENSE and the existing patent notice. No legal wording or
  ownership terms were changed.

## Evidence and integration

`v75-test-evaluation.md` describes tests, failure-injection checks and limits;
`test-results/` contains their outputs. `v75-file-changes.json` records changes
against the supplied folder, including SHA-256 hashes. The external integration
bundle contains the exact changed/new publication files, a deletion manifest,
and a hash-checking installer. It includes no .git directory or private archive.
The installer defaults to inspection and refuses conflicting newer files.
It accepts either the supplied-folder version or the baseline Git version for
modified/deleted files, and treats an already-updated file as a no-op.

Local backups live in `archive/pre_v75_localpauli_git/` and are not needed for
reproduction. The integration bundle does not transfer old ignored keys; the
listed obsolete tracked keys are explicitly removed with backup on application.
The earlier scientific audit report is retained as a historical report; paths
in it describe the earlier audit workspace. Current setup instructions are in
the root README, and current exported manuscript artifacts are under results/.
