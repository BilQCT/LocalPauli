# v75 tests and evaluation

Evaluated 2026-10-12 on macOS, using Python 3.12.2 with the exact package versions
in `numerics/requirements-tested.txt`, and the existing Julia 1.11.4 environment.
The scientific release candidate was copied to a standalone temporary directory
from tracked files plus intended nonignored additions. It had no sibling
manuscript, audit folder, historical archive, or ignored large dictionaries.

| Check | Outcome | What it establishes |
|---|---|---|
| Standalone full regeneration | PASS: 278 LPs, 4 exact certificates | The released numerical inputs and code suffice without private workspace files. Missing large dictionaries rebuild correctly. |
| Independent full-coordinate certificate verification | PASS: all 278 numerical and 4 exact cases | Checks operator membership, target reconstruction, full dual feasibility, objective/gap and symmetry averaging. |
| Separate --prepare-data invocation | PASS | Dictionary preparation and target/catalogue checks work independently of LP solving. |
| Comparison with earlier audited optima | PASS: maximum difference 4.441e-16 | All 278 values agree to nine decimals; three differ in last binary floating-point bits. |
| Manuscript export | PASS: 104 displayed entries, 31 changes | Preserves graph/decorations mapping; 2 small L3, 14 four-qubit L1 and 15 L2 cells change. |
| Installer on original Git baseline | PASS | Conflicts refused before writes, inspection is read-only, backups created, all final hashes match, repeat application is a no-op. |
| Reproduction after bundle integration | PASS | The updated baseline rebuilds omitted dictionaries, verifies all certificates and exports the expected 31 changes. |
| Fast Python regression suite | PASS: 4 tests with multiple cases | Detects both modulus mistakes, decoration-mapping errors, missing dictionaries, corrupt columns and state-convention errors. |
| Deliberately corrupted dual witness | REJECTED as expected | Adding 0.1 to a saved witness fails mathematical certificate verification. |
| Deliberately changed result after verification | REJECTED as expected | Export refuses stale verification rather than plotting unchecked results. |
| Julia regression suite | PASS: 37 assertions | State and dictionary conventions, LP optimality, saved-certificate round trip and stale-file rejection. |
| Julia legacy-helper suite | PASS: 15 assertions | Normalization, Pauli Y convention and gate-order agreement for n=2,3. |

Worst certificate errors: reconstruction 1.068e-13, dual violation 1.248e-13,
primal-dual gap 1.954e-14, reported-objective discrepancy 4.441e-16. All are
well below the 1e-8 acceptance threshold. The exact certificates are symbolic,
not merely high-precision floating-point comparisons.

Independent enumeration checks the 3/10/50 unrestricted and 2/6/24 independent
colored-graph catalogues through n=4. All 417 stored targets through n=5 are
checked from state vectors. No n=5 robustness claim is made. The 22/24 strict
four-qubit L2 improvement is verified from unrounded values.

## Evaluation of the tests

The tests do more than repeat expected numeric constants: primal witnesses
provide upper bounds, dual witnesses provide lower bounds, and they are checked
in the unreduced physical coordinates. Two independently written state builders
and the Julia implementation provide additional convention checks. Intentional
corruption is used to confirm that critical checks fail, not just pass on good
inputs. Rebuilt pools are compared as complete mathematical column sets; changes
in HDF5 serialization are reported without mistaking them for scientific changes.
The manuscript exporter verifies the hashes of the certified artifacts before
reading their results.

The first comparison of regenerated results intentionally required exact binary
float equality and detected three differences <=4.441e-16. Inspection established
that these were numerical roundoff with identical printed results and valid
primal/dual bounds. The recorded comparison uses a 1e-12 tolerance and also checks
nine-decimal agreement; no result was adjusted to hide that initial finding.

Logs are under `docs/test-results/`. Workspace-specific path prefixes in logs
are replaced with placeholders for publication; numerical outputs are unchanged.
The integration installer is separately exercised for inspection, application,
idempotence, and rejection of newer conflicting files; see its test log.

## Limits

Dependencies were already installed; a fresh network download/installation and
other OS/Python combinations were not evaluated. Broader ranges in requirements.txt
are compatibility suggestions, not a tested matrix. The Julia package manifest
is retained. Full n=4 GAP orbit regeneration was not repeated in this transfer;
the Python full-orbit checks and prior Julia n=2,3 agreement are documented in
the original numerical audit. This is not an audit of every manuscript theorem,
a legal review, or a comprehensive repository-history secret scan. The full
revised manuscript has not been compiled. No public release was pushed.
