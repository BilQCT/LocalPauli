"""Export a plain-text v75 revision guide and LaTeX values from verified results.

Original displayed values come from reference-v74.json. Indices are matched by graph and U,
never by assuming that stored catalogue IDs are publication indices.
"""
import csv
import hashlib
import itertools
import json
from pathlib import Path
from four_core import catalogue

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/manuscript_v75/manuscript'

def canonical(n,edges,U):
    return min((tuple(sorted(tuple(sorted((p[a],p[b]))) for a,b in edges)),tuple(sorted(p[i] for i in U)))
               for p in itertools.permutations(range(n)))

def main():
    reference=json.loads((ROOT/'numerics/reference-v74.json').read_text())
    report=json.loads((ROOT/'results/manuscript_v75/results.json').read_text())
    verify=json.loads((ROOT/'results/manuscript_v75/verification.json').read_text())
    assert verify['passed'] and verify['numerical_certificates']==278
    for name in ('results.json','certificates.npz','exact_small.json'):
        assert hashlib.sha256((ROOT/'results/manuscript_v75'/name).read_bytes()).hexdigest()==verify['verified_artifacts_sha256'][name], 'Run --verify-only again: '+name+' changed after verification'
    values={}
    for r in report['rows']:values[r['n'],canonical(r['n'],r['edges'],r['U']),r['kind']]=r['value']
    def lookup(n,edges,U,kind):return values[n,canonical(n,edges,U),kind]
    small=[(2,'P_2',[(0,1)],[]),(2,'P_2',[(0,1)],[0]),
           (3,'P_3',[(0,1),(1,2)],[]),(3,'P_3',[(0,1),(1,2)],[0]),
           (3,'P_3',[(0,1),(1,2)],[1]),(3,'P_3',[(0,1),(1,2)],[0,2]),
           (3,'K_3',[(0,1),(0,2),(1,2)],[]),(3,'K_3',[(0,1),(0,2),(1,2)],[0])]
    small_rows=[];comparisons=[]
    for index,(n,g,edges,U) in enumerate(small):
        decoration=r'\varnothing' if not U else r'\{'+','.join(str(i+1) for i in U)+r'\}'
        v=[lookup(n,edges,U,k) for k in ('S','C','D','L3')]
        small_rows.append(f'{n} & ${g}$ & ${decoration}$ & '+' & '.join(f'{x:.3f}' for x in v)+r' \\')
        for kind,old,new in zip(('S','C','D','L3'),reference['small'][index],v):
            comparisons.append({'table':'small','state':f'{g}, U={U}','kind':kind,'old':old,'new':new,'changed':abs(old-new)>0.0005})
    coords={k:[] for k in ('L1','L2','S')};four_rows=[]
    old_values=reference['four']
    for c in catalogue():
        edges=[[a-1,b-1] for a,b in c['edges']];U=[x-1 for x in c['U']]
        v={k:lookup(4,edges,U,k) for k in coords}
        four_rows.append({'index':c['index'],'graph':c['graph'],'U_one_based':c['U'],**v})
        for k in coords:
            coords[k].append(f'({c["index"]},{v[k]:.6f})')
            old=old_values[k][c['index']-1]
            comparisons.append({'table':'four','state':str(c['index']),'kind':k,'old':old,'new':v[k],'changed':abs(old-v[k])>0.0000005})
    plots=[]
    for k,style,label in [('L1','blue!55,draw=blue','L_1'),('L2','red!55,draw=red','L_2'),('S','orange!65,draw=brown','S')]:
        plots.extend([r'\addplot+[fill='+style+r'] coordinates {'+' '.join(coords[k])+'};',r'\addlegendentry{$\mathfrak R_{'+label+'}$}'])
    small_tex='\n'.join(small_rows)
    plot_tex='\n'.join(plots)
    definition=r'''Let $D$, $S$, and $C$ denote the deterministic operators, all stabilizer
projectors, and the CNC operators, respectively. For the numerical comparisons,
$M$ is the complete orbit of the $p=2$ max-weight seed under tensor products
of single-qubit Clifford conjugations. The seed has identity coefficient $1$,
zero coefficients at weights $1,\ldots,n-1$, and full-weight coefficients
$(-1)^{\gamma^{(2)}(u)}$, with $X,Y,Z$ identified with $u_i=0,1,2$.
Each product in $\gamma^{(2)}$ is evaluated modulo $3$, and the resulting
terms are added modulo $2$, as in Eq.~\eqref{eq:gamma p efficient}.
This gives $|M|=96$, $13{,}824$, and $165{,}888$ for $n=2,3,4$,
respectively. No additional seeds, tensor-product families, or independent
Bell outcome relabelings are included in $M$. The local phase spaces used
in our analysis are'''
    exact=r'''The two nontrivial three-qubit values are
\[
\mathfrak R_{L_3}(\rho_{P_3,\varnothing})
=\frac{5+3\sqrt2}{8}\approx1.155330086,
\qquad
\mathfrak R_{L_3}(\rho_{K_3,\varnothing})
=\frac12+\frac{9\sqrt2}{20}\approx1.136396103.
\]
Both are strictly below the CNC robustness
$1+\sqrt2/5\approx1.282842712$ of the product magic state.
The three-qubit values are certified by exact primal decompositions and
dual witnesses over $\mathbb Q(\sqrt2)$.'''
    method=r'''All computations use the expectation coordinates
$a_{Pj}=\operatorname{Tr}(PA_j)$ and $r_P=\operatorname{Tr}(P\rho)$,
in lexicographic $I,X,Y,Z$ order, and solve
\[
\min_{Aq=r}\|q\|_1
=\max_{|A^{\mathsf T}y|\leq\mathbf1}r^{\mathsf T}y.
\]
Every dictionary operator has identity coordinate $1$, so the identity
constraint enforces $\sum_jq_j=1$. Numerical primal reconstructions,
dual constraints, and primal--dual gaps were checked in the full Pauli
coordinates to tolerance $10^{-8}$. The four-qubit computation uses a
lossless average over local Clifford symmetries of each target; the
certificates were checked against every original dictionary operator.
These finite-dictionary values are upper bounds on the robustness over
the full efficiently representable locally closed phase space.'''
    gamma=r'''\gamma^{(p)}(u)
:=\left[\sum_{1\leq i_1<\cdots<i_p\leq n}
\left[u_{i_1}\cdots u_{i_p}\right]_3\right]_2,
\qquad [x]_m\in\{0,\ldots,m-1\}.'''
    changes=sum(r['changed'] for r in comparisons)
    counts={f'{table}/{k}':sum(r['changed'] for r in comparisons if r['table']==table and r['kind']==k)
            for table,ks in [('small',('S','C','D','L3')),('four',('L1','L2','S'))] for k in ks}
    assert len(comparisons)==104 and changes==31 and counts['four/L1']==14 and counts['four/L2']==15
    overview='\n'.join(f'{r["index"]:2d} {r["graph"]:3s} {str(r["U_one_based"]):10s} {r["L1"]:.9f} {r["L2"]:.9f} {r["S"]:.9f}' for r in four_rows)
    text=f'''V75 MANUSCRIPT REVISION GUIDE -- VERIFIED NUMERICAL CORRECTIONS
Reference: v74/main.tex supplied on 12 October 2026.
The reference manuscript has NOT been edited. Line numbers below refer to v74.

DECISION AND SCOPE
Retain the manuscript's seed: each monomial is multiplied in Z3; the terms
are summed in Z2. Use the complete proper local-Clifford orbit, preserving the
original computational choice of symmetry group. Hoban et al. provide the
general mixed-arithmetic framework and vertex criterion; they do not uniquely
specify this particular elementary-polynomial seed.

Of the 104 displayed robustness values, 31 change: two small-table L3 cells,
14 four-qubit L1 bars, and 15 four-qubit L2 bars. All other displayed values
agree at printed precision. Neither the graph drawings nor state numbering
changes. The 22/24 four-qubit L2 advantage remains true.

1. MAKE THE SEED UNAMBIGUOUS (near lines 1951-1983)
Keep the label eq:gamma p efficient. Replace the expression inside its
equation environment by the following explicit nested-modulus notation:

{gamma}

The existing sentence that each product is evaluated in Z3 and the final sum
is reduced modulo 2 is correct. The p=1 statement and incremental modulo-2
summation explanation remain correct. For the quadratic example use
  gamma^(2)(u) = [sum_(i<j) [u_i u_j]_3]_2.
The two-variable truth table remains unchanged. Do not insert a modulo-3
reduction around the whole sum.

1a. CORRECT THE PARITY CLAIM IMMEDIATELY AFTER THE TRUTH TABLE
The sentence claiming that the total truth-table parity is invariant under
party, input AND outcome relabelings is false. The displayed O5 seed has
two negative correlations, hence even parity. Flipping the outcome of the
first setting of the first party negates one entire row (three correlations),
giving five negative correlations and odd parity, within the SAME Bell orbit.
Thus parity cannot distinguish O4 from O5. Keep the displayed seed's parity
calculation if desired, but replace the subsequent two sentences with:

"This parity refers to the displayed representative and is not invariant
under outcome relabelings. The orbits $\\oO_4$ and $\\oO_5$ are distinct
max-weight orbits under the relabeling group used in
Table~\\ref{{tab:ns232_orbits_verbose}}."

This correction does not change either orbit size or any robustness value.

2. DEFINE THE NUMERICAL M (paragraph preceding L0, L1, L2, L3)
Replace the sentence beginning "Let D, S, C, and M denote ..." through
"The local phase spaces used in our analysis are" with:

{definition}

Keep the following displayed definitions and robustness inequalities.
S means all stabilizer projectors, including entangled stabilizers.
The Bell orbit is a different choice: its two-qubit size is 192 (the O5
table size), whereas the LC orbit used for robustness has size 96. These
counts describe different groups and are not a contradiction.

3. REPLACE THE SMALL-TABLE BODY (lines 2432-2441)
Keep the table header, labels, caption and graph catalogue. Insert these
eight rows, retaining horizontal rules between graph blocks if desired:

{small_tex}

Thus P3 changes 1.183 -> 1.155, and K3 changes 1.190 -> 1.136.
All S, C and L0 values, and the remaining six L3 values, are unchanged.

4. UPDATE THE REPEATED NUMERICAL PARAGRAPH (lines 2537-2540)
Replace the sentence beginning "Numerically, the L3 robustness values ..."
and its CNC comparison with the following (keep the subsequent strictness
conclusion):

{exact}

5. REPLACE ALL THREE FOUR-QUBIT PLOT SERIES (lines 2635-2640)
Keep the axis, caption, figure label, graph catalogue and numbering. Replace
the addplot/addlegendentry lines with this block. The S series is supplied
for completeness and has not changed:

{plot_tex}

Full values for checking (U uses the manuscript's one-based vertex labels):
ID Graph U          L1          L2          S
{overview}

6. ADD A REPRODUCIBILITY PARAGRAPH
Place after the phase-space definitions or near the numerical figure:

{method}

7. CLAIMS TO KEEP / SCOPE TO MAKE EXPLICIT
- The strict P3 and K3 separation from CNC survives and is stronger numerically.
- L2 improves over S in 22 of 24 four-qubit cases. States 5 and 10 remain
  stabilizer states with S=L2=1.
- The undecorated comparisons to product magic states remain unchanged in logic.
- "Exhaustive" means connected graphs with INDEPENDENT Hadamard decorations,
  up to color-preserving graph isomorphism: 2, 6, 24 cases for n=2,3,4.
  Unrestricted decorations instead give 3, 10, 50 cases; do not mix the counts.
- The robustness values quantify the stated additive-error estimation bounds;
  they do not alone prove a total-runtime improvement or an exact sampler
  for signed decompositions. Existing caveats on initialization remain relevant.
- The numerical audit does not constitute a new proof of the full simulation
  theorems. In particular it does not verify the unspecific sentence about
  "initial numerical results" at n=3 in the conclusion. Retain that sentence
  only if its separate evidence can be identified; otherwise remove it.

8. VERIFIED EVIDENCE AND FILES
- 278 newly solved numerical LPs, all with full-coordinate primal/dual checks.
- 4 newly verified exact certificates: P3/K3 in C and L3.
- All 417 saved target states through n=5 agree with independent state-vector
  construction (no claim that all n=5 robustness values were recomputed).
- Complete connected colored-graph catalogues checked independently through n=4.
- The original Julia orbit generator agrees with Python on all n=2 and n=3
  columns after the fix. Julia tests also validate n=4 seed and dictionary size.
- The two-qubit nonsignaling polytope was re-enumerated: 1,408 vertices in
  orbits of sizes 64, 576, 288, 288, 192, as in the manuscript.
- The two-qubit Pauli polytope was independently enumerated: 22,320 vertices,
  agreeing with the conclusion. The facet counts 36 and 60 agree with the
  corresponding local and full stabilizer dictionaries.
- docs/numerical-audit.md records the original audit; docs/v75-test-evaluation.md records standalone release testing.

Corrected data/certificates:
  LocalPauli/results/manuscript_v75/results.json
  LocalPauli/results/manuscript_v75/results.csv
  LocalPauli/results/manuscript_v75/certificates.npz
  LocalPauli/results/manuscript_v75/exact_small.json
  LocalPauli/results/manuscript_v75/verification.json
Machine-readable original-to-corrected comparison:
  results/manuscript_v75/manuscript/manuscript-value-changes.csv
Copyable LaTeX plot block:
  results/manuscript_v75/manuscript/v75-four-qubit-plots.tex
Copyable LaTeX small-table rows:
  results/manuscript_v75/manuscript/v75-small-table-rows.tex

The guide changes numerical definitions and claims only where needed. It does
not silently change the manuscript to preserve the old numerical values.
'''
    OUT.mkdir(exist_ok=True)
    (OUT/'v75-manuscript-changes.txt').write_text(text)
    (OUT/'v75-small-table-rows.tex').write_text(small_tex+'\n')
    (OUT/'v75-four-qubit-plots.tex').write_text(plot_tex+'\n')
    with (OUT/'manuscript-value-changes.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(comparisons[0]));writer.writeheader();writer.writerows(comparisons)
    # Replace the formerly stale repository table artifact with generated data.
    (ROOT/'robustness_tables.tex').write_text(r'''% Generated from verified results/manuscript_v75/results.json.
% Publication subset: independent Hadamard decorations, proper LC orbit M.
\begin{table}[t]
\centering
\begin{tabular}{cclrrrr}
\hline
$n$ & $G$ & $U$ & $\mathfrak R_S$ & $\mathfrak R_C$ & $\mathfrak R_{L_0}$ & $\mathfrak R_{L_3}$ \\
\hline
'''+small_tex+r'''
\hline
\end{tabular}
\caption{Corrected robustness for the independent-decoration two- and three-qubit catalogue.}
\end{table}

\begin{figure}[t]
\centering
\begin{tikzpicture}
\begin{axis}[ybar=0pt,bar width=4pt,width=\textwidth,height=7.2cm,
xmin=0.5,xmax=24.5,ymin=0,ymax=3.1,xtick={1,...,24},
xlabel={State Index},ylabel={Robustness},ymajorgrids=true,
legend columns=3,legend style={at={(0.5,1.02)},anchor=south,draw=none}]
'''+plot_tex+r'''
\end{axis}
\end{tikzpicture}
\caption{Corrected robustness for the 24 independent-decoration four-qubit states in the manuscript order.}
\end{figure}
''')
    print('Exported guide; comparison counts:',counts)

if __name__=='__main__':main()
