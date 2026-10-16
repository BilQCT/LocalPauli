"""Export tables and plot data from verified robustness results."""
import csv
import hashlib
import itertools
import json
from pathlib import Path
from four_core import catalogue

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/robustness/exports'

def canonical(n,edges,U):
    return min((tuple(sorted(tuple(sorted((p[a],p[b]))) for a,b in edges)),tuple(sorted(p[i] for i in U)))
               for p in itertools.permutations(range(n)))

def main():
    report=json.loads((ROOT/'results/robustness/results.json').read_text())
    verify=json.loads((ROOT/'results/robustness/verification.json').read_text())
    assert verify['passed'] and verify['numerical_certificates']==278
    for name in ('results.json','certificates.npz','exact_small.json'):
        assert hashlib.sha256((ROOT/'results/robustness'/name).read_bytes()).hexdigest()==verify['verified_artifacts_sha256'][name], 'Run --verify-only again: '+name+' changed after verification'
    values={}
    for r in report['rows']:values[r['n'],canonical(r['n'],r['edges'],r['U']),r['kind']]=r['value']
    def lookup(n,edges,U,kind):return values[n,canonical(n,edges,U),kind]
    small=[(2,'P_2',[(0,1)],[]),(2,'P_2',[(0,1)],[0]),
           (3,'P_3',[(0,1),(1,2)],[]),(3,'P_3',[(0,1),(1,2)],[0]),
           (3,'P_3',[(0,1),(1,2)],[1]),(3,'P_3',[(0,1),(1,2)],[0,2]),
           (3,'K_3',[(0,1),(0,2),(1,2)],[]),(3,'K_3',[(0,1),(0,2),(1,2)],[0])]
    small_rows=[]
    for index,(n,g,edges,U) in enumerate(small):
        decoration=r'\varnothing' if not U else r'\{'+','.join(str(i+1) for i in U)+r'\}'
        v=[lookup(n,edges,U,k) for k in ('S','C','D','L3')]
        small_rows.append(f'{n} & ${g}$ & ${decoration}$ & '+' & '.join(f'{x:.3f}' for x in v)+r' \\')
    coords={k:[] for k in ('L1','L2','S')};four_rows=[]
    for c in catalogue():
        edges=[[a-1,b-1] for a,b in c['edges']];U=[x-1 for x in c['U']]
        v={k:lookup(4,edges,U,k) for k in coords}
        four_rows.append({'index':c['index'],'graph':c['graph'],'U_one_based':c['U'],**v})
        for k in coords:
            coords[k].append(f'({c["index"]},{v[k]:.6f})')
    plots=[]
    for k,style,label in [('L1','blue!55,draw=blue','L_1'),('L2','red!55,draw=red','L_2'),('S','orange!65,draw=brown','S')]:
        plots.extend([r'\addplot+[fill='+style+r'] coordinates {'+' '.join(coords[k])+'};',r'\addlegendentry{$\mathfrak R_{'+label+'}$}'])
    small_tex='\n'.join(small_rows)
    plot_tex='\n'.join(plots)
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'small-table-rows.tex').write_text(small_tex+'\n')
    (OUT/'four-qubit-plots.tex').write_text(plot_tex+'\n')
    with (OUT/'four-qubit-values.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(four_rows[0]))
        writer.writeheader();writer.writerows(four_rows)
    # Replace the formerly stale repository table artifact with generated data.
    (ROOT/'robustness_tables.tex').write_text(r'''% Generated from verified results/robustness/results.json.
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
\caption{Robustness for the independent-decoration two- and three-qubit catalogue.}
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
\caption{Robustness for the 24 independent-decoration four-qubit states in the catalogue order.}
\end{figure}
''')
    print('Exported eight small-system rows and 24 four-qubit rows to',OUT)

if __name__=='__main__':main()
