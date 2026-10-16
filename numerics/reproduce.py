"""Regenerate corrected dictionaries and every n=2,3,4 catalogue robustness.

The publication subset has independent Hadamard decorations (2,6,24 cases).
Every LP has a saved primal and full-coordinate dual certificate. The default
is a complete regeneration; --verify-only validates saved certificates and data.
The standalone comparison fixture records the v74 reference hash; no manuscript is required.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('OMP_NUM_THREADS', '2')
import argparse
import csv
import hashlib
import itertools
import json
import platform
import shutil
from pathlib import Path
import time
import h5py
import numpy as np
import scipy
from scipy import sparse
import sympy
import small_core as sc
import four_core as fc
import exact_certificates as ec

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/manuscript_v75'
SEED = 'termwise_product_mod3_sum_mod2_v1'
CONVENTIONS = {'seed_convention':SEED,'orbit_convention':'proper_local_clifford',
               'pauli_convention':'lexicographic I,X,Y,Z; Tr(P A)',
               'state_family':'T^n H^U E(G) |+>^n'}

def read_matrix(name, n):
    with h5py.File(ROOT/f'data/{name}_{n}.h5') as f:
        a=f['data'][:].T
    assert a.shape[0]==4**n
    return a

def column_set(a):
    if sparse.issparse(a):a=a.toarray()
    assert np.all(a==np.rint(a))
    return {row.tobytes() for row in np.asarray(a.T,dtype=np.int8)}

def write_pool(path, matrix):
    """Write independently enumerated columns in Julia-compatible disk order."""
    a=matrix.toarray() if sparse.issparse(matrix) else matrix
    with h5py.File(path,'w') as f:
        f['data']=np.asarray(a.T,dtype=np.int8)
        f.attrs['pauli_convention']=CONVENTIONS['pauli_convention']
        f.attrs['trace_convention']='all identity coefficients +1'

def save_M(n, matrix, digits):
    a=matrix.toarray().astype(np.int8) if sparse.issparse(matrix) else matrix.astype(np.int8)
    path=ROOT/f'data/max_weight_operators_p_2_{n}.h5'
    seed=np.zeros(4**n,dtype=np.int8);seed[0]=1
    seed[np.all(digits!=0,axis=1)]=sc.gamma2_seed(n)
    # HDF5.jl reverses dimension order on disk; preserve the existing Julia API.
    with h5py.File(path,'w') as f:
        for k,v in CONVENTIONS.items():f.attrs[k]=v
        f['data']=a.T
        f['fiducial']=seed;f['n']=n;f['p']=2
        f['pauli_labels']=np.array([''.join('IXYZ'[d] for d in ds) for ds in digits],dtype=h5py.string_dtype())
    return path

def pool_small(n, write):
    digits,comm,phase=sc.pauli_data(n)
    levels=sc.isotropic_spaces(n,comm)
    c,_=sc.cnc_columns(n,levels,comm,phase)
    pools={'C':c,'D':sc.deterministic_columns(n),
           'S':sc.stabilizer_columns(n,levels,phase),
           'M':sc.maxweight_clifford_columns(n,digits)}
    for name,file in [('C','cnc_operators'),('D','deterministic'),('S','stabilizer_states')]:
        path=ROOT/f'data/{file}_{n}.h5'
        if write and not path.exists():
            write_pool(path,pools[name])
        stored=read_matrix(file,n)
        if name=='C' and write:
            # The historical n=3 file has 36,864 columns with trace -1.
            # Whole-column sign flips preserve the signed l1 optimum, but the
            # manuscript phase-space operators themselves must have trace +1.
            assert set(stored[0]).issubset({-1,1})
            assert column_set(pools[name])==column_set(stored*stored[[0]]),(name,n)
            path=ROOT/f'data/{file}_{n}.h5'
            backup=ROOT/'archive/pre_termwise_seed_20261012/data'/path.name
            backup.parent.mkdir(parents=True,exist_ok=True)
            if not backup.exists():shutil.copy2(path,backup)
            with h5py.File(path,'w') as f:
                f['data']=pools[name].toarray().astype(np.int8).T
                f.attrs['pauli_convention']=CONVENTIONS['pauli_convention']
                f.attrs['trace_convention']='all identity coefficients +1'
        else:
            assert column_set(pools[name])==column_set(stored),(name,n)
    if write:save_M(n,pools['M'],digits)
    assert column_set(pools['M'])==column_set(read_matrix('max_weight_operators_p_2',n))
    pools['L1']=sparse.hstack([pools['D'],pools['M']]).tocsc()
    pools['L2']=sparse.hstack([pools['S'],pools['D'],pools['M']]).tocsc()
    pools['L3']=sparse.hstack([pools['C'],pools['D'],pools['M']]).tocsc()
    return digits,comm,phase,pools

def pool_four(write):
    p=fc.build_pools()
    for name,file in [('D','deterministic'),('S','stabilizer_states')]:
        path=ROOT/f'data/{file}_4.h5'
        if write and not path.exists():write_pool(path,p[name])
        assert column_set(p[name])==column_set(read_matrix(file,4)),name
    m=np.zeros((256,len(p['M'])),dtype=np.int8);m[0]=1;m[p['full']]=p['M'].T
    if write:save_M(4,m,p['digits'])
    assert column_set(m)==column_set(read_matrix('max_weight_operators_p_2',4))
    return p

def independent(record):
    return not any(a in record['red'] and b in record['red'] for a,b in record['edges'])

def records(n):
    return json.loads((ROOT/f'data/bicolored_graphs_n{n}.json').read_text())

def case_four(record):
    return {'index':record['id'],'graph':'catalogue','edges':[[a+1,b+1] for a,b in record['edges']],
            'U':[a+1 for a in record['red']]}

def physical_target(n,edges,U):
    """Independent state-vector construction with direct Pauli matrix elements."""
    bitstrings=np.array(list(itertools.product(range(2),repeat=n)))
    ket=np.ones(2**n,dtype=complex)/np.sqrt(2**n)
    for a,b in edges:ket*=(-1.)**(bitstrings[:,a]*bitstrings[:,b])
    h=np.array([[1,1],[1,-1]])/np.sqrt(2)
    for q in U:
        gate=np.array([[1.]])
        for site in range(n):gate=np.kron(gate,h if site==q else np.eye(2))
        ket=gate@ket
    ket*=np.exp(1j*np.pi/4*bitstrings.sum(axis=1))
    ds=np.array(list(itertools.product(range(4),repeat=n)))
    x=(ds==1)|(ds==2);z=(ds==2)|(ds==3)
    mask=x@2**np.arange(n-1,-1,-1)
    coeff=(1j)**np.sum(ds==2,axis=1)
    values=coeff*np.sum(ket[np.arange(2**n)[None,:]^mask[:,None]].conj()*
                       (-1.)**(z.astype(int)@bitstrings.T)*ket[None,:],axis=1)
    assert np.max(abs(values.imag))<1e-12
    return values.real

def audit_catalogues():
    """Exhaustive labeled connected graphs, quotiented by color-preserving permutations."""
    result={}
    for n in (2,3,4):
        perms=list(itertools.permutations(range(n)))
        def canonical(edges,U):
            return min((tuple(sorted(tuple(sorted((p[a],p[b]))) for a,b in edges)),tuple(sorted(p[q] for q in U))) for p in perms)
        expected=set();pairs=list(itertools.combinations(range(n),2))
        for mask in range(2**len(pairs)):
            edges=[p for j,p in enumerate(pairs) if mask>>j&1]
            seen={0}
            while True:
                more=seen|{b for a,b in edges if a in seen}|{a for a,b in edges if b in seen}
                if more==seen:break
                seen=more
            if len(seen)!=n:continue
            for colors in range(2**n):
                U=[q for q in range(n) if colors>>q&1]
                expected.add(canonical(edges,U))
        actual=[canonical(x['edges'],x['red']) for x in records(n)]
        assert len(actual)==len(set(actual)) and set(actual)==expected
        result[str(n)]={'all_decorations':len(actual),'independent_decorations':sum(independent(x) for x in records(n))}
    assert [result[str(n)]['independent_decorations'] for n in (2,3,4)]==[2,6,24]
    # Verify the manuscript's 24 explicit graph drawings against the catalogue.
    keys=lambda edges,U:(tuple(sorted(tuple(sorted(e)) for e in edges)),tuple(sorted(U)))
    selected={keys(x['edges'],x['red']) for x in records(4) if independent(x)}
    assert selected=={keys([[a-1,b-1] for a,b in x['edges']],[a-1 for a in x['U']]) for x in fc.catalogue()}
    return result

def audit_targets():
    result={}
    for n in (2,3,4,5):
        with h5py.File(ROOT/f'data/bicolored_resource_projectors_pauli_n{n}.h5') as f:stored=f['data'][:]
        errors=[]
        for row,record in zip(stored,records(n),strict=True):
            rho=physical_target(n,record['edges'],record['red'])
            errors.append(float(np.max(abs(row-rho))))
            assert abs(rho[0]-1)<1e-12 and abs(rho@rho-2**n)<1e-11
        assert max(errors)<1e-12
        result[str(n)]={'targets_checked':len(errors),'max_coordinate_error':max(errors)}
    return result

def record_key(n,i,kind):return f'n{n}_s{i:03d}_{kind}'

def save_numeric(archive,key,result,matrix,target):
    selected=np.flatnonzero(abs(result.primal)>1e-10)
    cols=matrix[:,selected].toarray()
    weights=result.primal[selected]
    assert np.max(abs(cols@weights-target))<1e-8
    for name,value in {'seeds':cols,'weights':weights,'dual':result.dual,'rho':target,'pool_indices':selected}.items():archive[key+'_'+name]=value

def run():
    OUT.mkdir(parents=True,exist_ok=True)
    report={'schema_version':1,'conventions':CONVENTIONS,'rows':[],'pool_counts':{},
            'catalogues':audit_catalogues(),'targets':audit_targets(),
            'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sympy.__version__}}
    archive={};exact={}
    for n in (2,3):
        digits,comm,phase,pools=pool_small(n,True)
        counts={k:p.shape[1] for k,p in pools.items()};report['pool_counts'][str(n)]=counts
        for record in records(n):
            rho=physical_target(n,record['edges'],record['red'])
            assert np.max(abs(rho-sc.resource_coordinates(n,record['edges'],record['red'],digits)))<1e-12
            for kind in ('D','S','C','L1','L2','L3'):
                numeric=sc.solve_robustness(pools[kind],rho)
                key=record_key(n,record['id'],kind);save_numeric(archive,key,numeric,pools[kind],rho)
                report['rows'].append({'n':n,'record_id':record['id'],'edges':record['edges'],'U':record['red'],
                    'independent_U':independent(record),'kind':kind,'certificate_key':key,**numeric.summary()})
                if n==3 and not record['red'] and kind in ('C','L3'):
                    graph='P3' if len(record['edges'])==2 else 'K3'
                    case={'n':n,'graph':graph,'edges':record['edges'],'U':[]}
                    target=sc.exact_resource_coordinates(n,record['edges'],[],digits)
                    exact[graph+'_'+kind]=ec.exact_optimality_certificate(case,pools[kind],numeric,target,counts,pool=kind)
            print('Completed',n,record['id'],flush=True)
    p=pool_four(True);report['pool_counts']['4']={'D':p['D'].shape[1],'S':p['S'].shape[1],'M':len(p['M'])}
    for record in records(4):
        case=case_four(record);rho=physical_target(4,record['edges'],record['red'])
        assert np.max(abs(rho-fc.resource(case,p['digits'])))<1e-12
        W,generators=fc.twirl_basis(case,p['digits']);reduced=fc.reduced_pools(p,W,['clifford'])
        for kind in ('D','S','L1','L2'):
            result,cert=fc.solve(p,reduced,W,rho,kind,'clifford')
            key=record_key(4,record['id'],kind)
            for name,value in cert.items():archive[key+'_'+name]=value
            report['rows'].append({'n':4,'record_id':record['id'],'edges':record['edges'],'U':record['red'],
                'independent_U':independent(record),'kind':kind,'certificate_key':key,**result})
        print('Completed 4',record['id'],flush=True)
    report['input_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT/'data').glob('*.h5')) if any(s in p.name for s in ('max_weight','stabilizer_states','cnc_operators','deterministic_'))}
    report['manuscript_reference_sha256']=json.loads((ROOT/'numerics/reference-v74.json').read_text())['manuscript_sha256']
    (OUT/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    (OUT/'exact_small.json').write_text(json.dumps(exact,indent=2)+'\n')
    np.savez_compressed(OUT/'certificates.npz',**archive)
    fields=('n','record_id','independent_U','kind','value','primal_l1','primal_residual','dual_violation','certificate_key')
    with (OUT/'results.csv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(report['rows'])
    verify()

def verify():
    report=json.loads((OUT/'results.json').read_text());archive=np.load(OUT/'certificates.npz',allow_pickle=False)
    exact=json.loads((OUT/'exact_small.json').read_text());checks=[]
    # Historical byte hashes record the original serialization. Rebuilt HDF5
    # files can encode identical column sets with different bytes/metadata.
    # pool_small/pool_four independently check the complete mathematical pools.
    serialization_changes=[rel for rel,sha in report['input_sha256'].items()
        if not (ROOT/rel).exists() or hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=sha]
    for n in (2,3,4):
        if n<4:digits,comm,phase,pools=pool_small(n,False)
        else:pools=pool_four(False)
        for item in report['rows']:
            if item['n']!=n:continue
            k=item['certificate_key'];kind=item['kind']
            rho=physical_target(n,item['edges'],item['U'])
            seeds=archive[k+'_seeds'];weights=archive[k+'_weights'];dual=archive[k+'_dual']
            assert np.max(abs(archive[k+'_rho']-rho))<1e-12
            if n<4:
                expected=pools[kind][:,archive[k+'_pool_indices']].toarray()
                assert np.array_equal(seeds,expected)
                reconstructed=seeds@weights
                violation=max(0,float(np.max(abs(pools[kind].T@dual)))-1)
            else:
                families=archive[k+'_families'];indices=archive[k+'_pool_indices']
                allowed={'D':{0},'S':{1},'L1':{0,2},'L2':{0,1,2}}[kind]
                assert set(map(int,families))<=allowed
                assert np.array_equal(seeds,fc.original_seeds(pools,families,indices))
                case={'edges':[[a+1,b+1] for a,b in item['edges']],'U':[a+1 for a in item['U']]}
                W,generators=fc.twirl_basis(case,pools['digits'])
                assert np.array_equal(W,archive[k+'_W'])
                reconstructed=seeds@weights
                for ix,sign in generators:reconstructed=(reconstructed+sign*reconstructed[ix])/2
                violation=fc.full_dual_violation(pools,dual,kind,'clifford')
            metrics={'key':k,'residual':float(np.max(abs(reconstructed-rho))),
                'dual_violation':violation,'gap':float(abs(np.sum(abs(weights))-rho@dual)),
                'reported_error':float(abs(item['value']-rho@dual))}
            assert max(v for k,v in metrics.items() if k!='key')<1e-8,metrics
            checks.append(metrics)
        if n==3:
            members=column_set(pools['M'])
            for cert in exact.values():ec.verify_one(cert,digits,comm,phase,pools,members)
    # Nestedness and the reported 22/24 improvement are checked from unrounded values.
    groups={}
    for row in report['rows']:groups.setdefault((row['n'],row['record_id']),{})[row['kind']]=row['value']
    for key,v in groups.items():
        assert v['L2']<=min(v['S'],v['L1'])+1e-8 and v['L1']<=v['D']+1e-8
        if 'L3' in v:assert v['L3']<=min(v['L2'],v['C'])+1e-8 and v['C']<=v['S']+1e-8
    chosen=[x for x in records(4) if independent(x)]
    improved=[x['id'] for x in chosen if groups[4,x['id']]['L2']<groups[4,x['id']]['S']-1e-8]
    assert len(improved)==22
    result={'passed':True,'numerical_certificates':len(checks),'exact_small_certificates':len(exact),
        'four_independent_improvements':len(improved),'serialization_changes':serialization_changes,'worst_errors':{k:max(c[k] for c in checks) for k in ('residual','dual_violation','gap','reported_error')},'checks':checks}
    result['verified_artifacts_sha256']={name:hashlib.sha256((OUT/name).read_bytes()).hexdigest() for name in ('results.json','certificates.npz','exact_small.json')}
    (OUT/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k!='checks'},flush=True)

if __name__=='__main__':
    if not __debug__:raise RuntimeError('Do not use Python -O: assertions perform certificate checks.')
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--verify-only',action='store_true',help='Verify saved certificates; run --prepare-data first in a fresh checkout.')
    mode.add_argument('--prepare-data',action='store_true',help='Rebuild dictionaries and check saved targets/catalogues without solving LPs.')
    args=parser.parse_args()
    if args.prepare_data:
        OUT.mkdir(parents=True,exist_ok=True)
        for n in (2,3):pool_small(n,True)
        pool_four(True)
        audit={'catalogues':audit_catalogues(),'targets':audit_targets(),'passed':True}
        (OUT/'preparation.json').write_text(json.dumps(audit,indent=2)+'\n')
        print('Prepared independently enumerated dictionaries; target and catalogue checks passed.',flush=True)
    elif args.verify_only:verify()
    else:run()
