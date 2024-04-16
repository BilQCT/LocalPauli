using GAP
const g = GAP.Globals
const gjl = GAP.gap_to_julia
const jlg = GAP.julia_to_gap

#include("stab.jl");

# convert from default GAP convention to physics convention for symplectic form:
function form_conversion(a,n)
    b = [a[i] for i in 1:n];
    for i in 1:n; push!(b,a[2*n-i+1]);; end
    return b
end

# Generate symplectic group as permutation group:
function symplectic_perm_group(n)
    # Generate symplectic group:
    spn = g.SymplecticGroup(2*n,2)

    # Generate domain of symplectic group action:
    En = g.ShallowCopy(g.Elements(g.GF(2)^(2*n))); g.Sort(En);
    hom = g.ActionHomomorphism(spn,En)

    Fn = sort(all_dit_strings(2,2*n));
    Dn = [i for i in 1:length(Fn)];
    fdict = Dict(zip(Fn,Dn)); bdict = Dict(zip(Dn,Fn));
    
    return [hom(spn),fdict,bdict]
end

# Generate symplectic orbit:
function symplectic_orbit(n,spn,subset,fdict,bdict)
    subset = [form_conversion(a,n) for a in subset];
    # generate gap object of subset:
    dsubset = g.Set(jlg([g.Set(jlg([jlg(fdict[a]) for a in subset]))]))

    # generate orbits:
    dorbs = gjl(g.Orbits(spn,dsubset,g.OnSets)[1]); orbs = [];
    
    for dorb in dorbs
        orb = [ form_conversion(bdict[d],n) for d in dorb]
        push!(orbs,orb)
    end
    

    return orbs
end