using GAP
const g = GAP.Globals
const gjl = GAP.gap_to_julia
const jlg = GAP.julia_to_gap

include("pauli.jl");

# convert from default GAP convention to physics convention for symplectic form:
function form_conversion(a,n)
    b = [a[i] for i in 1:n];
    for i in 1:n; push!(b,a[2*n-i+1]);; end
    return b
end

"""
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
"""

mutable struct SympPerm
    """
    Description:

    Group:
    fdict:
    bdict:
    """
    Group::GapObj
    fdict::Dict{Vector{Int},Int}
    bdict::Dict{Int,Vector{Int}}

    function SympPerm(n::Int)
        # Generate symplectic group:
        spn = g.SymplecticGroup(2*n,2)

        # Generate domain of symplectic group action:
        En = g.ShallowCopy(g.Elements(g.GF(2)^(2*n))); g.Sort(En);
        hom = g.ActionHomomorphism(spn,En)

        Fn = sort(all_dit_strings(2,2*n));
        Dn = [i for i in 1:length(Fn)];

        Group = hom(spn)
        fdict = Dict(zip(Fn,Dn))
        bdict = Dict(zip(Dn,Fn))

        new(Group,fdict,bdict)
    end
end



mutable struct SympOrbit
    """
    Description:
    Set:
    Orbit:
    """
    Subset::Set{Vector{Int}}
    Orbit::Set{Set{Vector{Int}}}

    function SympOrbit(n::Int,SP::SympPerm,Subset::Set{Vector{Int}})
        Subset = [form_conversion(a,n) for a in Subset];
        G = SP.Group; fdict = SP.fdict; bdict = SP.bdict

        # generate gap object of subset:
        DomS = g.Set(jlg([g.Set(jlg([jlg(fdict[a]) for a in Subset]))]))

        # generate orbits:
        DomOrbs = gjl(g.Orbits(G,DomS,g.OnSets)[1]); Orbs = Vector{Set{Vector{Int}}}([]);
        
        for DomOrb in DomOrbs
            Orb = [ form_conversion(bdict[d],n) for d in DomOrb];
            push!(Orbs,Set(Orb))
        end
        
        Orbit = Set(Orbs);
        Subset = Set(Subset)

        new(Subset,Orbit)
    end
end




"""
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
"""