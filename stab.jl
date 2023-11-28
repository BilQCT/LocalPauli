# Stabilizer theory

function tensor_prod(A,n)
    An = A;
    for i in 1:(n-1); An = kron(An,A); end;
    return An
end

################################################################################################

function pauli(x)
    if typeof(length(x)//2) != Rational{Int64}
        return "Must be bit-string of even length."
    elseif typeof(length(x)//2) == Rational{Int64}
        # number of qubits:
        n = Int(length(x)//2);
        
        # Pauli Matrices:
        X = Matrix([0 1; 1 0]);
        Z = Matrix([1 0; 0 -1]);
        
        # x,z bit strings:
        a = x[1:n]; b = x[n+1:end];
        
        # initialize Pauli:
        if n == 1;
            XaZb = (X^a[1]*Z^b[1])
        else
            XaZb = (X^a[1]*Z^b[1]);
            for i in 2:n
                XaZb = kron(XaZb,(X^a[i]*Z^b[i]));
            end
        end
        
        # multiply by phase:
        arg_phase = (transpose(a)*b) % 4;
        phase = im^(arg_phase)
        
        return phase*XaZb
    end
end



################################################################################################

function sum_bitstrings(x,y)
    if length(x) != length(y)
        return "Bit-strings must be same length"
    else
        n = length(x);
        return [(x[i]+y[i])%2 for i in 1:n]
    end
end


################################################################################################


function omega(a,b)
    n = Int(length(a)/2);
    ax = a[1:n]; az =  a[n+1:end];
    bx = b[1:n]; bz =  b[n+1:end];
    
    # compute symplectic inner product:
    omega_ab = ((transpose(ax)*bz)+ (transpose(az)*bx)%2);
    
    return omega_ab
end


################################################################################################


function commuting_set(S)
    for a in S
        for b in S
            if omega(a,b) == 1;
                return false
            end
        end
    end
    return true
end



################################################################################################




function beta(a,b)
    if omega(a,b) == 1
        return "Must be commuting Paulis"
    else
        n = Int(length(a)/2);
        Ta = pauli(a); Tb = pauli(b); Tab = pauli(sum_bitstrings(a,b));
        Id_n = tensor_prod(I(2),n);
        
        # check beta:
        if Ta*Tb*Tab == Id_n
            return 0
        else
            return 1
        end
    end
end



################################################################################################



function canonical_generators(n)
    # initialize array:
    gens = [];
    
    for i in 1:n
        g = [0 for j in 1:(2*n)]; g[i] = 1;
        push!(gens,g);
    end
    
    return gens
end



################################################################################################




function canonical_maximal_isotropic(n)
    gens = canonical_generators(n);
    return isotropic(gens)
end


################################################################################################



function generate_elements(N,n,gens,iso)
    # if recursion parameter less than N:
    if n < N
        # number of old elements:
        m = length(iso); 
        
        # generator to add by:
        a = gens[n]
        
        # generate new elements:
        iso_prime = [sum_bitstrings(iso[i],a) for i in 1:m]
        
        # all current elements:
        iso = vcat(iso,iso_prime)
        
        # update recursion parameter:
        n += 1; generate_elements(N,n,gens,iso);
    else
        return collect(Set(iso))
    end   
end


################################################################################################



# input: array of 2n-bit strings; i.e., generators:
# output: all elements of isotropic subspace
function isotropic(gens)
    # remove possible redundant elements:
    gens = collect(Set(gens))

    if length(gens) == 1
        L = length(gens[1])
        Id = [0 for i in 1:L];
        return [Id,gens[1]];
    end
    
    # n: number of qubits; k: gens of stabilizer subspace:
    n = Int(length(gens[1])/2); k = length(gens);
    
    # initialize isotropic subspace:
    iso = copy(gens);
    
    # check if all elements commute:
    for i in 1:k
        a = gens[i];
        for j in 1:k
            if j != i
                b = gens[j];
                if omega(a,b) == 1
                    error("All elements must commute")
                end
            end
        end
    end
    
    # initialize recursion parameter:
    n = 1;
    
    # generate elements:
    return generate_elements(k,n,gens,iso)
end



################################################################################################


function gens_to_check_matrix(gens)
    m = length(gens); k = length(gens[1]);
    M = Array{Int64}(undef, 0, k);
    for i in 1:m
        gi = Vector{Int64}(gens[i]);
        M = vcat(M,transpose(gi));
    end
    return M
end



################################################################################################


function E(n)
    return all_dit_strings(2,2*n)
end


function bit_to_pauli_string(x)
    # raise error if not of correct length:
    if (length(x) % 2) != 0
        error("must have even length.");
    end;
    
    # number of qubits:
    n = Int64(length(x)/2);
    
    # initialize pauli string:
    s = "";
    
    for i in 1:n
        if (x[i] == 0) && (x[i+n] == 0)
            s = s*"i";
        elseif (x[i] == 1) && (x[i+n] == 0)
            s = s*"x";
        elseif (x[i] == 0) && (x[i+n] == 1)
            s = s*"z";
        elseif (x[i] == 1) && (x[i+n] == 1)
            s = s*"y";
        else
            error("not valid input");
        end;
    end;
    return s
end


function pauli_list(n)
    En = E(n);
    
    pn = [];
    for a in En
        push!(pn,bit_to_pauli_string(a));
    end
    
    perm = sortperm(pn);
    
    return En[perm], pn[perm];
end


################################################################################################


using Combinatorics

using Nemo

function find_gens(isotropic)
    # determine ring:
    Z2 = ResidueRing(ZZ, 2);
    
    # determine number of qubits:
    n = Int64(log2(length(isotropic)));
    
    # remove identity
    identity = [0 for i in 1:(2*n)];
    idx = findall(x->x != identity, isotropic)
    
    combs = combinations(isotropic[idx],n);
    for c in combs
        
        A = matrix(Z2, gens_to_check_matrix(c));
        r = rank(A);
        if r == n; return c; end;
    end
    
end


################################################################################################

using Combinatorics

function Gamma(gens)
    # generate isotropic subspace:
    N = length(gens);
    
    # all functions: gens -> Z2:
    Z = all_dit_strings(2,N);
    
    Gamma = [];
    
    for z in Z
        
        # initialize value assignment:
        gamma = Dict(zip(gens,z));
        
        # fix image of identity:
        identity = [0 for i in 1:length(gens[1])];
        gamma[identity] = 0;
        
        # generate all linear combinations of generators:
        for i in 2:N
            combs = collect(combinations(gens,i));

            for j in 1:length(combs)
                
                a = combs[j][1];
                s = gamma[a];
                
                for c in combs[j][2:end]
                    bta = beta(a,c)
                    a = sum_bitstrings(a,c);
                    s = (s + gamma[c]+bta) %2;
                end
                
                pair = Pair(a,s);
                if pair ∉ gamma; gamma[a] = s; end;
            end
        end
        push!(Gamma,gamma);
    end
    
    return Gamma;
end


######################################


# recursive loop to construct all dit strings:
function generate_dit_strings(d,N,n,s)
    if n < N
        Zd = [(i-1) for i in 1:d] ; dit_strings = [];
        for i in 1:length(s)
            for x in Zd
                #println(i)
                push!(dit_strings,push!(copy(s[i]),x));
            end
        end
        s = dit_strings; n = n+1;
        generate_dit_strings(d,N,n,s)
    else
        return s
    end
end

######################################

# input: outcomes (d) and number of generators (N) (N=dim(simplex))
function all_dit_strings(d,N)
    s = [[i-1] for i in 1:d]; n = 1;
    return generate_dit_strings(d,N,n,s)
end


################################################################################################


# Generating stabilizer states:

function stabilizer_coefficients(n,II)
    # generate pauli strings:
    En, pn = pauli_list(n);
    
    #println(En,"\n",pn,"\n")
    
    # number of paulis:
    N = length(En);
    
    # index list:
    idx = [i for i in 1:N];
    
    # mapping: En to idx:
    dict = Dict(zip(En,idx));
    
    #println(dict)
    
    A = Array{Int64}(undef, 0, N);
    
    for iso in II
        # generate valid value assignments
        gens = find_gens(iso);
        Gamma_iso = Gamma(gens);
    
        for gamma in Gamma_iso
            # initialize matrix
            a = zeros(Int64,1,N);
            
            # index placement in a:
            idx_iso = [dict[a] for a in collect(keys(gamma))];

            # coefficients:
            coefficients = [(-1)^s for s in collect(values(gamma))];

            a[idx_iso] .= coefficients;

            A = vcat(A,a);
        end
    end
    return A;
end