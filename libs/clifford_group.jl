using GAP

# read gap script:
@gap("Read(\"libs/clifford_group.g\");;");

@gap("LoadPackage(\"forms\");;");

function clifford_group(n)

    n_string = string(n);
    GAP.evalstr("n:= "*string(n)*";;");
    GAP.evalstr("cl"*n_string*"_pair := clifford_group(n);;");
    GAP.evalstr("nat"*n_string*" := cl"*n_string*"_pair[1];;");
    GAP.evalstr("cl"*n_string*" := cl"*n_string*"_pair[2];;");
end

function pauli_to_bit_strings(n)
    n_string = string(n);
    GAP.evalstr("n := "*n_string*";");
    GAP.evalstr("E"*n_string*" := all_bit_strings(2*n);");
    GAP.evalstr("P"*n_string*"_pair := list_bit_to_pauli_string(E"*n_string*");;");
    Pn_string = (GAP.gap_to_julia(GAP.evalstr("P"*n_string*"_pair[1];")));
    En = (GAP.gap_to_julia(GAP.evalstr("P"*n_string*"_pair[2];")));

    return Pn_string, En;
end


function generate_pauli_coefficients(pauli_keys,input_pauli,input_values)

    # generate dictionary between paulis and index 1,...,4^n:
    idx_list = [i for i in 1:length(pauli_keys)];
    dict = Dict(zip(pauli_keys,idx_list));

    pauli_idx = [dict[A] for A in input_pauli];

    C = [0 for i in 1:length(pauli_keys)];

    C[pauli_idx] .= input_values;

    return C
end


function generate_rational_coefficient_string(C)
    # convert to julia rational vector
    C = Vector{Rational{Int64}}(C);
    # convert to string and remove "Rational{Int64}:" characters:
    C_string = string(C)[16:end];

    # initialize output string:
    S = ""; 
    for i in 1:length(C_string)
        if (C_string[i] == '/') && (C_string[i+1] == '/')
            continue;
        end;
        S = S*C_string[i];
    end
    return S
end

# NOTE: should have previously generated clifford group:
function clifford_orbit(C,n)
    nn = string(n);

    # only works for rational inputs:
    CC = generate_rational_coefficient_string(C);

    C_orbit = GAP.evalstr("clifford_orbit(nat"*nn*",cl"*nn*","*CC*",E"*nn*")");

    return GAP.gap_to_julia(C_orbit);
end



# convert array to matrix:
# assume elements are vectors of same length:
function array_to_matrix(typ, arry)
    n = length(arry); d = length(arry[1]);
    # intialize matrix:
    M = Matrix{typ}(undef, 0, d);
    for i in 1:n
        M = vcat(M,transpose(arry[i]));
    end;
    return M
end


# convert array to polymake matrix:
# input: typ: pm.Rational, pm.Integer, etc.
# input: julia array:
function array_to_polymake_matrix(typ, arry)
    M = array_to_matrix(typ,arry);
    return pm.Matrix{typ}(M)
end


################################################################################################

# symplectic action on subspace:

function symplectic_orbit(n,subspace)
    if string(subspace)[1] == 'A'
        ss = "Set("*string(subspace)[4:end]*")";
    else
        ss = "Set("*string(subspace)*")";
    end

    nn = string(n);

    #input to GAP:
    GAP.evalstr("n:= "*nn*";;");
    GAP.evalstr("subspace := "*ss*";;");

    # call symplectic_orbit function:
    orbit = (GAP.evalstr("symp_orbit := symplectic_orbit_subspace(n,subspace);"));

    return GAP.gap_to_julia(orbit)
end