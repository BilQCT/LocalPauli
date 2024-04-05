using PlotlyJS

include("cnc.jl")


function create_3dim_table_with_given_3rd_party_local(omega::Set{Vector{Int}})
    
    header_paulis = ["III", "IXI", "IYI", "IZI", "", "IIX", "IXX", "IYX", "IZX"]

    cell_paulis = [
        ["XII", "YII", "ZII", "", "IIY", "XIY", "YIY", "ZIY"],
        ["XXI", "YXI", "ZXI", "", "IXY", "XXY", "YXY", "ZXY"],
        ["XYI", "YYI", "ZYI", "", "IYY", "XYY", "YYY", "ZYY"],
        ["XZI", "YZI", "ZZI", "", "IZY", "XZY", "YZY", "ZZY"],
        ["", "", "", "", "", "", "", ""],
        ["XIX", "YIX", "ZIX", "", "IIZ", "XIZ", "YIZ", "ZIZ"],
        ["XXX", "YXX", "ZXX", "", "IXZ", "XXZ", "YXZ", "ZXZ"],
        ["XYX", "YYX", "ZYX", "", "IYZ", "XYZ", "YYZ", "ZYZ"],
        ["XZX", "YZX", "ZZX", "", "IZZ", "XZZ", "YZZ", "ZZZ"]
    ]

    
    cell_paulis_values = [["<b>" * string(pauli) * "</b>" for pauli in row]  for row in cell_paulis]

    header_paulis_values = ["<b>" * string(pauli) * "</b>" for pauli in header_paulis]

    header_colors = []
    cell_colors = []

    in_set_color = "#51A664"
    out_set_color = "#BA4C55"

    blank_color = "white"

    for pauli_string in header_paulis
        if pauli_string == ""
            push!(header_colors, blank_color)
            continue
        end
        pauli = get_pauli_from_pauli_string(pauli_string)
        if pauli in omega
            push!(header_colors, in_set_color)
        else
            push!(header_colors, out_set_color)
        end
    end

    for col in cell_paulis
        col_colors = []
        for pauli_string in col
            if pauli_string == ""
                push!(col_colors, blank_color)
                continue
            end
            pauli = get_pauli_from_pauli_string(pauli_string)
            if pauli in omega
                push!(col_colors, in_set_color)
            else
                push!(col_colors, out_set_color)
            end
        end

        push!(cell_colors, col_colors)
        
    end

    line_colors = ["darkslategray" for _ in 1:8]
    line_colors[4] = "black"

    paulis_table = table(
        header=attr(
            values=header_paulis_values,
            line_color="darkslategray",
            fill_color=header_colors,
            align=["center"],
            font=attr(color="white", size=12)
        ),
        cells=attr(
            values=cell_paulis_values,
            line_color=["darkslategray"],
            fill_color=cell_colors,
            align=["center"],
            font=attr(color="white", size=12)
        ),
        columnwidth=[3, 3, 3, 3, 1, 3, 3, 3, 3]
    )

    return paulis_table

end

function draw_3dim_paulis(omega::Set{Vector{Int}}) 
    III_table = create_3dim_table_with_given_3rd_party_local(omega)
    p = plot(III_table)
    relayout!(p, showlegend=false, title_text="3-dimensional Pauli Operators")

    return p

end
