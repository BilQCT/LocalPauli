import re
from typing import Any, Dict, List, Tuple

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
from graphix import Circuit, GraphState
from graphix.clifford import Clifford
from graphix.command import C, E, M, N
from graphix.fundamentals import Axis
from graphix.graphsim import GraphState
from graphix.measurements import Measurement
from graphix.parameter import Placeholder
from graphix.pattern import Pattern


class MBPCCompiler:
    """
    Parses QASM circuits and compiles them into the Measurement-Based
    Pauli Computation (MBPC) model using the Graphix library.

    This compiler maps gate-based quantum circuits into optimized
    MBQC measurement patterns, specifically tailored for MBPC by performing
    Pauli measurements classically and reducing the resource state graph.

    Assumptions:
        1. Initial State: All logical qubits are implicitly initialized in |+>.
        2. Measurement Basis: QASM 'measure' operations are Pauli X-basis.
        3. T-Gate Parameterization: T-gates use symbolic `Placeholder`s.
        4. Strict Graph Reduction: Applies E1 & E2 to eliminate HS/HSZ unitaries.
        5. Theorem 3.1 Commutation: Commutes remaining S and Z decorations past
           the T-gates, absorbing them into X/Y measurement bases and classical
           outcome flips, leaving ONLY Hadamard (H) decorations on the graph.
    """

    def __init__(self):
        self.t_count = 0

    def parse_qasm_string(self, qasm_str: str) -> Circuit:
        """
        Parses a QASM string and builds a Graphix Circuit object.

        Args:
            qasm_str (str): A string containing OpenQASM instructions.

        Returns:
            Circuit: A Graphix Circuit object representing the parsed quantum circuit.
        """
        lines = qasm_str.strip().split("\n")
        num_qubits = 0

        # 1. Determine the number of qubits from the qreg instruction
        for line in lines:
            line = line.strip()
            if line.startswith("qreg"):
                match = re.search(r"qreg\s+\w+\[(\d+)\]", line)
                if match:
                    num_qubits = int(match.group(1))
                    break

        if num_qubits == 0:
            raise ValueError("Could not find 'qreg' definition in QASM.")

        circuit = Circuit(num_qubits)
        self.t_count = 0

        # 2. Parse gates sequentially
        for line in lines:
            # Strip comments and whitespace
            line = line.split("//")[0].strip()
            if not line or line.startswith(("OPENQASM", "include", "qreg", "creg")):
                continue

            # Match standard QASM gate syntax: gate_name q[0], q[1];
            match = re.match(r"([a-z]+)\s+(.*);", line)
            if not match:
                continue

            gate = match.group(1)
            args_str = match.group(2)

            # Extract all qubit indices using regex
            qubits = [int(idx) for idx in re.findall(r"\[(\d+)\]", args_str)]

            self._apply_gate(circuit, gate, qubits)

        return circuit

    def _apply_gate(self, circuit: Circuit, gate: str, qubits: List[int]) -> None:
        """
        Maps a parsed QASM gate identifier to the corresponding Graphix instruction.
        T-gates are parameterized to ensure they are not reduced during the Pauli reduction step.

        Args:
            circuit (Circuit): The Graphix Circuit object being constructed.
            gate (str): The name of the gate (e.g., 'h', 'cx', 't').
            qubits (List[int]): The list of target (and control) qubit indices.
        """
        if gate == "h":
            circuit.h(qubits[0])
        elif gate == "x":
            circuit.x(qubits[0])
        elif gate == "y":
            circuit.y(qubits[0])
        elif gate == "z":
            circuit.z(qubits[0])
        elif gate == "s":
            circuit.s(qubits[0])
        elif gate == "cx":
            circuit.cnot(control=qubits[0], target=qubits[1])
        elif gate == "cz":
            circuit.cz(qubits[0], qubits[1])
        elif gate == "t":
            # Parameterize the T-gate to protect it from Pauli reduction
            self.t_count += 1
            circuit.rz(qubits[0], Placeholder(f"T_{self.t_count}"))
        elif gate == "measure":
            # Assumption: MBQC natively utilizes X-basis measurements
            circuit.m(qubits[0], Axis.X)
        else:
            print(f"Warning: Unsupported gate '{gate}' skipped.")

    def compile_to_mbpc(self, circuit: Circuit) -> Tuple[Pattern, List[Dict[str, Any]]]:
        """
        Compiles the quantum circuit into a reduced MBPC measurement pattern.

        This method basically applies Theorem 3.1 from our paper.
        It uses the GraphState simulator to perform all Pauli measurements classically,
        effectively collapsing the graph down to only the non-Clifford (T-gate) nodes.

        This method leverages Theorem 3.1 of MBPC formalisms:
        It standardizes the pattern, locks the inputs to |+>, and uses the
        efficient GraphState simulator to eliminate all Pauli measurements,
        collapsing the graph down to only the non-Clifford (T-gate) nodes.

        Args:
            circuit (Circuit): The standard Graphix circuit to be compiled.

        Returns:
            Tuple[Pattern, List[Dict[str, Any]]]:
                - pattern: The fully optimized and reduced MBQC Pattern.
                - t_nodes_info: A list of dictionaries containing metadata for
                  each remaining non-Pauli measurement (node ID, target angle,
                  measurement plane, and feedforward adaptivity domains).
        """
        # 1. Transpile to base MBQC pattern and standardize
        pattern = circuit.transpile().pattern
        pattern.standardize()
        pattern.shift_signals()

        # 2. Enforce the input state as |+> (required for the reduction)
        pattern.remove_input_nodes()

        # 3. Graph reduction over all Pauli operators
        # This will absorb/eliminate all measurements except the Placeholder T-gates
        pattern.perform_pauli_measurements()

        t_nodes_info = []
        mbpc_cmds = []

        # 4. Separate Graph Topology (N, E, C) from Execution Commands (M, X, Z, S)
        nodes = []
        edges = []
        vops = {}
        other_cmds = []

        for cmd in list(pattern):
            if isinstance(cmd, N):
                nodes.append(cmd.node)
            elif isinstance(cmd, E):
                edges.append(cmd.nodes)
            elif isinstance(cmd, C):
                vops[cmd.node] = cmd.clifford
            else:
                other_cmds.append(cmd)

        # 5. Build GraphState simulator and apply Elliott et al. strict reductions
        g = GraphState(nodes=nodes, edges=edges, vops=vops)

        while True:
            reduced = False

            # Rule E1: Hollow nodes (H) cannot have loops (S)
            for n in list(g.nodes):
                if g.nodes[n].get("hollow", False) and g.nodes[n].get("loop", False):
                    g.equivalent_graph_e1(n)
                    reduced = True
                    break

            if reduced:
                continue

            # Rule E2: Hollow nodes (H) cannot be connected to other Hollow nodes
            for e in list(g.edges):
                u, v = e
                if g.nodes[u].get("hollow", False) and g.nodes[v].get("hollow", False):
                    g.equivalent_graph_e2(u, v)
                    reduced = True
                    break

            if not reduced:
                break  # The graph is fully reduced!

        # 6. Extract the fully reduced topology and local Clifford decorations
        reduced_nodes = list(g.nodes)
        reduced_edges = list(g.edges)
        reduced_vops = g.extract_vops()

        # 7. Rebuild the modified pattern
        new_cmds = []
        for n in reduced_nodes:
            new_cmds.append(N(node=n))
        for e in reduced_edges:
            new_cmds.append(E(nodes=e))

        t_nodes_info = []
        t_nodes_set = {
            cmd.node
            for cmd in other_cmds
            if isinstance(cmd, M) and cmd.measurement.try_to_pauli() is None
        }

        # Commute S and Z past T-gates, leaving only H decorations on the graph!
        for cmd in other_cmds:
            if isinstance(cmd, M):
                if cmd.measurement.try_to_pauli() is None:

                    cliff_obj = reduced_vops.get(cmd.node)
                    v_str = cliff_obj.name if cliff_obj is not None else "I"

                    s_i = 1 if "H" in v_str else 0
                    r_i = 1 if "S" in v_str else 0
                    u_i = 1 if "Z" in v_str else 0

                    # 1. THEOREM 3.1: Leave ONLY the Hadamard decoration
                    if s_i == 1:
                        new_cmds.append(C(node=cmd.node, clifford=Clifford.H))

                    # 2. Commute S past T: changes X measurement to Y measurement
                    mbpc_basis = "Y" if r_i == 1 else "X"

                    # 3. Commute Z past T: causes an unconditional measurement outcome flip
                    uncond_flip = True if u_i == 1 else False

                    t_nodes_info.append(
                        {
                            "node": cmd.node,
                            "angle": str(cmd.measurement.angle),
                            "plane": cmd.measurement.to_plane().name,
                            "h_decoration": bool(s_i),
                            "basis": mbpc_basis,
                            "unconditional_flip": uncond_flip,
                            "s_domain": set(cmd.s_domain),
                            "t_domain": set(cmd.t_domain),
                        }
                    )

                    # Mask the MBPC T-gate using the correctly commuted X or Y basis
                    meas_op = Measurement.Y if mbpc_basis == "Y" else Measurement.X
                    x_meas = M(
                        node=cmd.node,
                        measurement=meas_op,
                        s_domain=cmd.s_domain,
                        t_domain=cmd.t_domain,
                    )
                    new_cmds.append(x_meas)
                else:
                    new_cmds.append(cmd)
            else:
                new_cmds.append(cmd)

            # Restore any non-T-gate C commands (e.g., if there are unmeasured output qubits)
            for n, clifford in reduced_vops.items():
                if n not in t_nodes_set:
                    new_cmds.append(C(node=n, clifford=clifford))

        pattern.replace(new_cmds)
        return pattern, t_nodes_info


def print_mbpc_theorem_format(pattern, t_nodes_info):
    """
    Outputs the state and adaptivity exactly as formatted in
    Theorem 3.1 of the MBPC formalism, in a highly compact format.
    """
    nodes = []
    edges = []
    for cmd in pattern:
        if isinstance(cmd, N):
            nodes.append(cmd.node)
        elif isinstance(cmd, E):
            edges.append(cmd.nodes)

    print("============================================================")
    print(" MBPC RESOURCE STATE & ADAPTIVITY (Theorem 3.1 Normal Form)")
    print("============================================================")

    print(f"1. TOPOLOGY: Nodes (t={len(nodes)}): {nodes} | Edges: {edges}")

    # First, safely extract all the T-gate identifiers (e.g., T_1, T_2)
    t_labels = {}
    for info in t_nodes_info:
        angle_str = str(info["angle"])
        match = re.search(
            r"(?:[+-]?\s*\d*(?:\.\d+)?\s*\*\s*)?([A-Za-z0-9_]+)", angle_str
        )
        t_labels[info["node"]] = match.group(1) if match else "T"

    # Combine all H-decorations onto a single line with the T_k identifiers
    h_decs = [
        f"Node {info['node']} ({t_labels[info['node']]}): {'H' if info['h_decoration'] else 'I'}"
        for info in t_nodes_info
    ]
    print(f"2. H-DECORATIONS (V_i): {', '.join(h_decs)}")

    print("3. ADAPTIVE MEASUREMENTS:")
    for info in t_nodes_info:
        node = info["node"]
        basis = info["basis"]
        uncond_flip = "Yes" if info["unconditional_flip"] else "No"

        # Format empty sets cleanly as {} instead of set()
        s_dom = set(info["s_domain"]) if info["s_domain"] else "{}"
        t_dom = set(info["t_domain"]) if info["t_domain"] else "{}"

        # Grab the identifier we parsed earlier
        var_name = t_labels[node]

        # Print all the adaptivity rules compactly on one line
        print(
            f"   Node {node:<2}({var_name}): Basis = {basis} | s_dom = {s_dom} | t_dom = {t_dom} | UncondFlip = {uncond_flip}"
        )

    print("============================================================")


def draw_mbpc_graph(pattern, t_nodes_info):
    """
    Draws the reduced MBPC resource state graph using networkx.
    Color codes:
      - Red: Nodes with an Hadamard (H) decoration.
      - Blue: Nodes with No decoration (I).
    """
    G = nx.Graph()
    h_nodes = set()

    # 1. Build the graph topology from the Pattern commands
    for cmd in pattern:
        if isinstance(cmd, N):
            G.add_node(cmd.node)
        elif isinstance(cmd, E):
            G.add_edge(cmd.nodes[0], cmd.nodes[1])
        elif isinstance(cmd, C) and cmd.clifford == Clifford.H:
            h_nodes.add(cmd.node)

    # 2. Extract T-gate labels for pretty printing
    t_labels = {}
    for info in t_nodes_info:
        angle_str = info["angle"]
        match = re.search(
            r"(?:[+-]?\s*\d*(?:\.\d+)?\s*\*\s*)?([A-Za-z0-9_]+)", angle_str
        )
        var_name = match.group(1) if match else "T"
        t_labels[info["node"]] = var_name

    # 3. Assign colors and labels
    color_map = []
    labels = {}
    for node in G.nodes():
        # Assign Color based on H-decoration presence
        if node in h_nodes:
            color_map.append("#ff9999")  # Light Red for H
        else:
            color_map.append("#99ccff")  # Light Blue for I

        # Assign Label (Node ID + T_k identifier)
        if node in t_labels:
            labels[node] = f"Node {node}\n({t_labels[node]})"
        else:
            labels[node] = f"Node {node}\n(Output)"

    # 4. Draw the graph
    plt.figure(figsize=(8, 6))

    # spring_layout provides a nice force-directed graph spacing
    pos = nx.spring_layout(G, seed=42)

    nx.draw(
        G,
        pos,
        node_color=color_map,
        with_labels=True,
        labels=labels,
        node_size=2500,
        font_size=10,
        font_weight="bold",
        edge_color="#555555",
        width=2,
        edgecolors="black",
    )  # Adds a border around nodes

    # 5. Add a legend
    h_patch = mpatches.Patch(facecolor="#ff9999", edgecolor="black", label="H Gate (H)")
    i_patch = mpatches.Patch(
        facecolor="#99ccff", edgecolor="black", label="No H Gate (I)"
    )

    # Let matplotlib find the best empty corner
    plt.legend(handles=[h_patch, i_patch], loc="best", fontsize=11)

    # Add 10% whitespace padding around the graph
    plt.margins(0.10)

    plt.title("MBPC Reduced Resource State Topology", fontsize=14, fontweight="bold")
    plt.axis("off")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    qasm_example_3_2 = """
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[1];
    creg c[1];
    t q[0];     // Acts directly on implicit |+>
    h q[0];
    t q[0];
    measure q[0] -> c[0]; // Implicitly X-basis
    """

    print("Parsing QASM and compiling to MBPC...")
    compiler = MBPCCompiler()
    circuit = compiler.parse_qasm_string(qasm_example_3_2)

    # Run the compiler (which performs Pauli reduction, E1/E2 elimination, and S/Z commutation)
    pattern, t_nodes_info = compiler.compile_to_mbpc(circuit)

    # Print the strictly formatted Theorem 3.1 output
    print_mbpc_theorem_format(pattern, t_nodes_info)

    # Draw the graphical representation!
    print("Drawing the MBPC Resource Graph...")
    draw_mbpc_graph(pattern, t_nodes_info)
