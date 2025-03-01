
# Introduction

This This is the GitHub repository companion to [Okay, et al. 2024] "Classical simulation of universal measurement-based quantum computation using multipartite Bell scenarios" by Okay, Yücel, Ipek, [arXiv:2410.23734](https://arxiv.org/abs/2410.23734).

## Background

Among the central contributions of arXiv:2410.23734 is a probabilistic representation of Measurement-Based Pauli Computation (MBPC), a variant of Measurement-Based Quantum Computing (MBQC) that uses single-qubit Pauli measurements and is also universal; see e.g., [[Danos, et al. 2006]](https://arxiv.org/abs/0704.1263), [[Raussendorf, et al. 2017]](https://arxiv.org/abs/1511.08506).

Since MBPC is driven by single-qubit Pauli measurements, the computation does not initiate on a graph state, as this would only compute stabilizer circuits. Rather, we use a resource state, called a *magic cluster state*, which is a modification of a graph state that exhibits both entanglement and non-stabilizerness. Given a graph $G = (V,E)$ on $n$ vertices and a subset $U\subset V$, an $n$-qubit magic cluster state in the $n$-qubit Hilbert space $\mathcal{H}_n = (\mathbb{C}^2)^{\otimes n}$ is defined by

$$\left |\psi_{G,U} \right \rangle = E(G)T^{\otimes U} \left |+ \right \rangle^{\otimes n}\quad \text{where}\quad E(G) = \prod_{(i,j)\in E}CZ_{ij}$$

The central objects of our classical simulation algorithm are the Local Pauli (LP) polytope, convex geometric structures that enclose the set of quantum states; including magic cluste states. The $n$-qubit LP polytope, denoted $\text{LP}_n$, can be identified with $\text{NS}_n$, the non-signaling polytopes of the (n,3,2) multi-partite Bell scenario for all $n$.

An important structural property of LP polytopes is that they are closed under the tensor product operation familiar from quantum mechanics. In particular, for a vertex $A_1\in \text{LP}_{n_1}$ and vertex $A_2\in \text{LP}_{n_2}$ we have that $A_1\otimes A_2$ is a vertex of $\text{LP}_{n_1+n_2}$.

In [Okay, et al. 2024] a distinguished class of operators, called *locally closed (LC)* operators, was introduced. Two Pauli operators are said to *locally commute* if they commute on every tensor factor. A subset of Pauli operators is said to be locally closed if for every two locally commuting Pauli operators $A$ and $B$ in that subset, then their product $A\cdot B$ is also in the subset. For a single qubit $\text{LP}_1$ is a cube enclosing the Bloch sphere whose vertices are LC operators. The $n$-fold tensor product of $\text{LP}_1$ vertices is an LC operator that is always a vertex of $\text{LP}_n$ and can be identified with the deterministic vertices of $\text{NS}_n$; i.e., the vertices of the Bell polytope. Interestingly, the vertices of $\text{LP}_2$, the two-qubit LP polytope also consists only of LC operators. Thus the tensor product of $\text{LP}_1$ and $\text{LP}_2$ vertices is also locally closed. 

LC operators with efficient representation (in $n$) and update rules can be organized into a phase space that can be used for classical simulation based on sampling from a quasi-probability distribution. The LC phase space subsumes the phase space introduced in [[Raussendorf, et al. 2017]](https://arxiv.org/abs/1905.05374) based on closed non-contextual (CNC) sets, which itself generalizes the stabilizer phase space of (Howard and Campbell, 2016). Some examples of LC phase spaces include:

- $D$:$\quad$ Deterministic vertices; i.e., vertices of the $(n,3,2)$ Bell polytope.
- $L$:$\quad$ Tensor product of $\text{LP}_1$ and $\text{LP}_2$ vertices,
- $L'$:$\quad$ Deterministic vertices and the $n$-qubit CNC operators.

For any phase space we can define a robustness monotone:

$$\mathfrak{R}(\rho) := \min \left\{\sum_{\alpha\in V} |Q(\alpha)|~:~\rho = \sum_{\alpha\in V} Q(\alpha) A_\alpha \right\},$$

where $\alpha \in V$ labels the operators in the phase space.


## Code contents


In this GitHub repository we provide the $L_1$ and $L_2$ phase space operators for $n=3$ and deterministic vertices for $n=3,4,5$. We encode the phase spaces in matrices whose columns are labeled by $\alpha \in V$ and whose rows are Pauli coefficients of the corresponding operator $A_\alpha$. In particular, for some phase space operator $A_\alpha$ and a Pauli $P$, the matrix has components $M_{P,\alpha} = \text{Tr}(PA_\alpha) $.  We employ lexicographic ordering on the rows so that $I\otimes \cdots \otimes I$ precedes $I\otimes \cdots \otimes X$, which precedes $I\otimes \cdots \otimes Y$, and so on.

Our matrices have size:

**Three-qubits**

- $D$:$\quad$ $64\times 512$
- $L$:$\quad$ $64\times 27,008$
- $L'$:$\quad$ $64\times 71,648$


In the directory **libs** we have basic functionality for manipulating Pauli operators. The directory **local_robustness** has functions for computing the robustness for the $D$, $L_1$, and $L_2$ phase spaces, as well as the CNC phase space. We compute the robustness for two families of resource states

\begin{align}
\text{PBC~:}\quad &     \left |\psi_{G,U}^{PBC} \right \rangle = E(G)H^{\otimes U}T^{\otimes 3} \left |+++ \right \rangle,\notag\\
\text{MBQC:}\quad &     \left |\psi_{G,U}^{MBQC} \right \rangle =T^{\otimes 3} H^{\otimes U}E(G) \left |+++ \right \rangle,\notag
\end{align}

where $U$ runs over all possible subsets of $\{1,2,3\}$ and $G$ runs over all connected graphs on $3$ qubits. 


# Requirements:

Combinatorics   v1.0.2
GLPK            v1.2.1
HDF5            v0.17.2
JLD2            v0.4.53
JuMP            v1.23.0
YAML            v0.4.12

# Instructions

To create an environment for running LocalPauli, in Julia REPL use the following commands:

using Pkg
Pkg.generate("path/to/MyProject")
Pkg.activate("path/to/MyProject")

Pkg.add("YourPackage")


# Disclaimer on Intellectual Property

This repository contains material related to our arXiv preprint arXiv:2410.23734. Please note that the technology described herein is subject to a patent pending (application number: 18/925,447). By using, modifying, or distributing this material under the Apache License 2.0, you acknowledge that certain patent rights may apply. We request that users review the terms of the license and this notice to ensure compliance with our intellectual property rights.

If you have any questions regarding the scope of the patent or usage rights, please contact us at cihan.okay@bilkent.edu.tr.


