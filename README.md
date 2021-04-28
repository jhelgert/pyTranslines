
# pyTranslines

This is a python package to control networks of electrical transmission
lines. The dynamics on the lines are described by the telegrapher equations,
a 2x2 system of PDEs from the class of hyperbolic balance laws. There are
two controls inside the network:

- A so-called *outer control* to disable single lines inside the network
(or whole subgrids) at an arbitrary point in time.
- And a so-called *inflow control* that determines the power inflow into
the network at the source vertices.

The goal is to
minimize the sum of quadratic deviation between the given demand and the
delivered load to the consumers during a time interval `[0, T]`.

The packages uses a *first discretize, then optimize* ansatz to transform
the infinite dimensional optimization problem on function spaces into a
finite dimensional optimization problem. More precisely, the package
solves a mixed-integer quadratically constrained quadratic program (MIQCQP)
by Gurobi. Additionally, it contains a specialized heuristic coupled
with Gurobi and based on the CIAP decomposition to quickly generate feasible solutions. 
In particular for large problems the latter gives a significantly speed up.

## Install

First install [pyCIAP](https://github.com/jhelgert/pyCIAP). Then, clone this repo and run
``` bash
python3 setup.py install
```
inside the repo folder.

## Example

![](https://i.imgur.com/w7pE1iS.png)

Here, `u0` and `u1` are inflow controls at the source vertices 0 and 1 and
`Q9, Q10, Q11, Q12, Q13` are given demand functions at the sink (consumer)
vertices 9,10,11,12,13. The network contains two switches `s1` and `s2`
to disable lines (or a whole subgrid) inside the network. 

The goal is to
minimize the sum of quadratic deviation between the demand and the
delivered load to the consumers during a time interval `[0, T]`.

``` python
from pyTranslines import Translines

# vertices
V = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]

# arcs
A = [(0, 2), (1, 3), (3, 5), (2, 4), (2, 8), (4, 6), (5, 6),
     (6, 7), (8, 7), (2, 9), (8, 10), (3, 11), (7, 12), (7, 13)]

# producer / source vertices
producers = [0, 1]

# consumer / sink vertices
consumers = [9, 10, 11, 12, 13]

# Each config contains the disabled lines if the config is active.
# Example: 2 corresponds to the line A[2] = (3, 5)
configs = [(), (2, 3, 7), (8,), (2, 3, 7, 8)]

# Given consumer demand
demand = np.loadtxt('demand_dt_0_5.dat')

# end of time horizon [0, T]
T = 10

# Create our discretized control problem
Prob = Translines(V, A, producers, consumers, configs, demand)

# T time horizon, lr length of all network arcs, L, C, R, G physical parameters
Prob.set_parameters(T=T, lr=1.0, L=1.0, C=1.0, R=1.0e-3, G=2.0e-3)

# discretization step sizes
Prob.set_step_sizes(dt=0.5, dx=0.5)

# upper bounds for the power inflow at the source vertices
Prob.set_inflow_UB([120, 80])

# Set the configuration dwell times (as number of time steps, i.e. multiple of dt)
# for all configurations: 3 / dt = 3 / 0.5 = 1.5
Prob.set_min_dwell_times(min_up=3, min_down=3)

# Set the transmission disturbances on a specific line
Prob.set_disturbed_lines([13], [20.0], [500*(T/0.5 + 1)])

# bigM constraint used for linearizing the coupling conditions
Prob.set_BigM(150.0)

# Solve via Gurobi
Prob.solve(solver="Gurobi")

# Solve via the specialized heuristic
# Prob.solve(solver="GurobiCIAP")

Prob.plot_results()
```
