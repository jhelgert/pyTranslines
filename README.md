
# pyTranslines

...

## Example

``` python
from pyTranslines import Translines

# Vertices
V = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
# Arcs
A = [(0, 2), (1, 3), (3, 5), (2, 4), (2, 8), (4, 6), (5, 6),
     (6, 7), (8, 7), (2, 9), (8, 10), (3, 11), (7, 12), (7, 13)]
# producer / source vertices 
producers = [0, 1]
# consumer / sink vertices
consumers = [9, 10, 11, 12, 13]
# Each config contains the disabled lines if the config is active.
# Example: 2 corresponds to the line A[2] = (3, 5)
configs = [(), (2, 3, 7), (8,), (2, 3, 7, 8)]
# Given demand
demand = np.loadtxt('demand_dt_0_5.dat')

# end of time horizon [0, T]
T = 10
# Create our discretized control problem
Prob = Translines(V, A, producers, consumers, configs, demand)
Prob.set_parameters(T=T, lr=1.0, L=1.0, C=1.0, R=1.0e-3, G=2.0e-3)
Prob.set_step_sizes(dt=0.5, dx=0.5)
Prob.set_inflow_UB([120, 80])
# Set the configuration dwell times (as number of time steps, i.e. multiple of dt)
# for all configurations: 3 / dt = 3 / 0.5 = 1.5
Prob.set_min_dwell_times(min_up=3, min_down=3)
# Set the transmission disturbances on a specific line
Prob.set_disturbed_lines([13], [20.0], [500*(T/0.5 + 1)])
# bigM constraint used for linearizing the coupling conditions
Prob.set_BigM(150.0)

Prob.solve(solver="Gurobi")
Prob.plot_results()
```