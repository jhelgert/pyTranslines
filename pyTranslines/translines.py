import numpy as np
from gurobipy import Env, GRB, Model, quicksum as qsum
import matplotlib.pyplot as plt
from math import ceil
from pyCIAP import DSUR


class Translines():
    def __init__(self, V, A, producers, consumers, configs, demand) -> None:
        # network data
        self.V = V
        self.A = A
        self.producers = producers
        self.consumers = consumers
        # list of configurations, each config contains the disabled lines
        # for each switches configuration
        self.configs = configs
        self.noc = len(configs)
        # given demand for the consumers
        self.demand = demand
        # create the sets of outgoing/incoming arcs for each vertex
        self.__create_deltaInOut()
        # create the equal distribution matrices based on the given network
        self.__create_distribution_matrices()
        # create the set of the inner arcs, i.e. all arcs not in producers
        # or consumers
        self.__create_Vinner()
        # big M parameter for the linearization of the coupling conditions
        self.bigM = None

        self.min_down_dwell_times = None
        self.min_up_dwell_times = None
        self.transmission_disturbances = None

        # gurobi environment (only needed to disable outputs to the console)
        self.grb_env = Env()

    def set_parameters(self, T, lr, L, C, R, G) -> None:
        # time horizon
        self.T = T
        # arc length for all arcs
        self.lr = lr
        # physical constants for all arcs
        self.L = L
        self.C = C
        self.R = R
        self.G = G

    def set_BigM(self, bigM) -> None:
        self.bigM = bigM

    def set_step_sizes(self, dt, dx) -> None:
        assert dt <= np.sqrt(self.L*self.C) * \
            dx, "step sizes don't satisfy CFL condition"
        self.dt = dt
        self.dx = dx
        self.nt = int(self.T/dt) + 1
        self.nx = int(self.lr/dx) + 1
        self.timegrid = np.arange(0, self.T + dt, dt)

    def set_inflow_UB(self, ubs) -> None:
        ''' sets the upper bounds for the inflow controls '''
        assert len(ubs) == len(
            self.producers), "number of UBs doesn't match number of producers"
        self.inflow_UBs = np.asarray(ubs)

    def set_min_dwell_times(self, min_up, min_down) -> None:
        # Same dwell times for each configuration
        self.min_down_dwell_times = [min_down]*len(self.configs)
        self.min_up_dwell_times = [min_up]*len(self.configs)

    def set_disturbed_lines(self, lines, load_averages, thresholds) -> None:
        assert len(lines) == len(load_averages) == len(
            thresholds), "Number of disturbed lines doesn't match the len of the other arguments"
        self.transmission_disturbances = True
        self.disturbed_lines = lines
        self.xiPlus_mean = load_averages
        self.load_threshold = thresholds

    def __create_deltaInOut(self) -> None:
        self.delta_in = [[] for _ in self.V]
        self.delta_out = [[] for _ in self.V]
        for v, (i, j) in zip(self.V, self.A):
            self.delta_in[v] += [w for w,
                                 (k, l) in zip(self.V, self.A) if l == v]
            self.delta_out[v] += [w for w,
                                  (k, l) in zip(self.V, self.A) if k == v]

    def __create_Vinner(self) -> None:
        self.V_inner = [
            v for v in self.V if v not in self.producers+self.consumers]

    def __create_distribution_matrices(self) -> None:
        self.Dplus = np.zeros((len(self.configs), len(self.A), len(self.A)))
        self.Dminus = np.zeros((len(self.configs), len(self.A), len(self.A)))
        for config, off_arcs in enumerate(self.configs):
            for k, (i, j) in enumerate(self.A):
                inc1 = [r for r in self.delta_out[j] if r not in off_arcs]
                inc2 = [r for r in self.delta_in[i] if r not in off_arcs]
                if len(inc1) > 0:
                    self.Dplus[config, inc1, k] = 1.0 / len(inc1)
                if len(inc2) > 0:
                    self.Dminus[config, inc2, k] = 1.0 / len(inc2)

    def __build_MIQCQP(self, linearize_coupling=True) -> None:
        # Entries of the matrices B, i.e. all B_r are the same for all arcs
        # r in A
        b11 = 0.5 * (self.R/self.L + self.G/self.C)
        b21 = 0.5 * (self.R/self.L - self.G/self.C)
        b12 = b21
        b22 = b11

        # flux values from the characteristic formulation
        lambda_p = 1.0/np.sqrt(self.L*self.C)
        lambda_m = -1.0*lambda_p

        # Get rid of all the 'self' in front of the attributes for better
        # readability inside this method
        V = self.V
        A = self.A
        V_inner = self.V_inner
        consumers = self.consumers
        producers = self.producers
        configs = self.configs
        demand = self.demand
        Dplus = self.Dplus
        Dminus = self.Dminus
        delta_in = self.delta_in
        delta_out = self.delta_out
        nt = self.nt
        nx = self.nx
        dt = self.dt
        dx = self.dx
        bigM = self.bigM

        # Create the gurobi model
        self.miqcqp = Model()

        # binary outer control variables b[config, t]
        b = self.miqcqp.addVars(len(configs), nt, vtype='B')

        # continuous inflow control variables u[producer, t]
        u = self.miqcqp.addVars(producers, nt, lb=0.0)

        # set the upper bounds for the inflow controls
        if self.inflow_UBs is not None:
            for i, producer in enumerate(producers):
                for t in range(nt):
                    u[producer, t].ub = self.inflow_UBs[i]

        # continous characteristic variables xiPlus[r, x, t] and xiMinus[r, x, t]
        xiPlus = self.miqcqp.addVars(len(A), nx, nt)
        xiMinus = self.miqcqp.addVars(len(A), nx, nt)

        # Save the variables for later use
        self.binControl = b
        self.xiPlus = xiPlus
        self.xiMinus = xiMinus

        if linearize_coupling:
            # indices of the arcs affected by coupling constraints
            I1 = [r for r, (i, j) in enumerate(A) if j not in consumers]
            I2 = [r for r, (i, j) in enumerate(A) if i not in producers]

            # continuous helper variables hPlus[config, r, t] and hMinus[config, r, t]
            # Will be used for linearizing the coupling constraints
            hPlus = self.miqcqp.addVars(len(configs), I1, nt)
            hMinus = self.miqcqp.addVars(len(configs), I2, nt)

        # update the model to add all variables
        self.miqcqp.update()

        # set boundary values (consumer)
        for consumer in consumers:
            for r in delta_in[consumer]:
                for t in range(nt):
                    xiMinus[r, nx-1, t].lb = 0.0
                    xiMinus[r, nx-1, t].ub = 0.0

        # set inflow boundary conditions (producer)
        for producer in producers:
            for r in delta_out[producer]:
                for t in range(1, nt):
                    self.miqcqp.addConstr(
                        xiPlus[r, 0, t] == u[producer, t], name=f"inflow_boundary[{r},{t}]")

        # set initial values
        for r, _ in enumerate(A):
            # xiPlus[r, x0, t0] = 0
            xiPlus[r, 0, 0].lb = 0.0
            xiPlus[r, 0, 0].ub = 0.0
            # we start with an empty network:
            for x in range(nx):
                xiPlus[r, x, 0].lb = 0.0
                xiPlus[r, x, 0].ub = 0.0
                xiMinus[r, x, 0].lb = 0.0
                xiMinus[r, x, 0].ub = 0.0

        # build the objective expression
        expr = 0.0
        for s, consumer in enumerate(consumers):
            for t in range(self.nt):
                load = sum(xiPlus[r, nx-1, t] for r in delta_in[consumer])
                expr += self.dt * (demand[s, t] - load)**2
        self.miqcqp.setObjective(0.05*expr)

        # Network dynamics constraints (pde discretization)
        for r, _ in enumerate(A):
            for t in range(nt-1):
                for x in range(1, nx):
                    # right traveling characteristic
                    w1 = lambda_p*(dt/dx) * \
                        (xiPlus[r, x, t] - xiPlus[r, x-1, t])
                    w2 = -b11 * xiPlus[r, x, t]*dt - b12*xiMinus[r, x, t]*dt
                    self.miqcqp.addConstr(
                        xiPlus[r, x, t+1] == xiPlus[r, x, t] - w1 + w2, name=f"dynamics_left[{r},{x},{t}]")
                for x in range(nx-1):
                    # left-traveling characteristic
                    w1 = lambda_m * (dt/dx) * \
                        (xiMinus[r, x+1, t] - xiMinus[r, x, t])
                    w2 = -b12 * xiPlus[r, x, t]*dt - b11 * xiMinus[r, x, t]*dt
                    self.miqcqp.addConstr(
                        xiMinus[r, x, t+1] == xiMinus[r, x, t] - w1 + w2, name=f"dynamics_right[{r},{x},{t}]")

        # network coupling constraints (linearized)
        for v in V_inner:
            for t in range(1, nt):
                for r in delta_out[v]:
                    # build the rhs expression
                    rhs = 0.0
                    for c, _ in enumerate(configs):
                        for j in delta_in[v]:
                            if linearize_coupling:
                                self.miqcqp.addConstr(
                                    hPlus[c, j, t] <= b[c, t] * bigM, name=f"coupling1_lin_helper1[{c},{j},{t}]")
                                self.miqcqp.addConstr(
                                    hPlus[c, j, t] >= -1.0*b[c, t] * bigM, name=f"coupling1_lin_helper2[{c}, {j}, {t}]")
                                self.miqcqp.addConstr(
                                    hPlus[c, j, t] >= xiPlus[j, nx-1, t] - (1 - b[c, t]) * bigM, name=f"coupling1_lin_helper3[{c}, {j}, {t}]")
                                self.miqcqp.addConstr(
                                    hPlus[c, j, t] <= xiPlus[j, nx-1, t] + (1 - b[c, t]) * bigM, name=f"coupling1_lin_helper4[{c}, {j}, {t}]")
                                rhs += Dplus[c, r, j] * \
                                    lambda_p * hPlus[c, j, t]
                            else:
                                rhs += Dplus[c, r, j] * lambda_p * \
                                    b[c, t] * xiPlus[j, nx-1, t]
                    # add the constraint to the model
                    self.miqcqp.addConstr(
                        lambda_p * xiPlus[r, 0, t] == rhs, name=f"coupling1[{r}, {t}]")
                for r in delta_in[v]:
                    # build the rhs expression
                    rhs = 0
                    for c, _ in enumerate(configs):
                        for j in delta_out[v]:
                            if linearize_coupling:
                                self.miqcqp.addConstr(
                                    hMinus[c, j, t] <= b[c, t] * bigM, name=f"coupling2_lin_helper1[{c}, {j}, {t}]")
                                self.miqcqp.addConstr(
                                    hMinus[c, j, t] >= -1.0*b[c, t] * bigM, name=f"coupling2_lin_helper2[{c}, {j}, {t}]")
                                self.miqcqp.addConstr(
                                    hMinus[c, j, t] >= xiMinus[j, 0, t] - (1-b[c, t]) * bigM, name=f"coupling2_lin_helper3[{c}, {j}, {t}]")
                                self.miqcqp.addConstr(
                                    hMinus[c, j, t] <= xiMinus[j, 0, t] + (1-b[c, t]) * bigM, name=f"coupling2_lin_helper4[{c}, {j}, {t}]")
                                rhs += Dminus[c, r, j] * \
                                    lambda_m * hMinus[c, j, t]
                            else:
                                rhs += Dminus[c, r, j] * lambda_m * \
                                    b[c, t] * xiMinus[j, 0, t]
                    # add the constraint to the model
                    self.miqcqp.addConstr(
                        lambda_m * xiMinus[r, nx-1, t] == rhs, name=f"coupling2[{r}, {t}]")

        # SOS1 constraint (exactly one active config at each time)
        for t in range(nt):
            self.miqcqp.addConstr(qsum(b[c, t]
                                       for c, _ in enumerate(configs)) == 1, name=f"sos[{t}]")

        # dwell time constraints
        if self.min_up_dwell_times is not None:
            for c, _ in enumerate(configs):
                C_U = self.min_up_dwell_times[c]
                for k in range(nt - C_U):
                    # minimum up dwell time
                    lhs = qsum(b[c, t] for t in range(k+1, k + C_U + 1))
                    self.miqcqp.addConstr(
                        lhs >= C_U * (b[c, k+1] - b[c, k]), name=f"min_dwell_up_time[{c}, {k}]")
        if self.min_down_dwell_times is not None:
            for c, _ in enumerate(configs):
                C_D = self.min_down_dwell_times[c]
                for k in range(nt - C_D):
                    # minimum down dwell time
                    lhs = qsum(1 - b[c, t] for t in range(k+1, k + C_D + 1))
                    self.miqcqp.addConstr(
                        lhs >= C_D * (b[c, k] - b[c, k+1]), name=f"min_dwell_down_time[{c}, {k}]")

        # transmission disturbances
        if self.transmission_disturbances is not None:
            for i, line in enumerate(self.disturbed_lines):
                for x in range(self.nx):
                    lhs = qsum(
                        (xiPlus[line, x, t] - self.xiPlus_mean[i])**2 for t in range(nt))
                    self.miqcqp.addConstr(
                        lhs <= self.load_threshold[i], name=f"trans_disturbances[{i}, {x}]")

    def solve(self, solver="Gurobi", max_calls=50000) -> None:
        self.used_solver = solver
        # build the model
        self.__build_MIQCQP()
        # define a callback for the CIAP heuristic

        def CIAPcallback(model, where):
            if model._calls <= model._maxcalls and where == GRB.Callback.MIPNODE:
                # Get the relaxed control (solution of node relaxation)
                alpha = np.array([val for val in model.cbGetNodeRel(
                    model._binControl)]).reshape(self.noc, self.nt)
                # np.savetxt(f"logging/alpha_{model._calls:d}.txt", alpha)
                # Calculate a feasible control by DSUR
                beta = DSUR(alpha, self.dt, self.timegrid,
                            self.min_up_dwell_times[0], self.min_down_dwell_times[0])
                # Post the feasible control to Gurobi as heuristic solution
                model.cbSetSolution(model._binControl, beta.flatten())
                # Use the posted feasible control and try to compute a heuristic
                # solution (i.e. compute the remaining variable values)
                objval = model.cbUseSolution()
                # print(beta)
                # New Solution found: (objval == Inf, if the solution is infeasible
                # or not a new incumbent)
                if objval < GRB.INFINITY:
                    print(f"C                                 {objval:.4f}")
                model._calls += 1
            if where == GRB.Callback.MIPSOL:
                # Logging callback
                model._incObjVals.append(model.cbGet(GRB.Callback.MIPSOL_OBJ))
                model._incRuntime.append(model.cbGet(GRB.Callback.RUNTIME))
            if where == GRB.Callback.MESSAGE:
                # How have the solution been found?
                if (c := model.cbGet(GRB.Callback.MSG_STRING)[0]) in ['H', '*']:
                    model._incFoundBy += c

        def loggingCallback(model, where) -> None:
            if where == GRB.Callback.MIPSOL:
                # Logging callback
                model._incObjVals.append(model.cbGet(GRB.Callback.MIPSOL_OBJ))
                model._incRuntime.append(model.cbGet(GRB.Callback.RUNTIME))
            if where == GRB.Callback.MESSAGE:
                # How have the solution been found?
                if (c := model.cbGet(GRB.Callback.MSG_STRING)[0]) in ['H', '*']:
                    model._incFoundBy += c

        # Set attributes inside the gurobi model
        self.miqcqp._calls = 0
        self.miqcqp._maxcalls = max_calls
        self.miqcqp._binControl = [
            var for var in self.miqcqp.getVars() if var.vtype == "B"]
        self.miqcqp._incObjVals = []
        self.miqcqp._incFoundBy = ""
        self.miqcqp._incRuntime = []
        # call the solver
        if solver == "Gurobi":
            self.miqcqp.optimize(loggingCallback)
        if solver == "GurobiCIAP":
            self.miqcqp.optimize(CIAPcallback)

    def get_logging_info(self) -> None:
        return {
            "method": self.used_solver,
            "incObjVals": self.miqcqp._incObjVals,
            "incFoundBy": self.miqcqp._incFoundBy,
            "incRuntimes": self.miqcqp._incRuntime,
            "totalRuntime": self.miqcqp.Runtime
        }

    def plot_results(self) -> None:
        time = np.arange(0.0, self.T + self.dt, self.dt)
        # Plot the delivered loads
        fig, axes = plt.subplots(len(self.consumers), 1)
        for i, consumer in enumerate(self.consumers):
            for r in self.delta_in[consumer]:
                delivered_load = np.array(
                    [self.xiPlus[r, self.nx-1, t].X for t in range(self.nt)])
                axes[i].plot(time, self.demand[i, :self.nt],
                             time, delivered_load)
            if i < len(self.consumers) - 1:
                axes[i].set_xticks([])
        axes[-1].set_xlabel(r"$t$")
        fig.suptitle(
            r"Delivered load $\xi^{+}_r(x_{{\ell}_r}, t)$ and consumers demand $Q_r(t)$")
        fig.tight_layout()

        # Plot the outer controls
        fig, axes = plt.subplots(len(self.configs), 1)
        for c, _ in enumerate(self.configs):
            b = np.array([self.binControl[c, t].X for t in range(self.nt)])
            axes[c].step(time, b)
            axes[c].set_ylim([-0.125, 1.25])
            axes[c].set_yticks([0, 1.0])
            if c < len(self.configs) - 1:
                axes[c].set_xticks([])
        axes[-1].set_xlabel(r"$t$")
        fig.suptitle(r"Binary controls $b_c(t)$")
        fig.tight_layout()

        # Plot the inflow controls
        fig, axes = plt.subplots(len(self.producers), 1)
        for q, producer in enumerate(self.producers):
            inflow = np.array([self.xiPlus[producer, 0, t].X
                               for t in range(self.nt)])
            axes[q].step(time, inflow)
            if q < len(self.producers) - 1:
                axes[q].set_xticks([])
        axes[-1].set_xlabel(r"$t$")
        fig.suptitle(r"Inflow controls $u_q(t)$")
        fig.tight_layout()

        # plot the dynamics on the lines
        xiP, xiM, b = self.get_sols()

        fig, axes = plt.subplots(ceil(xiP.shape[0] // 2), 2)

        for r, ax in enumerate(axes.flatten()):
            if r < len(self.A) - 2:
                ax.set_xticks([])
            else:
                ax.set_xlabel(r"$t$")
            im = ax.pcolormesh(xiP[r], cmap="viridis",
                               label=r"$\xi_" + f"{r}" + r"$")
            # no y label for right column and only colorbars for the right col
            if (r + 1) % 2 == 0:
                ax.set_yticks([])
                fig.colorbar(im, ax=ax)
            else:
                ax.set_ylabel(r"$x$")
        fig.suptitle(
            r"Transmission characteristics $\xi^{+}_r$ on all network lines $r \in A$")
        plt.show()

    def getGrbModel(self, linearize_coupling=True) -> None:
        self.__build_MIQCQP()
        self.miqcqp.update()
        return self.miqcqp

    def get_sols(self):
        b = np.array([self.binControl[c, t].x for c, _ in enumerate(self.configs)
                      for t in range(self.nt)]).reshape(len(self.configs), self.nt)
        xiP = np.array([self.xiPlus[r, x, t].X for r, _ in enumerate(self.A) for x in range(
            self.nx) for t in range(self.nt)]).reshape(len(self.A), self.nx, self.nt)
        xiM = np.array([self.xiMinus[r, x, t].X for r, _ in enumerate(self.A) for x in range(
            self.nx) for t in range(self.nt)]).reshape(len(self.A), self.nx, self.nt)
        return xiP, xiM, b
