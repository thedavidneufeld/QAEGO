# QAEGO (Quantum Algorithms for Electrical Grid Optimization) Research
# University of Lethbridge
# David Neufeld 
# May 2023

import copy
import numpy as np
from numpy.linalg import inv
from HHL_conversion import hhl_helper


# Class for keeping bus, gen, and branch data
class PowerData:

    def __init__(self):
        self.buses = np.zeros((0, 6))
        self.gens = np.zeros((0, 6))
        self.branches = np.zeros((0, 3), dtype=complex)

    # pd is active load power
    # qd is reactive load power
    # v is voltage magnitude
    # pa is voltage phase angle
    # bustype is either 1:Slack, 2:PV, or 3(or any other value):PQ
    def addBus(self, busnum:int, pd:float, qd:float, v:float, pa:float, bustype:int):
        # Check if busnum already exists in buses
        in_buses = self.buses[:, 0] == busnum
        if in_buses.any():
            print('Bus', busnum, 'already exists!\n')
            print('Buses:\n', self.buses, '\n')
        else:
            # Check for valid bustype, if the value is anything 
            # other than 1, 2, or 3, set to 3
            if bustype < 1 or bustype > 3:
                bustype = 3
            # Add bus to buses
            newbus = np.array([busnum, pd, qd, v, pa, bustype])
            self.buses = np.vstack((self.buses, newbus))
            # Sort buses based on busnums
            self.buses = self.buses[self.buses[:, 0].argsort()]

    def removeBus(self, busnum:int):
        # Remove bus with busnum if it exists
        in_buses = self.buses[:, 0] != busnum
        if not in_buses.all():
            self.buses = self.buses[in_buses]
            # Remove any gens that might exist at busnum
            self.removeGen(busnum)
    
    # pg is injected active power
    # qg is injected reactive power
    # v is voltage magnitude
    # qmin is minimum value for qg
    # qmax is maximum value for qg
    def addGen(self, busnum:int, pg:float, qg:float, v:float, qmin:float, qmax:float):
        # Add to gen to gens if a bus exists at busnum
        in_buses = self.buses[:, 0] == busnum
        if in_buses.any():
            newgen = np.array([busnum, pg, qg, v, qmin, qmax])
            self.gens = np.vstack((self.gens, newgen))
            # Sort gens based on busnums
            self.gens = self.gens[self.gens[:, 0].argsort()]
        else:
            print('A gen cannot be added at bus', busnum, 'because that bus does not exist!\n')
            print('Buses:\n', self.buses, '\n')

    def removeGen(self, busnum:int):
        # Remove gen(s) with busnum if it/they exists
        in_gens = self.gens[:, 0] != busnum
        if not in_gens.all():
            self.gens = self.gens[in_gens]
    
    # fbus and tbus are the buses the branch is going between
    # LI is line impedance of the branch
    def addBranch(self, fbus:int, tbus:int, LI:complex):
        # Do not add a branch when fbus and tbus are the same
        if fbus == tbus:
            print('A branch cannot exist between the same bus!\n')
            return
        # Add branch to branches if fbus and tbus exist
        fin_bus = self.buses[:, 0] == fbus
        tin_bus = self.buses[:, 0] == tbus
        if fin_bus.any() and tin_bus.any():
            # If branch already exists, do not add it
            newbranch = np.array([fbus, tbus, LI])
            fin_branches = np.empty(0, dtype=bool)
            tin_branches = np.empty(0, dtype=bool)
            for b in self.branches:
                f = b[0:2] == fbus
                t = b[0:2] == tbus
                fin_branches = np.append(fin_branches, f.any())
                tin_branches = np.append(tin_branches, t.any())
            in_branches = fin_branches & tin_branches
            if in_branches.any():
                print('Branch from', fbus, 'to', tbus, 'already exists!\n')
                print('Branches:\n', self.branches, '\n')
                return
            self.branches = np.vstack((self.branches, newbranch))
            # Sort branches based on fbus and tbus
            self.branches = self.branches[self.branches[:, 1].argsort()]
            self.branches = self.branches[self.branches[:, 0].argsort()]
        else:
            print('At least one of buses', fbus, 'and', tbus, 'does not exist!')
            print('Buses:\n', self.buses, '\n')

    def removeBranch(self, fbus:int, tbus:int=-1):
        # If tbus == -1, remove all branches with an fbus or tbus value 
        # equivalent to the specified fbus value if they exist
        if tbus == -1:
            fin_bus = self.branches[:, 0] != fbus
            tin_bus = self.branches[:, 1] != fbus
            in_branches = fin_bus & tin_bus
            if not in_branches.all():
                self.branches = self.branches[in_branches]
        # Remove all branches with fbus and tbus if they exist
        else:
            fin_bus = self.branches[:, 0] != fbus
            tin_bus = self.branches[:, 1] != tbus
            in_branches = fin_bus | tin_bus
            if not in_branches.all():
                self.branches = self.branches[in_branches]
    

# Class for solving power flow using Newton-Raphson method
class NRPF:
    def __init__(self, method:str=None):
        self.set_solve_method(method)
    
    def set_solve_method(self, method:str):
        if method == 'classical' or method == 'hhl':
            self.solve_method = method
        elif method == None:
            self.solve_method = 'classical'
        else:
            self.solve_method = 'classical'
            print('Not a solve method, set to default (classical)')

    def NR(self, olddata:PowerData, ep):
            self._reset_values
            # Make a copy of the old data to retain its data
            data = copy.deepcopy(olddata)
            self._load_powerdata(data)
            # Value calculations
            self._create_GA_LA(data)
            self._create_Ybus(data)
            self._init_unknown()
            total_iterations = 0 # Total iterations of NR
            iterations = 0 # Iterations since limit check last failed
            # Perform NR method iteratively
            # Stop condition handled within the loop
            while(True):
                self._compute_powers(data)
                # If limit check fails, treat ith PV bus as PQ bus and recompute powers
                if not self._check_limits(data):
                    iterations = 0
                    self._compute_powers(data)
                self._compute_power_mismatches(data)
                self._compute_Jacobian(data)
                self._compute_increment(ep)
                self._compute_new_PAV(data)
                total_iterations += 1
                iterations += 1
                if self.solve_method == 'hhl':
                    # Compare HHL with Classical Method
                    print("Delta_PA difference at iteration ", total_iterations, ":\n", np.linalg.norm(self.Delta_PA_classic - self.Delta_PA), "\n")
                    print("Delta_V difference at iteration ", total_iterations, ":\n", np.linalg.norm(self.Delta_V_classic - self.Delta_V), "\n")
                # If convergence has been reached, exit loop
                if self._check_convergence(data, ep):
                    break
                self.V = self.V_new
                self.PA = self.PA_new
            self._load_values(data, olddata)
            self._calculate_injected(data)

            # Print Data
            print("\nTotal nmber of iterations: ", total_iterations)
            print("Iterations since last limit check failure: ", iterations)
            print("\nOld Voltages:\n", olddata.buses[:,3])
            print("\nNew Voltages:\n", data.buses[:,3])    
            print("\nOld Phase Angles:\n", olddata.buses[:,4])
            print("\nNew Phase Angles:\n", data.buses[:,4])
            print("\nOld Bus Values:\n", olddata.buses, "\n")
            print("New Bus Values (after NR):\n", data.buses, "\n")
            print("Old Gen Values:\n", olddata.gens, "\n")
            print("New Gen Values (after NR):\n", data.gens, "\n")

            # Return new data with updated values
            return data
    
    def _reset_values(self):
        self.P = None
        self.Q = None
        self.V = None
        self.PA = None
        self.Q_min = None
        self.Q_max = None
        self.LA = None
        self.GA = None
        self.Y = None
        self.G = None
        self.B = None
        self.P_new = None
        self.Q_new = None
        self.Delta_P = None
        self.Delta_Q = None
        self.Delta_PQ = None
        self.Delta_PAV = None
        self.Delta_PA = None
        self.Delta_PA_classic = None
        self.Delta_V = None
        self.Delta_V_classic = None
        self.PA_new = None
        self.V_new = None

    def _load_powerdata(self, data:PowerData):
        self.P = np.zeros(0)
        self.Q = np.zeros(0)
        self.V = np.zeros(0)
        self.PA = np.zeros(0)
        self.Q_min = np.zeros(0)
        self.Q_max = np.zeros(0)
        
        for b in data.buses:
            np.append(self.PA, b[4])
            # Check if bus is PQ or not
            if b[5] != 3:
                # Check if a generator is attatched to the bus
                has_gen = data.gens[:, 0] == b[0]
                if has_gen.any():
                    bus_gen = data.gens[has_gen]
                    # Add generator value for bus as its voltage magnitude
                    self.V = np.append(self.V, bus_gen[0][3])
                    # Add P and Q values
                    if b[5] == 2:
                        self.P = np.append(self.P, bus_gen[0][1])
                        self.Q = np.append(self.Q, bus_gen[0][2])
                        self.Q_min = np.append(self.Q_min, bus_gen[0][4])
                        self.Q_max = np.append(self.Q_max, bus_gen[0][5])
                    else:
                        self.P = np.append(self.P, -b[1])
                        self.Q = np.append(self.Q, -b[2])
                        self.Q_min = np.append(self.Q_min, 0)
                        self.Q_max = np.append(self.Q_max, 0)
                else:
                    raise Exception('All Slack and PV buses must have a gen! At least one gen missing.')
            else:
                # Add P, Q, and V values from bus
                self.P = np.append(self.P, -b[1])
                self.Q = np.append(self.Q, -b[2])
                self.V = np.append(self.V, b[3])
                self.Q_min = np.append(self.Q_min, 0)
                self.Q_max = np.append(self.Q_max, 0)
            # Add PA value from bus
            self.PA = np.append(self.PA, b[4])

    def _create_GA_LA(self, data:PowerData):
        # Line Admittance Matrix
        self.LA = np.zeros((data.buses.shape[0], data.buses.shape[0]), dtype = complex)
        for b in data.branches:
            # Branch b[0] and b[1] contain fbus and tbus information
            i = int(b[0].real)-1
            j = int(b[1].real)-1
            # Set LA for both buses to 1/branch
            self.LA[i][j] = 1/b[2]
            self.LA[j][i] = 1/b[2]

        # Compute Generator Admittance
        self.GA = np.zeros(data.buses.shape[0])

    def _create_Ybus(self, data:PowerData):
        self.Y = np.zeros((data.buses.shape[0], data.buses.shape[0]),  dtype = complex)
        for i in range(self.Y.shape[0]):
            self.Y[i][i] = self.Y[i][i] + self.GA[i]
            for j in range (self.Y.shape[1]):
                if i != j:
                    self.Y[i, j] = self.Y[i, j] - self.LA[i][j]
                    self.Y[i, i] = self.Y[i, i] + self.LA[i][j]
                    
        self.G = self.Y.real
        self.B = self.Y.imag

    def _init_unknown(self):
        # Flat start for bus voltages
        # Unknown values for bus volage will be set to 1
        for i in range(self.V.size):
            if self.V[i] == 0:
                self.V[i] = 1

    def _compute_powers(self, data:PowerData):
        self.P_new = np.zeros(data.buses.shape[0])
        self.Q_new = np.zeros(data.buses.shape[0])
        # We don't need to calculate for the slack bus (the first bus), so we start at the second bus
        for i in range(1, data.buses.shape[0]):
            for j in range (data.buses.shape[0]):
                self.P_new[i] = self.P_new[i]+self.V[i]*self.V[j]*(self.G[i][j]*np.cos(self.PA[i]-self.PA[j])+self.B[i][j]*np.sin(self.PA[i]-self.PA[j]))
                self.Q_new[i] = self.Q_new[i]+self.V[i]*self.V[j]*(self.G[i][j]*np.sin(self.PA[i]-self.PA[j])-self.B[i][j]*np.cos(self.PA[i]-self.PA[j]))

    def _check_limits(self, data:PowerData):
        for i in range(1, data.buses.shape[0]):
            if data.buses[i][5] == 2:
                # Check that Q for gen i is within limits
                if (self.Q_new[i] > self.Q_max[i]) or (self.Q_new[i] < self.Q_min[i]):
                    # Treat bus as PQ
                    data.buses[i][5] = 3
                    # Set bus's V to initial value
                    data.buses[i][3] = 1
                    self.V[i] = 1
                    return False
        # If the for loop is able to complete, all limit checks pass
        return True
    
    def _compute_power_mismatches(self, data:PowerData):
        self.Delta_P = np.zeros(0)
        self.Delta_Q = np.zeros(0)
        for i in range(1, data.buses.shape[0]):
            self.Delta_P = np.append(self.Delta_P, self.P[i] - self.P_new[i])
            if data.buses[i][5] == 3:
                self.Delta_Q = np.append(self.Delta_Q, self.Q[i] - self.Q_new[i])

    def _compute_Jacobian(self, data:PowerData):
        H = np.zeros((data.buses.shape[0]-1, data.buses.shape[0]-1))
        N = np.zeros((data.buses.shape[0]-1, data.buses.shape[0]-1))
        J = np.zeros((data.buses.shape[0]-1, data.buses.shape[0]-1))
        L = np.zeros((data.buses.shape[0]-1, data.buses.shape[0]-1))
        for i in range(1, data.buses.shape[0]):
            for j in range(1, data.buses.shape[0]):
                if i == j:
                    H[i-1][j-1] = -self.Q_new[i]-self.B[i][i]*(self.V[i]**2)
                    # data.buses[i][7] is Bus i's BusType, 2 is PV
                    if data.buses[i][5] != 2:
                        N[i-1][j-1] = self.P_new[i]+self.G[i][i]*(self.V[i]**2)
                        J[i-1][j-1] = self.P_new[i]-self.G[i][i]*(self.V[i]**2)
                        L[i-1][j-1] = self.Q_new[i]-self.B[i][i]*(self.V[i]**2)
                else:
                    H[i-1][j-1] = self.V[i]*self.V[j]*(self.G[i][j]*np.cos(self.PA[i]-self.PA[j])+self.B[i][j]*np.sin(self.PA[i]-self.PA[j]))
                    if data.buses[j][5] != 2:
                        N[i-1][j-1] = self.V[j]*self.V[i]*(self.G[i][j]*np.sin(self.PA[i]-self.PA[j])+self.B[i][j]*np.cos(self.PA[i]-self.PA[j]))
                    if data.buses[i][5] != 2:
                        J[i-1][j-1] = -self.V[i]*self.V[j]*(self.G[i][j]*np.cos(self.PA[i]-self.PA[j])+self.B[i][j]*np.sin(self.PA[i]-self.PA[j]))
                    if data.buses[i][5] != 2 and data.buses[j][5] != 2:
                        L[i-1][j-1] = self.V[j]*self.V[i]*(self.G[i][j]*np.sin(self.PA[i]-self.PA[j])-self.B[i][j]*np.cos(self.PA[i]-self.PA[j]))

        for i in reversed(range(1, data.buses.shape[0])):
            if data.buses[i][5] == 2:
                # If Bus i is PV, remove column for that bus from N
                N = np.delete(N, i-1, 1)
                # If Bus i is PV, remove row for that bus from J
                J = np.delete(J, i-1, 0)
                # If Bus i is PV, remove row and column for that bus from L
                L = np.delete(L, i-1, 1)
                L = np.delete(L, i-1, 0)

        self.Jac = np.block([[H, N], [J, L]])

    def _compute_increment(self, ep):
        self.Delta_PQ = np.append(self.Delta_P, self.Delta_Q)
        if self.solve_method == 'hhl':
            hhl = hhl_helper()
            self.Delta_PAV = hhl.run_HHL(self.Jac, self.Delta_PQ, ep)
            self.Delta_PAV_classic = np.dot(inv(self.Jac), self.Delta_PQ)
            self.Delta_PA_classic = self.Delta_PAV_classic[0:len(self.Delta_P)]
            self.Delta_V_classic = self.Delta_PAV_classic[len(self.Delta_P):]
        else:
            self.Delta_PAV = np.dot(inv(self.Jac), self.Delta_PQ)
        self.Delta_PA = self.Delta_PAV[0:len(self.Delta_P)]
        self.Delta_V = self.Delta_PAV[len(self.Delta_P):]

    def _compute_new_PAV(self, data:PowerData):
        self.PA_new = np.zeros(data.buses.shape[0])
        self.V_new = np.zeros(data.buses.shape[0])
        self.PA_new[0] = self.PA[0]
        self.V_new[0] = self.V[0]
        incr = 0
        for i in range(1, data.buses.shape[0]):
            self.PA_new[i] = self.PA[i] + self.Delta_PA[i-1]
            if data.buses[i][5] == 3:
                self.V_new[i] = self.V[i] + (self.Delta_V[incr])*self.V[i]
                incr = incr + 1
            else:
                self.V_new[i] = self.V[i]

    def _check_convergence(self, data:PowerData, ep):
        convergence_PA = np.zeros(data.buses.shape[0])
        convergence_V = np.zeros(data.buses.shape[0])
        for i in range(data.buses.shape[0]):
            convergence_PA[i] = abs(self.PA_new[i] - self.PA[i])
            convergence_V[i] = abs(self.V_new[i] - self.V[i])
        for i in range(data.buses.shape[0]):
            if (convergence_PA[i] > ep):
                return False
            # data.buses[i][7] is Bus i's BusType, 3 is PQ
            if data.buses[i][5] == 3:
                if(convergence_V[i] > ep):
                    return False
        # If this line is reached, all convergence levels are <= ep
        return True
    
    def _load_values(self, data:PowerData, olddata:PowerData):
        # Load calculated values into data.buses and change altered PV buses back to PV
        for i in range(data.buses.shape[0]):
            data.buses[i][4] = self.PA_new[i]
            data.buses[i][3] = self.V_new[i]
            if (data.buses[i][5] == 3) and (olddata.buses[i][5] == 2):
                data.buses[i][5] = 2

    def _calculate_injected(self, data:PowerData):
        # for loop temp modification
        for i in range(1, data.buses.shape[0]):
            if data.buses[i][5] == 2:
                for j in range(data.buses.shape[0]):
                    data.gens[i][2] = data.gens[i][2]+data.buses[i][3]*data.buses[j][3]*(self.G[i][j]*np.sin(data.buses[i][4]-data.buses[j][4])-self.B[i][j]*np.cos(data.buses[i][4]-data.buses[j][4]))
