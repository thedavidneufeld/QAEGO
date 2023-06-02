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
        # Make a copy of the old data to retain its data
        data = copy.deepcopy(olddata)
        P, Q, V, PA, Q_min, Q_max = self._load_powerdata(data)
        V0 = V
        PA0 = PA
        # Value calculations
        LA = self._calculate_LA(data)
        Y, G, B = self._create_Ybus(data, LA)
        self._init_unknown(data, V, PA)
        total_iterations = 0 # Total iterations of NR
        iterations = 0 # Iterations since limit check last failed
        # Perform NR method iteratively
        # Stop condition handled within the loop
        while(True):
            P_new, Q_new = self._compute_powers(data, V, PA, G, B)
            # If limit check fails, treat ith PV bus as PQ bus and recompute powers
            while(not self._check_limits(data, Q_new, Q_min, Q_max)):
                iterations = 0
                P, Q, V, PA, Q_min, Q_max = self._load_powerdata(data)
                self._init_unknown(data, V, PA)
                P_new, Q_new = self._compute_powers(data, V, PA, G, B)
            Delta_PQ = self._compute_power_mismatches(data, P, Q, P_new, Q_new)
            Jac = self._compute_Jacobian(data, V, PA, G, B, P_new, Q_new)
            Delta_PAV, Delta_PAV_classic = self._compute_increment(Delta_PQ, Jac, ep)
            PA_new, V_new = self._compute_new_PAV(data, V, PA, Delta_PAV)
            total_iterations += 1
            iterations += 1
            if self.solve_method == 'hhl':
                # Compare HHL with Classical Method
                print("Delta_PAV difference at iteration ", total_iterations, ":\n", np.linalg.norm(Delta_PAV_classic - Delta_PAV), "\n")
            # If convergence has been reached, exit loop
            if self._check_convergence(Delta_PAV, ep):
                break
            V = V_new
            PA = PA_new
        self._load_values(data, olddata, V_new, PA_new)

        # Print Data
        print("\nTotal nmber of iterations: ", total_iterations)
        print("Iterations since last limit check failure: ", iterations)
        print("\nOld Voltages:\n", V0)
        print("\nNew Voltages:\n", V_new)    
        print("\nOld Phase Angles:\n", PA0)
        print("\nNew Phase Angles:\n", PA_new)
        print("\nOld Bus Values:\n", olddata.buses, "\n")
        print("New Bus Values (after NR):\n", data.buses, "\n")

        # Return new data with updated values
        return data
    
    def NR_PC(self, olddata:PowerData, ep):
        # Make a copy of the old data to retain its data
        data = copy.deepcopy(olddata)
        P, Q, V, PA, Q_min, Q_max = self._load_powerdata(data)
        V0 = V
        PA0 = PA
        # For first iteration, set V_avg and PA_avg equal to V and PA
        V_avg = V
        PA_avg = PA
        # Value calculations
        LA = self._calculate_LA(data)
        Y, G, B = self._create_Ybus(data, LA)
        self._init_unknown(data, V, PA)
        total_iterations = 0 # Total iterations of NR
        iterations = 0 # Iterations since limit check last failed
        # Perform NR method iteratively
        # Stop condition handled within the loop
        while(True):
            P_new, Q_new = self._compute_powers(data, V, PA, G, B)
            # If limit check fails, treat ith PV bus as PQ bus and recompute powers
            while(not self._check_limits(data, Q_new, Q_min, Q_max)):
                iterations = 0
                P, Q, V, PA, Q_min, Q_max = self._load_powerdata(data)
                self._init_unknown(data, V, PA)
                P_new, Q_new = self._compute_powers(data)
            Delta_PQ = self._compute_power_mismatches(data, P, Q, P_new, Q_new)
            if iterations > 0:
                # Predicted values
                Delta_PAV_pred, x = self._compute_increment(Delta_PQ, Jac, ep)
                PA_pred, V_pred = self._compute_new_PAV(data, V, PA, Delta_PAV_pred)
                V_avg = 0.5*(V + V_pred)
                PA_avg = 0.5*(PA + PA_pred)
            # Corrected values
            Jac = self._compute_Jacobian(data, V_avg, PA_avg, G, B, P_new, Q_new)
            Delta_PAV, Delta_PAV_classic = self._compute_increment(Delta_PQ, Jac, ep)
            PA_new, V_new = self._compute_new_PAV(data, V, PA, Delta_PAV)
            total_iterations += 1
            iterations += 1
            if self.solve_method == 'hhl':
                # Compare HHL with Classical Method
                print("Delta_PAV difference at iteration ", total_iterations, ":\n", np.linalg.norm(Delta_PAV_classic - Delta_PAV), "\n")
            # If convergence has been reached, exit loop
            if self._check_convergence(Delta_PAV, ep):
                break
            V = V_new
            PA = PA_new
        self._load_values(data, olddata, V_new, PA_new)

        # Print Data
        print("\nTotal nmber of iterations: ", total_iterations)
        print("Iterations since last limit check failure: ", iterations)
        print("\nOld Voltages:\n", V0)
        print("\nNew Voltages:\n", V_new)    
        print("\nOld Phase Angles:\n", PA0)
        print("\nNew Phase Angles:\n", PA_new)
        print("\nOld Bus Values:\n", olddata.buses, "\n")
        print("New Bus Values (after NR):\n", data.buses, "\n")

        # Return new data with updated values
        return data
                
    def _load_powerdata(self, data:PowerData):
        P = np.zeros(0)
        Q = np.zeros(0)
        V = np.zeros(0)
        PA = np.zeros(0)
        Q_min = np.zeros(0)
        Q_max = np.zeros(0)
        
        for b in data.buses:
            # Add PA value from bus
            PA = np.append(PA, b[4])
            # Check if bus is PQ or not
            if b[5] != 3:
                # Check if a generator is attatched to the bus
                has_gen = data.gens[:, 0] == b[0]
                if has_gen.any():
                    bus_gen = data.gens[has_gen]
                    # Add generator value for bus as its voltage magnitude
                    V = np.append(V, bus_gen[0][3])
                    # Add P and Q values
                    if b[5] == 2:
                        P = np.append(P, bus_gen[0][1])
                        Q = np.append(Q, bus_gen[0][2])
                        Q_min = np.append(Q_min, bus_gen[0][4])
                        Q_max = np.append(Q_max, bus_gen[0][5])
                    else:
                        P = np.append(P, -b[1])
                        Q = np.append(Q, -b[2])
                        Q_min = np.append(Q_min, 0)
                        Q_max = np.append(Q_max, 0)
                else:
                    raise Exception('All Slack and PV buses must have a gen! At least one gen missing.')
            else:
                # Add P, Q, and V values from bus
                P = np.append(P, -b[1])
                Q = np.append(Q, -b[2])
                V = np.append(V, b[3])
                Q_min = np.append(Q_min, 0)
                Q_max = np.append(Q_max, 0)

        return P, Q, V, PA, Q_min, Q_max

    def _calculate_LA(self, data:PowerData):
        # Line Admittance Matrix
        LA = np.zeros((data.buses.shape[0], data.buses.shape[0]), dtype = complex)
        for b in data.branches:
            # Branch b[0] and b[1] contain fbus and tbus information
            i = int(b[0].real)-1
            j = int(b[1].real)-1
            # Set LA for both buses to 1/LI
            LA[i][j] = 1/b[2]
            LA[j][i] = 1/b[2]

        return LA

    def _create_Ybus(self, data:PowerData, LA):
        Y = np.zeros((data.buses.shape[0], data.buses.shape[0]),  dtype = complex)
        for i in range(Y.shape[0]):
            for j in range (Y.shape[1]):
                if i != j:
                    Y[i, j] = Y[i, j] - LA[i][j]
                    Y[i, i] = Y[i, i] + LA[i][j]
        G = Y.real
        B = Y.imag

        return Y, G, B

    def _init_unknown(self, data, V, PA):
        # Flat start for bus voltages and phase angles
        # Unknown values for bus volage will be set to 1
        # Unknown values for phase angles will be set to 0
        for i in range(1, data.buses.shape[0]):
            if (data.buses[i][5] == 2) or (data.buses[i][5] == 3):
                PA[i] = 0
                if data.buses[i][5] == 3:
                    V[i] = 1

    def _compute_powers(self, data:PowerData, V, PA, G, B):
        P_new = np.zeros(data.buses.shape[0])
        Q_new = np.zeros(data.buses.shape[0])
        # We don't need to calculate for the slack bus (the first bus), so we start at the second bus
        for i in range(1, data.buses.shape[0]):
            for j in range (data.buses.shape[0]):
                P_new[i] = P_new[i]+V[i]*V[j]*(G[i][j]*np.cos(PA[i]-PA[j])+B[i][j]*np.sin(PA[i]-PA[j]))
                Q_new[i] = Q_new[i]+V[i]*V[j]*(G[i][j]*np.sin(PA[i]-PA[j])-B[i][j]*np.cos(PA[i]-PA[j]))

        return P_new, Q_new

    def _check_limits(self, data:PowerData, Q_new, Q_min, Q_max):
        for i in range(1, data.buses.shape[0]):
            if data.buses[i][5] == 2:
                # Check that Q for gen i is within limits
                if (Q_new[i] > Q_max[i]) or (Q_new[i] < Q_min[i]):
                    # Treat bus as PQ
                    data.buses[i][5] = 3
                    return False
        # If the for loop is able to complete, all limit checks pass
        return True
    
    def _compute_power_mismatches(self, data:PowerData, P, Q, P_new, Q_new):
        Delta_P = np.zeros(0)
        Delta_Q = np.zeros(0)
        for i in range(1, data.buses.shape[0]):
            Delta_P = np.append(Delta_P, P[i] - P_new[i])
            if data.buses[i][5] == 3:
                Delta_Q = np.append(Delta_Q, Q[i] - Q_new[i])
        Delta_PQ = np.append(Delta_P, Delta_Q)

        return Delta_PQ

    def _compute_Jacobian(self, data:PowerData, V, PA, G, B, P_new, Q_new):
        H = np.zeros((data.buses.shape[0]-1, data.buses.shape[0]-1))
        N = np.zeros((data.buses.shape[0]-1, data.buses.shape[0]-1))
        J = np.zeros((data.buses.shape[0]-1, data.buses.shape[0]-1))
        L = np.zeros((data.buses.shape[0]-1, data.buses.shape[0]-1))
        for i in range(1, data.buses.shape[0]):
            for j in range(1, data.buses.shape[0]):
                if i == j:
                    H[i-1][j-1] = -Q_new[i]-B[i][i]*(V[i]**2)
                    # data.buses[i][5] is Bus i's BusType, 2 is PV
                    if data.buses[i][5] != 2:
                        N[i-1][j-1] = P_new[i]+G[i][i]*(V[i]**2)
                        J[i-1][j-1] = P_new[i]-G[i][i]*(V[i]**2)
                        L[i-1][j-1] = Q_new[i]-B[i][i]*(V[i]**2)
                else:
                    H[i-1][j-1] = V[i]*V[j]*(G[i][j]*np.cos(PA[i]-PA[j])+B[i][j]*np.sin(PA[i]-PA[j]))
                    if data.buses[j][5] != 2:
                        N[i-1][j-1] = V[j]*V[i]*(G[i][j]*np.sin(PA[i]-PA[j])+B[i][j]*np.cos(PA[i]-PA[j]))
                    if data.buses[i][5] != 2:
                        J[i-1][j-1] = -V[i]*V[j]*(G[i][j]*np.cos(PA[i]-PA[j])+B[i][j]*np.sin(PA[i]-PA[j]))
                    if data.buses[i][5] != 2 and data.buses[j][5] != 2:
                        L[i-1][j-1] = V[j]*V[i]*(G[i][j]*np.sin(PA[i]-PA[j])-B[i][j]*np.cos(PA[i]-PA[j]))

        for i in reversed(range(1, data.buses.shape[0])):
            if data.buses[i][5] == 2:
                # If Bus i is PV, remove column for that bus from N
                N = np.delete(N, i-1, 1)
                # If Bus i is PV, remove row for that bus from J
                J = np.delete(J, i-1, 0)
                # If Bus i is PV, remove row and column for that bus from L
                L = np.delete(L, i-1, 1)
                L = np.delete(L, i-1, 0)

        Jac = np.block([[H, N], [J, L]])

        return Jac

    def _compute_increment(self, Delta_PQ, Jac, ep):
        if self.solve_method == 'hhl':
            hhl = hhl_helper()
            Delta_PAV = hhl.run_HHL(Jac, Delta_PQ, ep)
            Delta_PAV_classic = np.dot(inv(Jac), Delta_PQ)
        else:
            Delta_PAV = np.dot(inv(Jac), Delta_PQ)
            Delta_PAV_classic = Delta_PAV

        return Delta_PAV, Delta_PAV_classic

    def _compute_new_PAV(self, data:PowerData, V, PA, Delta_PAV):
        PA_new = np.zeros(data.buses.shape[0])
        V_new = np.zeros(data.buses.shape[0])
        PA_new[0] = PA[0]
        V_new[0] = V[0]
        num_buses = data.buses.shape[0]
        incr = 0
        for i in range(1, num_buses):
            PA_new[i] = PA[i] + Delta_PAV[i-1]
            if data.buses[i][5] == 3:
                V_new[i] = V[i] + (Delta_PAV[num_buses-1+incr])*V[i]
                incr += 1
            else:
                V_new[i] = V[i]
            
        return PA_new, V_new

    def _check_convergence(self, Delta_PAV, ep):
        if max(Delta_PAV) > ep:
            return False
        else:
            return True
    
    def _load_values(self, data:PowerData, olddata:PowerData, V_new, PA_new):
        # Load calculated values into data.buses and change altered PV buses back to PV
        for i in range(data.buses.shape[0]):
            data.buses[i][4] = PA_new[i]
            data.buses[i][3] = V_new[i]
            if (data.buses[i][5] == 3) and (olddata.buses[i][5] == 2):
                data.buses[i][5] = 2