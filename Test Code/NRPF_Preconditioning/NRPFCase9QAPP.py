from PowerFlow import PowerData, NRPF
import matplotlib.pyplot as plt
from pathlib import Path
import time

def main():
    # Set directory
    dir = Path(__file__).parent.resolve()

    # Load case9Q data
    case9Q = PowerData()
    case9Q.loadcasedata(str(dir) + '/case9Q_simplified.m')

    # Print initial data
    print('Case 9Q:\n')
    print('Buses:\n', case9Q.buses, '\n')
    print('Gens\n', case9Q.gens, '\n')
    print('Branches\n', case9Q.branches, '\n')
    
    # Initialize Classical solver
    pf = NRPF()
    # Initialize HHL solver
    pf_hhl = NRPF('hhl')

    # Start timing Classical solve
    st_cl = time.time()

    # Solve Classically
    c9Q_pf = pf.NRAPP(case9Q, 1e-8)

    # Stop timing Classical solve and output time
    et_cl = time.time()
    cl_time = et_cl - st_cl
    print('Time to solve (Classical with APP):', cl_time, 'seconds\n')

    # Print Results
    print('Newton Raphson (APP) Solution:\n')
    print('\nIterations:', c9Q_pf['iterations'])
    print('Voltages each iteration:\n', c9Q_pf['voltages'])
    print('Phase Angles each iteration:\n', c9Q_pf['phase_angles'], '\n')

    # Save Plot
    plt.plot(c9Q_pf["dpav_difference"])
    plt.savefig("c9Q_pf_dpav_difference.png")
    plt.clf()

    # Start timing HHL solve
    st_hhl = time.time()

    # Solve with HHL
    c9Q_hhl = pf_hhl.NRAPP(case9Q, 1e-8)

    # Stop timing HHL solve and output time
    et_hhl = time.time()
    hhl_time = et_hhl - st_hhl
    print('Time to solve (HHL with APP):', hhl_time, 'seconds\n')

    # Print Results
    print('Newton Raphson (APP) with HHL Solution:\n')
    print('\nIterations:', c9Q_hhl['iterations'])
    print('Voltages each iteration:\n', c9Q_hhl['voltages'])
    print('Phase Angles each iteration:\n', c9Q_hhl['phase_angles'], '\n')

    # Save Plot
    plt.plot(c9Q_hhl["dpav_difference"])
    plt.savefig("c9Q_hhl_dpav_difference.png")

if __name__ == "__main__":
    main()