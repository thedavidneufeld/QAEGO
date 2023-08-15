from PowerFlow import PowerData, NRPF
import matplotlib.pyplot as plt
from pathlib import Path
import time

def main():
    # Set directory
    dir = Path(__file__).parent.resolve()

    # Load case3 data
    case3 = PowerData()
    case3.loadcasedata(str(dir) + '/case3.m')

    # Print initial data
    print('Case 3:\n')
    print('Buses:\n', case3.buses, '\n')
    print('Gens\n', case3.gens, '\n')
    print('Branches\n', case3.branches, '\n')
    
    # Initialize Classical solver
    pf = NRPF()
    # Initialize HHL solver
    pf_hhl = NRPF('hhl')

    # Start timing Classical solve
    st_cl = time.time()

    # Solve Classically
    c3_pf = pf.NR(case3, 1e-8)

    # Stop timing Classical solve and output time
    et_cl = time.time()
    cl_time = et_cl - st_cl
    print('Time to solve (Classical):', cl_time, 'seconds\n')

    # Print Results
    print('Newton Raphson Solution:\n')
    for i in range(c3_pf["iterations"]+1):
        print(f"Iteration {i}:")
        print("    Voltages:", c3_pf["voltages"][i])
        print("    Phase Angles:", c3_pf["phase_angles"][i])
    print('\n')

    # Save Plot
    plt.plot(c3_pf["dpav_difference"])
    plt.savefig("c3_pf_dpav_difference.png")
    plt.clf()

    # Start timing HHL solve
    st_hhl = time.time()

    # Solve with HHL
    c3_hhl = pf_hhl.NR(case3, 1e-8)

    # Stop timing HHL solve and output time
    et_hhl = time.time()
    hhl_time = et_hhl - st_hhl
    print('Time to solve (HHL):', hhl_time, 'seconds\n')

    # Print Results
    print('Newton Raphson with HHL Solution:\n')
    for i in range(c3_hhl["iterations"]+1):
        print(f"Iteration {i}:")
        print("    Voltages:", c3_hhl["voltages"][i])
        print("    Phase Angles:", c3_hhl["phase_angles"][i])
    print('\n')

    # Save Plot
    plt.plot(c3_hhl["dpav_difference"])
    plt.savefig("c3_hhl_dpav_difference.png")

if __name__ == "__main__":
    main()
