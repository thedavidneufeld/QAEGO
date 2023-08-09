from PowerFlow import PowerData, NRPF
import matplotlib.pyplot as plt
from pathlib import Path

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

    # Solve Classically
    c9Q_pf = pf.NR(case9Q, 1e-8)

    # Print Results
    print('Newton Raphson Solution:\n')
    for i in range(c9Q_pf["iterations"]+1):
        print(f"Iteration {i}:")
        print("    Voltages:", c9Q_pf["voltages"][i])
        print("    Phase Angles:", c9Q_pf["phase_angles"][i])
    print('\n')

    # Save Plot
    plt.plot(c9Q_pf["dpav_difference"])
    plt.savefig("c9Q_pf_dpav_difference.png")
    plt.clf()

    # Solve with HHL
    c9Q_hhl = pf_hhl.NR(case9Q, 1e-8)

    # Print Results
    print('Newton Raphson with HHL Solution:\n')
    for i in range(c9Q_hhl["iterations"]+1):
        print(f"Iteration {i}:")
        print("    Voltages:", c9Q_hhl["voltages"][i])
        print("    Phase Angles:", c9Q_hhl["phase_angles"][i])
    print('\n')

    # Save Plot
    plt.plot(c9Q_hhl["dpav_difference"])
    plt.savefig("c9Q_hhl_dpav_difference.png")

if __name__ == "__main__":
    main()