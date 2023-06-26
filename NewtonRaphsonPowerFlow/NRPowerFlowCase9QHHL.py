from PowerFlow import PowerData, NRPF

def main():
    # Load simplified case9Q data (found in case9Q_simplified.m)
    case9Q = PowerData()
    case9Q.addBus(1, 0, 0, 1, 0, 1)
    case9Q.addBus(2, 0, 0, 1, 0, 2)
    case9Q.addBus(3, 0, 0, 1, 0, 2)
    case9Q.addBus(4, 0, 0, 1, 0, 3)
    case9Q.addBus(5, 0.9, 0.3, 1, 0, 3)
    case9Q.addBus(6, 0, 0, 1, 0, 3)
    case9Q.addBus(7, 1, 0.35, 1, 0, 3)
    case9Q.addBus(8, 0, 0, 1, 0, 3)
    case9Q.addBus(9, 1.25, 0.5, 1, 0, 3)
    case9Q.addGen(1, 0, 0, 1, -3, 3)
    case9Q.addGen(2, 1.63, 0, 1, -3, 3)
    case9Q.addGen(3, 0.85, 0, 1, -3, 3)
    case9Q.addBranch(1, 4, complex(0, 0.0576))
    case9Q.addBranch(4, 5, complex(0.017, 0.092))
    case9Q.addBranch(5, 6, complex(0.039, 0.17))
    case9Q.addBranch(3, 6, complex(0, 0.0586))
    case9Q.addBranch(6, 7, complex(0.0119, 0.1008))
    case9Q.addBranch(7, 8, complex(0.0085, 0.072))
    case9Q.addBranch(8, 2, complex(0, 0.0625))
    case9Q.addBranch(8, 9, complex(0.032, 0.161))
    case9Q.addBranch(9, 4, complex(0.01, 0.085))
    print("TEST")
    # Initialize Classical solver
    pf = NRPF()
    # Initialize HHL solver
    pf_hhl = NRPF('hhl')

    # Solve Classically
    c9Q_pf = pf.NR(case9Q, 1e-8)
    # Solve with HHL
    c9Q_hhl = pf_hhl.NR(case9Q, 1e-8)

if __name__ == "__main__":
    main()