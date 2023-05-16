# QAEGO (Quantum Algorithms for Electrical Grid Optimization) Research
# University of Lethbridge
# David Neufeld 
# May 2023

import numpy as np
from numpy.linalg import norm
from math import ceil, log
# linear_solvers can be installed from
# https://github.com/anedumla/quantum_linear_solvers
from linear_solvers import NumPyLinearSolver, HHL
from qiskit import Aer
from qiskit.quantum_info import Statevector

class hhl_helper:

    def _make_nxn(self, matrix, vector):
        if matrix.shape[0] != matrix.shape[1]:
            # set n to the highest value between the rows and columns of matrix
            # by doing this, we can ensure that the matrix will be nxn and the vector will be nx1
            n = max(matrix.shape[0], matrix.shape[1])
            matrix = np.pad(matrix, ((0, n-matrix.shape[0]), (0, n-matrix.shape[1])), 'constant', constant_values = (0))
            vector.resize(1, n, refcheck=False)
            # create a matrix (lambda) that contains a small value to add to the matrix and vector
            # this will ensure that the matrix is invertible and still close in value to the original
            la = np.empty(n)
            la.fill(10**-6)
            la_mat = np.diag(la)
            matrix = matrix + la_mat
            vector = vector + la
        return matrix, vector
    
    def _make_2nx2n(self, matrix, vector):
        # this function assumes that the matrix is square
        # if matrix is not 2^n x 2^n, then set n to the next highest power of 2
        if not log(matrix.shape[0], 2).is_integer():
            n = 2**(ceil(log(matrix.shape[0], 2)))
            matrix = np.pad(matrix, ((0, n-matrix.shape[0]), (0, n-matrix.shape[1])), 'constant', constant_values = (0))
            vector.resize(1, n, refcheck=False)
            for i in range(n):
                if matrix[i][i] == 0:
                    matrix[i][i] = 1
        return matrix, vector
    
    def _make_hermitian(self, matrix, vector):
        # this function assumes that the matrix is square
        # if the matrix is not hermitian, make it hermitian and alter the vector accordingly
        matrix_H = np.matrix(matrix).H
        if not np.array_equal(matrix, matrix_H):
            zeros = np.zeros((matrix.shape[0], matrix.shape[0]))
            matrix = np.block([[zeros, matrix_H],
                        [matrix, zeros]])
            vector_conj = np.conj(vector)
            vector = np.append(vector_conj, vector)
        return matrix, vector

    def _hhl_compatible(self, matrix, vector):
        mat = matrix.copy()
        vec = vector.copy()
        mat, vec = self._make_nxn(mat,vec)
        mat, vec = self._make_2nx2n(mat, vec)
        mat, vec = self._make_hermitian(mat, vec)
        return mat, vec
    
    def _new_HHL(self, ep):
        backend = Aer.get_backend('aer_simulator')
        self.hhl = HHL(ep, quantum_instance=backend)

    def run_HHL(self, matrix, vector, ep):
        self._new_HHL(ep)
        mat, vec = self._hhl_compatible(matrix, vector)
        solution = self.hhl.solve(mat, vec)
        sv = Statevector(solution.state)
        data_start = int(sv.data.size / 2)
        data_end = data_start + vec.size
        sv_data = sv.data[data_start:data_end].real
        solution_norm = solution.euclidean_norm
        data = solution_norm * sv_data / norm(sv_data)
        vec_norm = norm(vec)
        final_data = (vec_norm * data)[0:vector.size]
        return final_data