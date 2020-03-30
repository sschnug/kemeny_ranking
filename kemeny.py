from __future__ import division
from timeit import default_timer as time
from cylp.cy import CyClpSimplex
from cylp.py.pivots import PositiveEdgePivot
import importlib
import itertools
import numpy as np
import scipy.sparse as sp
from utils import combs, perms
from reductions import extended_condorcet_simple


class KemenyRanking():
    def __init__(self, fp, verbose=True, condorcet_red=True):
        self.verbose = verbose
        self.condorcet_red = True
        self.parse_file(fp)

    def parse_file(self, fp):
        """ Reads and preprocesses input """
        # TODO add checks
        # TODO add specification
        if self.verbose:
            print('Parse input')

        with open(fp, encoding='utf8') as file:
            content = file.readlines()
            content = [x.strip() for x in content]                          # remove newlines
            content = [x.replace(':', '') for x in content]                 # remove ":"
            content = [np.array(x.split(), dtype=object) for x in content]  # split line into list
                                                                            # -> array

            raw_arr = np.array(content)
            self.voters_raw = raw_arr[:, 0]
            self.votes_raw = raw_arr[:, 1:]

            # Map to 0, N -> only votes!
            self.orig2id = {}
            self.id2orig = {}
            id_ = 0
            for i in np.unique(self.votes_raw):
                self.orig2id[i] = id_
                self.id2orig[id_] = i
                id_ += 1
            self.votes_arr = np.vectorize(self.orig2id.get)(self.votes_raw)

        if self.verbose:
            print('     ... finished')

            print('Problem statistics')
            print('  {} votes'.format(self.votes_arr.shape[0]))
            print('  {} candidates'.format(self.votes_arr.shape[1]))

    def build_Q(self):
        """ Creates incidence-matrix: form used in MIP-model """
        if self.verbose:
            print('Build incidence-matrix')

        N, n = self.votes_arr.shape                                              # N votes, n cands
        self.Q = np.zeros((n,n))
        for a,b in itertools.combinations(range(n), 2):
            a_pos = np.where(self.votes_arr == a)[1]
            b_pos = np.where(self.votes_arr == b)[1]
            plus = np.count_nonzero(a_pos < b_pos)
            minus = np.count_nonzero(a_pos > b_pos)
            self.Q[a,b] = plus
            self.Q[b,a] = minus

        if self.verbose:
            print('     ... finished')

    def solve_ilp(self):
        """ Solves problem exactly using MIP/ILP approach
            Used solver: CoinOR CBC
            Incidence-matrix Q holds complete information needed for opt-process
        """
        if self.verbose:
            print('Solve: build model')

        if self.condorcet_red:
            condorcet_red_mat = extended_condorcet_simple(self.votes_arr)

        n = self.Q.shape[0]
        x_n = n*n

        model = CyClpSimplex()                                           # MODEL
        x = model.addVariable('x', x_n, isInt=True)                      # VARS

        model.objective = self.Q.ravel()                                 # OBJ

        # x_ab = boolean (already int; need to constrain to [0,1])
        model += sp.eye(x_n) * x >= np.zeros(x_n)
        model += sp.eye(x_n) * x <= np.ones(x_n)

        idx = lambda i, j: np.ravel_multi_index((i, j), (n,n))

        # constraints for every pair
        start_time = time()
        n_pairwise_constr = n*(n-1)//2
        if self.verbose:
            print('  # pairwise constr: ', n_pairwise_constr)

        # Somewhat bloated just to get some vectorization / speed !
        combs_ = combs(range(n), 2)

        inds_a = np.ravel_multi_index(combs_.T, (n, n))
        inds_b = np.ravel_multi_index(combs_.T[::-1], (n, n))

        row_inds = np.tile(np.arange(n_pairwise_constr), 2)
        col_inds = np.hstack((inds_a, inds_b))

        pairwise_constraints = sp.coo_matrix((np.ones(n_pairwise_constr*2),
                                              (row_inds, col_inds)),
                                              shape=(n_pairwise_constr, n*n))
        end_time = time()
        if self.verbose:
            print("    Took {:.{prec}f} secs".format(end_time - start_time, prec=3))

        # and for every cycle of length 3
        start_time = time()
        n_triangle_constrs = n*(n-1)*(n-2)
        if self.verbose:
            print('  # triangle constr: ', n_triangle_constrs)

        # Somewhat bloated just to get some vectorization / speed !
        perms_ = perms(range(n), 3)

        inds_a = np.ravel_multi_index(perms_.T[(0,1), :], (n, n))
        inds_b = np.ravel_multi_index(perms_.T[(1,2), :], (n, n))
        inds_c = np.ravel_multi_index(perms_.T[(2,0), :], (n, n))

        row_inds = np.tile(np.arange(n_triangle_constrs), 3)
        col_inds = np.hstack((inds_a, inds_b, inds_c))

        triangle_constraints = sp.coo_matrix((np.ones(n_triangle_constrs*3),
                                              (row_inds, col_inds)),
                                              shape=(n_triangle_constrs, n*n))
        end_time = time()
        if self.verbose:
            print("    Took {:.{prec}f} secs".format(end_time - start_time, prec=3))


        model += pairwise_constraints * x == np.ones(n_pairwise_constr)
        model += triangle_constraints * x >= np.ones(n_triangle_constrs)

        if self.condorcet_red and condorcet_red_mat != None:
            I, J, V = sp.find(condorcet_red_mat)
            indices_pos = np.ravel_multi_index([J, I], (n,n))
            indices_neg = np.ravel_multi_index([I, J], (n,n))
            nnz = len(indices_pos)

            if self.verbose:
                print('  Extended Condorcet reductions: {} * 2 relations fixed'.format(nnz))

            lhs = sp.coo_matrix((np.ones(nnz*2),
                        (np.arange(nnz*2),
                         np.hstack((indices_pos, indices_neg)))),
                  shape=(nnz*2, n*n))
            rhs = np.hstack((np.ones(len(indices_pos)), np.zeros(len(indices_neg))))
            model += lhs * x == rhs

        cbcModel = model.getCbcModel()  # Clp -> Cbc model / LP -> MIP
        cbcModel.logLevel = self.verbose

        if self.verbose:
            print('Solve: run MIP\n')
        start_time = time()
        status = cbcModel.solve()           #-> "Call CbcMain. Solve the problem
                                            #   "using the same parameters used
                                            #   "by CbcSolver."
                                            # This deviates from cylp's docs which are sparse!
                                            # -> preprocessing will be used and is very important!
        end_time = time()
        if self.verbose:
            print("  CoinOR CBC used {:.{prec}f} secs".format(end_time - start_time, prec=3))

        x_sol = cbcModel.primalVariableSolution['x']
        self.obj_sol = cbcModel.objectiveValue
        x = np.array(x_sol).reshape((n, n)).round().astype(int)
        self.aggr_rank = np.argsort(x.sum(axis=0))[::-1]

    def postprocess(self):
        if self.verbose:
            print('Postprocessing')
        self.final_solution = np.vectorize(self.id2orig.get)(self.aggr_rank)
        if self.verbose:
            print('    ... finished')

    def print_sol(self):
        print('--------')
        print('SOLUTION')
        print('  objective: ', self.obj_sol)
        print('  aggregation: ')
        print(self.final_solution)
