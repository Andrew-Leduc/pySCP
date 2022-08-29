'''
b-matching implementation for emDyn, ported from matlab
'''
import numpy as np
import scipy.io       as sc_io
import scipy.sparse   as sc_sp
import scipy.optimize as sc_opt

from collections import namedtuple

import logging
import os, time

class graph_deg(object):
    def __init__(self, a, b):
        self.a = 1
        self.b = int( np.floor( self.a * a / b))


class graph_ext(object):
    def __init__(self, a, b, deg):
        self.rem = deg.a*a - b*deg.b
        self.ind = np.random.permutation(b)
        self.Inds = self.ind[:self.rem]


def load_mtl_synth_data():
    ''' load test-data + results from matlab to verify '''
    DATA_DIR       = '/home/dmaliout_work/Work_crd/Octave/octave_emdyn/Data'
    #TEMP_FILE_BMAT = 'dmm_temp_bmatching_data.mat'
    TEMP_FILE_BMAT = 'debug_bmat_data_p2.mat'

    fname_data = os.path.join(DATA_DIR,
                              TEMP_FILE_BMAT)
    mtl_data = sc_io.loadmat(fname_data)

    return mtl_data


class b_matching_solver(object):
    ''' class to run b-matching '''
    def __init__(self, num_row, num_col):
        ''' pre-compute some arrays to speed up the actual solver '''

        assert num_row >= num_col, "Transpose the matrix s.t. num_row >= num_col"

        self.num_row = num_row
        self.num_col = num_col

        a = max( num_row, num_col)
        b = min( num_row, num_col)

        self.a, self.b = a, b
        self.Nx = self.a * self.b

        #deg.a = 1;                       # Vertex Degree for Vertices a
        #deg.b = floor(  deg.a*a/b );     # Vertex Degree for Vertices b
        deg = graph_deg(a, b)

        #ext.rem = deg.a*a - b*deg.b;
        #ext.ind = randperm( b );
        #ext.Inds = ext.ind( 1:ext.rem );
        ext = graph_ext(a, b, deg)

        self.B = np.r_[deg.b*np.ones(b),\
                       deg.a*np.ones(a)]

        if  ext.rem >=1:
            raise NotImplementedError('rem > 0 not tested yet, debug first')
            B[ext.Inds] += 1

        # Populate the constraint matrix A:
        A_idxi_k = np.zeros(2*a*b)
        A_idxj_k = np.zeros(2*a*b)
        for i in range(b):
            idx_k = i*a + np.arange(a, dtype=int)
            A_idxi_k[idx_k] = i
            A_idxj_k[idx_k] = idx_k
        for i in range(a):
            idx_k = (a*b) + i*b + np.arange(b, dtype=int)
            rng_ab_step_a = np.arange(0, a*b, a, dtype=int)
            A_idxi_k[idx_k] = b+i
            A_idxj_k[idx_k] = i + rng_ab_step_a
        A_vals_k = np.ones(2*a*b)
        # create the sparse matrix from the start
        self.A = sc_sp.coo_matrix((A_vals_k, (A_idxi_k, A_idxj_k)), (a+b, a*b))

        l_bounds = np.zeros(self.Nx)
        u_bounds = np.ones (self.Nx)
        self.bounds   = np.array(list(zip(l_bounds, u_bounds)))


    def solve(self, mtx_ij):
        ''' run the LP solver for b-matching '''
        assert mtx_ij.shape == (self.num_row, self.num_col), "DEBUG, dim mismatch"
        # matlab and python flatten in rev order, hence .T
        obj_vec = mtx_ij.T.flatten()
        obj_vec = 10*obj_vec / np.max(np.abs(obj_vec))
        assert self.Nx == len(obj_vec), "DEBUG, dimension mismatch"

        res = sc_opt.linprog(-obj_vec, A_eq=self.A, b_eq=self.B, #A_ub=sp_A, b_ub=B,
                             bounds=self.bounds, method='highs')

        assert res.success, "Solver failed, status: %d.  DEBUG" % res.status

        # reshape into a matrix form (paying respect to matlab vs python mtx order)
        x_opt = res.x.copy()
        X_opt = np.reshape(x_opt, (self.b, self.a)).T

        #[inds.x, inds.y] = ind2sub( size(r),  find( x>1e-2 ))
        #[vals, ins] = sort( inds.x  )
        #inds.y = inds.y( ins )
        idx_j, idx_i = np.unravel_index( np.where(res.x)[0],\
                                         mtx_ij.T.shape )
        idx_srt = np.argsort(idx_i)
        idx_i, idx_j = idx_i[idx_srt], idx_j[idx_srt]
        #idx_opt = (idx_i, idx_j)
        IdxTuple = namedtuple('IdxTuple', ['idx_i', 'idx_j'])
        idx_opt = IdxTuple(idx_i, idx_j)

        return (x_opt, idx_opt, X_opt)

#@profile  # uncomment to profile with line_profiler / kernprof
def b_matching(mtx_ij, options={}):
    ''' port of emDyn b-matching '''

    num_row, num_col = mtx_ij.shape
    if num_row < num_col:
        mtx_ij = mtx_ij.T

    a = max( num_row, num_col)
    b = min( num_row, num_col)

    #deg.a = 1;                       # Vertex Degree for Vertices a
    #deg.b = floor(  deg.a*a/b );     # Vertex Degree for Vertices b
    deg = graph_deg(a, b)

    #ext.rem = deg.a*a - b*deg.b;
    #ext.ind = randperm( b );
    #ext.Inds = ext.ind( 1:ext.rem );
    ext = graph_ext(a, b, deg)

    B = np.r_[deg.b*np.ones(b),\
              deg.a*np.ones(a)]

    if  ext.rem >=1:
        #raise NotImplementedError('rem > 0 not tested yet, debug first')
        B[ext.Inds] += 1

    #A = np.zeros((a+b, a*b))
    #for i in range(b):
    #    A[i, (i-1)*a + (0:a)] = 1
    #for i in range(a):
    #    A[b+i,  (i-1) + (1 : a : (a*b))] = 1
    #logging.warn('Create A in a sparse-aware way from the start')
    #A = sparse(A)
    #sp_A = sc_sp.coo_matrix(A)

    # Populate the constraint matrix A:
    A_idxi_k = np.zeros(2*a*b)
    A_idxj_k = np.zeros(2*a*b)
    for i in range(b):
        idx_k = i*a + np.arange(a, dtype=int)
        A_idxi_k[idx_k] = i
        A_idxj_k[idx_k] = idx_k
    for i in range(a):
        idx_k = (a*b) + i*b + np.arange(b, dtype=int)
        rng_ab_step_a = np.arange(0, a*b, a, dtype=int)
        A_idxi_k[idx_k] = b+i
        A_idxj_k[idx_k] = i + rng_ab_step_a
    A_vals_k = np.ones(2*a*b)
    # create the sparse matrix from the start
    A = sc_sp.coo_matrix((A_vals_k, (A_idxi_k, A_idxj_k)), (a+b, a*b))

    ## Part 2:  invoke the LP solver (GLPK)
    use_sc_opt = True
    if use_sc_opt:
        Nx = a*b
        # matlab and python flatten in rev order, hence .T
        obj_vec = mtx_ij.T.flatten()
        obj_vec = 10*obj_vec / np.max(np.abs(obj_vec))
        assert Nx == len(obj_vec), "DEBUG"

        l_bounds = np.zeros(Nx)
        u_bounds = np.ones(Nx)
        bounds   = np.array(list(zip(l_bounds, u_bounds)))

        res = sc_opt.linprog(-obj_vec, A_eq=A, b_eq=B, #A_ub=sp_A, b_ub=B,
                             bounds=bounds, method='highs')

    else:  # port GLPK to use this
        sense = -1  #maximize / minimize (default)  1
        #ctype   = arrayfun(@(x)'S', zeros(a+b,1), 'UniformOutput', true )
        #vartype = cellfun( @(x)'C', cell(a*b,1),  'UniformOutput', true )
        ctype   = 'S' * len(a+b)
        vartype = 'C' * len(a*b)
        [x, fmin, status] = glpk( mtx_ij[:], A, B,\
                                np.zeros( mtx_ij[:].shape ),\
                                np.ones(  mtx_ij[:].shape ),\
                                ctype, vartype, sense
                            )

    assert res.success, "Solver failed, status: %d.  DEBUG" % res.status

    # reshape into a matrix form (paying respect to matlab vs python mtx order)
    x_opt = res.x.copy()
    X_opt = np.reshape(x_opt, (b, a)).T

    #[inds.x, inds.y] = ind2sub( size(r),  find( x>1e-2 ))
    #[vals, ins] = sort( inds.x  )
    #inds.y = inds.y( ins )
    idx_j, idx_i = np.unravel_index( np.where(res.x)[0], mtx_ij.T.shape )
    idx_srt = np.argsort(idx_i)
    idx_i, idx_j = idx_i[idx_srt], idx_j[idx_srt]
    #idx_opt = (idx_i, idx_j)
    IdxTuple = namedtuple('IdxTuple', ['idx_i', 'idx_j'])
    idx_opt = IdxTuple(idx_i, idx_j)

    return (x_opt, idx_opt, X_opt)


def run_b_matching():
    ''' '''
    num_row = 10
    num_col = 20

    mtx_ij = np.random.randn(num_row, num_col)
    out = b_matching(mtx_ij, options={})


def test_mtl_b_matching():
    ''' verify that we can match matlab synth results '''
    import pylab as plt, seaborn as sns

    mtl_data = load_mtl_synth_data()
    mtx_ij = mtl_data['mtx_ij']

    x, idx, X = b_matching(mtx_ij, options={})

    sns.heatmap(X); plt.show()
    print('Done.')

if __name__ == '__main__':
    #run_b_matching()
    test_mtl_b_matching()