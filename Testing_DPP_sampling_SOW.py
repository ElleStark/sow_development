# Applies fixed-size Determinantal Point Processing (k-DPP) sampling - method proposed in Calandriello et al., 2020:
# 10.48550/arxiv.2006.16947
# Purpose: obtain diverse subset of desired size from full-factorial SOW set
import sklearn.metrics.pairwise
from dppy.finite_dpps import FiniteDPP
import numpy as np
from numpy.random import randn
from sklearn.utils.extmath import row_norms
import matplotlib.pyplot as plt

# Function for computing likelihood kernel, in our case the similarity matrix.
# Code below based on sklearn implementation of an rbf kernel, which calculates similarity from Euclidean distances.
# Likelihood Evaluation Equation:
def eval_L_rbf(X, Y=None, gamma=None):
    #X, Y = sklearn.metrics.pairwise.check_pairwise_arrays(X, Y)
    X = np.atleast_2d(X)
    # Need to check what gamma parameter is for
    if gamma is None:
        gamma = 1.0 / X.shape[1]

    if Y is None:
        XX = row_norms(X, squared=True)[:, np.newaxis]
        K = -2 * (X.dot(X.T))
        K += XX
        K += XX
        K *= -gamma
        np.exp(K, K)  # exponentiate K in-place
        return K
    else:
        XX = row_norms(X, squared=True)[:, np.newaxis]
        YY = row_norms(Y, squared=True)[np.newaxis, :]
        K = -2 * (X.dot(Y.T))
        K += XX
        K += YY
        K *= -gamma
        np.exp(K, K)  # exponentiate K in-place
        return K


k = 300
test_data = randn(10000, 2)
np.linalg.matrix_rank(test_data)
#all_sow = np.loadtxt('data/full_factorial_sow_norm.csv', delimiter=',')
#print(all_sow[0:5,:])

#### Testing intermediate steps/ number of expected samples
def evaluate_L_diagonal(eval_L, X):
    """Helper function to evaluate a likelihood function on a set of points (i.e. compute the diagonal of the L matrix)"""
    diag_eval = getattr(eval_L, "diag", None)
    if callable(diag_eval):
        return diag_eval(X)
    else:
        # inspired by sklearn.gaussian_process.kernels.PairwiseKernel
        return np.apply_along_axis(eval_L, 1, X).ravel()

X_data = test_data
eval_L = 
diag_L = evaluate_L_diagonal(eval_L, X_data)
trace_L = diag_L.sum()

dict_bless = bless(X_data, eval_L, 1.0, rls_oversample_bless, rng,
                   nb_iter_bless=nb_iter_bless, verbose=verbose)

# Phase 1: use estimate RLS to sample the dict_dppvfx dictionary, i.e. the one used to construct A
# here theory says that to have high acceptance probability we need the oversampling factor to be ~deff^2
# but even with constant oversampling factor we seem to accept fast

D_A = reduce_lambda(X_data, eval_L, dict_bless, dict_bless.lam, rng, rls_oversample_parameter=rls_oversample_dppvfx)

# Phase 2: pre-compute L_hat, B_bar, l_i, det(I + L_hat), etc.
U_DD, S_DD, _ = np.linalg.svd(eval_L(D_A.X, D_A.X))
U_DD, S_root_inv_DD = stable_invert_root(U_DD, S_DD)
m = U_DD.shape[1]

E = S_root_inv_DD * U_DD.T

# The _T indicates that B_bar_T is the transpose of B_bar,
# we keep it that way for efficiency reasons
X_data = test_data
B_bar_T = E.dot(eval_L_rbf(D_A.X, X_data))
diag_L_hat = np.square(B_bar_T).sum(axis=0)
trace_L_hat = diag_L_hat.sum()

# While we have L_hat = B_bar_T.T * B_bar_T, we do not want to compute explicitly the (n x n) matrix
# instead we reason in terms of B_bar_T * B_bar_T.T which is a (m x m) matrix. We call this matrix A_mm.
# I_A_mm indicates I + A_mm (i.e. A_mm with identity added)
I_A_mm = B_bar_T.dot(B_bar_T.T)
I_A_mm[np.diag_indices(m)] += 1.0

# we now need to compute the l_i estimates using L_hat, it is more efficient to do it in terms of
# B_bar_T and I_A_mm
# in particular, we will use the diag(L - L_hat + L_hat(L_hat + I)^-1) estimator
# but we must first tune L to obtain a desired s
# we can use the fact the the non-zero eigenvalues of I + L_hat and I_A_mm are equal
eigvals, eigvec = np.linalg.eigh(I_A_mm)

if np.any(eigvals <= 1.0):
    raise ValueError('Some eigenvalues of L_hat are negative, this should never happen. '
                     'Minimum eig: {}'.format(np.min(eigvals - 1.0)))

natural_expected_size = trace_L - trace_L_hat + np.sum((eigvals - 1.0) / eigvals)

##### End test data

DPP = FiniteDPP('likelihood',
                **{'L_eval_X_data': (eval_L_rbf, test_data.T)})
DPP.flush_samples()
sample = DPP.sample_exact_k_dpp(size=k, mode='alpha', random_state=26, )
samples = test_data.T[sample_index, :]
samples = samples.T

plt.scatter(x=test_data[:, 0], y=test_data[:, 1])
plt.scatter(x=samples[:, 0], y=samples[:, 1], color='red')
plt.show()

sample = DPP.sample_exact_k_dpp(size=k, mode='alpha')
sow_subset = all_sow[sample, :]
#sow_subset = DPP.list_of_samples
print(sow_subset[0:5,:])
savetxt('initial_sow_set_kdpp.csv', sow_subset, delimiter=',')