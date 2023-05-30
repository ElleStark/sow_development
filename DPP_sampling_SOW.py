# Applies fixed-size Determinantal Point Processing (k-DPP) sampling - method proposed in Calandriello et al., 2020:
# 10.48550/arxiv.2006.16947
# Purpose: obtain diverse subset of desired size from full-factorial SOW set
import sklearn.metrics.pairwise
from dppy.finite_dpps import FiniteDPP
import numpy as np
from numpy.random import randn
from sklearn.utils.extmath import row_norms
import matplotlib.pyplot as plt
from sklearn.gaussian_process.kernels import RBF, PairwiseKernel, Kernel

# Class FastMixinKernel taken from https://github.com/LCSL/dpp-vfx/blob/master/exp_dppy_mnist_first_sample.py
class FastMixinKernel(Kernel):
    def __init__(self, gp_kernel, pairwise_kernel):
        self.gp_kernel = gp_kernel
        self.pairwise_kernel = pairwise_kernel

    def __call__(self, X, Y=None, **kwargs):
        return self.pairwise_kernel(X, Y, **kwargs)

    def diag(self, X):
        return self.gp_kernel.diag(X)

    def is_stationary(self):
        return self.gp_kernel.is_stationary()

# Function for computing likelihood kernel, in our case the similarity matrix.
# Code below based on sklearn implementation of an rbf kernel, which calculates similarity from Euclidean distances.
# Likelihood Evaluation Equation:
def eval_L_rbf(X, Y=None, gamma=None):
    # rbf kernel definition based on implementation in sklearn.metrics.pairwise.rbf_kernel

    X = np.atleast_2d(X)

    if gamma is None:
        gamma = 1.0 / (X.shape[1])

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


k = 20
test_data = randn(20000, 2)

# Implementation based on https://github.com/LCSL/dpp-vfx/blob/master/exp_dppy_mnist_first_sample.py
sigma = np.sqrt(3*test_data.shape[1])
dot_func = FastMixinKernel(
    RBF(sigma),
    PairwiseKernel(gamma=1/np.square(sigma), metric='rbf', pairwise_kernels_kwargs={'n_jobs': -2})
)



#all_sow = np.loadtxt('data/full_factorial_sow_norm.csv', delimiter=',')
#print(all_sow[0:5,:])

DPP = FiniteDPP('likelihood',
                **{'L_eval_X_data': (dot_func, test_data)})
DPP.flush_samples()
sample = DPP.sample_exact_k_dpp(size=k, mode='alpha', **{'rls_oversample_alphadpp': 4, 'rls_oversample_bless':4})
samples = test_data[sample, :]

plt.scatter(x=test_data[:, 0], y=test_data[:, 1])
plt.scatter(x=samples[:, 0], y=samples[:, 1], color='red')
plt.show()

sample = DPP.sample_exact_k_dpp(size=k, mode='alpha')
sow_subset = all_sow[sample, :]
#sow_subset = DPP.list_of_samples
print(sow_subset[0:5,:])
savetxt('initial_sow_set_kdpp.csv', sow_subset, delimiter=',')