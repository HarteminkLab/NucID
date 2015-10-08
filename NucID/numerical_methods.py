# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
import scipy.optimize as sopt
import scipy.integrate as sint

# <headingcell level=1>

# Normal log p.d.f

# <codecell>

def matrix_normal_logpdf(mat, mu, sd):
    '''calculate the log likelihood on a matrix of data, each column shares a mean and sd
    mat -- an N x J numpy matrix
    mu -- an 1xJ numpy matrix for mean
    sd -- an 1xJ numpy matrix for standard deviation
    '''
    sd2 = sd**2
    
    N = mat.shape[0]; J = mat.shape[1]
    A = -0.5 * N * J * np.log(2*np.pi)
    
    _B = np.log(sd)
    B = -N*_B.sum()
    
    C = -np.sum( (mat**2) / (2.0 * sd2) )
    D = np.sum(1.0 * mat * mu / sd2) # in case they are all integers
    
    _E = mu**2 / (2.0 * sd2)
    E = -N * np.sum(_E)
    
    return A + B + C + D + E

# <codecell>

def normal_logpdf(x, mu, sd):
    '''returns the log of a normal density
    mu -- mean of the normal distribution 
    sd -- standard deviation of the normal distribution
    '''
    
    return -0.5 * np.log(2*np.pi) - np.log(sd) - 0.5 / sd**2 * (x - mu)**2

# <headingcell level=1>

# Numerical integration

# <codecell>

def calculate_marginal_likelihood_quadature_integration_scalar(integrand_func, bounds, *args):
    '''
    Calculate the numerical integration of a likelihood function. Returns the log of the marginal likelihood.
    
    integrand_func -- the integrand to integrate, the first argument should be the variable to integrate. The last three 
    argument should be "adjustment", "log: True or False", "sign"
    
    bounds -- the upper and lower bound for integration
    '''
    
    # Find the max value of the integrand_func over the integration bounds
    opt_res = sopt.minimize_scalar(
                    integrand_func, bounds = bounds, 
                    args = tuple(list(args) + [0, True, -1]), method = 'bounded')
    
    
    adjustment = -1 * opt_res.fun
    
    int_res = sint.quad(integrand_func, bounds[0], bounds[1], 
                        args = tuple(list(args) + [adjustment, False, 1]), 
                        epsabs = 0, epsrel = 1e-5)
    
    return np.log(int_res[0]) + adjustment, int_res

