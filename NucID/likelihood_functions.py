
# coding: utf-8

# In[1]:

import numpy as np
import warnings
import numerical_methods as num
import parameters as para


# ##Parameterized mean curve as $e^a + expb \times (j - c)^2 + \lambda_j$

# where $expb = e^b$. $c$ is always set to 0 here.

# In[ ]:

def dnase_asinh_normalLogLik_quadmean_expa_with_oscillation(
        a, dnase_asinh, oscillation, expb = 0, c = 0, k = 1):
    '''parameterized mean curve as exp(a) + expb*(j - c)^2 + lambda_j
    dnase_asinh: a n x (2J+1) matrix
    '''
    
    J = (dnase_asinh.shape[1] - 1)/2
    j2 = ((np.arange(-J, J + 1) - c) ** 2).reshape(1, -1)
    mu =  np.exp(a) + expb * j2 + oscillation
    
    if(np.sum(mu < 0) > 0):
        print 'a mu is < 0'
        
    var = k * mu
    sd = np.sqrt(var)

    return num.matrix_normal_logpdf(dnase_asinh, mu, sd)


# ##Use the same likelihood function as above, but with $a$ integrated out

# In[ ]:

def dnase_asinh_normalLogLik_quadmean_expa_marginala_with_oscillation_integrand(
        a, dnase_asinh, oscillation, expb = 0, c = 0, k = 1, mu_a = 0, sd_a = 1,
        adjustment = 0, log = False, sign = 1):
    '''dnase_asinh: a n x (2J+1) matrix
    '''

    log_integrand1 = dnase_asinh_normalLogLik_quadmean_expa_with_oscillation(a, dnase_asinh, oscillation, 
                                                                             expb = expb, c = c, k = k)
    log_integrand2 = num.normal_logpdf(a, mu_a, sd_a)
    log_integrand = log_integrand1 + log_integrand2 
    
    
    if log:
        return sign * (log_integrand - adjustment)
    else:
        return sign * np.exp(log_integrand - adjustment)


    
def dnase_asinh_normalLogLik_quadmean_expa_marginal_overa_with_oscillation(
        dnase_asinh, oscillation,expb = 0, c = 0, k = 1, mu_a = 0, sd_a = 1):
    '''normal likelihood function for a asinh transformed DNase data window. 
    the mean curve is modelled as exp(a) + expb*(j - c)^2 + lambda_j, 
    where lambda_j is the detrended oscillation series. The parameter a is integrated over a normal
    prior with mean mu_a and standard deviation sd_a.
    
    dnase_asinh -- a nx(2J+1) matrix
    '''
    
    # integration range for a
    bounds = [mu_a - 5 * sd_a, mu_a + 5 * sd_a]
    
    if bounds[0] < np.log(np.max(-oscillation)):
        # make sure the mean is bounded larger than 0 because the variance is proportional to the mean
        bounds[0] = np.log(np.max(-oscillation)) 
        
    int_res = num.calculate_marginal_likelihood_quadature_integration_scalar(
                    dnase_asinh_normalLogLik_quadmean_expa_marginala_with_oscillation_integrand, bounds, 
                    dnase_asinh, oscillation, expb, c, k, mu_a, sd_a) 
    
    if int_res[1][0] / int_res[1][1] < 1e4:
        warnings.warn('integration relative error larger than 1e-4')
    
    return int_res[0]


## Use the same likelihood function as above, but with both $a$ and $b$ integrated out

# This requires the [cubature](https://github.com/saullocastro/cubature) package.

# In[ ]:

# integrand for the cubature package
def dnase_asinh_normalLogLik_quadmean_expab_marginalab_integrand(theta, dnase_asinh, adjustment = 0, c = 0, k = 1, 
                                                               mu_a = 0, sd_a = 1, mu_b = 0, sd_b = 1,
                                                               log = False):
    '''
    dnase_asinh: n x (2J+1) matrix
    theta = [a,b]
    '''

    aa = theta[0]; bb = theta[1]
    log_integrand1 = dnase_asinh_normalLogLik_quadmean_expab(theta, dnase_asinh, c = c, k = k)
    log_integrand2 = num.normal_logpdf(aa, mu_a, sd_a)
    log_integrand3 = num.normal_logpdf(bb, mu_b, sd_b)
    log_integrand = log_integrand1 + log_integrand2 + log_integrand3 - adjustment
    
    if log:
        return np.array([log_integrand])
    else:
        return np.array([np.exp(log_integrand)])


def dnase_asinh_normalLogLik_quadmean_expab_marginal_overab(dnase_asinh, c = 0, k = 1, 
                                                            mu_a = 0, sd_a = 1, mu_b = 0, sd_b = 1
                                                            ):
    
    import cubature as cuba
        
    # integration rectangle for a and b
    xmin = np.array([mu_a - 5*sd_a, mu_b - 5 * sd_b], float); xmax = np.array([mu_a + 5*sd_a, mu_b + 5*sd_b], float)
    
    int_res = cuba.cubature(2, dnase_asinh_normalLogLik_quadmean_expab_marginalab_integrand, xmin = xmin, xmax = xmax, 
                  adaptive = 'h', abserr = 0.0, relerr = 1e-5, args = tuple([dnase_asinh, 0, c, k, mu_a, sd_a, mu_b, sd_b]))
    
    
    if int_res[0][0] / int_res[1][0] < 1e4:
        warnings.warn('integration relative error larger than 1e-4')
    
    return np.log(int_res[0][0])


# ##Bayes factor

# A function to facilitate the calculation of Bayes factor, using the following:
#   * Default parameters from Crawford data in ***```parameters.py```***.
#   * The likelihood is ***```dnase_asinh_normalLogLik_quadmean_expa_marginal_overa_with_oscillation```*** (defined above). Mean curve is parameterized as $e^a + expb \times j^2 + \lambda_j$, where $a$ is integrated out, $expb$ is fixed at MLE and $\lambda_j$ is the detrended oscillation pattern.
#   

# In[4]:

def bayes_factor_marginal_over_a(dnase_asinh):
    
    likelihood = dnase_asinh_normalLogLik_quadmean_expa_marginal_overa_with_oscillation    
        
    nuc_likelihood = likelihood(
        dnase_asinh, oscillation = para.crawford_detrend_nuc_oscillation, expb = para.crawford_nuc_expb, 
        c = 0, k = para.crawford_nuc_k, mu_a = para.crawford_nuc_a_mu, sd_a = para.crawford_nuc_a_sd
        )
    
    background_likelihood = likelihood(
        dnase_asinh, oscillation = para.crawford_detrend_background_oscillation, expb = para.crawford_background_expb, 
        c = 0, k = para.crawford_background_k, mu_a = para.crawford_background_a_mu, sd_a = para.crawford_background_a_sd
        )
    
    return nuc_likelihood - background_likelihood

