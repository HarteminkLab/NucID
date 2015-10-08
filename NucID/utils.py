# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
import pandas as pd
import statsmodels.api as sm
from likelihood_functions import bayes_factor_marginal_over_a

# <codecell>

def moving_window_score(dnase, ignore_all_zero = False, score_function = bayes_factor_marginal_over_a, 
                        half_window_size = 73, **kwarg):
    '''given a DataFrame of DNase data, calculate the moving window score using the score function
    
    dnase -- a pandas DataFrame like:

            chr  coordinate  fwd_count  rev_count
        0    1           9          3          0
        1    1          10          4          0
        2    1          11         22          0
        3    1          12         10          0
        4    1          13         27          0

    score_function -- the function used to calculate the score. Should follow the following:
    
                    1. The function should be able to take a 1x(2J+1) array and return the score of it.
                    2. The function should always treat the array as a window from the forward strand. Reverse strand data
                    will be reversed before passed to the function. 
                    3. The first argument of the function should be the data array
                        
                    By default, the score_function is bayes_factor_marginal_over_a
                    
    half_window_size -- (2*half_window_size + 1) is the size of the window
    
    ignore_all_zero -- if the data in a window is all 0, do I ignore this window or not
    '''
    
    fwd = 'fwd_count'
    rev= 'rev_count'
    
    dnase_len = dnase.shape[0]
    window_size = 2*half_window_size + 1
    iterations = range(dnase_len - window_size + 1)
    
    score = np.zeros(len(iterations))
    
    for x in iterations:
        try:
            moving_data = np.vstack([dnase[fwd][x:(x + window_size)].values.reshape(1,-1), 
                                    dnase[rev][x:(x + window_size)].values[::-1].reshape(1,-1)])
            if ignore_all_zero and np.sum(moving_data == 0) == moving_data.size:
                score[x] = np.nan
            else:
                score[x] = score_function( moving_data , **kwarg)
                
        except:
            score[x] = np.nan
    
    
    return pd.DataFrame({
                'chr':dnase.chromo.iloc[0],
                'coordinate': dnase.coordinate.iloc[(half_window_size):(dnase.shape[0] - half_window_size)],
                'score':score
            })

# <codecell>

def lowess_smooth(series, fracs = None, verbose = False):
    '''smooth a series of values using LOWESS
    
    series -- a 1 dimensional array
    fracs -- the fraction of series used to calculate the smoothing curve. By default, a 20bp window is used for smoothing
    '''

    x = np.arange(len(series))
    
    if fracs is None:
        fracs = 20.0 / len(series)
        
    smoothed = sm.nonparametric.lowess(series, x, frac = fracs)[:, 1]
        
    return smoothed

