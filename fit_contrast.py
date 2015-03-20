import numpy as np
from scipy.misc import comb

def binomial_dist(k, n, p):
    return comb(n, k)*(p**k)*(1-p)**(n-k)

def log_likelihood(params, p_measured, theta, n = 100):
    C, phi, b = params
    p = (C/2)*(1 + np.sin(theta + phi)) + b
    log_likelihood = 0
    p_measured = np.array(p_measured)
    lkhood = np.sum(np.log(binomial_dist( (n*p_measured).round(), n, p)))
    if np.isnan(lkhood): return -np.inf
    else: return lkhood

def log_prior(params):
    # restrict the amplitude of the oscillation to 1
    C, phi, b = params
    if C < 0 or C > 1:
        return -np.inf
    elif b < 0 or b > 1:
        return -np.inf
    elif C + b > 1:
        return -np.inf
    else:
        return 0

def log_posterior(params, p_measured, theta):
    return log_prior(params) + log_likelihood(params, p_measured, theta)

def posterior(params, p_measured, theta):
    return np.exp(log_posterior(params, p_measured, theta))

def marginalized_posterior_coarse(p_measured, theta):
    '''
    Coarsely marginalize the posterior distribution over each variable. Find the max probability
    for each one, and find a range around this where the probability is greater than some cutoff
    '''
    
    bin_size = 0.1
    C_grid = np.arange(0., 1., bin_size)
    phi_grid = np.arange(0, 2*np.pi, bin_size)
    b_grid = np.arange(0., 1., bin_size)
    
    posterior_grid = np.zeros([len(C_grid), len(phi_grid), len(b_grid)])
    
    for i, c in enumerate(C_grid):
        for j, p in enumerate(phi_grid):
            for k, b in enumerate(b_grid):
                posterior_grid[i][j][k] = posterior([c, p, b], p_measured, theta)
    
    cutoff = 1e-3 # probability cutoff below which to ignore to save computation time
    # list of marginal distributions of [C, phi, b].
    # then take only the part of the marginal distribution where the probability is bigger than the cutoff
    marginals = [np.sum(posterior_grid,axis=(1,2)), np.sum(posterior_grid,axis=(0,2)), np.sum(posterior_grid,axis=(0,1)) ]
    clips = [np.where( marginal >= cutoff*np.max(marginal))[0] for marginal in marginals ]
    #plot(marginals[1])
    # This looks complicated. What it does:
    # evaluate each grid region (C, phi, b) at the endpoints of the clips, and add half a bin on each side for safety
    # then return the values
    clipped_regions = [ grid[[clp[0], clp[-1]]] + [-bin_size/2, bin_size/2] 
                       for clp, grid in zip(clips, [C_grid, phi_grid, b_grid]) ]
    return clipped_regions

def marginalized_posterior(p_measured, theta):
    '''
    Marginalize the posterior distribution over phi and b
    by computing the posterior distribution and summing over
    the phi and b axes. Normalize so that the sum of each posterior
    distribution is equal to 1.
    '''
    # region of interest for doing the bayesian thing
    roi = marginalized_posterior_coarse(p_measured, theta)
    C_grid, phi_grid, b_grid = [np.linspace(mn, mx, 30) for mn, mx in roi]    
    posterior_grid = np.zeros([len(C_grid), len(phi_grid), len(b_grid)])
    
    for i, c in enumerate(C_grid):
        for j, p in enumerate(phi_grid):
            for k, b in enumerate(b_grid):
                posterior_grid[i][j][k] = posterior([c, p, b], p_measured, theta)
    marginalized_posterior = [np.sum(posterior_grid,axis=(1,2)), np.sum(posterior_grid,axis=(0,2)), np.sum(posterior_grid,axis=(0,1))]
    marginalized_posterior = [post/sum(post) for post in marginalized_posterior] # normalize
    return [C_grid, phi_grid, b_grid], marginalized_posterior

def generate_errorbars(p_measured, theta):
    '''
    Right now this will only generate the error bars for the contrast.
    '''
    grids, post = marginalized_posterior(p_measured, theta)
    c_grid, phi_grid, b_grid = grids
    c_dist, phi_dist, b_dist = post
    
    cumdist = np.cumsum(c_dist)
    credible_region = np.where((cumdist>.05)&(cumdist<0.95))[0]
    min_index, max_index = credible_region[0], credible_region[-1]
    return c_grid[np.where(c_dist == np.max(c_dist))][0], (c_grid[min_index], c_grid[max_index])
    
