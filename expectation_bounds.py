import numpy as np
from scipy import integrate
from scipy.special import expi

def inv_cdf(u, lam):
    return -np.log(1-u)/lam

def sum_rate(x, y, snr_x, snr_y):
    return np.log2(1 + snr_x*x + snr_y*y)

def mac_rate(x, y, snr):
    s = 1/snr
    return np.log2(1 + x/(s+y))

def exp_times_expi(x):
    return np.exp(x)*expi(-x)

def lower_sum_rate(lam_x, lam_y, snr_x, snr_y):
    _param = (lam_x*lam_y)/(snr_x*lam_y+snr_y*lam_x)
    return -exp_times_expi(_param)/np.log(2)

def _single_upper_sum_rate(lam_x, lam_y, snr_x, snr_y):
    def cost_func(u, lam_x=1, lam_y=1, snr_x=1, snr_y=1):
        _u = 1.-u
        return sum_rate(inv_cdf(u, lam=lam_x), inv_cdf(_u, lam=lam_y),
                        snr_x=snr_x, snr_y=snr_y)
    expected = integrate.quad(cost_func, 0, 1, args=(lam_x, lam_y, snr_x, snr_y))
    return expected[0]

upper_sum_rate = np.vectorize(_single_upper_sum_rate)

def indep_sum_rate(lam_x, lam_y, snr_x, snr_y):
    rate = (exp_times_expi(lam_x/snr_x)*(snr_x-lam_x)-snr_x)/snr_x
    idx_uneq = np.where(lam_x/snr_x != lam_y/snr_y)
    if len(idx_uneq[0]) > 0:
        rate[idx_uneq] = (snr_x*lam_y*exp_times_expi(lam_x/snr_x)-snr_y*lam_x*exp_times_expi(lam_y/snr_y))/(snr_x*lam_y - snr_y*lam_x)
    return -rate/np.log(2)

def upper_mac_rate(lam_x, lam_y, snr_x, snr_y):
    def cost_func(u, s=1, lam_x=1, lam_y=1):
            if bound.startswith('max'):
                _u = 1.-u
            return func(inv_cdf(u, lam=lam_x), inv_cdf(_u, lam=lam_y), s=s)
    expected = integrate.quad(cost_func, 0, 1, args=(s,lam_x, lam_y))
    return expected[0]
