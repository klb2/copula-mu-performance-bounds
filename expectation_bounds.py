"""Basic functions for bounding the expected value of a function of two random
variables

This module contains different functions for calculating bounds on the expected
value of a function of two random variables.


Copyright (C) 2020 Karl-Ludwig Besser

This program is used in the article:
Karl-Ludwig Besser and Eduard Jorswieck, "Copula-Based Bounds for Multi-User
Communications", submitted to IEEE Communications Letters.

License:
This program is licensed under the GPLv3 license. If you in any way use this
code for research that results in publications, please cite our original
article listed above.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.
See the GNU General Public License for more details.

Author: Karl-Ludwig Besser, Technische Universität Braunschweig
"""

import numpy as np
from scipy import integrate
from scipy.special import expi

def inv_cdf(u, lam):  # Rayleigh fading
    return -np.log(1-u)/lam

def sum_rate(x, y, snr_x, snr_y):
    return np.log2(1 + snr_x*x + snr_y*y)

def mac_rate(x, y, snr):
    s = 1/snr
    return np.log2(1 + x/(s+y))

def sinr(x, y, s=1):
    return x/(s+y)


def exp_times_expi(x):
    return np.exp(x)*expi(-x)


### SUM RATE ###

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

#################



### MAC RATE ###

def _single_upper_mac_rate(lam_x, lam_y, snr):
    def cost_func(u, snr=1, lam_x=1, lam_y=1):
        _u = 1.-u
        return mac_rate(inv_cdf(u, lam=lam_x), inv_cdf(_u, lam=lam_y), snr=snr)
    expected = integrate.quad(cost_func, 0, 1, args=(snr,lam_x, lam_y))
    return expected[0]
upper_mac_rate = np.vectorize(_single_upper_mac_rate)

def lower_mac_rate(lam_x, lam_y, snr):
    s = 1./snr
    _part1 = np.exp(lam_y*s)*expi(-lam_y*s)
    _part2 = np.exp((lam_y*lam_x*s)/(lam_x+lam_y))*expi(-lam_x*lam_y*s/(lam_x+lam_y))
    return (_part1 - _part2)/np.log(2)

def indep_mac_rate(lam_x, lam_y, snr):
    s = 1./snr
    if lam_x == lam_y:
        lam = lam_x
        expect_xi = 1 + np.exp(lam*s)*(lam*s-1)*expi(-lam*s) + np.log(s)
        expect_psi = -np.exp(lam*s)*(expi(-lam*s)-np.exp(-lam*s)*np.log(s))
    else:
        expect_xi = -(np.exp(lam_x*s)*lam_y*expi(-lam_x*s)-np.exp(lam_y*s)*lam_x*expi(-lam_y*s)+np.log(s)*(lam_x-lam_y))/(lam_y-lam_x)
        expect_psi = -np.exp(lam_y*s)*(expi(-lam_y*s)-np.exp(-lam_y*s)*np.log(s))
    return (expect_xi - expect_psi)/np.log(2)

##############



### SINR ###

def lower_sinr(s=1, lam_x=1, lam_y=1):
    return (lam_y + lam_y**2*s*np.exp(lam_y*s)*expi(-lam_y*s))/lam_x

def _single_upper_sinr(s=1, lam_x=1, lam_y=1):
    def cost_func(u, s, lam_x, lam_y):
        _u = 1.-u
        return sinr(inv_cdf(u, lam=lam_x), inv_cdf(_u, lam=lam_y), s=s)
    expected = integrate.quad(cost_func, 0, 1, args=(s, lam_x, lam_y))
    return expected[0]
upper_sinr = np.vectorize(_single_upper_sinr)

def indep_sinr(s=1, lam_x=1, lam_y=1):
    expect_x = 1/lam_x
    expect_z = lam_y*np.exp(lam_y*s)*(-expi(-lam_y*s))  # for real x>0: Gamma(0, x)=-Ei(-x)
    return expect_x*expect_z

###############


def main(lam_x=1, lam_y=1):
    from probability_bounds import export_results
    def calc_mac_rate(snr_db, bound='min', lam_x=1, lam_y=1):
        snr_lin = 10**(snr_db/10.)
        if bound.startswith('ind'):
            return indep_mac_rate(lam_x=lam_x, lam_y=lam_y, snr=snr_lin)
        elif bound.startswith('min') or bound.startswith('low'):
            return lower_mac_rate(lam_x=lam_x, lam_y=lam_y, snr=snr_lin)
        elif bound.startswith('max') or bound.startswith('up'):
            return upper_mac_rate(lam_x=lam_x, lam_y=lam_y, snr=snr_lin)
        else:
            return NotImplemented
    snr_db = np.arange(-5, 21)
    keys = ['min', 'max', 'ind']
    results = {}
    for _key in keys:
        results[_key] = calc_mac_rate(snr_db, _key, lam_x=lam_x, lam_y=lam_y)
    results["snr"] = snr_db
    export_results(results, "expectation-mac-rate-lx{}-ly{}.dat".format(lam_x, lam_y))
    return results

if __name__ == "__main__":
    main()
