"""Computations of the outage probability for dependent Rayleigh fading
channels with a dependency model from literature.

This module contains different functions to calculate the outage probability of
Rayleigh fading channels with a specific dependency model from literature.
Monte Carlo simulations are performed as verification.


Copyright (C) 2020 Karl-Ludwig Besser

This program is used in the article:
Karl-Ludwig Besser and Eduard Jorswieck, "Copula-Based Bounds for Multi-User
Communications", submitted to IEEE Communications Letters

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
from scipy import stats
from scipy import special
from scipy import integrate
import matplotlib.pyplot as plt

from probability_bounds import export_results

def pdf_single_var(s, rho, t):
    return 1./(1-rho)*np.exp(-(s+rho*t)/(1-rho))*special.i0(2*np.sqrt(s*rho*t)/(1-rho))

def pdf_multi_var(s, rho, t, n_channels=2):
    _a_l = np.sqrt(n_channels*rho*t)
    _part1 = s**((n_channels-1)/2)/((1-rho)*_a_l**(n_channels-1))
    _part2 = np.exp(-(s+_a_l**2)/(1-rho))
    _part3 = special.iv(n_channels-1, 2*_a_l*np.sqrt(s)/(1-rho))
    return _part1*_part2*_part3

def cdf_multi_var(s, rho, t, n_channels=2):
    _a_l = np.sqrt(n_channels*rho*t)
    df = 2*n_channels
    _var = (1-rho)/2
    rv_ncx2 = stats.ncx2(df=df, nc=_a_l**2/_var, scale=_var)
    return rv_ncx2.cdf(s)

def integrand_cdf_multi_var(t, s, rho, n_channels=2):
    _a_l = np.sqrt(n_channels*rho*t)
    df = 2*n_channels
    _var = (1-rho)/2
    rv_ncx2 = stats.ncx2(df=df, nc=_a_l**2/_var, scale=_var)
    return rv_ncx2.cdf(s)*np.exp(-t)

def single_channel(N, rho, x0, y0):
    x = stats.norm(scale=np.sqrt(.5)).rvs(N)
    y = stats.norm(scale=np.sqrt(.5)).rvs(N)
    z = np.sqrt(1-rho)*x + np.sqrt(rho)*x0 + 1j*(np.sqrt(1-rho)*y+np.sqrt(rho)*y0)
    abs_z = np.abs(z)**2
    return abs_z

def copula_lower_bound(s, lam_x=1, lam_y=1):
    alpha_x = 1./lam_x
    alpha_y = 1./lam_y
    k = (alpha_x+alpha_y)*np.log(alpha_x+alpha_y) - alpha_x*np.log(alpha_x) - alpha_y*np.log(alpha_y)
    return np.maximum(1.-np.exp(-(s-k)/(alpha_x+alpha_y)), 0)

def copula_upper_bound(s, lam_x=1, lam_y=1):
    alpha_x = 1./lam_x
    alpha_y = 1./lam_y
    lam = 1./np.maximum(alpha_x, alpha_y)
    return stats.expon(scale=1./lam).cdf(s)

def outage_prob_integration(s, rho, n_channels=2, **kwargs):
    result = integrate.quad(integrand_cdf_multi_var, 0, np.inf,
                            args=(s, rho, n_channels))
    return result[0]

def outage_prob_monte_carlo(s, rho, n_channels, t):
    _cdf = cdf_multi_var(s, rho, t, n_channels)
    _prob = np.mean(_cdf)
    return _prob

def main(snr_db=0, rate=1):
    snr = 10**(snr_db/10.)
    s = (2**rate-1)/snr
    # Copula-based bounds
    upper = copula_upper_bound(s)
    lower = copula_lower_bound(s)
    print("Upper Bound: {}\nLower Bound: {}".format(upper, lower))

    # Correlation model from literature
    N = 50000
    n_channels = 2
    rho = np.linspace(0, 1, 20)
    functions = {"integration": outage_prob_integration, "mc": outage_prob_monte_carlo}
    results = {k: [] for k in functions.keys()}
    x0 = stats.norm(scale=np.sqrt(.5)).rvs(N)
    y0 = stats.norm(scale=np.sqrt(.5)).rvs(N)
    t = x0**2 + y0**2
    #t = stats.expon(scale=1).rvs(N)
    for _rho in rho:
        for name, _func in functions.items():
            if _rho == 1.:
                _prob = stats.expon(scale=2).cdf(s)
            else:
                _prob = _func(s=s, rho=_rho, n_channels=n_channels, t=t)
            print("Rho: {:.3f}\tProb: {:.7f}\tAlg: {}".format(_rho, _prob, name))
            results[name].append(_prob)
    results["rho"] = rho
    export_results(results, "rayleigh-comparison-corr-snr{}-rate{}.dat".format(snr_db, rate))
    plt.hlines(lower, min(rho), max(rho))
    plt.hlines(upper, min(rho), max(rho))
    for name, _results in results.items():
        plt.plot(rho, _results, label=name)
    plt.legend()
    plt.xlabel("Rho $\\rho$")
    plt.ylabel("Probability")
    plt.title("Two Rayleigh Fading Channels\nSNR={}dB, Rate={:.3f}, s=(2^rate-1)/snr={:.3f}".format(snr_db, rate, s))

def check_ncx2_dist(pdf=False):
    N = 10000
    x0 = 1.2
    y0 = .5
    rho = .9
    n_channels = 2
    abs_z = [single_channel(N, rho, x0, y0) for i in range(n_channels)]
    sum_abs_z = np.sum(abs_z, axis=0)
    s = np.linspace(0, np.max(sum_abs_z), 75)
    if pdf:
        #_pdf = pdf_single_var(s, rho, x0**2+y0**2)
        _pdf = pdf_multi_var(s, rho, x0**2+y0**2, n_channels)
        plt.hist(sum_abs_z, bins=50, density=True)
        plt.plot(s, _pdf, 'r-')
    else:
        _cdf = cdf_multi_var(s, rho, x0**2+y0**2, n_channels)
        plt.hist(sum_abs_z, bins=150, density=True, cumulative=True)
        plt.plot(s, _cdf, 'r-')

if __name__ == "__main__":
    main(snr_db=10, rate=1)
    #check_ncx2_dist()
    plt.show()
