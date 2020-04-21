"""Basic functions related to copulas.

This module contains different functions for basic copula calculations and
simulations.


Copyright (C) 2020 Karl-Ludwig Besser

This program is used in the article:
Karl-Ludwig Besser and Eduard Jorswieck, "Copula-Based Multi-User Performance
Bounds", submitted to IEEE Communications Letters.

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

class Copula:
    def __init__(self, dim=2, *args, **kwargs):
        self.dim = dim
    
    def _check_input(self, x, check_shape=True):
        x = np.array(x)
        if check_shape:
            pass
            #if np.shape(x)[-1] != self.dim:
            #    raise ValueError("The shape of the input has to match the dimension")
        if not np.all(np.where((0 <= x) & (x <= 1), True, False)):
            raise ValueError("All inputs need to be between 0 and 1.")
        return x
    
    def cdf(self, u):
        pass
    
    def rvs(self, size=1):
        pass

class CopulaFrechetUpper(Copula):
    def cdf(self, u):
        u = self._check_input(u)
        return np.min(u, axis=0)

    def rvs(self, size=1):
        u = np.random.rand(size, 1)
        return np.hstack([u]*self.dim)


class CopulaFrechetLower(Copula):
    def __init__(self, dim=2, *args, **kwargs):
        super().__init__(dim, *args, **kwargs)
        self.dim = 2
    
    def cdf(self, u):
        u = self._check_input(u, check_shape=False)
        return np.maximum(1 - self.dim + np.sum(u, axis=0), 0)

    def rvs(self, size=1):
        if self.dim != 2:
            raise ValueError("The sampling for this copula is currently only supported for 2 dimensions.")
        u = np.random.rand(size, 1)
        v = 1.-u
        return np.hstack((u, v))

    
class CopulaProduct(Copula):
    def cdf(self, u):
        u = self._check_input(u, check_shape=False)
        return np.prod(u, axis=0)
    
    def rvs(self, size=1):
        return np.random.rand(size, self.dim)
    

class CopulaGaussian(Copula):
    def __init__(self, dim=2, cov=None, *args, **kwargs):
        self.cov = cov
        super().__init__(dim=dim, *args, **kwargs)

    def cdf(self, u):
        u = np.array(u)
        u = self._check_input(u)
        is_meshgrid = (u.ndim == 3)
        if is_meshgrid:
            _shape = np.shape(u)[1:]
            u = np.reshape(u, (-1, int(u.size/len(u)))).T
        if isinstance(self.cov, (float, int)):
            cov = np.array([[1, self.cov], [self.cov, 1]])
        x = stats.norm.ppf(u)
        _cdf = stats.multivariate_normal.cdf(x, cov=cov, allow_singular=True)
        if is_meshgrid:
            _cdf = np.reshape(_cdf, _shape)
        return _cdf
    
    def rvs(self, size=1, cov=0):
        if self.cov is not None:
            cov = self.cov
        if isinstance(cov, (float, int)):
            cov = np.array([[1, cov], [cov, 1]])
        #mv_norm_dist = stats.multivariate_normal(cov=cov)
        gauss_samples = stats.multivariate_normal.rvs(size=size, cov=cov)
        uni_samples = stats.norm.cdf(gauss_samples)
        return uni_samples
