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
    
    def cdf(self, x):
        pass
    
    def rvs(self, size=1):
        pass

class CopulaFrechetUpper(Copula):
    def cdf(self, x):
        x = self._check_input(x)
        return np.min(x, axis=0)

    def rvs(self, size=1):
        u = np.random.rand(size, 1)
        return np.hstack([u]*self.dim)


class CopulaFrechetLower(Copula):
    def __init__(self, dim=2, *args, **kwargs):
        super().__init__(dim, *args, **kwargs)
        self.dim = 2
    
    def cdf(self, x):
        x = self._check_input(x, check_shape=False)
        return np.maximum(1 - self.dim + np.sum(x, axis=0), 0)

    def rvs(self, size=1):
        if self.dim != 2:
            raise ValueError("The sampling for this copula is currently only supported for 2 dimensions.")
        u = np.random.rand(size, 1)
        v = 1.-u
        return np.hstack((u, v))

    
class CopulaProduct(Copula):
    def cdf(self, x):
        x = self._check_input(x, check_shape=False)
        return np.prod(x, axis=0)
    
    def rvs(self, size=1):
        return np.random.rand(size, self.dim)
    

class CopulaGaussian(Copula):
    def __init__(self, dim=2, cov=None, *args, **kwargs):
        self.cov = cov
        super().__init__(dim=dim, *args, **kwargs)

    def cdf(self, x):
        return NotImplemented
    
    def rvs(self, size=1, cov=0):
        if self.cov is not None:
            cov = self.cov
        if isinstance(cov, (float, int)):
            cov = np.array([[1, cov], [cov, 1]])
        #mv_norm_dist = stats.multivariate_normal(cov=cov)
        gauss_samples = stats.multivariate_normal.rvs(size=size, cov=cov)
        uni_samples = stats.norm.cdf(gauss_samples)
        return uni_samples
