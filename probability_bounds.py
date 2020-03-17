import numpy as np
from scipy import stats

def export_results(results, filename):
    import pandas as pd
    data = pd.DataFrame.from_dict(results)
    data.to_csv(filename, sep="\t", index=False)

def copula_lower_sum_uniform(a, b, a_x=1, b_x=3, a_y=2, b_y=5, s=1):
    _s_min_lower = np.minimum(a_x+b_y, a_y+b_x)
    rv_s_lower = stats.uniform(loc=_s_min_lower, scale=b_x+b_y-_s_min_lower)
    t = rv_s_lower.cdf(s)
    c = np.minimum(a, b)
    idx_t = np.where(np.logical_and(a >= t, b >= t))
    c[idx_t] = np.maximum(a[idx_t] + b[idx_t] - 1, t)
    return c

def joint_pdf_lower_sum_uniform(a_x=1, b_x=3, a_y=2, b_y=5, s=6):
    n_samples = 50
    xlim = [0, 4]
    ylim = [0, 6]
    x, stepx = np.linspace(*xlim, num=n_samples, retstep=True)
    y, stepy = np.linspace(*ylim, num=n_samples, retstep=True)
    X, Y = np.meshgrid(x, y)
    rv_x = stats.uniform(loc=a_x, scale=b_x-a_x)
    rv_y = stats.uniform(loc=a_y, scale=b_y-a_y)
    marg_cdf_x = rv_x.cdf(X)
    marg_cdf_y = rv_y.cdf(Y)
    joint_cdf = copula_lower_sum_uniform(marg_cdf_x, marg_cdf_y,
                                         a_x=a_x, b_x=b_x, a_y=a_y, b_y=b_y, s=s)
    _gradx = np.gradient(joint_cdf, stepx, axis=0)
    joint_pdf = np.gradient(_gradx, stepy, axis=1)
    filename = "joint_pdf-sum_uniform_lower-X{}_{}-Y{}_{}-s{}.dat".format(
                a_x, b_x, a_y, b_y, s)
    print(np.min(joint_pdf), np.max(joint_pdf))
    results = {"X": X.ravel(), "Y": Y.ravel(), "pdf": joint_pdf.ravel()}
    export_results(results, filename)



def _g_upper_product_uniform(y, s, rv_x, rv_y):
    _part1 = np.minimum(rv_x.cdf(s/y)+rv_y.cdf(y), rv_y.cdf(s/y)+rv_x.cdf(y))
    return np.minimum(_part1, 1)

def upper_product_uniform(s, a_x=0, b_x=1, a_y=0, b_y=1):
    rv_x = stats.uniform(loc=a_x, scale=b_x-a_x)
    rv_y = stats.uniform(loc=a_y, scale=b_y-a_y)
    yopt = np.sqrt(s*(b_y-a_y)/(b_x-a_x))
    return np.minimum(_g_upper_product_uniform(yopt, s, rv_x, rv_y), np.minimum(rv_x.cdf(s/a_y), rv_y.cdf(s/a_x)))

def save_upper_product_uniform(**kwargs):
    _upper = upper_product_uniform(**kwargs)
    results = {"s": kwargs.get('s'), "upper": _upper}
    filename = "product_uniform_upper-X{a_x}_{b_x}-Y{a_y}_{b_y}.dat".format(**kwargs)
    export_results(results, filename)

if __name__ == "__main__":
    #joint_pdf_lower_sum_uniform()
    save_upper_product_uniform(s=np.linspace(0, 20, num=100),
                               a_x=1, b_x=3,
                               a_y=2, b_y=5)
