import numpy as np
from scipy.stats import norm

def wald_z(diff, se_diff):
    """
    Calculate Wald Z statistic and two-tailed p-value.
     diff: difference between two estimates
     se_diff: standard error of the difference
     Returns Z statistic and p-value
    """
    Z = diff / se_diff
    p = 2 * (1 - norm.cdf(abs(Z)))
    return Z, p

def holm_correction(pvals):
    """
    Holm step-down adjustment; returns adjusted p-values in original order.
    """
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])  # indices sorted by p
    adjusted = [0.0] * m
    prev = 0.0
    for rank, i in enumerate(order):
        k = m - rank
        adj = min(1.0, pvals[i] * k)
        adj = max(adj, prev)  # enforce monotonicity
        adjusted[i] = adj
        prev = adj
    return adjusted

def prop_err_division(res_num, res_denom):
    """
    Calculate propagated error a division of two quantities with
    associated standard errors.
    res_num: (value, std_error) of the numerator
    res_denom: (value, std_error) of the denominator
    """
    return np.sqrt(((1/res_denom[0]) * res_num[1])**2 +
                   ((res_num[0] / (res_denom[0]**2)) * res_denom[1])**2)
