# to root
import sys
import os
import time
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.ratio_bounding_mp import ratio_bounding_mp
from mpmath import mp, mpf, log, exp, loggamma, fsum, pi
from rpy2.robjects.packages import importr
from src.utils.utils import logdiffexp
from tabulate import tabulate

# function of each term in log scale
def f(theta: tuple, k: int):
    return (mpf(k)) * theta[0] - theta[1] * loggamma(mpf(k)+1)

# dcom from compoisson
def z_dcom(lamb, nu, log_error = 0.001):
    def log_sum(x, y):
        if (x == - mp.inf):
            return y
        elif (y == - mp.inf):
            return x
        elif (x > y):
            return (x + log(1 + exp(y - x)))
        else:
            return (y + log(1 + exp(x - y)))

    lamb = mpf(lamb)
    nu = mpf(nu)
    log_error = mpf(log_error)

    if (lamb < 0 or nu < 0):
        raise ValueError("Invalid arguments, only defined for lambda >= 0, nu >= 0")

    z = - mp.inf
    z_last = mpf(0)
    j = 0

    while(abs(z - z_last) > log_error):
        z_last = z
        z = log_sum(z, j * log(lamb) - nu * loggamma(mpf(j)+1))
        j = j+1
    return z

# dcomp from CompGLM
def z_dcomp(lamb, nu, sumTo = 100):
    sum = mpf(1)
    factorial = mpf("1")
    lamb_pow = mpf(lamb)
    for i in range(1, sumTo+1):
        factorial *= i
        sum += lamb_pow / mp.power(factorial, nu)
        lamb_pow *= lamb
    return sum

# Gaunt's approximation with 4 terms
def z_gaunt(log_lambda, nu):
    log_lambda = mpf(log_lambda)
    nu = mpf(nu)
    nu2 = nu**2
    log_common = log(nu) + log_lambda/nu
    resids = [log(mpf(0))] * 4
    lcte = (nu * exp(log_lambda/nu)) - ((nu - 1)/(2*nu)*log_lambda + (nu - 1)/2*log(2*pi) + mpf('0.5')*log(nu))
    c1 = (nu2 - 1) / 24
    c2 = (nu2 - 1) / 1152 * (nu2 + 23)
    c3 = (nu2 - 1) / 414720 * (5*nu2**2 - 298*nu2 + 11237)
    resids[0] = 1
    resids[1] = c1 * exp(-1 * log_common)
    resids[2] = c2 * exp(-2 * log_common)
    resids[3] = c3 * exp(-3 * log_common)
    return lcte + log(fsum(resids))

def time_mean(func, *args, runs=3):
    """Run func(*args) 3 times and return the mean execution time (seconds)."""
    times = []
    for _ in range(runs):
        start = time.perf_counter()
        func(*args)
        end = time.perf_counter()
        times.append(end - start)
    return (sum(times) / len(times)) * 1000

if __name__ == "__main__":
    mp.dps = 200

    mu = [mpf(4), mpf(4), mpf(4), mpf(4), mpf(4)]
    nu = [mpf('0.1'), mpf('0.3'), mpf('1'), mpf('3'), mpf('30')]
    lamb = [mu[i]**nu[i] for i in range(len(mu))]

    loglamb = [log(x) for x in lamb]
    initial_k = 0
    M = len(nu) * [10**5]
    error = mpf(2)**mpf(-52)

    comp_reg = importr('COMPoissonReg')

    data = []
    for i in range(len(lamb)):
        # Measure time for each method (mean of 3 runs)
        t_rbp_loose = time_mean(ratio_bounding_mp, f, (loglamb[i], nu[i]), M[i], mpf(0), 10**(-4), initial_k)
        t_rbp_strict = time_mean(ratio_bounding_mp, f, (loglamb[i], nu[i]), M[i], mpf(0), 10**(-16), initial_k)
        t_dcmp = time_mean(comp_reg.dcmp, 0, float(lamb[i]), float(nu[i]))
        t_dcom = time_mean(z_dcom, float(lamb[i]), float(nu[i]))
        t_dcomp = time_mean(z_dcomp, float(lamb[i]), float(nu[i]))
        t_gaunt = time_mean(z_gaunt, float(loglamb[i]), float(nu[i]))

        data.append([
            f"mu={mu[i]} | nu={nu[i]}",
            t_rbp_loose,
            t_rbp_strict,
            t_dcmp,
            t_dcom,
            t_dcomp,
            t_gaunt
        ])

    headers = ["", "2.2x10^-4|RBP", "2.2x10^-16|RBP", "COMPoissonReg", "compoisson", "CompGLM", "Gaunt's Approx."]
    print(tabulate(data, headers, tablefmt="fancy_grid", floatfmt=".6f"))
