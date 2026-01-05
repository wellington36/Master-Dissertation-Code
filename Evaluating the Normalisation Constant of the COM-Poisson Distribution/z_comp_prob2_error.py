# to root
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from src.ratio_bounding_mp import ratio_bounding_mp

from mpmath import mp, mpf, log, exp, loggamma, fsum, pi, factorial
from rpy2.robjects.packages import importr
from src.utils.utils import logdiffexp
from tabulate import tabulate
import math


# function of each term in log scale
def f(theta: tuple, k: int):
    """
    terms of the normalization contant series

    theta   : (log_lambda, nu)
    k       : k-th term
    """
    return (mpf(k)) * theta[0] - theta[1] * loggamma(mpf(k)+1)

# dcom from compoisson (https://cran.r-project.org/src/contrib/Archive/compoisson)
# https://cran.r-project.org/web/packages/compoisson/index.html
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

def z_dcomp(lamb, nu, sumTo = 100):
    sum = mpf(1)
    factorial = mpf("1")
    lamb_pow = mpf(lamb)

    for i in range(1, sumTo+1):
        factorial *= i
        sum += lamb_pow / mp.power(factorial, nu)
        lamb_pow *= lamb
    
    return sum

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


if __name__ == "__main__":
    mp.dps = 200

    prob = 2

    mu = [mpf(4), mpf(4), mpf(4), mpf(4), mpf(4)]
    nu = [mpf('0.1'), mpf('0.3'), mpf('1'), mpf('3'), mpf('30')]
    lamb = [mu[i]**nu[i] for i in range(len(mu))]

    loglamb = [log(x) for x in lamb]
    initial_k = 0
    M = len(nu) * [10**5]
    gold_error = mpf(2)**mpf(-64)

    
    error = 10**(-4)
    error_m4 = []
    for i in range(len(lamb)):
        gold_value = ratio_bounding_mp(f, (loglamb[i], nu[i]), M[i], mpf(0), gold_error, initial_k=initial_k)[1]
        bp_value = ratio_bounding_mp(f, (loglamb[i], nu[i]), M[i], mpf(0), error, initial_k=initial_k)[1]

        prob2_gold_value = prob * loglamb[i] - nu[i] * log(factorial(prob)) - gold_value
        prob2_bp_value = prob * loglamb[i] - nu[i] * log(factorial(prob)) - bp_value


        error_m4.append([exp(logdiffexp(prob2_bp_value, prob2_gold_value))])
        #error_m4.append([exp(prob2_bp_value)])

    
    error = 10**(-16)
    error_m16 = []
    for i in range(len(lamb)):
        gold_value = ratio_bounding_mp(f, (loglamb[i], nu[i]), M[i], mpf(0), gold_error, initial_k=initial_k)[1]
        bp_value = ratio_bounding_mp(f, (loglamb[i], nu[i]), M[i], mpf(0), error, initial_k=initial_k)[1]

        prob2_gold_value = prob * loglamb[i] - nu[i] * log(factorial(prob)) - gold_value
        prob2_bp_value = prob * loglamb[i] - nu[i] * log(factorial(prob)) - bp_value

        error_m16.append([exp(logdiffexp(prob2_bp_value, prob2_gold_value))])
        #error_m16.append([exp(prob2_bp_value)])
    
    # Libraries
    comp_reg = importr('COMPoissonReg')

    def dcmp_in_log_scale(x, lambda_, nu):
        # Call the dcmp function with log=TRUE in R
        result = comp_reg.dcmp(x, lambda_, nu)
        return float(result[0])  # Return the first result as a Python float

    libraries = []
    for i in range(len(lamb)):
        gold_value = ratio_bounding_mp(f, (loglamb[i], nu[i]), M[i], mpf(0), gold_error, initial_k=initial_k)[1]
        prob2_gold_value = prob * loglamb[i] - nu[i] * log(factorial(prob)) - gold_value

        z_dcom_value = z_dcom(float(lamb[i]), float(nu[i]))
        z_dcomp_value = log(z_dcomp(float(lamb[i]), float(nu[i])))
        z_gaunt_value = z_gaunt(log(float(lamb[i])), float(nu[i]))

        prob2_dcmp = log(list(comp_reg.dcmp(prob, float(lamb[i]), float(nu[i])))[0])
        prob2_dcom = prob * loglamb[i] - nu[i] * log(factorial(prob)) - z_dcom_value
        prob2_dcomp = prob * loglamb[i] - nu[i] * log(factorial(prob)) - z_dcomp_value
        prob2_gaunt = prob * loglamb[i] - nu[i] * log(factorial(prob)) - z_gaunt_value

        libraries.append([exp(logdiffexp(prob2_gold_value, p)) for p in [prob2_dcmp,
                                                                         prob2_dcom,
                                                                         prob2_dcomp,
                                                                         prob2_gaunt]])
        #libraries.append([exp(prob2_dcmp), exp(prob2_dcom), exp(prob2_dcomp), exp(prob2_gaunt)])

    
    # Organize in a table
    data = []
    for idx, (a, b, c) in enumerate(zip(error_m4, error_m16, libraries), start=1):
        data.append([f"mu={float(mu[idx-1]):.2f} | nu={float(nu[idx-1]):.2f}", a[0], b[0], c[0], c[1], c[2], c[3]])

    headers = ["", "2.2x10^-4|RBP", "2.2x10^-16|RBP", "COMPoissonReg", "compoisson", "CompGLM", "Gaunt's Approx."]

    print(tabulate(data, headers, tablefmt="fancy_grid"))
    