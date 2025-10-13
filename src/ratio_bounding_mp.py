from mpmath import mpf, log, expm1
from src.utils.utils import logsumexp

def ratio_bounding_mp(f, theta, M, L, eps, initial_k):
    k = initial_k
    L = mpf(L)
    eps = mpf(eps)
    leps = log(eps)
    bucket_size = 5
    converged = False

    def term(k): return mpf(f(theta, k))

    # Monotonicity test
    if L == 0:
        is_decreasing = True
    else:
        is_decreasing = term(M) - term(M-1) > log(L)

    # M sufficiency
    if ((is_decreasing and term(M) - log(-expm1(term(M) - term(M-1))) >= log(2) + leps)
        or (not is_decreasing and term(M) + log(L) - log(1 - L) >= log(2) + leps)):
        raise ValueError("It is not possible to reach the stopping criterion with the given M.")

    log_terms = [log(mpf(0))] * bucket_size
    bucket_iten = 1
    log_terms[0] = term(k)
    log_terms[1] = term(k+1)
    k += 1

    prev_term = log_terms[0]
    current_term = log_terms[1]
    log_sum_total = log(mpf(0))

    while not converged and k < M+initial_k:
        bucket_iten += 1
        k += 1
        prev_term = current_term
        current_term = term(k)
        log_terms[bucket_iten] = current_term

        if bucket_iten + 1 >= bucket_size:
            if (current_term <= prev_term and
                ((is_decreasing and current_term - log(-expm1(current_term - prev_term)) <= log(2) + leps) or 
                (not is_decreasing and current_term + log(L) - log(1 - L) <= log(2) + leps))):
                converged = True

            log_sum_total = logsumexp([log_sum_total] + log_terms[:bucket_iten+1])
            bucket_size *= 2
            log_terms = [log(mpf(0))] * bucket_size
            bucket_iten = -1  # restart fresh
    
    Bound1 = current_term + log(L) - log(1 - L)
    Bound2 = current_term - log(-expm1(current_term - prev_term))
    result = logsumexp([log_sum_total, Bound1 - log(2), Bound2 - log(2)])
    return (k - initial_k, result)
