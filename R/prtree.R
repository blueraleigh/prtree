prtree = function(x, m, lambda, sigma, maxit)
{
    n = nrow(x)
    d = ncol(x)
    m = min(m, n)
    b = matrix(0L, m, m)
    r = matrix(0, n, m)
    if (m < n)
        v = kmeans(x, x[sample(1:n, m), ])$centers
    else
        v = 1*x
    x = t(x)
    v = t(v)
    converged = .Call(
        C_prtree,
        as.integer(maxit),
        as.numeric(lambda),
        as.numeric(sigma),
        x,
        v,
        r,
        b)
    structure(
        list(r=r, v=t(v), b=b, x=t(x), converged=converged),
        class="prtree"
    )
}