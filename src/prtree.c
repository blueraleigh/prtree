#include "prtree.h"

/*
** Compute a soft clustering of the data by the vertices of the
** principal tree.
**
** [in]     sigma: The bandwith used for computing membership coefficients.
** [in]     x:     A d-by-n matrix of data coordinates.
** [in]     v:     A d-by-m matrix of vertex coordinates.
** [in,out] r:     A n-by-m matrix of membership coefficients.
** [in,out] c:     A m-length vector containing the column sums of r.
** [in]     z:     A n-length vector for temporary row sums.
**
** Returns the mean quantization error of the principal tree.
*/
static double softmax(
    int d
    , int n
    , int m
    , double sigma
    , const double *x
    , const double *v 
    , double *r
    , double *c
    , double *z)
{
    int i;
    int j;
    int k;
    double e = 0;
    double maxrij;
    /*
    for (i = 0; i < n; ++i)
    {
        z[i] = 0;
        for (j = 0; j < m; ++j)
        {
            c[j] = 0;
            z[i] += r[i+j*n] = exp(-dist2(d, x+i*d, v+j*d) / sigma);
        }
    }
    */

    for (i = 0; i < n; ++i)
    {
        z[i] = R_NegInf;
        for (j = 0; j < m; ++j)
        {
            c[j] = 0;
            r[i+j*n] = -dist2(d, x+i*d, v+j*d) / sigma;
            if (r[i+j*n] > z[i])
                z[i] = r[i+j*n];
        }
    }
    for (i = 0; i < n; ++i)
    {
        maxrij = z[i];
        z[i] = 0;
        for (j = 0; j < m; ++j)
            z[i] += r[i+j*n] = exp(r[i+j*n] - maxrij);
    }

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
            c[j] += r[i+j*n] /= z[i];
    }
    
    // finally, compute the quantization error
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            if (r[i+j*n] > 0)
                e += r[i+j*n]*(dist2(d, x+i*d, v+j*d) + sigma*log(r[i+j*n]));
        }
    }
    return e / n;
}


/* Update the embedding coordinates */
static void embed(
    int d
    , int n
    , int m
    , double lam
    , const double *x
    , double *v
    , double *l
    , double *r
    , double *c
    , double *xr
    , int lwork
    , double *work
    , int *ipiv)
{
    int i;
    int info;
    double one = 1;
    double zero = 0;
    
    // A*B where A=x, B=r
    F77_CALL(dgemm)("N","N",&d,&m,&n,&one,x,&d,r,&n,&zero,xr,&d);

    // augment l
    for (i = 0; i < m*m; ++i)
        l[i] *= lam;
    for (i = 0; i < m; ++i)
        l[i+i*m] += c[i];

    // invert l (which is symmetric)
    F77_CALL(dsytrf)("L",&m,l,&m,ipiv,work,&lwork,&info);
    F77_CALL(dsytri)("L",&m,l,&m,ipiv,work,&info);

    // A*B where A=xr, B=l (now the symmetric inverse of the original l)
    F77_CALL(dsymm)("R","L",&d,&m,&one,l,&m,xr,&d,&zero,v,&d);
}


/*
** Compute a principal tree.
**
** [in]     maxit:  The max number of iterations.
** [in]     lambda: Regularization parameter on length of the tree.
** [in]     sigma:  The bandwith used for soft clustering.
** [in]     x:      A d-by-n matrix of data coordinates.
** [in,out] v:      A d-by-m matrix of vertex coordinates.
** [in,out] r:      A n-by-m matrix of membership coefficients.
** [in,out] b:      A m-by-m adjaceny matrix encoding the principal tree.
**
** Returns 1 if the algorithm converved, 0 otherwise.
*/
static int fit(
    int d,
    int n,
    int m,
    int maxit,
    double lambda,
    double sigma,
    const double *x,
    double *v,
    double *r,
    int *b
)
{
    int i;
    int converged;

    // BLAS/LAPACK workspace
    int info;
    int lwork = -1;
    int *ipiv = malloc(m * sizeof(int));
    double *work = malloc(1 * sizeof(double));
    F77_CALL(dsytrf)("L",&m,0,&m,0,work,&lwork,&info);
    lwork = *work;
    work = realloc(work, lwork * sizeof(double));

    // distance work space
    int *dij_idx = malloc(ndist(m) * sizeof(int));
    double *dij = malloc(ndist(m) * sizeof(double));

    // mst work space
    int *parent = malloc(m * sizeof(int));
    int *size = malloc(m * sizeof(int));

    // laplacian work space
    double *l = malloc(m * m * sizeof(double));

    // softmax work space
    double *c = malloc(m * sizeof(double));
    double *z = malloc(n * sizeof(double));

    // embed workspace
    double *xr = malloc(d * m * sizeof(double));

    // sum of squared edge lengths of the mst
    double len;
    // empirical mean quantization error
    double mqe;

    double score0;
    double score1;
    dist(d, m, v, dij, dij_idx);
    len = mst(m, dij, dij_idx, b, l, parent, size);
    mqe = softmax(d, n, m, sigma, x, v, r, c, z);
    score0 = mqe + lambda * len / 2;
    Rprintf("%-14s %-14s\n", "Iter", "Score");
    Rprintf("%-14d %-14f\n", 0, score0);
    for (i = 0, converged = 0; i < maxit && !converged; ++i)
    {
        embed(d, n, m, lambda, x, v, l, r, c, xr, lwork, work, ipiv);
        dist(d, m, v, dij, dij_idx);
        len = mst(m, dij, dij_idx, b, l, parent, size);
        mqe = softmax(d, n, m, sigma, x, v, r, c, z);
        score1 = mqe + lambda * len / 2;
        converged = ((score0 - score1) / score0) < 0.0001;
        score0 = score1;
        Rprintf("%-14d %-14f\n", i+1, score0);
    }

    free(ipiv);
    free(work);
    free(dij);
    free(dij_idx);
    free(parent);
    free(size);
    free(l);
    free(c);
    free(z);
    free(xr);

    return converged;
}


SEXP C_prtree(
    SEXP maxit,
    SEXP lambda,
    SEXP sigma,
    SEXP x,
    SEXP v,
    SEXP r,
    SEXP b
)
{
    int d = INTEGER(getAttrib(x, R_DimSymbol))[0];
    int n = INTEGER(getAttrib(x, R_DimSymbol))[1];
    int m = INTEGER(getAttrib(v, R_DimSymbol))[1];
    int converged = fit(d, n, m, *INTEGER(maxit), 
        *REAL(lambda), *REAL(sigma), REAL(x), 
        REAL(v), REAL(r), INTEGER(b));
    return ScalarInteger(converged);
}
