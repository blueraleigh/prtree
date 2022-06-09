#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>

// Given a flat index (k) into the lower triangle of a n-by-n matrix
// stored in column major order, return the corresponding row (i) and 
// column (j) indices. the diagonal is not included as part of the
// lower triangle. e.g.,
//
//   - - - - -    0 -> (1, 0)  5 -> (3, 1)
//   0 - - - -    1 -> (2, 0)  6 -> (4, 1)
//   1 4 - - -    2 -> (3, 0)  7 -> (3, 2)
//   2 5 7 - -    3 -> (4, 0)  8 -> (4, 2)
//   3 6 8 9 -    4 -> (2, 1)  9 -> (4, 3)
//
// Reference:
//
// Modified from Algorithm 3 of the following reference
//   https://hal.archives-ouvertes.fr/hal-02047514/document
//
// The linked reference uses a 1-based index and includes the diagonal. 
// Here we use a 0-based index and exclude the diagonal.
//
// See also: https://en.wikipedia.org/wiki/Triangular_number
static void k2ij(int n, int k, int *i, int *j)
{
    ++k;
    --n;
    int kp = ( n * (n + 1) ) / 2 - k;
    int p = ( sqrt(1 + 8 * kp) - 1 ) / 2;
    *i = k - n*(n-1)/2 + p*(p+1)/2;
    *j = n - p - 1;
}


// return the number of distance calculations
// needed to compare n objects
static int ndist(int n)
{
    return n * (n-1) / 2;
}


// return the squared distance between two n-length vectors
static double dist2(int n, double *x, double *y)
{
    double d2 = 0;
    for (--n; n >= 0; --n)
        d2 += (x[n] - y[n]) * (x[n] - y[n]);
    return d2;
}


// compute squared distances between all columns of the d-by-n matrix x
// and sort them in increasing order
static void dist(int d, int n, double *x, double *dij, int *dij_idx)
{
    int k;
    int i;
    int j;
    int nd = ndist(n);
    for (k = 0; k < nd; ++k)
    {
        k2ij(n, k, &i, &j);
        dij[k] = dist2(d, x+i*d, x+j*d);
        dij_idx[k] = k;
    }
    rsort_with_index(dij, dij_idx, nd); // R_ext/Utils.h
}


// find the root of x
static int mst_find(int *parent, int x)
{
    while (parent[x] != x)
    {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}


static int mst_union(int *parent, int *size, int x, int y)
{
    x = mst_find(parent, x);
    y = mst_find(parent, y);
    if (x != y)
    {
        // form an edge between x and y
        if (size[x] < size[y])
        {
            int z = x;
            x = y;
            y = z; 
        }
        parent[y] = x;
        size[x] += size[y];
        return 1;
    }
    // x and y share a root -- connecting them
    // by an edge would make a cycle (disallowed)
    return 0;
}


// compute the minimum spanning tree b (and return its length)
static double mst(
    int n, double *dij, int *dij_idx, int *b, int *parent, int *size)
{
    int k;
    int i;
    int j;
    int nd = ndist(n);
    double len = 0;
    
    // initialize the forest
    for (k = 0; k < n; ++k)
    {
        size[k] = 1;
        parent[k] = k;
    }

    for (k = 0; k < nd; ++k)
    {
        k2ij(n, dij_idx[k], &i, &j);

        if (mst_union(parent, size, i, j))
        {
            b[i+j*n] = b[j+i*n] = 1;
            len += dij[k];
        }
        else
        {
            b[i+j*n] = b[j+i*n] = 0;
        }

        b[i+i*n] = b[j+j*n] = 0;
    }

    return len;
}


// compute the laplacian of the minimum spanning tree b
static void lapl(int n, int *b, double *l)
{
    int k;
    int i;
    int j;
    int nd = ndist(n);
    memset(l, 0, n*n*sizeof(double));
    for (k = 0; k < nd; ++k)
    {
        k2ij(n, k, &i, &j);
        if (b[i+j*n])
        {
            l[i+j*n] = -1;
            l[j+i*n] = -1;
            l[i+i*n] += 1;
            l[j+j*n] += 1;
        }
    }
}


// sigma is the bandwith used for computing membership coefficients
// x is a d-by-n matrix of data coordinates
// v is a d-by-m matrix of vertex coordinates
// r is a n-by-m matrix of membership coefficients
// c is a m-length vector containing the column sums of r
// z is a n-length vector for temporary row sums
static double softmax(int d, int n, int m, double sigma, double *x, double *v, 
    double *r, double *c, double *z)
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


// update the embedding coordinates
static void embed(int d, int n, int m, double lam, 
    double *x, double *v, double *l, double *r, double *c, 
    double *xr, int lwork, double *work, int *ipiv)
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

    // invert l (symmetric)
    F77_CALL(dsytrf)("L",&m,l,&m,ipiv,work,&lwork,&info);
    F77_CALL(dsytri)("L",&m,l,&m,ipiv,work,&info);

    // A*B where A=xr, B=l(inverse,symmetric)
    F77_CALL(dsymm)("R","L",&d,&m,&one,l,&m,xr,&d,&zero,v,&d);
}


static int fit(
    int d,
    int n,
    int m,
    int maxit,
    double lambda,
    double sigma,
    double *x,  // d-by-n matrix of data coordinates
    double *v,  // d-by-m matrix of initial vertex coordinates
    double *r,  // n-by-m matrix of membership probabilities
    int *b      // m-by-m adjacency matrix with mst
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
    double mse;

    // score to minimize
    double score0;
    double score1;
    dist(d, m, v, dij, dij_idx);
    len = mst(m, dij, dij_idx, b, parent, size);
    lapl(m, b, l);
    mse = softmax(d, n, m, sigma, x, v, r, c, z);
    score0 = mse + lambda * len / 2;
    Rprintf("%-14s %-14s\n", "Iter", "Score");
    Rprintf("%-14d %-14f\n", 0, len+mse);
    for (i = 0, converged = 0; i < maxit && !converged; ++i)
    {
        embed(d, n, m, lambda, x, v, l, r, c, xr, lwork, work, ipiv);
        dist(d, m, v, dij, dij_idx);
        len = mst(m, dij, dij_idx, b, parent, size);
        lapl(m, b, l);
        mse = softmax(d, n, m, sigma, x, v, r, c, z);
        score1 = mse + lambda * len / 2;
        converged = ((score0 - score1) / score0) < 0.0001;
        score0 = score1;
        Rprintf("%-14d %-14f\n", i+1, len+mse);
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
    SEXP x,  // d-by-n matrix of data coordinates
    SEXP v,  // d-by-m matrix of initial vertex coordinates
    SEXP r,  // n-by-m matrix of membership probabilities
    SEXP b   // m-by-m adjacency matrix with mst
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


/*
** Force-directed graph layout as described in:
**
** Graph Drawing by Stress Majorization
** Emden R. Gansner, Yehuda Koren and Stephen North
** In Proceedings 12th Symposium on Graph Drawing (GD), 
** pages 239–250, 2004
**
** https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.439.26&rep=rep1&type=pdf
** 
*/

static double stress(int d, int n, double *x, double *dij)
{
    int i;
    int j;
    int k;
    int nd = ndist(n);
    double d0;
    double d1;
    double s = 0;
    for (k = 0; k < nd; ++k)
    {
        k2ij(n, k, &i, &j);
        d0 = dij[k];
        d1 = sqrt(dist2(d, x+i*d, x+j*d));
        s += ( (d1 - d0) * (d1 - d0) ) / (d0 * d0);
    }
    return s;
}


static void update_layout(int d, int n, double *Lw, double *Lx, double *x,
    double *z)
{

    int i;
    int info;
    int nminus1 = n - 1;
    double one = 1;
    double zero = 0;
    
    // A*B where A=Lx, B=t(x)
    // Lx is n-by-n
    // x is d-by-n
    F77_CALL(dgemm)("N","T",&n,&d,&n,&one,Lx,&n,x,&d,&zero,z,&n);

    F77_CALL(dpotrs)("L",&nminus1,&d,Lw,&n,z,&n,&info);

    for (i = 0; i < nminus1; ++i)
    {
        x[0+i*d] = z[i+0*n];
        x[1+i*d] = z[i+1*n];
    }
}

static void dijkstra_(int d, int m, int *b, 
    double *x, double *dij, int *notvisited, int s)
{
    int v;
    int u = s;
    int nQ = m - 1;
    double alt;

    for (v = 0; v < m; ++v)
    {
        notvisited[v] = 1;
        dij[s+v*m] = R_PosInf;
    }
    dij[s+s*m] = 0;

    do 
    {
        notvisited[u] = 0;
        for (v = 0; v < m; ++v)
        {
            if (b[u+v*m] && notvisited[v] && dij[s+u*m] < R_PosInf)
            {
                alt = dij[s+u*m] + sqrt(dist2(d, x+u*d, x+v*d));
                if (alt < dij[s+v*m])
                    dij[s+v*m] = alt;
            }
        }
        alt = R_PosInf;
        for (v = 0; v < m; ++v)
        {
            if (notvisited[v] && dij[s+v*m] <= alt)
            {
                u = v;
                alt = dij[s+v*m];
            }
        }
    } while (nQ-- > 0);
}

static double *dijkstra(int d, int m, int *b, double *x)
{
    int i;
    int j;
    int k;
    double *dij = malloc(m * m * sizeof(double));
    double *D = malloc((m*(m-1)/2) * sizeof(double));
    int *notvisited = malloc(m * sizeof(*notvisited));
    
    for (i = 0; i < m; ++i)
        dijkstra_(d, m, b, x, dij, notvisited, i);

    for (j = 0, k = 0; j < m; ++j)
    {
        for (i = j+1; i < m; ++i)
            D[k++] = dij[i+j*m];
    }

    free(dij);
    free(notvisited);
    return D;
}

SEXP C_prtree_layout(SEXP layout, SEXP b, SEXP v)
{
    int i;
    int j;
    int k;
    int d = INTEGER(getAttrib(v, R_DimSymbol))[0];
    int m = INTEGER(getAttrib(v, R_DimSymbol))[1];
    int mminus1 = m - 1;

    int nd = ndist(m);

    double d0;
    double d1;
    double s0;
    double s1;
    double delta;
    double eps = 0.0001;
    
    // calculate shortest path distances between all vertices
    double *dij = dijkstra(d, m, INTEGER(b), REAL(v));

    // layout coordinates
    double *x = REAL(layout);
    d = 2;

    // laplacian matrices
    double *Lw = calloc(m * m, sizeof(double));
    double *Lx = calloc(m * m, sizeof(double));

    double *z = malloc(m * d * sizeof(double));

    for (k = 0; k < nd; ++k)
    {
        k2ij(m, k, &i, &j);
        Lw[i+j*m] = -1 / (dij[k] * dij[k]);
        Lw[j+i*m] = -1 / (dij[k] * dij[k]);
        Lw[i+i*m] += 1 / (dij[k] * dij[k]);
        Lw[j+j*m] += 1 / (dij[k] * dij[k]);
    }

    // cholesky factorization (dropping the last row and column)
    int info;
    F77_CALL(dpotrf)("L", &mminus1, Lw, &m, &info);

    s0 = stress(d, m, x, dij);
    int s = 0;
    Rprintf("Stress %d: %14f\n", s++, s0);
    do {

        memset(Lx, 0, m*m*sizeof(double));
        for (k = 0; k < nd; ++k)
        {
            k2ij(m, k, &i, &j);
            d0 = dij[k];
            d1 = sqrt(dist2(d, x+i*d, x+j*d));
            if (d1 != 0)
            {
                Lx[i+j*m] = -(1/d0)*(1/d1);
                Lx[j+i*m] = -(1/d0)*(1/d1);
                Lx[i+i*m] += (1/d0)*(1/d1);
                Lx[j+j*m] += (1/d0)*(1/d1);
            }
        }
        update_layout(d, m, Lw, Lx, x, z);
        s1 = stress(d, m, x, dij);
        delta = (s0 - s1) / s0;
        s0 = s1;
        Rprintf("Stress %d: %14f\n", s++, s1);
    } while (delta > eps);

    free(dij);
    free(Lw);
    free(Lx);
    free(z);

    return R_NilValue;
}

