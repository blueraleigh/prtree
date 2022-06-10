#include "prtree.h"


static double stress(int d, int n, const double *x, const double *dij)
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

/* Update the layout coordinates */
static void layout(
    int d
    , int n
    , const double *Lw
    , double *Lx
    , double *x
    , double *z)
{
    int i;
    int info;
    int nminus1 = n - 1;
    double one = 1;
    double zero = 0;
    
    // x is d-by-n
    // Lx is n-by-n
    // A*B where A=Lx, B=t(x)
    F77_CALL(dgemm)("N","T",&n,&d,&n,&one,Lx,&n,x,&d,&zero,z,&n);

    // Lw is nminus1-by-nminus1 but the
    // leading storage dimension is n
    F77_CALL(dpotrs)("L",&nminus1,&d,Lw,&n,z,&n,&info);

    // the last row of x is fixed (0,0)
    for (i = 0; i < nminus1; ++i)
    {
        // record the new layout coordinates
        x[0+i*d] = z[i+0*n];
        x[1+i*d] = z[i+1*n];
        // reset diagonal of Lx for next update
        Lx[i+i*n] = 0;
    }
    Lx[nminus1+nminus1*n] = 0;
}


/*
** Compute the force-directed graph layout as described in:
**
** Graph Drawing by Stress Majorization
** Emden R. Gansner, Yehuda Koren and Stephen North
** In Proceedings 12th Symposium on Graph Drawing (GD), 
** pages 239–250, 2004
**
** https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.439.26&rep=rep1&type=pdf
** 
** [in,out] X: Planar (x,y) coordinates for each vertex.
** [in]     b: Binary adjaceny matrix.
** [in]     v: Multidimensional embedding coordinates for each vertex.
*/
SEXP C_prtree_layout(SEXP X, SEXP b, SEXP v)
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
    double *x = REAL(X);
    d = 2;

    // laplacian matrices
    double *Lw = calloc(m * m, sizeof(double));
    double *Lx = calloc(m * m, sizeof(double));

    // LAPACK workspace
    double *z = malloc(m * d * sizeof(double));

    for (k = 0; k < nd; ++k)
    {
        k2ij(m, k, &i, &j);
        Lw[i+j*m] = -1 / (dij[k] * dij[k]);
        Lw[j+i*m] = -1 / (dij[k] * dij[k]);
        Lw[i+i*m] += 1 / (dij[k] * dij[k]);
        Lw[j+j*m] += 1 / (dij[k] * dij[k]);
    }

    int info;
    // compute the cholesky factorization of Lw
    // after removing the last row and column.
    // this ensures that Lw is positive definite.
    F77_CALL(dpotrf)("L", &mminus1, Lw, &m, &info);

    int s = 0;
    s0 = stress(d, m, x, dij);
    Rprintf("Stress %d: %14f\n", s++, s0);
    do {

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
        layout(d, m, Lw, Lx, x, z); // also resets diagonal of Lx to 0
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
