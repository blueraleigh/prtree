#include "prtree.h"

/*
** Return the number of distance calculations needed 
** to compare n objects.
*/
int ndist(int n)
{
    return n * (n-1) / 2;
}


/* 
** Return the squared distance between two n-length vectors 
*/
double dist2(int n, const double *x, const double *y)
{
    double d2 = 0;
    for (--n; n >= 0; --n)
        d2 += (x[n] - y[n]) * (x[n] - y[n]);
    return d2;
}


/* 
** Compute squared distances between all columns of the d-by-n
** matrix x and sort them in increasing order.
*/
void dist(int d, int n, const double *x, double *dij, int *dij_idx)
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
