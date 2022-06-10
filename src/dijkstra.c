#include "prtree.h"

static void dijkstra_(
    int d
    , int m
    , const int *b
    , const double *x
    , double *dij
    , int *notvisited
    , int s)
{
    int v;
    int u = s;
    int n = m - 1;
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
    } while (n-- > 0);
}

/*
** Compute shortest graph distances between all vertices using
** Dijkstra's algorithm.
**
** [in] d: Coordinate dimension
** [in] m: Number of vertices
** [in] b: A m-by-m binary adjaceny matrix
** [in] x: A d-by-m matrix of vertex coordinates
**
** Returns a vector of pairwise distances that must be free'd later.
** The distances are the lower triangular portion of the full
** m-by-m distance matrix stored in column-major order.
*/
double *dijkstra(int d, int m, const int *b, const double *x)
{
    int i;
    int j;
    int k;
    double *Dij = malloc(m * m * sizeof(*Dij));
    double *dij = malloc((m*(m-1)/2) * sizeof(*dij));
    int *notvisited = malloc(m * sizeof(*notvisited));
    
    for (i = 0; i < m; ++i)
        dijkstra_(d, m, b, x, Dij, notvisited, i);

    for (j = 0, k = 0; j < m; ++j)
    {
        for (i = j+1; i < m; ++i)
            dij[k++] = Dij[i+j*m];
    }

    free(Dij);
    free(notvisited);

    return dij;
}
