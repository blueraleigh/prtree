#include "prtree.h"

/*
** Return the root vertex of x's forest.
*/
static int mst_find(int *parent, int x)
{
    while (parent[x] != x)
    {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}


/*
** Form the union of the two forests to which vertices
** x and y belong. Do nothing if they belong to the same
** forest. Returns 1 if a union operation was performed,
** 0 otherwise.
*/
static int mst_union(int *parent, int *size, int x, int y)
{
    x = mst_find(parent, x);
    y = mst_find(parent, y);
    if (x != y)
    {
        // x and y are currently in different forests,
        // so form an edge between them.
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
    // x and y already belong to the same forest.
    // connecting them by an edge would form a cycle.
    return 0;
}


/*
** Compute the minimum spanning tree using Kruskal's algorithm
**
** [in]     n:           The number of vertices.
** [in]     dij:         Pairwise distances between all vertices sorted in
**                         increasing order.
** [in]     dij_idx:     Permutation vector recording the original indices
**                         of the sorted distance vector.
** [in,out] b:           Binary adjaceny matrix encoding the minimum 
**                         spanning tree.
** [in,out] l:           The Laplacian matrix of b.
** [in]     parent,size: Disjoint set data structures for computing the 
**                         spanning tree.
**
** Returns the length of the minimum spanning tree.
*/                
double mst(
    int n
    , const double *dij
    , const int *dij_idx
    , int *b
    , double *l
    , int *parent
    , int *size)
{
    int k;
    int i;
    int j;
    int nd = ndist(n);
    double len = 0;
    
    // initialize the forest
    for (i = 0; i < n; ++i)
    {
        size[i] = 1;
        parent[i] = i;
        b[i+i*n] = 0;
        l[i+i*n] = 0;
    }

    for (k = 0; k < nd; ++k)
    {
        k2ij(n, dij_idx[k], &i, &j);

        if (mst_union(parent, size, i, j))
        {
            b[i+j*n] = b[j+i*n] = 1;
            l[i+j*n] = l[j+i*n] = -1;
            l[j+j*n] += 1;
            l[i+i*n] += 1;
            len += dij[k];
        }
        else
        {
            b[i+j*n] = b[j+i*n] = 0;
            l[i+j*n] = l[j+i*n] = 0;
        }
    }

    return len;
}
