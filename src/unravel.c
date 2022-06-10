#include "prtree.h"

/*
** Given a flat index (k) into the lower triangle of a n-by-n matrix
** stored in column major order, return the corresponding row and 
** column indices (i,j). The diagonal is not included as part of the
** lower triangle. e.g.,
**
**   - - - - -    0 -> (1, 0)  5 -> (3, 1)
**   0 - - - -    1 -> (2, 0)  6 -> (4, 1)
**   1 4 - - -    2 -> (3, 0)  7 -> (3, 2)
**   2 5 7 - -    3 -> (4, 0)  8 -> (4, 2)
**   3 6 8 9 -    4 -> (2, 1)  9 -> (4, 3)
**
** Reference:
**
** Modified from Algorithm 3 of the following reference
**   https://hal.archives-ouvertes.fr/hal-02047514/document
**
** The linked reference uses a 1-based index and includes the diagonal. 
** Here we use a 0-based index and exclude the diagonal.
**
** See also: https://en.wikipedia.org/wiki/Triangular_number
*/
void k2ij(int n, int k, int *i, int *j)
{
    ++k;
    --n;
    int kp = ( n * (n + 1) ) / 2 - k;
    int p = ( sqrt(1 + 8 * kp) - 1 ) / 2;
    *i = k - n*(n-1)/2 + p*(p+1)/2;
    *j = n - p - 1;
}