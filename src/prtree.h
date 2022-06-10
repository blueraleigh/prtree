#ifndef PRTREE_H
#define PRTREE_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>

double mst(int,const double*,const int*,int*,double*,int*,int*);
void k2ij(int,int,int*,int*);
double *dijkstra(int, int, const int*, const double*);

int ndist(int);
double dist2(int,const double*,const double*);
void dist(int,int,const double*,double*,int*);

#endif