#pragma once

#include <stdlib.h>

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

class DLLEXPORT CSorter
{
public:

  CSorter(int n, double *x, double *y, double *lp, double *lv, double *iv=NULL, double *jv=NULL);
  CSorter(int n, double *x, double *y, double cx, double cy, double *v = 0);
  CSorter(int n, double *t, double *k);
  CSorter(int n, double *x, double *y, double *v, double thresh);
  CSorter(double *t, double *x, double *y, double *d, int n);  // different order used to differentiate from 1st type
  
  ~CSorter();

protected:

  int num;
  int *sortindices;
  double *swapper;
};
