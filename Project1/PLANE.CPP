#include "StdAfx.h"
#include "svd.h"
#include "plane.h"

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

int calc_plane(int n,double **x, double *pt, double *vc)
{
  int i;
  double s[3],e[3];

  if (n <= 2)
    return 0;

  CVector work(n);
  CMatrix A(n, 3);

  pt[0] = pt[1] = pt[2] = 0.0;  // centroid
  for (i=0; i<n; i++)
  {
    pt[0] += x[i][0];
    pt[1] += x[i][1];
    pt[2] += x[i][2];
  }
  pt[0] /= n;
  pt[1] /= n;
  pt[2] /= n;
  for (i=0; i<n; i++)
  {
    A.m[i][0] = x[i][0] - pt[0];
    A.m[i][1] = x[i][1] - pt[1];
    A.m[i][2] = x[i][2] - pt[2];
  }

  /* singular value decomposition (values sorted in decreasing order) */
  dsvdc(A.m, n, 3, s, e, NULL, A.m, work.v, 1, 1);

  for (i = 0; i < 3; i++)   // normal is singular vector
    vc[i] = A.m[i][2];      // of smallest singular value

  // make vector point "up"

  if (vc[2] < 0.0)
  {
    vc[0] *= -1.0;
    vc[1] *= -1.0;
    vc[2] *= -1.0;
  }

  return 1;
}
