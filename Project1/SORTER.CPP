#include "StdAfx.h"
#include "Sorter.h"

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

static double gThresh;
static double *sortxp, *sortyp, sortp[2], sortv[2];

static int compfuncs(const void *a, const void *b);
static int compfunce(const void *a, const void *b);
static int compfunct(const void *a, const void *b);
static int compfuncf(const void *a, const void *b);

// sort based on projection of points onto a line

CSorter::CSorter(int n, double *x, double *y, double *lp, double *lv, double *iv, double *jv)
{
  sortp[0] = lp[0];
  sortp[1] = lp[1];
  sortv[0] = lv[0];
  sortv[1] = lv[1];

  num = n;
  if (num < 2)
    return;

  sortindices = new int[num];
  swapper = new double[num];
  sortxp = x;
  sortyp = y;

  int i;
  for (i=0; i<num; i++)
    sortindices[i] = i;

  qsort((char *) sortindices, num, sizeof(int), compfuncs);

  for (i=0; i<num; i++)
    swapper[i] = x[sortindices[i]];
  for (i=0; i<num; i++)
    x[i] = swapper[i];

  for (i=0; i<num; i++)
    swapper[i] = y[sortindices[i]];
  for (i=0; i<num; i++)
    y[i] = swapper[i];
    
  if (iv)
  {
    for (i=0; i<num; i++)
      swapper[i] = iv[sortindices[i]];
    for (i=0; i<num; i++)
      iv[i] = swapper[i];
  }
  
  if (jv)
  {
    for (i=0; i<num; i++)
      swapper[i] = jv[sortindices[i]];
    for (i=0; i<num; i++)
      jv[i] = swapper[i];
  }
}

// sort based on radial position about a centroid

CSorter::CSorter(int n, double *x, double *y, double cx, double cy, double *v)
{
  sortp[0] = cx;
  sortp[1] = cy;

  num = n;
  if (num < 2)
    return;

  sortindices = new int[num];
  swapper = new double[num];
  sortxp = x;
  sortyp = y;

  int i;
  for (i=0; i<num; i++)
    sortindices[i] = i;

  qsort((char *) sortindices, num, sizeof(int), compfunce);

  for (i=0; i<num; i++)
    swapper[i] = x[sortindices[i]];
  for (i=0; i<num; i++)
    x[i] = swapper[i];

  for (i=0; i<num; i++)
    swapper[i] = y[sortindices[i]];
  for (i=0; i<num; i++)
    y[i] = swapper[i];

  // find largest opening, rotate arrays so start and end points at gap

  double maxgap = 0.0;
  int g=0;
  for (i=0; i<num; i++)
  {
    int j = i-1;
    if (i == 0)
      j = num - 1;
    double gap = _hypot(x[i] - x[j], y[i] - y[j]);
    if (i == 0 || gap > maxgap)
    {
      g = i;
      maxgap = gap;
    }
  }

  if (v)
  {
    for (i=0; i<num; i++)
    {
      int j = i-1;
      if (i == 0)
        j = num - 1;
      double vi[2], vj[2];
      vi[0] = x[i] - cx;
      vi[1] = y[i] - cy;
      vj[0] = x[j] - cx;
      vj[1] = y[j] - cy;
      if ((dot(vi,v) > 0.0) && (dot(vj,v) > 0.0) && (curl(vj, v) > 0.0) && (curl(v, vi) > 0.0))
      {
        //bugout(0, _T("9 %f %f %f %f"), cx, cy, x[i], y[i]);
        //bugout(0, _T("9 %f %f %f %f"), cx, cy, x[j], y[j]);
        g = i;
        break;
      }
    }
  }

  if (g > 0)  // shift the points
  {
    for (i=0; i<num; i++)
      swapper[i] = x[(i+g)%num];
    for (i=0; i<num; i++)
      x[i] = swapper[i];

    for (i=0; i<num; i++)
      swapper[i] = y[(i+g)%num];
    for (i=0; i<num; i++)
      y[i] = swapper[i];
  }
}

// sort based on t value of points projected onto nominal curves, k values are sorted

CSorter::CSorter(int n, double *t, double *k)
{
  num = n;
  if (num < 2)
    return;

  sortindices = new int[num];
  swapper = new double[num];
  sortxp = t;
  sortyp = k;

  int i;
  for (i=0; i<num; i++)
    sortindices[i] = i;

  qsort((char *) sortindices, num, sizeof(int), compfunct);

  for (i=0; i<num; i++)
    swapper[i] = t[sortindices[i]];
  for (i=0; i<num; i++)
    t[i] = swapper[i];

  for (i=0; i<num; i++)
    swapper[i] = k[sortindices[i]];
  for (i=0; i<num; i++)
    k[i] = swapper[i];
}

// sort designed to separate two branches that begin at the first node

CSorter::CSorter(int n, double *x, double *y, double *v, double thresh)
{
  sortp[0] = x[0];
  sortp[1] = y[0];
  sortv[0] = v[0];
  sortv[1] = v[1];
  
  gThresh = thresh;

  num = n;
  if (num < 2)
    return;

  sortindices = new int[num];
  swapper = new double[num];
  sortxp = x;
  sortyp = y;

  int i;
  for (i=0; i<num; i++)
    sortindices[i] = i;

  qsort((char *) sortindices, num, sizeof(int), compfuncf);

  for (i=0; i<num; i++)
    swapper[i] = x[sortindices[i]];
  for (i=0; i<num; i++)
    x[i] = swapper[i];

  for (i=0; i<num; i++)
    swapper[i] = y[sortindices[i]];
  for (i=0; i<num; i++)
    y[i] = swapper[i];
}

// sort based on t value of points projected onto nominal curves, x and y values are sorted

CSorter::CSorter(double *t, double *x, double *y, double *d, int n)
{
  num = n;
  if (num < 2)
    return;

  sortindices = new int[num];
  swapper = new double[num];
  sortxp = t;

  int i;
  for (i=0; i<num; i++)
    sortindices[i] = i;

  qsort((char *) sortindices, num, sizeof(int), compfunct);

  for (i=0; i<num; i++)
    swapper[i] = t[sortindices[i]];
  for (i=0; i<num; i++)
    t[i] = swapper[i];

  for (i=0; i<num; i++)
    swapper[i] = x[sortindices[i]];
  for (i=0; i<num; i++)
    x[i] = swapper[i];

  for (i=0; i<num; i++)
    swapper[i] = y[sortindices[i]];
  for (i=0; i<num; i++)
    y[i] = swapper[i];

  for (i=0; i<num; i++)
    swapper[i] = d[sortindices[i]];
  for (i=0; i<num; i++)
    d[i] = swapper[i];
}

CSorter::~CSorter()
{
  if (sortindices)
    delete [] sortindices;

  if (swapper)
    delete [] swapper;
}


// compare function for measured curves based on t points on nominal curves.

int compfunct(const void *a, const void *b)
{
  double da = sortxp[*(int *)a];  // sortxp should be nomt
  double db = sortxp[*(int *)b];
  if (da > db)
    return 1;
  if (da < db)
    return -1;
  return 0;
}

// compare function for CV and CC curves and mean camber curves.
// sorts based on projections onto a line
// sortp is a point on the line.
// sortv is the line vector

int compfuncs(const void *a, const void *b)
{
  double pt[2];
  pt[0] = sortxp[*(int *)a];
  pt[1] = sortyp[*(int *)a];
  double d1 = projdist(sortp, sortv, pt);
  pt[0] = sortxp[*(int *)b];
  pt[1] = sortyp[*(int *)b];
  double d2 = projdist(sortp, sortv, pt);

  if (d1 < d2)
    return -1;
  if (d1 > d2)
    return 1;
  return 0;
}

// compare function for LE and TE curves.
// sortp is approximately the centroid of the nominal curve.

int compfunce(const void *a, const void *b)
{
  double x1 = sortxp[*(int *)a];
  double y1 = sortyp[*(int *)a];
  double x2 = sortxp[*(int *)b];
  double y2 = sortyp[*(int *)b];
  double a1 = atan2(y1 - sortp[1], x1 - sortp[0]);
  double a2 = atan2(y2 - sortp[1], x2 - sortp[0]);

  if (a1 < a2)
    return -1;
  if (a1 > a2)
    return 1;
  return 0;
}

// compare function for separating forks as part of voronoi mean camber point calculations
// sorts based on projections onto a line
// sortp is a point on the line.
// sortv is the line vector

int compfuncf(const void *a, const void *b)
{
  double vec[2];
  
  vec[0] = sortxp[*(int *)a] - sortp[0];
  vec[1] = sortyp[*(int *)a] - sortp[1];
  double d1 = normalize(vec, vec);
  if (d1 > 0.0)
    if (dot(vec, sortv) < gThresh)  // not on main line
      d1 += 1.0e6;

  vec[0] = sortxp[*(int *)b] - sortp[0];
  vec[1] = sortyp[*(int *)b] - sortp[1];
  double d2 = normalize(vec, vec);
  if (d2 > 0.0)
    if (dot(vec, sortv) < gThresh)  // not on main line
      d2 += 1.0e6;

  if (d1 < d2)
    return -1;
  if (d1 > d2)
    return 1;
  return 0;
}
