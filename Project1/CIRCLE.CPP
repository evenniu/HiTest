#include "stdafx.h"
#include "circle.h"

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

// during reporting, it can be useful to retrieve information about the curve;
// I call these annotations
void CCircle::addDoubleAnnotation(const char* /*key*/, double /*value*/)
{
}
double CCircle::getDoubleAnnotation(const char* /*key*/) const
{
  return std::numeric_limits<double>::quiet_NaN();
}
void CCircle::addBoolAnnotation(const char* /*key*/, bool /*value*/)
{
}
bool CCircle::getBoolAnnotation(const char* /*key*/) const
{
  return false;
}

CCircle::CCircle(double *c, double *s, double *e, double r)
{
  m_isEnglish = true;  // NEED to set this
  int i;
  for (i=0; i<3; i++)
    center[i] = c[i];

  if (r <= 0.0)  // actually is a point
  {
    a0 = a1 = 0.0;
    rad = 0.0;
  }
  else
  {
    a0 = atan2(s[1] - c[1], s[0] - c[0]);
    a1 = atan2(e[1] - c[1], e[0] - c[0]);
    if (a1 < a0)
      a1 += 2*M_PI;

    if (a1 - a0 < 0.0001)
      a1 = a0 + 2*M_PI;

    rad = r;
  }

  storesize = sizeof(long) + 10L*sizeof(double);
  Min[0] = Min[1] = 1.0;
  Max[0] = Max[1] = 0.0;
  good = 1;
}

CCircle::CCircle(int n, double *x, double *y)  // calculate from set of points
{
  good = 0;
  m_isEnglish = true;  // NEED to set this

  if (!FitCircle(n, x, y))
    return;
  a0 = 0.0;
  a1 = 2*M_PI;
  storesize = sizeof(long) + 10L*sizeof(double);
  Min[0] = Min[1] = 1.0;
  Max[0] = Max[1] = 0.0;

  good = 1;
}

CCircle::~CCircle()
{}

int CCircle::Extreme(double *v, double *et, double *ep, double, double, int, bool)
{
  if (rad <= 0.0)  // should not be here
  {
    *et = T0();
    CalcPoint(ep, *et);
    return 1;
  }

  *et = atan2(v[1], v[0]);  // not likely here either (circle)
  CalcPoint(ep,*et);
  return 1;
}

double CCircle::ClosestPoint(double *tgt, double *bestxyz, double *bestt, double *tangent, double t0, double t1, int)
{
  if (t0 == 0.0 && t1 == 0.0)
  {
    t0 = T0();
    t1 = T1();
  }

  if (rad <= 0.0)
  {
    bestxyz[0] = X();
    bestxyz[1] = Y();

    if (bestt)
      *bestt = 0.0;

    if (tangent)
    {
      tangent[0] = 1.0;
      tangent[1] = 0.0;
    }

    return _hypot(bestxyz[0] - tgt[0], bestxyz[1] - tgt[1]);
  }

  double t = atan2(tgt[1] - Y(), tgt[0] - X());
  if (t < t0)
    t += 2*M_PI;
  if  (t >= t0 && t <= t1)
  {
    CalcPoint(bestxyz, t, tangent);
    if (bestt)
      *bestt = t;
    return _hypot(bestxyz[0] - tgt[0], bestxyz[1] - tgt[1]);
  }

  double d0,d1,b0[3],tn0[3],b1[3],tn1[3];

  CalcPoint(b0, t0, tn0);
  d0 = _hypot(b0[0] - tgt[0], b0[1] - tgt[1]);

  CalcPoint(b1, t1, tn1);
  d1 = _hypot(b1[0] - tgt[0], b1[1] - tgt[1]);

  if (d0 < d1)
  {
    if (bestt)
      *bestt = t0;
    if (tangent)
    {
      tangent[0] = tn0[0];
      tangent[1] = tn0[1];
    }
    return d0;
  }
  else
  {
    if (bestt)
      *bestt = t1;
    if (tangent)
    {
      tangent[0] = tn1[0];
      tangent[1] = tn1[1];
    }
    return d1;
  }
}

int CCircle::LineIntersect(double *xy, double *ij, double *sol, double t0,
                           double t1, double *tsol, double *dist, bool /*isRay*/, bool /*useSol*/)
{
  if (rad <= 0.0)
    return 0;

  if (t0 == 0.0 && t1 == 0.0)
  {
    t0 = T0();
    t1 = T1();
  }

  double cx = center[0];
  double cy = center[1];

  // calculate coefficients of quadratic (ax*x + bx + c = 0) that solves

  double a = ij[0]*ij[0] + ij[1]*ij[1];
  double b = 2.0*(ij[0]*(xy[0] - cx) + ij[1]*(xy[1] - cy));
  double c = (xy[0] - cx)*(xy[0] - cx) + (xy[1] - cy)*(xy[1] - cy) - rad*rad;

  double radical = b*b - 4*a*c;

  if (radical < 0.0)   // no solutions
    return 0;

  if (fabs(radical) < 1.0e-6)  // one solution
  {
    double t = (-b + sqrt(radical))/2/a;
    double dum[3];
    sol[0] = xy[0] + t*ij[0];
    sol[1] = xy[1] + t*ij[1];
    if (ClosestPoint(sol, dum, tsol, 0, t0, t1) > 1.0e-4) // not in domain
      return 0;
    if (dist)
      *dist = _hypot(xy[0] - sol[0], xy[1] - sol[1]);
    return 1;
  }
  else                              // two solutions
  {
    double bt0,bt1,d0,d1,s0[3],s1[3],dum[3];

#pragma warning(suppress : 4458)
    double a0 = (-b + sqrt(radical))/2/a;
    s0[0] = xy[0] + a0*ij[0];
    s0[1] = xy[1] + a0*ij[1];
    d0 = ClosestPoint(s0, dum, &bt0, 0, t0, t1);

#pragma warning(suppress : 4458)
    double a1 = (-b - sqrt(radical))/2/a;
    s1[0] = xy[0] + a1*ij[0];
    s1[1] = xy[1] + a1*ij[1];
    d1 = ClosestPoint(s1, dum, &bt1, 0, t0, t1);

    if (d0 > 1.0e-4 && d1 > 1.0e-4) // intersections not in domain
      return (0);

    if (d0 < d1)
    {
      sol[0] = s0[0];
      sol[1] = s0[1];
      if (tsol)
        *tsol = bt0;
    }
    else
    {
      sol[0] = s1[0];
      sol[1] = s1[1];
      if (tsol)
        *tsol = bt1;
    }

    if (dist)
      *dist = _hypot(xy[0] - sol[0], xy[1] - sol[1]);

    return 1;
  }
}

void CCircle::Extent()
{
  if (Min[0] <= Max[0] && Min[1] <= Max[1])
    return;
  Min[0] = center[0] - rad;
  Min[1] = center[1] - rad;

  Max[0] = center[0] + rad;
  Max[1] = center[1] + rad;
}

void CCircle::Extent(double *min, double *max)
{
  Extent();

  min[0] = Min[0];
  min[1] = Min[1];

  max[0] = Max[0];
  max[1] = Max[1];
}

int CCircle::CalcPoint(double *xyz, double a, double *tan, double *curv) const 
{
  xyz[0] = center[0] + rad*cos(a);
  xyz[1] = center[1] + rad*sin(a);
  //xyz[2] = center[2];
  if (m_palign)
    m_palign->MeasToBest(xyz, 1, xyz);

  if (tan)
  {
    if (rad > 0.0)
    {
      tan[0] = -rad*sin(a);
      tan[1] =  rad*cos(a);
    }
    else
    {
      tan[0] = 0.0;
      tan[1] = 0.0;
    }
    if (m_palign)
      m_palign->MeasToBest(tan, 0, tan);
  }

  if (curv)
  {
    if (rad > 0.0)
    {
      curv[0] = -rad*cos(a);
      curv[1] = -rad*sin(a);
    }
    else
    {
      curv[0] = 0.0;
      curv[1] = 0.0;
    }
    if (m_palign)
      m_palign->MeasToBest(curv, 0, curv);
  }

  return 1;
}

void CCircle::Write(FILE *fp)
{
  short type;
  if (rad > 0.0)
    type = CIRCLE_TYPE;
  else
    type = POINT_TYPE;
  fwrite(&type,sizeof(short),1,fp);
  fwrite(&m_extreme,sizeof(double),1,fp);

  fwrite(center,sizeof(double),3,fp);
  fwrite(&a0,sizeof(double),1,fp);
  fwrite(&a1,sizeof(double),1,fp);
  fwrite(&rad,sizeof(double),1,fp);
  fwrite(Min,sizeof(double),2,fp);
  fwrite(Max,sizeof(double),2,fp);
  fwrite(&seqno,sizeof(long),1,fp);
}

// FitCircle uses the Gauss-Newton method to estimate the center
// and radius of a circle

int CCircle::FitCircle(int npts, double *xp, double *yp)
{
  if (!CircInit(npts, xp, yp))
    return 0;

  int i;
  double *f, *h, *fak, **x, **J;

  int nparms = 3;

  if (npts < nparms)
    return 0;

  /* Allocate local variables. */
  x = new double *[npts];
  J = new double *[npts];
  for (i=0; i<npts; i++)
  {
    x[i] = new double[2];
    J[i] = new double[3];
  }
  f = new double[npts];
  h = new double[npts];
  fak = new double[npts];

  double EPS = 1.e-12;
  int MAXITS = 30;

  int iter;
  for (iter = 0; iter < MAXITS; iter++)
  {
    double temp, grad[4];

    for (i = 0; i < npts; i++)    // Compute x and f.
    {
      x[i][0] = xp[i] - center[0];
      x[i][1] = yp[i] - center[1];

      fak[i] = temp = norm(2,x[i]);
      f[i] = rad - temp;
    }

    for (i = 0; i < npts; i++)    // Compute J
    {
      double d;
      if ((d = fak[i]) != 0.0)
        temp = -1. / d;
      else
        temp = 0.;
      J[i][0] = x[i][0] * temp;
      J[i][1] = x[i][1] * temp;
      J[i][2] = -1.;
    }

    // Compute scaled gradient of objective function

    multATx(npts, nparms, J, f, grad);

    if (norm(nparms,grad) / npts < EPS)   // Check for convergence
      break;

    // Solve Jh = f in the least squares sense
    if (!LSnormal(npts, nparms, J, f, h))
    {
      for (i=0; i<npts; i++)
      {
        delete [] x[i];
        delete [] J[i];
      }

      delete [] x;
      delete [] J;
      delete [] f;
      delete [] h;
      delete [] fak;

      return 0;
    }

    center[0] += h[0]; // Update center
    center[1] += h[1];

    rad += h[2]; // Update rad
  }


  // Compute the least squares error. The residuals are in f.
  ls_err = norm(npts,f) / sqrt((double) npts);

  // Clean up

  for (i=0; i<npts; i++)
  {
    delete [] x[i];
    delete [] J[i];
  }

  delete [] x;
  delete [] J;
  delete [] f;
  delete [] h;
  delete [] fak;

  return 1;
}

/*---------------------------------------------------------------*/
/*  Get initial values for circle routine.                       */
/*  This only works well if the data are well distributed.       */
/*---------------------------------------------------------------*/

int CCircle::CircInit(int n, double *x, double *y)
{
  int i;

  // initialize center

  center[0] = center[1] = center[2] = 0.0;
  for (i=0; i<n; i++)
  {
    center[0] += x[i];
    center[1] += y[i];
  }
  center[0] /= n;
  center[1] /= n;

  // initialize radius
  for (rad = 0., i = 0; i < n; i++)
    rad += _hypot(x[i] - center[0], y[i] - center[1]);
  rad /= n;

  return 1;
}

double CCircle::norm(int n, double *dx)
{
  register int i;
  static double cutlo = 8.232e-11, cuthi = 1.304e19;
  double sum = 0., hitst, xmax;

  if (n < 1)
    return 0.;

  if (n == 1)
    return fabs(*dx);

  if (n == 2)
    return _hypot(dx[0],dx[1]);

  hitst = cuthi / (double) n;
  i = 0;

  // Zero Sum
  while (i < n)
  {
    if (*dx != 0.0)
      break;
    i++;
    dx++;
  }

  if (i >= n)
    return (0.0);

START:
  if (fabs(*dx) > cutlo)
  {
    for (; i < n; i++, dx++)   // Loop over elements
    {
      if (fabs(*dx) > hitst)
        goto GOT_LARGE;
      sum += (*dx) * (*dx);
    }
    sum = sqrt(sum);
    return (sum);              // Sum completed normally
  }
  else
  {                             // Small sum prepare for phase 2
    xmax = fabs(*dx);
    dx++;
    i++;
    sum += 1.0;
    for (; i < n; i++, dx++)
    {
      if (fabs(*dx) > cutlo) // Got normal element. Rescale and process
      {
        sum = (sum * xmax) * xmax;
        goto START;
      }
      if (fabs(*dx) > xmax)
      {
        sum = 1.0 + sum * (xmax / (*dx)) * (xmax / (*dx));
        xmax = fabs(*dx);
        continue;
      }
      sum += ((*dx) / xmax) * ((*dx) / xmax);
    }
    return (xmax * sqrt(sum));
  }

GOT_LARGE:
  sum = 1.0 + (sum / (*dx)) / (*dx);   // Rescale and process
  xmax = fabs(*dx);
  dx++;
  i++;
  for (; i < n; i++, dx++)
  {
    if (fabs(*dx) > xmax)
    {
      sum = 1.0 + sum * (xmax / (*dx)) * (xmax / (*dx));
      xmax = fabs(*dx);
      continue;
    }
    sum += ((*dx) / xmax) * ((*dx) / xmax);
  }
  return (xmax * sqrt(sum));    /* End of small sum. */

}

void CCircle::multATx(int rows, int cols, double **A, double *x, double *y)
{
  int i, j;
  double tmp;

  for (i = 0; i < cols; i++)
  {
    tmp = 0.;
    for (j = 0; j < rows; j++)
      tmp += A[j][i] * x[j];
    y[i] = tmp;
  }
}

/*--------------------------------------------------------------------*/
/*  Solve the least squares problem min ||Ax-b|| using the method of  */
/*  normal equations.                                                 */
/*  Note: A must have rank n.                                         */
/*--------------------------------------------------------------------*/

int CCircle::LSnormal(int m, int n, double **A, double *b, double *x)
{
  int i, j, k;
  double sum, **C;

  /* Compute d = A'b. */
  multATx(m,n,A,b,x);

  C = new double *[n];
  for (i=0; i<n; i++)
    C[i] = new double[n];

  /* Compute the lower triangle of A'A. */
  for (i = 0; i < n; i++)
  {
    for (j = 0; j <= i; j++)
    {
      sum = 0.;
      for (k = 0; k < m; k++)
        sum += A[k][i] * A[k][j];
      C[i][j] = sum;
    }
  }

  int rv = LDLT(n,C); // Compute the LDL' decomposition of C

  if (rv)
    rv = LDLTsolve(n, C, x); // Solve the system Cx = d

  for (i=0; i<n; i++)
    delete [] C[i];
  delete [] C;

  return rv;
}

//   Calculate the LDL' factorization A = L*D*L' of a symmetric, positive
//   definite n X n matrix.  L is unit lower triangular and D is diagonal
//   with positive diagonal elements.
// Input:
//    A   - real array of dimension n X n containing the matrix to be
//          factored.  It is assumed that A is symmetric and positive
//          definite.
//    n   - the size of the input matrix A.
// Output:
//    A   - real array of dimension n X n containing L in its lower
//          triangle and D on its diagonal.

int CCircle::LDLT(int n, double **A)
{
  int i, j, l;
  double vj, sum, *v;

  v = new double[n];

  for (j = 0; j < n; j++)
  {
    for (i = 0; i < j; i++)
      v[i] = A[j][i] * A[i][i];
    sum = A[j][j];
    for (i = 0; i < j; i++)
      sum -= A[j][i] * v[i];
    A[j][j] = v[j] = sum;
    for (i = j+1; i < n; i++)
    {
      sum = 0.;
      for (l = 0; l < j; l++)
        sum += A[i][l] * v[l];
      if ((vj = v[j]) == 0.)
      {
        delete [] v;
        return 0;
      }
      A[i][j] = (A[i][j] - sum) / vj;
    }
  }

  delete [] v;
  return 1;
}

//  Solve the system Ax = LDL'x = b.  Overwrite b with x.
//  Note that L and D are stored in A:  D overwrites the diagonal of
//  A and L overwrites the lower triangle of A (except the diagonal).
//  The diagonal of L is all ones.

int CCircle::LDLTsolve(int n, double **A, double *b)
{
  int i;

  lower_backsub(n, A, b);
  for (i = 0; i < n; i++)
  {
    if (fabs(A[i][i]) < 1.0e-15)
      return 0;

    b[i] /= A[i][i];
  }

  upper_backsub(n, A, b);

  return 1;
}

//  Solve the lower triangular system Lz = b.  Overwrite b with z.
//  Note that L[i][i] = 1 for all i.

void CCircle::lower_backsub(int n, double **L, double *b)
{
  int i, j;
  double sum;

  for (i = 0; i < n; i++)
  {
    for (sum = 0., j = 0; j < i; j++)
      sum += L[i][j] * b[j];
    b[i] -= sum;
  }
}

//  Solve the upper triangular system L'x = y.  Overwrite y with x.
//  Note that L[i][i] = 1 for all i.

void CCircle::upper_backsub(int n, double **L, double *y)
{
  int i, j;
  double sum;

  for (i = n-1; i >= 0; i--)
  {
    sum = 0.;
    for (j = i+1; j < n; j++)
      sum += L[j][i] * y[j];
    y[i] -= sum;
  }
}
