#include "StdAfx.h"
#include <stdio.h>
#include <memory.h>
#include <string.h>
#include <time.h>
#include "CIRCLE.H"
#include <assert.h>
#include "SUBCURVE.H"



#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

#define minDebugLevel 0

#ifndef _DEBUG
void bugout(int, const WCHAR*, ...)
{
#else
void bugout(int level, const WCHAR* format, ...)
{
  if(level > minDebugLevel)
    return;

  va_list args;

  WCHAR out_buf[300];
  va_start(args, format);
  //_vstprintf(out_buf, format, args);
  vswprintf_s(out_buf, 300, format, args);
  va_end(args);

  FILE* fout;
  _wfopen_s(&fout, L"C:/Users/Public/Documents/Hitest.out", L"at");
  if(fout)
  {
    fwprintf(fout, L"%s\n", out_buf);
    fclose(fout);
  }
#endif
}

void ang_to_dms(double ang, wchar_t *buf)
{
  wchar_t sign('+');
  if (ang < 0)
    sign = wchar_t('-');
  ang = fabs(ang);
  int deg = (int)ang;
  ang = ang - deg;
  ang *= 60;
  int min = (int)ang;
  ang = ang - min;
  ang *= 60;
  int sec = (int)(ang + 0.5);
  if (sec == 60)
  {
    sec = 0;
    min ++;
    if (min == 60)
    {
      min = 0;
      deg++;
    }
  }
  
  swprintf_s(buf, 40, L"%c%d:%02d:%02d",sign,deg,min,sec);
}



// This method returns the angle between two vectors
// The returned angle is in the range [0, pi]
double AngleBetweenVectors2d(const double * vec1, const double* vec2)
{
  return atan2(
    norm_cross2d(vec1, vec2),
    dot(vec1, vec2)
    );
}
double AngleBetweenVectors3d(const double * vec1, const double* vec2)
{
  return atan2(
    norm_cross3d(vec1, vec2),
    dot(vec1, vec2, 3)
    );
}



/*--------------------------------------------------------*/
/*  DDOT  Form the dot product of two vectors.            */
/*  INPUT                                                 */
/*    n     Number of elements to sum.                    */
/*    x     Address of first element of x vector.         */
/*    y     Address of first element of y vector.         */
/*  OUPUT                                                 */
/*    ddot  Dot product x and y.                          */
/*--------------------------------------------------------*/

double ddot(int n, double *x, double *y)
{
  int i;
  double dot = 0.;

  if (n < 1)
    return 0.;
  for (i = 0; i < n; i++)
    dot += x[i] * y[i];

  return dot;

}  /* end ddot */


/*
* Return the current knot the parameter u is less than or equal to.
* Find this "breakpoint" allows the evaluation routines to concentrate on
* only those control points actually effecting the curve around u.]
*
*      m   is the index of last control point (i.e., numpts - 1)
*      k   is the order of the curve (or surface direction)
*      kv  is the knot vector ([0..m+k-1]) to find the break point in.
*/

int FindBreakPoint(double u, double * kv, int m, int k)
{
  int l = k-1;
  int r = m+k;
  int half = (r - l)/2;
  int in = l + half;

  while (half)
  {
    if (u < kv[in])
      r = in;
    else
      l = in;
    half = (r - l)/2;
    in = l + half;
  }

  return in;
}

/*---------------------------------------------------------------*/
/*  NORM  Compute the Euclidean norm of dx while being very      */
/*        careful of destructive underflow and overflow.         */
/*  INPUT                                                        */
/*    n       Number of elements to use.                         */
/*    dx      Address of first element of x vector.              */
/*  OUTPUT                                                       */
/*    Euclidean norm of dx.                                      */
/*---------------------------------------------------------------*/

double l2norm(const double *dx, int n)
{
  register int i;
  static double cutlo = 8.232e-11, cuthi = 1.304e19;
  double sum = 0., hitst, xmax;

  if (n <= 2)
  {
    if (n == 2)
    {
      double xv = fabs(dx[0]);
      double yv = fabs(dx[1]);
      if (xv > yv)
      {
        double h = yv/xv;
        return xv*sqrt(1.0 + h*h);
      }
      else if (yv)
      {
        double h = xv/yv;
        return yv*sqrt(1.0 + h*h);
      }
      else
        return 0;
    }
    else if (n == 1)
      return fabs(*dx);
    else
      return 0.;
  }

  hitst = cuthi / (double) n;
  i = 0;

  /* Zero Sum. */
  while (i < n)
  {
    if (*dx != 0.0)
      break;
    i++;
    dx++;
  }
  if (i >= n)
    return (sum);

START:
  if (fabs(*dx) > cutlo)
  {
    for (; i < n; i++, dx++)   /* Loop over elements. */
    {
      if (fabs(*dx) > hitst)
        goto GOT_LARGE;
      sum += (*dx) * (*dx);
    }
    sum = sqrt(sum);
    return (sum);      /* Sum completed normally. */
  }
  else   /* Small sum prepare for phase 2. */
  {
    xmax = fabs(*dx);
    dx++;
    i++;
    sum += 1.0;
    for (; i < n; i++, dx++)
    {
      if (fabs(*dx) > cutlo)   /* Got normal element. Rescale and process. */
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
  sum = 1.0 + (sum / (*dx)) / (*dx);   /* Rescale and process. */
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

double sqnorm(double *dx, int n)
{
  register int i;
  static double cutlo = 8.232e-11, cuthi = 1.304e19;
  double sum = 0.0, hitst, xmax, tmp, temp;

  if (n == 0)
    return 0.;

  hitst = cuthi / (double) n;
  i = 0;

  /* Zero sum */
  while (i < n && *dx == 0.)
  {
    ++dx;
    ++i;
  }
  if (i >= n)
    return (sum);

  if ((tmp = fabs(*dx)) > cutlo)
  {
START:
    for (; i < n; ++dx, ++i)
    {               /* Loop over elements. */
      if ((tmp = fabs(*dx)) > hitst)
      {
        sum = 1.0 + (sum / tmp) / tmp;  /* Rescale and process. */
        xmax = tmp;
        ++dx;
        ++i;
        for (; i < n; ++dx, ++i)
        {
          if ((tmp = fabs(*dx)) > xmax)
          {
            sum = 1.0 + sum * (temp = xmax / tmp, temp * temp);
            xmax = tmp;
            continue;
          }
          sum += (temp = tmp / xmax, temp * temp);
        }
        return ((sum * xmax) * xmax);
      }
      sum += tmp * tmp;
    }
    return (sum);   /* Sum completed normally. */
  }
  else
  {                       /* Small sum prepare for phase 2. */
    xmax = tmp;
    ++dx;
    ++i;
    sum += 1.0;
    for (; i < n; ++dx, ++i)
    {
      if (fabs(*dx) > cutlo)
      {
        /* Got normal element. Rescale and process. */
        sum = (sum * xmax) * xmax;
        goto START;
      }
      if ((tmp = fabs(*dx)) > xmax)
      {
        sum = 1.0 + sum * (temp = xmax / tmp, temp * temp);
        xmax = tmp;
        continue;
      }
      sum += (temp = tmp / xmax, temp * temp);
    }
    return ((sum * xmax) * xmax);
  }
}

double dot(const double *u, const double *v, int ndim)
{
  double rv = 0.0;

  int i;
  for (i=0; i<ndim; i++)
    rv += u[i]*v[i];

  return rv;
}


// Calculate the norm of the cross product
double norm_cross2d(const double v1[2], const double v2[2])
{
  return abs(v1[0] * v2[1] - v1[1] * v2[0]);
}
double norm_cross3d(const double v1[3], const double v2[3])
{
  double crossproduct[3];
  cross(v1, v2, crossproduct);
  return l2norm(crossproduct, 3);
}


// Calculate cross product, vo = v1 X v2.  For 3d vectors only
void cross(const double v1[3], const double v2[3], double vo[3])
{
  double vt[3];

  vt[0] = v1[1] * v2[2] - v1[2] * v2[1];
  vt[1] = v1[2] * v2[0] - v1[0] * v2[2];
  vo[2] = v1[0] * v2[1] - v1[1] * v2[0];
  vo[0] = vt[0];
  vo[1] = vt[1];
}

void cross_product(double *U, double *V, double *N)
{
  N[0] = U[1] * V[2] - V[1] * U[2];
  N[1] = -U[0] * V[2] + V[0] * U[2];
  N[2] = U[0] * V[1] - V[0] * U[1];
}


double clip(double value, double low, double high)
{
  assert(low < high);
  if(value < low) return low;
  if(value > high) return high;
  return value;
}


double normalize(double *vo, const double *vi, int ndim)
{
  int i;
  double temp = l2norm(vi, ndim);

  if (temp > 1.0e-15)
  {
    for (i=0; i<ndim; i++)
      vo[i] = vi[i] / temp;
  }
  else
  {
    for (i=0; i<ndim; i++)
      vo[i] = vi[i];
  }

  return(temp);
}

//  SGIVENS:  Construct a Givens plane rotation.

void SGivens(double *da, double *db, double *c, double *s)
{
  double r, scale, z, rho, d1, d2;

  rho = *db;
  if (fabs(*da) > fabs(*db))
    rho = *da;
  scale = fabs(*da) + fabs(*db);
  if (scale == 0.)
  {
    *c = 1.;
    *s = 0.;
    r = 0.;
  }
  else
  {
    d1 = *da / scale;
    d2 = *db / scale;
    r = scale * sqrt(d1 * d1 + d2 * d2);
    r = (rho < 0. ? -1. : 1.) * r;
    *c = *da / r;
    *s = *db / r;
  }
  z = 1.;
  if (fabs(*da) > fabs(*db))
    z = *s;
  if (fabs(*db) >= fabs(*da) && *c != 0.)
    z = 1. / *c;
  *da = r;
  *db = z;
}

double curl(double *v1,double *v2)
{
  return v1[0] * v2[1] - v1[1] * v2[0];
}

// intersect two lines


int int_line_line(double x1, double y1, double i1, double j1, double x2,
                  double y2, double i2, double j2, double *sx, double *sy)
{
  double den = i2*j1 - i1*j2;
  if (fabs(den) < 1.0e-12)
    return 0;

  double num = i2*(y2-y1) - j2*(x2-x1);

  double alpha = num/den;

  *sx = x1 + alpha*i1;
  *sy = y1 + alpha*j1;

  return 1;
}

int int_line_line(double *p1, double *v1, double *p2, double *v2, double *sol)
{
  return int_line_line(p1[0], p1[1], v1[0], v1[1], p2[0], p2[1], v2[0], v2[1], &sol[0], &sol[1]);
}

// fit parabola to points (x0,y0), (x1,y1) and (x2,y2).  f(x) = aX^2 + bX + c

int parabola(double *x, double *y, double *a, double *b, double *c)
{
  double den;

  den = x[0]*x[0]*x[1] + x[0]*x[2]*x[2] + x[1]*x[1]*x[2] -
    x[1]*x[2]*x[2] - x[0]*x[1]*x[1] - x[0]*x[0]*x[2];

  if (fabs(den) < 1.0e-8)  // should only happen if x is duplicated
    return 0;

  *a = y[0]*x[1] + x[0]*y[2] + y[1]*x[2] - x[1]*y[2] - x[2]*y[0] - x[0]*y[1];
  *a /= den;

  *b = x[0]*x[0]*y[1] + y[0]*x[2]*x[2] + x[1]*x[1]*y[2] -
    y[1]*x[2]*x[2] - y[2]*x[0]*x[0] - x[1]*x[1]*y[0];
  *b /= den;

  *c = x[0]*x[0]*x[1]*y[2] + x[0]*y[1]*x[2]*x[2] + y[0]*x[1]*x[1]*x[2] -
    y[0]*x[1]*x[2]*x[2] - y[1]*x[2]*x[0]*x[0] - y[2]*x[1]*x[1]*x[0];
  *c /= den;

  return 1;
}

// intersect parabola and line

int int_para_line(double a, double b, double c,           // parabola definition
                  double x, double y, double i, double j, // line definition
                  double *rx, double *ry)                 // return values
{
  double aa,bb,cc,root;

  if (fabs(y - a*x*x - b*x - c) < 1.0e-8)  // (x,y) on parabola
  {
    *rx = x;
    *ry = y;
  }

  aa = a*i*i;
  bb = 2*a*x*i + b*i - j;
  cc = a*x*x + b*x + c - y;
  root = bb*bb - 4*aa*cc;

  if (fabs(aa) < 1.0e-12)  // just a line
  {
    if (fabs(bb) < 1.0e-12)  // but line is parallel to x axis
      return 0; // no intersections

    double ta = -cc/bb;
    *rx = x+ta*i;
    *ry = y+ta*j;
    return 1;
  }

  if (root < 0.0)
    return 0; // no intersections

  double a1, a2, x1, y1, x2, y2;

  a1 = (-bb + sqrt(root))/2.0/aa;
  x1 = x+a1*i;
  y1 = y+a1*j;

  a2 = (-bb - sqrt(root))/2.0/aa;
  x2 = x+a2*i;
  y2 = y+a2*j;

  if (_hypot(x1-x, y1-y) < _hypot(x2-x, y2-y))  // choose closest to (x,y)
  {
    *rx = x1;
    *ry = y1;
  }
  else
  {
    *rx = x2;
    *ry = y2;
  }

  return 1;
}

// intersect parabola and circle centered at (0,0) of radius r closest to (x,y)

int int_para_circ(double a, double b, double c,   // parabola definition
                  double x, double y, double r,   // seed point and radius
                  double *rx, double *ry)         // return values
{
  // Newton-Raphson to find zeroes of f(x) = ax^2 + bx + c +/- sqrt(r - x^2)
  // Use the zero where rx is closest to input x, which should be a decent
  // guess.

  double x1 = x, y1, dx, dy;

  int iter;
  for (iter=0; iter<10; iter++)   // check top half of circle
  {
    if (x1 < -r + 1.0e-6)           // not a good candidate
    {
      x1 = -r;
      break;
    }

    if (x1 > r - 1.0e-6)
    {
      x1 = r;
      break;
    }

    y1 = a*x1*x1 + b*x1 + c - sqrt(r*r - x1*x1);
    dy = 2*a*x1 + b + x1/sqrt(r*r - x1*x1);

    if (fabs(dy) < 1.0e-6)  // function at an extrema, bump it off
      dy = 1.0e-6;

    dx = -y1/dy;

    if (fabs(dx) < 1.0e-6)
      break;

    x1 += dx;
  }
  y1 = a*x1*x1 + b*x1 + c;

  double x2 = x, y2;

  for (iter=0; iter<10; iter++)   // check bottom half of circle
  {
    if (x2 < -r + 1.0e-6)           // not a good candidate
    {
      x2 = -r;
      break;
    }

    if (x2 > r - 1.0e-6)
    {
      x2 = r;
      break;
    }

    y2 = a*x2*x2 + b*x2 + c + sqrt(r*r - x2*x2);
    dy = 2*a*x2 + b - x2/sqrt(r*r - x2*x2);

    if (fabs(dy) < 1.0e-6)
      dy = 1.0e-6;

    dx = -y2/dy;

    if (fabs(dx) < 1.0e-6)
      break;

    x2 += dx;
  }
  y2 = a*x2*x2 + b*x2 + c;

  if (_hypot(x1 - x,y1 - y) < _hypot(x2 - x,y2 - y))
  {
    *rx = x1;
    *ry = y1;
  }
  else
  {
    *rx = x2;
    *ry = y2;
  }
  return 1;
}

// construct vector normal to parabola at a given x, point toward input (ni,nj)

void para_norm(double a, double b, double /*c*/, double x,
               double *ni, double *nj)
{
  double vi = 2*a*x + b;
  double vj = -1.0;

  double den = _hypot(vi,vj);  // length for normalizing vector

#ifdef CARE_ABOUT_DIRECTION   // don't so far

  if (vi*(*ni) + vj*(*nj))
    < 0.0)
  {
    *ni = -vi/den;
    *nj = -vj/den;
  }
  else
  {
    *ni = vi/den;
    *nj = vj/den;
  }
#else
  *ni = vi/den;
  *nj = vj/den;
#endif
}

void rotate(int n, double theta, double *xi, double *yi, double *xo, double *yo)
{
  double s = sin(theta);
  double c = cos(theta);
  int i;
  for (i = 0; i < n; i++)
  {
    double temp = xi[i]*c - yi[i]*s;   // allow rotating into same arrays
    yo[i] = xi[i]*s + yi[i]*c;
    xo[i] = temp;
  }
}


/*--------------------------------------------------------------*/
/*  Use Kahan's summation formula to compute weighted `. */
/*--------------------------------------------------------------*/

void centroid(int n, int ndim, double **X, double *X0, double *w)
{
  int i, j;
  long double c;
  double s, t, wbar;

  if (w)
  {
    double cw, tw;
    wbar = w[0];
    cw = 0.0;
    for (i = 1; i < n; i++)
    {
      cw += w[i];
      tw = wbar + cw;
      cw -= (tw - wbar);
      wbar = tw;
    }
  }
  else
    wbar = (double)n;

  for (j = 0; j < ndim; j++)
  {
    s = w ? w[0]*X[0][j] : X[0][j];
    c = 0.0;
    for (i = 1; i < n; i++)
    {
      c += (w ? w[i]*X[i][j] : X[i][j]);
      t = s + (double)c;
      c -= (t - s);
      s = t;
    }
    X0[j] = s/wbar;
  }
}

// about the same as above...

double wcentroid(int n, int ndim, double **X, double *w, double *X0)
{
  int i, j;
  long double c, cw;
  double s, t, tw, wbar;

  wbar = w[0];
  cw = 0.;
  for (i = 1; i < n; i++)
  {
    cw += w[i];
    tw = wbar + (double)cw;
    cw -= (tw - wbar);
    wbar = tw;
  }
  for (j = 0; j < ndim; j++)
  {
    s = w[0] * X[0][j];
    c = 0.;
    for (i = 1; i < n; i++)
    {
      c += w[i] * X[i][j];
      t = s + (double)c;
      c -= (t - s);
      s = t;
    }
    X0[j] = s / wbar;
  }
  return wbar;

}  /* end wcentroid */



/*----------------------------------------------------------------------*/
/*  Translate a set of points by a given vector u.                      */
/*  P is an npoints X ndim matrix containing the points.                */
/*----------------------------------------------------------------------*/

void translate_all(int npoints, int ndim, double *u, double **P)
{
  int i, j;
  double uj;

  for (j = 0; j < ndim; j++)
  {
    if ((uj = u[j]) != 0.)
    {
      for (i = 0; i < npoints; i++)
        P[i][j] += uj;
    }
  }
}


/*----------------------------------------------------------------------*/
/*  Same as translate_all except it subtracts instead of adding.        */
/*----------------------------------------------------------------------*/

void shift_all(int npoints, int ndim, double *u, double **P)
{
  int i, j;
  double uj;

  for (j = 0; j < ndim; j++)
  {
    if ((uj = u[j]) != 0.)
    {
      for (i = 0; i < npoints; i++)
        P[i][j] -= uj;
    }
  }
}

/*-------------------------------------------------------------------*/
/*  Compute the signed area of the parallelogram formed by the three */
/*  points (x[k],y[k]), (x[j],y[j]), and (x[i],y[i]).                */
/*-------------------------------------------------------------------*/

double sarea(double *x, double *y, int k, int j, int i)
{
  return (x[k] - x[i]) * (y[j] - y[i]) - (y[k] - y[i]) * (x[j] - x[i]);
}

/*-------------------------------------------------------------------*/
/*  Compute the vector area of the parallelogram formed by the three */
/*  points (x[k],y[k],z[k]), (x[j],y[j],z[j]), and (x[i],y[i],z[i]). */
/*-------------------------------------------------------------------*/

double varea(double *x, double *y, double *z, int k, int j, int i, double *N)
{
  double xki, xji, yki, yji, zki, zji;

  xki = x[k] - x[i];
  xji = x[j] - x[i];
  yki = y[k] - y[i];
  yji = y[j] - y[i];
  zki = z[k] - z[i];
  zji = z[j] - z[i];
  N[0] = yki * zji - zki * yji;
  N[1] = -xki * zji + zki * xji;
  N[2] = xki * yji - yki * xji;

  return l2norm(N,3);
}

/*------------------------------------------------------------*/
/*  Solve the tridiagonal system LUx = r.                     */
/*------------------------------------------------------------*/

void solvetrid(int n, double *a, double *d, double *b, double *r, double *x)
{
  double *y = new double[n];

  /* Solve Ly = r. */

  y[0] = r[0];
  int i;
  for (i = 1; i < n; i++)
    y[i] = r[i] - b[i-1] * y[i-1];

  /* Solve Ux = y. */
  x[n-1] = y[n-1] / d[n-1];
  for (i = n-2; i >= 0; i--)
    x[i] = (y[i] - a[i] * x[i+1]) / d[i];

  delete [] y;
}

/*------------------------------------------------------------*/
/*  Compute the LU decomposition of a tridiagonal matrix.     */
/*  a is superdiagonal, d is diagonal, and b is subdiagonal.  */
/*------------------------------------------------------------*/

void tridlu(int n, double *a, double *d, double *b)
{
  int k;

  for (k = 0; k < n-1; k++)
  {
    b[k] /= d[k];
    d[k+1] -= b[k] * a[k];
  }
}

/*-------------------------------------------------------------------*/
/*  Find the one positive zero of the cubic polynomial               */
/*  p(t) = at^3 - bt^2 + bt - c, a > 0, c > 0.                       */
/*-------------------------------------------------------------------*/

double solve_cubic(double a, double b, double c, int *flag)
{
  double u, q, r, s1, s2, t;
  static double third = 1. / 3.;

  *flag = 0;
  if (a <= 0.)
    return 0.;
  u =  -(b / a) * third;
  q = u * (u + 1.);
  r = u * u * (u + 1.5) - .5 * c / a;
  t = r * r - q * q * q;
  if (t < 0.)
  {
    if (t < -1.e-2)
    {
      *flag = -1;
      return 0.;
    }
    t = 0.;
  }
  s1 = -SGN(r) * pow(fabs(r) + sqrt(t), third);
  if (s1 == 0.)
    s2 = 0.;
  else
    s2 = q / s1;

  return s1 + s2 + u;
}

/*---------------------------------------------------------------*/
/*  Compute Euclidean distance between two points while          */
/*  being very careful of destructive underflow and overflow.    */
/*---------------------------------------------------------------*/

double dist(double *dx, double *dy, int n)
{
  register int i;
  static double cutlo = 8.232e-11, cuthi = 1.304e19;
  double sum = 0.0, hitst, xmax, tmp;

  if (n <= 2)
  {
    if (n == 2)
    {
      double xv = fabs(dx[0] - dy[0]);
      double yv = fabs(dx[1] - dy[1]);
      if (xv > yv)
      {
        double h = yv/xv;
        return xv*sqrt(1.0 + h*h);
      }
      else if (yv)
      {
        double h = xv/yv;
        return yv*sqrt(1.0 + h*h);
      }
      else
        return 0;
    }
    else if (n == 1)
      return fabs(*dx - *dy);
    else
      return 0.;
  }

  hitst = cuthi / (double) n;
  i = 0;

  /* Zero Sum. */
  while (i < n && *dx == *dy)
  {
    ++dx;
    ++dy;
    ++i;
  }
  if (i >= n)
    return (sum);

  if ((tmp = fabs(*dx - *dy)) > cutlo)
  {
START:
    for (; i < n; ++dx, ++dy, ++i)
    {               /* Loop over elements. */
      if ((tmp = fabs(*dx - *dy)) > hitst)
      {
        sum = 1.0 + (sum / tmp) / tmp;  /* Rescale and process. */
        xmax = tmp;
        ++dx;
        ++dy;
        ++i;
        for (; i < n; ++dx, ++dy, ++i)
        {
          if ((tmp = fabs(*dx - *dy)) > xmax)
          {
            sum = 1.0 + sum * (xmax / tmp)*(xmax / tmp);
            xmax = tmp;
            continue;
          }
          sum += (tmp / xmax)*(tmp / xmax);
        }
        return (sqrt(sum) * xmax);
      }
      sum += tmp * tmp;
    }
    return sqrt(sum);       /* Sum completed normally. */
  }
  else
  {                       /* Small sum prepare for phase 2. */
    xmax = tmp;
    ++dx;
    ++dy;
    ++i;
    sum += 1.0;
    for (; i < n; ++dx, ++dy, ++i)
    {
      if (fabs(*dx - *dy) > cutlo)
      {       /* Got normal element. Rescale and process. */
        sum = (sum * xmax) * xmax;
        goto START;
      }
      if ((tmp = fabs(*dx - *dy)) > xmax)
      {
        sum = 1.0 + sum * (xmax / tmp)*(xmax / tmp);
        xmax = tmp;
        continue;
      }
      sum += (tmp / xmax)*(tmp / xmax);
    }
    return sqrt(sum) * xmax;
  }
}

/*------------------------------------------------------------------*/
/*   Compute cosine of angle between V and vector P1P2.             */
/*------------------------------------------------------------------*/

double scalar(double *p1, double *p2, double *V, int ndim)
{
  double Q[3];
  int i;

  for (i = 0; i < ndim; i++)
    Q[i] = p2[i] - p1[i];

  return -dot(V, Q, ndim) / (1.e-10 + l2norm(Q, ndim) * l2norm(V, ndim));
}

/*-------------------------------------------------------------*/
/*  Project the point x onto the simplex                       */
/*  l <= x[0] <= ... <= x[n-1] <= u.                           */
/*-------------------------------------------------------------*/

void simplex(int n, double *x, double l, double u)
{
  int i, j, k=0, kl, ku;
  double t=0.0;

  kl = 0;
  while (kl < n && x[kl] <= l)
    x[kl++] = l;

  if (kl == n)
    return;

  ku = n - 1;
  while (ku >= 0 && x[ku] >= u)
    x[ku--] = u;

  if (ku < 0)
    return;

  for (j = kl+1; j <= ku; )
  {
    for (; j <= ku; j++)
    {
      if (x[j-1] > x[j])
      {
        t = 0.5 * (x[j-1] + x[j]);
        x[j-1] = x[j] = t;
        j++;
        k = 2;
        break;
      }
    }
    for (; j <= ku && t > x[j]; j++)
    {
      t =  (k * t + x[j]) / (double) (k + 1);
      for (i = j; i >= j-k; i--)
        x[i] = t;
      k++;
    }
  }
}

/* Zero a double vector. */

void dzero(int n, double *x)
{
  memset(x, 0, n * sizeof(double));
}

/* Copy lower triangular A to lower triangular B when both are allocated */
/* in contiguous memory. */

void lMcopy(int n, double **A, double **B)
{
  size_t size = (n * (n + 1)) >> 1;
  memcpy(B[0], A[0], size * sizeof(double));
}

/*--------------------------------------------------------------------*/
/*  Solve the lower triangular system Lz = b.  Overwrite b with z.    */
/*  Note that L[i][i] = 1 for all i.                                  */
/*--------------------------------------------------------------------*/

void lower_backsub(int n, double **L, double *b)
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

/*--------------------------------------------------------------------*/
/*  Solve the system Ax = LDL'x = b.  Overwrite b with x.             */
/*  Note that L and D are stored in A:  D overwrites the diagonal of  */
/*  A and L overwrites the lower triangle of A (except the            */
/*  diagonal).  The diagonal of L is all ones.                        */
/*--------------------------------------------------------------------*/

int LDLTsolve(int n, double **A, double *b)
{
  int i;

  if (LDLT(n,A))
    return 1;
  lower_backsub(n, A, b);
  for (i = 0; i < n; i++)
    b[i] /= A[i][i];
  upper_backsub(n, A, b);

  return 0;
}

/*-------------------------------------------------------------------------*/
/*  Calculate the LDL' factorization A = L*D*L' of a symmetric, positive   */
/*  definite n X n matrix.  L is unit lower triangular and D is diagonal   */
/*  with positive diagonal elements.                                       */
/*  INPUT                                                                  */
/*    A  Array of dimension n X n containing the matrix to be factored.    */
/*       A must be symmetric and positive definite.                        */
/*    n  Size of the input matrix A.                                       */
/*  OUTPUT                                                                 */
/*    A  Array of dimension n X n containing L in its lower triangle and   */
/*       D on its diagonal.                                                */
/*-------------------------------------------------------------------------*/

int LDLT(int n, double **A)
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
        return 1;
      A[i][j] = (A[i][j] - sum) / vj;
    }
  }

  delete [] v;

  return 0;
}

/*--------------------------------------------------------------------*/
/*  Solve the upper triangular system L'x = y.  Overwrite y with x.   */
/*  Note that L[i][i] = 1 for all i.                                  */
/*--------------------------------------------------------------------*/

void upper_backsub(int n, double **L, double *y)
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

/*--------------------------------------------------------------------*/
/*  Multiply matrix A times vector x and store in y.                  */
/*--------------------------------------------------------------------*/

void multAx(int rows, int cols, double **A, double *x, double *y)
{
  int i, j;
  double tmp;

  for (i = 0; i < rows; i++)
  {
    tmp = 0.;
    for (j = 0; j < cols; j++)
      tmp += A[i][j] * x[j];
    y[i] = tmp;
  }
}

void insertion_sort_points(int n, double *t, double **x)
{
  int i, j, k;
  double v, P[3];

  for (i = 1; i < n; i++)
  {
    v = t[i];
    j = i;
    for (k = 0; k < 3; k++)
      P[k] = x[i][k];
    while (j && t[j-1] > v)
    {
      t[j] = t[j-1];
      for (k = 0; k < 3; k++)
        x[j][k] = x[j-1][k];
      j--;
    }
    t[j] = v;
    for (k = 0; k < 3; k++)
      x[j][k] = P[k];
  }
}

void CompensatePoint(int how, double prad, double zval, double *ball, double *ijk, double *comped)
{
  // coming in ijk[2] is kvalue, ijk[0] and ijk[1] must be adjusted
  // how: 1=probe comp XY,
  //      2=probe comp Z,
  //      4=slide back to plane,

  normalize(ijk, ijk, 2);

  double rest = sqrt(1.0 - ijk[2]*ijk[2]);  // length left over for I and J.

  ijk[0] *= rest;
  ijk[1] *= rest;  // now a unit vector
  //bugout(0, _T("%f,%f,%f,%f,%f,%f BALL CENTER XYZ IJK"), ball[0], ball[1], ball[2], ijk[0], ijk[1], ijk[2]);

  //if (how == 4)
  //  bugout(0, _T("10 %f %f M"), ball[0], ball[1]);

  // determine touch point

  comped[0] = ball[0];
  comped[1] = ball[1];
  comped[2] = ball[2];

  if (how & 1)
  {
    comped[0] -= prad*ijk[0];
    comped[1] -= prad*ijk[1];
  }


  if (how & 2)
    comped[2] -= prad*ijk[2];

  double save[3]; save[0] = comped[0]; save[1] = comped[1]; save[2] = comped[2];
  //bugout(0, _T("%f,%f,%f CONTACT POINT"), comped[0], comped[1], comped[2]);

  // push into plane

  if (how & 4)
  {
    double ij[3];
    normalize(ij, ijk, 2);

    double theta = asin(ijk[2]);
    double adxy = (comped[2]-zval)*tan(theta);
    //bugout(0, _T("Expected Full comp is %f"), prad/cos(theta));

    comped[0] += adxy*ij[0];
    comped[1] += adxy*ij[1];
    comped[2] = zval;
  }

  //if (how == 4)
  //{
  //  bugout(0, _T("10 %f %f P"), comped[0], comped[1]);
  //  bugout(0, _T("9 %f %f %f %f (%f)"), ball[0], ball[1], comped[0], comped[1], dist(ball, comped));
  //}
  //bugout(0, _T("%f,%f,%f In Plane"), comped[0], comped[1], comped[2]);

  //bugout(0, _T("7 %f %f ball"), ball[0], ball[1]);
  //bugout(0, _T("7 %f %f comped"), comped[0], comped[1]);
  //bugout(0, _T("9 %f %f %f %f"), ball[0], ball[1], comped[0], comped[1]);
  //bugout(0, _T("7 %f %f half"), save[0], save[1]);
  //bugout(0, _T("8 %f %f %f %f"), ball[0], ball[1], save[0], save[1]);
  //bugout(0, _T("9 %f %f %f %f"), save[0], save[1], comped[0], comped[1]);
}

// distance between a point and a line

double ptlinedist(double *lp, double *lv, double *pt, double *proj)
{
  double h = dist(lp, pt);
  double ij[2], v[2];
  normalize(ij, lv);
  v[0] = pt[0] - lp[0];
  v[1] = pt[1] - lp[1];
  double l = dot(ij, v);
  double d = h*h - l*l;

  if (proj)
  {
    proj[0] = lp[0] + l*ij[0];
    proj[1] = lp[1] + l*ij[1];
  }

  if (d <= 0.0)
    return 0.0;
  return sqrt(d);
}

// signed projected distance (distance from line point).
// optionally return point projection.

double projdist(double *lp, const double *lv, double *pt, double *proj)
{
  double ij[2], v[2];
  normalize(ij, lv);
  v[0] = pt[0] - lp[0];
  v[1] = pt[1] - lp[1];
  double l = dot(ij, v);
  if (proj)
  {
    proj[0] = lp[0] + l*ij[0];
    proj[1] = lp[1] + l*ij[1];
  }

  return l;
}

// Does line segment p1-p2 intersect with p3-p4?

int intersect(double *p1, double *p2, double *p3, double *p4)
{
  double denom = (p2[0]-p1[0])*(p3[1]-p4[1]) - (p2[1]-p1[1])*(p3[0]-p4[0]);
  if (denom == 0)
    return 0;

  double numer = (p3[0]-p1[0])*(p3[1]-p4[1]) - (p3[0]-p4[0])*(p3[1]-p1[1]);
  numer /= denom;
  if (numer < 0.0 || numer > 1.0)
    return 0;

  numer = (p2[0]-p1[0])*(p3[1]-p1[1]) - (p3[0]-p1[0])*(p2[1]-p1[1]);
  numer /= denom;
  if (numer < 0.0 || numer > 1.0)
    return 0;

  return 1;
}


double radius_3points(double *p1, double *p2, double *p3)
{
  double x, y, m1[2], m2[2], v1[2], v2[2];

  m1[0] = 0.5*(p1[0] + p2[0]);
  m1[1] = 0.5*(p1[1] + p2[1]);
  v1[0] = p1[1] - p2[1];
  v1[1] = p2[0] - p1[0];
  normalize(v1, v1);

  m2[0] = 0.5*(p2[0] + p3[0]);
  m2[1] = 0.5*(p2[1] + p3[1]);
  v2[0] = p2[1] - p3[1];
  v2[1] = p3[0] - p2[0];
  normalize(v2, v2);

  if (!int_line_line(m1[0], m1[1], v1[0], v1[1], m2[0], m2[1], v2[0], v2[1],
    &x, &y))
    return 1.0e20;

  return _hypot(m1[0]-x, m1[1]-y);
}

// This function determines if the ray that starts at rayStart with direction
// rayDir intersects the line segment that joins the points p1 and p2.

int int_ray_segment(double *rayStart, double *rayDir, double *p1, double *p2,
                    double *inter)
{
  // used when there is an intersection

  double rayEnd[2];  // push it way out and use segment/segment intersection.
  rayEnd[0] = rayStart[0] + 1000*rayDir[0];
  rayEnd[1] = rayStart[1] + 1000*rayDir[1];

  // construct vectors from rayStart to the line segment end points.

  double dir1[2], dir2[2];
  dir1[0] = p1[0] - rayStart[0];
  dir1[1] = p1[1] - rayStart[1];
  dir2[0] = p2[0] - rayStart[0];
  dir2[1] = p2[1] - rayStart[1];

  // these are positive when the the associated end is on plus side of the ray
  double dot1 = dot(dir1, rayDir);
  double dot2 = dot(dir2, rayDir);

  // if both end points are wrong side of ray, the whole segment is
  if (dot1 <= 0.0 && dot2 <= 0.0)
    return 0;

  // examine curl between dir1 and raydir and between raydir and dir2.
  // If these are the same sign then the ray intersects the segment.
  // If the first curl is zero return 1.  If the second curl is zero return
  // 0.  This is unlikely, but avoids counting twice in loops.

  // these will be positive if the angle sweep between the two vectors is CCW.
  double curl1 = curl(dir1, rayDir);
  double curl2 = curl(rayDir, dir2);

  if (curl1 == 0.0)
  {
    if (inter)
    {
      inter[0] = p1[0];  // is the first point
      inter[1] = p1[1];
    }
    return 1;
  }

  if (curl2 == 0.0)
    return 0;

  if (curl1*curl2 < 0.0)
    return 0;

  // At this point we know there is an intersection.  If dot1 and dot2 are
  // non-negative, the intersection must be on the positive side of the ray.

  if (dot1 >= 0.0 && dot2 >= 0.0)
  {
    if (inter)
    {
      double dir[2];
      dir[0] = p2[0] - p1[0];
      dir[1] = p2[1] - p1[1];
      normalize(dir, dir);
      return int_line_line(rayStart[0], rayStart[1], rayDir[0], rayDir[1],
        p1[0], p1[1], dir[0], dir[1], &inter[0], &inter[1]);
    }
    else
      return 1;
  }

  // if we are still here, the ray intersects the segment and the segment
  // straddles the boundary of the ray half plane this isn't very likely as
  // this function is used, but we deal with it anyway.  Need to find the
  // intersection point and see if it is on the positive side.  Note the
  // intersection could have been done straight away and examined, but
  // the code above should be fast and we aren't concerned what the
  // intersection point actually is.


  if (inter)
  {
    double dir[2];
    dir[0] = p2[0] - p1[0];
    dir[1] = p2[1] - p1[1];
    normalize(dir, dir);
    return int_line_line(rayStart[0], rayStart[1], rayDir[0], rayDir[1],
      p1[0], p1[1], dir[0], dir[1], &inter[0], &inter[1]);
  }
  else
    return intersect(rayStart, rayEnd, p1, p2);
}

// Routine to invent points to add a phantom end to a partial section.
// n1 and n2 are the indices where the phanom points need to be created.
// nold is the size of the arrays x and y that contain the existing points.
// n is the number of phantom points to create.  nx and ny are the arrays
// that will contain the phantom points.  Calling routine will reassemble
// the points

int InventEnd(int n1, int n2, int nold, double *x, double *y, int n, double *nx, double *ny, double mult,
              double *newstart, double *newend)
{
  //bugout(3, _T("InventEnd called"));
  if (n < 8 || nold < 4)
    return 0;

  double d, length, delta=0, offset, v1[2], v2[2];

  int loopMax = (n - 2)/4;

  if (mult <= 0.0)
  {
    // try making delta the same for each side of corner.  having some trouble with bulging.

    v1[0] = x[n1] - x[(n1 + nold - 1)%nold];
    v1[1] = y[n1] - y[(n1 + nold - 1)%nold];
    length = normalize(v1, v1);
    delta = 0.4*length/loopMax;
    //bugout(0, _T("length %f delta %f"), length, delta);

    v1[0] = x[n2] - x[(n2 + nold + 1)%nold];
    v1[1] = y[n2] - y[(n2 + nold + 1)%nold];
    length = normalize(v1, v1);

    d = 0.4*length/loopMax;
    if (d < delta)
      delta = d;
    //bugout(0, _T("length %f d %f delta %f"), length, d, delta);

    v1[0] = x[n2] - x[n1];     // vector to close off end
    v1[1] = y[n2] - y[n1];
    length = normalize(v1, v1);

    d = 0.5*length/(n-2*loopMax-1);
    if (d < delta)
      delta = d;
    //bugout(0, _T("length %f d %f delta %f"), length, d, delta);
  }

  //bugout(0, _T("9 %f %f %f %f invent begin"), x[(n1 + nold - 1)%nold], y[(n1 + nold - 1)%nold], x[n1], y[n1]);
  
  // points to follow n1;

  // use same delta on both sides...

  v1[0] = x[n2] - x[(n2 + nold + 1)%nold];
  v1[1] = y[n2] - y[(n2 + nold + 1)%nold];
  double len2 = normalize(v1, v1);


  v1[0] = x[n1] - x[(n1 + nold - 1)%nold];
  v1[1] = y[n1] - y[(n1 + nold - 1)%nold];
  length = normalize(v1, v1);

  if (len2 < length)
    length = len2;

  if (mult <= 0.0)
  {
    //delta = 0.2*length/loopMax;
    if (newstart)
    {
      newstart[0] = x[n1] - (loopMax + 1)*delta*v1[0];
      newstart[1] = y[n1] - (loopMax + 1)*delta*v1[1];
    }

    for (int i=0; i <= loopMax; i++)
    {
      nx[i] = x[n1] - (loopMax - i)*delta*v1[0];
      ny[i] = y[n1] - (loopMax - i)*delta*v1[1];
    }
  }
  else
  {
    delta = mult*length;

    for (int i=0; i <= loopMax; i++)
    {
      nx[i] = x[n1] + (i + 1)*delta*v1[0];
      ny[i] = y[n1] + (i + 1)*delta*v1[1];
    }
  }

  // points to precede n2;

  v1[0] = x[n2] - x[(n2 + nold + 1)%nold];
  v1[1] = y[n2] - y[(n2 + nold + 1)%nold];
  length = normalize(v1, v1);

  //bugout(0, _T("9 %f %f %f %f invent begin"), x[(n2 + nold + 1)%nold], y[(n2 + nold + 1)%nold], x[n2], y[n2]);

  if (mult <= 0.0)
  {
    //delta = 0.2*length/loopMax;
    if (newend)
    {
      newend[0] = x[n2] - (loopMax+1)*delta*v1[0];
      newend[1] = y[n2] - (loopMax+1)*delta*v1[1];
    }

    for (int i=0; i<= loopMax; i++)
    {
      nx[n-i-1] = x[n2] - (loopMax - i)*delta*v1[0];
      ny[n-i-1] = y[n2] - (loopMax - i)*delta*v1[1];
    }
  }
  else
  {
    delta = mult*length;
    for (int i=0; i<= loopMax; i++)
    {
      nx[n-i-1] = x[n2] + (i+1)*delta*v1[0];
      ny[n-i-1] = y[n2] + (i+1)*delta*v1[1];
    }
  }

  v2[0] = nx[n-loopMax-1] - nx[loopMax];     // vector to close off end
  v2[1] = ny[n-loopMax-1] - ny[loopMax];
  length = normalize(v2, v2);

  int i;
  if (mult <= 0.0)
  {
    //delta = 0.5*length/(n-2*loopMax-1);

    for (i=loopMax+1; i<n-loopMax-1; i++)
    {
      if (i < n/2)
      {
        nx[i] = nx[loopMax] + (i-loopMax)*delta*v2[0];
        ny[i] = ny[loopMax] + (i-loopMax)*delta*v2[1];
      }
      else if (n%2 && i == n/2)  // if odd number of pionts use midpint
      {
        nx[i] = 0.5*(x[n1] + x[n2]);
        ny[i] = 0.5*(y[n1] + y[n2]);
      }
      else
      {
        nx[i] = nx[n-loopMax-1] - (n-loopMax-1-i)*delta*v2[0];
        ny[i] = ny[n-loopMax-1] - (n-loopMax-1-i)*delta*v2[1];
      }
    }
  }
  else
  {
    offset = length/(n - 2*loopMax-1);

    for (i=loopMax+1; i<n-loopMax-1; i++)
    {
      nx[i] = nx[loopMax] + (i-loopMax)*offset*v2[0] + delta*v1[0];
      ny[i] = ny[loopMax] + (i-loopMax)*offset*v2[1] + delta*v1[1];
    }
  }
  //for (i=1; i<n; i++) bugout(3, _T("8 %f %f %f %f"), nx[i-1], ny[i-1], nx[i], ny[i]);
  //for (i=0; i<n; i++) bugout(0, _T("10 %f %f INV%d"), nx[i], ny[i], i);
  //if (newstart) bugout(0, _T("10 %f %f S"), newstart[0], newstart[1]);
  //if (newend) bugout(0, _T("10 %f %f E"), newend[0], newend[1]);

  return 1;
}

// is angle between the two line segments less than eps (radians)?

int Colinear(double *p1, double *p2, double *p3, double eps)
{
  double v1[2], v2[2];

  v1[0] = p2[0] - p1[0];
  v1[1] = p2[1] - p1[1];
  normalize(v1, v1);

  v2[0] = p3[0] - p1[0];
  v2[1] = p3[1] - p2[1];
  normalize(v2, v2);

  double angle = acos(dot(v1, v2));

  if (angle < eps)
    return 1;

  return 0;
}

// intersect circle and line segment

void int_circ_line(double *cir, double r, double *p1, double *p2, double *sol)
{
  double d = dist(cir, p1);
  double v[2];
  v[0] = p2[0] - p1[0];
  v[1] = p2[1] - p1[1];
  normalize(v, v);

  if (d < 1.0e-6)  // simple problem
  {
    sol[0] = cir[0] + r*v[0];
    sol[1] = cir[1] + r*v[1];
    return;
  }

  // more complex option
  double cir_p1[2];
  cir_p1[0] = p1[0] - cir[0];
  cir_p1[1] = p1[1] - cir[1];

  // setup quadratic

  double A = 1.0;
  double B = 2*dot(cir_p1, v);
  double C = d*d - r*r;

  double radical = B*B - 4*A*C;
  if (radical <= 0)   // something is goofy
  {
    sol[0] = p1[0];  // return something
    sol[1] = p1[1];
    return;
  }

  double alpha, alpha2;

  alpha = (-B + sqrt(radical))/(2*A);
  alpha2 = (-B - sqrt(radical))/(2*A);

  // want smallest positve solution, or closest to zero if both negative
  if (alpha < 0.0 && alpha2 > 0.0)
    alpha = alpha2;
  else if (alpha*alpha2 > 0.0 && fabs(alpha2) < fabs(alpha))
    alpha = alpha2;

  sol[0] = p1[0] + alpha*v[0];
  sol[1] = p1[1] + alpha*v[1];
}

bool CalculateLine(CMatrix *points, double *center, double *theta)
{
  int i;

  double inTheta = *theta;  // angle passed in should be approximate solution, used to set direction of line

  int numpts = points->R();

  if (numpts < 2)
    return false;

  center[0] = center[1] = 0.0;
  for(i=0 ; i<numpts ; i++)   // calc centroid
  {
    center[0] += points->m[i][0];
    center[1] += points->m[i][1];
  }
  center[0] /= numpts;
  center[1] /= numpts;

  // compute the sums used in the regression
  double sum1 = 0.0, sum2 = 0.0;
  double xp,yp;
  for(i=0; i<numpts ; i++)
  {
    xp = points->m[i][0] - center[0];
    yp = points->m[i][1] - center[1];
    sum1 += xp*yp;
    sum2 += xp*xp - yp*yp;
  }

  // angle of line for (x,y)-points
  if (fabs(sum2) < 1.0e-12)
    return false;

  *theta = 0.5*atan(2.0*sum1/sum2); 

  double s1 = sin(*theta);
  double c1 = cos(*theta);
  double s2 = sin(*theta + M_PI_2);
  double c2 = cos(*theta + M_PI_2);

  for(i=0,sum1=0.0,sum2=0.0 ; i<numpts ; i++)
  {
    xp = points->m[i][0] - center[0];
    yp = points->m[i][1] - center[1];
    double temp = xp*s1 - yp*c1;
    sum1 += temp*temp;
    temp = xp*s2 - yp*c2;
    sum2 += temp*temp;
  }

  *theta = (sum1 > sum2 ? *theta + M_PI_2 : *theta);
  if (*theta < 0.0)
    *theta += 2.0*M_PI;

  if (fabs(inTheta) >= 0 && fabs(inTheta) < 2.0*M_PI)  // if theta passed in is in [0, 2*M_PI) do direction check
  {
    if (fabs(*theta - inTheta) > M_PI)  // pointing in wrong direction
    {
      *theta += M_PI;

      if (*theta > 2.0*M_PI)  // put back into [0, 2.0*M_PI) if necessary
        *theta -= 2.0*M_PI;
    }
  }

  return true;
}

// following is functional rather than orthogonal line regression
// SLOPE  =  ( n*S(x*y) - Sx*Sy ) / ( n*S(x*x) ?(Sx)*(Sx) )

bool CalculateFunctionalLine(CMatrix *points, double *center, double *theta)
{
  int i;

  int numpts = points->R();

  if (numpts < 2)
    return false;

  center[0] = center[1] = 0.0;
  for(i=0 ; i<numpts ; i++)   // calc centroid
  {
    center[0] += points->m[i][0];
    center[1] += points->m[i][1];
  }
  center[0] /= numpts;
  center[1] /= numpts;

  // compute the sums used in the regression
  double sumXY = 0.0, sumX = 0.0, sumY=0.0, sumXX = 0.0;
  double xp,yp;
  for(i=0; i<numpts ; i++)
  {
    xp = points->m[i][0]/* - center[0]*/;
    yp = points->m[i][1]/* - center[1]*/;
    sumX += xp;
    sumY += yp;
    sumXY += xp*yp;
    sumXX += xp*xp;
  }

  if (fabs(numpts*sumXX - sumX*sumX) < 1.0e-12)
    return false;

  double slope = (numpts*sumXY - sumX*sumY) / (numpts*sumXX - sumX*sumX);

  *theta = atan(slope); 

  return true;
}

// check a stream of points, removing any that might cause a loop in the data.

// new method
int fix_data(int n, double *pX, double *pY, double *pK, int closed, double howClose, bool checkEnds, bool showDebug)
{
  if (n < 3)
    return n;

  int i, j;
  //first loop, look for points that are too close

  if (howClose < 1.0e-04)  // use 0.00001 if comes in as zero
    howClose = 1.0e-10;
  else
    howClose = howClose*howClose;  // square, so that we don't need to take square root in all the checks
  
  for (i=1; i<n; i++)
  {
    if ((pX[i]-pX[i-1])*(pX[i]-pX[i-1]) + (pY[i]-pY[i-1])*(pY[i]-pY[i-1]) > howClose)
      continue;

   /* if (showDebug)
      bugout(0, L"10 %f %f TC too close removal", pX[i], pY[i]);*/
    
    // point i is too close to i-1 so remove it.
    for (j=i; j < n-1; j++)
    {
      pX[j] = pX[j+1];
      pY[j] = pY[j+1];
      if (pK != NULL)
        pK[j] = pK[j+1];
    }
    n--;
    i--; // so this index gets checked again in next iteration
  }
  if (n < 3)
    return n;

  double nrm1, nrm2, nrm3, c, c2, v1[2], v2[2], v3[2];

  // check that first point isn't a hook
  if (checkEnds || !closed)
  {
    v1[0] = pX[1] - pX[0];  // first line segment
    v1[1] = pY[1] - pY[0];
    v2[0] = pX[2] - pX[1];  // second line segment
    v2[1] = pY[2] - pY[1];

    nrm1 = l2norm(v1, 2);
    nrm2 = l2norm(v2, 2);
    if (nrm1 * nrm2 != 0.)
      c = dot(v1, v2, 2) / (nrm1 * nrm2); // Cosine of angle
    else
      c = -1.0;

    if (c < 0.0)
    {
      /*if (showDebug)
        bugout(0, L"10 %f %f HS hook start removal", pX[0], pY[0]);*/

      for (j = 0; j < n-1; j++)
      {
        pX[j] = pX[j+1];
        pY[j] = pY[j+1];
        if (pK != NULL)
          pK[j] = pK[j+1];
      }
      n--;
    }
 
     if (n < 3)
      return n;

    // Check for hook at end of data.

    v1[0] = pX[n-3] - pX[n-2];
    v1[1] = pY[n-3] - pY[n-2];
    v2[0] = pX[n-2] - pX[n-1];
    v2[1] = pY[n-2] - pY[n-1];

    nrm1 = l2norm(v1,2);
    nrm2 = l2norm(v2,2);
    if (nrm1 * nrm2 != 0.)
      c = dot(v1, v2, 2) / (nrm1 * nrm2); // Cosine of angle
    else
      c = -1.0;

    if (c < 0.0)
    {
     /* if (showDebug)
        bugout(0, L"10 %f %f HE Hook End removal", pX[n-1], pY[n-1]);*/
      n--;
    }
  }
  
  if (n < 3)
    return n;

  v1[0] = pX[1] - pX[0];  // first line segment
  v1[1] = pY[1] - pY[0];

  for (i=1; i<n; i++)
  {
    if (i == n-1 && !closed) // done
      break;
      
    // segment before current point
    if (i > 1)  // for i == 0, v1 computed before loop.
    {
      v1[0] = v2[0];  // copy segment 2 to be the new segment 1
      v1[1] = v2[1];
    }
    
    // segment after current point
    v2[0] = pX[(i+1)%n] - pX[i];  // line segment after current point
    v2[1] = pY[(i+1)%n] - pY[i];
    
    nrm1 = l2norm(v1, 2);
    nrm2 = l2norm(v2, 2);
    if (nrm1 * nrm2 != 0.)
      c = dot(v1, v2, 2) / (nrm1 * nrm2); // Cosine of angle
    else
      c = -1.0;

    if (c > 0.0) // okay
      continue;

      //bugout(0, L"10 %f %f Q hook?", pX[(i+1)%n], pY[(i+1)%n]);
      
    // if we are here, then the two segments that include the current point form an angle that is sharper than 90 degrees.
    // we need to look at the next segment to see if we do indeed have a single sharp point.  If we do, then we want to keep the point.
    // but if the direction changes again soon, then we want to throw this poin1t away.

    int removeIndex = -1;

    // looking for zig zag shape.  when the probe catches on the part and then springs back, it may take a few points before we see the zig zag
    // check the next ten points to see if we see it.  if we have a smaller number of points, don't look so far ahead.
    
    int ahead = n/10;
    if (ahead > 10)
      ahead = 10;
      
    for (j=0; j<ahead; j++)
    {
      v3[0] = pX[(i+j+2)%n] - pX[(i+1)%n];  // line segment after next point
      v3[1] = pY[(i+j+2)%n] - pY[(i+1)%n];
    
      nrm3 = l2norm(v3, 2);
      if (nrm2 * nrm3 != 0.)
        c2 = dot(v2, v3, 2) / (nrm2 * nrm3); // Cosine of angle
      else
        c2 = -1.0;
    
      if (c2 < 0.0)  // zig zag action
      {
        removeIndex = (i + 1)%n;
        break;
      }
    }

    if (removeIndex < 0)
      continue;
  
      //bugout(0, L"10 %f %f SP=%d spike removal", pX[removeIndex], pY[removeIndex], j);
          
    // point i is too close to i-1 so remove it.
    for (j=removeIndex; j < n-1; j++)
    {
      pX[j] = pX[j+1];
      pY[j] = pY[j+1];
    }

    v2[0] = pX[(i+1)%n] - pX[i];  // compute this again with new next point so will be good for next iteration
    v2[1] = pY[(i+1)%n] - pY[i];

    n--;
    i--; // so this index gets checked again in next iteration
  }


  return n;
}

// new method
int fix_data_new(int n, double* pX, double* pY, double* pK, int closed, double howClose, bool checkEnds, bool showDebug,
                 double theat)
{

  if(n < 3)
    return n;

  int i, j;
  // first loop, look for points that are too close

  if(howClose < 1.0e-04) // use 0.00001 if comes in as zero
    howClose = 1.0e-10;
  else
    howClose = howClose * howClose; // square, so that we don't need to take square root in all the checks

  for(i = 1; i < n; i++)
  {
    if((pX[i] - pX[i - 1]) * (pX[i] - pX[i - 1]) + (pY[i] - pY[i - 1]) * (pY[i] - pY[i - 1]) > howClose)
      continue;

    //if(showDebug)
    //  bugout(0, L"fix_data_new:cur Point(%f ,%f), pre Point(%f ,%f) dist=%lf, howClose=%lf,TC too close removal", pX[i], pY[i], pX[i - 1],
    //         pY[i - 1], ((pX[i] - pX[i - 1]) * (pX[i] - pX[i - 1]) + (pY[i] - pY[i - 1]) * (pY[i] - pY[i - 1])),
    //         howClose);

    // point i is too close to i-1 so remove it.
    for(j = i; j < n - 1; j++)
    {
      pX[j] = pX[j + 1];
      pY[j] = pY[j + 1];
      if(pK != NULL)
        pK[j] = pK[j + 1];
    }
    n--;
    i--; // so this index gets checked again in next iteration
  }
  if(n < 3)
    return n;

  double nrm1, nrm2, nrm3, c, c2, v1[2], v2[2], v3[2];

  // check that first point isn't a hook
  if(checkEnds || !closed)
  {
    v1[0] = pX[1] - pX[0]; // first line segment
    v1[1] = pY[1] - pY[0];
    v2[0] = pX[2] - pX[1]; // second line segment
    v2[1] = pY[2] - pY[1];

    nrm1 = l2norm(v1, 2);
    nrm2 = l2norm(v2, 2);
    if(nrm1 * nrm2 != 0.)
      c = dot(v1, v2, 2) / (nrm1 * nrm2); // Cosine of angle
    else
      c = -1.0;

    if(c < 0.0) // nyc
    {
      /*if(showDebug)
        bugout(0, L"10 %f %f HS hook start removal", pX[0], pY[0]);*/

      for(j = 0; j < n - 1; j++)
      {
        pX[j] = pX[j + 1];
        pY[j] = pY[j + 1];
        if(pK != NULL)
          pK[j] = pK[j + 1];
      }
      n--;
    }

    if(n < 3)
      return n;

    // Check for hook at end of data.

    v1[0] = pX[n - 3] - pX[n - 2];
    v1[1] = pY[n - 3] - pY[n - 2];
    v2[0] = pX[n - 2] - pX[n - 1];
    v2[1] = pY[n - 2] - pY[n - 1];

    nrm1 = l2norm(v1, 2);
    nrm2 = l2norm(v2, 2);
    if(nrm1 * nrm2 != 0.)
      c = dot(v1, v2, 2) / (nrm1 * nrm2); // Cosine of angle
    else
      c = -1.0;

    if(c < 0.0) // nyc
    {
      /*if(showDebug)
        bugout(0, L"10 %f %f HE Hook End removal", pX[n - 1], pY[n - 1]);*/
      n--;
    }
  }

  if(n < 3)
    return n;

  v1[0] = pX[1] - pX[0]; // first line segment
  v1[1] = pY[1] - pY[0];

  for(i = 1; i < n; i++)
  {
    double point_spacing = 0;
    point_spacing = sqrt((pX[i] - pX[i - 1]) * (pX[i] - pX[i - 1]) + (pY[i] - pY[i - 1]) * (pY[i] - pY[i - 1]));
    //bugout(0, L"fix_data_new:cur Point(%f ,%f), pre Point(%f ,%f), Point dist=%lf,i=%d ", pX[i],pY[i], pX[i - 1], pY[i - 1],
    //       point_spacing,i);
    if(i == n - 1 && !closed) // done
      break;

    // segment before current point
    if(i > 1) // for i == 0, v1 computed before loop.
    {
      v1[0] = v2[0]; // copy segment 2 to be the new segment 1
      v1[1] = v2[1];
    }

    // segment after current point
    v2[0] = pX[(i + 1) % n] - pX[i]; // line segment after current point
    v2[1] = pY[(i + 1) % n] - pY[i];

    nrm1 = l2norm(v1, 2);
    nrm2 = l2norm(v2, 2);
    if(nrm1 * nrm2 != 0.)
      c = dot(v1, v2, 2) / (nrm1 * nrm2); // Cosine of angle
    else
      c = -1.0;

    if(c > theat) // okay niu
      continue;
    double nrm4, nrm5, v4[2], v5[2], c4;
    bool Deal = false;
    for(j = 0; j < 5; j++)
    {
      v4[0] = pX[(i + 1) % n] - pX[(i - 1) % n];
      v4[1] = pY[(i + 1) % n] - pY[(i - 1) % n];
      v5[0] = pX[(i + j + 2) % n] - pX[(i + 1) % n];
      v5[1] = pY[(i + j + 2) % n] - pY[(i + 1) % n];
      nrm4 = l2norm(v4, 2);
      nrm5 = l2norm(v5, 2);
      if(nrm4 * nrm5 != 0.)
        c4 = dot(v4, v5, 2) / (nrm4 * nrm5); // Cosine of angle
      else
        c4 = -1.0;
      if(c4 > 0.0) //
      {
        /*for(int k = i; k < n - 1; k++)
        { 
          pX[k] = pX[k + 1];
          pY[k] = pY[k + 1];
        }*/
        //bugout(0, L"fix_data_new:do avg point( %f ,%f) index=%d ****", pX[i], pY[i], i); 
        pX[i] = (pX[i - 1] + pX[i + 1])/2;
        pY[i] = (pY[i - 1] + pY[i + 1])/2;
        v2[0] = pX[(i + 1) % n] - pX[i]; // compute this again with new next point so will be good for next iteration
        v2[1] = pY[(i + 1) % n] - pY[i];

       // n--;
       // i--; // so this index gets checked again in next iteration

        Deal = true;
        break;
      }
    }
    if(Deal)
    {
      continue;
    }

    //if(showDebug)
    //  bugout(0, L"10 %f %f Q hook****", pX[(i + 1) % n], pY[(i + 1) % n]);

    // if we are here, then the two segments that include the current point form an angle that is sharper than 90
    // degrees. we need to look at the next segment to see if we do indeed have a single sharp point.  If we do, then we
    // want to keep the point. but if the direction changes again soon, then we want to throw this poin1t away.

    int removeIndex = -1;

    // looking for zig zag shape.  when the probe catches on the part and then springs back, it may take a few points
    // before we see the zig zag check the next ten points to see if we see it.  if we have a smaller number of points,
    // don't look so far ahead.
    int max_X_index = -1;
    int max_Y_index = -1;
    int min_X_index = -1;
    int min_Y_index = -1;
    double maxXval = 1.0e-5;
    double maxYval = 1.0e-5;
    double minXval = 1.0e5;
    double minYval = 1.0e5;
    for(int m = 0; m < n; m++)
    {
      if(pX[m] > maxXval)
      {
        maxXval = pX[m];
        max_X_index = m;
      }
      if(pX[m] < minXval)
      {
        minXval = pX[m];
        min_X_index = m;
      }
      if(pY[m] > maxYval)
      {
        maxYval = pY[m];
        max_Y_index = m;
      }
      if(pY[m] < minYval)
      {
        minYval = pY[m];
        min_Y_index = m;
      }
    }
    int ahead = n / 10;
    if(ahead > 4)//4 or 10
      ahead = 4;

    for(j = 0; j < ahead; j++)
    {
      if(pX[(i + j + 2) % n] == maxXval || pY[(i + j + 2) % n] == maxXval || pY[(i + j + 2) % n] == minXval ||
         pY[(i + j + 2) % n] == minYval)
      {
        break;
      }

      v3[0] = pX[(i + j + 2) % n] - pX[(i + 1) % n]; // line segment after next point
      v3[1] = pY[(i + j + 2) % n] - pY[(i + 1) % n];

      nrm3 = l2norm(v3, 2);
      if(nrm2 * nrm3 != 0.)
        c2 = dot(v2, v3, 2) / (nrm2 * nrm3); // Cosine of angle
      else
        c2 = -1.0;

      if(c2 < theat) // zig zag action  nyc
      {
        removeIndex = (i + 1) % n;
     /*   bugout(0, L"c2(%lf) < theat(%lf),Vector angle=%lf,removeIndex=%d", c2, theat, acos(theat) * 180.0 / M_PI,
               removeIndex);*/
        break;
      }
    }

    if(removeIndex < 0)
      continue;

   /* if(showDebug)
      bugout(0, L"will remove( %f ,%f) SP=%d spike removal ****nyc", pX[removeIndex], pY[removeIndex], j);*/

    // point i is too close to i-1 so remove it.
    for(j = removeIndex; j < n - 1; j++)
    {
      pX[j] = pX[j + 1];
      pY[j] = pY[j + 1];
    }
    

     //pX[removeIndex] = (pX[removeIndex - 1] + pX[removeIndex + 1])/2;
     //pY[removeIndex] = (pY[removeIndex - 1] + pY[removeIndex + 1])/2;
    //bugout(0, L"do avg point( %f ,%f) SP=%d ****nyc", pX[removeIndex], pY[removeIndex], j); 
    v2[0] = pX[(i + 1) % n] - pX[i]; // compute this again with new next point so will be good for next iteration
    v2[1] = pY[(i + 1) % n] - pY[i];

    n--;
    i--; // so this index gets checked again in next iteration
  }
 // bugout(0, L"fix_data: will return, n=%d", n);
  return n;
}
/* old method
int fix_data(int n, double *pX, double *pY, double *pK, int closed, double howClose, bool checkEnds, bool showDebug)
{
  int i, j;
  double P[3], P1[3], c, nrm1, nrm2;

  if (n < 3)
    return n;

  // first check for points that are too close

  for (i = 0; i < n; i++)
  {
    if (i == n-1 && !closed)
      break;

    P[0] = pX[(i+1)%n] - pX[i];
    P[1] = pY[(i+1)%n] - pY[i];

    nrm1 = l2norm(P,2);

    if (nrm1 < howClose)    // point i is too close to point i+1
    {
      if (showDebug)
        bugout(0, _T("10 %f %f TC too close removal"), pX[(i+1)%n], pY[(i+1)%n]);
      for (j = i+1; j < n-1; j++)  // remove point i+1
      {
        pX[j] = pX[j+1];
        pY[j] = pY[j+1];
        if (pK != NULL)
          pK[j] = pK[j+1];
      }
      n--;
      i--;   // check i again, will this cause problems?
    }
  }

  // Check for hook at start of data.

  P[0] = pX[0] - pX[1];
  P[1] = pY[0] - pY[1];

  P1[0] = pX[2] - pX[1];
  P1[1] = pY[2] - pY[1];

  nrm1 = l2norm(P,2);
  nrm2 = l2norm(P1,2);
  if (nrm1 * nrm2 != 0.)
    c = dot(P,P1,2) / (nrm1 * nrm2); // Cosine of angle
  else
    c = -1.0;

  if (c > 0.0 && checkEnds)
  {
    if (showDebug)
      bugout(0, _T("10 %f %f HS hook start removal"), pX[0], pY[0]);
    for (i = 1; i < n; i++)
    {
      pX[i-1] = pX[i];
      pY[i-1] = pY[i];
      if (pK != NULL)
        pK[i-1] = pK[i];
    }
    n--;
  }

  if (n < 3)
    return n;

  // Check for hook at end of data.

  P[0] = pX[n-1] - pX[n-2];
  P[1] = pY[n-1] - pY[n-2];

  P1[0] = pX[n-3] - pX[n-2];
  P1[1] = pY[n-3] - pY[n-2];

  nrm1 = l2norm(P,2);
  nrm2 = l2norm(P1,2);
  if (nrm1 * nrm2 != 0.)
    c = dot(P,P1,2) / (nrm1 * nrm2); // Cosine of angle
  else
    c = -1.0;

  if (c > 0.0 && checkEnds)
  {
    if (showDebug)
      bugout(0, _T("10 %f %f HE Hook End removal"), pX[n-1], pY[n-1]);
    n--;
  }

  // Delete spikes and close points

  for (i = 1; i < n-1; i++)
  {
    P[0] = pX[i+1] - pX[i];
    P[1] = pY[i+1] - pY[i];

    P1[0] = pX[i-1] - pX[i];
    P1[1] = pY[i-1] - pY[i];

    nrm1 = l2norm(P,2);
    nrm2 = l2norm(P1,2);

    if (nrm1 * nrm2 != 0.0)
      c = dot(P,P1,2) / (nrm1 * nrm2); // Cosine of angle
    else
      c = 1.0;


    if (c >= 0.0)  // Angle between 0 and pi/2
    {
      if (showDebug)
      {
        bugout(0, _T("8 %f %f %f %f c=%f"), pX[i], pY[i], pX[i-1], pY[i-1], c);
        bugout(0, _T("9 %f %f %f %f c=%f"), pX[i], pY[i], pX[i+1], pY[i+1], c);
      }

      double saveC = c;
      // Make sure this won't cause a bunch of points to get removed because it is really a sharp bend?

      if ((i > 1) && (i < n-1) && (nrm1*nrm2 != 0))
      {
        bool reversing = false;
        // first check curls, see if we are reversing direction more than once...

        if (i < n-2)
        {
          double v1[2], v2[2], v3[3];

          v1[0] = pX[i] - pX[i-1];
          v1[1] = pY[i] - pY[i-1];
          
          v2[0] = pX[i+1] - pX[i];
          v2[1] = pY[i+1] - pY[i];
          
          v3[0] = pX[i+2] - pX[i+1];
          v3[1] = pY[i+2] - pY[i+1];
          
          normalize(v1, v1);
          normalize(v2, v2);
          normalize(v3, v3);

          double curl1 = curl(v1, v2);
          double curl2 = curl(v2, v3);

          if (dot(v2, v3) < 0.5  && curl1*curl2 < 0.0)
            reversing = true;
        }

        // if this is removed will next one be removed next time around?

        if (!reversing)
        {
          P[0] = pX[i+1] - pX[i-1];
          P[1] = pY[i+1] - pY[i-1];
          P1[0] = pX[i-2] - pX[i-1];
          P1[1] = pY[i-2] - pY[i-1];

          nrm1 = l2norm(P,2);
          nrm2 = l2norm(P1,2);
          if (nrm1 * nrm2 != 0.0)
            c = dot(P,P1,2) / (nrm1 * nrm2); // Cosine of angle
          else
            c = -1.0;

          if (c > 0.0)   // next would also be bad, leave this alone.
            continue;
        }
      }

      // delete knot at i
      if (showDebug)
        bugout(0, _T("10 %f %f SP=%.1f spike removal"), pX[i], pY[i], saveC);

      for (j = i; j < n-1; j++)
      {
        pX[j] = pX[j+1];
        pY[j] = pY[j+1];
        if (pK != NULL)
          pK[j] = pK[j+1];
      }
      n--;
      i--;   // check i again, will this cause problems?
    }
  }

  return n;
}
*/

// 3D line/plane intersection. used in probe compensation, not used, didn't make any diffeence

bool int_plane_line(double *pvec, double pd, double *lpt, double *lvec, double *sol)
{
  // Compute the distance of a point on the line to the plane

  double temp1 = dot(pvec, lpt, 3) - pd;

  // Get the hypotenuse of the triangle.  Must use unit vectors to get the right answer.

  double temp2 = dot(pvec, lvec, 3);
  if( fabs(temp2) < 1.0e-8)
    return false;

  double lamda = temp1 / temp2;
  sol[0] = lpt[0] - lamda * lvec[0];
  sol[1] = lpt[1] - lamda * lvec[1];
  sol[2] = lpt[2] - lamda * lvec[2];

  //double cd = dot(sol, pvec, 3);
  //bugout(0, _T("intersect check pd %f check %f difference %f"), pd, cd, cd-pd);

  return true;
}

bool PointInPolygon(int num, double *vx, double *vy, double testx, double testy)
{
  int i, j;
  bool c = false;
  for (i = 0, j = num-1; i < num; j = i++)
  {
    if ( ((vy[i]>testy) != (vy[j]>testy)) && (testx < (vx[j]-vx[i]) * (testy-vy[i]) / (vy[j]-vy[i]) + vx[i]) )
       c = !c;
  }
  return c;
}

// this is used to see if start of LEC is same as end of CVC and tests like that.  The values should be either very close, or not close at all
// so the test value doesn't need to be all that small

bool ParametersEqual(double period, double p1, double p2)
{
  if (fabs(p1 - p2) < 0.0001)
    return true;

  double diff = 0.0;  // compute difference between two as a positive number
  if (p1 > p2)
    diff = fmod(p1 - p2, period);
  else
    diff = fmod(p2 - p1, period);

  if (diff < 0.0001)
    return true;

  return false;
}

// check to see if parameter p is in [r1,r2] modulo period
bool ParameterInRange(double p, double r1, double r2, double period)
{
  // 1 is probably sufficient, but using 5 to  make sure we have a positive values for first parameter of fmod
  double dp = fmod(p + 5*period - r1, period);  // difference between parameter and start of range
  double dw = r2 - r1;  // width of range
  if (dp < dw)
    return true;

  return false;
}

DLLEXPORT double signedAngleBetween2DVectors(const double* a, const double* b)
{
  double cross_component = a[0]*b[1] - a[1]*b[0];
  double dot_component = dot(a, b);
  return atan2(cross_component, dot_component);
}


// This function returns radians_unknown
// with an altered version with some multiple of (2*pi) to make the number
// as close as possible to radians_good
DLLEXPORT double UnwrapAngle(double radians_good, double radians_towrap)
{
  return radians_towrap + round((radians_good-radians_towrap) / M_2PI) * M_2PI;
}


// This function takes two termination-points and a center
// It alters the termination points adding/subtracting the period to them
// and figures out which one is the start and which is the end
std::pair<double, double> FitRangeToPoint(double termination1, double termination2, double x, double period)
{
  // Figure out the start
  double start1 = termination1 + floor((x - termination1) / period) * period;
  double start2 = termination2 + floor((x - termination2) / period) * period;
  double start = (start1 > start2) ? start1 : start2;
  // Figure out the end
  double end1 = termination1 + ceil((x - termination1) / period) * period;
  double end2 = termination2 + ceil((x - termination2) / period) * period;
  double end = (end1 < end2) ? end1 : end2;
  // all done
  return std::pair<double, double>(start, end);
}



/*
// This function returns a half-piece of a given section
// You give it a pointer to a Section object, and you tell it CVC or CCC
// and it gives you a CSubCurve for half the curve, first nominal and second measured
std::unique_ptr<CSubCurve> MakeNominalHalfCurve(const CSection *in_section, unsigned int which_piece)
{
  // check validity
  assert(which_piece == CVC || which_piece == CCC);

  // get the boundaries of the piece; identify a t-point on the desired half
  const double t0 = in_section->NomPart(which_piece)->T0();
  const double t1 = in_section->NomPart(which_piece)->T1();
  const double period = in_section->NomCurve()->Period();
  const double nom_t_included = 0.5 * (t0 + t1);
  const auto distanceToRounded = [](double x) -> double { return std::abs(x - round(x)); };

  // make sure at most one edge is partial
  if(in_section->LEType() == EDGE_PARTIAL && in_section->TEType() == EDGE_PARTIAL)
  {
    return std::make_unique<CSubCurve>(in_section->NomCurve(), t0, t1, period);
  }

  // identify the t-points of the nominal endpoints
  double nom_t_leading = in_section->NomPart(LEC)->Extreme();
  double nom_t_trailing = in_section->NomPart(TEC)->Extreme();

  // change things if edges are partial
  if(in_section->LEType() == EDGE_PARTIAL)
  {
    // pick the one farther away from the trailing point
    const double distance0 = distanceToRounded((t0 - nom_t_trailing) / period);
    const double distance1 = distanceToRounded((t1 - nom_t_trailing) / period);
    nom_t_leading = (distance0 > distance1) ? t0 : t1;
  }
  else if(in_section->TEType() == EDGE_PARTIAL)
  {
    // pick the one farther away from the leading point
    const double distance0 = distanceToRounded((t0 - nom_t_leading) / period);
    const double distance1 = distanceToRounded((t1 - nom_t_leading) / period);
    nom_t_trailing = (distance0 > distance1) ? t0 : t1;
  }

  // identify a t-point on the desired half
  // create the desired nominal t-range
  std::pair<double, double> nom_t_range = FitRangeToPoint(nom_t_leading, nom_t_trailing, nom_t_included, period);

  // build the subcurves
  return std::make_unique<CSubCurve>(in_section->NomCurve(), nom_t_range.first, nom_t_range.second, period);
}
std::unique_ptr<CSubCurve> MakeMeasuredHalfCurve(const CSection *in_section, unsigned int which_piece)
{
  // check validity
  assert(which_piece == CVC || which_piece == CCC);

  // get the boundaries of the piece; identify a t-point on the desired half
  const double t0 = in_section->MeaPart(which_piece)->T0();
  const double t1 = in_section->MeaPart(which_piece)->T1();
  const double period = in_section->MeaCurve()->Period();
  const double meas_t_included = 0.5 * (t0 + t1);
  const auto distanceToRounded = [](double x) -> double { return std::abs(x - round(x)); };

  // make sure at most one edge is partial
  if(in_section->LEType() == EDGE_PARTIAL && in_section->TEType() == EDGE_PARTIAL)
  {
    return std::make_unique<CSubCurve>(in_section->MeaCurve(), t0, t1, period);
  }

  // identify the t-points of the nominal endpoints
  double meas_t_leading = in_section->MeaPart(LEC)->Extreme();
  double meas_t_trailing = in_section->MeaPart(TEC)->Extreme();

  // change things if edges are partial
  if(in_section->LEType() == EDGE_PARTIAL)
  {
    // pick the one farther away from the trailing point
    const double distance0 = distanceToRounded((t0 - meas_t_trailing) / period);
    const double distance1 = distanceToRounded((t1 - meas_t_trailing) / period);
    meas_t_leading = (distance0 > distance1) ? t0 : t1;
  }
  else if(in_section->TEType() == EDGE_PARTIAL)
  {
    // pick the one farther away from the leading point
    const double distance0 = distanceToRounded((t0 - meas_t_leading) / period);
    const double distance1 = distanceToRounded((t1 - meas_t_leading) / period);
    meas_t_trailing = (distance0 > distance1) ? t0 : t1;
  }

  // identify a t-point on the desired half
  // create the desired nominal t-range
  std::pair<double, double> meas_t_range = FitRangeToPoint(meas_t_leading, meas_t_trailing, meas_t_included, period);

  // build the subcurves
  return std::make_unique<CSubCurve>(in_section->MeaCurve(), meas_t_range.first, meas_t_range.second, period);
}
 */

// compute 2D angle between two vectors

double SignedAngle(double *v1, double *v2)
{
  double nv1[2], nv2[2];
  
  double l1 = normalize(nv1, v1);
  double l2 = normalize(nv2, v2);
  if (fabs(l1) < 1.0e-5 || fabs(l2) < 1.0e-5)
    return 0.0;

  double dt = dot(nv1, nv2);
  
  if (dt >= 1.0)
    return 0.0;
  
  if (dt <= -1.0)
    return M_PI;
  
  double angle = acos(dt);
  if (curl(v1, v2) < 0)
    angle *= -1.0;
    
  return angle;
}
 
int TrimTails(int n1, int n2, double *x, double *y, double *z, double *vi, double *vj,
              int *end1, int *end2, int types, CCompRecordArray *cra, bool isEnglish)
{
  //int n1new = RemoveLoops(n1, x, y, z, vi, vj);
  //int n2new = RemoveLoops(n1, &x[n1], &y[n1], &z[n1], &vi[n1], &vj[n1]);

  //bugout(0, _T("in TrimTails types %x n1 %d n2 %d"), types, n1, n2);
  CMatrix sc1(n1,5);
  CMatrix sc2(n2,5);

  int parallel = 0;  // test if scans are essentially parralel or antiparallel
  double v1[2], v2[2];
  v1[0] = x[n1-1] - x[0];
  v1[1] = y[n1-1] - y[0];
  v2[0] = x[n1+n2-1] - x[n1];
  v2[1] = y[n1+n2-1] - y[n1];
  if (dot(v1,v2) > 0.0)
    parallel = 1;

  int showPoints = false;
  //if (IsCustomerID(L"HOW"))
  //  showPoints = true;

  int i;
  for (i=0; i<n1; i++)
  {
    sc1.m[i][0] = x[i];
    sc1.m[i][1] = y[i];
    sc1.m[i][2] = z[i];
    if (vi && vj)
    {
      sc1.m[i][3] = vi[i];
      sc1.m[i][4] = vj[i];
    }
    //bugout(0, _T("7 %f %f scan1 %d"), x[i], y[i], i);
    if (cra && showPoints)
    {
      CompRecord *compRecord = new CompRecord(10, x[i], y[i], 0.0, 0.0, i+1);
      cra->push_back(compRecord);
    }
    if (i > 0 && cra)
    {
      CompRecord *compRecord = new CompRecord(8, x[i-1], y[i-1], x[i], y[i], 0);
      cra->push_back(compRecord);
    }
  }
  for (i=0; i<n2; i++)
  {
    int j = parallel ? n2 - 1 - i : i;
    sc2.m[j][0] = x[n1 + i];
    sc2.m[j][1] = y[n1 + i];
    sc2.m[j][2] = z[n1 + i];
    if (vi && vj)
    {
      sc2.m[j][3] = vi[n1 + i];
      sc2.m[j][4] = vj[n1 + i];
    }
    //bugout(0, _T("7 %f %f scan2 %d"), x[n1+i], y[n1+i], j);

    if (cra && showPoints)
    {
      CompRecord *compRecord = new CompRecord(10, x[n1+i], y[n1+i], 0.0, 0.0, n1+i+1);
      cra->push_back(compRecord);
    }
    if (i > 0 && cra)
    {
      CompRecord *compRecord = new CompRecord( 8, x[n1+i-1], y[n1+i-1], x[n1+i], y[n1+i], 0);
      cra->push_back(compRecord);
    }
  }

  int s1, e1, s2, e2;

  int missed = 0, missingOk = 0;
  if (types & 4)  // le partial
    missingOk = 1;
  if (types & 16) // te partial
    missingOk++;

  if (!TrimEnd(n1, sc1, n2, sc2, &s1, &e2, 1, isEnglish))
  {
    s1 = 0;
    e2 = n2-1;
    missed++;

    // if this is a partial section:
    // make sure that points 0, 1, 2, ...  of sc1 don't have any wild points
    // make sure that points n2-1, n2-2, n2-3, ... of sc2 don't have any wild points
    if (missingOk)
    {
      int lim = n1/10;
      if (lim < 6)
        lim = 6;

      int newInter = -1;
      for (i=1; i<lim; i++)
      {
        double vv1[2], vv2[2];
        vv1[0] = sc1.m[i][0] - sc1.m[i-1][0];
        vv1[1] = sc1.m[i][1] - sc1.m[i-1][1];
        vv2[0] = sc1.m[i+1][0] - sc1.m[i][0];
        vv2[1] = sc1.m[i+1][1] - sc1.m[i][1];
        normalize(vv1, vv1);
        normalize(vv2, vv2);
        if (dot(vv1,vv2) < 0)
          newInter = i+1;
      }

      if (newInter > 0)
        s1 = newInter;

      lim = n2/10;
      if (lim < 6)
        lim = 6;

      newInter = -1;
      for (i=1; i<lim; i++)
      {
        int j = n2 - i - 1;

        double vv1[2], vv2[2];
        vv1[0] = sc2.m[j][0] - sc2.m[j-1][0];
        vv1[1] = sc2.m[j][1] - sc2.m[j-1][1];
        vv2[0] = sc2.m[j+1][0] - sc2.m[j][0];
        vv2[1] = sc2.m[j+1][1] - sc2.m[j][1];
        normalize(vv1, vv1);
        normalize(vv2, vv2);
        if (dot(vv1, vv2) < 0)
          newInter = j-1;
      }

      if (newInter > 0)
        e2 = newInter;
    }
  }

  if (!TrimEnd(n2, sc2, n1, sc1, &s2, &e1, 1, isEnglish))
  {
    s2 = 0;
    e1 = n1-1;
    missed++;

    // make sure that points 0, 1, 2, ...  of sc2 don't have any wild points
    // make sure that points n1-1, n1-2, n1-3, ... of sc1 don't have any wild points
    // if this is a partial section:
    if (missingOk)
    {
      int lim = n2/10;
      if (lim < 6)
        lim = 6;

      int newInter = -1;
      for (i=1; i<lim; i++)
      {
        double vv1[2], vv2[2];
        vv1[0] = sc2.m[i][0] - sc2.m[i-1][0];
        vv1[1] = sc2.m[i][1] - sc2.m[i-1][1];
        vv2[0] = sc2.m[i+1][0] - sc2.m[i][0];
        vv2[1] = sc2.m[i+1][1] - sc2.m[i][1];
        normalize(vv1, vv1);
        normalize(vv2, vv2);
        if (dot(vv1, vv2) < 0)
          newInter = i+1;
      }

      if (newInter > 0)
        s2 = newInter;

      lim = n1/10;
      if (lim < 6)
        lim = 6;

      newInter = -1;
      for (i=1; i<lim; i++)
      {
        int j = n1 - i - 1;

        double vv1[2], vv2[2];
        v1[0] = sc1.m[j][0] - sc1.m[j-1][0];
        v1[1] = sc1.m[j][1] - sc1.m[j-1][1];
        v2[0] = sc1.m[j+1][0] - sc1.m[j][0];
        v2[1] = sc1.m[j+1][1] - sc1.m[j][1];
        normalize(vv1, vv1);
        normalize(vv2, vv2);
        if (dot(vv1, vv2) < 0)
          newInter = j-1;
      }

      if (newInter > 0)
        e1 = newInter;
    }
  }
  //bugout(0, _T("n1 %d s1 %d e1 %d"), n1, s1, e1);
  //bugout(0, _T("n2 %d s2 %d e2 %d"), n2, s2, e2);

  if (missed > missingOk)
  {
    return 0;
  }

  // reorder and stuff back into original arrays

  // find point halfway between s2 and e2
  int mid = (s2+e2)/2;  // default
  double halfdist = 0.5*dist(sc2.m[s2], sc2.m[e2]);
  double mindiff = 1.0e20;
  for (i=s2 + 1; i<e2; i++)
  {
    double d = dist(sc2.m[i], sc2.m[s2]);
    double diff = fabs(d - halfdist);
    if (diff < mindiff)
    {
      mindiff = diff;
      mid = i;
    }
  }

  int j = 0;
  for (i=mid; i<=e2; i++)
  {
    x[j] = sc2.m[i][0];
    y[j] = sc2.m[i][1];
    z[j] = sc2.m[i][2];
    if (vi && vj)
    {
      vi[j] = sc2.m[i][3];
      vj[j] = sc2.m[i][4];
    }
    if (j && _hypot(x[j] - x[j-1], y[j] - y[j-1]) < 0.0002)
      continue;  // skip it
    //bugout(0, _T("7 %f %f %d trimmed"), x[j], y[j], j);
    j++;
  }
  *end1 = j-1;
  for (i=s1; i<=e1; i++)
  {
    x[j] = sc1.m[i][0];
    y[j] = sc1.m[i][1];
    z[j] = sc1.m[i][2];
    if (vi && vj)
    {
      vi[j] = sc1.m[i][3];
      vj[j] = sc1.m[i][4];
    }
    if (j && _hypot(x[j] - x[j-1], y[j] - y[j-1]) < 0.0002)
      continue;  // skip it
    //bugout(0, _T("7 %f %f %d trimmed"), x[j], y[j], j);
    j++;
  }
  *end2 = j-1;
  for (i=s2; i<mid; i++)
  {
    x[j] = sc2.m[i][0];
    y[j] = sc2.m[i][1];
    z[j] = sc2.m[i][2];
    if (vi && vj)
    {
      vi[j] = sc2.m[i][3];
      vj[j] = sc2.m[i][4];
    }
    if (j && _hypot(x[j] - x[j-1], y[j] - y[j-1]) < 0.0002)
      continue;  // skip it
    //bugout(0, _T("7 %f %f %d trimmed"), x[j], y[j], j);
    j++;
  }

  if (j && _hypot(x[0] - x[j-1], y[0] - y[j-1]) < 0.0002)
    j--;

  return j;
}


// used to trim overlapping scans of uncompensated points

int TrimEnd(int n1, CMatrix &sc1, int n2, CMatrix &sc2, int *s1, int *e2, int useInter, bool isEnglish)
{
  // find line segments intersect
  //bugout(0, _T("n1 %d sc1 %d %d n2 %d sc2 %d %d useInter %d"), n1, sc1.R(), sc1.C(), n2, sc2.R(), sc2.C(), useInter);

  int i,j=0;
  for (i=0; i<n1/2; i++)  // big hammer, should be quicker way
  {
    //bugout(0,_T("8 %f %f %f %f curve1 pt %d"), sc1.m[i][0], sc1.m[i][1], sc1.m[i+1][0], sc1.m[i+1][1], i);
    for (j=n2-2; j>n2/2; j--)
    {
      //if (i == 0) bugout(0, _T("9 %f %f %f %f curve2 pt %d"), sc2.m[j][0], sc2.m[j][1], sc2.m[j+1][0], sc2.m[j+1][1], j);
      if (intersect(sc1.m[i], sc1.m[i+1], sc2.m[j], sc2.m[j+1]))
        break;
    }

    if (j > n2/2)
    {
      //bugout(0, _T("8 %f %f %f %f INTERSECTION"), sc1.m[i][0], sc1.m[i][1], sc1.m[i+1][0], sc1.m[i+1][1]);
      //bugout(0, _T("9 %f %f %f %f INTERSECTION"), sc2.m[j][0], sc2.m[j][1], sc2.m[j+1][0], sc2.m[j+1][1]);
      //bugout(0, _T("would use %d %d"), i+1, j);
      break;
    }
  }


  if (i == n1/2)  // no intersections
  {
    //bugout(0, _T("No Intersection"));
    double d, v1[2], v2[2];
    d  = 0.25*(dist(sc1.m[0],sc1.m[1]) + dist(sc1.m[1],sc1.m[2]) +
      dist(sc2.m[n2-1],sc2.m[n2-2]) + dist(sc2.m[n2-2],sc2.m[n2-3]));

    // find first candidate for continuation of curve 1

    v1[0] = sc1.m[0][0] - sc1.m[1][0];
    v1[1] = sc1.m[0][1] - sc1.m[1][1];
    //normalize(v1,v1);
    for (j=n2-1; j>n2/2; j--)
    {
      v2[0] = sc2.m[j][0] - sc1.m[0][0];
      v2[1] = sc2.m[j][1] - sc1.m[0][1];
      //normalize(v2,v2);
      if (dot(v1,v2) > 0)   // point is on correct side
        break;
    }

    if (j > n2/2)  // found a candidate, is is close enough?
    {
      double dt = dist(sc1.m[0], sc2.m[j]);
      if (dt < 2.5*d)
      {
        *s1 = 0;
        if (dt < 0.3*d) // too close, use next point
          *e2 = --j;
        else            // use this point
          *e2 = j;

        return 1;
      }
    }

    // find first candidate for continuation of curve 2

    v1[0] = sc2.m[n2-1][0] - sc2.m[n2-2][0];
    v1[1] = sc2.m[n2-1][1] - sc2.m[n2-2][1];
    //normalize(v1,v1);
    for (i=0; i<n1/2; i++)  // big hammer, should be quicker way
    {
      v2[0] = sc1.m[i][0] - sc2.m[n2-1][0];
      v2[1] = sc1.m[i][1] - sc2.m[n2-1][1];
      //normalize(v2,v2);
      if (dot(v1,v2) > 0)   // point is on correct side
        break;
    }

    if (i < n1/2)  // found a candidate, is is close enough?
    {
      double dt = dist(sc2.m[n2-1], sc1.m[i]);
      if (dt < 2.5*d)
      {
        *e2 = n2-1;
        if (dt < 0.3*d) // too close, use next point
          *s1 = ++i;
        else            // use this point
          *s1 = i;

        return 1;
      }
    }

    // are the two curves parallel, look for close encounters.


    int bi=0,bj=0;
    double minD = 1.0e20;

    for (i=0; i<n1/2; i++)
    {
      double dd, lv[2];
      lv[0] = sc1.m[i+1][0] - sc1.m[i][0];
      lv[1] = sc1.m[i+1][1] - sc1.m[i][1];
      double length = normalize(lv, lv);

      for (j=n2-2; j>n2/2; j--)
      {
        dd = projdist(sc1.m[i], lv, sc2.m[j]);
        if (dd >= 0 && dd <= length)
        {
          dd = ptlinedist(sc1.m[i], lv, sc2.m[j]);
          if (dd < minD)
          {
            minD = dd;
            bi = i;
            bj = j;
          }
        }

        dd = projdist(sc1.m[i], lv, sc2.m[j+1]);
        if (dd >= 0 && dd <= length)
        {
          dd = ptlinedist(sc1.m[i], lv, sc2.m[j+1]);
          if (dd < minD)
          {
            minD = dd;
            bi = i;
            bj = j;
          }
        }
      }
    }

    double threshold = isEnglish ? 0.01 : 0.25;//0.005  0.125  dgj

    if (minD < threshold)
    {
      i = bi;
      j = bj;
    }
    else
    {
      //bugout(0, _T("Total Failure"));
      return 0;
    }
  }

  // don't use circle algoritm to better locate actual intersection if useInter is set or if nott enough points on both tails
  if (useInter || i > n1 - 5 || j < 3)
  {
    *s1 = i+1;
    *e2 = j;
    return 1;
  }

  // fit a circle to 6 points

  int k;
  double x[6], y[6];
  for (k=0; k<3; k++)
  {
    x[k] = sc1.m[i+2+k][0];
    y[k] = sc1.m[i+2+k][1];
    x[3+k] = sc2.m[j-1-k][0];
    y[3+k] = sc2.m[j-1-k][1];
    //bugout(3, "7 %f %f", x[k], y[k]);
    //bugout(3, "7 %f %f", x[3+k], y[3+k]);
  }

  CCircle cc(6, x, y);
  if (!cc.Valid())
  {
    *s1 = i+1;
    *e2 = j;
    return 1;
  }

  double r, ctr[2];
  ctr[0] = cc.X();
  ctr[1] = cc.Y();
  r = cc.Rad();

  double d1 = dist(sc1.m[i], sc1.m[i+1]); // distance between adjacent points
  double d2 = dist(sc2.m[j], sc2.m[j+1]);
  double dave = 0.5*(d1 + d2);
  double gap = dist(sc1.m[i+1], sc2.m[j]);

  int diff = j - i;         // this should stay constant later

  if (gap > 1.6*dave)       // gap too wide, include i or j+1
  {
    diff++;
  }
  else if (gap < 0.5*dave)  // gap too narow, include
    diff--;

  int newi = i-3;
  int testj = j-3;
  if (newi < 0)
  {
    testj -= newi;
    newi = 0;
  }

  // where does scan 1 get closer to circle
  while (newi < n1-1 && testj < n2-1)
  {
    double ss1 = fabs(r - dist(ctr, sc1.m[newi]));
    double ss2 = fabs(r - dist(ctr, sc2.m[testj]));
    if (ss1 < ss2)  // i is closer to being on circle
      break;
    newi++;
    testj++;
  }

  if (newi == n1 || testj == n2-1)
  {
    *s1 = i+1;
    *e2 = j;
    return 1;
  }

  int newj = j+4;
  int testi = i+4;
  if (newj > n2-1)
  {
    int m = newj - n2 + 1;
    testi -= m;
    newj = n2-1;
  }

  // where does scan 2 get closer to circle
  while (newj > 0 && testi > 0)
  {
    double ss1 = fabs(r - dist(ctr, sc1.m[testi]));
    double ss2 = fabs(r - dist(ctr, sc2.m[newj]));
    if (ss2 < ss1)  // j is closer to being on circle
      break;
    newj--;
    testi--;
  }

  if (testi == 0 || newj == 0)
  {
    *s1 = i+1;
    *e2 = j;
    return 1;
  }

  int newdiff = newj - newi;

  while (newdiff < diff)  // need to pick up points
  {
    // add newi or newj+1
    double ss1 = fabs(r - dist(ctr, sc1.m[newi]));
    double ss2 = fabs(r - dist(ctr, sc2.m[newj+1]));
    if (ss1 < ss2)  // i is closer to being on circle
      newi--;
    else
      newj++;
    newdiff++;
  }

  while (newdiff > diff)  // need to kill up points
  {
    // kill newi+1 or newj
    double ss1 = fabs(r - dist(ctr, sc1.m[newi+1]));
    double ss2 = fabs(r - dist(ctr, sc2.m[newj]));
    if (ss1 > ss2)  // i is further from being on circle
      newi++;
    else
      newj--;
    newdiff--;
  }

  *s1 = newi+1;
  *e2 = newj;
  return 1;
}


/*
int RemoveLoops(int n, double *x, double *y, double *z, double *vi, double *vj)
{
if (n < 20)
return n;

double sumAngle = 0;
double v1[2], v2[2];
v1[0] = x[1] - x[0];
v1[1] = y[1] - y[0];
normalize(v1, v1);
int i;
for (i=2; i<n; i++)
{
double angle;
if (i % 2 == 0)
{
v2[0] = x[i] - x[i-1];
v2[1] = y[i] - y[i-1];
normalize(v2, v2);
angle = asin(curl(v1, v2));
}
else
{
v1[0] = x[i] - x[i-1];
v1[1] = y[i] - y[i-1];
normalize(v1, v1);
angle = asin(curl(v2, v1));
}
sumAngle += angle;
if (sumAngle > 2.0*M_PI)
{
int xxx = 999;
}
}

return n;
}
*/

__time64_t ConvertDateToTime64(int yr, int mo, int da, int hr, int mn, int sc)
{
  struct tm th;
  th.tm_year = yr;
  th.tm_mon = mo;
  th.tm_mday = da;
  th.tm_hour = hr;
  th.tm_min = mn;
  th.tm_sec = sc;

  return _mktime64(&th);  
}
