#include "StdAfx.h"
#include "OdrSpline.h"
#include "LevMarq.h"

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

#define SQRT // sqrt // For centripetal parameterization

ODRSpline::ODRSpline(int np, double *px, double *py, double *pk, bool isEnglish, int cl, int tn,
                     double th)
{
  m_english = isEnglish;
  m_good = 0;

  m_npts = np;

  if (pk)
    m_ndim = 3;
  else
    m_ndim = 2;
  m_nctl = 0;
  m_closed = cl;
  m_tension = tn;
  m_thin_tol = th;

  m_pcoeffs = 0;
  m_pbeta = 0;
  m_scale_fac = 0.0;   // set this to 1.0 to disable scaling

  m_pknots = 0;
  m_pcx = 0;
  m_pcy = 0;
  m_pck = 0;
  m_pa = 0;
  m_pb = 0;
  m_pd = 0;
  m_puu = 0;
  m_pzz = 0;

  m_pxm = 0;

  if (m_npts < 4 || !px || !py)
    return;

  m_pxm = new CMatrix(m_npts, m_ndim);

  int i,j;
  for(i=j=0; i<m_npts; i++,j++)
  {
    // check that same point is entered twice or more in a row

    if (j > 1)
    {
      if (fabs(px[j] - px[j-1]) < 1.0e-6 && fabs(py[j] - py[j-1]) < 1.0e-6)
      {
        i--;     // reuse this slot;
        m_npts--;  // one less point
        continue;
      }
    }

    m_pxm->m[i][0] = px[j];
    m_pxm->m[i][1] = py[j];
    if (m_ndim == 3)
      m_pxm->m[i][2] = pk[j];
  }
  //bugout(5, _T("in ODRSpline after duplicate point check, m_npts=%d"), m_npts);

  if (m_npts < 4)
    return;

  transform_points(m_npts, m_pxm->m, m_X0);

  double *x = new double[m_npts + 2];
  double *y = new double[m_npts + 2];
  double *k = pk ? new double[m_npts + 2] : 0;

  if (m_closed)
  {
    x[0] = x[m_npts] = m_pxm->m[0][0];
    y[0] = y[m_npts] = m_pxm->m[0][1];
    if (pk)
      k[0] = k[m_npts] = m_pxm->m[0][2];
    for (i=1; i<m_npts; i++)
    {
      x[i] = m_pxm->m[i][0];
      y[i] = m_pxm->m[i][1];
      //bugout(5, _T("9 %f %f %f %f"), m_X0[0]+x[i-1], m_X0[1]+y[i-1], m_X0[0]+x[i], m_X0[1]+y[i]);
      if (pk)
        k[i] = m_pxm->m[i][2];
    }
  }
  else
  {
    for (i=0; i<m_npts; i++)
    {
      x[i+1] = m_pxm->m[i][0];
      y[i+1] = m_pxm->m[i][1];
      if (pk)
        k[i+1] = m_pxm->m[i][2];
    }
  }

  int *loc_parms = new int[m_npts+1];
  for (i = 0; i <= m_npts; i++)
    loc_parms[i] = 0;
  m_thin_tol /= m_scale_fac;
  if (m_closed)
    m_nctl = thin(m_npts+1, x, y, k, loc_parms) + 1;
  else
  {
    m_nctl = thin(m_npts, &x[1], &y[1], k ? &k[1] : 0, loc_parms) + 1;
    x[0] = x[1];
    y[0] = y[1];
    x[m_nctl] = x[m_nctl-1];
    y[m_nctl] = y[m_nctl-1];
    if (pk)
    {
      k[0] = k[1];
      k[m_nctl] = k[m_nctl-1];
    }
  }

  if (m_nctl < 4)
  {
    delete [] x;
    delete [] y;
    if (k)
      delete [] k;
    delete [] loc_parms;
    return;
  }

  m_pknots = new double[m_nctl+5];

  int ncoeffs, nshape;

  ncoeffs = 2*(m_nctl+1);

  if (m_closed)
    nshape = 2*(m_nctl-2);
  else
    nshape = ncoeffs;
  m_pbeta = new double[nshape];
  if (m_closed)
    m_pcoeffs = new double[2*ncoeffs];
  else
    m_pcoeffs = m_pbeta;

  m_pcx = &m_pcoeffs[0];
  m_pcy = &m_pcoeffs[m_nctl+1];

  if (pk)
    m_pck = new double[m_nctl+1];

  // Get initial values for the spline parameters by interpolation at
  // the knots.

  double c;

  if (m_closed)
  {
    if (!closed_curve(m_nctl, x, y, k, m_pknots+3, m_pcx, m_pcy, m_pck, 1))
    {
      delete [] x;
      delete [] y;
      if (k)
        delete [] k;
      delete [] loc_parms;
      return;
    }

    // Continue the knots periodically.
    c = m_pknots[m_nctl+1] - m_pknots[3];
    for (i=0; i<3; i++)
    {
      m_pknots[2-i] = m_pknots[m_nctl-i] - c;
      m_pknots[m_nctl+2+i] = m_pknots[4+i] + c;
    }
  }
  else if (!init_spline(m_nctl, x, y, k, m_pknots, m_pcx, m_pcy, m_pck))
  {
    delete [] x;
    delete [] y;
    if (k)
      delete [] k;
    delete [] loc_parms;
    return;
  }

  delete [] x;
  delete [] y;
  if (k)
    delete [] k;

  delete [] loc_parms;  // remove these if next line not WHAT_IF

  // need to unscale and shift back control points
  for (i = 0; i <= m_nctl; i++)
  {
    m_pcx[i] = m_scale_fac * m_pcx[i] + m_X0[0];
    m_pcy[i] = m_scale_fac * m_pcy[i] + m_X0[1];
  }

  //bugout(5, _T("odrspline finished"));

  m_good = 1;
}



ODRSpline::~ODRSpline()
{
  if (m_pknots)
    delete [] m_pknots;
  if (m_pa)
    delete [] m_pa;
  if (m_pb)
    delete [] m_pb;
  if (m_pd)
    delete [] m_pd;
  if (m_puu)
    delete [] m_puu;
  if (m_pzz)
    delete [] m_pzz;

  if (m_pck)  // don't delete m_pcx and m_pcy, they aren't allocated
    delete [] m_pck;

  if (m_pxm)
    delete m_pxm;

  if (m_closed && m_pcoeffs)
    delete [] m_pcoeffs;
  if (m_pbeta)
    delete [] m_pbeta;
}

// Transform a set of points.  The data are shifted to centroid coordinates.

int ODRSpline::transform_points(int n, double **x, double *X0)
{
  // Shift data to centroid coordinates...

  centroid(n, 2, x, X0);
  shift_all(n, 2, X0, x);

  /* Scale the data only if m_scale_fac is not 1. */
  int i,j;
  double data_max, temp;
  if (m_scale_fac != 1.)
  {
    for (data_max = 0., i = 0; i < n; i++)
      for (j = 0; j < 2; j++)
        if ((temp = fabs(x[i][j])) > data_max)
          data_max = temp;
    if (data_max > 2.)
    {
      m_scale_fac = data_max / 2.;
      for (i = 0; i < n; i++)
        for (j = 0; j < 2; j++)
          x[i][j] /= m_scale_fac;
    }
    else
      m_scale_fac = 1.;
  }

  return 1;
}

//  Remove knots from a piecewise linear spline without causing the resulting
//  spline to deviate from the original one by more than a given tolerance.
//  The measure of deviation is area per unit length.  Since this measure can
//  fail to detect narrow peaks, code is included to find peaks that are not
//  already chosen as knots.  For smooth curves where such peaks may be
//  considered to be noise, this code may be bypassed or deleted.  There is
//  also code for finding corners, i.e., knots where there is a sharp turn.

int ODRSpline::thin(int n, double *x, double *y, double *kv, int *loc_parms)
{
  int i, j, k, d;
  double L, Lmax, f, ff, df, D, dev;
  double nrm1, nrm2, t;
  double tmp, nrmN, P[3], P1[3];

  if (m_thin_tol <= 0.)   // don't thin
  {
    for (i = 0; i < n; i++)
      loc_parms[i] = i;
    return n;  // this bypasses checks for reversing points, etc.
  }

  f = 0;
  D = 0.;
  P[0] = x[1] - x[0];
  P[1] = y[1] - y[0];
  Lmax = l2norm(P,2);

  k = 0;  /* First point triplet is (0,1,2). */
  i = 2;
  d = 1;
  loc_parms[0] = 0;

  while (i < n)
  {
    df = sarea(x, y, k, i-1, i);
    nrmN = fabs(df);
    ff = f + df;

    P[0] = x[k] - x[i];
    P[1] = y[k] - y[i];
    L = l2norm(P,2);

    /* Compute relative weight for point i-1 */
    /* as a corner and update max.           */
    if (L)
      dev = nrmN / L;  /* possible alternative: nrmN */
    else
      dev = 0.;

    if (dev > D)
    {
      D = dev;
      d = i-1;
    }

    /* Check for peak points */
    if (L >= Lmax)
      Lmax = L;
    else  /* Possible peak */
    {
      P[0] = x[i-1] - x[i];
      P[1] = y[i-1] - y[i];
      tmp = l2norm(P,2);
      if (tmp > 0.1 * m_thin_tol)
      {
        k++;
        x[k] = x[i-1];
        y[k] = y[i-1];
        if (kv)
          kv[k] = kv[i-1];
        loc_parms[k] = i-1;
        f = 0.;
        D = 0.;
        i++;
        Lmax = tmp;
        continue;
      }
    }
    /* Check for knots */
    nrm1 = fabs(ff);
    nrm2 = fabs(df);
    if ((nrm1 <= m_thin_tol * L) || (nrm2 <= m_thin_tol * 0.001))  /* Exclude some peaks? */
    {
      /* Don't include point i-1 */
      i++;
      f = ff;
    }
    else   /* Include a point in the list of knots. The     */
    {      /* point included is the most recent local max   */
      /* in terms of deviation from the linear spline. */
      k++;
      x[k] = x[d];
      y[k] = y[d];
      if (kv)
        kv[k] = kv[d];
      loc_parms[k] = d;
      f = 0;
      D = 0.;
      d++;

      i = d+1; /* Start over at last knot. */
      if (i < n)
      {
        P[0] = x[i] - x[d];
        P[1] = y[i] - y[d];
        Lmax = l2norm(P,2);
      }
    }
  }
  /* Include the last point. */
  k++;
  x[k] = x[n-1];
  y[k] = y[n-1];
  if (kv)
    kv[k] = kv[n-1];
  loc_parms[k] = n-1;

  // Identify corners (these are deleted).

  i = 0;
  while (1)
  {
    i++;
    if (i >= k)
      break;

    P[0] = x[i+1] - x[i];
    P[1] = y[i+1] - y[i];
    P1[0] = x[i-1] - x[i];
    P1[1] = y[i-1] - y[i];
    nrm1 = l2norm(P,2);
    nrm2 = l2norm(P1,2);
    if (nrm1 * nrm2 != 0.0)
      t = dot(P,P1,2) / (nrm1 * nrm2); /* Cosine of angle. */
    else
      t = -1.0;
    if (t >= 0.)     /* Angle between 0 and pi/2. */
    {
      /* Delete knot at x[i] */
      for (j = i; j < k; j++)
      {
        x[j] = x[j+1];
        y[j] = y[j+1];
        if (kv)
          kv[j] = kv[j+1];
        loc_parms[j] = loc_parms[j+1];
      }
      k--;
    }
  }

  return k+1;
}

//  Compute the chord length or centripetal parametrization and set the knots
//  for a closed spline curve.
//  Input:  Data P[0],...,P[m-3] as arrays x[0:m-3] and y[0:m-3]
//  Output: knots t[0],...,t[m-2].

int ODRSpline::closed_knots(int m, double *x, double *y, double *t)
{
  int i;
  double P[3], temp;

  t[0] = 0.;
  for (i = 1; i <= m-3; i++)
  {
    P[0] = x[i] - x[i-1];
    P[1] = y[i] - y[i-1];
    t[i] = t[i-1] + SQRT(l2norm(P,2));
  }
  P[0] = x[0] - x[m-3];
  P[1] = y[0] - y[m-3];
  t[m-2] = temp = t[m-3] + SQRT(l2norm(P,2));

  if (!temp && m > 3) // All knot parameters are equal.
    return 0;

  if (temp < 1. && m > 3)
  {
    for (i = 0; i < m-2; i++)
      t[i] /= temp;
    t[m-2] = 1.;
  }

  return 1;
}

//  Set up the spline interpolation matrix (d is diagonal, a is above, b is
//  below) using periodic endpoint conditions.  The matrix is cyclic
//  tridiagonal.  The distinct knots are t[0],...,t[m-2]; they are implicitly
//  continued periodically to the left and right.

int ODRSpline::setcmatrix(int m, double *a, double *d, double *b, double *t)
{
  int i;
  double delim1, delim2, deli, delip1, del, del1, del2;

  if (m <= 3) // Not enough data to compute a closed curve.
    return 0;

  deli = t[1] - t[0];
  if (deli == 0.) // Parameters not distinct.
    return 0;

  delim1 = t[m-2] - t[m-3];
  delim2 = t[m-3] - t[m-4];

  for (i = 0; i < m-3; i++)
  {
    delip1 = t[i+2] - t[i+1];
    del = delim1 + deli;
    del1 = delim2 + delim1;
    del2 = deli + delip1;
    a[i] = delim1 * delim1 / (delim1 + del2) / del;
    b[i] = deli * deli / (del1 + deli) / del;
    d[i] = (deli * del1 / (del1 + deli) +
      delim1 * del2 / (delim1 + del2)) / del;
    delim2 = delim1;
    delim1 = deli;
    deli = delip1;
  }

  // i = m - 3
  delip1 = t[1] - t[0];
  del = delim1 + deli;
  del1 = delim2 + delim1;
  del2 = deli + delip1;
  a[m-3] = delim1 * delim1 / (delim1 + del2) / del;
  b[m-3] = deli * deli / (del1 + deli) / del;
  d[m-3] = (deli * del1 / (del1 + deli) +
    delim1 * del2 / (delim1 + del2)) / del;

  return 1;
}

//  Solve a cyclic tridiagonal linear system.
//  Make sure the spline variables have been allocated!

int ODRSpline::solve_ctrid(int m, double *a, double *d, double *b, double *r,
                           double *x)
{
  int i;
  double alpha, beta, gamma, fac;

  alpha = a[m-3];
  beta = b[0];

  // Set up the diagonal of the modified tridiagonal system.
  gamma = -d[0];
  d[0] -= gamma;
  d[m-3] -= alpha * beta / gamma;

  // Factor the modified system.
  for (i = 1; i < m-2; i++)
  {
    b[i] /= d[i-1];
    d[i] -= b[i] * a[i-1];
  }

  // Solve Ax = r.
  solvetrid(m-2, a, d, &b[1], r, x);

  // Set up m_puu.
  m_puu[0] = gamma;
  m_puu[m-3] = alpha;
  for (i = 1; i < m-3; i++)
    m_puu[i] = 0.;

  // Solve Azz = u.
  solvetrid(m-2, a, d, &b[1], m_puu, m_pzz);

  // Compute the solution vector.
  fac = (x[0] + beta * x[m-3] / gamma) /
    (1. + m_pzz[0] + beta * m_pzz[m-3] / gamma);
  for (i = 0; i < m-2; i++)
    x[i] -= fac * m_pzz[i];

  return 1;
}

//  Compute a 2D or 3D closed cubic spline curve interpolant If kts = 1 knots
//  are computed; otherwise they must be supplied.
//  Input:  Data P[0],...,P[m-3] as arrays x[0:m-3] and y[0:m-3].

int ODRSpline::closed_curve(int m, double *x, double *y, double *k, double *t,
                            double *cx, double *cy, double *ck, int kts)
{
  if (m < 3)
    return 0;

  if (m_pa)
    delete [] m_pa;
  m_pa = new double[m - 2];

  if (m_pb)
    delete [] m_pb;
  m_pb = new double[m - 2];

  if (m_pd)
    delete [] m_pd;
  m_pd = new double[m - 2];

  if (m_puu)
    delete [] m_puu;
  m_puu = new double[m - 2];

  if (m_pzz)
    delete [] m_pzz;
  m_pzz = new double[m - 2];

  // Compute chord length or centripetal parametrization.

  if (kts)
    if (!closed_knots(m, x, y, t))
      return 0;

  if (m > 3)
  {
    // Set up the cyclic tridiagonal system and solve for cx.
    if (!setcmatrix(m, m_pa, m_pd, m_pb, t))
      return 0;
    if (!solve_ctrid(m, m_pa, m_pd, m_pb, x, &cx[1]))
      return 0;

    // Set up the cyclic tridiagonal system and solve for cy.
    if (!setcmatrix(m, m_pa, m_pd, m_pb, t))
      return 0;
    if (!solve_ctrid(m, m_pa, m_pd, m_pb, y, &cy[1]))
      return 0;

    if (k && ck)
    {
      // Set up the cyclic tridiagonal system and solve for ck.
      if (!setcmatrix(m, m_pa, m_pd, m_pb, t))
        return 0;
      if (!solve_ctrid(m, m_pa, m_pd, m_pb, k, &ck[1]))
        return 0;
    }
  }
  else
  {
    cx[1] = x[m-3];
    cy[1] = y[m-3];
    if (k && ck)
      ck[1] = k[m-3];
  }

  // Continue the control points periodically.
  cx[0] = cx[m-2];
  cy[0] = cy[m-2];
  if (k && ck)
    ck[0] = ck[m-2];

  int i;
  for (i = 1; i < 3; i++)
  {
    cx[m-2+i] = cx[i];
    cy[m-2+i] = cy[i];
    if (k && ck)
      ck[m-2+i] = ck[i];
  }

  return 1;
}

//  Compute the chord length or centripetal parametrization and set the knots.
//  Following Boehm, we leave out the second and second to last data point in
//  computing the parametrization.
//  Input:  Data P[0],...,P[m] as arrays x[0:m] and y[0:m].
//  Output: knots t[0],...,t[m+4].

int ODRSpline::chord(int m, double *x, double *y, double *t)
{
  int i, j;
  double P[3], temp;

  // Set the first 4 knots.
  t[0] = t[1] = t[2] = t[3] = 0.;

  P[0] = x[2] - x[0];
  P[1] = y[2] - y[0];
  t[4] = SQRT(l2norm(P,2));

  for (i = 2; i < m-2; i++)
  {
    P[0] = x[i+1] - x[i];
    P[1] = y[i+1] - y[i];
    t[i+3] = t[i+2] + SQRT(l2norm(P,2));
  }
  P[0] = x[m] - x[m-2];
  P[1] = y[m] - y[m-2];
  t[m+1] = t[m] + SQRT(l2norm(P,2));

  if (t[m+1] == 0.) // All knot parameters are zero.
    return 0;

  if ((temp = t[m+1]) < 1.)
    for (j = 4; j <= m+1; j++)
      t[j] /= temp;

  // Set the last 4 knots.
  t[m+4] = t[m+3] = t[m+2] = t[m+1];

  return 1;
}

//  Set up the triadiagonal matrix (d is diagonal, a is above, b is below).

int ODRSpline::setmatrix(int L, double *a, double *d, double *b, double *t)
{
  int i;
  double delim1, delim2, deli, delip1, del, del1, del2;

  if (L < 1)
    return 0;

  d[0] = d[L] = 1.;
  a[0] = b[L-1] = 0.;
  if (L == 1)
    return 1;

  // L >= 2
  delim1 = 0.;
  deli = t[4] - t[3];
  delip1 = t[5] - t[4];
  if (deli == 0. || delip1 == 0.) // parameters not distinct.
    return 0;

  for (i = 1; i < L-1; i++)
  {
    delim2 = delim1;
    delim1 = deli;
    deli = delip1;
    delip1 = t[i+5] - t[i+4];
    del = delim1 + deli;
    del1 = delim2 + delim1;
    del2 = deli + delip1;
    d[i] = (deli * del1 / (del1 + deli) + delim1 * del2 /
      (delim1 + del2)) / del;
    a[i] = delim1 * delim1 / (delim1 + del2) / del;
    b[i-1] = deli * deli / (del1 + deli) / del;
  }

  // i = L - 1
  delim2 = delim1;
  delim1 = deli;
  deli = delip1;
  del = delim1 + deli;
  del1 = delim2 + delim1;
  d[L-1] = (deli * del1 / (del1 + deli) + delim1 * deli /
    (delim1 + deli)) / del;
  a[L-1] = delim1 * delim1 / (delim1 + deli) / del;
  b[L-2] = deli * deli / (del1 + deli) / del;

  return 1;
}

//  Set the right hand side of the system using Bessel conditions at the
//  endpoints.  Given the data in P[0],..,P[m], with P[1] and P[m-1] initially
//  empty, compute values for P[1] and P[m-1] so the resulting cubic B-spline
//  curve is C^2.  This routine needs to be called once for each coordinate.

int ODRSpline::setrhs(int m, double *t, double *P)
{
  if (m < 3)
    return 0;

  if (m == 3)
  {
    P[1] = (2. * P[0] + P[3]) / 3.;
    P[2] = (2. * P[3] + P[0]) / 3.;
    return 1;
  }

  double t1, t2, tmp;

  if (m == 4)
  {
    t1 = (t[5] - t[4]) / (t[5] - t[3]);
    t2 = 1. - t1;
    if (t2 == 0.)
    {
      tmp = .5 *(P[2] - t1 * t1 * P[0]) / t1;
      P[1] = (2. * (t1 * P[0] + tmp) + P[0]) / 3.;
    }
    else
    {
      tmp = .5 * (P[2] - t1 * t1 * P[0] - t2 * t2 * P[4]) / (t1 * t2);
      P[1] = (2. * (t1 * P[0] + t2 * tmp) + P[0]) / 3.;
    }
    t1 = (t[5] - t[4]) / (t[5] - t[3]);
    t2 = 1. - t1;
    tmp = .5 * (P[2] - t1 * t1 * P[0] - t2 * t2 * P[4]) / (t1 * t2);
    P[3] = (2. * (t1 * tmp + t2 * P[4]) + P[4]) / 3.;

    return 1;
  }

  // m > 4 case

  t1 = (t[5] - t[4]) / (t[5] - t[3]);
  t2 = 1. - t1;
  if (t2 == 0.)
  {
    tmp = .5 *(P[2] - t1 * t1 * P[0]) / t1;
    P[1] = (2. * (t1 * P[0] + tmp) + P[0]) / 3.;
  }
  else
  {
    tmp = .5 * (P[2] - t1 * t1 * P[0] - t2 * t2 * P[3]) / (t1 * t2);
    P[1] = (2. * (t1 * P[0] + t2 * tmp) + P[0]) / 3.;
  }
  t1 = (t[m+1] - t[m]) / (t[m+1] - t[m-1]);
  t2 = 1. - t1;
  tmp = .5 *(P[m-2] - t1*t1*P[m-3] - t2*t2*P[m]) / (t1 * t2);
  P[m-1] = (2. * (t1 * tmp + t2 * P[m]) + P[m]) / 3.;

  return 1;
}

//  Compute an initial 3d cubic spline curve.

int ODRSpline::init_spline(int m, double *x, double *y, double *k, double *t,
                           double *cx, double *cy, double *ck)
{
  if (m <= 2) // Not enough data to compute a spline curve.
    return 0;

  if (m_pa)
    delete [] m_pa;
  m_pa = new double[m-2];

  if (m_pb)
    delete [] m_pb;
  m_pb = new double[m-2];

  if (m_pd)
    delete [] m_pd;
  m_pd = new double[m-1];

  // Compute chord length or centripetal parametrization.
  if (!chord(m, x, y, t))
    return 0;

  // Set up the tridiagonal linear system and factor it.
  if (!setmatrix(m-2, m_pa, m_pd, m_pb, t))
    return 0;
  tridlu(m-1, m_pa, m_pd, m_pb);

  // Set up the right hand side for the x-coordinate and solve for cx.
  if (!setrhs(m, t, x))
    return 0;
  solvetrid(m-1, m_pa, m_pd, m_pb, &x[1], &cx[1]);

  // Set up the right hand side for the y-coordinate and solve for cy.
  if (!setrhs(m, t, y))
    return 0;
  solvetrid(m-1, m_pa, m_pd, m_pb, &y[1], &cy[1]);

  // Set the first two and last two control points.
  cx[0] = x[0];
  cy[0] = y[0];
  cx[m] = x[m];
  cy[m] = y[m];

  if (k && ck)
  {
    // Set up the right hand side for the k-coordinate and solve for ck.
    if (!setrhs(m, t, k))
      return 0;
    solvetrid(m-1, m_pa, m_pd, m_pb, &k[1], &ck[1]);

    ck[0] = k[0];
    ck[m] = k[m];
  }

  return 1;
}