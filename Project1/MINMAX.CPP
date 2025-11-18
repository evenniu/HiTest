#include "stdafx.h"

#include "minmax.h"
//#include "automath.h"

/*
int numpts;                // total number of inputs 
double cc;
double spin;
double transx;
double transy;
static double scale_fac;

double shrink;

IVEC mia1 = {0,NULL};      // misc. integer array 1 
DVEC wt = {0,NULL};        // weight tolerance       // Brian Shier
double maxtol;             // largest weight or sum of squared weights 
DVEC mda1 = {0,NULL};      // misc. double array 1 
DVEC mda2 = {0,NULL};      // misc. double array 2 
DVEC mda3 = {0,NULL};      // misc. double array 3 
DVEC mda4 = {0,NULL};      // misc. double array 4 
DVEC mda5 = {0,NULL};      // misc. double array 5 
DVEC mda6 = {0,NULL};      // misc. double array 6 
DVEC tottol = {0,NULL};    // total tolerance 
DVEC zone = {0,NULL};      // tolerance used in algorithms 
DVEC omega = {0,NULL};     // weight multipliers 
DMAT noms = {0,0,NULL};    // Nominals 
DMAT vals = {0,0,NULL};    // Measured values 
DMAT infs = {0,0,NULL};    // Solution 
DMAT uvw = {0,0,NULL};     // rectangular: dx,dy,dz 
DMAT rtp = {0,0,NULL};     // spherical: rho,theta,phi 
DMAT mdda1 = {0,0,NULL};   // misc double double array 1 
ALIGNREC cmmalign;

*/
double BIGTOL = 1.0e-12;

CMinMax::CMinMax(int np)
{
  numpts = np;
  tol = 1.0e-7;

  alloc_arrays();
}

CMinMax::~CMinMax()
{
  free_arrays();
}


void CMinMax::free_arrays()
{
  freeivector(&mia1);
  //freedvector(&wt);    // Brian Shier

  freedmatrix(&noms);
  freedmatrix(&vals);
  freedmatrix(&infs);
  freedmatrix(&uvw);
  freedmatrix(&rtp);
  freedmatrix(&mdda1);

  //freedvector(&tottol);
  freedvector(&zone);
  freedvector(&omega);

  freedvector(&mda1);
  freedvector(&mda2);
  freedvector(&mda3);
  freedvector(&mda4);
  freedvector(&mda5);
  freedvector(&mda6);

  int i;
  for (i = 0; i < 3; i++)
    delete [] cmmalign.mat[i];
}

int CMinMax::alloc_arrays()
{
  ivector(numpts,&mia1);
  //dvector(numpts,&wt);     // Brian Shier

  dmatrix(numpts,3,&noms);
  dmatrix(numpts,3,&vals);
  dmatrix(numpts,3,&infs);
  dmatrix(numpts,3,&uvw);
  dmatrix(numpts,3,&rtp);
  dmatrix(numpts,3,&mdda1);

  //dvector(numpts,&tottol);
  dvector(numpts,&zone);
  dvector(numpts,&omega);

  dvector(numpts,&mda1);
  dvector(numpts,&mda2);
  dvector(numpts,&mda3);
  dvector(numpts,&mda4);
  dvector(numpts,&mda5);
  dvector(numpts,&mda6);

  int i;
  for (i = 0; i < 3; i++)  // allocate alignment record rows
    cmmalign.mat[i] = new double[3];

  return(1);
}

void CMinMax::calc_align(double *dx,double *dy,double *dtheta)
{
  int i, m;
  double cxi, cyi, cxv, cyv, delta, rad;

  // calculate centroids 
  for (i = 0, cxi = cyi = cxv = cyv = 0.0; i < numpts; i++)
  {
    cxi += infs.m[i][0];
    cyi += infs.m[i][1];
    cxv += vals.m[i][0];
    cyv += vals.m[i][1];
  }
  cxi /= numpts;
  cyi /= numpts;
  cxv /= numpts;
  cyv /= numpts;

  // determine rotation angle using point which is farthest from centroid 
  // angle is from original alignment to best fit alignment 
  delta = 0.0;
  m = 0;
  for (i = 0; i < numpts; i++)
  {
    rad = _hypot(vals.m[i][0] - cxv,vals.m[i][1] - cyv);
    if (rad > delta)
    {
      delta = rad;
      m = i;
    }
  }
  if (delta == 0.0)
    *dtheta = 0.0;
  else
  {
    delta *= delta;
    *dtheta = atan2(CURL(vals.m[m][0] - cxv,vals.m[m][1] - cyv,
      infs.m[m][0] - cxi,infs.m[m][1] - cyi) / delta,
      ((vals.m[m][0] - cxv) * (infs.m[m][0] - cxi) +
      (vals.m[m][1] - cyv) * (infs.m[m][1] - cyi)) / delta);
  }

  double c = cos(*dtheta);
  double s = sin(*dtheta);
  double tcxv = c * cxv - s * cyv;  // Rotate the vector
  double tcyv = s * cxv + c * cyv;  // (cxv,cyv) by *dtheta.
  *dx = cxi - tcxv;
  *dy = cyi - tcyv;
}

// Scale the weights and set the scaling factor maxtol.
void CMinMax::scale_weights()
{
  int i;

  for (maxtol = 0.0, i = 0; i < numpts; i++)
  {
    if (omega.v[i] > maxtol)
      maxtol = omega.v[i];
  }
  if (maxtol == 0.0)
    maxtol = 1.0;

  if (maxtol != 1.0)
  {
    double sqrt_maxtol = sqrt(maxtol);
    for (i = 0; i < numpts; i++)
    {
      omega.v[i] /= maxtol;
      zone.v[i] *= sqrt_maxtol;
    }
  }
}

/*
void CMinMax::calc_omegas(int nUseTols)
{
int i;

for (i = 0; i < numpts; i++)
{
omega.v[i] = SQR(wt.v[i]);
}

// scale the weights and set maxtol
scale_weights();
}
*/

bool CMinMax::MinMaxFit(bool translateOnly/*=false*/)
{
  int i;

  if (numpts < 2)
    return false;


  // Weights should have already been calculated at this point!
  mmfit(translateOnly);         // generic min-max 

  calc_align(&transx,&transy,&spin);

  cmmalign.orig[0] = transx;
  cmmalign.orig[1] = transy;
  cmmalign.mat[0][0] = cos(spin);
  cmmalign.mat[1][1] = cos(spin);
  cmmalign.mat[0][1] = -sin(spin);
  cmmalign.mat[1][0] = sin(spin);
  cmmalign.mat[0][2] = cmmalign.mat[1][2] = 0.0;
  cmmalign.mat[2][0] = cmmalign.mat[2][1] = 0.0;
  cmmalign.mat[2][2] = 1.0;

  for (i = 0; i < numpts; i++)
  {
    vals.m[i][2] = 0.0;
    infs.m[i][2] = 0.0;
    noms.m[i][2] = 0.0;
  }

  for (i = 0; i < numpts; i++)
  {
    uvw.m[i][0] = noms.m[i][0] - infs.m[i][0];
    uvw.m[i][1] = noms.m[i][1] - infs.m[i][1];
    uvw.m[i][2] = 0.0;
  }

  /* provide method to pass back alignment

  // fill in matrix to pass alignment back
  if (pdMat)
  {
  for (i = 0; i < 4; i++)
  {
  for (int j = 0; j < 3; j++ )
  {
  if (!i)
  pdMat[i][j] = -cmmalign.orig[j];
  else
  pdMat[i][j] = cmmalign.mat[i-1][j];
  }
  }
  }

  */

  return true;
}

/* min max fit for rotation only (about pivot point) */

void CMinMax::mmfit(bool translateOnly)
{
  int i;
  double d, maxdev = 0.0, tmp;

  //lsfit(); done in blade

  for (i = 0; i < numpts;i++)
  {
    uvw.m[i][0] = noms.m[i][0] - infs.m[i][0];
    uvw.m[i][1] = noms.m[i][1] - infs.m[i][1];
    if ((tmp = _hypot(uvw.m[i][0],uvw.m[i][1]) / zone.v[i]) > maxdev)
      maxdev = tmp;
  }

  maxdev *= maxtol;

  cc = 0.0;
  for (i = 0; i < 40; i++)
  {
    d = trans(&infs,&uvw);
    if (translateOnly)  // no rotation, so jump out after first call to trans
      break;

    d = rotop(&infs,&uvw);

    if (fabs(d-maxdev) < tol)
      break;
    maxdev = d;
  }
}

void CMinMax::TransferIn(CMatrix *bnoms, CMatrix *bvals, CVector *bweights, CAlignment *balign)
{
  int i,j;
  for (i=0; i< numpts; i++)
  {
    noms.m[i][0] = bnoms->m[i][0];
    noms.m[i][1] = bnoms->m[i][1];
    vals.m[i][0] = infs.m[i][0] = bvals->m[i][0];
    vals.m[i][1] = infs.m[i][1] = bvals->m[i][1];

    omega.v[i] = bweights->v[i]*bweights->v[i];
    if (fabs(bweights->v[i]) < 1.0e-5)
      zone.v[i] = 10000.0;
    else
      zone.v[i] = 1.0/bweights->v[i];

    //bugout(0, _T("TransferIn %d  N %f %f M %f %f V %f %f W %f"), i,m_noms.m[i][0],m_noms.m[i][1],m_vals.m[i][0],m_vals.m[i][1],m_ijks.m[i][0],m_ijks.m[i][1],m_omega.v[i]);
  }

  // scale the weights and set maxtol
  scale_weights();

  for (i=0; i<2; i++)
    for (j=0; j<2; j++)
    {
      cmmalign.mat[i][j] = balign->m_mat[i][j];
      //bugout(0, _T(" mat[%d][%d] = %f"), i, j, m_cmmalign.mat[i][j]);
    }

    cmmalign.orig[0] = balign->m_borig[0];
    cmmalign.orig[1] = balign->m_borig[1];
    //bugout(0, _T(" orig= %f %f"), m_cmmalign.orig[0], m_cmmalign.orig[1]);
}

void CMinMax::TransferOut(CMatrix *binfs, CAlignment *balign)
{
  int i,j;

  for (i=0; i< numpts; i++)
  {
    binfs->m[i][0] = infs.m[i][0];
    binfs->m[i][1] = infs.m[i][1];
  }

  balign->m_borig[0] = cmmalign.orig[0];
  balign->m_borig[1] = cmmalign.orig[1];
  balign->m_morig[0] = 0.0;
  balign->m_morig[1] = 0.0;

  for (i=0; i<2; i++)
    for (j=0; j<2; j++)
      balign->m_mat[i][j] = cmmalign.mat[i][j];
}



void CMinMax::afrotate(double theta, DMAT* xyo, DMAT* xyi)
{
  int i;
  double s, c, xc, ys, temp;

  /* this may look a little funny, but each rotation uses 3 multiplications
  rather than 4. */

  s = sin(theta);
  c = cos(theta);
  for(i = 0; i < numpts; i++)
  {
    xc = xyi->m[i][0] * c;
    ys = xyi->m[i][1] * s;
    temp = xc - ys;
    xyo->m[i][1] = (xyi->m[i][0] + xyi->m[i][1]) * (c + s) - xc - ys;
    xyo->m[i][0] = temp;
  }
}

double CMinMax::SQR(double val) /* was a macro but caused problems in XENIX */
{
  return (val * val);
}
/*--------------------------------------------------------------------*/
/* Dynamically allocate a matrix as an array of pointers to pointers  */
/* to rows.                                                           */
/*--------------------------------------------------------------------*/

void CMinMax::dmatrix(int row, int col, DMAT* a)
{
  int i;

  if(row <= a->dim1 && col <= a->dim2) /* already large enough */
  {
    /* calling function may be expecting matrix to be zeroed out */
    for(i = 0; i < row; i++)
      memset((LPTSTR)a->m[i], 0, col * sizeof(double));
    return;
  }

  if(a->dim1)       /* malloced before */
    freedmatrix(a); /* kill smaller one */

  a->dim1 = row;
  a->dim2 = col;

  a->m = new double*[row];
  //VERIFY(a->m);

  for(i = 0; i < row; i++)
  {
    a->m[i] = new double[col];
    //VERIFY(a->m[i]);
  }
}

/*--------------------------------------------------------------------*/
/* Dynamically allocate a double vector.                              */
/*--------------------------------------------------------------------*/

void CMinMax::dvector(int n, DVEC* a)
{
  if(n <= a->dim) /* already large enough */
  {
    memset((LPTSTR)a->v, 0, n * sizeof(double));
    return;
  }

  if(a->dim)        /* malloced before */
    freedvector(a); /* kill smaller one */

  a->dim = n;
  a->v = new double[n];

  //VERIFY(a->v);
}

/*--------------------------------------------------------------------*/
/* Dynamically allocate an integer vector.                            */
/*--------------------------------------------------------------------*/
void CMinMax::ivector(int n, IVEC* a)
{
  if(n <= a->dim) /* already large enough */
  {
    memset((LPTSTR)a->v, 0, n * sizeof(int));
    return;
  }

  if(a->dim)        /* malloced before */
    freeivector(a); /* kill smaller one */

  a->dim = n;
  a->v = new int[n];

  //VERIFY(a->v);
}

/*--------------------------------------------------------------------*/
/* Free an allocated matrix of doubles                                */
/*--------------------------------------------------------------------*/

void CMinMax::freedmatrix(DMAT* a)
{
  int i;

  if(a == NULL || a->m == NULL)
    return;

  for(i = 0; i < a->dim1; i++)
    delete[] a->m[i];
  delete[] a->m;

  a->dim1 = a->dim2 = 0;
  a->m = NULL;
}

/*--------------------------------------------------------------------*/
/* Free an allocated vector of doubles                                */
/*--------------------------------------------------------------------*/

void CMinMax::freedvector(DVEC* a)
{
  if(a == NULL || a->v == NULL)
    return;

  delete[] a->v;

  a->dim = 0;
  a->v = NULL;
}

/*--------------------------------------------------------------------*/
/* Free an allocated vector of ints                                   */
/*--------------------------------------------------------------------*/

void CMinMax::freeivector(IVEC* a)
{
  if(a == NULL || a->v == NULL)
    return;

  delete[] a->v;

  a->dim = 0;
  a->v = NULL;
}

double CMinMax::ivmax(int* imax, double* r)
{
  double ri, rmax;
  int i;

  for(*imax = 0, rmax = r[0], i = 1; i < numpts; i++)
    if((ri = r[i]) > rmax)
    {
      rmax = ri;
      *imax = i;
    }
  return (rmax);
}

double CMinMax::ivmin(int* imin, double* r)
{
  double ri, rmin;
  int i;

  for(*imin = 0, rmin = r[0], i = 1; i < numpts; i++)
    if((ri = r[i]) < rmin)
    {
      rmin = ri;
      *imin = i;
    }
  return (rmin);
}

double CMinMax::dsign(double a)
{
  if(a > 0.0)
    return (1.0);
  else if(a < 0.0)
    return (-1.0);
  else
    return (0.0);
}

/*---------------------------------------------------------------*/
/*  Copy matrix A into matrix B.                                 */
/*  Use Mcopy if they are allocated in contiguous memory.        */
/*  Do not use Mcopy to copy a submatrix if the submatrix has a  */
/*  different number of columns as the original matrix!          */
/*---------------------------------------------------------------*/

void CMinMax::mcopy(int rows, int cols, double** A, double** B)
{
  int i;

  for(i = 0; i < rows; i++)
    memcpy(B[i], A[i], cols * sizeof(double));
}

// start of rotop.cpp code

/*********************************************************************/
/*  "setmxmn"  Find angles associated with max and min coordinates.  */
/*      If more than one, choose smallest such angle.                */
/*********************************************************************/

void CMinMax::setmxmn(double* amax, double* amin, int* imax, int* imin, DMAT* xy, DMAT* zw, double* r)
{
  double a, dmax, dmin, xmax, ymax, zmax, wmax, denom;
  double amx, amn, xmin, ymin, zmin, wmin;
  int i;

  dmax = ivmax(imax, r);
  dmin = ivmin(imin, r);

  *amax = M_PI;
  *amin = -M_PI;

  for(i = 0; i < numpts; i++)
  {
    if(fabs(r[i] - dmax) * maxtol < tol)
    {
      xmax = xy->m[i][0];
      ymax = xy->m[i][1];
      zmax = zw->m[i][0];
      wmax = zw->m[i][1];

      denom = _hypot(xmax, ymax) * _hypot(zmax, wmax);
      if(denom < BIGTOL)
      {
        if(fabs(dmax) > fabs(dmin))
          *amax = 0.0;
        else
          *amax = .5 * dsign(dmax) * M_PI;
        *imax = i;
        break;
      }

      a = (xmax * zmax + ymax * wmax) / denom;
      if(fabs(a) > 1.0)
        a = dsign(a);

      amx = ACOS(a);
      if(xmax * wmax - ymax * zmax < 0.0)
        amx = -amx;

      if(fabs(amx) < fabs(*amax))
      {
        *amax = amx;
        *imax = i;
      }
    }
  }

  for(i = 0; i < numpts; i++)
  {
    if(fabs(r[i] - dmin) * maxtol < tol)
    {
      xmin = xy->m[i][0];
      ymin = xy->m[i][1];
      zmin = zw->m[i][0];
      wmin = zw->m[i][1];
      denom = _hypot(xmin, ymin) * _hypot(zmin, wmin);
      if(denom < BIGTOL)
      {
        if(fabs(dmin) > fabs(dmax))
          *amin = 0.0;
        else
          *amin = .5 * dsign(dmin) * M_PI;
        *imin = i;
        break;
      }
      a = (xmin * zmin + ymin * wmin) / denom;
      if(fabs(a) > 1.0)
        a = dsign(a);
      amn = ACOS(a);
      if(xmin * wmin - ymin * zmin < 0.0)
        amn = -amn;
      if(fabs(amn) < fabs(*amin))
      {
        *amin = amn;
        *imin = i;
      }
    }
  }
}

/**************************************************************************/
/*  "rotop":                                                              */
/*  Rotate until max distance between points is minimized.                */
/**************************************************************************/

double CMinMax::rotop(DMAT* xy, DMAT* uv)
{
  int its = 0, i, j, j1, j2, imin, imax;
  double theta_max, theta_min, theta, theta1, x0, y0, fac = 0.5;
  double dmax, dmin, dm, cm;
  double x1, y1, x2, y2, w1, z1, w2, z2, c1, s1, rotx21, rotx22;

  rotpt(xy, uv, &x0, &y0); /* get good rotation point */

  for(i = 0; i < numpts; i++) /* shift origin */
  {
    xy->m[i][0] -= x0;
    xy->m[i][1] -= y0;
    noms.m[i][0] -= x0;
    noms.m[i][1] -= y0;
  }

  cmdists(mda1.v, xy);      /* compute distances and signs */
  dmax = ivmax(&j, mda1.v); /* compute max and min         */
  dmin = ivmin(&j, mda1.v);

  while(fabs(dmax + dmin) * maxtol >= tol)
  {
    dm = MAX(fabs(dmax), fabs(dmin)); /* maximal deviation */

    if(its++ > 20)
      break;

    setmxmn(&theta_max, &theta_min, &imax, &imin, xy, &noms, mda1.v);

    if(((fabs(theta_max) < tol * maxtol) && (fabs(dmax) == dm)) ||
       ((fabs(theta_min) < tol * maxtol) && (fabs(dmin) == dm)))
      break;

    if(fabs(theta_max - theta_min) < BIGTOL)
    {
      afrotate(theta_max, xy, xy);
      break;
    }

    if(dsign(dmax) + dsign(dmin) != 0.0)
    {
      if(dmin >= 0.0)
        afrotate(MIN(theta_max, theta_min) + BIGTOL, xy, xy);
      else
        afrotate(MAX(theta_max, theta_min) - BIGTOL, xy, xy);
      cmdists(mda1.v, xy);      /* compute distances and signs */
      dmax = ivmax(&j, mda1.v); /* compute max and min         */
      dmin = ivmin(&j, mda1.v);
    }
    else
    {
      if(dmax + dmin >= 0.0) /* rotate in positive direction */
      {
        j1 = imax;
        j2 = imin;
        theta1 = theta_max - .1 * BIGTOL;
      }
      else /* rotate in negative direction */
      {
        j1 = imin;
        j2 = imax;
        theta1 = theta_min + .1 * BIGTOL;
      }

      x1 = xy->m[j1][0];
      y1 = xy->m[j1][1];
      x2 = xy->m[j2][0];
      y2 = xy->m[j2][1];
      w1 = noms.m[j1][1];
      z1 = noms.m[j1][0];
      w2 = noms.m[j2][1];
      z2 = noms.m[j2][0];

      c1 = cos(theta1);
      s1 = sin(theta1);

      rotx21 = c1 * x2 - s1 * y2;
      rotx22 = s1 * x2 + c1 * y2;

      /* theta is the maximal allowable rotation */
      if(fabs(_hypot(x1, y1) - _hypot(z1, w1)) / zone.v[j1] < _hypot(z2 - rotx21, w2 - rotx22) / zone.v[j2])
        theta = brent(0.0, theta1, x1, y1, x2, y2, z1, w1, z2, w2, zone.v[j1], zone.v[j2]);
      else
        theta = theta1;

      cm = dm;

      /* we've found a maximal rotation theta; now rotate */
      /* so that all maximal vectors are reduced.           */
      for(i = 0; i < 30; i++)
      {
        afrotate(theta, &mdda1, xy); /* try rotating by theta */
        cmdists(mda1.v, &mdda1);

        dmax = ivmax(&j, mda1.v);
        dmin = ivmin(&j, mda1.v);

        dm = MAX(fabs(dmax), fabs(dmin));
        if((dmax * dmin >= 0.0) || (dm >= cm)) /* did it work?      */
          theta = fac * theta;                 /* no, reduce theta */
        else
          break; /* yes, quit loop   */
      }          /* end for i */

      mcopy(numpts, 3, mdda1.m, xy->m);
    }
  }

  for(i = 0; i < numpts; i++) /* shift back */
  {
    noms.m[i][0] += x0;
    noms.m[i][1] += y0;
    xy->m[i][0] += x0;
    xy->m[i][1] += y0;
    uv->m[i][0] = noms.m[i][0] - xy->m[i][0];
    uv->m[i][1] = noms.m[i][1] - xy->m[i][1];
  }

  return (maxtol * MAX(fabs(dmax), fabs(dmin)));
}

// end of rotop.cpp code

// start of trans.cpp code

/*****************************************************************************/
/*                                                                           */
/* "fmaxv"  find the outermost vectors in the active constraint set.         */
/*                                                                           */
/*****************************************************************************/

double CMinMax::fmaxv(int* imin, int* imax, DMAT* uv, int j, int* I)
{
  int i;
  double ui, vi, uj, vj, amax, amin, argu, a;

  amax = 0.0;
  amin = 0.0;
  *imax = j;
  *imin = j;
  uj = uv->m[j][0];
  vj = uv->m[j][1];
  for(i = 0; i < numpts; i++)
  {
    if((I[i]) && (i != j))
    {
      ui = uv->m[i][0];
      vi = uv->m[i][1];
      argu = (ui * uj + vi * vj) / (_hypot(ui, vi) * _hypot(uj, vj));
      if(fabs(argu) > 1.0)
        argu = dsign(argu);
      a = ACOS(argu);
      if(uj * vi - vj * ui < 0)
        a = -a;
      if(a < amin - .5 * BIGTOL)
      {
        amin = a;
        *imin = i;
      }
      if(a > amax + .5 * BIGTOL)
      {
        amax = a;
        *imax = i;
      }
    }
  }
  return (amax - amin);
}

/*****************************************************************************/
/*                                                                           */
/* "findmin":  minimize the maximum deviation under translation in           */
/*             direction (v1,v2).  Returns t.                                */
/*****************************************************************************/

double CMinMax::findmin(DMAT* uv, double v1, double v2)
{
  int i, k, kwillbe, m, *I;
  double vdk, vdkw, wk, wkw = 0.0, *vdot, *r;
  double t, tmaybe, t1, t2, d;
  double A, B, C, Q, disc;

  I = mia1.v;
  vdot = mda1.v;
  r = mda2.v;

  abdists(r, &k, I, &m, uv); /* k is first index with r(k) = dmax. */

  d = SQR2(v1, v2);

  for(i = 0; i < numpts; i++) /* -2*vdot(i) = p'_i(0) */
    vdot[i] = omega.v[i] * (v1 * uv->m[i][0] + v2 * uv->m[i][1]);

  vdk = vdot[k];                  /* If more than one max vector, choose */
  for(i = k + 1; i < numpts; i++) /* the one with the largest derivative */
    if((I[i]) && (vdot[i] < vdk))
    {
      k = i;
      vdk = vdot[k];
    }
  wk = omega.v[k];
  t = 0.0;

  while(t < MIN(1.0, vdk / (wk * d)) - BIGTOL)
  {
    tmaybe = MIN(1.0, vdk / (wk * d));
    kwillbe = k;
    vdkw = vdk;

    for(i = 0; i < numpts; i++)
    {
      /* intersection of p_k and p_i */
      A = (omega.v[i] - wk) * d;
      B = 2.0 * (vdk - vdot[i]);
      C = SQR(r[i]) - SQR(r[k]);
      disc = B * B - 4.0 * A * C;
      if(disc >= 0.0)
      {
        t1 = t2 = t;
        Q = -.5 * (B + dsign(B) * sqrt(disc));
        if(Q != 0.0)
        {
          t1 = C / Q;
          if(A != 0.0)
            t2 = Q / A;
        }

        if((omega.v[i] - wk) * d * t1 > vdot[i] - vdk)
        {
          if((t1 > t + BIGTOL) && (t1 < tmaybe - BIGTOL))
          {
            tmaybe = t1;
            kwillbe = i;
            vdkw = vdot[i];
            wkw = omega.v[k];
          }
          if((t1 == tmaybe) && ((omega.v[i] - wkw) * d * t1 > vdot[i] - vdkw))
          {
            kwillbe = i;
            vdkw = vdot[i];
            wkw = omega.v[i];
          }
        }
        else if((omega.v[i] - wk) * d * t2 > vdot[i] - vdk)
        {
          if((t2 > t + BIGTOL) && (t2 < tmaybe - BIGTOL))
          {
            tmaybe = t2;
            kwillbe = i;
            vdkw = vdot[i];
            wkw = omega.v[k];
          }
          if((t2 == tmaybe) && ((omega.v[i] - wkw) * d * t2 > vdot[i] - vdkw))
          {
            kwillbe = i;
            vdkw = vdot[i];
            wkw = omega.v[i];
          }
        }
      } /* end if disc */
    }   /* end for i */

    t = tmaybe;
    k = kwillbe;
    vdk = vdot[k];
    wk = omega.v[k];
  }

  return (t);
}

/*****************************************************************************/
/*                                                                           */
/* "trans":                                                                  */
/*     Translate until max distance between points is minimized.             */
/*     Returns new points and maximal distance.                              */
/*****************************************************************************/

double CMinMax::trans(DMAT* xy, DMAT* uv)
{
  int i, imin, imax, j, m, *I;
  double t, v1, v2, sqd = 0.0, spread, *r;

  I = mia1.v;
  r = mda1.v;

  bool bUntilDone = true;
  while(bUntilDone)
  { /*  first check conditions for completion */
    sqd = abdists(r, &j, I, &m, uv);
    if(sqd < BIGTOL)
      break;

    spread = fmaxv(&imin, &imax, uv, j, I);

    if(spread > M_PI - BIGTOL)
      break;

    v1 = (zone.v[imin] * uv->m[imax][0] + zone.v[imax] * uv->m[imin][0]) / (zone.v[imin] + zone.v[imax]);
    v2 = (zone.v[imin] * uv->m[imax][1] + zone.v[imax] * uv->m[imin][1]) / (zone.v[imin] + zone.v[imax]);

    if(_hypot(v1, v2) < BIGTOL)
      break;

    t = findmin(uv, v1, v2);

    if(t == 0.0)
      break;

    for(i = 0; i < numpts; i++)
    {
      uv->m[i][0] -= (t * v1);
      uv->m[i][1] -= (t * v2);
    }
  }

  for(i = 0; i < numpts; i++)
  {
    xy->m[i][0] = noms.m[i][0] - uv->m[i][0];
    xy->m[i][1] = noms.m[i][1] - uv->m[i][1];
  }

  return (maxtol * sqd);
}

// end of trans.cpp code

// start of rotpt.cpp code

static double xmid, ymid;

double CMinMax::tsign(double a)
{
  if(a > BIGTOL)
    return (1.0);
  else if(a < -BIGTOL)
    return (-1.0);
  else
    return (0.0);
}

/*********************************************************************/
/*  "intrsct": Find intersection point of two lines.                 */
/*             Intersection point is (xi,yi)+ti*(ui,vi), i=1,2.      */
/*********************************************************************/

int CMinMax::intrsct(double* x0, double* y0, double* l1, double* l2)
{
  double det1, t1;

  *x0 = 0.0;
  *y0 = 0.0;

  det1 = l1[3] * l2[2] - l1[2] * l2[3];

  if(fabs(det1) < tol)
    return (0); /* no intersection point found */

  t1 = (l2[2] * (l2[1] - l1[1]) - l2[3] * (l2[0] - l1[0])) / det1;
  *x0 = l1[0] + t1 * l1[2];
  *y0 = l1[1] + t1 * l1[3];

  if(_hypot(*x0 - xmid, *y0 - ymid) > 2.5 * cc)
    return (0);
  else
    return (1);
}

double CMinMax::fmaxu(int m, double* au, double* av, int j)
{
  int i, ii;
  double ui, vi, uj, vj, amax, amin, argu, a;

  amax = 0.0;
  amin = 0.0;
  uj = au[j];
  vj = av[j];
  for(i = 0; i < 2; i++)
  {
    ii = (j + i + 1) % m;
    ui = au[ii];
    vi = av[ii];
    argu = (ui * uj + vi * vj) / (_hypot(ui, vi) * _hypot(uj, vj));
    if(fabs(argu) > 1.0)
      argu = dsign(argu);
    a = ACOS(argu);
    if(uj * vi - vj * ui < 0)
      a = -a;
    if(a < amin - .5 * BIGTOL)
      amin = a;
    if(a > amax + .5 * BIGTOL)
      amax = a;
  }
  return (amax - amin);
}

/****************************************************************************/
/*  "rotpt":                                                                */
/*   Find a point about which rotation reduces all maximal vectors.         */
/****************************************************************************/

void CMinMax::rotpt(DMAT* xy, DMAT* uv, double* x0, double* y0)
{
  int i, j, k, m, *I;
  int mm, jj, kk, p, flag, n_flags, flags[10], sflag[10];
  double *au, *av, *ax, *ay, *ae, *aw;
  double x1 = 0.0, x11, y1 = 0.0, y11, fac, temp, sqd, ebar;
  double d, ds, xp, yp, xcm, ycm, denom;
  double xt[10], yt[10], xpt[10], ypt[10], c[5];
  double l[5][5], lp[5][5], ltemp[5], sign[10];

  I = mia1.v;
  au = mda1.v;
  av = mda2.v;
  ax = mda3.v;
  ay = mda4.v;
  ae = mda5.v;
  aw = mda6.v;

  if(cc == 0.0) /* only do this once per dataset! */
  {
    for(xmid = ymid = 0.0, i = 0; i < numpts; i++)
    { /* compute centroid of nominal points */
      xmid += noms.m[i][0];
      ymid += noms.m[i][1];
    }
    xmid /= numpts;
    ymid /= numpts;

    for(i = 0; i < numpts; i++) /* compute typical scale factor */
    {
      temp = MAX(fabs(noms.m[i][0] - xmid), fabs(noms.m[i][1] - ymid));
      if(temp > cc)
        cc = temp;
    }
  }

  sqd = abdists(au, &k, I, &m, uv);

  /* (more liberal estimate of the maximal vectors) */

  for(m = i = 0; i < numpts; i++)
    m += (I[i] = (fabs(au[i] - sqd) * maxtol < 10.0 * tol));

  for(i = j = 0; i < numpts; i++)
  {
    if(I[i])
    {
      ae[j] = zone.v[i];
      aw[j] = omega.v[i];
      au[j] = uv->m[i][0];
      av[j] = uv->m[i][1];
      ax[j] = xy->m[i][0];
      ay[j] = xy->m[i][1];
      j++;
    }
  }

  /* default: compute a weighted centroid of the measured points  */
  /* corresponding to maximal vectors (optimal for m=2).   */
  for(ebar = (*x0) = (*y0) = 0.0, i = 0; i < m; i++)
  {
    *x0 += ax[i] / ae[i];
    *y0 += ay[i] / ae[i];
    ebar += 1.0 / ae[i];
  }
  *x0 /= ebar;
  *y0 /= ebar;

  if((m < 3) || (m > 5))
    return;

  for(i = 0; i < m; i++)
  {
    au[i] *= aw[i];
    av[i] *= aw[i];
  }

  if(m == 3) /* compute rotation point that maximizes the */
  {          /* magnitude of the weighted minimum curl.   */
    denom = (au[1] - au[0]) * (av[2] - av[1]) - (av[1] - av[0]) * (au[2] - au[1]);
    if(fabs(denom) >= .2 * SQR(sqd))
    {
      if(fabs(denom) > 1.e-8)
      {
        for(i = 0; i < m; i++)
          c[i] = ax[i] * av[i] - ay[i] * au[i];
        *x0 = ((c[2] - c[1]) * (au[1] - au[0]) - (c[1] - c[0]) * (au[2] - au[1])) / denom;
        *y0 = ((c[2] - c[1]) * (av[1] - av[0]) - (c[1] - c[0]) * (av[2] - av[1])) / denom;
      }
      else
        *x0 = *y0 = 0;

      return;
    }
  }

  if(m == 4) /* do same for some subset of three maximal vectors */
  {
    for(i = 0; i < m; i++)
      c[i] = (ax[i] * av[i]) - (ay[i] * au[i]);
    ax = mda3.v;

    for(j = 0; j < m - 1; j++) /* skip one */
    {
      if(fmaxu(m, au, av, j) > M_PI)
      {
        denom =
            (au[j + 1] - au[j]) * (av[(j + 2) % m] - av[j + 1]) - (av[j + 1] - av[j]) * (au[(j + 2) % m] - au[j + 1]);
        if(fabs(denom) >= .2 * SQR(sqd))
        {
          if(fabs(denom) > 1.e-8)
          {
            x1 = ((c[(j + 2) % m] - c[j + 1]) * (au[j + 1] - au[j]) -
                  (c[j + 1] - c[j]) * (au[(j + 2) % m] - au[j + 1])) /
                 denom;
            y1 = ((c[(j + 2) % m] - c[j + 1]) * (av[j + 1] - av[j]) -
                  (c[j + 1] - c[j]) * (av[(j + 2) % m] - av[j + 1])) /
                 denom;
          }
          else
            x1 = y1 = 0;
        }
        d = CURL(au[j], av[j], ax[j] - x1, ay[j] - y1); /* next 2 are equal */
        ds = tsign(d);
        d = fabs(d);
        k = (j + 3) % m;
        if((tsign((temp = CURL(au[k], av[k], ax[k] - x1, ay[k] - y1))) == ds) && (fabs(temp) >= d)) /* optimal point
                                                                                                       found */
        {
          *x0 = x1;
          *y0 = y1;
          return;
        }
      }
    }
  }

  for(i = 0; i < m; i++)
  {
    l[i][0] = ax[i]; /* lines thru maximal vectors */
    l[i][1] = ay[i];
    l[i][2] = au[i];
    l[i][3] = av[i];

    lp[i][0] = ax[i]; /* lines normal to maximal vectors */
    lp[i][1] = ay[i];
    lp[i][2] = av[i];
    lp[i][3] = -au[i];
  }

  mm = m * (m - 1) / 2; /* number of intersection points      */

  for(n_flags = k = i = 0; i < m - 1; i++) /* first find the intersection points */
  {
    for(p = 0, j = i + 1; j < m; j++, k++)
    {
      if((flags[k] = intrsct(&xpt[k], &ypt[k], l[i], l[j])) != 0)
      {
        n_flags++;
        while((p % m == i) || (p % m == j))
          p++;
        jj = p % m;
        sign[k] = tsign(au[jj] * (ay[jj] - ypt[k]) - av[jj] * (ax[jj] - xpt[k]));
      }
      else
        sign[k] = 0.0;
    }
  }

  for(k = 0; k < mm; k++) /* identify intersection points on boundary */
    if(flags[k])          /* of active region.                        */
      for(j = 0; j < m; j++)
        if(((ds = tsign(CURL(au[j], av[j], ax[j] - xpt[k], ay[j] - ypt[k]))) != 0) && (ds != sign[k]))
        {
          flags[k] = -1;
          break;
        }

  for(i = 0; i < m; i++) /*  find good rotation point for point i */
  {
    sflag[i] = 0;
    xt[i] = yt[i] = 0.0;
    for(n_flags = 0, d = BIGTOL, j = 0; j < m - 1; j++)
    { /* find intersection of lp[i] with other */
      /* lines, and choose one on boundary     */
      /* of active region with largest curl    */
      jj = (i + j + 1) % m;
      if((flag = (intrsct(&xp, &yp, lp[i], l[jj]))) != 0)
      {
        temp = CURL(au[i], av[i], ax[i] - xp, ay[i] - yp);
        ds = tsign(temp);
        temp = fabs(temp);
        for(k = 1; k < m; k++) /* is point on boundary of */
        /* active region?          */
        {
          kk = (jj + k) % m;
          if((kk != i) && (tsign(CURL(au[kk], av[kk], ax[kk] - xp, ay[kk] - yp)) != ds))
          {
            flag = 0;
            break;
          }
        } /* end for k */

        if((flag) && (temp > d))
        {
          n_flags++;
          d = temp;
          xt[i] = xp;
          yt[i] = yp;
          sflag[i] = 1;
        }
      }
    }

    for(k = 0; k < mm; k++)
    { /* choose vertex of active     */
      /* region with largest curl    */
      /* compared to previous result */
      temp = fabs(CURL(au[i], av[i], ax[i] - xpt[k], ay[i] - ypt[k]));
      if((flags[k] == 1) && (temp > d))
      {
        n_flags++;
        d = temp;
        xt[i] = xpt[k];
        yt[i] = ypt[k];
        sflag[i] = 1;
      }
    } /* end for k */

    if(n_flags) /* any intersections found here?        */
      continue; /* yes ... continue                     */
    /* otherwise, find point in active region */
    /* (or on boundary)                       */

    for(xcm = ycm = 0.0, kk = k = 0; k < mm; k++)
      if(flags[k] == 1)
      {
        xcm += xpt[k];
        ycm += ypt[k];
        kk++;
      }

    if(kk == 0)
      return;

    xcm /= kk;
    ycm /= kk;

    ltemp[0] = xcm;
    ltemp[1] = ycm;
    ltemp[2] = av[i];
    ltemp[3] = -au[i];
    intrsct(&x11, &y11, ltemp, l[i]);
    for(fac = .5, j = 0; j < 30; j++, fac = fac / 2.0)
    {
      x1 = x11 + fac * cc * ltemp[2] / sqd;
      y1 = y11 + fac * cc * ltemp[3] / sqd;
      ds = tsign(CURL(au[0], av[0], ax[0] - x1, ay[0] - y1));

      for(flag = 1, k = 1; k < m; k++) /* check if in active region */
        if(tsign(CURL(au[k], av[k], ax[k] - x1, ay[k] - y1)) != ds)
        {
          flag = 0;
          break;
        }

      if(flag)
      {
        xt[i] = x1;
        yt[i] = y1;
        sflag[i] = 1;
        break;
      }

      x1 = x11 - fac * cc * ltemp[2] / sqd;
      y1 = y11 - fac * cc * ltemp[3] / sqd;
      ds = tsign(CURL(au[0], av[0], ax[0] - x1, ay[0] - y1));

      for(flag = 1, k = 1; k < m; k++) /* check if in active region */
        if(tsign(CURL(au[k], av[k], ax[k] - x1, ay[k] - y1)) != ds)
        {
          flag = 0;
          break;
        }

      if(flag)
      {
        xt[i] = x1;
        yt[i] = y1;
        sflag[i] = 1;
        break;
      }
    }
  }

  for(*x0 = 0.0, *y0 = 0.0, mm = k = 0; k < m; k++)
  {
    if(!sflag[k])
      continue;

    *x0 += xt[k];
    *y0 += yt[k];
    mm++;
  }

  if(mm)
  {
    *x0 /= mm;
    *y0 /= mm;
  }
}

// end of rotpt.cpp code

// start of cmdists.cpp code

/*****************************************************************************/
/*  "cmdists":                                                               */
/*  Compute distances between pairs of nominal and measured points           */
/*  and attach a sign depending on a rotation about (x0,y0).                 */
/*****************************************************************************/

void CMinMax::cmdists(double* r, DMAT* xy)
{
  int i;
  double ui, vi, ri;

  for(i = 0; i < numpts; i++)
  {
    ui = noms.m[i][0] - xy->m[i][0];
    vi = noms.m[i][1] - xy->m[i][1];
    ri = _hypot(ui, vi) / zone.v[i];
    r[i] = ((xy->m[i][0] * vi - xy->m[i][1] * ui < 0.0) ? -ri : ri);
  }
}

/*****************************************************************************/
/*  "abdists":                                                               */
/*  Compute distances between pairs of nominal and measured points           */
/*****************************************************************************/

double CMinMax::abdists(double* r, int* k, int* I, int* m, DMAT* uv)
{
  double ri;
  double sqd = 0.0;
  int i;

  for(*k = 1, i = 0; i < numpts; i++)
  {
    ri = _hypot(uv->m[i][0], uv->m[i][1]) / zone.v[i];

    if(ri > sqd + BIGTOL)
    {
      *k = i;
      sqd = ri;
    }
    r[i] = ri;
  }
  for(*m = 0, i = 0; i < numpts; i++)
    *m += (I[i] = (fabs(r[i] - sqd) * maxtol < tol));

  return (sqd);
}

// end of cmdists.cpp code

// start of brent.cpp code

double CMinMax::func(double theta, double x1, double y1, double x2, double y2, double z1, double w1, double z2,
                     double w2, double e1, double e2)
{
  double s, c;
  double rotx11, rotx12, rotx21, rotx22;

  s = sin(theta);
  c = cos(theta);
  rotx11 = c * x1 - s * y1;
  rotx12 = s * x1 + c * y1;
  rotx21 = c * x2 - s * y2;
  rotx22 = s * x2 + c * y2;
  return (_hypot(z1 - rotx11, w1 - rotx12) / e1 - _hypot(z2 - rotx21, w2 - rotx22) / e2);
}

double CMinMax::brent(double a, double b, double x1, double y1, double x2, double y2, double z1, double w1, double z2,
                      double w2, double e1, double e2)
{
  double fa;        /* f(a)                             */
  double fb;        /* f(b)                             */
  double fc;        /* f(c)                             */
  double prev_step; /* Distance from the last but one   */
  /* to the last approximation        */
  double tol_act; /* Actual tolerance                 */
  double p;       /* Interpolation step is calcu-     */
  double q;       /* lated in the form p/q; divi-     */
  /* sion operations is delayed       */
  /* until the last moment            */
  double new_step; /* Step at this iteration           */
  double e = 0.0, c = 0.0;
  double r, s;

  fa = func(a, x1, y1, x2, y2, z1, w1, z2, w2, e1, e2);
  fb = func(b, x1, y1, x2, y2, z1, w1, z2, w2, e1, e2);
  fc = fb;
  while(fb != 0.0) /* main iteration loop */
  {
    if((fb > 0.0) == (fc > 0.0))
    {
      c = a;
      fc = fa;
      prev_step = b - a;
      e = prev_step;
    }
    if(fabs(fc) < fabs(fb))
    {
      a = b;
      b = c;
      c = a;   /* Swap data for b to be the    */
      fa = fb; /* best approximation           */
      fb = fc;
      fc = fa;
    }
    new_step = 0.5 * (c - b); /* test for convergence */
    tol_act = 2.0 * BIGTOL * MAX(fabs(b), 1.0);
    if(fabs(new_step) <= tol_act || fb == 0.0)
      break; /* acceptable approximation was found */

    /* decide between bisection and interpolation */
    if(fabs(e) < tol_act || fabs(fa) <= fabs(fb))
    {
      prev_step = new_step; /* do bisection */
      e = new_step;
    }
    else
    {
      s = fb / fa; /* interpolation */
      if(a == c)   /* If we have only two distinct   */
      {            /* points, do linear interpolation */
        p = 2.0 * new_step * s;
        q = 1.0 - s;
      }
      else
      { /* do inverse quadratic interpolation */
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * new_step * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }

      if(p > 0.0) /* p was calculated with the op-*/
        q = -q;   /* posite sign; make p positive   */
      else        /* and assign possible minus to   */
        p = -p;   /* q            */

      /*  interpolated point acceptable */
      if(2.0 * p < 3.0 * new_step * q - fabs(tol_act * q) && p < fabs(0.5 * e * q))
      {
        prev_step = p / q;
        e = prev_step;
      }
      else
      {
        prev_step = new_step;
        e = new_step;
      }
    }

    a = b; /* move on to next point */
    fa = fb;
    if(fabs(prev_step) > tol_act)
      b += prev_step;
    else if(b > c)
      b -= tol_act;
    else
      b += tol_act;

    fb = func(b, x1, y1, x2, y2, z1, w1, z2, w2, e1, e2);
  }
  return (b);
}
