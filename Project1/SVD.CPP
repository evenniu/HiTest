#include "StdAfx.h"
/*************************************************************************/
/*  This file contains a routine for computing the SVD of a matrix.      */
/*  Based on a Fortran program of G. Stewart.                            */
/*************************************************************************/

#include "svd.h"

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

/*----------------------------------------------------------------------*/
/*  DSVDC                                                               */
/*     dsvdc is a routine to reduce a double precision n X p matrix A   */
/*     by orthogonal transformations U and V to diagonal form.  The     */
/*     diagonal elements s[i] are the singular values of A.  The        */
/*     columns of U are the corresponding left singular vectors,        */
/*     and the columns of V the right singular vectors.                 */
/*  Input:                                                              */
/*     double A[n][p] A contains the matrix whose singular value        */
/*                    decomposition is to be computed.  A is            */
/*                    destroyed by dsvdc.                               */
/*     int n          n is the number of rows of the matrix A.          */
/*     int p          p is the number of columns of the matrix A.       */
/*     double work[n] scratch array.                                    */
/*     int job        job controls the computation of the singular      */
/*                    vectors.  It has the decimal expansion ab         */
/*                    with the following meaning                        */
/*                        a = 0     Do not compute the left singular    */
/*                                  vectors.                            */
/*                        a = 1     Return the n left singular vectors  */
/*                                  in U.                               */
/*                        a >= 2    Return the first min(n,p) singular  */
/*                                  vectors in U.                       */
/*                        b = 0     Do not compute the right singular   */
/*                                  vectors.                            */
/*                        b = 1     Return the right singular vectors   */
/*                                  in V.                               */
/*     int sort       If sort != 0 the singular values are sorted;      */
/*                    otherwise they are not.                           */
/*  Output:                                                             */
/*     double s[mm],  mm = min(n+1,p).                                  */
/*                    The first min(n,p) entries of s contain the       */
/*                    singular values of A arranged in descending       */
/*                    order of magnitude.                               */
/*     double e[p]    e ordinarily contains zeros.  However, see        */
/*                    the discussion of info for exceptions.            */
/*     double U[n][k] If joba = 1 then k = n,                           */
/*                    if joba >= 2 then k = min(n,p).                   */
/*                    U contains the matrix of left singular vectors.   */
/*                    U is not referenced if job = 0.  If n <= p        */
/*                    or if job = 2, then U may be identified with A    */
/*                    in the function call.                             */
/*     double V[p][p] V contains the matrix of right singular vectors.  */
/*                    V is not referenced if job = 0.  If p <= n,       */
/*                    then V may be identified with A in the            */
/*                    function call.                                    */
/*     int info       The return value.                                 */
/*                    The singular values (and their corresponding      */
/*                    singular vectors) s[info],...,s[m-1]              */
/*                    are correct (here m=min(n,p)).  Thus, if          */
/*                    info = 0, all the singular values and their       */
/*                    vectors are correct.  In any event, the matrix    */
/*                    B = U'*A*V is the bidiagonal matrix               */
/*                    with the elements of s on its diagonal and the    */
/*                    elements of e on its superdiagonal (U'            */
/*                    is the transpose of U).  Thus the singular        */
/*                    values of A and B are the same.                   */
/*----------------------------------------------------------------------*/

int dsvdc(double **A, int n, int p, double *s, double *e, double **U,
          double **V, double *work, int job, int sort)
{
  int i, j, k, l, mm, lu, ll, kk, ls=0, i1;
  int nct, ncu, lls, nrt;
  int maxit, getU, getV, kj, flag, jobu, iter;
  double d1, test, b, c, f, g, t, scale, sc, shift, temp, sum;
  double cs, t1, ztest, el, sl, sm, sn, ek, sk, dx, dy;

  maxit = 30;

  /* Decide what to compute. */
  getU = 0;
  getV = 0;
  jobu = job % 100 / 10;
  if (jobu > 1)
    ncu = MIN(n, p);
  else
    ncu = n;
  if (jobu != 0)
    getU = 1;
  if (job % 10 != 0)
    getV = 1;

  /* Reduce A to bidiagonal form, storing the diagonal elements */
  /* in s and the superdiagonal elements in e.                           */
  nct = MIN(n-1, p);
  i1 = MIN(p-2, n);
  nrt = MAX(0, i1);
  lu = MAX(nct, nrt);

  for (l = 0; l < lu; l++)
  {
    if (l < nct)
    {
      /* Compute the transformation for the l-th column and */
      /* place the l-th diagonal in s[l:n-1].                      */
      for (sum = 0., k = l; k < n; k++)
        sum += (temp = A[k][l], temp * temp);
      s[l] = sqrt(sum);
      if (s[l] != 0.)
      {
        if (A[l][l] != 0.)
          s[l] = DSIGN2(s[l], A[l][l]);
        for (sc = s[l], k = l; k < n; k++)
          A[k][l] /= sc;
        A[l][l] += 1.;
        s[l] = -s[l];
      }
    }
    for (j = l+1; j < p; j++)
    {
      if (l < nct && s[l] != 0.)
      {
        /* Apply the transformation. */
        for (t = 0., k = l; k < n; k++)
          t += A[k][l] * A[k][j];
        t /= -A[l][l];
        for (i = l; i < n; i++)
          A[i][j] += t * A[i][l];
      }
      /* Place the l-th row of A into e for the                  */
      /* subsequent calculation of the row transformation. */
      e[j] = A[l][j];
    }
    if (getU && l < nct)
    {
      /* Place the transformation in U for subsequent back */
      /* multiplication.                                                    */
      for (i = l; i < n; i++)
        U[i][l] = A[i][l];
    }
    if (l < nrt)
    {
      /* Compute the l-th row transformation and place the */
      /* l-th super-diagonal in e(l).                                 */
      for (sum = 0., k = l+1; k < p; k++)
        sum += (temp = e[k], temp * temp);
      e[l] = sqrt(sum);
      if (e[l] != 0.)
      {
        if (e[l+1] != 0.)
          e[l] = DSIGN2(e[l], e[l+1]);
        for (sc = e[l], k = l+1; k < p; k++)
          e[k] /= sc;
        e[l+1] += 1.;
        e[l] = -e[l];
        /* Apply the transformation. */
        for (i = l+1; i < n; i++)
        {
          work[i] = 0.;
          for (j = l+1; j < p; j++)
            work[i] += e[j] * A[i][j];
        }
        for (i = l+1; i < n; i++)
          for (j = l+1; j < p; j++)
            A[i][j] -= e[j]/e[l+1] * work[i];
      }
      if (getV)
      {
        /* Place the transformation in V for subsequent */
        /* back multiplication.                                     */
        for (i = l+1; i < p; i++)
          V[i][l] = e[i];
      }
    }
  }
  /* Set up the final bidiagonal matrix of order mm. */
  mm = MIN(p, n+1);
  if (nct < p)
    s[nct] = A[nct][nct];
  if (mm > n)
    s[mm-1] = 0.;
  if (nrt < mm-1)
    e[nrt] = A[nrt][mm-1];
  e[mm-1] = 0.;
  /*   if required, generate U. */
  if (getU)
  {
    for (j = nct; j < ncu; j++)
    {
      for (i = 0; i < n; i++)
        U[i][j] = 0.;
      U[j][j] = 1.;
    }
    for (ll = 0; ll < nct; ll++)
    {
      l = nct - ll;
      if (s[l-1] == 0.)
      {
        for (i = 0; i < n; i++)
          U[i][l-1] = 0.;
        U[l-1][l-1] = 1.;
      }
      else
      {
        for (j = l; j < ncu; j++)
        {
          for (t = 0., k = l-1; k < n; k++)
            t += U[k][l-1] * U[k][j];
          t /= -U[l-1][l-1];
          for (i = l-1; i < n; i++)
            U[i][j] += t * U[i][l-1];
        }
        for (i = l-1; i < n; i++)
          U[i][l-1] *= -1.;
        U[l-1][l-1] += 1.;
        for (i = 0; i < l-1; i++)
          U[i][l-1] = 0.;
      }
    }
  }
  /* If required, generate V. */
  if (getV)
  {
    for (ll = 0; ll < p; ll++)
    {
      l = p - ll;
      if (l <= nrt && e[l-1] != 0.)
      {
        for (j = l; j < p; j++)
        {
          for (t = 0., k = l; k < p; k++)
            t += V[k][l-1] * V[k][j];
          t /= -V[l][l-1];
          for (i = l; i < p; i++)
            V[i][j] += t * V[i][l-1];
        }
      }
      for (i = 0; i < p; i++)
        V[i][l-1] = 0.;
      V[l-1][l-1] = 1.;
    }
  }
  /* Main iteration loop for the singular values. */
  kk = mm;
  iter = 0;
  while (mm > 0)
  {
    /* If too many iterations have been performed, set */
    /* flag and return.                                                */
    if (iter >= maxit)
      return mm-1;
    for (ll = 0; ll < mm; ll++)
    {
      l = mm - ll - 1;
      if (l == 0)
        break;
      test = fabs(s[l-1]) + fabs(s[l]);
      ztest = test + fabs(e[l-1]);
      if (ztest == test)
      {
        e[l-1] = 0.;
        break;
      }
    }
    if (l == mm - 1)
      flag = 4;
    else
    {
      for (lls = l; lls < mm+1; lls++)
      {
        ls = mm - lls + l;
        if (ls == l)
          break;
        test = 0.;
        if (ls != mm)
          test += fabs(e[ls-1]);
        if (ls != l + 1)
          test += fabs(e[ls-2]);
        ztest = test + fabs(s[ls-1]);
        if (ztest == test)
        {
          s[ls-1] = 0.;
          break;
        }
      }
      if (ls == l)
        flag = 3;
      else
      {
        if (ls == mm)
          flag = 1;
        else
        {
          flag = 2;
          l = ls;
        }
      }
    }
    l++;
    /* Perform the task indicated by flag. */
    switch (flag)
    {
    case 1:   /* deflate negligible s[mm]. */
      f = e[mm - 2];
      e[mm - 2] = 0.;
      for (kj = l; kj < mm; kj++)
      {
        k = mm - kj + l - 1;
        t1 = s[k-1];
        SGivens(&t1, &f, &cs, &sn);
        s[k-1] = t1;
        if (k != l)
        {
          f = -sn * e[k - 2];
          e[k - 2] = cs * e[k - 2];
        }
        if (getV)
        {
          for (i = 0; i < p; i++)
          {
            dx = V[i][k-1];
            dy = V[i][mm-1];
            V[i][mm-1] = cs * dy - sn * dx;
            V[i][k-1] = cs * dx + sn * dy;
          }
        }
      }
      continue;
    case 2:   /*   split at negligible s[l]. */
      f = e[l - 2];
      e[l - 2] = 0.;
      for (k = l-1; k < mm; k++)
      {
        t1 = s[k];
        SGivens(&t1, &f, &cs, &sn);
        s[k] = t1;
        f = -sn * e[k];
        e[k] = cs * e[k];
        if (getU)
        {
          for (i = 0; i < n; i++)
          {
            dx = U[i][k];
            dy = U[i][l-2];
            U[i][l-2] = cs * dy - sn * dx;
            U[i][k] = cs * dx + sn * dy;
          }
        }
      }
      continue;
    case 3:   /*   perform one QR step. */
      /*   calculate the shift. */
      d1 = MAX(fabs(s[mm-1]), fabs(s[mm - 2]));
      d1 = MAX(d1, fabs(e[mm - 2]));
      d1 = MAX(d1, fabs(s[l-1]));
      scale = MAX(d1, fabs(e[l-1]));
      sm = s[mm-1] / scale;
      sk = s[mm - 2] / scale;
      ek = e[mm - 2] / scale;
      sl = s[l-1] / scale;
      el = e[l-1] / scale;
      b = ((sk + sm) * (sk - sm) + ek * ek) / 2.;
      c = (d1 = sm * ek, d1 * d1);
      shift = 0.;
      if (b != 0. || c != 0.)
      {
        shift = sqrt(b * b + c);
        if (b < 0.)
          shift = -shift;
        shift = c / (b + shift);
      }
      f = (sl + sm) * (sl - sm) + shift;
      g = sl * el;
      /* Chase zeros. */
      for (k = l; k < mm; k++)
      {
        SGivens(&f, &g, &cs, &sn);
        if (k != l)
          e[k - 2] = f;
        f = cs * s[k-1] + sn * e[k-1];
        e[k-1] = cs * e[k-1] - sn * s[k-1];
        g = sn * s[k];
        s[k] = cs * s[k];
        if (getV)
        {
          for (i = 0; i < p; i++)
          {
            dx = V[i][k-1];
            dy = V[i][k];
            V[i][k] = cs * dy - sn * dx;
            V[i][k-1] = cs * dx + sn * dy;
          }
        }
        SGivens(&f, &g, &cs, &sn);
        s[k-1] = f;
        f = cs * e[k-1] + sn * s[k];
        s[k] = -sn * e[k-1] + cs * s[k];
        g = sn * e[k];
        e[k] = cs * e[k];
        if (getU && k < ncu)         /* was n */
        {
          for (i = 0; i < n; i++)
          {
            dx = U[i][k-1];
            dy = U[i][k];
            U[i][k] = cs * dy - sn * dx;
            U[i][k-1] = cs * dx + sn * dy;
          }
        }
      }
      e[mm - 2] = f;
      iter++;
      continue;
    case 4:   /* Convergence. */
      /* Make the singular value positive. */
      if (s[l-1] < 0.)
      {
        s[l-1] = -s[l-1];
        if (getV)
          for (i = 0; i < p; i++)
            V[i][l-1] *= -1.;
      }
      /* If required, order the singular values   */
      /* and the corresponding singular vectors. */
      if (sort)
      {
        while (l != kk && s[l-1] < s[l])
        {
          t = s[l-1];
          s[l-1] = s[l];
          s[l] = t;
          if (getV && l < p)
          {
            for (i = 0; i < p; i++)
            {
              temp = V[i][l-1];
              V[i][l-1] = V[i][l];
              V[i][l] = temp;
            }
          }
          if (getU && l < ncu)   /* was n */
          {
            for (i = 0; i < n; i++)
            {
              temp = U[i][l-1];
              U[i][l-1] = U[i][l];
              U[i][l] = temp;
            }
          }
          l++;
        }
      }
      iter = 0;
      mm--;
    }
  }
  return 0;
}
