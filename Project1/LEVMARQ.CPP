#include "StdAfx.h"
#include "LevMarq.h"
#include <memory.h>

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

#define LM_TOL 0.75
#define AZERO  1.0e-10
#define LM_DIV 2.0

CLevMarq::CLevMarq(int np)
{
  npp = np;

  c = 0;
  A0 = 0;
  ASave = 0;

  Allocate();
}

CLevMarq::~CLevMarq()
{
  Destroy();
}

int CLevMarq::Allocate()
{
  if (npp < 2)
    return 0;

  c = new double[npp];

  int i;
  int s = (npp*(npp + 1))/2;

  A0 = new double *[npp];
  A0[0] = new double[s];
  for (i=1; i<npp; i++)
    A0[i] = A0[i-1] + i;
  for (i=0; i<s; i++)
    A0[0][i] = 0.0;

  ASave = new double *[npp];
  ASave[0] = new double[s];
  for (i=1; i<npp; i++)
    ASave[i] = ASave[i-1] + i;
  for (i=0; i<s; i++)
    ASave[0][i] = 0.0;

  return 1;
}

void CLevMarq::Destroy()
{
  if (c)
  {
    delete [] c;
    c = 0;
  }

  if (A0)
  {
    delete [] A0[0];
    delete [] A0;
    A0 = 0;
  }

  if (ASave)
  {
    delete [] ASave[0];
    delete [] ASave;
    ASave = 0;
  }
}

/*-------------------------------------------------------------*/
/*  LEVMQ  Solve the Levenberg-Marquardt subproblem.           */
/*  INPUT                                                      */
/*    n        Number of equations.                            */
/*    np       Number of parameters.                           */
/*    A        System matrix.                                  */
/*    b        Right-hand side (actually, A^Tb).               */
/*    alpha    Starting parameter vector.                      */
/*    delta    Starting trust region radius.                   */
/*  OUTPUT                                                     */
/*    alpha    Solution vector.                                */
/*    bnorm    Norm of solution.                               */
/*    lambda   Approximate Lagrange multiplier for delta.      */
/*    linres   Squared error.                                  */
/*-------------------------------------------------------------*/

int CLevMarq::levmq(int n, int np, double **A, double *b, double *alpha, double delta, double *bnorm,
                    double *lambda, double *linres, int update, int first)
{
  int i, j, k, its, flag = 0;
  double temp, dlambda, delta0, delta1, fac = 2.;
  static double lambda0, lambda1, l0, u0;

  if (np != npp)
  {
    npp = np;
    Destroy();
    Allocate();
    first = 1;
  }

  if (update)
  {
    /* Compute lower part of A^TA. */
    for (j = 0; j < np; j++)
    {
      for (i = j; i < np; i++)
      {
        for (temp = 0., k = 0; k < n; k++)
          temp += A[k][i] * A[k][j];
        ASave[i][j] = temp;
      }
    }
  }

  if (first)
  {
    lambda0 = l0 = AZERO;
    lambda1 = u0 = 1.e4;
  }

  delta0 = .5 * (1. - LM_TOL) * delta;
  delta1 = .5 * (1. + LM_TOL) * delta;

  /* Modify lambda until ||alpha(lambda)|| ~= delta. */
  for (its = 0; its < 30; its++)
  {
    /* Compute bnorm = ||alpha(lambda)||. */
    memcpy(alpha, b, np*sizeof(double));
    lMcopy(np, ASave, A0);
    for (i = 0; i < np; i++)
      A0[i][i] += *lambda;
    if (LDLTsolve(np, A0, alpha))
    {
      //bugout(0, L"Matrix singular.\n"); // Shouldn't happen since lambda > 0.
      return -9;
    }

    *bnorm = l2norm(alpha, np);

    if (*lambda <= AZERO && *bnorm <= delta)
    {
      flag = 0;
      *lambda = AZERO;
      break;
    }

    if (*bnorm < LM_TOL * delta)
    {
      if (flag == 1)
        fac /= LM_DIV;
      lambda0 = *lambda / (fac + 1.);
      if (lambda1 <= lambda0)
        lambda1 = u0;
      flag = -1;
    }
    else if (*bnorm > delta)
    {
      if (flag == -1)
        fac /= LM_DIV;
      lambda1 = *lambda * (fac + 1);
      if (lambda1 <= lambda0)
        lambda0 = l0;
      flag = 1;
    }
    else
    {
      flag = 1;
      break;
    }

    if (fabs(*bnorm - delta1) > delta0)
    {
      /* Adjust lambda by Newton's method. */
      memcpy(c, alpha, np*sizeof(double));
      lower_backsub(np, A0, c);
      for (i = 0; i < np; i++)
      {
        if ((temp = A0[i][i]) > 0.)
          c[i] /= sqrt(temp);
      }

      dlambda = (temp = *bnorm / l2norm(c,np), temp * temp) *
        (*bnorm / delta1 - 1.);
      *lambda += dlambda;

      if (*lambda < l0)
        *lambda = l0;
      if (*lambda > u0)
        *lambda = u0;
    }
  }

 /* if (its >= 30)
    bugout(0, L"Warning: no convergence in levmq.\n");*/

  for (*linres = 0., k = 0; k < np; k++)
  {
    temp = alpha[k];
    *linres -= temp * (b[k] + *lambda * temp);
  }

  return flag;
}
