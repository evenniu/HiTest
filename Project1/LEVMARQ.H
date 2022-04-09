#pragma once

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

class DLLEXPORT CLevMarq
{
public:

  CLevMarq(int np);
  ~CLevMarq();
  int Allocate();
  void Destroy();

  int npp;

  double *c;
  double **A0;
  double **ASave;

  int levmq(int n, int np, double **A, double *b, double *alpha,
    double delta, double *bnorm, double *lambda, double *linres,
    int update, int first);
};
