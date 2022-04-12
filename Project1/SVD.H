#pragma once

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

DLLEXPORT int dsvdc(double **A, int n, int p, double *s, double *e, double **U, double **V, double *work, int job, int sort);
