#pragma once
#include <memory>

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define DSIGN(x) ((x) > 0.0 ? 1.0 : ((x) < 0.0 ? -1.0 : 0.0))
#define DSIGN2(x,y) ((y) > 0.0 ? fabs(x) : ((y) < 0.0 ? -fabs(x) : 0.0))
#define SQR2(x,y)  (((x)*(x)) + ((y)*(y)))
#define CURL(x1,y1,x2,y2)  (((x1)*(y2)) - ((x2)*(y1)))
#define SGN(x) ((x) > 0.? 1. : ((x) < 0.? -1. : 0.))
#define ACOS(X) acos((X) > 1.0 ? 1.0 : ((X) < -1.0 ? -1.0 : (X)))

#define AINT -12.  // Parameters for the minimization routine fmin.
#define BINT  6.
#define FTOL  1.e-6
#define ITMAX 100
#define M_SQRT1_2  0.70710678118655
#define M_SQRT2    1.4142135623730950488016887242097
#define M_E         2.71828182845904523536
#define M_LOG2E     1.44269504088896340736
#define M_LOG10E    0.434294481903251827651
#define M_LN2       0.693147180559945309417
#define M_LN10      2.30258509299404568402
#define M_2PI       6.2831853071795864769253
#define M_PI        3.14159265358979323846
#define M_PI_2      1.57079632679489661923
#define M_PI_4      0.785398163397448309616
#define M_1_PI      0.318309886183790671538
#define M_2_PI      0.636619772367581343076
#define M_1_SQRTPI  0.564189583547756286948
#define M_2_SQRTPI  1.12837916709551257390
#define M_SQRT_2    0.707106781186547524401
class CMatrix;

struct CompRecord
{
  short type;
  double xy1[2];
  double xy2[2];
  int index;
  CompRecord()
  {
    type = 0;
    index = 0;
    xy1[0] = xy1[1] = xy2[0] = xy2[1];
  }
  CompRecord(short typ, double x1, double y1, double x2, double y2, int ind)
  {
    type = typ;
    xy1[0] = x1;
    xy1[1] = y1;
    xy2[0] = x2;
    xy2[1] = y2;
    index = ind;
  }
};

struct MapRecord
{
  double v[6];
  MapRecord(double v0, double v1, double v2, double v3, double v4, double v5)
  {
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
    v[3] = v3;
    v[4] = v4;
    v[5] = v5;
  }
};

typedef vector < CompRecord * > CCompRecordArray;
typedef vector < MapRecord * > CMapRecordArray;
//typedef CArray < DumpRecord *, DumpRecord * > CDumpRecordArray; convert to std::vector


DLLEXPORT extern void ang_to_dms(double angle, wchar_t *buf);

// These methods return the angle between two vectors
// The returned angle is in the range [0, pi]
DLLEXPORT double AngleBetweenVectors2d(const double * vec1, const double* vec2);
DLLEXPORT double AngleBetweenVectors3d(const double * vec1, const double* vec2);


DLLEXPORT double clip(double value, double low, double high);

DLLEXPORT double ddot(int n, double *x, double *y);
DLLEXPORT int FindBreakPoint(double, double *, int, int);
DLLEXPORT double l2norm(const double *, int ndim=2);
DLLEXPORT double sqnorm(double *dx, int n=2);
DLLEXPORT double dot(const double *, const double *, int ndim=2);
DLLEXPORT double normalize(double * vo, const double * vi, int ndim=2);
DLLEXPORT void SGivens(double *, double *, double *, double *);
DLLEXPORT double curl(double *, double *);

DLLEXPORT int int_line_line(double, double, double, double, double, double, double, double, double *, double *);
DLLEXPORT int int_line_line(double *p1, double *v1, double *p2, double *l2, double *sol);
DLLEXPORT int parabola(double *, double *, double *, double *, double *);
DLLEXPORT int int_para_line(double, double, double, double, double, double, double, double *, double *);
DLLEXPORT int int_para_circ(double, double, double, double, double, double, double *, double *);
DLLEXPORT void para_norm(double, double, double, double, double *, double *);

DLLEXPORT void rotate(int, double, double *, double *, double *, double *);
DLLEXPORT void multAx(int rows, int cols, double **A, double *x, double *y);
DLLEXPORT void centroid(int n, int ndim, double **X, double *X0, double *w = 0);
DLLEXPORT double wcentroid(int n, int ndim, double **X, double *w, double *X0);
DLLEXPORT void translate_all(int npoints, int ndim, double *u, double **P);
DLLEXPORT void shift_all(int npoints, int ndim, double *u, double **P);
DLLEXPORT double sarea(double *x, double *y, int k, int j, int i);
DLLEXPORT double varea(double *x, double *y, double *z, int k, int j, int i, double *N);
DLLEXPORT void solvetrid(int n, double *a, double *d, double *b, double *r, double *x);
DLLEXPORT void tridlu(int n, double *a, double *d, double *b);
DLLEXPORT double solve_cubic(double a, double b, double c, int *flag);
DLLEXPORT double dist(double *dx, double *dy, int n=2);
DLLEXPORT double scalar(double *p1, double *p2, double *V, int ndim=2);
DLLEXPORT void simplex(int n, double *x, double l, double u);
DLLEXPORT void dzero(int n, double *x);
DLLEXPORT void lMcopy(int n, double **A, double **B);
DLLEXPORT void lower_backsub(int n, double **L, double *b);
DLLEXPORT int LDLTsolve(int n, double **A, double *b);
DLLEXPORT int LDLT(int n, double **A);
DLLEXPORT void upper_backsub(int n, double **L, double *y);
DLLEXPORT void insertion_sort_points(int n, double *t, double **x);
DLLEXPORT void CompensatePoint(int how, double prad, double zval, double *ball, double *ijk, double *comped);
DLLEXPORT double ptlinedist(double *lp, double *lv, double *pt, double *proj = 0);
DLLEXPORT double projdist(double *lp, const double *lv, double *pt, double *proj = 0);
DLLEXPORT int intersect(double *p1, double *p2, double *p3, double *p4);
DLLEXPORT double radius_3points(double *p1, double *p2, double *p3);
DLLEXPORT int int_ray_segment(double *rayStart, double *rayDir, double *p1, double *p2, double *inter = 0);
DLLEXPORT int InventEnd(int n1, int n2, int nold, double *x, double *y, int n, double *nx, double *ny, double mult=1.0,
                     double *newstart=0, double *newend=0);
DLLEXPORT int SupplementEnd(double *x, double *y, int n, double *nx, double *ny);

DLLEXPORT double norm_cross2d(const double v1[2], const double v2[2]);
DLLEXPORT double norm_cross3d(const double v1[3], const double v2[3]);
DLLEXPORT void cross(const double v1[3], const double v2[3], double vo[3]);
DLLEXPORT void cross_product(double *U, double *V, double *N);
DLLEXPORT void int_circ_line(double *cir, double r, double *p1, double *p2, double *sol);
DLLEXPORT int Colinear(double *p1, double *p2, double *p3, double eps);

DLLEXPORT bool CalculateLine(CMatrix *points, double *center, double *theta);
DLLEXPORT bool CalculateFunctionalLine(CMatrix *points, double *center, double *theta);

DLLEXPORT int fix_data(int n, double *pX, double *pY, double *pK, int closed, double howClose, bool checkEnds, bool showDebug=false);
DLLEXPORT int fix_data_new(int n, double* pX, double* pY, double* pK, int closed, double howClose, bool checkEnds,
                           bool showDebug = false, double theat = 0);

DLLEXPORT bool int_plane_line(double *pvec, double pd, double *lpt, double *lvec, double *sol);

DLLEXPORT bool PointInPolygon(int num, double *vx, double *vy, double testx, double testy);
DLLEXPORT double SignedAngle(double *v1, double *v2);
DLLEXPORT bool ParametersEqual(double period, double p1, double p2);
DLLEXPORT bool ParameterInRange(double p, double r1, double r2, double period);
DLLEXPORT double signedAngleBetween2DVectors(const double* a, const double* b);
DLLEXPORT double UnwrapAngle(double radians_good, double radians_towrap);
std::pair<double, double> FitRangeToPoint(double termination1, double termination2, double x, double period);


// This function returns a half-piece of a given section
// You give it a pointer to a Section object, and you tell it CVC or CCC
// and it gives you a CSubCurve for half the curve, first nominal and second measured
class CSubCurve;
class CSection;
std::unique_ptr<CSubCurve> MakeNominalHalfCurve(const CSection *in_section, unsigned int which_piece);
std::unique_ptr<CSubCurve> MakeMeasuredHalfCurve(const CSection *in_section, unsigned int which_piece);


//DLLEXPORT int RemoveLoops(int n, double *x, double *y, double *z, double *vi, double *vj);
DLLEXPORT int TrimEnd(int n1, CMatrix &sc1, int n2, CMatrix &sc2, int *s1, int *e2, int useInter, bool isEnglish);
DLLEXPORT int TrimTails(int n1, int n2, double *x, double *y, double *z, double *vi, double *vj,
                        int *end1, int *end2, int types, CCompRecordArray *cra, bool isEnglish);

DLLEXPORT __time64_t ConvertDateToTime64(int yr, int mo, int da, int hr, int mn, int sc);
