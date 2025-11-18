#ifndef AUTOFIT_H
#define AUTOFIT_H

typedef struct dmat
{
  int dim1,dim2;
  double **m;
  dmat()
  {
    dim1 = dim2 = 0;
    m = NULL;
  }
}
DMAT;

typedef struct dvec
{
  int dim;
  double *v;
  dvec()
  {
    dim = 0;
    v = NULL;
  }
}
DVEC;

typedef struct ivec
{
  int dim;
  int *v;
  ivec()
  {
    dim = 0;
    v = NULL;
  }
}
IVEC;

typedef struct ALIGNREC   /* AUTOFIT */
{
  double *mat[3];   /* these need to be malloced */
  double orig[3];
}
ALIGNREC;

class CMinMax
{
public:
  CMinMax(int np);
  ~CMinMax();


  int numpts;

  IVEC mia1;                 /* misc integer array 1 */
  //DVEC wt;                   /* weight tolerance */      // Brian Shier
  double transx;             /* two dimensional major axis translation */
  double transy;             /* two dimensional minor axis translation */
  double spin;               /* two dimensional rotation */
  double maxtol;             /* largest tolerance */
  double cc;                 /* approximate size of part */
  double tol;                /* alogithm tolerances */
  DVEC    mda1;              /* misc double array 1 */
  DVEC    mda2;              /* misc double array 2 */
  DVEC    mda3;              /* misc double array 3 */
  DVEC    mda4;              /* misc double array 4 */
  DVEC    mda5;              /* misc double array 5 */
  DVEC    mda6;              /* misc double array 6 */
  //DVEC    tottol;            /* total tolerance */
  DVEC    zone;              /* tolerance used by algoritms */
  DVEC    omega;             /* weight multipliers */
  DMAT    noms;              /* Nominals */
  DMAT    vals;              /* Measured values */
  DMAT    infs;              /* solution */
  DMAT    uvw;
  DMAT    rtp;
  DMAT    mdda1;             /* misc double double array 1 */
  ALIGNREC cmmalign;

  bool MinMaxFit(bool translateOnly=false);
  void TransferIn(CMatrix *noms, CMatrix *vals, CVector *weights, CAlignment *align);
  void TransferOut(CMatrix *infs, CAlignment *align);
  void free_arrays(void);
  int alloc_arrays(void);
  int autofit(void);
  void calc_omegas(int nUseTols, int nFitType, BOOL bTolWeights=TRUE);
  void mmfit(bool translateOnly);
  double afcentroid(int, int, double **, double *, double *);
  void freeivector(IVEC *);
  void freedvector(DVEC *);
  void freedmatrix(DMAT *);
  void ivector(int, IVEC *);
  void dmatrix(int, int, DMAT *);
  void dvector(int, DVEC *);
  void dcopy(int n, double *x, double *y);
  void icopy(int n, int *x, int *y);
  void mcopy(int rows, int cols, double **A, double **B);
  double dsign(double a);
  double tsign(double a);
  void setmxmn(double *amax,double *amin,int *imax,int *imin, DMAT *xy, DMAT *zw, double *r);
  double fmaxv(int *imin,int *imax,DMAT *uv,int j,int *I);
  double SQR(double);
  void afrotate(double theta, DMAT *xyo, DMAT *xyi);
  int intrsct(double *x0,double *y0,double *l1,double *l2);
  double fmaxu(int m,double *au,double *av,int j);
  void calc_align(double *dx,double *dy,double *dtheta);
  void scale_weights();
  void calc_omegas(int nUseTols);

  double trans(DMAT *, DMAT *);
  double rotop(DMAT *, DMAT *);
  void rotpt(DMAT *,DMAT *,double *,double *);

  double findmin(DMAT *uv,double v1,double v2);
  void cmdists(double *,DMAT *);
  double abdists(double *, int *, int *, int *, DMAT *);
  double func(double, double, double, double, double, double, double, double, double, double, double);
  double brent(double, double, double, double, double, double, double, double, double, double, double, double);
  double ivmax(int *imax,double *r);
  double ivmin(int *imin,double *r);
};

#endif
