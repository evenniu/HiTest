#pragma once

#include <stdio.h>

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

class DLLEXPORT CNurbCurve : public CCurve
{
public:
  CNurbCurve(FILE *fp, bool isEnglish);
  //CNurbCurve(int fd, bool isEnglish);
  CNurbCurve(int n, double *x, double *y, double *z = 0, bool isEnglish = true, int cl=0, int tn=1, int fx=1,
             double tc=0.0, CDumpRecordArray *splDmp=NULL, int index1=-1, int index2=-1, double mult=0.0, bool showDebug=false);
  CNurbCurve(int n, CCurve *ccc, CCurve *cvc, CCurve *lec = 0, CCurve *tec = 0, int useextreme = 0, double *nomStart = NULL);
  CNurbCurve(int n, double lextr, double textr, double ler, double ter, CCurve *whole, int useextremes);
  ~CNurbCurve();

  virtual int Type() const override
  {
    return NURB_TYPE;
  }

  void Knot(int i, double v)
  {
    m_pknot[i] = v;
  }
  void Weight(int i, double w)
  {
    m_ppoints->m[i][0] = w;
    m_ppoints->m[i][1] = w;
    m_ppoints->m[i][2] = w;
#ifdef DIM3

    m_ppoints->m[i][3] = w;
#endif

  }
  void Control(int i, double *v, double *k=0)
  {
    m_ppoints->m[i][0] *= v[0];
    m_ppoints->m[i][1] *= v[1];
#ifdef DIM3

    m_ppoints->m[i][2] *= v[2];
#endif

    if (k && m_ppointk)
      m_ppointk[i] = *k;
  }
  int CalcPoint(double *xyz, double a, double *tan = 0, double *curv = 0) const override;
  int Valid()
  {
    return m_good;
  }
  void Valid(int v)
  {
    m_good = v;
  }
  void Domain(double s, double b)
  {
    m_t0=s;
    m_t1=b;
  }
  void Closed(int c)
  {
    m_closed = (short)c;
  }
  void Poly(int p)
  {
    m_poly = (short)p;
  }

  int KNum()
  {
    return (int)m_kNum;
  }
  int Closed()
  {
    return m_closed;
  }
  void SetT0(double t)
  {
    m_t0 = t;
  }
  void SetT1(double t)
  {
    m_t1 = t;
  }
  double T0() const override
  {
    return m_t0;
  }
  double T1() const override
  {
    return m_t1;
  }
  double Knot(int i)
  {
    return m_pknot[i];
  }
  long StoreSize() //WWW make unsigned long
  {
    return sizeof(short) + sizeof(double) + m_storesize;
  }
  void Write(FILE *fp);
  void Extent(double *min, double *max);
  void Extent();

  double TotalLength(int numPoints, double t0, double t1);
  void Walk(double &t, double d, double *pt, int numPoints);
  void Debug(int numPoints);

  void CenterOfMass(double *CM, bool isGEI);
  double ComputeArea(double t0 = 0.0, double t1 = 0.0);

  bool SlidePerUnitWindow(double dbLength, int index, int &index2, double t0, double &t1);
  double arc_inverse(double length, double t0, int &index);
  double compute_length(double t0 = 0.0, double t1 = 0.0);
  double deriv_func(double t);


//protected:
public:  // need to get the data members from outside, and don't want to write access functions for each one.

  void MeanCamber(int n, CCurve *ccc, CCurve *cvc, CCurve *lec = 0, CCurve *tec = 0, int useextreme = 0, double *nomStart = NULL);

  int Allocate(int havek);  // allocate memory
  int FindSpan(double u) const;
  void DersBasisFuns(int i, double u, int n=0) const;

  int m_nSpan, m_pSpan;// used internally in FindSpan
  int m_mid;           // used internally in FindSpan
  int m_good;          // set 0 if new's fail
  short m_Num;         // size of control point array
  short m_Order;       // order of curve
  short m_kNum;        // size of knot array
  short m_closed;      // closed
  short m_poly;        // if 1 is polynomial, not rational
  long m_storesize;    // bytes needed to write to binary file
  double *m_pknot;     // knot array
  CMatrix *m_pbf;      // basis functions and derivatives
  CMatrix *m_ppoints;  // control points (x*w,y*w,z*w,w)
  double *m_ppointk;   // control points for k component of normal vector
  double *m_pleft;     // used by DersBasisFuns
  double *m_pright;    // used by DersBasisFuns
  CMatrix *m_pndu;     // used by DersBasisFuns
  CMatrix *m_pa;       // used by DersBasisFuns
  double m_t0,m_t1;    // domain
  double m_min[2];
  double m_max[2];

  double GaussQuad(int integrand, double a, double b);
  double moment_x_func(double u);
  double moment_y_func(double u);
  double area_func(double u);

  // during reporting, it can be useful to retrieve information about the curve;
  // I call these annotations
  virtual void addDoubleAnnotation(const char* key, double value) override final;
  virtual double getDoubleAnnotation(const char* key) const override final;
  virtual void addBoolAnnotation(const char* key, bool value) override final;
  virtual bool getBoolAnnotation(const char* key) const override final;
};
