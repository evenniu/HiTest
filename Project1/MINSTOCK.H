#pragma once

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

class DLLEXPORT CMinStock
{
public:
  CMinStock(void);
  ~CMinStock(void);

  int m_numChecks; //number of checks.

  CMatrix *m_nomCV;
  CMatrix *m_nomCC;
  CMatrix *m_meaCV;
  CMatrix *m_meaCC;
  double *m_stock;
  double m_pTol[4];
  double *m_nomT;
  double m_mTol[4];
  CCurve *m_nom;
  CCurve *m_mea;
  CCurve *m_nMCL;
  CCurve *m_mMCL;
  int *m_side;
  bool m_calculated[4];

  CAlignment *m_firstAlign;
  CAlignment *m_alignCV;
  CAlignment *m_alignCC;
  CAlignment *m_alignLE;
  CAlignment *m_alignTE;

  double m_nomleext;
  double m_nomteext;
  double m_nomLength;
  double m_nomLEChd[2];
  double m_nomLECtr[2];
  double m_nomTEChd[2];
  double m_nomTECtr[2];
  double m_mclNomLERatio;
  double m_mclNomTERatio;
  double m_nomPeriod;

  double m_mealeext;
  double m_meateext;
  double m_meaLEChd[2];
  double m_meaLECtr[2];
  double m_meaTEChd[2];
  double m_meaTECtr[2];
  double m_mclMeaLERatio;
  double m_mclMeaTERatio;
  double m_meaPeriod;

  double m_deltaLE;
  double m_delta;
  double m_deltaTE;

  double m_leoff;
  double m_teoff;

  double m_minStock;
  double m_minStockLE;
  double m_minStockTE;
  int m_minIndex;
  int m_minIndexLE;
  int m_minIndexTE;

  double m_maxStock;
  double m_maxStockLE;
  double m_maxStockTE;
  int m_maxIndex;
  int m_maxIndexLE;
  int m_maxIndexTE;

  bool Initialize(double dLE, double d, double dTE);
  void SetCurves(CCurve *nom, CCurve *mea, CCurve *nMCL, CCurve *mMCL,
    double nomleext, double nomteext, double mealeext, double meateext, double leoff, double teoff);
  bool Calculate(int side, CAlignment *align, double startVal, double stopVal);
  bool ReCalculate(int side, CAlignment *align);

  bool GetNomCV(int i, double *xy);
  bool GetNomCC(int i, double *xy);
  bool GetMeaCV(int i, double *xy);
  bool GetMeaCC(int i, double *xy);
  double GetStock(int i);
  void SetPTol(int side, double tol);
  void SetMTol(int side, double tol);
  double GetPTol(int side);
  double GetMTol(int side);

  double CalcStock(int i, double t, double *nomCV, double *nomCC, double *meaCV, double *meaCC, double *nomMC);
};
