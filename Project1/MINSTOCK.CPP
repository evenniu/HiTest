#include "stdafx.h"
#include "minstock.h"

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

CMinStock::CMinStock(void)
{
  m_nom = NULL;;
  m_mea = NULL;;
  m_nMCL = NULL;;
  m_mMCL = NULL;;

  m_numChecks = 0;

  m_nomCV = NULL;
  m_nomCC = NULL;
  m_meaCV = NULL;
  m_meaCC = NULL;
  m_stock = NULL;
  m_nomT = NULL;
  m_side = NULL;

  m_alignCV = NULL;
  m_alignCC = NULL;
  m_alignLE = NULL;
  m_alignTE = NULL;

}

CMinStock::~CMinStock(void)
{
  if (m_nomCV)
    delete m_nomCV;

  if (m_nomCC)
    delete m_nomCC;

  if (m_meaCV)
    delete m_meaCV;

  if (m_meaCC)
    delete m_meaCC;

  if (m_stock)
    delete [] m_stock;

  if (m_nomT)
    delete [] m_nomT;

  if (m_side)
    delete [] m_side;
}

bool CMinStock::Initialize(double dLE, double d, double dTE)
{
  if (d < 0.001)
    return false;

  if (dLE < 0.0001)
    dLE = d;

  if (dTE < 0.0001)
    dTE = d;

  m_delta = d;
  m_deltaLE = dLE;
  m_deltaTE = dTE;

  if (!m_nom || !m_mea || !m_nMCL || !m_mMCL)
    return false;

  m_mclNomLERatio = 0.001;// myGetProfileDouble(L"NomMclLEBackoff", 0.001);
  m_mclNomTERatio = 0.001;//myGetProfileDouble(L"NomMclTEBackoff", 0.001);
  m_mclMeaLERatio = 0.004;// myGetProfileDouble(L"MeaMclLEBackoff", 0.004);
  m_mclMeaTERatio = 0.004;//myGetProfileDouble(L"MeaMclTEBackoff", 0.004);
  
  double dt = m_nMCL->T1() - m_nMCL->T0();

  m_nom->CalcPoint(m_nomLEChd, m_nomleext);
  m_nMCL->CalcPoint(m_nomLECtr, m_nMCL->T0() + m_mclNomLERatio*dt);

  m_nom->CalcPoint(m_nomTEChd, m_nomteext);
  m_nMCL->CalcPoint(m_nomTECtr, m_nMCL->T1() - m_mclNomTERatio*dt);

  m_nomLength = dist(m_nomLEChd, m_nomLECtr) +
    m_nMCL->TotalLength(1000, m_nMCL->T0() + m_mclNomLERatio*dt, m_nMCL->T1() - m_mclNomTERatio*dt) +
    dist(m_nomTECtr, m_nomTEChd);

  int n;

  // compute how many points will be needed. add some extras just in case.

  n = (int) (m_leoff/m_deltaLE + (m_nomLength - m_leoff - m_teoff)/m_delta + m_teoff/m_deltaTE) + 10;

  m_nomCV = new CMatrix(n, 2);
  if (!m_nomCV)
    return false;

  m_nomCC = new CMatrix(n, 2);
  if (!m_nomCC)
    return false;

  m_meaCV = new CMatrix(n, 2);
  if (!m_meaCV)
    return false;

  m_meaCC = new CMatrix(n, 2);
  if (!m_meaCC)
    return false;

  m_stock = new double[n];
  if (!m_stock)
    return false;

  m_nomT = new double[n];
  if (!m_nomT)
    return false;

  m_side = new int[n];
  if (!m_side)
    return false;

  m_numChecks = n;

  m_pTol[0] = m_pTol[1] = m_pTol[2] = m_pTol[3] = 0.0;
  m_mTol[0] = m_mTol[1] = m_mTol[2] = m_mTol[3] = 0.0;
  m_calculated[0] = m_calculated[1] = m_calculated[2] = m_calculated[3] = false;

  return true;
}

void CMinStock::SetCurves(CCurve *nom, CCurve *mea, CCurve *nMCL, CCurve *mMCL,
                          double nomleext, double nomteext, double mealeext, double meateext, double leoff, double teoff)
{
  m_nom = nom;
  m_mea = mea;
  m_nMCL = nMCL;
  m_mMCL = mMCL;

  m_nomPeriod = nom->Period();
  m_meaPeriod = mea->Period();

  m_nomleext = nomleext;
  m_nomteext = nomteext;
  m_mealeext = mealeext;
  m_meateext = meateext;

  m_leoff = leoff;
  m_teoff = teoff;
}


bool CMinStock::GetNomCV(int i, double *xy)
{
  if (i < 0 || i >= m_numChecks)
    return false;

  xy[0] = m_nomCV->m[i][0];
  xy[1] = m_nomCV->m[i][1];

  return true;
}

bool CMinStock::GetNomCC(int i, double *xy)
{
  if (i < 0 || i >= m_numChecks)
    return false;

  xy[0] = m_nomCC->m[i][0];
  xy[1] = m_nomCC->m[i][1];

  return true;
}

bool CMinStock::GetMeaCV(int i, double *xy)
{
  if (i < 0 || i >= m_numChecks)
    return false;

  xy[0] = m_meaCV->m[i][0];
  xy[1] = m_meaCV->m[i][1];

  return true;
}

bool CMinStock::GetMeaCC(int i, double *xy)
{
  if (i < 0 || i >= m_numChecks)
    return false;

  xy[0] = m_meaCC->m[i][0];
  xy[1] = m_meaCC->m[i][1];

  return true;
}

double CMinStock::GetStock(int i)
{
  if (i >= 0 && i < m_numChecks)
    return m_stock[i];
  else
    return 0.0;
}

void CMinStock::SetPTol(int i, double tol)
{
  if (i >= 0 && i < 4)
    m_pTol[i] = tol;
}

void CMinStock::SetMTol(int i, double tol)
{
  if (i >= 0 && i < 4)
    m_mTol[i] = tol;
}

double CMinStock::GetPTol(int i)
{
  return (i >= 0 && i < 4) ? m_pTol[i] : 0.0;
}

double CMinStock::GetMTol(int i)
{
  return (i >= 0 && i < 4) ? m_mTol[i] : 0.0;
}

bool CMinStock::Calculate(int side, CAlignment *align, double startVal, double stopVal)
{
  if (m_numChecks < 1)
    return false;

  m_calculated[side] = true;

  m_mea->Align(align);
  m_mMCL->Align(align);

  m_firstAlign = align;
  m_alignCV = align;
  m_alignCC = align;
  m_alignLE = align;
  m_alignTE = align;

  double dt = m_mMCL->T1() - m_mMCL->T0();

  m_mea->CalcPoint(m_meaLEChd, m_mealeext);
  m_mMCL->CalcPoint(m_meaLECtr, m_mMCL->T0() + m_mclMeaLERatio*dt);

  m_mea->CalcPoint(m_meaTEChd, m_meateext);
  m_mMCL->CalcPoint(m_meaTECtr, m_mMCL->T1() - m_mclMeaTERatio*dt);

  if (startVal > 0.0)
    dt = startVal;
  else
    dt = 2*m_deltaLE;

  if (stopVal  <= 0.0)
    stopVal = 2*m_deltaTE;

  double nomCV[2], nomCC[2], meaCV[2], meaCC[2], nomMC[2];
  m_minStock = m_minStockLE = m_minStockTE = 1.0e20;
  m_maxStock = m_maxStockLE = m_maxStockTE = -1.0e20;
  m_minIndex = m_minIndexLE = m_minIndexTE = -1;
  m_maxIndex = m_maxIndexLE = m_maxIndexTE = -1;

  int i;
  bool jumpOut = false;
  for (i=0; i < m_numChecks; i++)
  {
    CalcStock(i, dt, nomCV, nomCC, meaCV, meaCC, nomMC);
    m_nomT[i] = dt;

    m_nomCV->m[i][0] = nomCV[0];
    m_nomCV->m[i][1] = nomCV[1];
    m_nomCC->m[i][0] = nomCC[0];
    m_nomCC->m[i][1] = nomCC[1];
    m_meaCV->m[i][0] = meaCV[0];
    m_meaCV->m[i][1] = meaCV[1];
    m_meaCC->m[i][0] = meaCC[0];
    m_meaCC->m[i][1] = meaCC[1];

    m_stock[i] = dist(meaCV, meaCC) - dist(nomCV, nomCC);

    double dLE = dist(m_nomLEChd, nomMC);
    double dTE = dist(m_nomTEChd, nomMC);
    if (dLE < m_leoff)
    {
      if (m_stock[i] < m_minStockLE)
      {
        m_minStockLE = m_stock[i];
        m_minIndexLE = i;
      }
      if (m_stock[i] > m_maxStockLE)
      {
        m_maxStockLE = m_stock[i];
        m_maxIndexLE = i;
      }
      m_side[i] = LEC;
    }
    else if (dTE < m_teoff)
    {
      if (m_stock[i] < m_minStockTE)
      {
        m_minStockTE = m_stock[i];
        m_minIndexTE = i;
      }
      if (m_stock[i] > m_maxStockTE)
      {
        m_maxStockTE = m_stock[i];
        m_maxIndexTE = i;
      }
      m_side[i] = TEC;
    }
    else
    {
      if (m_stock[i] < m_minStock)
      {
        m_minStock = m_stock[i];
        m_minIndex = i;
      }
      if (m_stock[i] > m_maxStock)
      {
        m_maxStock = m_stock[i];
        m_maxIndex = i;
      }
      m_side[i] = CVC;
    }

    if (jumpOut)
      break;

    if (dt < m_leoff)
      dt += m_deltaLE;
    else if (dt < m_nomLength - m_teoff)
    {
      if (dt + m_delta > m_nomLength - m_teoff)
        dt = m_nomLength - m_teoff + 0.5*m_deltaTE;
      else
        dt += m_delta;
    }
    else
      dt += m_deltaTE;

    if (dt > m_nomLength - stopVal)
    {
      dt = m_nomLength - stopVal;
      jumpOut = true;
    }
  }

  if (i < m_numChecks)
    m_numChecks = i+1;

  m_mea->Align(NULL);
  m_mMCL->Align(NULL);

  return false;
}


bool CMinStock::ReCalculate(int side, CAlignment *align)
{
  if (side != LEC && side != TEC)
    return false;

  if (align == m_firstAlign)
    return true;

  if (m_calculated[side])
    return true;

  m_calculated[side] = true;

  if (side == CVC || side == CCC)
  {
    m_alignCV = align;
    m_alignCC = align;
  }
  if (side == LEC)
    m_alignLE = align;
  if (side == TEC)
    m_alignTE = align;

  m_mea->Align(align);
  m_mMCL->Align(align);

  double dt = m_mMCL->T1() - m_mMCL->T0();

  m_mea->CalcPoint(m_meaLEChd, m_mealeext);
  m_mMCL->CalcPoint(m_meaLECtr, m_mMCL->T0() + m_mclMeaLERatio*dt);

  m_mea->CalcPoint(m_meaTEChd, m_meateext);
  m_mMCL->CalcPoint(m_meaTECtr, m_mMCL->T1() - m_mclMeaTERatio*dt);

  double nomCV[2], nomCC[2], meaCV[2], meaCC[2], nomMC[2];

  if (side == LEC)
  {
    m_minStockLE = 1.0e20;
    m_maxStockLE = -1.0e20;
    m_minIndexLE = -1;
    m_maxIndexLE = -1;
  }
  else if (side == TEC)
  {
    m_minStockTE = 1.0e20;
    m_maxStockTE = -1.0e20;
    m_minIndexTE = -1;
    m_maxIndexTE = -1;
  }
  else
  {
    m_minStock = 1.0e20;
    m_maxStock = -1.0e20;
    m_minIndex = -1;
    m_maxIndex = -1;
  }

  int i;
  for (i=0; i<m_numChecks; i++)
  {
    if (m_side[i] != side)
      continue;

    CalcStock(i, m_nomT[i], nomCV, nomCC, meaCV, meaCC, nomMC);

    m_nomCV->m[i][0] = nomCV[0];
    m_nomCV->m[i][1] = nomCV[1];
    m_nomCC->m[i][0] = nomCC[0];
    m_nomCC->m[i][1] = nomCC[1];
    m_meaCV->m[i][0] = meaCV[0];
    m_meaCV->m[i][1] = meaCV[1];
    m_meaCC->m[i][0] = meaCC[0];
    m_meaCC->m[i][1] = meaCC[1];

    m_stock[i] = dist(meaCV, meaCC) - dist(nomCV, nomCC);

    if (side == LEC)
    {
      if (m_stock[i] < m_minStockLE)
      {
        m_minStockLE = m_stock[i];
        m_minIndexLE = i;
      }
      if (m_stock[i] > m_maxStockLE)
      {
        m_maxStockLE = m_stock[i];
        m_maxIndexLE = i;
      }
      m_side[i] = LEC;
    }
    else if (side == TEC)
    {
      if (m_stock[i] < m_minStockTE)
      {
        m_minStockTE = m_stock[i];
        m_minIndexTE = i;
      }
      if (m_stock[i] > m_maxStockTE)
      {
        m_maxStockTE = m_stock[i];
        m_maxIndexTE = i;
      }
      m_side[i] = TEC;
    }
    else
    {
      if (m_stock[i] < m_minStock)
      {
        m_minStock = m_stock[i];
        m_minIndex = i;
      }
      if (m_stock[i] > m_maxStock)
      {
        m_maxStock = m_stock[i];
        m_maxIndex = i;
      }
      m_side[i] = CVC;
    }
  }

  m_mea->Align(NULL);
  m_mMCL->Align(NULL);

  return true;
}

double CMinStock::CalcStock(int /*i*/, double t, double *nomCV, double *nomCC, double *meaCV, double *meaCC, double *nomMC)
{
  if (m_numChecks < 1)
    return 0.0;

  double nd1 = dist(m_nomLEChd, m_nomLECtr);
  double nd3 = dist(m_nomTEChd, m_nomTECtr);
  double nd2 = m_nomLength - nd1 - nd3;

  if (t < 0 || t > m_nomLength)
    return 0.0;

  double dt, dir[2];
  double tStart, tEnd, dAdj, tAdj;

  if (t < nd1)  // on the first straight part
  {
    dir[0] = m_nomLECtr[0] - m_nomLEChd[0];
    dir[1] = m_nomLECtr[1] - m_nomLEChd[1];
    normalize(dir, dir);
    nomMC[0] = m_nomLEChd[0] + t*dir[0];
    nomMC[1] = m_nomLEChd[1] + t*dir[1];
  }
  else if (t < nd1 + nd2)  // on the MCL
  {
    dt = m_nMCL->T1() - m_nMCL->T0();
    tStart = m_nMCL->T0() + m_mclNomLERatio*dt;
    tEnd = m_nMCL->T1() - m_mclNomLERatio*dt;
    dAdj = (t - nd1);
    tAdj = tStart + dAdj*(tEnd - tStart)/nd2;
    m_nMCL->CalcPoint(nomMC, tAdj, dir);
  }
  else  // on the last straight part
  {
    dir[0] = m_nomTECtr[0] - m_nomTEChd[0];
    dir[1] = m_nomTECtr[1] - m_nomTEChd[1];
    normalize(dir, dir);
    nomMC[0] = m_nomTEChd[0] + (m_nomLength - t)*dir[0];
    nomMC[1] = m_nomTEChd[1] + (m_nomLength - t)*dir[1];
  }

  double nrm[2], cvt0, cvt1, cct0, cct1;  // doesn't really matter in all of this which side is CCV and which is CVX
  nrm[0] = dir[1];
  nrm[1] = -dir[0];
  cvt0 = m_nomleext;
  cvt1 = m_nomteext;
  if (cvt1 < cvt0)
    cvt1 += m_nomPeriod;
  cct0 = m_nomteext;
  cct1 = m_nomleext;
  if (cct1 < cct0)
    cct1 += m_nomPeriod;
  m_nom->LineIntersect(nomMC, nrm, nomCV, cvt0, cvt1);
  m_nom->LineIntersect(nomMC, nrm, nomCC, cct0, cct1);

  double bestT, meaMC[2];

  dir[0] = m_meaLECtr[0] - m_meaLEChd[0];
  dir[1] = m_meaLECtr[1] - m_meaLEChd[1];
  dt = normalize(dir, dir);
  bestT =  projdist(m_meaLEChd, dir, nomMC, meaMC);
  if (bestT > dt)
  {
    dir[0] = m_meaTEChd[0] - m_meaTECtr[0];
    dir[1] = m_meaTEChd[1] - m_meaTECtr[1];
    dt = normalize(dir, dir);
    bestT =  projdist(m_meaTECtr, dir, nomMC, meaMC);
    if (bestT < 0)
    {
      dt = m_mMCL->T1() - m_mMCL->T0();
      m_mMCL->ClosestPoint(nomMC, meaMC, &bestT, dir, m_mMCL->T0() + m_mclMeaLERatio*dt, m_mMCL->T1() - m_mclMeaTERatio*dt);
    }
  }

  nrm[0] = dir[1];
  nrm[1] = -dir[0];
  cvt0 = m_mealeext;
  cvt1 = m_meateext;
  if (cvt1 < cvt0)
    cvt1 += m_meaPeriod;
  cct0 = m_meateext;
  cct1 = m_mealeext;
  if (cct1 < cct0)
    cct1 += m_meaPeriod;
  m_mea->LineIntersect(meaMC, nrm, meaCV, cvt0, cvt1);
  m_mea->LineIntersect(meaMC, nrm, meaCC, cct0, cct1);

  return 0.0;
}
