#include "StdAfx.h"
#include "bestfit.h"
#include <assert.h>
#include <Eigen/Dense>

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

CBestFit::CBestFit(int np) : m_numpts(np), m_vals(0), m_infs(0), m_noms(0), m_ijks(0), m_omega(0), /*m_b(0), m_Q(0), m_Q1(0),
m_xbar(0), m_gcurve(0),*/ m_btxpt(0), m_btypt(0), m_outtol(0), m_bestPartOf(0), m_vartol(0), m_LEvartol(0), m_t(0)
{
  if (np > 0)
  {
    m_vals = new CMatrix(np,2);
    m_infs = new CMatrix(np,2);
    m_noms = new CMatrix(np,2);
    m_ijks = new CMatrix(np,2);
    m_t = new CVector(np);
    m_omega = new CVector(np);
    //m_b = new CVector(np);
    //m_xbar = new CVector(np);

    m_btxpt = new double [np];
    m_btypt = new double [np];
    m_outtol = new double [np];
    m_bestPartOf = new int [np];
    m_valWasUsedInFit = new bool[np];
  }

  m_reportFit = 0;
  m_totalBad = 0;
  m_totalChecked = 0;

  //m_Q = new CMatrix(2, 2);
  //m_Q1 = new CMatrix(2, 2);
}

CBestFit::~CBestFit()
{
  if (m_vals)
    delete m_vals;

  if (m_infs)
    delete m_infs;

  if (m_noms)
    delete m_noms;

  if (m_ijks)
    delete m_ijks;

  if (m_omega)
    delete m_omega;

  if (m_t)
    delete m_t;

  //if (m_b)
  //  delete m_b;

  //if (m_Q)
  //  delete m_Q;

  //if (m_Q1)
  //  delete m_Q1;

  //if (m_xbar)
  //  delete m_xbar;

  if (m_outtol)
    delete [] m_outtol;

  if (m_btxpt)
    delete [] m_btxpt;

  if (m_btypt)
    delete [] m_btypt;

  if (m_bestPartOf)
    delete []m_bestPartOf;

  if (m_vartol)
    delete m_vartol;

  if (m_LEvartol)
    delete m_LEvartol;

  if(m_valWasUsedInFit)
    delete[] m_valWasUsedInFit;
}

// LeastSquaresFit: Find rigid body motion (rotation and translation) that minimizes the sum of the squared distances
//                  between nominal points and measured points

bool CBestFit::LeastSquaresFit(CFitParams &fp, double *nompt, double *meapt)
{
  int i,j;

  Eigen::Map<const Eigen::VectorXd> eigenOmega(m_omega->v, m_numpts);
  /*if(eigenOmega.squaredNorm() == 0.0)
  {
    throw std::logic_error("No weight left at all!");
  }*/

  CVector Xmean(2);
  CVector Zmean(2);

  if (nompt && meapt)
  {
    Xmean.v[0] = meapt[0];   // set translation points
    Xmean.v[1] = meapt[1];
    Zmean.v[0] = nompt[0];
    Zmean.v[1] = nompt[1];
  }
  else
  {
    switch (fp.tranfit)
    {
    case 0:  // full translation
      // get the weighted centroids
      centroid(m_numpts, 2, m_vals->m, Xmean.v, m_omega->v);
      centroid(m_numpts, 2, m_noms->m, Zmean.v, m_omega->v);
      if(!std::isfinite(Xmean.v[0]) || !std::isfinite(Xmean.v[1]))
      {
        throw std::logic_error("Non-numbers crept into Xmean!");
      }
      if(!std::isfinite(Zmean.v[0]) || !std::isfinite(Zmean.v[1]))
      {
        throw std::logic_error("Non-numbers crept into Zmean!");
      }
      break;
    case 1:  // no translation
      Xmean.v[0] = 0.0;
      Xmean.v[1] = 0.0;
      Zmean.v[0] = 0.0;
      Zmean.v[1] = 0.0;
      break;
    case 2:  // in tol
      centroid(m_numpts, 2, m_vals->m, Xmean.v, m_omega->v);
      centroid(m_numpts, 2, m_noms->m, Zmean.v, m_omega->v);
      if(!std::isfinite(Xmean.v[0]) || !std::isfinite(Xmean.v[1]))
      {
        throw std::logic_error("Non-numbers crept into Xmean!");
      }
      if(!std::isfinite(Zmean.v[0]) || !std::isfinite(Zmean.v[1]))
      {
        throw std::logic_error("Non-numbers crept into Zmean!");
      }

      // rectangular tolerance zones

      for (i=0; i<2; i++)
      {
        /* backwards???*/
        if (Xmean.v[i] - Zmean.v[i] > fp.tranPTol[i])
          Xmean.v[i] =  Zmean.v[i] + fp.tranPTol[i];
        else if (Xmean.v[i] - Zmean.v[i] < fp.tranMTol[i])
          Xmean.v[i] =  Zmean.v[i] + fp.tranMTol[i];
        /**/
        if (Xmean.v[i] - Zmean.v[i] > -fp.tranMTol[i])
          Xmean.v[i] =  Zmean.v[i] - fp.tranMTol[i];
        else if (Xmean.v[i] - Zmean.v[i] < -fp.tranPTol[i])
          Xmean.v[i] =  Zmean.v[i] - fp.tranPTol[i];

        if(!std::isfinite(Xmean.v[0]) || !std::isfinite(Xmean.v[1]))
        {
          throw std::logic_error("Non-numbers crept into Xmean!");
        }
        if(!std::isfinite(Zmean.v[0]) || !std::isfinite(Zmean.v[1]))
        {
          throw std::logic_error("Non-numbers crept into Zmean!");
        }
      }

      break;
    }
  }

  double Theta = 0.0;

  switch (fp.rotfit)
  {
  case 0:  // full
  case 2:  // intol
    {
      CMatrix x(m_numpts,2);
      CMatrix z(m_numpts,2);

      for (i = 0; i < m_numpts; i++)
      {
        for (j = 0; j < 2; j++)
        {
          x.m[i][j] = m_vals->m[i][j] - Xmean.v[j];
          z.m[i][j] = m_noms->m[i][j] - Zmean.v[j];
        }
      }

      double num = 0.0, den = 0.0;

      for (i=0; i<m_numpts; i++)
      {
        num += m_omega->v[i]*(z.m[i][1]*x.m[i][0]-z.m[i][0]*x.m[i][1]);
        den += m_omega->v[i]*(z.m[i][0]*x.m[i][0]+z.m[i][1]*x.m[i][1]);
      }

      if (fabs(den) < 1.0e-12)
        return false;

      Theta = atan(num/den);
      if (fp.rotfit == 2) // intol
      {
       // int angleunits = myGetProfileInt(L"AngleFormat", 1);
        int angleunits = 1;
        //double ptol = fabs(M_PI*fp.rotPTol/180);
        //double mtol = -fabs(M_PI*fp.rotMTol/180);
        double ptol = M_PI*fp.rotPTol/180;
        double mtol = M_PI*fp.rotMTol/180;

        if (angleunits == 0)
        {
          ptol /= 60;
          mtol /= 60;
        }

        if (Theta > ptol)
          Theta = ptol;
        else if (Theta < mtol)
          Theta = mtol;
      }
    }
    break;
  }

  m_align.m_borig[0] = Zmean.v[0];
  m_align.m_borig[1] = Zmean.v[1];
  m_align.m_morig[0] = Xmean.v[0];   // set translation points
  m_align.m_morig[1] = Xmean.v[1];

  m_align.m_mat[0][0] = m_align.m_mat[1][1] = cos(Theta);
  m_align.m_mat[1][0] = sin(Theta);
  m_align.m_mat[0][1] = -m_align.m_mat[1][0];

  for (i=0; i<m_numpts; i++)
  {
    double meas[2],bfit[2];

    meas[0] = m_vals->m[i][0];
    meas[1] = m_vals->m[i][1];

    m_align.MeasToBest(meas,1,bfit);

    if(!std::isfinite(bfit[0]) || !std::isfinite(bfit[1]))
    {
      throw std::logic_error("Non-numbers crept into bfit!");
    }
    m_infs->m[i][0] = bfit[0];
    m_infs->m[i][1] = bfit[1];
  }

  return true;
}

bool CBestFit::Omega()                  // fill m_omega array
{
  for (int i=0; i<m_numpts; i++)
    m_omega->v[i] = 1.0;

  return true;
}

bool CBestFit::Omega(int i, double v)  // fill m_omega value
{
  if (i < 0  || i >= m_numpts)
    return false;

  m_omega->v[i] = v;

  return true;
}

void CBestFit::PutVal(int i, double x, double y)   // copy vector to m_vals array
{
  if (i < 0 || i >= m_numpts)
    return;

  m_vals->m[i][0] = x;
  m_vals->m[i][1] = y;
}

bool CBestFit::GetVal(int i, double *v)    // copy vector from m_vals array
{
  if (i < 0 || i >= m_numpts)
    return false;

  v[0] = m_vals->m[i][0];
  v[1] = m_vals->m[i][1];

  return true;
}

bool CBestFit::GetInf(int i, double *v)    // copy vector from m_infs array
{
  if (i < 0 || i >= m_numpts)
    return false;

  v[0] = m_infs->m[i][0];
  v[1] = m_infs->m[i][1];

  return true;
}

void CBestFit::PutInf(int i, double x, double y)   // copy vector to m_infs array
{
  if (i < 0 || i >= m_numpts)
    return;

  if(!std::isfinite(x) || !std::isfinite(y))
  {
    throw std::logic_error("Non-numbers crept into (x, y)!");
  }
  m_infs->m[i][0] = x;
  m_infs->m[i][1] = y;
}


void CBestFit::PutNom(int i, double x, double y)   // copy vector to m_noms array
{
  if (i < 0 || i >= m_numpts)
    return;

  m_noms->m[i][0] = x;
  m_noms->m[i][1] = y;
}

void CBestFit::PutT(int i, double t)   // copy vector to m_noms array
{
  if (i < 0 || i >= m_numpts)
    return;

  m_t->v[i] = t;
}

bool CBestFit::GetNom(int i, double *v)    // copy vector from m_noms array
{
  if (i < 0 || i >= m_numpts)
    return false;

  v[0] = m_noms->m[i][0];
  v[1] = m_noms->m[i][1];

  return true;
}

bool CBestFit::GetT(int i, double *t)    // copy vector from m_noms array
{
  if (i < 0 || i >= m_numpts)
    return false;

  *t = m_t->v[i];

  return true;
}

void CBestFit::PutVec(int i, double x, double y)   // copy vector to m_noms array
{
  if (i < 0 || i >= m_numpts)
    return;

  m_ijks->m[i][0] = x;
  m_ijks->m[i][1] = y;
}

bool CBestFit::GetVec(int i, double *v)    // copy vector from m_ijks array
{
  if (i < 0 || i >= m_numpts)
    return false;

  v[0] = m_ijks->m[i][0];
  v[1] = m_ijks->m[i][1];

  return true;
}

bool CBestFit::Align(CAlignment *a)      // copy alignment out
{
  if (!a)
    return false;

  for (int i=0; i<2; i++)
  {
    a->m_morig[i] = m_align.m_morig[i];
    a->m_borig[i] = m_align.m_borig[i];
    for (int j=0; j<2; j++)
      a->m_mat[i][j] = m_align.m_mat[i][j];
  }

  return true;
}

bool CBestFit::ImportAlign(CAlignment *a)      // copy alignment in
{
  if (!a)
    return false;

  for (int i=0; i<2; i++)
  {
    m_align.m_morig[i] = a->m_morig[i];
    m_align.m_borig[i] = a->m_borig[i];
    for (int j=0; j<2; j++)
      m_align.m_mat[i][j] = a->m_mat[i][j];
  }

  return true;
}

// function to find global (hopefully) minimum of a function, zfunc
// zfunc is a pointer which should be set before calling

#define GOLDEN 0.6180339887499   // golden ratio

int CBestFit::ApplyVarTol(double *minOut, double *maxOut)
{
  *minOut = *maxOut = 0.0;

  if (!m_vartol)
    return -1;

  int numout = 0;
  for (int p=0; p<m_numpts; p++)
  {
    m_btxpt[p] = m_btypt[p] = m_outtol[p] = 0.0;

    double d, mea[2], tol[2];

    mea[0] = m_infs->m[p][0];  // m_bxpt[p];
    mea[1] = m_infs->m[p][1];  // m_bypt[p];

    if (m_vartol->OutTol(mea, tol, &d, p))
    {
      if (d < *minOut)
        *minOut = d;
      if (d > *maxOut)
        *maxOut = d;
      m_btxpt[p] = tol[0];
      m_btypt[p] = tol[1];
      m_outtol[p] = d;
      numout++;
    }
  }

  return numout;
}

int CBestFit::ReturnVarTol(double *minOut, double *maxOut)
{
  *minOut = *maxOut = 0.0;

  if (!m_vartol)
    return -1;

  int numout = 0;
  for (int p=0; p<m_numpts; p++)
  {
    double d = m_outtol[p];

    if (d != 0.0)
    {
      if (d < *minOut)
        *minOut = d;
      if (d > *maxOut)
        *maxOut = d;
      numout++;
    }
  }

  return numout;
}

int CBestFit::ApplyLEVarTol(double *minOut, double *maxOut)
{
  *minOut = *maxOut = 0.0;

  if (!m_LEvartol)
    return -1;

  int numout = 0;
  for (int p=0; p<m_numpts; p++)
  {
    m_btxpt[p] = m_btypt[p] = m_outtol[p] = 0.0;

    if (m_bestPartOf[p] != LEC)
      continue;

    double d, mea[2], tol[2];

    mea[0] = m_infs->m[p][0];  // m_bxpt[p];
    mea[1] = m_infs->m[p][1];  // m_bypt[p];

    if (m_LEvartol->OutTol(mea, tol, &d, p))
    {
      if (d < *minOut)
        *minOut = d;
      if (d > *maxOut)
        *maxOut = d;
      m_btxpt[p] = tol[0];
      m_btypt[p] = tol[1];
      m_outtol[p] = d;
      numout++;
    }
  }

  return numout;
}

int CBestFit::ReturnLEVarTol(double *minOut, double *maxOut)
{
  *minOut = *maxOut = 0.0;

  if (!m_LEvartol)
    return -1;

  int numout = 0;
  for (int p=0; p<m_numpts; p++)
  {
    if (m_bestPartOf[p] != LEC)
      continue;

    double d = m_outtol[p];

    if (d != 0.0)
    {
      if (d < *minOut)
        *minOut = d;
      if (d > *maxOut)
        *maxOut = d;
      numout++;
    }
  }

  return numout;
}

bool CBestFit::ReturnFit(double *x, double *y, double *theta, int method)
{
  if (method < 0)
    method = m_reportFit;
  else if (method > 1)
    method = 0;

  if (method == 0)  // old method - relative to pivot point
  {
    *x     = m_align.m_borig[0] - m_align.m_morig[0];
    *y     = m_align.m_borig[1] - m_align.m_morig[1];
    *theta = 180.0*asin(m_align.m_mat[1][0])/M_PI;  // note this is degrees
  }
  else              // must be new method - relative to stack point
  {
    //*x     = m_align.m_borig[0] - m_align.m_morig[0];
    //*y     = m_align.m_borig[1] - m_align.m_morig[1];
    //*theta = 180.0*asin(m_align.m_mat[1][0])/M_PI;  // note this is degrees
    bugout(0, L"%f %f %f would be returned by method 0", *x, *y, *theta);

    double xy[2] = {0.0, 0.0};
    m_align.BestToMeas(xy, 1, xy);
    *x = -xy[0];                           // minus to agree with how has been done in the past
    *y = -xy[1];

    int bi = -1;
    double maxD = -1.0;
    int i;
    for (i=0; i<m_numpts; i++)
    {
      double d = dist(m_vals->m[i], xy);
      if (d > maxD)
      {
        bi = i;
        maxD = d;
      }
    }
    if (bi < 0)
      return false;

    double vm[2], vn[2];

    vm[0] = m_vals->m[bi][0] - xy[0];
    vm[1] = m_vals->m[bi][1] - xy[1];
    vn[0] = m_infs->m[bi][0];
    vn[1] = m_infs->m[bi][1];

    normalize(vm, vm);
    normalize(vn, vn);

    *theta = 180.0*acos(dot(vm, vn))/M_PI;
    if (curl(vn, vm) > 0.0)                 // backwards to agree with way done in past
      *theta *= -1.0;
    //bugout(0, _T("%f %f %f is returned by method 1"), *x, *y, *theta);
  }

  return true;
}

CFitParams::CFitParams()
{
  for (int i=0; i<4; i++)
  {
    fitcurve[i] = 1;
    weightcurve[i] = 1;
  }
  algorithm = BestFitAlgorithm::LeastSquares;
  rotfit = 0;
  tranfit = 0;
  fitToMiddleOfZone = 0;
  usenominals = 0;
  rotMTol = 0.0;
  rotPTol = 0.0;
  profileMTol = 0.0;
  profilePTol = 0.0;
  pivot = 1;
  tranPTol[0] = tranPTol[1] = 0.0;
  tranMTol[0] = tranMTol[1] = 0.0;
  lepercent = 5.0;
  tepercent = 95.0;
  leoff1 = 0.0;
  leoff2 = 0.0;
  teoff1 = 0.0;
  teoff2 = 0.0;
  forgeAngle = 0.0;

  complexEdgeZone = nullptr;
  chordZone = nullptr;
}

// routines used by "Manipulate" view

void CBestFit::XYDelta(double *xy)
{
  CAlignment *align = GetAlign();
  align->m_morig[0] -= xy[0];
  align->m_morig[1] -= xy[1];
  Align(align);

  int i;
  for(i=0; i<m_numpts; i++)
  {
    double pxy[2];
    pxy[0] = m_vals->m[i][0];
    pxy[1] = m_vals->m[i][1];
    align->MeasToBest(pxy, 1, pxy);
    PutInf(i, pxy[0], pxy[1]);
  }
}

void CBestFit::RotDelta(double angle)
{
  CAlignment *align = GetAlign();
  align->Rotate(2, angle);
  Align(align);

  int i;
  for(i=0; i<m_numpts; i++)
  {
    double xy[2];
    xy[0] = m_vals->m[i][0];
    xy[1] = m_vals->m[i][1];
    align->MeasToBest(xy, 1, xy);
    PutInf(i, xy[0], xy[1]);
  }
}