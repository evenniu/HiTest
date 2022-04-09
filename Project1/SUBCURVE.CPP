#include "StdAfx.h"
#include <stdio.h>
#include "SubCurve.h"

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

// during reporting, it can be useful to retrieve information about the curve;
// I call these annotations
void CSubCurve::addDoubleAnnotation(const char* /*key*/, double /*value*/)
{
}
double CSubCurve::getDoubleAnnotation(const char* /*key*/) const
{
  return std::numeric_limits<double>::quiet_NaN();
}
void CSubCurve::addBoolAnnotation(const char* /*key*/, bool /*value*/)
{
}
bool CSubCurve::getBoolAnnotation(const char* /*key*/) const
{
  return false;
}

CSubCurve::CSubCurve(CCurve *c, double s, double e, double p) : CCurve()
{
  m_psb = c;
  m_st0 = s;
  m_st1 = e;
  m_period = p;
  m_isEnglish = c->IsEnglish();

  m_min[0] = m_min[1] = 1.0;
  m_max[0] = m_max[1] = 0.0;
}

/*
CSubCurve::CSubCurve(FILE *, int)
{
m_min[0] = m_min[1] = 1.0;
m_max[0] = m_max[1] = 0.0;
}
*/
int CSubCurve::CalcPoint(double *xyz, double t, double *tan, double *curv) const 
{
  if (!m_psb)
    return 0;

  if (t > m_st1)
    t = m_st1;
  if (t < m_st0)
    t = m_st0;

  if (t > m_psb->T1())
    t -= m_period;
  if (t < m_psb->T0())
    t += m_period;

  int rv = m_psb->CalcPoint(xyz, t, tan, curv);
  if (rv && m_palign)
  {
    m_palign->MeasToBest(xyz, 1, xyz);
    if (tan)
      m_palign->MeasToBest(tan, 0, tan);
    if (curv)
      m_palign->MeasToBest(curv, 0, curv);
  }

  return rv;
}

void CSubCurve::Write(FILE * /*fp*/)
{}

void CSubCurve::Extent()
{
  if (m_min[0] <= m_max[0] && m_min[1] <= m_max[1])
    return;

  int i,j;
  double t,td,xyz[3];

  td = 0.01*(m_st1 - m_st0);

  for (i=j=0, t=m_st0; i<100; i++, t += td)
  {
    if (CalcPoint(xyz, t))
    {
      if (j++ == 0)
      {
        m_min[0] = m_max[0] = xyz[0];
        m_min[1] = m_max[1] = xyz[1];
      }
      else
      {
        if (xyz[0] < m_min[0])
          m_min[0] = xyz[0];
        if (xyz[1] < m_min[1])
          m_min[1] = xyz[1];

        if (xyz[0] > m_max[0])
          m_max[0] = xyz[0];
        if (xyz[1] > m_max[1])
          m_max[1] = xyz[1];
      }
    }
  }
}

void CSubCurve::Extent(double *min, double *max)
{
  Extent();

  min[0] = m_min[0];
  min[1] = m_min[1];
  max[0] = m_max[0];
  max[1] = m_max[1];
}

int CSubCurve::Extreme(double *ij, double *t, double *p, double t0, double t1, int save, bool endOK)
{
  if (t0 == 0.0 && t1 == 0.0)
  {
    t0 = m_st0;
    t1 = m_st1;
  }

  if (m_psb)
  {
    int rv = m_psb->Extreme(ij, t, p, t0, t1, save, endOK);
    if (save == 1)
      m_extreme = *t;
    return rv;
  }

  return 0;
}
