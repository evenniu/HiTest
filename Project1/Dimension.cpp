#include "stdafx.h"
#include "Dimension.h"

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

CDimension::CDimension()
{
  m_type = m_oldtype = DIM_SKIP;
  m_default = m_olddefault = true;
  m_nom = -1.0e20;
  m_ptol = m_mtol = m_usl = m_lsl = 0.0;
  m_methodOverride = -1;
  m_RoadTHCKRotateOverride = 0.0;
  m_RoadTHCKRotateXaxisOverride = 0.0;
  m_RoadTHCKoffsetOverride = 0.0;
  for(int i = 0; i < 50;i++)
  {
    m_header1OverrideArr[i] = { '\0' };
    m_header2OverrideArr[i] = { '\0' };
  }

}

CDimension::~CDimension()
{
}

void CDimension::SetTolerances(double mtol, double ptol, bool isSpecLimits)
{
  if (isSpecLimits)
  {
    m_lsl = mtol;
    m_usl = ptol;
    m_type = DIM_SPEC;
  }
  else
  {
    m_mtol = mtol;
    m_ptol = ptol;
    m_type = DIM_TOL;
  }
}

void CDimension::Set()
{
  if (m_nom < -1.0e19)
    return;

  if (m_type == DIM_TOL)
  {
    m_usl = m_nom + m_ptol;
    m_lsl = m_nom + m_mtol;
  }
  else if (m_type == DIM_SPEC)
  {
    m_ptol = m_usl - m_nom;
    m_mtol = m_lsl - m_nom;
  }
}

bool CDimension::isDefined()
{
  switch(m_type)
  {
  case DIM_TOL:
  case DIM_SPEC:
    break;
  case DIM_NOMTOL:
  case DIM_NOMSPEC:
    if (m_nom < -1.0e19)
      return false;
    break;
  default:
    return false;
  }

  return true;
}

void CDimension::Copy(const CDimension& d, bool makeDefault)
{
  m_type = d.m_type;
  m_oldtype = m_type;
  m_default = makeDefault ? true : d.m_default;
  m_olddefault = m_default;

  m_nom = d.m_nom;
  m_ptol = d.m_ptol;
  m_mtol = d.m_mtol;
  m_usl = d.m_usl;
  m_lsl = d.m_lsl;
  m_methodOverride = d.m_methodOverride;
}
