// ToleranceFile.cpp: implementation of the CToleranceFile class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "ToleranceSection.h"
#include <Eigen/Dense>

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

ComplexZone::ComplexZone()
{
  Init();
}

void ComplexZone::Init()
{
  numberOfZoneRegions = 0;
  interpolationType = ComplexZoneInterpolationType::SMOOTH;
  Eigen::Map<Eigen::VectorXd>(zoneDistances, MAX_SPECIAL_ZONE_REGIONS).setZero();
  Eigen::Map<Eigen::VectorXd>(zoneMinusTolerances, MAX_SPECIAL_ZONE_REGIONS).setZero();
  Eigen::Map<Eigen::VectorXd>(zonePlusTolerances, MAX_SPECIAL_ZONE_REGIONS).setZero();
}

void ComplexZone::Copy(const ComplexZone& inputZone)
{
  numberOfZoneRegions = inputZone.numberOfZoneRegions;
  interpolationType = inputZone.interpolationType;
  Eigen::Map<Eigen::VectorXd>(zoneDistances, MAX_SPECIAL_ZONE_REGIONS)
      .lazyAssign(Eigen::Map<const Eigen::VectorXd>(inputZone.zoneDistances, MAX_SPECIAL_ZONE_REGIONS));
  Eigen::Map<Eigen::VectorXd>(zoneMinusTolerances, MAX_SPECIAL_ZONE_REGIONS)
      .lazyAssign(
          Eigen::Map<const Eigen::VectorXd>(inputZone.zoneMinusTolerances, MAX_SPECIAL_ZONE_REGIONS));
  Eigen::Map<Eigen::VectorXd>(zonePlusTolerances, MAX_SPECIAL_ZONE_REGIONS)
      .lazyAssign(
          Eigen::Map<const Eigen::VectorXd>(inputZone.zonePlusTolerances, MAX_SPECIAL_ZONE_REGIONS));
}

CToleranceSection::CToleranceSection()
{
  Init();
}

CToleranceSection::~CToleranceSection()
{
}

void CToleranceSection::Init()
{
  m_name[0] = 0;

  for(int i = 0; i < 10; i++)
  {
    m_leOffset[i] = m_teOffset[i] = -1.0e20;
    m_leFlatSize[i] = m_teFlatSize[i] = -1.0e20;
    m_leAngle[i] = m_teAngle[i] = -1.0e20;
    m_lePoints[i][0] = m_lePoints[i][1] =  m_tePoints[i][0] = m_tePoints[i][1] = -1;
    m_complexEdgeZone[i].Init();
    m_chordZone[i].Init();

    m_leCuppingOffset[i] = -1.0e20;
    m_teCuppingOffset[i] = -1.0e20;
    m_teCuppingOffsetB[i] = -1.0e20;
    m_teCuppingOffsetC[i] = -1.0e20;
  }

  m_leChange = m_teChange = -1.0e20;
  m_leChange2 = m_teChange2 = -1.0e20;
  m_leAngOffset[0] = m_leAngOffset[1] = m_teAngOffset[0] = m_teAngOffset[1] = -1.0e20;
  m_concavePositionReferenceAngle = -1.0e20;
  m_leRadii[0] = m_leRadii[1] = -1.0e20;
  m_teRadii[0] = m_teRadii[1] = -1.0e20;
  m_camber2016_offsets[0] = m_camber2016_offsets[1] = std::numeric_limits<double>::quiet_NaN();
  m_leVarTol[0] = m_leVarTol[1] = m_teVarTol[0] = m_teVarTol[1] = -1.0e20;
  m_midVarTol[0] = m_midVarTol[1] = m_chordVarTol[0] = m_chordVarTol[1] = -1.0e20;
  m_gageParam[0] = m_gageParam[1] = -1.0e20;
  m_extremeAng[0] = m_extremeAng[1] = m_extremeAng[2] = m_extremeAng[3] = -1.0e20;
  m_stockDist = m_stockDistLE = m_stockDistTE = -1.0e20;
  m_waveWidth = -1.0e20;

  m_leZoneMult = m_teZoneMult = 1.0;

  m_changeOverride[0] = m_changeOverride[1] = -1;

  m_inspectionSection = false;

  for(int i = 0; i < ParamArraySize; i++) // set to false, so default section doesn't get set to true
    m_default[i] = false;       // copy will initial to true for other sections

  for(int i = 0; i < MAXFITS; i++)
    m_exclude[i] = 0;

  for(int i = 0; i < CalcArraySize; i++)
  {
    m_offset[i] = 0.0;
  }
}

void CToleranceSection::Copy(const CToleranceSection* tol, bool makeDefault)
{
  for(int i = 0; i < 10; i++)
  {
    m_leOffset[i] = tol->m_leOffset[i];
    m_teOffset[i] = tol->m_teOffset[i];
    m_leFlatSize[i] = tol->m_leFlatSize[i];
    m_teFlatSize[i] = tol->m_teFlatSize[i];
    m_leAngle[i] = tol->m_leAngle[i];
    m_teAngle[i] = tol->m_teAngle[i];
    m_lePoints[i][0] = tol->m_lePoints[i][0];
    m_lePoints[i][1] = tol->m_lePoints[i][1];
    m_tePoints[i][0] = tol->m_tePoints[i][0];
    m_tePoints[i][1] = tol->m_tePoints[i][1];
    m_complexEdgeZone[i].Copy(tol->m_complexEdgeZone[i]);
    m_chordZone[i].Copy(tol->m_chordZone[i]);

    m_leCuppingOffset[i] = tol->m_leCuppingOffset[i];
    m_teCuppingOffset[i] = tol->m_teCuppingOffset[i];
    m_teCuppingOffsetB[i] = tol->m_teCuppingOffsetB[i];
    m_teCuppingOffsetC[i] = tol->m_teCuppingOffsetC[i];
  }

  for(int i = 0; i < MAXFITS; i++)
    m_exclude[i] = tol->m_exclude[i];

  m_leAngOffset[0] = tol->m_leAngOffset[0];
  m_leAngOffset[1] = tol->m_leAngOffset[1];
  m_teAngOffset[0] = tol->m_teAngOffset[0];
  m_teAngOffset[1] = tol->m_teAngOffset[1];
  m_concavePositionReferenceAngle = tol->m_concavePositionReferenceAngle;
  m_leRadii[0] = tol->m_leRadii[0];
  m_leRadii[1] = tol->m_leRadii[1];
  m_teRadii[0] = tol->m_teRadii[0];
  m_teRadii[1] = tol->m_teRadii[1];
  m_camber2016_offsets[0] = tol->m_camber2016_offsets[0];
  m_camber2016_offsets[1] = tol->m_camber2016_offsets[1];
  m_leVarTol[0] = tol->m_leVarTol[0];
  m_leVarTol[1] = tol->m_leVarTol[1];
  m_teVarTol[0] = tol->m_teVarTol[0];
  m_teVarTol[1] = tol->m_teVarTol[1];
  m_midVarTol[0] = tol->m_midVarTol[0];
  m_midVarTol[1] = tol->m_midVarTol[1];
  m_chordVarTol[0] = tol->m_chordVarTol[0];
  m_chordVarTol[1] = tol->m_chordVarTol[1];
  m_leChange = tol->m_leChange;
  m_teChange2 = tol->m_teChange2;
  m_leChange2 = tol->m_leChange2;
  m_teChange = tol->m_teChange;
  m_changeOverride[0] = tol->m_changeOverride[0];
  m_changeOverride[1] = tol->m_changeOverride[1];

  m_gageParam[0] = tol->m_gageParam[0];
  m_gageParam[1] = tol->m_gageParam[1];
  m_extremeAng[0] = tol->m_extremeAng[0];
  m_extremeAng[1] = tol->m_extremeAng[1];
  m_extremeAng[2] = tol->m_extremeAng[2];
  m_extremeAng[3] = tol->m_extremeAng[3];

  m_stockDist = tol->m_stockDist;
  m_stockDistLE = tol->m_stockDistLE;
  m_stockDistTE = tol->m_stockDistTE;

  m_waveWidth = tol->m_waveWidth;

  m_leZoneMult = tol->m_leZoneMult;
  m_teZoneMult = tol->m_teZoneMult;

  m_inspectionSection = makeDefault ? false : tol->m_inspectionSection;

  for (int i = 0; i < CalcArraySize; i++)
  {
    m_dim[i].Copy(tol->m_dim[i], makeDefault);
    m_dim[i].m_RoadTHCKRotateOverride = tol->m_dim[i].m_RoadTHCKRotateOverride;//nyc 2021-11-08
    m_dim[i].m_RoadTHCKRotateXaxisOverride = tol->m_dim[i].m_RoadTHCKRotateXaxisOverride;
    m_dim[i].m_RoadTHCKoffsetOverride = tol->m_dim[i].m_RoadTHCKoffsetOverride;
    wcscpy_s(m_dim[i].m_header1OverrideArr, tol->m_dim[i].m_header1OverrideArr);
    wcscpy_s(m_dim[i].m_header2OverrideArr, tol->m_dim[i].m_header2OverrideArr);
    

    m_offset[i] = tol->m_offset[i];
  }
  for(int i = 0; i < ParamArraySize; i++)
    m_default[i] = makeDefault ? true : tol->m_default[i];
}

void CToleranceSection::SetName(wchar_t *name)
{
  wcscpy_s(m_name, name);
}

// Originally I used m_leAngOffset instead of creating a new parameter.  Use that value if they haven't defined m_leRadii

void CToleranceSection::GetLERadii(double *r1, double *r2)
{
  if (m_leRadii[0] > -1.0e19 && m_leRadii[1] > -1.0e19)
  {
    *r1 = m_leRadii[0];
    *r2 = m_leRadii[1];
  }
  else
  {
    *r1 = m_leAngOffset[0];
    *r2 = m_leAngOffset[1];
  }

  if (*r1 > *r2)  // swap them
  {
    double t = *r1;
    *r1 = *r2;
    *r2 = t;
  }
}

void CToleranceSection::GetTERadii(double *r1, double *r2)
{
  if (m_teRadii[0] > -1.0e19 && m_teRadii[1] > -1.0e19)
  {
    *r1 = m_teRadii[0];
    *r2 = m_teRadii[1];
  }
  else
  {
    *r1 = m_teAngOffset[0];
    *r2 = m_teAngOffset[1];
  }

  if (*r1 > *r2)  // swap them
  {
    double t = *r1;
    *r1 = *r2;
    *r2 = t;
  }
}
