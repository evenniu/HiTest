// ToleranceSection.h: interface for the CToleranceSection class.

#pragma once

#include "Dimension.h"

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

enum class ComplexZoneInterpolationType : int32_t
{
  CONSTANT,
  LINEAR,
  SMOOTH
};

struct DLLEXPORT ComplexZone
{
  int numberOfZoneRegions;
  ComplexZoneInterpolationType interpolationType;
  double zoneDistances[MAX_SPECIAL_ZONE_REGIONS];
  double zoneMinusTolerances[MAX_SPECIAL_ZONE_REGIONS];
  double zonePlusTolerances[MAX_SPECIAL_ZONE_REGIONS];

  ComplexZone();
  ~ComplexZone()
  {
  }
  void Init();
  void Copy(const ComplexZone& inputZone);
};


class DLLEXPORT CToleranceSection
{
public:

  void SetName(wchar_t *name);

  wchar_t m_name[50];     // Section name

  bool m_inspectionSection;

  int m_lePoints[10][2], m_tePoints[10][2];
  double m_leOffset[10], m_teOffset[10], m_leChange, m_teChange, m_leChange2, m_teChange2;
  double m_leFlatSize[10], m_teFlatSize[10];
  double m_leAngle[10], m_teAngle[10];  // for ANGLEOFFSET method for widths.
  double m_gageParam[2];
  double m_leAngOffset[2], m_teAngOffset[2];
  double m_leRadii[2], m_teRadii[2];
  double m_camber2016_offsets[2];
  double m_leVarTol[2], m_teVarTol[2], m_midVarTol[2], m_chordVarTol[2];
  double m_leZoneMult, m_teZoneMult;
  double m_extremeAng[4];
  double m_stockDistLE, m_stockDist, m_stockDistTE;
  double m_waveWidth;
  int m_exclude[MAXFITS];
  int m_changeOverride[2];
  ComplexZone m_complexEdgeZone[10];
  ComplexZone m_chordZone[10];
  double m_leCuppingOffset[10], m_teCuppingOffset[10], m_teCuppingOffsetB[10], m_teCuppingOffsetC[10];
  double m_concavePositionReferenceAngle;

  CDimension m_dim[CalcArraySize];
  double m_offset[CalcArraySize];
  bool m_default[ParamArraySize];  // from default sections?

  CToleranceSection();
  ~CToleranceSection();
  void Init();
  void Copy(const CToleranceSection* tol, bool makeDefault = false);
  void GetLERadii(double *r1, double *r2);
  void GetTERadii(double *r1, double *r2);
};
