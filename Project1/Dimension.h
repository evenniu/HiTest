// Dimension.h: interface for the CDimension class.

#pragma once

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

class DLLEXPORT CDimension
{
public:
  int m_type;      // DIM_xxx
  int m_oldtype;   // just used by editor
  bool m_default;  // from default section
  bool m_olddefault;  // used by editor

  double m_nom, m_ptol, m_mtol, m_usl, m_lsl;
  int m_methodOverride;      // defines which method to use for some calcs.
  wchar_t m_header1OverrideArr[50]; // defines which method to use for some calcs.
  wchar_t m_header2OverrideArr[50]; // defines which method to use for some calcs.
 
  //vector<string> m_header1OverrideVec;
  //vector<string> m_header2OverrideVec;
  double m_RoadTHCKRotateOverride;
  double m_RoadTHCKRotateXaxisOverride;
  double m_RoadTHCKoffsetOverride;
  CDimension();
  ~CDimension();
  void Set();                    // build spec limits from tols or vice versa
  bool isDefined();
  void SetTolerances(double mtol, double ptol, bool isSpecLimits);

  void Copy(const CDimension& d, bool makeDefault = false);
};

