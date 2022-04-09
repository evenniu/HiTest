#pragma once

#include "Alignment.h"
#define NS_OK              0
#define NS_BADNOSEINDICES  1
#define NS_BADTAILINDICES  2

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

namespace Hexagon
{
namespace Blade
{
struct MeanCamberCurveParameters2016;
}
}

class DLLEXPORT CNominalSection
{
public:

  wchar_t m_name[MAXBUFSZ];

  bool m_english;
  bool m_startCC;  // only used during calculation
  bool m_finishLE;
  bool m_finishTE;
  int m_convexSide;

  int m_npts;
  int m_errorCode;
  int m_ni1, m_ni2, m_ti1, m_ti2;
  int m_bladetype;
  int m_havek;
  int m_haveTol;
  int m_inose;
  int m_itail;
  const Hexagon::Blade::MeanCamberCurveParameters2016 *m_mclParams;
  int m_noseindex[2];
  int m_tailindex[2];
  int m_noseforce;
  int m_tailforce;
  int m_leSquare;
  int m_teSquare;
  int m_nPhantomLE;
  int m_nPhantomTE;
  int m_skewReport;

  short m_skewed;
  short m_letype;
  short m_tetype;

  double m_zval;
  double m_period;
  double m_nose[3];
  double m_tail[3];
  double m_t0[4];
  double m_t1[4];
  double m_pitch[3];
  double m_change[2];
  double m_skew[3];
  double m_skewOrigin[3];
  double m_ler;  // used in calculation only
  double m_ter;
  double m_tnose;
  double m_ttail;
  double m_pxLE[19], m_pyLE[19];   // phantom points for partial LE
  double m_pxTE[19], m_pyTE[19];   // phantom points for partial LE
  double m_endKLE[4], m_endKTE[4];
  double m_origLE[2][2];  // may need to restore original points to get written to math file
  double m_origTE[2][2];

  double *m_ox;  // original ordering
  double *m_oy;
  double *m_oz;
  double *m_oi;
  double *m_oj;
  double *m_ok;
  double *m_omtol;
  double *m_optol;

  double *m_xv;
  double *m_yv;
  double *m_kv;
  int m_tolsegCount;            //公差段个数
  //std::vector<string> m_tolsegNameVec; //公差段名称集合
  int *m_start_point;
  int *m_end_point;
  int* m_curvature_coef;
  double *m_pTol_start;
  double *m_pTol_end;
  double *m_mTol_start;
  double *m_mTol_end;
  BOOL m_tolCCW;
  CAlignment *m_align;

  CCurve *m_whole;
  CCurve *m_mc;
  CCurve *m_cvc;
  CCurve *m_ccc;
  CCurve *m_lec;
  CCurve *m_tec;

  CBladeError *m_error;     // last error holder
  double m_areaStartPoint[3][2]; //开口区域起始点坐标
  double m_areaEndPoint[3][2];  //开口区域终止点坐标
  int m_areaIndex[3][2];  //开口区域起始-终止点索引
  int m_areaCount;  //开口区域个数,最大3个
  CNominalSection();
  ~CNominalSection();

  int ApplySkew(double* zaxis, double* origin);
  void SetSize(int n, int needk, int needTol, bool english);
  int Calculate(double** m, double* ns, double* tl, double* leoff, double* teoff, int useextr, bool* vecWarn,
                const Hexagon::Blade::MeanCamberCurveParameters2016* mclParams);
  int Calculate(double** m, int ni1, int ni2, int ti1, int ti2, int nosetype, int tailtype, double* leoff,
                double* teoff, int useextr, const Hexagon::Blade::MeanCamberCurveParameters2016* mclParams);
  void MakeTolSegStore(int num);

private:

  void AdjustSubCurves(CCurve** tcvc, CCurve** tccc, CCurve** tlec, CCurve** ttec);
  void DetermineNormalNoseTail(double** m, double* ns, double* tl, double* leoff, double* teoff);
  void NormalReorderPoints(double** m, bool* vecWarn);
  int BuildNormalSubCurves(double* ns, double* tl, double leoff, double teoff, int forceoffsets, double* leo,
                           double* teo, int useextr, const Hexagon::Blade::MeanCamberCurveParameters2016* mclParams);
  void InitialNormalNoseTailTValues(double leoff, double teoff);
  void RefineNormalNoseTailTValues(double leoff, double teoff, int forceoffsets, double* leo, double* teo);
  int BuildNormalMeanCamber(double* ns, double* tl, int useextr,
                            const Hexagon::Blade::MeanCamberCurveParameters2016* mclParams);
  int BuildSquarePartialMeanCamber(int useextr, const Hexagon::Blade::MeanCamberCurveParameters2016* mclParams);
  int DetermineSquarePartialEndPoints(double** m, int ntype, int ttype, double* tooClose);
  int FillArrayComputeSplineSquarePartial(double** m, double* leoff, double* teoff, double tooClose);
  int DetermineTransitionSquarePartial(double leoff, double teoff);
};
