#include "StdAfx.h"
#include "NominalSection.h"
#include "Circle.h"
#include "Nurb.h"
#include "SubCurve.h"
#include "SmallestCircle.h"
#include "Sorter.h"
#include "MeanCamberCurve.h"
#include "CurvePolygon.h"
#include "HermiteCurve.h"

using namespace Hexagon;
using namespace Blade;


#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

CNominalSection::CNominalSection()
{
  m_name[0] = 0;

  m_convexSide = 0;

  m_npts = 0;
  m_letype = EDGE_NORMAL;
  m_tetype = EDGE_NORMAL;
  m_skewed = 0;
  m_skewReport = 0;

  m_zval = 0.0;
  m_period = -1.0;
  m_pitch[2] = -1.0;

  m_xv = 0;
  m_yv = 0;
  m_kv = 0;
  m_ox = 0;
  m_oy = 0;
  m_oz = 0;
  m_oi = 0;
  m_oj = 0;
  m_ok = 0;

  m_mclParams = nullptr;
  m_whole = 0;
  m_mc = 0;
  m_cvc = 0;
  m_ccc = 0;
  m_lec = 0;
  m_tec = 0;

  m_align = 0;
  m_omtol = NULL;
  m_optol = NULL;
  m_tolCCW = false;
  //m_error->ClearErrors();
  //m_tolsegNameVec = new std::vector<string>();
  m_tolsegCount = 0;
  m_start_point = NULL;
  m_end_point = NULL;
  m_curvature_coef = NULL;
  m_pTol_start = NULL;
  m_pTol_end = NULL;
  m_mTol_start = NULL;
  m_mTol_end = NULL;
  m_areaCount = 0;
  for(int i = 0; i < 3; i++)
  {
    m_areaStartPoint[i][0] = 0.0;
    m_areaStartPoint[i][1] = 0.0;
    m_areaEndPoint[i][0] = 0.0;
    m_areaEndPoint[i][1] = 0.0;
    m_areaIndex[i][0] = 0;
    m_areaIndex[i][1] = 0;
  }
}

CNominalSection::~CNominalSection()
{
  if(m_error != nullptr)
  {
    delete m_error;
    m_error = nullptr;
  }
  if (m_xv)
    delete [] m_xv;
  if (m_yv)
    delete [] m_yv;
  if (m_kv)
    delete [] m_kv;

  if (m_ox)
    delete [] m_ox;
  if (m_oy)
    delete [] m_oy;
  if (m_oz)
    delete [] m_oz;
  if (m_oi)
    delete [] m_oi;
  if (m_oj)
    delete [] m_oj;
  if (m_ok)
    delete [] m_ok;

  if (m_omtol)
    delete [] m_omtol;
  if (m_optol)
    delete [] m_optol;

  if(m_mclParams)
    delete m_mclParams;

  if (m_whole)
    delete m_whole;
  if (m_mc)
    delete m_mc;
  if (m_cvc)
    delete m_cvc;
  if (m_ccc)
    delete m_ccc;
  if (m_lec)
    delete m_lec;
  if (m_tec)
    delete m_tec;

  if (m_align)
    delete m_align;
   if(m_start_point)
    delete[] m_start_point;
   if(m_end_point)
    delete[] m_end_point;  
   
   if(m_curvature_coef)
     delete[] m_curvature_coef;
   if(m_pTol_start)
   delete [] m_pTol_start;
   if(m_pTol_end)
    delete[] m_pTol_end;
   if(m_mTol_start)
    delete[] m_mTol_start;
   if(m_mTol_end)
    delete[] m_mTol_end;
}

void CNominalSection::SetSize(int n, int needk, int needTol, bool english)
{
  m_npts = n;
  m_havek = needk;
  m_haveTol = needTol;
  m_english = english;
}

int CNominalSection::ApplySkew(double *zaxis, double *origin)
{
  m_skewed = 1;
  m_skew[0] = zaxis[0];
  m_skew[1] = zaxis[1];
  m_skew[2] = zaxis[2];
  m_skewOrigin[0] = origin[0];
  m_skewOrigin[1] = origin[1];
  m_skewOrigin[2] = origin[2];

  m_align = new CAlignment(3);
  m_align->MatrixFromAxis(m_skew);

  m_align->m_morig[0] = origin[0];
  m_align->m_morig[1] = origin[1];
  m_align->m_morig[2] = origin[2];

  m_align->m_borig[0] = origin[0];
  m_align->m_borig[1] = origin[1];
  m_align->m_borig[2] = origin[2];

  return 1;
}

void CNominalSection::DetermineNormalNoseTail(double **m, double *ns, double *tl, double *leoff, double *teoff)
{
  int i;
  
  if (ns && tl)   // m_nose and m_tail points are passed in
  {
    m_nose[0] = ns[0];
    m_nose[1] = ns[1];
    m_tail[0] = tl[0];
    m_tail[1] = tl[1];
    double mind = 1.0e20;
    for (i=0; i<m_npts; i++)
    {
      double d = dist(m_nose, m[i]);
      if (d < mind)
      {
        mind = d;
        m_inose = i;
      }
    }
    mind = 1.0e20;
    for (i=0; i<m_npts; i++)
    {
      double d = dist(m_tail, m[i]);
      if (d < mind)
      {
        mind = d;
        m_itail = i;
      }
    }
  }
  else if (ns)   // m_nose is specified but not m_tail
  {
    m_nose[0] = ns[0];
    m_nose[1] = ns[1];

    double bestd = 1.0e20;
    for (i=0; i<m_npts; i++)
    {
      double d = dist(m_nose, m[i]);
      if (d < bestd)
      {
        bestd = d;
        m_inose = i;
      }
    }

    bestd = 1.0e20;
    int besti = -1;
    int lim1 = (m_inose + 3*m_npts/8)%m_npts;
    double *rads = new double[m_npts/4];
    int *inds = new int[m_npts/4];
    for (i=0; i<m_npts/4; i++)
    {
      int j = (lim1 + i)%m_npts;
      inds[i] = j+1;
      rads[i] = radius_3points(m[j], m[(j+1)%m_npts], m[(j+2)%m_npts]);
      //bugout(0, _T("10 %f %f %d=%f"), m[j][0], m[j][1], j+1, rads[i]);
      if (rads[i] < bestd)
      {
        bestd = rads[i];
        besti = i;
        m_itail = j+1;
      }
    }

    // look for start and end indices that "close" to the smallest radius  say 10 percent;

    int sind = -1, eind = -1;
    for (i=0; i<m_npts/4; i++)
    {
      if (rads[i]/rads[besti] < 1.1)
      {
        if (sind < 0)
          sind = inds[i];
        eind = inds[i];
      }
    }

    if (eind < sind)
      eind += m_npts;

    m_itail = (sind + eind)/2;
    if (m_itail >= m_npts)
      m_itail -= m_npts;

    //bugout(0, _T(" s %d e %d b %d"), sind, eind, m_itail);
    delete [] rads;
    delete [] inds;

    m_tail[0] = m[m_itail][0];
    m_tail[1] = m[m_itail][1];
  }
  else
  {
    double *rad = new double[m_npts];
    double bestd = 1.0e20;
    for (i=0; i<m_npts; i++)
    {
      rad[i] = radius_3points(m[i], m[(i+1)%m_npts], m[(i+2)%m_npts]);
      if (rad[i] < bestd)
      {
        bestd = rad[i];
        m_itail = i;
      }
    }

    bestd = 1.0e20;
    int lim1 = (m_itail + 3*m_npts/8)%m_npts;
    for (i=0; i<m_npts/4; i++)
    {
      int j = (lim1 + i)%m_npts;
      if (rad[j] < bestd)
      {
        bestd = rad[j];
        m_inose = j;
      }
    }

    delete [] rad;

    m_nose[0] = m[m_inose][0];
    m_nose[1] = m[m_inose][1];
    m_tail[0] = m[m_itail][0];
    m_tail[1] = m[m_itail][1];
  }

  //bugout(0, _T("10 %f %f N %d first guess"), m_nose[0], m_nose[1], m_inose);
  //bugout(0, _T("10 %f %f T %d first guess"), m_tail[0], m_tail[1], m_itail);
  // Step 3. Fit circles and check size to determine which is m_nose and m_tail
  //         only swap if m_nose was not defined

  double cx[5], cy[5];

  int j;
  for (i=0; i<5; i++)
  {
    j = (m_inose + m_npts-2 + i) % m_npts;
    cx[i] = m[j][0];
    cy[i] = m[j][1];
  }
  CCircle c0(5, cx, cy);
  m_ler = c0.Rad();

  for (i=0; i<5; i++)
  {
    j = (m_itail + m_npts-2 + i) % m_npts;
    cx[i] = m[j][0];
    cy[i] = m[j][1];
  }
  CCircle c1(5, cx, cy);
  m_ter = c1.Rad();

  if (!ns)
  {
    if (m_ler < m_ter)  // want m_inose to be leading edge
    {
      j = m_inose;
      m_inose = m_itail;
      m_itail = j;
      m_ler = c1.Rad();
      m_ter = c0.Rad();
    }

    m_nose[0] = m[m_inose][0];
    m_nose[1] = m[m_inose][1];
    m_tail[0] = m[m_itail][0];
    m_tail[1] = m[m_itail][1];
  }
  
  // adjust end offsets if needed
  double clen = dist(m_nose, m_tail);
  if (m_ler < 0.25*clen && *leoff < 3.0*m_ler)
    *leoff = 3.0*m_ler;
  if (m_ter < 0.25*clen && *teoff < 3.0*m_ter)
    *teoff = 3.0*m_ter;
  
  //bugout(4, _T("leoff %f teoff %f"), leoff, teoff);
  //if (leo) bugout(4, _T("leo %f"), *leo); else bugout(4, _T("leo undefined"));
  //if (teo) bugout(4, _T("teo %f"), *teo); else bugout(4, _T("teo undefined"));
}

void CNominalSection::NormalReorderPoints(double **m, bool *vecWarn)
{
  // changing the log here, from what it used to be, need to search out to find a point in the + direction and in the - direction
  // that is about 5% out from the LE toward the TE
  // was using projdist, but this didn't work  for a blade with a large amount of curve

  double len = dist(m[m_inose], m[m_itail]);
  int i, j;
  int jp = -1, jm = -1;
  for (i=1; i<m_npts-1; i++)
  {
    j = (m_inose + i)%m_npts;
    double d = dist(m[m_inose], m[j]);
    if (d > 0.05*len)
    {
      jp = j;
      break;
    }
  }
  for (i=1; i<m_npts-1; i++)
  {
    j = (m_inose + m_npts - i)%m_npts;
    double d = dist(m[m_inose], m[j]);
    if (d > 0.05*len)
    {
      jm = j;
      break;
    }
  }

  double v1[2], v2[2];
  j = (m_inose + m_npts/8)%m_npts;
  if (jp >= 0 && jm >= 0)
    j = jp;
  v1[0] = m[j][0] - m[m_inose][0];
  v1[1] = m[j][1] - m[m_inose][1];

  j = (m_inose + m_npts - m_npts/8)%m_npts;
  if (jp >= 0 && jm >= 0)
    j = jm;
  v2[0] = m[j][0] - m[m_inose][0];
  v2[1] = m[j][1] - m[m_inose][1];

  int cw = 0;
  if (curl(v1, v2) < 0.0) // points are clock wise
    cw = 1;

  // Step 5. Fill array, ordered CCW, point 0 in middle of side before LE

  m_xv = new double[m_npts];
  m_yv = new double[m_npts];
  m_ox = new double[m_npts];
  m_oy = new double[m_npts];
  m_oz = new double[m_npts];
  if (m_havek)
  {
    m_kv = new double[m_npts];
    m_oi = new double[m_npts];
    m_oj = new double[m_npts];
    m_ok = new double[m_npts];
  }
  if (m_haveTol)
  {
    m_omtol = new double[m_npts];
    m_optol = new double[m_npts];
  }


  for (i=0; i<m_npts; i++)
  {
    //if (cw)
     j = i% m_npts;
    //else
    //  j = (m_npts - zeroindex + i) % m_npts;

    //thisNormal[0] = m[i][3];
    //thisNormal[1] = m[i][4];
    //thisNormal[2] = m[i][5];
    //normalize(thisNormal, thisNormal);

    //if (i > 0)
    //{
    //  if (dot(thisNormal, lastNormal) < 0.0) // more than 90 change between adjacent points
    //  {
    //    //bugout(0,L"10 %lf %lf VECWARN", m[i][0], m[i][1]);
    //    *vecWarn = true;
    //  }
    //}

    //lastNormal[0] = thisNormal[0];
    //lastNormal[1] = thisNormal[1];
    //lastNormal[2] = thisNormal[2];

    m_xv[j] = m[i][0];
    m_yv[j] = m[i][1];
    //bugout(0, L"10 %f %f %d was %d", m_xv[j], m_yv[j], j, i);
    if (m_havek)
      m_kv[j] = m[i][5];

    m_ox[j] = m[j][0];
    m_oy[j] = m[j][1];
    m_oz[j] = m[j][2];

    if (m_havek)
    {
      m_oi[j] = m[j][3];
      m_oj[j] = m[j][4];
      m_ok[j] = m[j][5];
    }

    if (m_haveTol)
    {
      m_omtol[j] = m[j][6];
      m_optol[j] = m[j][7];
    }
  }
}

void CNominalSection::InitialNormalNoseTailTValues(double leoff, double teoff)
{
  double ltan[3], ttan[3], bp1[2], bp2[2];
  m_whole->ClosestPoint(m_nose, bp1, &m_tnose, ltan, 0.0, m_period/2.0, 400);
  m_whole->ClosestPoint(m_tail, bp2, &m_ttail, ttan, m_period/2.0, m_period, 400);

  // determine reasonable break off for LE and CV/CC transition

  m_whole->CircIntersect(m_nose, 1.5*leoff, bp1, 0.0, m_tnose, &m_t0[LEC]);
  m_whole->CircIntersect(m_nose, 1.5*leoff, bp2, m_tnose, m_period/2, &m_t1[LEC]);
  
  //bugout(0, _T("8 %f %f %f %f pass 0 le start-end"),bp1[0],bp1[1],bp2[0],bp2[1]);
  //bugout(0, _T("pass 0 calculated m_t0[LEC] %f m_t1[LEC] %f m_period %f"), m_t0[LEC], m_t1[LEC], m_period);
  //bugout(0, _T("10 %f %f LE1_%d"), bp1[0], bp1[1], pass);
  //bugout(0, _T("10 %f %f LE2_%d"), bp2[0], bp2[1], pass);

  // determine reasonable break off for TE and CV/CC transition

  m_whole->CircIntersect(m_tail, 1.5*teoff, bp1, m_period/2, m_ttail, &m_t0[TEC]);
  m_whole->CircIntersect(m_tail, 1.5*teoff, bp2, m_ttail, m_period, &m_t1[TEC]);
  //bugout(0, _T("8 %f %f %f %f pass 0 te start-end"),bp1[0],bp1[1],bp2[0],bp2[1]);
  //bugout(0, _T("10 %f %f TE1_%d"), bp1[0], bp1[1], pass);
  //bugout(0, _T("10 %f %f TE2_%d"), bp2[0], bp2[1], pass);
}

void CNominalSection::RefineNormalNoseTailTValues(double leoff, double teoff, int forceoffsets, double* leo,
                                                  double* teo)
{
  double lijk[3], tijk[3], mclle[2], mclte[2];

  double mclNomRatioLE = 0.001;// myGetProfileDouble(L"NomMclLEBackoff", 0.001);
  double mclNomRatioTE = 0.001;// myGetProfileDouble(L"NomMclTEBackoff", 0.001);

  double deltaLE = mclNomRatioLE * (m_mc->T1() - m_mc->T0());
  double deltaTE = mclNomRatioTE * (m_mc->T1() - m_mc->T0());

  m_mc->CalcPoint(mclle, m_mc->T0() + deltaLE, lijk);
  normalize(lijk, lijk);
  lijk[0] *= -1.0;
  lijk[1] *= -1.0;
  m_lec->NewLineIntersect(mclle, lijk, m_nose, 0.0, 0.0, &m_tnose, NULL, true);

  m_mc->CalcPoint(mclte, m_mc->T1() - deltaTE, tijk);
  // tijk[0] *= -1.0;
  // tijk[1] *= -1.0;
  normalize(tijk, tijk);
  m_tec->NewLineIntersect(mclte, tijk, m_tail, 0.0, 0.0, &m_ttail, NULL, true);

  if(forceoffsets)
  {
    leoff = *leo; // set back to what they were
    teoff = *teo;
  }

  double bcp[3], tan[3], nrm[3], bp1[2], bp2[2];

  // determine requested break off for LE and CV/CC transition

  if(!m_mc->CircIntersect(m_nose, leoff, bcp, 0.0, 0.0, 0, tan))
  {
    tan[0] = mclle[0] - m_nose[0];
    tan[1] = mclle[1] - m_nose[1];
    normalize(tan, tan);
    bcp[0] = m_nose[0] + leoff * tan[0];
    bcp[1] = m_nose[1] + leoff * tan[1];
  }
  nrm[0] = tan[1];
  nrm[1] = -tan[0];
  normalize(nrm, nrm);
  m_whole->NewLineIntersect(bcp, nrm, bp1, 0.0, m_tnose, &m_t0[LEC]);
  m_whole->NewLineIntersect(bcp, nrm, bp2, m_tnose, m_period / 2, &m_t1[LEC]);
  // bugout(0, _T("8 %f %f %f %f pass 1 le start-end"),bp1[0],bp1[1],bp2[0],bp2[1]);
  // bugout(0, _T("pass 1 calculated m_t0[LEC] %f m_t1[LEC] %f m_period %f"), m_t0[LEC], m_t1[LEC], m_period);

  // determine requested break off for TE and CV/CC transition

  if(!m_mc->CircIntersect(m_tail, teoff, bcp, 0.0, 0.0, 0, tan))
  {
    tan[0] = mclte[0] - m_tail[0];
    tan[1] = mclte[1] - m_tail[1];
    normalize(tan, tan);
    bcp[0] = m_tail[0] + teoff * tan[0];
    bcp[1] = m_tail[1] + teoff * tan[1];
  }

  nrm[0] = tan[1];
  nrm[1] = -tan[0];
  normalize(nrm, nrm);
  m_whole->NewLineIntersect(bcp, nrm, bp1, m_period / 2, m_ttail, &m_t0[TEC]);
  m_whole->NewLineIntersect(bcp, nrm, bp2, m_ttail, m_period, &m_t1[TEC]);
  // bugout(0, _T("8 %f %f %f %f pass 1 te start-end"),bp1[0],bp1[1],bp2[0],bp2[1]);

  if(m_convexSide > 0 && m_convexSide < 5)
  {
    double pt1[2], vec[2], pt2[2];

    pt1[0] = m_xv[0];
    pt1[1] = m_yv[0];

    vec[0] = vec[1] = 0;
    if(m_convexSide == 1)
      vec[0] = 1.0;
    else if(m_convexSide == 2)
      vec[0] = -1.0;
    else if(m_convexSide == 3)
      vec[1] = 1.0;
    else if(m_convexSide == 4)
      vec[1] = -1.0;

    if(m_whole->NewLineIntersect(pt1, vec, pt2, m_t1[LEC], m_t0[TEC]))
    {
      if(projdist(pt1, vec, pt2) < 0.0) // the first point is on the Convex side
        m_startCC = false;
      else
        m_startCC = true;
    }
  }
}

int CNominalSection::BuildNormalMeanCamber(double* ns, double* tl, int useextr,
                                           const MeanCamberCurveParameters2016* mclParams)
{
  if(m_mc)
  {
    delete m_mc;
    m_mc = nullptr;
    delete m_mclParams;
    m_mclParams = nullptr;
  }

  // bugout(0, L"ler %f ter %f leoff %f teoff %f leo %f teo %f",m_ler,m_ter,leoff,teoff,leo?*leo:0.0,teo?*teo:0.0);
  // bugout(0, L"m_lec m_t0 %f m_t1 %f m_period %f", m_lec->T0(), m_lec->T1(), m_period);

  int numNomMCLPoints = 75;// myGetProfileInt(L"NominalMCLPoints", 75);

  if(ns || tl)
  {
    // m_lec and m_tec sometimes aren't big enough to create a good MCL

    CCurve* tcvc = m_cvc;
    CCurve* tccc = m_ccc;
    CCurve* tlec = m_lec;
    CCurve* ttec = m_tec;

    AdjustSubCurves(&tcvc, &tccc, &tlec, &ttec);

    // bugout(0, _T("before non radius mcl"));
    if(mclParams && ((mclParams->noseBackoff && mclParams->noseBackoff->backoffDistance > 0) ||
                     (mclParams->tailBackoff && mclParams->tailBackoff->backoffDistance > 0)))
    {
      MeanCamberCurveParameters2016 newParams = *mclParams;
      newParams.wholeCurve = m_whole;
      newParams.noseStart = m_lec->T0();
      newParams.noseEnd = m_lec->T1();
      newParams.tailStart = m_tec->T0();
      newParams.tailEnd = m_tec->T1();
      try
      {
        m_mc = Hexagon::Blade::createMeanCamberCurve2016(newParams, m_english);
        m_mclParams = new MeanCamberCurveParameters2016(newParams);
        return 1;
      }
      catch(...)
      {
        std::wstring errorSection = L"nominal section " + std::wstring(m_name);
        ErrorStruct es(BE_MEANCAMBERFAILED, errorSection.c_str());
        m_error->AddError(&es);
        return 0;
      }
    }
    else
    {
      m_mc = new CNurbCurve(numNomMCLPoints, tccc, tcvc, tlec, ttec, useextr);
      delete tcvc;
      delete tccc;
      delete tlec;
      delete ttec;
    }
  }
  else
  {
    // bugout(0, _T("before radius mcl"));
    if(mclParams && ((mclParams->noseBackoff && mclParams->noseBackoff->backoffDistance > 0) ||
                     (mclParams->tailBackoff && mclParams->tailBackoff->backoffDistance > 0)))
    {
      MeanCamberCurveParameters2016 newParams = *mclParams;
      newParams.wholeCurve = m_whole;
      newParams.noseStart = m_lec->T0();
      newParams.noseEnd = m_lec->T1();
      newParams.tailStart = m_tec->T0();
      newParams.tailEnd = m_tec->T1();
      try
      {
        m_mc = Hexagon::Blade::createMeanCamberCurve2016(newParams, m_english);
        m_mclParams = new MeanCamberCurveParameters2016(newParams);
        return 1;
      }
      catch(...)
      {
        std::wstring errorSection = L"nominal section " + std::wstring(m_name);
        ErrorStruct es(BE_MEANCAMBERFAILED, errorSection.c_str());
        m_error->AddError(&es);
        return 0;
      }
    }
    else
    {
      m_mc = new CNurbCurve(numNomMCLPoints, m_tnose, m_ttail, m_ler, m_ter, m_whole, 0);
    }
  }
  m_mc->IsEnglish(m_english);

  if(!m_mc->Valid())
  {
    delete m_mc;
    m_mc = 0;
    delete m_mclParams;
    m_mclParams = nullptr;
    return 0;
  }

  return 1;
}

int CNominalSection::BuildNormalSubCurves(double* ns, double* tl, double leoff, double teoff, int forceoffsets,
                                          double* leo, double* teo, int useextr,
                                          const MeanCamberCurveParameters2016* mclParams)
{
  // determine extreme t values for m_nose and m_tail points

  double xxx[2];
  m_whole->CalcPoint(xxx, m_whole->T0());
  bugout(0, L"BuildNormalSubCurves (%f ,%f) T0(%f)", xxx[0], xxx[1], m_whole->T0());
  m_whole->CalcPoint(xxx, m_period/2.0);
  bugout(0, L"BuildNormalSubCurves (%f ,%f) P/2(%f)", xxx[0], xxx[1], m_period / 2.0);
  m_whole->CalcPoint(xxx, m_period );
  bugout(0, L"BuildNormalSubCurves (%f ,%f)  leoff=%f teoff=%f, T1(%f)", xxx[0], xxx[1], leoff, teoff, m_whole->T1());

  // first, just find the closest points to the specified nose and tail.
  m_tnose=0.0;
  m_ttail=0.0;
    
  for (int pass = 0; pass < 2; pass++)
  {
    if (pass == 0)
    {
      InitialNormalNoseTailTValues(leoff, teoff);
    }
    else
    {
     //RefineNormalNoseTailTValues(leoff, teoff, forceoffsets, leo, teo);
    }

  /*  if (m_startCC)
    {
      m_t0[CVC] = m_t1[LEC];
      m_t1[CVC] = m_t0[TEC];
      m_t0[CCC] = m_t1[TEC];
      m_t1[CCC] = m_t0[LEC];
    }
    else
    {
      m_t0[CVC] = m_t1[TEC];
      m_t1[CVC] = m_t0[LEC];
      m_t0[CCC] = m_t1[LEC];
      m_t1[CCC] = m_t0[TEC];
    }*/

    //int i;
    //for (i=0; i<4; i++)
    //  if (m_t0[i] > m_t1[i])
    //    m_t1[i] += m_period;
    m_t0[LEC] = m_whole->T0();
    m_t1[LEC] = m_whole->T1();
    if (m_lec)
      delete m_lec;
   // 


  
   m_t0[LEC] = m_whole->T0();
   m_t1[LEC] = m_whole->T1();
   m_lec = new CSubCurve(m_whole, m_t0[LEC], m_t1[LEC], m_period);
    m_lec->Extreme(m_whole->T0());
   m_t0[TEC] = m_period / 2;
   m_t1[TEC] = m_whole->T1();
    if(m_tec)
      delete m_tec;
    m_tec = new CSubCurve(m_whole, m_t0[TEC], m_t1[TEC], m_period/2);
    m_tec->Extreme(m_whole->T1());
    /* 
    if (m_cvc)
      delete m_cvc;
    m_cvc = new CSubCurve(m_whole, m_t0[CVC], m_t1[CVC], m_period);

    if (m_ccc)
      delete m_ccc;
    m_ccc = new CSubCurve(m_whole, m_t0[CCC], m_t1[CCC], m_period);
    
    */

    // Calc Pitch Circle

    m_pitch[0] = m_pitch[1] = m_pitch[2] = -1.0;
    bugout(0, L"BuildNormalSubCurves: create SubCurve LEC [%f %f] CVC [%f %f]  TEC [%f %f] CCC [%f %f] P", m_t0[LEC], m_t1[LEC],
           m_t0[CVC], m_t1[CVC], m_t0[TEC], m_t1[TEC], m_t0[CCC], m_t1[CCC]);
    //auto thicknessResult = Hexagon::Blade::maxThickness(*m_whole, *m_cvc, *m_ccc, 2048, 1e-5);
    //m_pitch[2] = thicknessResult.thickness;

    //int rv;
    //rv = BuildNormalMeanCamber(ns, tl, useextr, mclParams);

    //thicknessResult = Hexagon::Blade::maxThickness(*m_whole, *m_cvc, *m_ccc, 2048, 1e-5);
    //m_pitch[0] = thicknessResult.thicknessCenter[0];
    //m_pitch[1] = thicknessResult.thicknessCenter[1];
    //m_pitch[2] = thicknessResult.thickness;

    //double pp[2];
    //m_cvc->CalcPoint(pp, m_cvc->T0());
    //m_cvc->CalcPoint(pp, m_cvc->T1());
    //m_ccc->CalcPoint(pp, m_ccc->T0());
    //m_ccc->CalcPoint(pp, m_ccc->T1());
  }

 // RefineNormalNoseTailTValues(leoff, teoff, forceoffsets, leo, teo);

  return 1;
}

// normal end sections

int CNominalSection::Calculate(double** m, double* ns, double* tl, double* leo, double* teo, int useextr, bool* vecWarn,
                               const MeanCamberCurveParameters2016* mclParams)
{
  bugout(0, L"Calculating section %s, leoff=%f teoff=%f", m_name, *leo, *teo);
  int forceoffsets = 0;
  double leoff = 0.0, teoff = 0.0;
  if (leo && teo)
  {
    leoff = *leo;
    teoff = *teo;
    forceoffsets = 1;
  }

  // If last point same as first, don't use it.

  if (dist(m[0], m[m_npts-1]) < 1.0e-6)
    m_npts--;

  // Determine indices of furthest points, m_nose and m_tail candidates

  m_ler = 0.0;
  m_ter = 0.0;
  DetermineNormalNoseTail(m, ns, tl, &leoff, &teoff);
  // Check to see if points are defined CW or CCW, determine good starting point (away from the ends)
  
  NormalReorderPoints(m, vecWarn);


  // Determine if point 0 is on CC or CV side.

  double chord[2], v1[2], v2[2];
  chord[0] = m[m_itail][0] - m[m_inose][0];
  chord[1] = m[m_itail][1] - m[m_inose][1];
  normalize(chord, chord);
  v1[0] = m_xv[0];
  v1[1] = m_yv[0];
  v2[0] = m_xv[m_npts/2];
  v2[1] = m_yv[m_npts/2];

  // Note: this may be overridden after spline is computed
  m_startCC = false;
  if (ptlinedist(m[0], chord, v1) < ptlinedist(m[m_npts/2], chord, v2))
    m_startCC = true;

  // Build whole spline curve
  int nomfixdat = FALSE; //myGetProfileInt(L"NominalRemove", 0) == 0 ? FALSE : TRUE;
  int nomtension = FALSE;// myGetProfileInt(L"NominalTension", 0) == 0 ? FALSE : TRUE;
  bugout(0, L"Calculate: new CNurbCurve m_npts(%d)", m_npts);
 // nomfixdat = true;
  /*for(int k = 0; k < m_npts; k++)
  {
    bugout(0, L"Calculating %lf %lf", m_xv[k], m_yv[k]);
  }*/
  m_whole = new CNurbCurve(m_npts, m_xv, m_yv, m_kv, m_english, 0, nomtension, nomfixdat);

  // set m_period of closed curve
  m_period = m_whole->T1() - m_whole->T0();
  bugout(0, L"Calculating curve T0=%lf ,T1=%lf, leoff(%lf) ,teoff(%lf)", m_whole->T0(), m_whole->T1(), leoff, teoff);
  // determine where cvc, ccc, lec and tec transition.  Compute mean camber line
  int rv = BuildNormalSubCurves(ns, tl, leoff, teoff, forceoffsets, leo, teo, useextr, mclParams);

  bugout(0, L"Calculate Done with section %s", m_name);
  return rv;
}

// sections with at least one square/partial end.

int CNominalSection::DetermineSquarePartialEndPoints(double **m, int ntype, int ttype, double *tooClose)
{
  m_errorCode = NS_OK;

  m_letype = (short)ntype;
  m_tetype = (short)ttype;
  //bugout(10, L"Calculating section %s", m_name);
  //bugout(0, L"%s m_letype %d m_tetype %d", m_name, m_letype, m_tetype);
  //bugout(0, L"before m_ni1 %d m_ni2 %d m_ti1 %d m_ti2 %d", m_ni1, m_ni2, m_ti1, m_ti2);

  m_leSquare = 1;  // square
  m_teSquare = 1;

  m_nPhantomLE = 0;
  m_nPhantomTE = 0;

  // do indices need to be swapped?

  m_finishLE = false;

  if (m_letype != EDGE_NORMAL)
  {
    int n = abs(m_ni1 - m_ni2);
    if (n == 0)
    {
      m_errorCode = NS_BADNOSEINDICES;

      ErrorStruct es(BE_BADNOSEINDICES, m_name);
      m_error->AddError(&es);
      return 0;
    }

    int swap = 0;

    if (n < m_npts/2 && m_ni1 > m_ni2)   // a smaller piece, this should cross 0
      swap = 1;
    else if (n > m_npts/2 && m_ni1 < m_ni2)  // larger piece should cross 0
      swap = 1;

    if (swap)
    {
      n = m_ni1;
      m_ni1 = m_ni2;
      m_ni2 = n;
    }

    m_origLE[0][0] = m[m_ni1][0];
    m_origLE[0][1] = m[m_ni1][1];
    m_origLE[1][0] = m[m_ni2][0];
    m_origLE[1][1] = m[m_ni2][1];

    if (m_ni2 - m_ni1 == 1 || (m_ni1 == m_npts-1 && m_ni2 == 0))
      m_finishLE = true;
  }
  //bugout(0, _T("10 %lf %lf ON0"), m_origLE[0][0], m_origLE[0][1]);
  //bugout(0, _T("10 %lf %lf ON1"), m_origLE[1][0], m_origLE[1][1]);

  m_finishTE = false;
  if (m_tetype != EDGE_NORMAL)
  {
    int n = abs(m_ti1 - m_ti2);
    if (n == 0)
    {
      m_errorCode = NS_BADTAILINDICES;


      ErrorStruct es(BE_BADTAILINDICES, m_name);
      m_error->AddError(&es);

      return 0;
    }

    int swap = 0;

    if (n < m_npts/2 && m_ti1 > m_ti2)   // a smaller piece, this should cross 0
      swap = 1;
    else if (n > m_npts/2 && m_ti1 < m_ti2)  // larger piece should cross 0
      swap = 1;

    if (swap)
    {
      n = m_ti1;
      m_ti1 = m_ti2;
      m_ti2 = n;
    }

    m_origTE[0][0] = m[m_ti1][0];
    m_origTE[0][1] = m[m_ti1][1];
    m_origTE[1][0] = m[m_ti2][0];
    m_origTE[1][1] = m[m_ti2][1];

    if (m_ti2 - m_ti1 == 1 || (m_ti1 == m_npts-1 && m_ti2 == 0))
      m_finishTE = true;
  }
  //bugout(0, _T("10 %lf %lf OT0"), m_origTE[0][0], m_origTE[0][1]);
  //bugout(0, _T("10 %lf %lf OT1"), m_origTE[1][0], m_origTE[1][1]);

  //bugout(3, _T("after m_ni1 %d m_ni2 %d m_ti1 %d m_ti2 %d"), m_ni1, m_ni2, m_ti1, m_ti2);

  if (m_letype == EDGE_PARTIAL)
  {
    if (m_ni1 == m_npts-1 && m_ni2 == 0)
      m_leSquare = 2;
    else if ((m_ni1+1)%m_npts == m_ni2)
      m_leSquare = 2;
  }

  // but when points are invented, we don't want them removed, so a value will then be applied.

  //if (m_leSquare == 2)  not true if points are not adjacent
  //  m_finishLE = true;

  if (m_tetype == EDGE_PARTIAL)
  {
    if (m_ti1 == m_npts-1 && m_ti2 == 0)
      m_teSquare = 2;
    else if ((m_ti1+1)%m_npts == m_ti2)
      m_teSquare = 2;
  }

  //if (m_teSquare == 2)  not true if they end points are not adjacent
  //  m_finishTE = true;
  //bugout(3, _T("m_leSquare %d m_teSquare %d m_finishLE %d m_finishTE %d"), m_leSquare, m_teSquare, m_finishLE, m_finishTE);

  if (m_finishLE)
  {
    m_nPhantomLE = 19;

    double x[4], y[4];
    x[0] = m[(m_ni1+m_npts-1)%m_npts][0];
    y[0] = m[(m_ni1+m_npts-1)%m_npts][1];
    m_endKLE[0] = m[(m_ni1+m_npts-1)%m_npts][5];
    x[1] = m[m_ni1][0];
    y[1] = m[m_ni1][1];
    m_endKLE[1] = m[m_ni1][5];
    x[2] = m[m_ni2][0];
    y[2] = m[m_ni2][1];
    m_endKLE[2] = m[m_ni2][5];
    x[3] = m[(m_ni2+1)%m_npts][0];
    y[3] = m[(m_ni2+1)%m_npts][1];
    m_endKLE[3] = m[(m_ni2+1)%m_npts][5];

    if (m_letype == EDGE_PARTIAL)
    {
      //adjust points 0 and 3 slightly to help with the normal vectors for those points
      double slide[2];
      slide[0] = 0.05*(x[1] - x[0]);
      slide[1] = 0.05*(y[1] - y[0]);
      x[1] += slide[0];
      y[1] += slide[1];
      slide[0] = 0.05*(x[2] - x[3]);
      slide[1] = 0.05*(y[2] - y[3]);
      x[2] += slide[0];
      y[2] += slide[1];
    }
    else
    {
      //zzz may need to mess with square ends that have two corners and one end point
    }

    double newstart[2], newend[2];
    InventEnd(1,2,4,x,y,19,m_pxLE,m_pyLE,-1.0,newstart,newend);
    double minD = 1.0;
    for (int p = 1; p<19; p++)
    {
      double d = _hypot(m_pxLE[p] - m_pxLE[p-1], m_pyLE[p] - m_pyLE[p-1]);
      if (p == 1 || d < minD)
      {
        minD = d;
      }
    }
    *tooClose = 0.9*minD;

    m[m_ni1][0] = newstart[0];
    m[m_ni1][1] = newstart[1];
    m[m_ni2][0] = newend[0];
    m[m_ni2][1] = newend[1];
  }

  //bugout(3, _T("m_ti1 %d m_ti2 %d npts %d"), m_ti1, m_ti2, m_npts);
  if (m_finishTE)
  {
    m_nPhantomTE = 19;

    double x[4], y[4];
    x[0] = m[(m_ti1+m_npts-1)%m_npts][0];
    y[0] = m[(m_ti1+m_npts-1)%m_npts][1];
    m_endKTE[0] = m[(m_ti1+m_npts-1)%m_npts][5];
    x[1] = m[m_ti1][0];
    y[1] = m[m_ti1][1];
    m_endKTE[1] = m[m_ti1][5];
    x[2] = m[m_ti2][0];
    y[2] = m[m_ti2][1];
    m_endKTE[2] = m[m_ti2][5];
    x[3] = m[(m_ti2+1)%m_npts][0];
    y[3] = m[(m_ti2+1)%m_npts][1];
    m_endKTE[3] = m[(m_ti2+1)%m_npts][5];

    if (m_tetype == EDGE_PARTIAL)
    {
      //adjust points 0 and 3 slightly to help with the normal vectors for those points
      double slide[2];
      slide[0] = 0.05*(x[1] - x[0]);
      slide[1] = 0.05*(y[1] - y[0]);
      x[1] += slide[0];
      y[1] += slide[1];
      slide[0] = 0.05*(x[2] - x[3]);
      slide[1] = 0.05*(y[2] - y[3]);
      x[2] += slide[0];
      y[2] += slide[1];
    }
    else
    {
      //zzz may need to mess with square ends that have two corners and one end point
    }

    //bugout(3, _T("before InventEnd"));
    double newstart[2], newend[2];
    InventEnd(1,2,4,x,y,19,m_pxTE,m_pyTE,-1.0,newstart,newend);
    double minD = 1.0;
    for (int p = 1; p<19; p++)
    {
      double d = _hypot(m_pxTE[p] - m_pxTE[p-1], m_pyTE[p] - m_pyTE[p-1]);
      if (p == 1 || d < minD)
      {
        minD = d;
      }
    }
    *tooClose = 0.9*minD;

    m[m_ti1][0] = newstart[0];
    m[m_ti1][1] = newstart[1];
    m[m_ti2][0] = newend[0];
    m[m_ti2][1] = newend[1];
    //bugout(3, _T("after InventEnd"));
  }
  //bugout(3, _T("nPhantom %d m_leSquare %d m_teSquare %d"), nPhantom, m_leSquare, m_teSquare);
  
  return 1;
}

int CNominalSection::FillArrayComputeSplineSquarePartial(double **m, double *leoff, double *teoff, double tooClose)
{
  int i, j;
  // Step 1. if last point same as first, don't use it.

  if (dist(m[0], m[m_npts-1]) < 1.0e-6)
    m_npts--;

  if (m_ni1 < 0 || m_ni2 < 0)
    m_leSquare = 0;

  if (m_ti1 < 0 || m_ti2 < 0)
    m_teSquare = 0;

  //bugout(3, _T("m_leSquare %d m_teSquare %d"), m_leSquare, m_teSquare);

  // Step 2. Determine indices of furthest points, m_nose and m_tail candidates


  if (m_leSquare)
  {
    // m_inose should be okay for partial sections, decides cw/ccw & cv/cc
    
    if (m_ni1 > m_ni2) // compute mid index when points are wrapping around.
      m_inose = ((m_ni1 + m_ni2 + m_npts)/2)%m_npts;
    else
      m_inose = (m_ni1 + m_ni2)/2;  // should be good enough

    if (m_finishLE) // use phantom point
    {
      m_nose[0] = m_pxLE[9];
      m_nose[1] = m_pyLE[9];
    }
    else
    {
      m_nose[0] = m[m_inose][0];
      m_nose[1] = m[m_inose][1];
    }
    //bugout(0, _T("10 %lf %lf Nose m_ni1 %d m_ni2 %d m_inose %d"), m_nose[0], m_nose[1], m_ni1, m_ni2, m_inose);
  }

  if (m_teSquare)
  {
    // m_itail should be okay for partial sections, decides cw/ccw & cv/cc

    if (m_ti1 > m_ti2)  // compute mid index when points are wrapping around.
      m_itail = ((m_ti1 + m_ti2 + m_npts)/2)%m_npts;
    else
      m_itail = (m_ti1 + m_ti2)/2;

    if (m_finishTE) // use phantom point
    {
      m_tail[0] = m_pxTE[9];
      m_tail[1] = m_pyTE[9];
    }
    else
    {
      m_tail[0] = m[m_itail][0];
      m_tail[1] = m[m_itail][1];
    }
    //bugout(0, _T("10 %lf %lf Tail m_ni1 %d m_ni2 %d m_itail %d"), m_tail[0], m_tail[1], m_ti1, m_ti2, m_itail);
  }

  // this routine is called only when at least one end is square or partial
  if (!m_leSquare)
  {
    double bestd = 1.0e20;
    int lim1 = (m_itail + 3*m_npts/8)%m_npts;
    for (i=0; i<m_npts/4; i++)
    {
      int jj = (lim1 + i)%m_npts;
      double rad = radius_3points(m[jj], m[(jj+1)%m_npts], m[(jj+2)%m_npts]);
      if (rad < bestd)
      {
        bestd = rad;
        m_inose = jj;
      }
    }
    if (*leoff < 2.0*bestd)
      *leoff = 2.0*bestd;

    if (m_ni1 >= 0)  // above was to make sure leoff is reasonable,  set inose back if it was set.
    {
      m_inose = m_ni1;
    }

    m_nose[0] = m[m_inose][0];
    m_nose[1] = m[m_inose][1];
  }

  if (!m_teSquare)
  {
    double bestd = 1.0e20;
    int lim1 = (m_inose + 3*m_npts/8)%m_npts;
    for (i=0; i<m_npts/4; i++)
    {
      int jj = (lim1 + i)%m_npts;
      double rad = radius_3points(m[jj], m[(jj+1)%m_npts], m[(jj+2)%m_npts]);
      if (rad < bestd)
      {
        bestd = rad;
        m_itail = jj;
      }
    }
    if (*teoff < 2.0*bestd)
      *teoff = 2.0*bestd;

    if (m_ti1 >= 0)  // above was to make sure that teoff is okay, if set, set itail back
    {
      m_itail = m_ti1;
    }

    m_tail[0] = m[m_itail][0];
    m_tail[1] = m[m_itail][1];
  }

  //bugout(0, _T("10 %f %f m_nose"), m_nose[0], m_nose[1]);
  //bugout(0, _T("10 %f %f m_tail"), m_tail[0], m_tail[1]);
  
  // Step 4. Check to see if points are defined CW or CCW

  double v1[2], v2[2];
  j = (m_inose + m_npts/8)%m_npts;
  v1[0] = m[j][0] - m[m_inose][0];
  v1[1] = m[j][1] - m[m_inose][1];
  j = (m_inose + m_npts - m_npts/8)%m_npts;
  v2[0] = m[j][0] - m[m_inose][0];
  v2[1] = m[j][1] - m[m_inose][1];

  int cw = 0;
  if (curl(v1, v2) < 0.0) // points are clock wise
    cw = 1;

  // Step 5. Fill array, ordered CCW, p3oint 0 in middle of side before LE

  m_xv = new double[m_npts + m_nPhantomLE + m_nPhantomTE];
  m_yv = new double[m_npts + m_nPhantomLE + m_nPhantomTE];
  m_ox = new double[m_npts + m_nPhantomLE + m_nPhantomTE];
  m_oy = new double[m_npts + m_nPhantomLE + m_nPhantomTE];
  m_oz = new double[m_npts + m_nPhantomLE + m_nPhantomTE];
  if (m_havek)
  {
    m_kv = new double[m_npts + m_nPhantomLE + m_nPhantomTE];
    m_oi = new double[m_npts + m_nPhantomLE + m_nPhantomTE];
    m_oj = new double[m_npts + m_nPhantomLE + m_nPhantomTE];
    m_ok = new double[m_npts + m_nPhantomLE + m_nPhantomTE];
  }
  if (m_haveTol)
  {
    m_omtol = new double[m_npts + m_nPhantomLE + m_nPhantomTE];
    m_optol = new double[m_npts + m_nPhantomLE + m_nPhantomTE];
  }

  // try to start curve halfway between m_nose and m_tail
  // want the m_nose point to be around the quarter point
  // old logic did not work for partial sections.

  double halfdist = 0.5*dist(m_nose, m_tail);
  double closest = 1.0e20;

  int izero=0;
  for (i=0; i<m_npts/2; i++)
  {
    int jj;
    double xy[2];
    if (cw)  // look for point closest to halfdist in [m_inose, m_inose + m_npts/2]
      jj = (m_inose + i)%m_npts;
    else     // look for point closest to halfdist in [m_inose - m_npts/2, m_inose]
      jj = (m_npts + m_inose - i)%m_npts;

    xy[0] = m[jj][0];
    xy[1] = m[jj][1];
    double d = fabs(dist(m_nose, xy) - halfdist);
    if (d < closest)
    {
      closest = d;
      izero = jj;
    }
  }
  //bugout(3, _T("izero %d cw %d closest %f"), izero, cw, closest);

  int ani1 = -1, ani2 = -1, ati1 = -1, ati2 = -1;

  for (i=0; i<m_npts; i++)
  {
    if (cw)
      j = (m_npts + izero - i) % m_npts;
    else
      j = (izero + i) % m_npts;

    m_xv[i] = m[j][0];
    m_yv[i] = m[j][1];
    if (m_havek)
      m_kv[i] = m[j][5];

    if (m_ni1 >= 0)
      if (_hypot(m[m_ni1][0] - m_xv[i], m[m_ni1][1] - m_yv[i]) < 1.0e-6)
        ani1 = i;

    if (m_ni2 >= 0 && _hypot(m[m_ni2][0] - m_xv[i], m[m_ni2][1] - m_yv[i]) < 1.0e-6)
      ani2 = i;

    if (m_ti1 >= 0)
      if (_hypot(m[m_ti1][0] - m_xv[i], m[m_ti1][1] - m_yv[i]) < 1.0e-6)
        ati1 = i;

    if (m_ti2 >= 0 && _hypot(m[m_ti2][0] - m_xv[i], m[m_ti2][1] - m_yv[i]) < 1.0e-6)
      ati2 = i;

    m_ox[i] = m[i][0];
    m_oy[i] = m[i][1];
    m_oz[i] = m[i][2];
    if (m_havek)
    {
      m_oi[i] = m[i][3];
      m_oj[i] = m[i][4];
      m_ok[i] = m[i][5];
    }
    if (m_haveTol)
    {
      m_omtol[i] = m[i][6];
      m_optol[i] = m[i][7];
    }
  }

  // Step 6. Determine if point 0 is on CC or CV side.

  double chord[2];
  chord[0] = m[m_itail][0] - m[m_inose][0];
  chord[1] = m[m_itail][1] - m[m_inose][1];
  normalize(chord, chord);
  v1[0] = m_xv[0];
  v1[1] = m_yv[0];
  v2[0] = m_xv[m_npts/2];
  v2[1] = m_yv[m_npts/2];

  // this may be overridden after spline is computed
  m_startCC = false;
  if (ptlinedist(m[0], chord, v1) < ptlinedist(m[m_npts/2], chord, v2))
    m_startCC = true;

  // Step 7. Build splines

  if (m_nPhantomLE > 0)  // need to add the phantom points
  {
    int insLE = ani1;
    if  (ani2 > ani1)
      insLE = ani2;

    //bugout(3, _T("inserting points at %d"), insLE);

    for (i = m_npts+18; i>insLE; i--)
    {
      int jj = (i - 19 + m_npts)%m_npts;
      m_xv[i] = m_xv[jj];
      m_yv[i] = m_yv[jj];
      if (m_havek)
        m_kv[i] = m_kv[jj];

      if (m_ti1 >= 0)
        if (_hypot(m[m_ti1][0] - m_xv[i], m[m_ti1][1] - m_yv[i]) < 1.0e-6)
          ati1 = i;

      if (m_ti2 >= 0)
        if (_hypot(m[m_ti2][0] - m_xv[i], m[m_ti2][1] - m_yv[i]) < 1.0e-6)
          ati2 = i;
    }

    double chordKLE = 0.0;  // don't remember what the future intent of this is.
    //bugout(3, _T("insLE %d chordK %f endK %f %f %f %f"), insLE, chordK, endK[0], endK[1], endK[2], endK[3]);
    for (i=0; i<19; i++)
    {
      if (cw)
        j = 18 - i;
      else
        j = i;

      m_xv[insLE + j] = m_pxLE[i];
      m_yv[insLE + j] = m_pyLE[i];

      if (m_havek)
      {
        double ramp, start;
        int step;

        if (i < 4)
        {
          ramp = 0.2*(m_endKLE[1] - m_endKLE[0]);
          start = m_endKLE[1] - ramp;
          step = i;
        }
        else if (i == 4)
        {
          m_kv[insLE + j] = m_endKLE[1];
          continue;
        }
        else if (i < 9)
        {
          ramp = chordKLE - m_endKLE[1];
          start = m_endKLE[1] + 0.25*ramp;
          step = i - 5;
        }
        else if (i == 9)
        {
          m_kv[insLE + j] = chordKLE;
          continue;
        }
        else if (i < 14)
        {
          ramp = m_endKLE[2] - chordKLE;
          start = chordKLE;
          step = i - 10;
        }
        else if (i == 14)
        {
          m_kv[insLE + j] = m_endKLE[2];
          continue;
        }
        else
        {
          ramp = 0.2*(m_endKLE[3] - m_endKLE[2]);
          start = m_endKLE[2] + 0.25*ramp;
          step = i;
        }
        m_kv[insLE + j] = start + 0.25*step*ramp;
      }
    }
  }

  int numptsLE = m_npts + m_nPhantomLE;

  if (m_nPhantomTE > 0)  // need to add the phantom points
  {
    int insTE = ati1;
    if  (ati2 > ati1)
      insTE = ati2;

    for (i = numptsLE+18; i>insTE; i--)
    {
      int jj = (i - 19 + numptsLE)%numptsLE;
      m_xv[i] = m_xv[jj];
      m_yv[i] = m_yv[jj];
      if (m_havek)
        m_kv[i] = m_kv[jj];
    }

    double chordKTE = 0.0;  // don't remember what the future intent of this is.
    //bugout(3, _T("insTE %d chordK %f endK %f %f %f %f"), insTE, chordK, endK[0], endK[1], endK[2], endK[3]);
    for (i=0; i<19; i++)
    {
      if (cw)
        j = 18 - i;
      else
        j = i;

      m_xv[insTE + j] = m_pxTE[i];
      m_yv[insTE + j] = m_pyTE[i];

      if (m_havek)
      {
        double ramp, start;
        int step;

        if (i < 4)
        {
          ramp = 0.2*(m_endKTE[1] - m_endKTE[0]);
          start = m_endKTE[1] - ramp;
          step = i;
        }
        else if (i == 4)
        {
          m_kv[insTE + j] = m_endKTE[1];
          continue;
        }
        else if (i < 9)
        {
          ramp = chordKTE - m_endKTE[1];
          start = m_endKTE[1] + 0.25*ramp;
          step = i - 5;
        }
        else if (i == 9)
        {
          m_kv[insTE + j] = chordKTE;
          continue;
        }
        else if (i < 14)
        {
          ramp = m_endKTE[2] - chordKTE;
          start = chordKTE;
          step = i - 10;
        }
        else if (i == 14)
        {
          m_kv[insTE + j] = m_endKTE[2];
          continue;
        }
        else
        {
          ramp = 0.2*(m_endKTE[3] - m_endKTE[2]);
          start = m_endKTE[2] + 0.25*ramp;
          step = i;
        }
        m_kv[insTE + j] = start + 0.25*step*ramp;
      }
    }
  }

  int nomfixdat = FALSE;//myGetProfileInt(L"NominalRemove", 0) == 0 ? FALSE : TRUE;
  int nomtension = FALSE;// myGetProfileInt(L"NominalTension", 0) == 0 ? FALSE : TRUE;

  //for (i=0; i<m_npts+nPhantom; i++) bugout(3, _T("i=%d k=%f"), i, m_kv[i]);
  //for (i=0; i<m_npts + m_nPhantomLE + m_nPhantomTE; i++) bugout(0, _T("10 %f %f %d whole curve point"), m_xv[i], m_yv[i], i);

  tooClose = 0.000001;  // this was set too large for partial sections.  For now, just make this a small number so no points are thrown out.
  
  m_whole = new CNurbCurve(m_npts + m_nPhantomLE + m_nPhantomTE, m_xv, m_yv, m_kv, m_english, 1, nomtension, nomfixdat, tooClose);
  
  //bugout(0, m_name);
  //m_whole->BugOut(600);

  // set m_period of closed curve

  m_period = m_whole->T1() - m_whole->T0();

  // restore endpoints that may have been moved by InventEnd, these arrays are written to the math file

  if (m_finishLE)
  {
    m_ox[m_ni1] = m_origLE[0][0];
    m_oy[m_ni1] = m_origLE[0][1];
    m_ox[m_ni2] = m_origLE[1][0];
    m_oy[m_ni2] = m_origLE[1][1];
  }

  if (m_finishTE)
  {
    m_ox[m_ti1] = m_origTE[0][0];
    m_oy[m_ti1] = m_origTE[0][1];
    m_ox[m_ti2] = m_origTE[1][0];
    m_oy[m_ti2] = m_origTE[1][1];
  }

  return 1;
}

int CNominalSection::DetermineTransitionSquarePartial(double leoff, double teoff)
{
  // determine reasonable break off for LE and CV/CC transition

  double xyz[2];
  if (m_leSquare == 2)  // partial end
  {
    //bugout(0, _T("10 %f %f n0"), m_origLE[0][0], m_origLE[0][1]);
    //bugout(0, _T("10 %f %f n1"), m_origLE[1][0], m_origLE[1][1]);
    //bugout(0, _T("10 %f %f ne"), m_nose[0], m_nose[1]);
    if (m_whole->ClosestPoint(m_origLE[0], xyz, &m_t0[LEC], 0, 0.0, 0.0, 500) > 0.001)
      m_whole->ClosestPoint(m_origLE[0], xyz, &m_t0[LEC], 0, 0.0, 0.0, 1500);

    if (m_whole->ClosestPoint(m_origLE[1], xyz, &m_t1[LEC], 0, 0.0, 0.0, 500) > 0.001)
      m_whole->ClosestPoint(m_origLE[1], xyz, &m_t1[LEC], 0, 0.0, 0.0, 1500);

    if (m_whole->ClosestPoint(m_nose, xyz, &m_tnose, 0, 0.0, 0.0, 500) > 0.001)
      m_whole->ClosestPoint(m_nose, xyz, &m_tnose, 0, 0.0, 0.0, 1500);

    if (m_t0[LEC] > m_t1[LEC])
    {
      double tmp = m_t0[LEC];
      m_t0[LEC] = m_t1[LEC];
      m_t1[LEC] = tmp;
    }
  }
  else if (m_leSquare == 1)  // square end
  {
    double vec[2], tanv[2], nrmv[2], bp[2], bp2[2];

    m_whole->NewClosestPoint(m_nose, bp, &m_tnose, tanv, 0.0, m_period/2.0, 800);
    vec[0] = m_tail[0] - m_nose[0];
    vec[1] = m_tail[1] - m_nose[1];
    nrmv[0] = tanv[1];
    nrmv[1] = -tanv[0];
    normalize(nrmv, nrmv);
    if (dot(nrmv, vec) < 0.0)
    {
      nrmv[0] *= -1.0;
      nrmv[1] *= -1.0;
    }
    bp[0] = m_nose[0] + leoff*nrmv[0];
    bp[1] = m_nose[1] + leoff*nrmv[1];

    m_whole->NewLineIntersect(bp, tanv, bp2, 0.0, m_tnose, &m_t0[LEC]);
    //bugout(0, L"10 %lf %lf L1", bp2[0], bp2[1]);
    m_whole->NewLineIntersect(bp, tanv, bp2, m_tnose, m_period/2, &m_t1[LEC]);
    //bugout(0, L"10 %lf %lf L2", bp2[0], bp2[1]);
  }
  else  // normal end
  {
    double bp[2];

    m_whole->NewClosestPoint(m_nose, bp, &m_tnose, 0, 0.0, m_period/2.0, 800);
    //bugout(0, L"10 %lf %lf NP", bp[0], bp[1]);
    m_whole->CircIntersect(m_nose, leoff, bp, 0.0, m_tnose, &m_t0[LEC]);
    //bugout(0, L"10 %lf %lf L1", bp[0], bp[1]);
    m_whole->CircIntersect(m_nose, leoff, bp, m_tnose, m_period/2, &m_t1[LEC]);
    //bugout(0, L"10 %lf %lf L2", bp[0], bp[1]);
  }

  // determine reasonable break off for TE and CV/CC transition

  if (m_teSquare == 2) // partial end
  {
    //bugout(0, _T("10 %f %f t0"), m_origTE[0][0], m_origTE[0][1]);
    //bugout(0, _T("10 %f %f t1"), m_origTE[1][0], m_origTE[1][1]);
    //bugout(0, _T("10 %f %f te"), m_tail[0], m_tail[1]);
    if (m_whole->ClosestPoint(m_origTE[0], xyz, &m_t0[TEC], 0, 0.0, 0.0, 500) > 0.001)
      m_whole->ClosestPoint(m_origTE[0], xyz, &m_t0[TEC], 0, 0.0, 0.0, 1500);

    if (m_whole->ClosestPoint(m_origTE[1], xyz, &m_t1[TEC], 0, 0.0, 0.0, 500) > 0.001)
      m_whole->ClosestPoint(m_origTE[1], xyz, &m_t1[TEC], 0, 0.0, 0.0, 1500);

    if (m_whole->ClosestPoint(m_tail, xyz, &m_ttail, 0, 0.0, 0.0, 500) > 0.001)
      m_whole->ClosestPoint(m_tail, xyz, &m_ttail, 0, 0.0, 0.0, 1500);

    if (m_t0[TEC] > m_t1[TEC])
    {
      double tmp = m_t0[TEC];
      m_t0[TEC] = m_t1[TEC];
      m_t1[TEC] = tmp;
    }
  }
  else if (m_teSquare == 1) // square end
  {
    double vec[2], tanv[2], nrmv[2], bp[2], bp2[2];

    //bugout(0, L"10 %lf %lf m_tail", m_tail[0], m_tail[1]);
    m_whole->NewClosestPoint(m_tail, bp, &m_ttail, tanv, m_period/2.0, m_period, 800);
    //bugout(0, L"10 %lf %lf TP", bp[0], bp[1]);

    vec[0] = m_nose[0] - m_tail[0];
    vec[1] = m_nose[1] - m_tail[1];
    nrmv[0] = tanv[1];
    nrmv[1] = -tanv[0];
    normalize(nrmv, nrmv);
    if (dot(nrmv, vec) < 0.0)
    {
      nrmv[0] *= -1.0;
      nrmv[1] *= -1.0;
    }
    bp[0] = m_tail[0] + teoff*nrmv[0];
    bp[1] = m_tail[1] + teoff*nrmv[1];

    m_whole->NewLineIntersect(bp, tanv, bp2, m_period/2, m_ttail, &m_t0[TEC]);
    //bugout(0, L"10 %lf %lf T1", bp2[0], bp2[1]);
    m_whole->NewLineIntersect(bp, tanv, bp2, m_ttail, m_period, &m_t1[TEC]);
    //bugout(0, L"10 %lf %lf T2", bp2[0], bp2[1]);
  }
  else  // normal end
  {
    double bp[2];

    m_whole->NewClosestPoint(m_tail, bp, &m_ttail, 0, m_period/2.0, m_period, 800);
    //bugout(0, L"10 %lf %lf TP", bp[0], bp[1]);
    m_whole->CircIntersect(m_tail, teoff, bp, m_period/2, m_ttail, &m_t0[TEC]);
    //bugout(0, L"10 %lf %lf T1", bp[0], bp[1]);
    m_whole->CircIntersect(m_tail, teoff, bp, m_ttail, m_period, &m_t1[TEC]);
    //bugout(0, L"10 %lf %lf T2", bp[0], bp[1]);
  }

  if (m_convexSide > 0 && m_convexSide < 5)
  {
    double pt1[2], vec[2], pt2[2];

    pt1[0] = m_xv[0];
    pt1[1] = m_yv[0];

    vec[0] = vec[1] = 0;
    if (m_convexSide == 1)
      vec[0] = 1.0;
    else if (m_convexSide == 2)
      vec[0] = -1.0;
    else if (m_convexSide == 3)
      vec[1] = 1.0;
    else if (m_convexSide == 4)
      vec[1] = -1.0;

    if (m_whole->LineIntersect(pt1, vec, pt2, m_t1[LEC], m_t0[TEC]))
    {
      if (projdist(pt1, vec, pt2) < 0.0)  // the first point is on the Convex side
        m_startCC = false;
      else
        m_startCC = true;
    }
  }

  if (m_startCC)
  {
    m_t0[CVC] = m_t1[LEC];
    m_t1[CVC] = m_t0[TEC];
    m_t0[CCC] = m_t1[TEC];
    m_t1[CCC] = m_t0[LEC];
  }
  else
  {
    m_t0[CVC] = m_t1[TEC];
    m_t1[CVC] = m_t0[LEC];
    m_t0[CCC] = m_t1[LEC];
    m_t1[CCC] = m_t0[TEC];
  }

  int i;
  for (i=0; i<4; i++)
    if (m_t0[i] > m_t1[i])
      m_t1[i] += m_period;
  
  //double zzz[2];
  //m_whole->CalcPoint(zzz, m_t0[LEC]);    
  //bugout(0, _T("10 %lf %lf L0"), zzz[0], zzz[1]);
  //m_whole->CalcPoint(zzz, m_t1[LEC]);    
  //bugout(0, _T("10 %lf %lf L1"), zzz[0], zzz[1]);
  //m_whole->CalcPoint(zzz, m_t0[TEC]);    
  //bugout(0, _T("10 %lf %lf T0"), zzz[0], zzz[1]);
  //m_whole->CalcPoint(zzz, m_t0[TEC]);    
  //bugout(0, _T("10 %lf %lf T1"), zzz[0], zzz[1]);

  m_lec = new CSubCurve(m_whole, m_t0[LEC], m_t1[LEC], m_period);
  m_lec->Extreme(m_tnose);

  m_tec = new CSubCurve(m_whole, m_t0[TEC], m_t1[TEC], m_period);
  m_tec->Extreme(m_ttail);

  m_cvc = new CSubCurve(m_whole, m_t0[CVC], m_t1[CVC], m_period);

  m_ccc = new CSubCurve(m_whole, m_t0[CCC], m_t1[CCC], m_period);

  return 1;
}

int CNominalSection::BuildSquarePartialMeanCamber(int useextr, const MeanCamberCurveParameters2016* mclParams)
{
  if(m_mc)
  {
    delete m_mc;
    m_mc = nullptr;
    delete m_mclParams;
    m_mclParams = nullptr;
  }

  if(mclParams && ((mclParams->noseBackoff && mclParams->noseBackoff->backoffDistance > 0) ||
                   (mclParams->tailBackoff && mclParams->tailBackoff->backoffDistance > 0)))
  {
    MeanCamberCurveParameters2016 newParams = *mclParams;
    newParams.wholeCurve = m_whole;
    newParams.noseStart = m_lec->T0();
    newParams.noseEnd = m_lec->T1();
    newParams.tailStart = m_tec->T0();
    newParams.tailEnd = m_tec->T1();
    try
    {
      m_mc = Hexagon::Blade::createMeanCamberCurve2016(newParams, m_english);
      m_mclParams = new MeanCamberCurveParameters2016(newParams);
      return 1;
    }
    catch(...)
    {
      std::wstring errorSection = L"nominal section " + std::wstring(m_name);
      ErrorStruct es(BE_MEANCAMBERFAILED, errorSection.c_str());
      m_error->AddError(&es);
      return 0;
    }
  }

  // build mean camber

  CCurve *uselec = m_lec;
  CCurve *usetec = m_tec;

  //bugout(0, L"LEC");
  //m_lec->BugOut(200);
  //bugout(0, L"TEC");
  //m_tec->BugOut(200);

  if (m_leSquare == 2 || (m_noseforce && m_leSquare == 1) || m_ni1 < 0)
  {
    double xyz[3];
    xyz[0] = m_nose[0];
    xyz[1] = m_nose[1];
    xyz[2] = 0.0;
    if (m_leSquare == 2 || m_ni1 < 0)
      uselec = new CCircle(xyz, 0, 0, -1.0);
    else
    {
      double cp[2];
      cp[0] = xyz[0] + 1.0;
      cp[1] = xyz[1];
      uselec = new CCircle(xyz, cp, cp, 1.0);
    }
    useextr &= ~1;
  }

  if (m_teSquare == 2 || (m_tailforce && m_teSquare == 1) || m_ti1 < 0)
  {
    double xyz[3];
    xyz[0] = m_tail[0];
    xyz[1] = m_tail[1];
    xyz[2] = 0.0;
    if (m_teSquare == 2 || m_ti1 < 0)
      usetec = new CCircle(xyz, 0, 0, -1.0);
    else
    {
      double cp[2];
      cp[0] = xyz[0] + 1.0;
      cp[1] = xyz[1];
      usetec = new CCircle(xyz, cp, cp, 1.0);
    }

    useextr &= ~2;
  }

  int numNomMCLPoints = 75;// myGetProfileInt(L"NominalMCLPoints", 75);
  
  //bugout(0, L"CVC");
  //m_cvc->BugOut(200);
  //bugout(0, L"CCC");
  //m_ccc->BugOut(200);
  //bugout(0, L"LEC");
  //uselec->BugOut(200);
  //bugout(0, L"TEC");
  //usetec->BugOut(200);

  m_mc = new CNurbCurve(numNomMCLPoints, m_cvc, m_ccc, uselec, usetec, useextr);
  m_mc->IsEnglish(m_english);

  if (uselec && uselec != m_lec)
    delete uselec;
  if (usetec && usetec != m_tec)
    delete usetec;
  if (!m_mc->Valid())
  {
    delete m_mc;
    m_mc = 0;
    delete m_mclParams;
    m_mclParams = nullptr;
    return 0;
  }

  return 1;
}
// For sections with at least one square and/or partial ends

int CNominalSection::Calculate(double** m, int ni1, int ni2, int ti1, int ti2, int ntype, int ttype, double* leo,
                               double* teo, int useextr, const MeanCamberCurveParameters2016* mclParams)
{
  m_ni1 = ni1;
  m_ni2 = ni2;
  m_ti1 = ti1;
  m_ti2 = ti2;
  
  /*
  if (0 <= m_ni1 && 0 <= m_ni2)
  {
    bugout(0, _T("10 %lf %lf N1Orig"), m[m_ni1][0], m[m_ni1][1]);
    bugout(0, _T("10 %lf %lf N2Orig"), m[m_ni2][0], m[m_ni2][1]);
  }

  if (0 <= m_ni1 && 0 <= m_ni2)
  {
    bugout(0, _T("10 %lf %lf T1Orig"), m[m_ti1][0], m[m_ti1][1]);
    bugout(0, _T("10 %lf %lf T2Orig"), m[m_ti2][0], m[m_ti2][1]);
  }
  */
  
  double tooClose = 0.0;  // this is used when spline is constructed, zero will apply something normal,

  if (DetermineSquarePartialEndPoints(m, ntype, ttype, &tooClose) == 0)
    return 0;

  int forceoffsets = 0;
  double leoff = 0.0, teoff = 0.0;
  if (leo && teo)
  {
    leoff = *leo;
    teoff = *teo;
    forceoffsets = 1;   // only for normal end if there is one...
  }

  if (FillArrayComputeSplineSquarePartial(m, &leoff, &teoff, tooClose) == 0)
    return 0;

  if (forceoffsets)
  {
    leoff = *leo; // set back to what they were
    teoff = *teo;
  }

  //bugout(3, _T("leoff %f teoff %f"), leoff, teoff);
  if (DetermineTransitionSquarePartial(leoff, teoff) == 0)
    return 0;
    
  // determine extreme t values for m_nose and m_tail points

  m_tnose = 0.0;
  m_ttail = 0.0;

  // Calc Pitch Circle

  int rv;
  rv = BuildSquarePartialMeanCamber(useextr, mclParams);

  if(m_mc == NULL || !m_mc->Valid())
    return 0;

  auto thicknessResult = Hexagon::Blade::maxThickness(*m_whole, *m_cvc, *m_ccc, 2048, 1e-5);
  m_pitch[0] = thicknessResult.thicknessCenter[0];
  m_pitch[1] = thicknessResult.thicknessCenter[1];
  m_pitch[2] = thicknessResult.thickness;

  // m_cvc->CalcPoint(pxyz, m_cvc->T0());
  // bugout(0, _T("10 %f %f S1"), pxyz[0], pxyz[1]);
  // m_cvc->CalcPoint(pxyz, m_cvc->T1());
  // bugout(0, _T("10 %f %f E1"), pxyz[0], pxyz[1]);
  // m_ccc->CalcPoint(pxyz, m_ccc->T0());
  // bugout(0, _T("10 %f %f S1"), pxyz[0], pxyz[1]);
  // m_ccc->CalcPoint(pxyz, m_ccc->T1());
  // bugout(0, _T("10 %f %f E1"), pxyz[0], pxyz[1]);

  // refine m_nose and m_tail calculations

  double lijk[3], tijk[3], xyz[2];

  double mclNomRatioLE = 0.001;//myGetProfileDouble(L"NomMclLEBackoff", 0.001);
  double mclNomRatioTE = 0.001;// myGetProfileDouble(L"NomMclTEBackoff", 0.001);
  if(dynamic_cast<Hexagon::Blade::HermiteOpenCurve*>(m_mc))
  {
    mclNomRatioLE = mclNomRatioTE = 1e-8;
  }

  double deltaLE = mclNomRatioLE*(m_mc->T1() - m_mc->T0());
  double deltaTE = mclNomRatioTE*(m_mc->T1() - m_mc->T0());

  m_mc->CalcPoint(xyz, m_mc->T0() + deltaLE, lijk);
  normalize(lijk, lijk);
  lijk[0] *= -1.0;
  lijk[1] *= -1.0;
  m_lec->NewLineIntersect(xyz, lijk, m_nose, 0.0, 0.0, &m_tnose, NULL, true);

  m_mc->CalcPoint(xyz, m_mc->T1() - deltaTE, tijk);
  //tijk[0] *= -1.0;
  //tijk[1] *= -1.0;
  normalize(tijk, tijk);
  m_tec->NewLineIntersect(xyz, tijk, m_tail, 0.0, 0.0, &m_ttail, NULL, true);

  m_lec->Extreme(m_tnose);
  m_tec->Extreme(m_ttail);

  //bugout(0, _T("7 %f %f nomsect final m_nose"), m_nose[0], m_nose[1]);
  //bugout(0, _T("7 %f %f nomsect final m_tail"), m_tail[0], m_tail[1]);

  return 1;
}

void CNominalSection::AdjustSubCurves(CCurve **tcvc, CCurve **tccc, CCurve **tlec, CCurve **ttec)
{
  double nose[2], tail[2], mid[2], nrm[2], cvmid[2], ccmid[2], xy[2];
  double d, d1, d2;
  double t0[4] = {0.0, 0.0, 0.0, 0.0};
  double t1[4] = {0.0, 0.0, 0.0, 0.0};

  m_whole->CalcPoint(nose, m_tnose);
  m_whole->CalcPoint(tail, m_ttail);

  //bugout(0, _T("10 %lf %lf N"), nose[0], nose[1]);
  //bugout(0, _T("10 %lf %lf T"), tail[0], tail[1]);

  mid[0] = 0.5*(nose[0] + tail[0]);
  mid[1] = 0.5*(nose[1] + tail[1]);
  //bugout(0, _T("10 %lf %lf Mid"), mid[0], mid[1]);

  nrm[0] = nose[1] - tail[1];
  nrm[1] = tail[0] - nose[0];
  normalize(nrm, nrm);

  m_cvc->LineIntersect(mid, nrm, cvmid, 0.0, 0.0);
  m_ccc->LineIntersect(mid, nrm, ccmid, 0.0, 0.0);
  //bugout(0, _T("10 %lf %lf CVMid"), cvmid[0], cvmid[1]);
  //bugout(0, _T("10 %lf %lf CCMid"), ccmid[0], ccmid[1]);

  d = d1 = dist(nose, cvmid);
  d2 = dist(nose, ccmid);
  if (d2 > d1)
    d = d2;

  m_whole->CircIntersect(m_nose, 0.5*d, xy, 0.0, m_tnose, &t0[LEC]);
  m_whole->CircIntersect(m_nose, 0.5*d, xy, m_tnose, m_period/2, &t1[LEC]);

  d = d1 = dist(tail, cvmid);
  d2 = dist(tail, ccmid);
  if (d2 > d1)
    d = d2;

  m_whole->CircIntersect(m_tail, 0.75*d, xy, 0.0, m_tnose, &t0[LEC]);
  m_whole->CircIntersect(m_tail, 0.75*d, xy, m_tnose, m_period/2, &t1[LEC]);

  m_whole->CircIntersect(m_tail, 0.75*d, xy, m_period/2, m_ttail, &t0[TEC]);
  m_whole->CircIntersect(m_tail, 0.75*d, xy, m_ttail, m_period, &t1[TEC]);

  if (m_startCC)
  {
    t0[CVC] = t1[LEC];
    t1[CVC] = t0[TEC];
    t0[CCC] = t1[TEC];
    t1[CCC] = t0[LEC];
  }
  else
  {
    t0[CVC] = t1[TEC];
    t1[CVC] = t0[LEC];
    t0[CCC] = t1[LEC];
    t1[CCC] = t0[TEC];
  }

  *tcvc = new CSubCurve(m_whole, m_t0[CVC], m_t1[CVC], m_period);
  *tccc = new CSubCurve(m_whole, m_t0[CCC], m_t1[CCC], m_period);
  *tlec = new CSubCurve(m_whole, m_t0[LEC], m_t1[LEC], m_period);
  *ttec = new CSubCurve(m_whole, m_t0[TEC], m_t1[TEC], m_period);
}


void CNominalSection::MakeTolSegStore(int num)
{
  m_tolsegCount = num;
  m_start_point = new int[num];
  m_end_point = new int[num];
  m_curvature_coef = new int[num];
  m_pTol_start = new double[num];
  m_pTol_end = new double[num];
  m_mTol_start = new double[num];
  m_mTol_end = new double[num];
}