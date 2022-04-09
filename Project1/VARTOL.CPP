#include "stdafx.h"
#include "vartol.h"

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

CVariableTol::CVariableTol()
{
  for (int i=0; i<2; i++)
  {
    m_x[i] = m_y[i] = NULL;
    m_tol[i] = NULL;
    m_letol[i] = m_stacktol[i] = m_tetol[i] = m_chordtol[i] = -1.0e20;
  }

  m_nom = m_mcc = NULL;
  m_leoff = m_teoff = -1.0e20;
  m_tnose = m_ttail = -1.0e20;
  m_minouttol = m_maxouttol = 0.0;

  m_teSquare = false;
}

CVariableTol::~CVariableTol()
{
  for (int i=0; i<2; i++)
  {
    if (m_x[i])
      delete m_x[i];
    if (m_y[i])
      delete m_y[i];
    if (m_tol[i])
      delete m_tol[i];
  }
}

int CVariableTol::CalcTol(int nle, int nte, int nside)
{
  if (!m_nom || !m_mcc || m_leoff < -1.0e10 || m_letol[0] < -1.0e10)
    return 0;

  m_numPoints = nle + nte + 2*nside;
  m_x[0] = new double[m_numPoints];
  m_y[0] = new double[m_numPoints];
  m_x[1] = new double[m_numPoints];
  m_y[1] = new double[m_numPoints];

  // calculate all the important stuff on the nominal curve

  double period = m_nom->Period();

  // calculate coordinates of nose and tail point

  double pt_nose[2], pt_tail[2];
  m_nom->CalcPoint(pt_nose, m_tnose);
  m_nom->CalcPoint(pt_tail, m_ttail);
  //bugout(2, _T("10 %f %f N"), pt_nose[0], pt_nose[1]);
  //bugout(2, _T("10 %f %f T"), pt_tail[0], pt_tail[1]);

  // calculate transition points on mean camber line

  double pt_lemc[2], vc_ltan[2], vc_lnrm[2];
  //bugout(2, _T("about to m_mcc CircIntersect for nose m_leoff %f period %f"), m_leoff, m_nom->Period());
  if (!m_mcc->CircIntersect(pt_nose, m_leoff, pt_lemc, 0.0, 0.0, 0, vc_ltan))
  {
    double pt[2], vc[2];
    double delta = 0.004*(m_mcc->T1() - m_mcc->T0());
    m_mcc->CalcPoint(pt, m_mcc->T0() + delta, vc_ltan);
    vc[0] = pt_nose[0] - pt[0];
    vc[1] = pt_nose[1] - pt[1];
    double off = normalize(vc, vc) - m_leoff;
    pt_lemc[0] = pt[0] + off*vc[0];
    pt_lemc[1] = pt[1] + off*vc[1];
  }

  vc_lnrm[0] = vc_ltan[1];
  vc_lnrm[1] = -vc_ltan[0];
  normalize(vc_lnrm, vc_lnrm);

  double pt_temc[2], vc_ttan[2], vc_tnrm[2];

  // make sure that m_teoff is big enough, if not just change it.
  // Cutoff is 1.5 the distance between mcc end and tail minus the lower chordtol

  m_mcc->CalcPoint(pt_temc, m_mcc->T1());
  double d = dist(pt_temc, pt_tail);
  if (m_teoff < 1.5*d-m_chordtol[0])
    m_teoff = 1.5*d-m_chordtol[0];

  //bugout(2, _T("about to m_mcc CircIntersect for tail m_teoff %f"), m_teoff);
  if (!m_mcc->CircIntersect(pt_tail, m_teoff, pt_temc, 0.0, 0.0, 0, vc_ttan))
  {
    double pt[2], vc[2];
    double delta = 0.004*(m_mcc->T1() - m_mcc->T0());
    m_mcc->CalcPoint(pt, m_mcc->T1() - delta, vc_ttan);
    vc[0] = pt_tail[0] - pt[0];
    vc[1] = pt_tail[1] - pt[1];
    double off = normalize(vc, vc) - m_leoff;
    pt_temc[0] = pt[0] + off*vc[0];
    pt_temc[1] = pt[1] + off*vc[1];
  }

  vc_tnrm[0] = vc_ttan[1];
  vc_tnrm[1] = -vc_ttan[0];
  normalize(vc_tnrm, vc_tnrm);

  // calculate transition points on the nominal curve

  double pt_le0[2], t_le0, pt_le1[2], t_le1;

  //bugout(2, _T("about to m_nom lineIntersect for 1st nose"));
  if (!m_nom->LineIntersect(pt_lemc, vc_lnrm, pt_le0, 0.0, m_tnose, &t_le0))
    return 0;
  //bugout(2, _T("about to m_nom lineIntersect for 2nd nose"));
  if (!m_nom->LineIntersect(pt_lemc, vc_lnrm, pt_le1, m_tnose, 0.5*period, &t_le1))
    return 0;

  double pt_te0[2], t_te0, pt_te1[2], t_te1;

  //bugout(2, _T("about to m_nom lineIntersect for 1st tail"));
  if (!m_nom->LineIntersect(pt_temc, vc_tnrm, pt_te0, 0.5*period, m_ttail, &t_te0))
    return 0;
  //bugout(2, _T("about to m_nom lineIntersect for 2nd tail"));
  if (!m_nom->LineIntersect(pt_temc, vc_tnrm, pt_te1, m_ttail, period, &t_te1))
    return 0;

  //bugout(0, _T("10 %f %f le0 t %f"), pt_le0[0], pt_le0[1], t_le0);
  //bugout(0, _T("10 %f %f le1 t %f"), pt_le1[0], pt_le1[1], t_le1);
  //bugout(0, _T("10 %f %f te0 t %f"), pt_te0[0], pt_te0[1], t_te0);
  //bugout(0, _T("10 %f %f te1 t %f"), pt_te1[0], pt_te1[1], t_te1);
  // calculate transition points near stack point (intersection with y-axis or x-axis)


  double origin[2], axis[2], pt_side0[2], t_side0, pt_side1[2], t_side1;
  origin[0] = origin[1] = axis[0] = axis[1] = 0.0;

  if (fabs(pt_nose[1] - pt_tail[1]) > fabs(pt_nose[0] - pt_tail[0]))
    axis[0] = 1.0;
  else
    axis[1] = 1.0;

  // had problem where stack point was beyond the blade for some reason.  adjust here.

  origin[0] = 0.5*(pt_nose[0] + pt_tail[0]);
  origin[1] = 0.5*(pt_nose[1] + pt_tail[1]);

//bugout(0, _T("10 %f %f N"), pt_nose[0], pt_nose[1]);
//bugout(0, _T("10 %f %f T"), pt_tail[0], pt_tail[1]);
//bugout(0, _T("8 %f %f %f %f T"),  0.0, 0.0, axis[0]*0.2, axis[1]*0.2);

  //bugout(2, _T("about to m_nom lineIntersect for 1st stack"));
  if (!m_nom->LineIntersect(origin, axis, pt_side0, m_ttail, m_tnose + period, &t_side0))
    return 0;
  //bugout(2, _T("about to m_nom lineIntersect for 2nd stack"));
  if (!m_nom->LineIntersect(origin, axis, pt_side1, m_tnose, m_ttail, &t_side1))
    return 0;
  //bugout(0, _T("10 %f %f side0 t %f"), pt_side0[0], pt_side0[1], t_side0);
  //bugout(0, _T("10 %f %f side1 t %f"), pt_side1[0], pt_side1[1], t_side1);

  int n1 = nside/2;
  int n2 = nside - n1;

  int n = m_numPoints - n1; // point counter, 0 will be in middle of a side
  Ramp(n, n1,  t_te1,   m_tetol,    t_side0, m_stacktol, -1);
  n = 0;
  Ramp(n, n2,  t_side0, m_stacktol, t_le0,   m_letol, 0);
  Ramp(n, nle, t_le0,   m_letol,    t_le1,   m_letol, 0);
  Ramp(n, n1,  t_le1,   m_letol,    t_side1, m_stacktol, 0);
  Ramp(n, n2,  t_side1, m_stacktol, t_te0,   m_tetol, 1);
  //bugout(2, _T("about to stretch"));

  if (!Stretch(n, nte, t_te0, t_te1, m_tetol, pt_temc, pt_tail, m_chordtol))
    return 0;

  for (int s=0; s<2; s++)
  {
    m_tol[s] = new CNurbCurve(m_numPoints, m_x[s], m_y[s], 0, m_nom->IsEnglish(), 1, 0, 1, 0.0, 0);
    if (!m_tol[s] || !m_tol[s]->Valid())
      return 0;
  }

  //bugout(2, _T("all is okay"));
  return 1;
}

int CVariableTol::Ramp(int &n, int num, double t0, double *tol0, double t1, double *tol1, int bias)
{
  //bugout(0, _T("Ramp t0 %f t1 %f period %f (before)"), t0, t1, m_nom->Period());
  if (t1 < t0)
    t1 += m_nom->Period();
  //bugout(0, _T("Ramp t0 %f t1 %f period %f (after)"), t0, t1, m_nom->Period());

  // bias -1, means put more points at the start, bias 1 means put more at the end, bias 0 means evenly space

  // using num instead of num - 1, should leave the last point unfilled.
  // this should be picked up when the other segment is filled.
  // avoiding duplicates here.

  int nb = num/5;                     // number of points in bias zone
  int nc;                             // index of point to change spacing of t values
  double tw = t1 - t0;                // width of t's to cover this segment
  double twb = tw/20;                 // width of t's to cover bias segment
  double dt1, dt2;                    // spacing for first and second parts
  double mtolw = (tol1[0] - tol0[0]); // width of minus tol to cover this segment
  double ptolw = (tol1[1] - tol0[1]); // width of plus tol to cover this segment
  double mtolwb = mtolw/20;           // width of minus tol to cover bias segment
  double ptolwb = ptolw/20;           // width of plus tol to cover bias segment
  double dmtol1, dmtol2;              // amount to change minus tolerance by each time
  double dptol1, dptol2;              // amount to change plus tolerance by each time

  if (bias < 0)  // extra points at start
  {
    nc = nb;
    dt1 = twb/nb;
    dt2 = (tw - twb)/(num - nb);
    dmtol1 = mtolwb/(nb);
    dmtol2 = (mtolw - mtolwb)/(num - nb);
    dptol1 = ptolwb/(nb);
    dptol2 = (ptolw - ptolwb)/(num - nb);
  }
  else if (bias > 0)
  {
    nc = num - nb;
    dt1 = (tw - twb)/(num - nb);
    dt2 = twb/nb;
    dmtol1 = (mtolw - mtolwb)/(num - nb);
    dmtol2 = mtolwb/(nb);
    dptol1 = (ptolw - ptolwb)/(num - nb);
    dptol2 = ptolwb/(nb);
  }
  else
  {
    nc = num;
    dt1 = dt2 = tw/(num - 1);
    dmtol1 = dmtol2 = mtolw/num;
    dptol1 = dptol2 = ptolw/num;
  }

  double dmtol=dmtol1, dptol=dptol1, dt=dt1, t = t0;
  double mtol = tol0[0];
  double ptol = tol0[1];

  for (int i=0; i<num; i++)
  {
    if (i == nc)
    {
      dmtol = dmtol2;
      dptol = dptol2;
      dt = dt2;
    }
    double xy[2], tan[2], nrm[2];
    m_nom->CalcPoint(xy, t, tan);
    nrm[0] = tan[1];
    nrm[1] = -tan[0];
    normalize(nrm, nrm);
    m_x[0][n] = xy[0] + mtol*nrm[0];
    m_y[0][n] = xy[1] + mtol*nrm[1];
    m_x[1][n] = xy[0] + ptol*nrm[0];
    m_y[1][n] = xy[1] + ptol*nrm[1];

    n++;
    t += dt;
    mtol += dmtol;
    ptol += dptol;
  }

  return 1;
}



int CVariableTol::Stretch(int &n, int num, double t0, double t1, double *tol, double *mcpt, double *end, double *tolend)
{
  if (t1 < t0)
    t1 += m_nom->Period();

  double mcvec[2];  // mcl point to te vector
  mcvec[0] = end[0] - mcpt[0];
  mcvec[1] = end[1] - mcpt[1];
  double mcdist = normalize(mcvec, mcvec);
  if (mcdist < 1.0e-6)
    return 0;

  //bugout(0, _T("10 %f %f MC"), mcpt[0], mcpt[1]);
  //bugout(0, _T("10 %f %f E"), end[0], end[1]);
  //bugout(0, _T("mcdist %f tol %f %f tolend %f %f"), mcdist, tol[0], tol[1], tolend[0], tolend[1]);

  double dt = (t1 - t0)/num;

  double lasttan[2];
  double dum[2];
 // double mclNomRatioTE = myGetProfileDouble(L"NomMclTEBackoff", 0.001);
  double mclNomRatioTE = 0.001;
  double nudgeTE = 5*mclNomRatioTE*(m_mcc->T1() - m_mcc->T0());
  m_mcc->CalcPoint(dum, m_mcc->T1() - nudgeTE, lasttan);
  normalize(lasttan, lasttan);

  double epsilon = 0.005*(m_mcc->T1() - m_mcc->T0());

  int nsave = n;
  for (int l=0; l<2; l++)
  {
    n = nsave;

    for (int i=0; i<num; i++)
    {
      double t = t0 + i*dt;
      double xy[2], bxy[2], bt, tanv[2], nrm[2];

      // evaluate nominal point and determine normal vector
      m_nom->CalcPoint(xy, t, tanv);
      nrm[0] = tanv[1];
      nrm[1] = -tanv[0];
      normalize(nrm, nrm);

      // adjust for normal tolerance
      m_x[l][n] = xy[0] + tol[l]*nrm[0];
      m_y[l][n] = xy[1] + tol[l]*nrm[1];
      //bugout(0, _T("8 %f %f %f %f N l=%d i=%d"), xy[0], xy[1], m_x[l][n], m_y[l][n], l, i);

      // find closest point on mcl or on segment joining mcl termination and end point
      // only care about the tangent vector

      if (m_mcc->ClosestPoint(xy, bxy, &bt, tanv) && m_mcc->T1() - bt > epsilon)
      {
        normalize(tanv, tanv);
        //bugout(0, _T("closest point okay bt %f t0 %f t1 %f"), bt, m_mcc->T0(), m_mcc->T1());
        //bugout(0, _T("9 %f %f %f %f B %d l=%d"), xy[0], xy[1], bxy[0], bxy[1], i, l);
        //bugout(0, _T("8 %f %f %f %f"), bxy[0], bxy[1], bxy[0]+0.04*tanv[0], bxy[1]+0.04*tanv[1]);
      }
      else
      {
        //bugout(0, _T("closest point failed, using lasttan"));
        tanv[0] = lasttan[0];
        tanv[1] = lasttan[1];
        //bugout(0, _T("8 %f %f %f %f"), xy[0], xy[1], xy[0]+0.04*tanv[0], xy[1]+0.04*tanv[1]);
      }

      double pdist = projdist(mcpt, mcvec, xy);
      double ratio = pdist/mcdist;

      if (ratio < 0.0)
        ratio = 0.0;
      if (ratio > 1.0)
        ratio = 1.0;
      //bugout(0, _T("pdist %f ratio %f"), pdist, ratio);

      // now strech/compress

      //xy[0] = m_x[l][n];
      //xy[1] = m_y[l][n];

      m_x[l][n] += ratio*tolend[l]*tanv[0];
      m_y[l][n] += ratio*tolend[l]*tanv[1];

      if (ratio > 0.98 && m_teSquare)
      {
        double toend[2];
        toend[0] = xy[0] - end[0];
        toend[1] = xy[1] - end[1];
        double edist = normalize(toend, toend);
        ratio = 2.0*edist/ m_teWidth;
        if (ratio < 0.0)
          ratio = 0.0;
        if (ratio > 1.0)
          ratio = 1.0;
        m_x[l][n] = xy[0] + ratio*tol[l]*toend[0];
        m_y[l][n] = xy[1] + ratio*tol[l]*toend[1];
        m_x[l][n] += tolend[l]*lasttan[0];
        m_y[l][n] += tolend[l]*lasttan[1];
      }
      //bugout(0, _T("9 %f %f %f %f %d Stretch %s"), xy[0], xy[1], m_x[l][n], m_y[l][n], n, l ? CString("Plus") : CString("Minus"));
      //bugout(0, _T("10 %f %f %d Stretch %s"), m_x[l][n], m_y[l][n], n, l ? CString("Plus") : CString("Minus"));
      n++;
    }
  }

  return 1;
}

//
//int CVariableTol::CalcTol(CSection *sect)
//{
//  m_numPoints = sect->GetNumNomPoints();
//  if (m_numPoints < 20)
//    return 0;
//
//  CCurve *nomCurve = sect->NomCurve();
//  if (!nomCurve)
//    return 0;
//
//  m_x[0] = new double[m_numPoints];
//  m_y[0] = new double[m_numPoints];
//  m_x[1] = new double[m_numPoints];
//  m_y[1] = new double[m_numPoints];
//
//  int i, j, nseed = 800;
//  double bt, t0=0.0, t1=0.0;
//  double mid[2]={0.0, 0.0};
//  for (i=0; i<m_numPoints; i++)
//  {
//    double xy[3], mtol, ptol;
//    if (!sect->GetNomPoint(i, xy) || !sect->GetTols(i, &mtol, &ptol))
//      return 0;
//
//    double d, bxy[2], tanv[2], nrm[2];
//    d = nomCurve->ClosestPoint(xy, bxy, &bt, tanv, 0.0, 0.0, nseed);
//    if (d > 0.01)
//      nomCurve->ClosestPoint(xy, bxy, &bt, tanv, 0.0, 0.0, 1000);
//
//    nseed = -1;
//    if (i == 0)
//      t0 = bt;
//    else if (i == 1)
//      t1 = bt;
//
//    mid[0] += xy[0];
//    mid[1] += xy[1];
//
//    nrm[0] = tanv[1];
//    nrm[1] = -tanv[0];
//    normalize(nrm, nrm);
//
//    if (mtol == 0.0 && ptol == 0.0)  // look behind and look ahead, find most recent and closest future points with valid tolerances. 
//    {                                // Interpolate tolerances based on t values
//      double mtol0=0.0, ptol0=0.0, mtol1=0.0, ptol1=0.0, bt0, bt1, xy0[3], xy1[3];  // 3 because GetNomPoint returns k value with x,y.
//      for (j = 1; j<6; j++)  // find last good point
//      {
//        int k = (i + m_numPoints - j) % m_numPoints;
//        if (!sect->GetNomPoint(k, xy0) || !sect->GetTols(k, &mtol0, &ptol0))
//          return 0;
//        if (mtol0 < ptol0)
//          break;
//      }
//      if (j == 5) // didn't find one
//        return 0;
//
//      for (j = 1; j<6; j++)  // find next good point
//      {
//        int k = (i + j) % m_numPoints;
//        if (!sect->GetNomPoint(k, xy1) || !sect->GetTols(k, &mtol1, &ptol1))
//          return 0;
//        if (mtol1 < ptol1)
//          break;
//      }
//      if (j == 5) // didn't find one
//        return 0;
//
//      if (mtol0 == mtol1 && ptol0 == ptol1)
//      {
//        mtol = mtol0;
//        ptol = ptol0;
//      }
//      else
//      {
//        bt0 = bt1 = bt; // need this if using nseed.
//        // shouldn't have go too far a round a corner if tol has changed.
//        // may need to change nseed here if assumption is wrong. 
//        d = nomCurve->ClosestPoint(xy0, bxy, &bt0, NULL, 0.0, 0.0, nseed);
//        if (d > 0.01)
//          nomCurve->ClosestPoint(xy0, bxy, &bt0, NULL, 0.0, 0.0, 1000);
//
//        d = nomCurve->ClosestPoint(xy1, bxy, &bt1, NULL, 0.0, 0.0, nseed);
//        if (d > 0.01)
//          nomCurve->ClosestPoint(xy1, bxy, &bt1, NULL, 0.0, 0.0, 1000);
//
//        double period = nomCurve->Period();
//        if (bt1 - bt0 > period/2)
//          bt1 -= period;
//        else if(bt0 - bt1 > period/2)
//          bt0 -= period;
//
//        if (bt < bt0 && bt < bt1)
//          bt += period;
//        else if (bt > bt0 && bt > bt1)
//          bt -= period;
//
//        mtol = mtol0 + (bt - bt0)*(mtol1 - mtol0)/(bt1 - bt0);
//        ptol = ptol0 + (bt - bt0)*(ptol1 - ptol0)/(bt1 - bt0);
//      }
//    }
//
//    m_x[0][i] = xy[0]+nrm[0]*mtol;
//    m_y[0][i] = xy[1]+nrm[1]*mtol;
//
//    m_x[1][i] = xy[0]+nrm[0]*ptol;
//    m_y[1][i] = xy[1]+nrm[1]*ptol;
//
//    //bugout(0, _T("10 %f %f %d"), m_x[1][i], m_y[1][i], i);
//    //bugout(0, _T("9 %f %f %f %f %d"), m_x[0][i], m_y[0][i], m_x[1][i], m_y[1][i], i);
//  }
//
//  mid[0] /= m_numPoints;  // this is to push the start/end of the spline away from LE or TE if it is near
//  // otherwise may missing crossing points on ID tol.
//  mid[1] /= m_numPoints;
//  int bi=0;
//  double dmin = 1.0e10;
//  for (i=0; i<m_numPoints; i++)
//  {
//    double xy[3];
//    sect->GetNomPoint(i, xy);
//    double d = dist(xy, mid);
//    if (d < dmin)
//    {
//      dmin = d;
//      bi = i;
//    }
//  }
//
//  // check to see that first two points have ascending spline parameters, if not reverse the points before 
//  bool reverse = false;  // if the points in the nominal file are defined clockwise, will need to reverse the points prior to splining
//
//  // make sure both spline parameters are in [0,period)
//  if (t0 < 0.0)
//    t0 += nomCurve->Period();
//  if (t1 < 0.0)
//    t1 += nomCurve->Period();
//
//  if (t0 >= nomCurve->Period())
//    t0 -= nomCurve->Period();
//  if (t1 >= nomCurve->Period())
//    t1 -= nomCurve->Period();
//
//  if (t1 - t0 > 0.5*nomCurve->Period())  // points straddled the start/end point and points are going in the wrong direction
//    reverse = true;
//  else if (t1 - t0 < -0.5*nomCurve->Period())  // points straddled the start/end point and points are going in the right direction
//    reverse = false;
//  else if (t1 - t0 < 0.0)  // not straddling, in the wrong order
//    reverse = true;
//  else
//    reverse = false;
//
//  double *tmp = new double[m_numPoints];
//
//  for (i=0; i<m_numPoints; i++)
//  {
//    j = reverse?  m_numPoints - 1 - i : i;
//    tmp[j] = m_x[0][(bi+i)%m_numPoints];
//  }
//  for (i=0; i<m_numPoints; i++)
//    m_x[0][i] = tmp[i];
//
//  for (i=0; i<m_numPoints; i++)
//  {
//    j = reverse?  m_numPoints - 1 - i : i;
//    tmp[j] = m_y[0][(bi+i)%m_numPoints];
//  }
//  for (i=0; i<m_numPoints; i++)
//    m_y[0][i] = tmp[i];
//
//  for (i=0; i<m_numPoints; i++)
//  {
//    j = reverse?  m_numPoints - 1 - i : i;
//    tmp[j] = m_x[1][(bi+i)%m_numPoints];
//  }
//  for (i=0; i<m_numPoints; i++)
//    m_x[1][i] = tmp[i];
//
//  for (i=0; i<m_numPoints; i++)
//  {
//    j = reverse?  m_numPoints - 1 - i : i;
//    tmp[j] = m_y[1][(bi+i)%m_numPoints];
//  }
//  for (i=0; i<m_numPoints; i++)
//    m_y[1][i] = tmp[i];
//
//  for (int s=0; s<2; s++)
//  {
//    m_tol[s] = new CNurbCurve(m_numPoints, m_x[s], m_y[s], 0, nomCurve->IsEnglish(), 1, 0, 1, 0.0, 0);
//    if (!m_tol[s] || !m_tol[s]->Valid())
//      return 0;
//  }
//
//  return m_numPoints;
//}
//

int CVariableTol::CalcTol(int nlo, double *xlo, double *ylo, int nhi, double *xhi, double *yhi)
{
  if (nlo < 20 || nhi < 20)
    return 0;

  // need to chane order of points?
  bool swap = false;
  double xy[2],pt[2],t0,t1;
  xy[0] = xlo[0];
  xy[1] = ylo[0];
  m_nom->ClosestPoint(xy, pt, &t0, NULL, 0.0, 0.0, 700);

  t1 = t0;
  xy[0] = xlo[2];
  xy[1] = ylo[2];
  m_nom->ClosestPoint(xy, pt, &t1, NULL, 0.0, 0.0, -1);
  if (t1 < t0)
  {
    if (fabs(t1 - t0) < 0.5* m_nom->Period())  // not straddling the periodic point
      swap = true;
  }
  else
  {
    if (fabs(t1 - t0) > 0.5* m_nom->Period())  // is straddling the periodic point
      swap = true;
  }

  if (swap)
  {
    int h = nlo/2;
    int i;
    for (i=0; i<h; i++)
    {
      int j = nlo - 1 - i;
      t0 = xlo[i];
      xlo[i] = xlo[j];
      xlo[j] = t0;
      t0 = ylo[i];
      ylo[i] = ylo[j];
      ylo[j] = t0;
    }
  }

  swap = false;
  xy[0] = xhi[0];
  xy[1] = yhi[0];
  m_nom->ClosestPoint(xy, pt, &t0, NULL, 0.0, 0.0, 700);

  t1 = t0;
  xy[0] = xhi[2];
  xy[1] = yhi[2];
  m_nom->ClosestPoint(xy, pt, &t1, NULL, 0.0, 0.0, -1);
  if (t1 < t0)
  {
    if (fabs(t1 - t0) < 0.5* m_nom->Period())  // not straddling the periodic point
      swap = true;
  }
  else
  {
    if (fabs(t1 - t0) > 0.5* m_nom->Period())  // is straddling the periodic point
      swap = true;
  }

  if (swap)
  {
    int h = nhi/2;
    int i;
    for (i=0; i<h; i++)
    {
      int j = nhi - 1 - i;
      t0 = xhi[i];
      xhi[i] = xhi[j];
      xhi[j] = t0;
      t0 = yhi[i];
      yhi[i] = yhi[j];
      yhi[j] = t0;
    }
  }

  CNurbCurve *loCurve = new CNurbCurve(nlo, xlo, ylo, 0, m_nom->IsEnglish(), 1, 0, 0, 0.0, 0);
  if (!loCurve || !loCurve->Valid())
  {
    if (loCurve)
      delete loCurve;
    return 0;
  }

  CNurbCurve *hiCurve = new CNurbCurve(nhi, xhi, yhi, 0, m_nom->IsEnglish(), 1, 0, 0, 0.0, 0);
  if (!hiCurve || !hiCurve->Valid())
  {
    if (hiCurve)
      delete hiCurve;
    delete loCurve;
    return 0;
  }

  if (loCurve->compute_length() > hiCurve->compute_length())
  {
    //Can't happen, they must be in wrong order, so swap
    m_tol[0] = hiCurve;
    m_tol[1] = loCurve;
  }
  else
  {
    m_tol[0] = loCurve;
    m_tol[1] = hiCurve;
  }

  return nlo + nhi;  // anything besides 0 will do
}


int CVariableTol::OutTol(double *pt_act, double *pt_tol, double *outtol, int /*n*/)
{
  double d, vc_tan[2], vc_nrm[2], vc_ij[2];

  static double bt0, bt1;

  //int ns = n == 0 ? 600 : -1;  // not working??
  int ns = 2000;

  // check plus tol curve first, this is better behaved, info can help minus tol check
  d = m_tol[1]->ClosestPoint(pt_act, pt_tol, &bt1, vc_tan, 0.0, 0.0, ns);
  vc_nrm[0] = vc_tan[1];
  vc_nrm[1] = -vc_tan[0];
  vc_ij[0] = pt_act[0] - pt_tol[0];
  vc_ij[1] = pt_act[1] - pt_tol[1];
  normalize(vc_nrm, vc_nrm);
  normalize(vc_ij, vc_ij);
  if (dot(vc_ij, vc_nrm) > 0.0)  // outside upper tolerance curve
  {
    //bugout(0, _T("9 %f %f %f %f OUTPLUSTOL"), pt_act[0], pt_act[1], pt_tol[0], pt_tol[1]);
    //bugout(0, _T("10 %f %f T"), pt_tol[0], pt_tol[1]);
    //bugout(0, _T("Nrm %f %f tol->act %f %f"), vc_nrm[0], vc_nrm[1], vc_ij[0], vc_ij[1]);
    *outtol = d;
    if (*outtol > m_maxouttol)
      m_maxouttol = *outtol;
    return 1;
  }

  //d = m_tol[0]->ClosestPoint(pt_act, pt_tol, &bt0, vc_tan, 0.0, 0.0, ns);  try using direction determined above

  d = m_tol[0]->ClosestNominal(pt_act, vc_nrm, pt_tol, &bt0, vc_tan, 0., 0., ns);
  if (d > 1000.0)  // not usenorm of ClosestNominal failed
  {
    ns = 2000;
    d = m_tol[0]->ClosestPoint(pt_act, pt_tol, &bt0, vc_tan, 0., 0., ns);
  }
  else
    ns = -1;

  vc_nrm[0] = vc_tan[1];
  vc_nrm[1] = -vc_tan[0];
  vc_ij[0] = pt_act[0] - pt_tol[0];
  vc_ij[1] = pt_act[1] - pt_tol[1];
  normalize(vc_nrm, vc_nrm);
  normalize(vc_ij, vc_ij);
  if (dot(vc_ij, vc_nrm) < 0.0)  // inside lower tolerance curve
  {
    //bugout(0, _T("9 %f %f %f %f OUTMINUSTOL"), pt_act[0], pt_act[1], pt_tol[0], pt_tol[1]);
    *outtol = -d;
    if (*outtol < m_minouttol)
      m_minouttol = *outtol;
    return 1;
  }

  return 0;
}
