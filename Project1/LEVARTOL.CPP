#include "stdafx.h"
#include "levartol.h"

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

CLEVariableTol::CLEVariableTol()
{
  for (int i=0; i<2; i++)
  {
    m_x[i] = m_y[i] = NULL;
    m_tol[i] = NULL;
    m_letol[i] = -1.0e20;
  }

  m_nom = NULL;
  m_t0 = m_t1 = m_tnose = -1.0e20;
  m_minouttol = m_maxouttol = 0.0;
}

CLEVariableTol::~CLEVariableTol()
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

int CLEVariableTol::CalcTol(int nle)
{
  if (!m_nom || m_letol[0] < -1.0e10)
    return 0;

  m_numPoints = nle + 6;  //add 3 more points on each side.

  m_x[0] = new double[m_numPoints];
  m_y[0] = new double[m_numPoints];
  m_x[1] = new double[m_numPoints];
  m_y[1] = new double[m_numPoints];

  int n1 = nle/2;
  int n2 = nle - n1;

  double delta = (m_tnose - m_t0)/n1;
  double t = m_t0 - 3*delta;
  int i;
  for (i=0; i<3; i++)
  {
    double xy[2], tan[2], nrm[2];
    m_nom->CalcPoint(xy, t, tan);
    nrm[0] = tan[1];
    nrm[1] = -tan[0];
    normalize(nrm, nrm);
    m_x[0][i] = xy[0] + m_stacktol[0]*nrm[0];
    m_y[0][i] = xy[1] + m_stacktol[0]*nrm[1];
    m_x[1][i] = xy[0] + m_stacktol[1]*nrm[0];
    m_y[1][i] = xy[1] + m_stacktol[1]*nrm[1];
    t += delta;
  }

  delta = (m_t1 - m_tnose)/n2;
  t = m_t1;
  for (i=0; i<3; i++)
  {
    double xy[2], tan[2], nrm[2];
    m_nom->CalcPoint(xy, t, tan);
    nrm[0] = tan[1];
    nrm[1] = -tan[0];
    normalize(nrm, nrm);
    m_x[0][nle + 3 + i] = xy[0] + m_stacktol[0]*nrm[0];
    m_y[0][nle + 3 + i] = xy[1] + m_stacktol[0]*nrm[1];
    m_x[1][nle + 3 + i] = xy[0] + m_stacktol[1]*nrm[0];
    m_y[1][nle + 3 + i] = xy[1] + m_stacktol[1]*nrm[1];
    t += delta;
  }

  int n = 3;
  Ramp(n, n1, m_t0,    m_stacktol, m_tnose, m_letol);
  Ramp(n, n2, m_tnose, m_letol,    m_t1,    m_stacktol);

  for (int s=0; s<2; s++)
  {
    //for(int i=1;i<m_numPoints;i++)bugout(0, _T("9 %f %f %f %f tolseg s=%d i=%d"),m_x[s][i-1],m_y[s][i-1],m_x[s][i],m_y[s][i],s,i);
    m_tol[s] = new CNurbCurve(m_numPoints, m_x[s], m_y[s], 0, m_nom->IsEnglish(), 0, 0, 1, 0.0, 0);
    //bugout(2, _T("about to check m_tol[%d] and valid"), s);
    if (!m_tol[s] || !m_tol[s]->Valid())
      return 0;
  }

  return 1;
}

int CLEVariableTol::Ramp(int &n, int num, double t0, double *tol0, double t1, double *tol1)
{
  if (t1 < t0)
    t1 += m_nom->Period();
  //bugout(0, _T("Ramp t0 %f t1 %f period %f (after)"), t0, t1, m_nom->Period());

  // using num instead of num - 1, should leave the last point unfilled.
  // this should be picked up when the other segment is filled.
  // avoiding duplicates here.

  double tw = t1 - t0;                // width of t's to cover this segment
  double dt = tw/(num - 1);           // spacing between t values
  double mtolw = (tol1[0] - tol0[0]); // width of minus tol to cover this segment
  double ptolw = (tol1[1] - tol0[1]); // width of plus tol to cover this segment
  double dmtol = mtolw/num;           // amount to change minus tolerance by each time
  double dptol = ptolw/num;           // amount to change plus tolerance by each time

  double mtol = tol0[0];
  double ptol = tol0[1];
  double t = t0;

  for (int i=0; i<num; i++)
  {
    double xy[2], tan[2], nrm[2];
    m_nom->CalcPoint(xy, t, tan);
    nrm[0] = tan[1];
    nrm[1] = -tan[0];
    normalize(nrm, nrm);
    m_x[0][n] = xy[0] + mtol*nrm[0];
    m_y[0][n] = xy[1] + mtol*nrm[1];
    m_x[1][n] = xy[0] + ptol*nrm[0];
    m_y[1][n] = xy[1] + ptol*nrm[1];
    //bugout(0, _T("10 %f %f %d Ramp %s"), m_x[0][n], m_y[0][n], n, CString("Minus"));
    //bugout(0, _T("10 %f %f %d Ramp %s"), m_x[1][n], m_y[1][n], n, CString("Plus"));

    n++;
    t += dt;
    mtol += dmtol;
    ptol += dptol;
  }

  return 1;
}

int CLEVariableTol::OutTol(double *pt_act, double *pt_tol, double *outtol, int /*n*/)
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
  if (dot(vc_ij, vc_nrm) > 0.0)  // outside lower tolerance curve
  {
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
    *outtol = -d;
    if (*outtol < m_minouttol)
      m_minouttol = *outtol;
    return 1;
  }

  return 0;
}
