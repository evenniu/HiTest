#include "StdAfx.h"
#include "Nurb.h"
#include "Circle.h"
#include "OdrSpline.h"
#include "Sorter.h"
#include "SubCurve.h"
#include <io.h>
#include <share.h>
#include <fcntl.h>

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

#ifdef DIM3
#define PSIZE 4
#else
#define PSIZE 3
#endif

// during reporting, it can be useful to retrieve information about the curve;
// I call these annotations
void CNurbCurve::addDoubleAnnotation(const char* /*key*/, double /*value*/)
{
}
double CNurbCurve::getDoubleAnnotation(const char* /*key*/) const
{
  return std::numeric_limits<double>::quiet_NaN();
}
void CNurbCurve::addBoolAnnotation(const char* /*key*/, bool /*value*/)
{
}
bool CNurbCurve::getBoolAnnotation(const char* /*key*/) const
{
  return false;
}

int CNurbCurve::Allocate(int havek)  // allocate arrays
{
  m_kNum = (short)(m_Num + m_Order);
  m_nSpan = m_kNum - m_Order - 1;
  m_pSpan = m_Order - 1;
  m_mid = -1;
  m_good = 0;

  m_pknot = new double[(int)m_kNum];
  if (m_pknot == NULL)
    return 0;

  m_pbf = new CMatrix((int)m_Order, (int)m_Order);
  if (m_pbf == NULL)
    return 0;

  m_pndu = new CMatrix((int)m_Order, (int)m_Order);
  if (m_pndu == NULL)
    return 0;

  m_pa = new CMatrix((int)m_Order, (int)m_Order);
  if (m_pa == NULL)
    return 0;

  m_pleft = new double[(int)m_Order];
  if (m_pleft == NULL)
    return 0;

  m_pright = new double[(int)m_Order];
  if (m_pright == NULL)
    return 0;

  m_ppoints = new CMatrix((int)m_Num, PSIZE);
  if (m_ppoints == NULL)
    return 0;

  if (havek)
  {
    m_ppointk = new double[(int)m_Num];
    if (m_ppointk == NULL)
      return 0;
  }

  m_storesize = 5L*sizeof(short) + sizeof(long) + 6L*sizeof(double);
  m_storesize += (long)m_kNum*sizeof(double);    // knots
  m_storesize += (long)m_Num*PSIZE*sizeof(double);  // control points
  if (havek)
    m_storesize += (long)m_Num*sizeof(double);  // control points
  m_good = 1;

  m_min[0] = m_min[1] = 1.0;
  m_max[0] = m_max[1] = 0.0;

  return 1;
}

// constructor used by IGES file read

// constructor used by math file read

/* using int fd version instead */
CNurbCurve::CNurbCurve(FILE *fp, bool isEnglish) : CCurve()
{
  m_isEnglish = isEnglish;
  bool havek = false;
  short curveType;
  m_good = 0;
  m_pknot = 0;
  m_ppoints = 0;
  m_ppointk = 0;
  m_pbf = 0;
  m_pndu = 0;
  m_pa = 0;
  m_pleft = 0;
  m_pright = 0;

  fread(&curveType, sizeof(short), 1, fp);
  if (curveType >= 1000)
  {
    havek = true;
    curveType -= 1000;
  }

  if (curveType != NURB_TYPE)
    return;

  fread(&m_extreme,sizeof(double), 1, fp);
  fread(&m_Num,sizeof(short), 1, fp);
  fread(&m_Order,sizeof(short), 1, fp);

  if (!Allocate(havek))
    return;

  m_good = 0;   // still may fail

  fread(&m_kNum,sizeof(short), 1, fp);
  fread(&m_closed,sizeof(short), 1, fp);
  fread(&m_poly,sizeof(short), 1, fp);
  fread(&m_t0,sizeof(double), 1, fp);
  fread(&m_t1,sizeof(double), 1, fp);
  fread(m_min,sizeof(double), 2, fp);
  fread(m_max,sizeof(double), 2, fp);

  if (!fread(m_pknot,sizeof(double), m_kNum, fp))
    return;

  int i;
  for(i=0;i<m_Num;i++)
    if (!fread(m_ppoints->m[i],sizeof(double), PSIZE, fp))
      return;

  if (havek)
  {
    if (!fread(m_ppointk,sizeof(double),m_Num, fp))
      return;
  }

  if (m_closed)
    m_period = m_t1 - m_t0;

  m_good = 1;
}

// constructor used to derive curve from a set of points

CNurbCurve::CNurbCurve(int n, double *x, double *y, double *k, bool isEnglish, int cl, int tn, int fx,
                       double tc, CDumpRecordArray *splDmp, int indexLE, int indexTE, double mult, bool showDebug) : CCurve()
{
  m_isEnglish = isEnglish;
  m_good = 0;
  m_pknot = 0;
  m_ppoints = 0;
  m_ppointk = 0;
  m_pbf = 0;
  m_pndu = 0;
  m_pa = 0;
  m_pleft = 0;
  m_pright = 0;

  double *xval = NULL, *yval = NULL, *kval = NULL;  //points will be copied to these arrays before spline fit

  int index1 = -1, index2 = -1;
  if (indexLE > 0 && indexTE <= 0)
  {
    // only LE is partial
    index1 = indexLE;
  }
  else if (indexTE > 0 && indexLE <= 0)
  {
    // only TE is partial
    index1 = indexTE;
  }
  else if (indexLE > 0 && indexTE > 0)
  {
    // both ends are partial
    if (indexLE < indexTE)
    {
      index1 = indexLE;
      index2 = indexTE;
    }
    else
    {
      index1 = indexTE;
      index2 = indexLE;
    }
  }
  else
  {
    // neither end is partial
  }

  if (index1 < 1 && index2 <1)
  {
    // no phantom points need be created, fix all data at once
    if (fx)
    {
      int nlast = n;

      if (tc <= 0.0)
        tc = isEnglish ? 0.0005 : 0.0125; // about .0005 inches
      else
        tc = isEnglish ? tc : 25.4*tc;

      while (1)
      {
       // n = fix_data(n, x, y, k, cl, tc, cl ? false : true, showDebug);
        n = fix_data_new(n, x, y, k, cl, tc, cl ? false : true, showDebug, 0.1); // 0.707 

        if (n == nlast)  // break when no points were removed
          break;

        nlast = n;
      }
    }

    xval = new double[n];
    yval = new double[n];
    if (k)
      kval = new double[n];

    int i;
    for (i=0; i<n; i++)
    {
      xval[i] = x[i];
      yval[i] = y[i];
      if (k)
        kval[i] = k[i];

      if (splDmp)
      {
        DumpRecord *drec = new DumpRecord(xval[i], yval[i]);
        splDmp->push_back(drec);
      }
    }
  }
  else if (index2 > 0)
  {
    // both ends are partial, need to create two sets of phantom points
    // fix the data in three pieces
    int n1 = index1+1;
    int n2 = index2 - index1;
    int n3 = n - index2 - 1;

    int start2 = index1+1;
    int start3 = index2+1;

    if (fx)
    {
      int nlast = n1;

      if (tc <= 0.0)
        tc = isEnglish ? 0.0005 : 0.0125; // about .0005 inches
      else
        tc = isEnglish ? tc : 25.4*tc;

      while (1)
      {
        n1 = fix_data(n1, x, y, k, 0, tc, false, showDebug);

        if (n1 == nlast)  // break when no points were removed
          break;

        nlast = n1;
      }

      nlast = n2;

      while (1)
      {
        n2 = fix_data(n2, &x[start2], &y[start2], k ? &k[start2] : NULL, 0, tc, false, showDebug);

        if (n2 == nlast)  // break when no points were removed
          break;

        nlast = n2;
      }

      nlast = n3;

      while (1)
      {
        n3 = fix_data(n3, &x[start3], &y[start3], k ? &k[start3] : NULL, 0, tc, false, showDebug);

        if (n3 == nlast)  // break when no points were removed
          break;

        nlast = n3;
      }
    }

    n = n1 + n2 + n3 + 20;

    xval = new double[n];
    yval = new double[n];
    if (k)
      kval = new double[n];

    double tx[4],ty[4], px1[10], py1[10], px2[10], py2[10];

    tx[0] = x[n1 - 2];
    ty[0] = y[n1 - 2];
    tx[1] = x[n1 - 1];
    ty[1] = y[n1 - 1];

    tx[2] = x[start2];
    ty[2] = y[start2];
    tx[3] = x[start2 + 1];
    ty[3] = y[start2 + 1];

    InventEnd(1,2,4,tx,ty,10,px1,py1,mult);

    tx[0] = x[n1 + n2 - 2];
    ty[0] = y[n1 + n2 - 2];
    tx[1] = x[n1 + n2 - 1];
    ty[1] = y[n1 + n2 - 1];

    tx[2] = x[start3];
    ty[2] = y[start3];
    tx[3] = x[start3 + 1];
    ty[3] = y[start3 + 1];

    InventEnd(1,2,4,tx,ty,10,px2,py2,mult);

    int i, j=0;
    for (i=0; i<n1; i++)
    {
      xval[j] = x[i];
      yval[j] = y[i];
      if (k)
        kval[j] = k[i];

      if (splDmp)
      {
        DumpRecord *drec = new DumpRecord(xval[j], yval[j]);
        splDmp->push_back(drec);
      }
      j++;
    }

    for (i=0; i<10; i++)
    {
      xval[j] = px1[i];
      yval[j] = py1[i];
      if (k)
        kval[j] = i < 5 ? k[n1-1] : k[start2];
      j++;
    }

    for (i=0; i<n2; i++)
    {
      xval[j] = x[start2 + i];
      yval[j] = y[start2 + i];
      if (k)
        kval[j] = k[start2 +i];

      if (splDmp)
      {
        DumpRecord *drec = new DumpRecord(xval[j], yval[j]);
        splDmp->push_back(drec);
      }
      j++;
    }

    for (i=0; i<10; i++)
    {
      xval[j] = px2[i];
      yval[j] = py2[i];
      if (k)
        kval[j] = i < 5 ? k[n1+n2-1] : k[start3];
      j++;
    }

    for (i=0; i<n3; i++)
    {
      xval[j] = x[start3 + i];
      yval[j] = y[start3 + i];
      if (k)
        kval[j] = k[start3 + i];

      if (splDmp)
      {
        DumpRecord *drec = new DumpRecord(xval[j], yval[j]);
        splDmp->push_back(drec);
      }
      j++;
    }
  }
  else
  {
    // there is one partial end, fix the data in two pieces, all up to index point first, then all after
    int n1 = index1+1;
    int n2 = n - n1;
    int start = index1+1;

    if (fx)
    {
      int nlast = n1;

      if (tc <= 0.0)
        tc = isEnglish ? 0.0005 : 0.0125; // about .0005 inches
      else
        tc = isEnglish ? tc : 25.4*tc;

      while (1)
      {
        n1 = fix_data(n1, x, y, k, 0, tc, false, showDebug);

        if (n1 == nlast)  // break when no points were removed
          break;

        nlast = n1;
      }

      nlast = n2;

      while (1)
      {
        n2 = fix_data(n2, &x[start], &y[start], k ? &k[start] : NULL, 0, tc, false, showDebug);

        if (n2 == nlast)  // break when no points were removed
          break;

        nlast = n2;
      }
    }

    n = n1 + n2 + 10;

    xval = new double[n];
    yval = new double[n];
    if (k)
      kval = new double[n];

    double tx[4],ty[4], px[10], py[10];

    tx[0] = x[n1 - 2];
    ty[0] = y[n1 - 2];
    tx[1] = x[n1 - 1];
    ty[1] = y[n1 - 1];

    tx[2] = x[start];
    ty[2] = y[start];
    tx[3] = x[start + 1];
    ty[3] = y[start + 1];

    InventEnd(1,2,4,tx,ty,10,px,py,mult);

    int i, j=0;
    for (i=0; i<n1; i++)
    {
      xval[j] = x[i];
      yval[j] = y[i];
      if (k)
        kval[j] = k[i];

      if (splDmp)
      {
        DumpRecord *drec = new DumpRecord(xval[j], yval[j]);
        splDmp->push_back(drec);
      }
      j++;
    }

    for (i=0; i<10; i++)
    {
      xval[j] = px[i];
      yval[j] = py[i];
      if (k)
        kval[j] = i < 5 ? k[n1-1] : k[start];
      j++;
    }

    for (i=0; i<n2; i++)
    {
      xval[j] = x[start + i];
      yval[j] = y[start + i];
      if (k)
        kval[j] = k[start + i];

      if (splDmp)
      {
        DumpRecord *drec = new DumpRecord(xval[j], yval[j]);
        splDmp->push_back(drec);
      }
      j++;
    }
  }

  ODRSpline *os = new ODRSpline(n, xval, yval, kval, m_isEnglish, cl, tn, 0.0);

  delete [] xval;
  delete [] yval;
  if (k)
    delete [] kval;

  if (!os || !os->isGood())
    return;

  m_Num = (short)(os->GetNctl() + 1);
  m_Order = 4;
  m_kNum = (short)(m_Num + m_Order);

  if (!Allocate(k != 0))
    return;

  m_good = 1;

  m_min[0] = m_min[1] = 1.0;
  m_max[0] = m_max[1] = 0.0;

  m_closed = (short)cl;
  m_poly = 1;

  int i;
  for (i=0; i < m_kNum; ++i)
    m_pknot[i] = os->GetKnots()[i];

  for (i=0; i < m_Num; ++i)
  {
    m_ppoints->m[i][0] = os->GetCx()[i];
    m_ppoints->m[i][1] = os->GetCy()[i];
    m_ppoints->m[i][2] = 1.0;
    if (k)
      m_ppointk[i] = os->GetCk()[i];
  }

  m_t0 = os->GetKnots()[3];
  m_t1 = os->GetKnots()[os->GetNctl()+1];
  if (m_closed)
    m_period = m_t1 - m_t0;

  delete os;
}

void CNurbCurve::MeanCamber(int n, CCurve *ccc, CCurve *cvc, CCurve *lec, CCurve *tec, int useextreme, double *nomStart)
{
  n = n - 40;

  bugout(0, L"START MeanCamber n=%d", n);
  int i, ng, nle = 0,nte = 0;
  double ts,te,td,p1[3],p2[3];
  double lex[20],ley[20],tex[20],tey[20];

  m_kNum = 0;  // zero all this out, because it may not call Alloc
  m_good = 0;
  m_pknot = 0;
  m_ppoints = 0;
  m_ppointk = 0;
  m_pbf = 0;
  m_pndu = 0;
  m_pa = 0;
  m_pleft = 0;
  m_pright = 0;

  if (!ccc || !cvc)
  {
    return;
  }

  if (lec)
  {
    if (lec->Type() == POINT_TYPE)
      nle = 0;  // may not be on curve
    else if (lec->Type() == CIRCLE_TYPE)
    {
      lex[0] = ((CCircle *)lec)->X();
      ley[0] = ((CCircle *)lec)->Y();
      nle = 1;
    }
    else
    {
      if (lec->Valid())
        nle = lec->EndCamber(20,lex,ley,useextreme & 1);
    }
  }
  
  if (tec)
  {
    if (tec->Type() == POINT_TYPE)
      nte = 0;  // may not be on curve
    else if (tec->Type() == CIRCLE_TYPE)
    {
      tex[0] = ((CCircle *)tec)->X();
      tey[0] = ((CCircle *)tec)->Y();
      nte = 1;
    }
    else
    {
      if (tec->Valid())
        nte = tec->EndCamber(20,tex,tey,useextreme & 2);
    }
  }
 
  m_extreme = -1.0;
  m_palign = 0;

  if (n + nle + nte < 3)
  {
    return;
  }

  double *x = new double[n + nle + nte];

  if (!x)
  {
    return;
  }

  double *y = new double[n + nle + nte];

  if (!y)
  {
    delete [] x;
    return;
  }

  for (i = 0; i < nle; i++)
  {
    x[i] = lex[nle - i - 1];  // reverse sort, closest to edge first
    y[i] = ley[nle - i - 1];
  }

  // find end points that will work
  // note not called with s and e anywhere, consider yanking???

  ts = ccc->T0();
  te = ccc->T1();

  if (nle > 0)
  {
    ccc->CalcPoint(p1,ts);
    ccc->CalcPoint(p2,te);

    if (_hypot(lex[nle - 1] - p1[0],ley[nle - 1] - p1[1]) >
        _hypot(lex[nle - 1] - p2[0],ley[nle - 1] - p2[1]))
    {
      td = ts;
      ts = te;
      te = td;
    }
  }
  else if (nte > 0)
  {
    ccc->CalcPoint(p1,ts);
    ccc->CalcPoint(p2,te);

    if (_hypot(tex[0] - p2[0],tey[0] - p2[1]) >
        _hypot(tex[0] - p1[0],tey[0] - p1[1]))
    {
      td = ts;
      ts = te;
      te = td;
    }
  }
  else if (lec)
  {
    double ns[2];
    lec->CalcPoint(ns,lec->T0());
    ccc->CalcPoint(p1,ts);
    ccc->CalcPoint(p2,te);

    if (dist(ns, p1) > dist(ns, p2))
    {
      td = ts;
      ts = te;
      te = td;
    }
  }
  else if (tec)
  {
    double tl[2];
    tec->CalcPoint(tl,tec->T0());
    ccc->CalcPoint(p1,ts);
    ccc->CalcPoint(p2,te);

    if (dist(tl, p2) > dist(tl, p1))
    {
      td = ts;
      ts = te;
      te = td;
    }
  }
 else if (nomStart != NULL) // if we drop down to this case, we must be dealing with a section that has both ends partial
 {
   ccc->CalcPoint(p1,ts);
   ccc->CalcPoint(p2,te);
   if (dist(nomStart, p1) > dist(nomStart, p2))
   {
     td = ts;
     ts = te;
     te = td;
   }
 }

  int n1 = n/2;
  int n2 = n - n1;

  double thalf = 0.5*(ts + te);
  double tdelta = (te - ts)/(n-1);

  double *tx = new double[n];
  double *ty = new double[n];
  int ng1 = CamberPoints(n1, ccc, cvc, tx, ty, thalf, ts);

  for (i=0; i<ng1; i++)
  {
    x[nle + i] = tx[ng1-1-i];
    y[nle + i] = ty[ng1-1-i];
  }
 
  ng = nle + ng1;

  int ng2 = CamberPoints(n2, ccc, cvc, tx, ty, thalf + tdelta, te);

  for (i=0; i<ng2; i++)
  {
    x[ng + i] = tx[i];
    y[ng + i] = ty[i];
  }
  ng += ng2;

  delete [] tx;
  delete [] ty;

  if (ng - nle < 4)
  {
    delete [] x;
    delete [] y;
    return;
  }


  for (i = 0; i < nte; i++)
  {
    x[ng + i] = tex[i];
    y[ng + i] = tey[i];
  }

  ng += nte;
  //bugout(0, _T("ng = %d"), ng);

  if (ng < 4)
  {
    delete [] x;
    delete [] y;
    //bugout(0, L"mmc failed ng is %d", ng);
    return;
  }

  // There is a case where the section is severely curved that camber points are found that are actually outside of the section.
  // check to see if there are any and remove

  double len = _hypot(x[ng-1]-x[0],y[ng-1]-y[0]); // distance between first and last points

  // check first and second point to make sure first doesn't seem to be bad, if it is remove it
  {
    double d = _hypot(x[1] - x[0], y[1] - y[0]); // distance first two points
    if(d > 0.2 * len)
    {
      //double d2 = _hypot(x[1] - x[2], y[1] - y[2]); // distance second and third point looks like first point is bad
      {
        for(i = 1; i < ng; i++)
        {
          x[i] = x[i + 1];
          y[i] = y[i + 1];
        }
        ng--;
      }
    }
  }

  // first point should be okay, now loop through points looking for a wild one.
  for(i = 1; i < ng; i++)
  {
    double d = _hypot(x[i]-x[i-1],y[i]-y[i-1]); // distance this point and the previous
    double d2 = i+1<ng ? _hypot(x[i+1]-x[i-1], y[i+1]-y[i-1]) : 0.0; // distance between the next point and the previous point
    if (d > 0.2*len && d2 < 0.2*len)
    {
      int j;
      for (j=i; j<ng; j++)
      {
        x[j] = x[j+1];
        y[j] = y[j+1];
      }
      i--;
      ng--;
    }
  }

  // Sort points.  There are two extreme cases, very straight or very curved.
  // Points may be out of order from the algorithm, or because the le and te
  // points are on the wrong end.
  // First step is to find furthest points and create a line.
  // Second step is to find point that is farthest from this line.
  // if the line segents joining the exreme points and furthest point are
  // fairly co-linear, then do a linear sort, otherwise do a radial sort.

  // oops, found the bug, keep this in mind though.

  double lp[2], ep[2], lv[2];

  lp[0] = x[0];
  lp[1] = y[0];
  ep[0] = x[ng-1];
  ep[1] = y[ng-1];

  if (lec && lec->Type() == POINT_TYPE)  // use it in the sort
  {
    //lec->CalcPoint(lp, 0.0);  // messes up check for swapVector below???
  }
  //lv[0] = x[ng/2] - lp[0];   // for blades with big hook
  //lv[1] = y[ng/2] - lp[1];
  lv[0] = x[3*ng/4] - lp[0];   // for blades with big hook and hook on te
  lv[1] = y[3*ng/4] - lp[1];
  normalize(lv, lv);

  bool swapVector = false;

  if (lec)
  {
    double pt[2];
    lec->CalcPoint(pt, lec->T0());
    if (_hypot(pt[0]-x[0],pt[1]-y[0]) > _hypot(pt[0]-x[ng-1], pt[1]-y[ng-1]))
      swapVector = true;
  }
  else if (tec)
  {
    double pt[2];
    tec->CalcPoint(pt, tec->T0());
    if (_hypot(pt[0]-x[0],pt[1]-y[0]) < _hypot(pt[0]-x[ng-1], pt[1]-y[ng-1]))
      swapVector = true;
  }

  if (swapVector)
  {
    lv[0] *= -1.0;
    lv[1] *= -1.0;
    lp[0] = x[ng-1];
    lp[1] = y[ng-1];
  }

  // If this is a really bowed blade, we should use a circle to sort...
  // find point that projects closest onto the line that joins the start and end points
  double bestdist = 1.0e10;
  double halfdist = 0.5*dist(lp, ep);
  double dropdist = 0.0;
  double vec[2];
  vec[0] = ep[0] - lp[0];
  vec[1] = ep[1] - lp[1];
  normalize(vec, vec);

  int bi = 0;
  for (i=1; i<ng-1; i++)
  {
    double pt[2], prj[2];
    pt[0] = x[i];
    pt[1] = y[i];
    double thisdist = fabs(projdist(lp, vec, pt, prj) - halfdist);
    if (thisdist < bestdist)
    {
      bi = i;
      bestdist = thisdist;
      dropdist = dist(pt, prj);
    }
  }

  //bugout(9, L"10 %f %f S", lp[0], lp[1]);
  //bugout(9, L"10 %f %f E", ep[0], ep[1]);
  //bugout(9, L"10 %f %f M", x[bi], y[bi]);
  //bugout(9, L"halfdist %f dropdist %f ratio %f", halfdist, dropdist, dropdist/halfdist);

  if (dropdist/halfdist > 0.35)  // use a circle sort
  {
    double xx[3], yy[3];
    xx[0] = lp[0];
    xx[1] = x[bi];
    xx[2] = ep[0];
    yy[0] = lp[1];
    yy[1] = y[bi];
    yy[2] = ep[1];
    CCircle circ(3, xx, yy);

    //bugout(9, L"10 %f %f Ctr", circ.X(), circ.Y());
    CSorter sp(ng, x, y, circ.X(), circ.Y());
  }
  else
  {
    //bugout(9, L"9 %f %f %f %f sort line", lp[0], lp[1], lp[0]+lv[0], lp[1]+lv[1]);
    //bugout(9, L"10 %f %f S", lp[0], lp[1]);
    CSorter sp(ng, x, y, lp, lv);
  }

  int slide = 0;   // remove close points

  for (i=1; i<ng; i++)
  {
    if (slide > 0)  // some weeded out already
    {
      x[i] = x[i+slide];
      y[i] = y[i+slide];
    }

    // we could check change of direction in dx and dy here too.

    if (_hypot(x[i] - x[i-1], y[i] - y[i-1]) < 1.0e-4)
    {
      slide++;
      ng--;
      i--;  // in case 3 or more close ones in a row
    }
  }

  // if this was a circle sort, it might have reversed the order.  need to check this out, and reverse if it did

  double d0 = _hypot(lp[0] - x[0], lp[1] - y[0]);
  double dn = _hypot(lp[0] - x[ng - 1], lp[1] - y[ng-1]);
  if (dn < d0)
  {
    // need to reverse
    double temp;
    for (i=0; i < ng/2; i++)
    {
      temp = x[i];
      x[i] = x[ng-i-1];
      x[ng-i-1] = temp;

      temp = y[i];
      y[i] = y[ng-i-1];
      y[ng-i-1] = temp;
    }
  }

  // should we use tension, fix_dat and/or thin_tol

  //make .5 percent of total length the too close cutoff.
  double tooClose = 0.005*_hypot(x[0]-x[ng-1], y[0]-y[ng-1]);
  // NurbCurve interprets tooClose as inches.
  if (!ccc->IsEnglish())
    tooClose /= 25.4;

  //bugout(0, _T(" final ng = %d tooClose = %f"), ng, tooClose);
  //for (i=0; i<ng; i++)
  //  bugout(0, _T("10 %f %f %d"), x[i], y[i], i);

  int nlast = ng;

  while (1)
  {
    ng = fix_data(ng, x, y, NULL, 0, tooClose, true);

    //bugout(5, _T("in ODRSpline after fix_data loop, m_npts=%d"), m_npts);
    if (ng == nlast)
      break;
    nlast = ng;
  }

  ODRSpline *os = new ODRSpline(ng, x, y, 0, ccc->IsEnglish(), 0, 0, 0.0);
  if (!os || !os->isGood())
  {
    delete [] x;
    delete [] y;
    if (os)
      delete os;
    return;
  }

  m_Num = (short)(os->GetNctl() + 1);
  m_Order = 4;
  m_kNum = (short)(m_Num + m_Order);

  if (!Allocate(0))
  {
    return;
  }

  m_good = 1;

  m_min[0] = m_min[1] = 1.0;
  m_max[0] = m_max[1] = 0.0;

  m_closed = 0;
  m_poly = 1;

  for (i=0; i < m_kNum; ++i)
    m_pknot[i] = os->GetKnots()[i];

  for (i=0; i < m_Num; ++i)
  {
    m_ppoints->m[i][0] = os->GetCx()[i];
    m_ppoints->m[i][1] = os->GetCy()[i];
    m_ppoints->m[i][2] = 1.0;
  }

  m_t0 = os->GetKnots()[3];
  m_t1 = os->GetKnots()[os->GetNctl()+1];

  delete os;
  delete [] x;
  delete [] y;
}

// constructor to derive mean camber curve from component curves

CNurbCurve::CNurbCurve(int n, CCurve *ccc, CCurve *cvc, CCurve *lec, CCurve *tec, int useextreme, double *nomStart) : CCurve()
{
  m_isEnglish = ccc->IsEnglish();
  MeanCamber(n, ccc, cvc, lec, tec, useextreme, nomStart);
}

// constructor to derive mean camber curve from a closed curve

CNurbCurve::CNurbCurve(int n, double lextr, double textr, double ler, double ter,
                       CCurve *whole, int useextremes)
{
  m_isEnglish = whole->IsEnglish();
  double nose[2], tail[2], p[2], t[4];

  whole->CalcPoint(nose, lextr);
  whole->CalcPoint(tail, textr);

  double period = whole->Period();

  whole->CircIntersect(nose, 3.0*ler, p, 0.0, lextr, &t[0]);
  whole->CircIntersect(nose, 3.0*ler, p, lextr, 0.5*period, &t[1]);
  whole->CircIntersect(tail, 4.0*ter, p, 0.5*period, textr, &t[2]);
  whole->CircIntersect(tail, 4.0*ter, p, textr, period, &t[3]);

  CSubCurve *lec = new CSubCurve(whole, t[0], t[1], period);
  CSubCurve *ccc = new CSubCurve(whole, t[1], t[2], period);
  CSubCurve *tec = new CSubCurve(whole, t[2], t[3], period);
  CSubCurve *cvc = new CSubCurve(whole, t[3], t[0] + period, period);
  bugout(0, L"New CNurbCurve N=%d, cal MeanCamber***",n);
  //MeanCamber(n, ccc, cvc, lec, tec, useextremes);//会出现问题 先临时注释

  delete lec;
  delete ccc;
  delete tec;
  delete cvc;
}

CNurbCurve::~CNurbCurve()
{
  if (m_ppoints)
    delete m_ppoints;

  if (m_ppointk)
    delete []m_ppointk;

  if (m_pknot)
    delete []m_pknot;

  if (m_pbf)
    delete m_pbf;

  if (m_pndu)
    delete m_pndu;

  if (m_pa)
    delete m_pa;

  if (m_pleft)
    delete []m_pleft;

  if (m_pright)
    delete []m_pright;
}

int CNurbCurve::FindSpan(double u) const
{
  if(m_mid > 0 && (u >= m_pknot[m_mid]) && (u < m_pknot[m_mid + 1]))
    return m_mid;

  if(u >= m_pknot[m_nSpan + 1])
  {
    // need this if knots okay?

    while((m_pknot[m_nSpan + 1] == m_pknot[m_nSpan]) && (m_nSpan > m_pSpan))
      const_cast<CNurbCurve*>(this)->m_nSpan--;
    return m_nSpan;
  }

  int low = m_pSpan;
  int high = m_nSpan+1;
  int mid = (low + high)/2;

  while(1)
  {
    if(u < m_pknot[mid])
      high = mid;
    else if(u >= m_pknot[mid + 1])
      low = mid;
    else
      return mid;
    mid = (low + high) / 2;
  }

  //return mid; can't get here
}


void CNurbCurve::DersBasisFuns(int i, double u, int n) const
{
  int j,r;

  m_pndu->m[0][0] = 1.0;
  for (j=1;j<=m_pSpan;j++)
  {
    m_pleft[j] = u - m_pknot[i+1-j];
    m_pright[j] = m_pknot[i+j] - u;
    double saved = 0.0;
    for (r=0; r<j; r++)
    {
      m_pndu->m[j][r] = m_pright[r+1] + m_pleft[j-r];
      double temp = m_pndu->m[r][j-1]/m_pndu->m[j][r];

      m_pndu->m[r][j] = saved + m_pright[r+1]*temp;
      saved = m_pleft[j-r]*temp;
    }
    m_pndu->m[j][j] = saved;
  }

  for (j=0;j<=m_pSpan;j++)
    m_pbf->m[0][j] = m_pndu->m[j][m_pSpan];

  for (r=0; r<=m_pSpan; r++)
  {
    int s1 = 0;
    int s2 = 1;

    m_pa->m[0][0] = 1.0;
    int k;
    for(k=1; k<=n; k++)
    {
      double d = 0.0;
      int rk = r-k;
      int pk = m_pSpan-k;

      if (r >= k)
      {
        m_pa->m[s2][0] = m_pa->m[s1][0]/m_pndu->m[pk+1][rk];
        d = m_pa->m[s2][0]*m_pndu->m[rk][pk];
      }

      int j1,j2;

      if (rk >= -1)
        j1 = 1;
      else
        j1 = -rk;

      if (r-1 <= pk)
        j2 = k-1;
      else
        j2 = m_pSpan-r;

      for(j=j1; j<=j2; j++)
      {
        m_pa->m[s2][j] = (m_pa->m[s1][j] - m_pa->m[s1][j-1])/m_pndu->m[pk+1][rk+j];
        d += m_pa->m[s2][j]*m_pndu->m[rk+j][pk];
      }
      if (r <= pk)
      {
        m_pa->m[s2][k] = -m_pa->m[s1][k-1]/m_pndu->m[pk+1][r];
        d += m_pa->m[s2][k]*m_pndu->m[r][pk];
      }
      m_pbf->m[k][r] = d;
      j=s1;
      s1=s2;
      s2=j;
    }
  }

  r = m_pSpan;

  int k;
  for (k=1; k<=n; k++)
  {
    for(j=0; j<=m_pSpan; j++)
      m_pbf->m[k][j] *= r;
    r *= (m_pSpan-k);
  }
}

int CNurbCurve::CalcPoint(double *p, double t, double *tan, double *curv) const 
{
  int i, j, k;
  double ck[4][PSIZE];
  double ep = (m_t1 - m_t0)*1.0e-16;

  // Evaluate non-uniform basis functions

  if (m_period > 0.0)
  {
    if (t < m_t0)
      t += m_period;
    else if (t > m_t1)
      t -= m_period;
  }


  if (t < m_t0 + ep)
    t = m_t0 + ep;

  if (t > m_t1 - ep)
    t = m_t1 - ep;

  int span = FindSpan(t);

  DersBasisFuns(span, t, curv ? 2 : 1);  // calc basis funcs and derivatives

  int du = 0;
  if (tan)
    du = 1;
  if (curv)
    du = 2;

  int first = span - m_Order + 1;
  int numax = m_poly ? PSIZE : PSIZE-1;

  for (k=0; k<=du; k++)         // derivative loop
  {
    for(i=0;i<numax;i++)       // axis loop
    {
      ck[k][i] = 0.0;
      for(j=0;j<m_Order;j++)    // degree loop
        ck[k][i] += m_pbf->m[k][j]*m_ppoints->m[first + j][i];
    }

    if (!m_poly)
      for(i=0;i<PSIZE-1;i++)
        ck[k][i] /= ck[k][PSIZE-1];
  }

  p[0] = ck[0][0];
  p[1] = ck[0][1];

  if (m_palign)
    m_palign->MeasToBest(p, 1, p);

  if (tan)
  {
    tan[0] = ck[1][0];
    tan[1] = ck[1][1];
    if (m_palign)
      m_palign->MeasToBest(tan, 0, tan);
  }

  if (curv)
  {
    curv[0] = ck[2][0];
    curv[1] = ck[2][1];
    if (m_palign)
      m_palign->MeasToBest(curv, 0, curv);
  }

  return 1;
}

void CNurbCurve::Extent()
{
  if (m_min[0] <= m_max[0] && m_min[1] <= m_max[1])
    return;

  int i,j;
  double t,td,xyz[3];

  td = 0.01*(m_t1 - m_t0);

  for (i=j=0, t=m_t0; i<100; i++, t += td)
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

void CNurbCurve::Extent(double *min, double *max)
{
  Extent();

  min[0] = m_min[0];
  min[1] = m_min[1];
  max[0] = m_max[0];
  max[1] = m_max[1];
}

void CNurbCurve::Write(FILE *fp)
{
  int i;

  short curveType = NURB_TYPE;
  if (m_ppointk)
    curveType += (short)1000;
  fwrite(&curveType, sizeof(short), 1, fp);
  fwrite(&m_extreme,sizeof(double), 1, fp);

  fwrite(&m_Num,sizeof(short), 1, fp);
  fwrite(&m_Order,sizeof(short), 1, fp);
  fwrite(&m_kNum,sizeof(short), 1, fp);
  fwrite(&m_closed,sizeof(short), 1, fp);
  fwrite(&m_poly,sizeof(short), 1, fp);

  fwrite(&m_t0,sizeof(double), 1, fp);
  fwrite(&m_t1,sizeof(double), 1, fp);
  fwrite(m_min,sizeof(double), 2, fp);
  fwrite(m_max,sizeof(double), 2, fp);

  fwrite(m_pknot,sizeof(double), m_kNum, fp);

  for(i=0;i<m_Num;i++)
  {
    fwrite(m_ppoints->m[i],sizeof(double), PSIZE, fp);
  }

  if (m_ppointk)
  {
    fwrite(m_ppointk,sizeof(double), m_Num, fp);
  }
}

double CNurbCurve::TotalLength(int numPoints, double t0, double t1)
{
  // big hammer approximation

  if (numPoints < 500)
    numPoints = 500;

  double rv = 0.0;
  double delta = (t1 - t0)/(numPoints - 1);
  double p1[2], p2[2];

  for (int i=0; i<numPoints; i++)
  {
    CalcPoint(p2, t0 + i*delta);
    if (i > 0)
      rv += dist(p1, p2);
    p1[0] = p2[0];
    p1[1] = p2[1];
  }

  return rv;
}

void CNurbCurve::Walk(double &t, double d, double *pt, int numPoints)
{
  double length = 0.0, lastlength = 0.0;
  double delta = (T1() - T0())/(numPoints - 1);
  double p1[2], p2[2], t0 = t;
  p2[0] = p2[1] = 0.0;

  for (int i=0; i<numPoints; i++)
  {
    CalcPoint(p2, t0 + i*delta);
    if (i > 0)
    {
      length += dist(p1, p2);
      if (length >= d)  // overshot it
      {
        double ratio = (d - lastlength)/(length - lastlength);
        t = t0 + (i-1+ratio)*delta;
        CalcPoint(pt, t);
        return;
      }
      lastlength = length;
    }
    p1[0] = p2[0];
    p1[1] = p2[1];
  }

  // should never be here.
  pt[0] = p2[0];
  pt[1] = p2[1];
}

void CNurbCurve::Debug(int numPoints)
{
  double delta = (T1() - T0())/(numPoints - 1);
  double p1[2], p2[2], t0 = T0();

  for (int i=0; i<numPoints; i++)
  {
    CalcPoint(p2, t0 + i*delta);
    if (i > 0)
    {
    }
    p1[0] = p2[0];
    p1[1] = p2[1];
  }
}

#define MOMENT_X_FUNC  0
#define MOMENT_Y_FUNC  1
#define AREA_FUNC      2
#define LENGTH_FUNC    3

// Methods used to compute center of gravity.


/*-------------------------------------------------------------------*/
/*  Compute the product y(u) * y(u) * x'(u).                         */
/*-------------------------------------------------------------------*/

double CNurbCurve::moment_x_func(double u)
{
  double dFdu[2], F[2];
  CalcPoint(F, u, dFdu);
  return F[1] * F[1] * dFdu[0];
}

/*-------------------------------------------------------------------*/
/*  Compute the product x(u) * x(u) * y'(u).                         */
/*-------------------------------------------------------------------*/

double CNurbCurve::moment_y_func(double u)
{
  double dFdu[2], F[2];
  CalcPoint(F, u, dFdu);
  return F[0] * F[0] * dFdu[1];
}

/*-------------------------------------------------------------------*/
/*  Compute the product x(u) * y'(u) - y(u) * x'(u).                 */
/*-------------------------------------------------------------------*/

double CNurbCurve::area_func(double u)
{
  double dFdu[2], F[2];
  CalcPoint(F, u, dFdu);
  return F[0] * dFdu[1] - F[1] * dFdu[0];
}

/*-------------------------------------------------------------------*/
/*  Gaussian quadrature: Compute integral from a to b.               */
/*-------------------------------------------------------------------*/

double CNurbCurve::GaussQuad(int integrand, double a, double b)
{
  int i;
  double s, r, u, sum;
  static double nodes[5] = {.148874338981631211, .433395394129247191,
    .679409568299024406, .865063366688984511,
    .97390652851717172};
  static double weights[5] = {.29552422471475287, .269266719309996355,
    .219086362515982044, .149451349150580593,
    .0666713443086881376};

  if (a == b)
    return 0.;

  s = .5 * (b - a);
  r = .5 * (b + a);
  for (sum = 0., i = 0; i < 5; i++)
  {
    u = s * nodes[i];

    double v1=0.0, v2=0.0;

    switch (integrand)
    {
    case MOMENT_X_FUNC:
      v1 = moment_x_func(r + u);
      v2 = moment_x_func(r - u);
      break;
    case MOMENT_Y_FUNC:
      v1 = moment_y_func(r + u);
      v2 = moment_y_func(r - u);
      break;
    case AREA_FUNC:
      v1 = area_func(r + u);
      v2 = area_func(r - u);
      break;
    case LENGTH_FUNC:
      v1 = deriv_func(r + u);
      v2 = deriv_func(r - u);
    }

    sum += weights[i] * (v1 +v2);
  }

  return s * sum;

}

/*-------------------------------------------------------------------*/
/*  Compute the area enclosed by a 2D spline curve using Gaussian    */
/*  quadrature.                                                      */
/*-------------------------------------------------------------------*/

double CNurbCurve::ComputeArea(double t0, double t1)
{
  if (t0 == 0.0 && t1 == 0.0)
  {
    t0 = T0();
    t1 = T1();
  }

  int j1 = FindSpan(t0);
  int j2 = FindSpan(t1);

  if (j1 == j2)
    return GaussQuad(AREA_FUNC, t0, t1);

  double sum = GaussQuad(AREA_FUNC, t0, m_pknot[j1+1]);

  int i;
  for (i = j1 + 1; i < j2; i++)
  {
    sum += GaussQuad(AREA_FUNC, m_pknot[i], m_pknot[i+1]);
  }
  sum += GaussQuad(AREA_FUNC, m_pknot[j2], t1);

  return sum;
}

/*-------------------------------------------------------------------*/
/*  Compute the center of mass of the area enclosed by a 2D spline   */
/*  curve using Gaussian quadrature.                                 */
/*-------------------------------------------------------------------*/

void CNurbCurve::CenterOfMass(double *CM, bool isGEI)
{
  double t0 = T0();
  double t1 = T1();

  double A = ComputeArea(t0, t1);
  int j1 = FindSpan(t0);
  int j2 = FindSpan(t1);

  if (j1 == j2)
  {
    CM[0] = GaussQuad(MOMENT_X_FUNC, t0, t1) / A;
    CM[1] = GaussQuad(MOMENT_Y_FUNC, t0, t1) / A;
    return;
  }

  double Mx = GaussQuad(MOMENT_X_FUNC, t0, m_pknot[j1+1]);
  double My = GaussQuad(MOMENT_Y_FUNC, t0, m_pknot[j1+1]);
  int i;
  for (i = j1 + 1; i < j2; i++)
  {
    Mx += GaussQuad(MOMENT_X_FUNC, m_pknot[i], m_pknot[i+1]);
    My += GaussQuad(MOMENT_Y_FUNC, m_pknot[i], m_pknot[i+1]);
  }
  Mx += GaussQuad(MOMENT_X_FUNC, m_pknot[j2], t1);
  My += GaussQuad(MOMENT_Y_FUNC, m_pknot[j2], t1);

  // change to agree with reports that we have this wrong.
  if (isGEI)
  {
    CM[0] = Mx / A;
    CM[1] = My / A;
  }
  else
  {
    CM[1] = Mx / A;
    CM[0] = -My / A;
  }
}

bool CNurbCurve::SlidePerUnitWindow(double dbLength, int index, int &index2, double t0, double &t1)
{
  if (m_Num <= 0 || dbLength <= 0 || !m_ppoints || !m_pknot) 
    return false;

  if (m_Num - index <= 0)
    return false;

  t1 = arc_inverse(dbLength, t0, index2);
  index2 += index -2;

  return false;
}


//  Find the parameter t1 such that the arclength of a spline curve from t0 to t1 is equal to a given value length.

double CNurbCurve::arc_inverse(double length, double t0, int &index)
{
  int i;
  double f, fprime, delta, eps = 1.e-4;

  // should have total length of curve so can make a decent guess at starting tau

  double t1 = t0 + 0.1*m_period;

  for (i = 0; i < 30; i++)
  {
    fprime = deriv_func(t1);
    if (fprime == 0.)
      break;
    f = compute_length(t0, t1) - length;
    delta = f / fprime;
    t1 -= delta;
    if (fabs(delta) < eps)
      break;
  }

  index = FindSpan(t1);

  return t1;
}

//  Compute the norm of the tangent vector to a spline at t.

double CNurbCurve::deriv_func(double t)
{
  double p[2],  v[2];

  CalcPoint(p, t, v);

  return normalize(v, v);
}

//  Compute the length of arc from t[0] to u of a 2D or 3D spline curve using Gaussian quadrature.

double CNurbCurve::compute_length(double t0, double t1)
{
  if (t0 == 0.0 && t1 == 0.0)
  {
    t0 = T0();
    t1 = T1();
  }

  int j1 = FindSpan(t0);
  int j2 = FindSpan(t1);

  if (j1 == j2)
    return GaussQuad(LENGTH_FUNC, t0, t1);
   
  double mult = 1.0; 
  if (j1 > j2)
  {
    int si = j1;
    j1 = j2;
    j2 = si;
    double sd = t0;
    t0 = t1;
    t1 = sd;
    mult = -1.0;
  }

  double sum = GaussQuad(LENGTH_FUNC, t0, m_pknot[j1+1]);

  int i;
  for (i = j1 + 1; i < j2; i++)
  {
    sum += GaussQuad(LENGTH_FUNC, m_pknot[i], m_pknot[i+1]);
  }
  sum += GaussQuad(LENGTH_FUNC, m_pknot[j2], t1);

  return mult*sum;
}

