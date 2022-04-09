#include "StdAfx.h"
#include <Eigen/Dense>
#include "CurvePolygon.h"

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG


#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif


#include <HexagonGDT/ScalarOptimization.h>


/* method to compute radius of curvature
*
* double der1[2], der2[2], rc[2];
* nomenc->CalcPoint(p, ta, der1, der2);
* double num = der1[0]*der1[0] + der1[1]*der1[1];
* double den = der1[0]*der2[1] - der1[1]*der2[0];
* if (fabs(den) > 1.0e-9)
* {
*    rc[0] = -der1[1]*num/den;
*    rc[1] =  der1[0]*num/den;
* }
*/

// Some day let curve know measurement units so can set iteration tolerances
CCurve::CCurve()
{
  m_extreme = -1.0;
  m_period = -1.0;
  //m_palign = NULL;
} // default constructor

CCurve::~CCurve()
{
}

// find extreme point on curve (LE or TE) releative to a normal vector, v

int CCurve::Extreme(double *v, double *et, double *ep, double pt0, double pt1, int save, bool endOK)
{
  if (pt0 == 0.0 && pt1 == 0.0)
  {
    pt0 = T0();
    pt1 = T1();
  }

  // kludge, this is needed for chord length for blades with a lot of camber
  // but it cause problems for some blades in mean camber line calculations
  // so save is set to -1 for chordlengh calculations

  if (Type() != SUBCURVE_TYPE && m_period > 0.0 && save == -1)
  {
    pt0 -= m_period/8;
    pt1 += m_period/8;
  }

  // Note that Circles, Points and Subcurves overload this function

  int i;
  double bt0=0.0,bt1=0.0,td,xyz[2],tv[2],nv[2],de,d0 = 0.0,d1 = 0.0;
  double ep0[2];

  // for blades with extreme curve, there maybe be more than one initial solution.
  // if there are then select the one that is "more extreme".

  //normalize(nv,v);
  const int numberOfStartingPointsForSearching = 1000;
  td = (pt1 - pt0)/(numberOfStartingPointsForSearching - 1);
  bt1 = pt0 - td;
  int numSol = 0;
  double sd0[numberOfStartingPointsForSearching], sbt0[numberOfStartingPointsForSearching], sbt1[numberOfStartingPointsForSearching];
  for (i=0; i<numberOfStartingPointsForSearching; i++)  // find bracketing t's
  {
    bt0 = bt1;
    bt1 += td;
    CalcPoint(xyz, bt1, tv);
    if (i == 0)
    {
      ep0[0] = xyz[0];
      ep0[1] = xyz[1];
    }

    //normalize(tv,tv);
    d0 = d1;
   // d1 = dot(nv,tv);

    if (fabs(d1) < 1.0e-6) // landed right on it!
    {
      *et = bt1;
      if (save == 1)
        m_extreme = bt1;
      if (ep)
      {
        ep[0] = xyz[0];
        ep[1] = xyz[1];
      }
      return 1;
    }

    if (i == 0)
      continue;

    if (d0*d1 < 0.0)  // straddling it
    {
      sbt0[numSol] = bt0;  // save the solution, but look for more.
      sbt1[numSol] = bt1;
      sd0[numSol] = d0;
      numSol++;
    }
  }

  if (numSol == 0)   // no bracketing interval found
  {
    if (endOK)
    {
      //if (projdist(ep0, nv, tv) > projdist(xyz, nv, tv))  // point zero was more extreme
      //{
      //  *et = pt0;
      //  ep[0] = ep0[0];
      //  ep[1] = ep0[1];
      //}
      //else
      //{
      //  *et = pt1;
      //  ep[0] = xyz[0];
      //  ep[1] = xyz[1];
      //}
      return 1;
    }
    return 0;
  }
  else if (numSol == 1)  // only one bracket found, use that
  {
    bt0 = sbt0[0];
    bt1 = sbt1[0];
    d0 = sd0[0];
  }
  else  // multiple brackets found
  {
    int bi = 0;
    double bd = -1.0e10;

    for (i=0; i<numSol; i++)
    {
      double test1[2], test2[2];
      CalcPoint(test1, sbt0[i]);
      CalcPoint(test2, sbt1[i]);

      /*double d = projdist(xyz, v, test1);
      if (d > bd)
      {
        bi = i;
        bd = d;
      }

      d = projdist(xyz, v, test2);
      if (d > bd)
      {
        bi = i;
        bd = d;
      }*/
    }

    bt0 = sbt0[bi];
    bt1 = sbt1[bi];
    d0 = sd0[bi];
  }

  while (bt1 - bt0 > 1.0e-5)  // binary search
  {
    *et = 0.5*(bt0 + bt1);
    CalcPoint(xyz, *et, tv);
    //normalize(tv,tv);
    de = dot(nv,tv);
    if (fabs(de) < 1.0e-4) // close enough
    {
      if (save == 1)
        m_extreme = *et;
      if (ep)
      {
        ep[0] = xyz[0];
        ep[1] = xyz[1];
      }
      return 1;
    }
    if (d0*de < 0.0)   // solution in left half
    {
      bt1 = *et;
      // d1 = de; not used???  problem here???
    }
    else               // solution in right half
    {
      bt0 = *et;
      d0 = de;
    }
  }

  // shouldn't be here, but go with whatever
  if (ep)
  {
    ep[0] = xyz[0];
    ep[1] = xyz[1];
    if (save == 1)
      m_extreme = *et;
  }
  return 1;
}

int CCurve::MatchTan(double *tv, double *et, double *ep, double pt0, double pt1, bool debug)
{
  if (pt0 == 0.0 && pt1 == 0.0)
  {
    pt0 = T0();
    pt1 = T1();
  }

  int i;
  double bt0=0.0,bt1=0.0,td=0.0,xy[2],tanv[2],nv[2],de = 0.0,d0 = 0.0,d1 = 0.0;

  nv[0] = tv[1];
  nv[1] = -tv[0];
  //normalize(nv,nv);  // nv is unit vector normal to curve
  td = (pt1 - pt0)/15.0;
  bt1 = pt0 - td;
  for (i=0; i<16; i++)  // find bracketing t's
  {
    bt0 = bt1;
    bt1 += td;
    CalcPoint(xy, bt1, tanv);
   // normalize(tanv,tanv);
    d0 = d1;
   // d1 = dot(nv,tanv);
    if (dot(tv,tanv) < 0.0)  // pointing in wrong direction
      continue;

    if (fabs(d1) < 1.0e-6)   // landed right on it!!
    {
      *et = bt1;
      if (ep)
      {
        ep[0] = xy[0];
        ep[1] = xy[1];
      }
      if (debug)
     //   bugout(0, L"MatchTan landed right on it!");
      return 1;
    }

    if (i == 0)
      continue;

    if (d0*d1 < 0.0)         // straddling it
      break;
  }

  if (i == 16)   // not there
  {
    if (debug)
    //  bugout(0, L"MatchTan did not find bracket!");
    return 0;
  }

  while (bt1 - bt0 > 1.0e-7)  // binary search
  {
    *et = 0.5*(bt0 + bt1);
    CalcPoint(xy, *et, tanv);
  //  normalize(tanv,tanv);
 //   de = dot(nv,tanv);
    if (fabs(de) < 1.0e-5) // close enough
    {
      if (ep)
      {
        ep[0] = xy[0];
        ep[1] = xy[1];
      }

    //  if (debug)
   //     bugout(0, L"MatchTan found it.  dot product %f!", de);

      return 1;
    }
    if (d0*de < 0.0)   // solution in left half
    {
      bt1 = *et;
      // d1 = de; not used???  problem here???
    }
    else               // solution in right half
    {
      bt0 = *et;
      d0 = de;
    }
  }

  // shouldn't be here, but go with whatever
  if (ep)
  {
    ep[0] = xy[0];
    ep[1] = xy[1];
  }

  if (debug)
//    bugout(0, L"MatchTan fell through cracks!  Last de was %lf", de);
  return 1;
}

#define OFF_EDGE       1
#define SMALL_DEN      2
#define NO_CONVERGE    3

int ACamberPoint(CCurve *a, CCurve *b, double *cp, double *p1, double *p2,
                 double *r1, double ta, double tb0, double tb1, double *tb)
{
 
  return NO_CONVERGE;
}

int CamberPoints(int n, CCurve *a, CCurve *b, double *x, double *y, double ta0, double ta1, double tb0, double tb1,
                 int useextreme)
{
  int type = a->Type();
  if (type == CIRCLE_TYPE || type == POINT_TYPE)  // should not be here
    return 0;

  if (n < 1 || !a || !b)
    return 0;

  if (ta0 == 0.0 && ta1 == 0.0)
  {
    ta0 = a->T0();
    ta1 = a->T1();
  }

  if (tb0 == 0.0 && tb1 == 0.0)
  {
    tb0 = b->T0();
    tb1 = b->T1();
  }

  int samecurve = a == b;

  //bugout(8, _T("START CamberPoints in EndCamber n=%d ta0 %f ta1 %f tb0 %f tb1 %f"), n, ta0, ta1, tb0, tb1);
  double textr=0.0,ts,te,td,r,cp[3],p1[3],p2[3],extr[3];

  bool extrIsGood = false;

  if (samecurve && useextreme)
  {
    a->CalcPoint(extr, a->Extreme());
    textr = a->Extreme();
    extrIsGood = true;
  }

  if (samecurve)
  {
    double xxx[2], yyy[2];
    a->CalcPoint(xxx,ta0);
    b->CalcPoint(yyy,tb1);
    a->CalcPoint(xxx,ta1);
    b->CalcPoint(yyy,tb0);
  }
  else
  {
    double xxx[2], yyy[2];
    a->CalcPoint(xxx,ta0);
    a->CalcPoint(yyy,ta1);
    b->CalcPoint(xxx,tb0);
    b->CalcPoint(yyy,tb1);
  }
  ts = ta0;
  te = ta1;

  // if need to go backwards, ta0 and ta1 should be reversed in calling routine

  if (n < 2)
  {
    td = te - ts;
    ts += 0.9*td;
  }
  else
    td = (te - ts)/(n - 1);

  int i, ng, retry_end = 0, retry_start = 0;
  double t = ts - td;
  double nudge = 50.0;
  double lastr = -1.0;
  for (i=0,ng=0; i<n; i++)
  {
    t += td;
    int rv = ACamberPoint(a, b, cp, p1, p2, &r, t, tb0, tb1);
    switch (rv)
    {
    case 0:  // good
      // x[ng] = cp[0];  done below
      // y[ng] = cp[1];
      // ng++;
      // retry_end = 0;  // do we want to try for more?
      break;
    case OFF_EDGE:  // trying to drive off the edge
      {
        //bugout(8, L"OFF_EDGE");
        if (i == 0)   // missed first point, readjust on the fly
        {
          if (retry_start%5 == 1)
            nudge /= 2.0;
          if (retry_start == 50)
            return 0;
          ts += td/nudge;  // gosh I hope I was close, or this may be slow
          td = (te - ts)/(n - 1);
          t = ts - td;
          i--;
          retry_start++;
        }
        else if (ng > 0)   // missing at end
        {
          if (retry_end++ > 3)
            break;  // give up and don't save the point
          t -= 1.1*td;  // move back a tad
          td *= 0.90;
          i--;
        }

        continue;
      }
    default:        // problems, but not driving of edge
      //bugout(8, L"10 %f %f OE", cp[0], cp[1]);
      
      double d1 = dist(cp, p1);
      double d2 = dist(cp, p2);
      
      if (d1 < 1.0e-5)
        continue;
      
      d1 = d2/d1;  // expect this to be close to 1.
      if (d1 < 0.99 || d1 > 1.01)
        continue;      
      
      //bugout(8, _T("OTHER ERROR"));
      //continue;
      break;    // go with it, wavy surface?
    }

    if (retry_end > 4)
      break;

    if (samecurve && ng > 0)  // attempt to terminate at circle center
    {
      double len, nv[3];

      if (ng > 1)  // this will be at least the third point
      {
        double l1, l2, v1[2], v2[2];
        v1[0] = x[ng-1] - x[ng-2];
        v1[1] = y[ng-1] - y[ng-2];
        v2[0] = cp[0] - x[ng-1];
        v2[1] = cp[1] - y[ng-1];
       /* l1 = normalize(v1,v1);
        l2 = normalize(v2,v2);*/

      }

      if (useextreme)
      {
        // already calculated extr above
        nv[0] = extr[0] - cp[0];
        nv[1] = extr[1] - cp[1];
        //len = normalize(nv,nv);
        if (len < r)
        {
          x[ng] = cp[0] + (len - r)*nv[0];
          y[ng] = cp[1] + (len - r)*nv[1];

          ng++;
          lastr = r;
          break;
        }
      }
      else
      {
        if (extrIsGood && ng > 1)
        {
          if (dist(p1, p2) < 0.25*r) // the cv and cc points are too close
            break;

          nv[0] = cp[0] - x[ng-2]; // skipping a point in case too close
          nv[1] = cp[1] - y[ng-2];
          double vec[2];
          vec[0] = extr[0] - cp[0];
          vec[1] = extr[1] - cp[1];
          if (dot(vec, nv) < 0.0)  // the point is on other side of curve
            break;
        }
        else
        {
          nv[0] = p2[1] - p1[1];
          nv[1] = p1[0] - p2[0];
          if (((cp[0] - x[ng-1])*nv[0] + (cp[1] - y[ng-1])*nv[1]) < 0.0)
          {
            nv[0] *= -1;  // point towards end
            nv[1] *= -1;
          }
        }

        normalize(nv,nv);

        if (a->LineIntersect(cp,nv,extr,0.0, 0.0,&textr))
        {
          len = _hypot(extr[0] - cp[0],extr[1] - cp[1]);
          if (len < r)
            break;  // don't use this point, adjust at bottom

          int nleft = n - ng - 1;
          if (nleft > 1)
          {
            double gap = 0.5*(tb0 - ta1);
            te = textr - gap; // adjust extreme parameter approximation
            tb0 = textr + gap;
            td = (te - t)/(nleft - 1);
          }
          extrIsGood = true;
        }
      }
    }
    x[ng] = cp[0];
    y[ng] = cp[1];

    ng++;
    lastr = r;
    if (retry_end)  // backed up to pick up end point, so quit
      break;
  }

  if (samecurve && ng > 1 && extrIsGood)
  {
    double len = _hypot(extr[0] - x[ng-1],extr[1] - y[ng-1]);

    if (len - lastr > 1.0e-5) // need to adjust last point
    {
      for (int j = 0; j < 5; j++)
      {
        double vec[2];
        vec[0] = x[ng-1] - extr[0];
        vec[1] = y[ng-1] - extr[1];
        normalize(vec, vec);
        x[ng-1] = extr[0] + lastr*vec[0];
        y[ng-1] = extr[1] + lastr*vec[1];

        double r1, r2, p1_temp[2], p2_temp[2];
        cp[0] = x[ng-1];
        cp[1] = y[ng-1];
        r1 = a->ClosestPoint(cp, p1_temp, 0, 0, ts, textr, 50);
        r2 = b->ClosestPoint(cp, p2_temp, 0, 0, textr, tb1, 50);
        lastr = 0.5*(r1 + r2);
        len = _hypot(extr[0] - x[ng-1],extr[1] - y[ng-1]);
        if (len - lastr < 1.0e-5)
          break;
      }
    }

    if (useextreme)
    {
      // project last point onto line that connect prev point to extr

      double lp[2], lv[2], pt[2], proj[2];

      lp[0] = x[ng-2];
      lp[1] = y[ng-2];
      lv[0] = extr[0] - lp[0];
      lv[1] = extr[1] - lp[1];
      pt[0] = x[ng-1];
      pt[1] = y[ng-1];
      if (projdist(lp, lv, pt, proj) != 0.0)
      {
        x[ng-1] = proj[0];
        y[ng-1] = proj[1];
       // bugout(8, L"10 %f %f F", x[ng-1], y[ng-1]);
      }
    }
  }

  return ng;
}

int CCurve::EndCamber(int n, double *x, double *y, int useextreme)
{
  int i;
  int type = Type();
  //bugout(8, "EndCamber useextreme %d", useextreme);

  if (type == CIRCLE_TYPE || type == POINT_TYPE)  // should not be here
    return 0;

  if (n < 1 || !x || !y || !Valid())
    return 0;

  double pt0,pt1,te;
  double saveExtreme = m_extreme;  // need to store this for calculation, not sure if we need to revert it, but will to be safe.

  pt0 = T0();
  pt1 = T1();

  if (useextreme && m_extreme != -1.0)
  {
    te = m_extreme;
  }
  else
  {
    double endvec[3];
    double dir[2],mid[2],end[2],p0[3],p1[3];
    CalcPoint(p0, T0());
    CalcPoint(p1, T1());
    
    //bugout(0, _T("8 %f %f %f %f EndCamber T0 and T1 points"),p0[0],p0[1],p1[0],p1[1]);
    //bugout(0, _T("t0 %f t1 %f period %f"), T0(), T1(), m_period);
    
    CalcPoint(end,0.5*(T0() + T1()));
    endvec[0] = p0[1] - p1[1];  // assumption is curve is U-shaped
    endvec[1] = p1[0] - p0[0];
    normalize(endvec,endvec);

    // make sure endvec points out

    mid[0] = 0.5*(p0[0] + p1[0]);
    mid[1] = 0.5*(p0[1] + p1[1]);
    dir[0] = end[0] - mid[0];
    dir[1] = end[0] - mid[0];

    if (dot(dir,endvec) < 0.0)
    {
      endvec[0] *= -1.0;
      endvec[1] *= -1.0;
    }

    //if (!Extreme(endvec, &te, end, pt0, pt1, 0))
    if (!Extreme(endvec, &te, end, pt0, pt1, 1))  // save extreme point parameter
      return 0;
  }

  double td = .05*(pt1 - pt0);

  int rv = CamberPoints(n, this, this, x, y, pt0, te - td, te + td, pt1, useextreme);
    
  //for(i=0; i<rv; i++)
  //  bugout(0, _T("10 %lf %lf xy%d"), x[i], y[i], i);

  if (rv < n/2)  // maybe would be better to try from other side?
  {
    double *tryx = new double[n];
    double *tryy = new double[n];

    int rv2 = CamberPoints(n, this, this, tryx, tryy, pt1, te + td, te - td, pt0, useextreme);

    //for(i=0; i<rv2; i++)
    //  bugout(0, _T("10 %lf %lf txy%d"), x[i], y[i], i);

    if (rv2 > rv)
    {
      for (i=0; i<rv2; i++)
      {
        x[i] = tryx[i];
        y[i] = tryy[i];
      }
      rv = rv2;
    }

    delete [] tryx;
    delete [] tryy;
  }

  m_extreme = saveExtreme;  // restore extreme parameter

  return rv;
}

void CCurve::BugOut(int numPoints, double t0, double t1)
{
  double t, td, p1[2], p2[2];

  td = -1.0;

  if (t0 == 0.0 && t1 == 0.0)
  {
    t0 = T0();
    t1 = T1();
  }

  t = t0;
  td = (t1 - t0)/numPoints;

  CalcPoint(p1,t0);
  p2[0] = p2[1] = 0.9;

  int k;
  for(k=1; k<=numPoints; k++)
  {
    t += td;
    if (k % 2 == 1)
    {
      CalcPoint(p2,t);
      //bugout(0, L"8 %f %f %f %f", p1[0], p1[1], p2[0], p2[1]);
    }
    else
    {
      CalcPoint(p1,t);
      //bugout(0, L"9 %f %f %f %f", p1[0], p1[1], p2[0], p2[1]);
    }
  }
}

bool CCurve::isPeriodic() const
{
  return Type() != SUBCURVE_TYPE && m_period > 0.0;
}

void CCurve::parametricBounds(double* const outBounds) const
{
  outBounds[0] = T0();
  outBounds[1] = T1();
}
void CCurve::evaluate(const double* const t, const ptrdiff_t n, double* const outPoints,
                      double* const outDerivative, double* const outSecondDerivative) const
{
  CalcPoints(outPoints, t, n, outDerivative, outSecondDerivative);
}

double CCurve::t0() const
{
  return T0();
}

double CCurve::t1() const
{
  return T1();
}

double CCurve::period() const
{
  return m_period;
}
/// <summary>
/// 
/// </summary>
/// <param name="tgt">input</param>
/// <param name="bestxyz">output</param>
/// <param name="bestt">output</param>
/// <param name="tangent">output</param>
/// <param name="pt0"></param>
/// <param name="pt1"></param>
/// <param name="nseed"></param>
/// <returns></returns>
double CCurve::ClosestPoint(double *tgt, double *bestxyz, double *bestt, double *tangent, double pt0, double pt1, int nseed)
{
  // NOTE: subcurves, circles and points overload this function
  if (pt0 == 0.0 && pt1 == 0.0)
  {
    pt0 = T0();
    pt1 = T1();
    if (Type() != SUBCURVE_TYPE && m_period > 0.0)
    {
      pt0 -= m_period/8;
      pt1 += m_period/8;
      if (nseed < 0 && bestt)
      {
        if (*bestt >= T1())
          *bestt -= m_period;
        else if (*bestt < T0())
          *bestt += m_period;
      }
    }
  }

  int i;
  double t[31], d[31],dotdt[31],dott[31];
  double td,dv[3],tv[3];

  if (nseed < 0 && bestt && *bestt >= pt0 && *bestt <= pt1)
  {
    i = 0;
    t[i] = *bestt;   // seed algorithm with *bestt
    td = 0.99*(pt1 - pt0);
  }
  else
  {
    if (nseed < 4)
      nseed = 10;

    td = 0.99*(pt1 - pt0)/nseed;
    d[0] = 1.0e20;

    for(i=0;i<nseed;i++)      // test at nseed points to seed search
    {
      double ti = pt0 + 0.5*td + i*td;  // seeds inside domain

      CalcPoint(dv, ti);     // no tv here is okay, findtol does request tan

      dv[0] -= tgt[0];
      dv[1] -= tgt[1];

      double di = dv[0]*dv[0] + dv[1]*dv[1];

      if (di < d[0])
      {
        d[0] = di;
        t[0] = ti;
      }
    }
  }

  int smallest = 0;
  if (fabs(pt1 - pt0) < 1.0e-08)
  {
    CalcPoint(bestxyz, t[0], tv);

    dv[0] = bestxyz[0] - tgt[0];
    dv[1] = bestxyz[1] - tgt[1];
    i = 30;
  }
  else
    i = 0;

  double dt;

  int osc = 0;

  // search for closest point, break out when sufficiently close
  for( ; i<30; i++)
  {
    CalcPoint(bestxyz, t[i], tv);

    dv[0] = tgt[0] - bestxyz[0];
    dv[1] = tgt[1] - bestxyz[1];

    d[i] = dv[0]*dv[0] + dv[1]*dv[1];     // square of distance
    if (d[i] < d[smallest])
      smallest = i;
    dott[i] = tv[0]*tv[0] + tv[1]*tv[1];  // square of tangent length
    dotdt[i] = dv[0]*tv[0] + dv[1]*tv[1]; // dot product

    int estimate = 0;

    if (i > 0 && d[i] > d[i-1])  // this guess is worse than last
      estimate = 1;
    else
    {
      dt = dotdt[i]/dott[i];

      if (fabs(dt/(pt1 - pt0)) < 1.0e-07)
        break;

      if (t[i] == pt0 && dt < 0.0)  // trying to go off the end
        break;

      if (t[i] == pt1 && dt > 0.0)  // trying to go off the end
        break;

      t[i+1] = t[i] + dt;

      if (i - osc > 1)
      {
        double tmean = 0.25*(t[i-2] + t[i-1] + t[i] + t[i+1]);
        double tdev[4];
        tdev[0] = t[i-2] - tmean;
        tdev[1] = t[i-1] - tmean;
        tdev[2] = t[i] - tmean;
        tdev[3] = t[i+1] - tmean;
        if (tdev[0]*tdev[2] > 0.0 && tdev[1]*tdev[3] > 0.0) // oscillating
        {
          osc = i+1; // try average of last four, or ave of i-2 and i-3?
          t[i+1] = tmean;
        }
      }
      else if (fabs(dt) > td)
      {
        dt = dt > 0 ? 0.75*td : -0.75*td;
        t[i+1] = t[i] + dt;
      }

      if (t[i+1] < pt0)
        t[i+1] = pt0;
      if (t[i+1] > pt1)
        t[i+1] = pt1;
    }

    if (estimate)
    {
      double ang0 = 0.0, ang1 = 0.0;
      if (d[i-1] > 1.0e-16 && dott[i-1] > 1.0e-16)
      {
        ang0 = dotdt[i-1]/sqrt(d[i-1])/sqrt(dott[i-1]);
        ang0 = acos(ang0) - 0.5*M_PI;
      }
      if (d[i] > 1.0e-16 && dott[i] > 1.0e-16)
      {
        ang1 = dotdt[i]/sqrt(d[i])/sqrt(dott[i]);
        ang1 = acos(ang1) - 0.5*M_PI;
      }

      if (fabs(ang0 - ang1) < 1.0e-8)     // flat, take average
        t[i+1] = 0.5*(t[i] + t[i-1]);
      else                                // interpolate
        t[i+1] = (ang0*t[i] - ang1*t[i-1])/(ang0 - ang1);
    }
  }


  if (i == 30)   // fell out, tangent etc, actually from i = 29
  {
    CalcPoint(bestxyz, t[smallest], tv);
    i = smallest;
  }

  if (bestt)
    *bestt = t[i];

  if (tangent)
  {
    tangent[0] = tv[0];
    tangent[1] = tv[1];
  }

  if (Type() != SUBCURVE_TYPE && m_period > 0.0 && bestt)
  {
    if (*bestt >= T1())
      *bestt -= m_period;
    else if (*bestt < T0())
      *bestt += m_period;
  }

  return _hypot(dv[0], dv[1]);
}


double CCurve::NewClosestPoint(const double* tgt, double* bestxyz, double* bestt, double* tangent, double pt0,
                               double pt1, int nseed)
{
  return Hexagon::Blade::closestPoint(*this, tgt, bestxyz, bestt, tangent, pt0, pt1, nseed);
}

int CCurve::NewLineIntersect(double* xy, double* ij, double* sol, double t0 /*= 0.0*/, double t1 /*= 0.0*/,
                             double* tsol /*= 0*/, double* dist /*= 0*/, bool isRay /*= false*/,
                             bool useSol /*= false*/, const ptrdiff_t numberOfSeeds /*= 500*/) const
{
  return Hexagon::Blade::lineIntersection(*this, xy, ij, sol, t0, t1, tsol, dist, isRay, useSol, numberOfSeeds);
}
/// <summary>
/// 
/// </summary>
/// <param name="tgt">input</param>
/// <param name="norm">input</param>
/// <param name="bestxyz">output</param>
/// <param name="bestt">output</param>
/// <param name="tangent">output</param>
/// <param name="pt0"></param>
/// <param name="pt1"></param>
/// <param name="nseed"></param>
/// <param name="match"></param>
/// <returns></returns>
double CCurve::ClosestNominal(double *tgt, double *norm, double *bestxyz, double *bestt, double *tangent, double pt0,
                              double pt1, int nseed, int match)
{
  if (pt0 == 0.0 && pt1 == 0.0)
  {
    pt0 = T0();
    pt1 = T1();
    if (Type() != SUBCURVE_TYPE && m_period > 0.0)
    {
      pt0 -= m_period/8;
      pt1 += m_period/8;
    }
  }

  int i, badmatch = 0;
  double t[31], d[31],dotdt[31],dott[31];
  double td=0.0,dv[3],tv[3],norm2[2];

  bool determineSeed = true;
  d[1] = 1.0e20;
  t[1] = bestt ? *bestt : pt0;

  if (nseed < 0 && bestt && *bestt >= pt0 && *bestt <= pt1)  // go with input besst as the seed
  {
    i = 0;
    t[i] = *bestt;   // seed algorithm with *bestt
    td = 0.99*(pt1 - pt0);

    CalcPoint(dv, t[i], tv);

    norm2[0] = tv[1];
    norm2[1] = -tv[0];
    normalize(norm2, norm2);

    if (dot(norm, norm2) > 0.707)  // look for within 45 degrees
      determineSeed = false;  // let's go with this
    else
      nseed = 400;  // better look again
  }

  if (determineSeed)
  {
    if (nseed < 4)
      nseed = 10;

    td = 0.99*(pt1 - pt0)/nseed;
    d[0] = 1.0e20;

    for(i=0;i<nseed;i++)      // test at nseed points to seed search
    {
      double ti = pt0 + 0.5*td + i*td;  // seeds inside domain

      CalcPoint(dv, ti, tv);
      dv[0] -= tgt[0];
      dv[1] -= tgt[1];
      double di = dv[0]*dv[0] + dv[1]*dv[1];

      if (di < d[0])
      {
        d[1] = di;
        t[1] = ti;

        norm2[0] = tv[1];
        norm2[1] = -tv[0];
        normalize(norm2, norm2);

        double dp = dot(norm, norm2);

        if (dp > 0.250)  // look for within 75 degrees
        {
          d[0] = di;
          t[0] = ti;
        }
      }
    }

    if (d[0] > 1.0e19)
    {
      d[0] = d[1];
      t[0] = t[1];
      badmatch = 1;
    }
  }

  int smallest = 0;
  if (fabs(pt1 - pt0) < 1.0e-08)
  {
    //bugout(0, _T("RIGHT ON"));
    CalcPoint(bestxyz, t[0], tv);

    dv[0] = bestxyz[0] - tgt[0];
    dv[1] = bestxyz[1] - tgt[1];
    i = 30;
  }
  else
    i = 0;

  double dt;

  int osc = 0, bounce = 0;

  // search for closest point, break out when sufficiently close
  for( ; i<30; i++)
  {
    CalcPoint(bestxyz, t[i], tv);

    dv[0] = tgt[0] - bestxyz[0];
    dv[1] = tgt[1] - bestxyz[1];

    d[i] = dv[0]*dv[0] + dv[1]*dv[1];     // square of distance
    if (d[i] < d[smallest])
      smallest = i;
    dott[i] = tv[0]*tv[0] + tv[1]*tv[1];  // square of tangent length
    dotdt[i] = dv[0]*tv[0] + dv[1]*tv[1]; // dot product

    int estimate = 0;

    if (i > 0 && d[i] > d[i-1])  // this guess is worse than last
    {
      estimate = 1;
      if (++bounce > 5)
      {
        //bugout(0, _T("Oscillating"));
        badmatch = 1;
        break;
      }
    }
    else
    {
      dt = dotdt[i]/dott[i];

      if (fabs(dt/(pt1 - pt0)) < 1.0e-07)
        break;

      if (t[i] == pt0 && dt < 0.0)  // trying to go off the end
      {
        //bugout(0, _T("OFF BEGINNING"));
        badmatch = 1;
        break;
      }

      if (t[i] == pt1 && dt > 0.0)  // trying to go off the end
      {
        //bugout(0, _T("OFF END"));
        badmatch = 1;
        break;
      }

      t[i+1] = t[i] + dt;

      if (i - osc > 1)
      {
        double tmean = 0.25*(t[i-2] + t[i-1] + t[i] + t[i+1]);
        double tdev[4];
        tdev[0] = t[i-2] - tmean;
        tdev[1] = t[i-1] - tmean;
        tdev[2] = t[i] - tmean;
        tdev[3] = t[i+1] - tmean;
        if (tdev[0]*tdev[2] > 0.0 && tdev[1]*tdev[3] > 0.0) // oscillating
        {
          osc = i+1; // try average of last four, or ave of i-2 and i-3?
          t[i+1] = tmean;
        }
      }
      else if (fabs(dt) > td)
      {
        dt = dt > 0 ? 0.75*td : -0.75*td;
        t[i+1] = t[i] + dt;
      }

      if (t[i+1] < pt0)
        t[i+1] = pt0;
      if (t[i+1] > pt1)
        t[i+1] = pt1;
    }

    if (estimate)
    {
      double ang0 = 0.0, ang1 = 0.0;
      if (d[i-1] > 1.0e-16 && dott[i-1] > 1.0e-16)
      {
        ang0 = dotdt[i-1]/sqrt(d[i-1])/sqrt(dott[i-1]);
        ang0 = acos(ang0) - 0.5*M_PI;
      }
      if (d[i] > 1.0e-16 && dott[i] > 1.0e-16)
      {
        ang1 = dotdt[i]/sqrt(d[i])/sqrt(dott[i]);
        ang1 = acos(ang1) - 0.5*M_PI;
      }

      if (fabs(ang0 - ang1) < 1.0e-8)     // flat, take average
        t[i+1] = 0.5*(t[i] + t[i-1]);
      else                                // interpolate
        t[i+1] = (ang0*t[i] - ang1*t[i-1])/(ang0 - ang1);
    }
  }

  //bugout(0, _T(" i %d smallest %d"), i, smallest);
  if (i == 30)   // fell out, tangent etc, actually from i = 29
  {
    CalcPoint(bestxyz, t[smallest], tv);
    //bugout(0, _T("FELL OUT SMALLEST"));
    i = smallest;
  }

  if (bestt)
    *bestt = t[i];

  if (match)
  {
    norm2[0] = tv[1];
    norm2[1] = -tv[0];
    normalize(norm2, norm2);

    if (dot(norm, norm2) < 0.866)  // if not within 30 degrees try match tan
    {
      tv[0] = -norm[1];
      tv[1] = norm[0];
      double tsol, sol[2];
      if (MatchTan(tv, &tsol, sol, pt0, pt1, false))
      {
        CalcPoint(bestxyz, tsol, tv);
        if (bestt)
          *bestt = tsol;
      }
    }
  }

  if (tangent)
  {
    tangent[0] = tv[0];
    tangent[1] = tv[1];
  }

  norm2[0] = tv[1];
  norm2[1] = -tv[0];
  normalize(norm2, norm2);

  if (fabs(dot(norm, norm2)) < 0.866)  // if not within 30 degrees kick out an error
    badmatch = 1;

  if (!determineSeed && badmatch)
  {
    //bugout(3, _T("9 %f %f %f %f recall here"), tgt[0], tgt[1], bestxyz[0], bestxyz[1]);
    //bugout(3, _T("10 %f %f A"), tgt[0], tgt[1]);
  }


  return badmatch*10000.0 + _hypot(dv[0], dv[1]);
}

int CCurve::LineIntersect(double *xy, double *ij, double *sol, double pt0, double pt1,
                          double *tsol, double *dist, bool isRay, bool useSol)
{
  if (pt0 == 0.0 && pt1 == 0.0)
  {
    pt0 = T0();
    pt1 = T1();
  }

  int i;

  if (dist)
    *dist = 10000.0;

  // search for bracket where curl changes sign
  // first try is middle third

  double ta, tb, tc, curla, curlb, curlc, delta, point[2], point2[2];

  ta = pt0 + (pt1 - pt0)/3.0;
  CalcPoint(point, ta);
  curla = (point[0] - xy[0])*ij[1] - (point[1] - xy[1])*ij[0];
  tb = pt1 - (pt1 - pt0)/3.0;
  CalcPoint(point2, tb);
  curlb = (point2[0] - xy[0])*ij[1] - (point2[1] - xy[1])*ij[0];

  bool seedGood = true;

  if (curla*curlb > 0)  // no sign change, need to look closer
    seedGood = false;

  if (isRay)        // if this is really a ray intersection, make sure its right direction
    seedGood = false;

  if (useSol && *tsol > pt0 && *tsol < pt1)
  {
    double width = (pt1 - pt0)/20;

    ta = *tsol - width;
    CalcPoint(point, ta);
    curla = (point[0] - xy[0])*ij[1] - (point[1] - xy[1])*ij[0];

    tb = *tsol + width;
    CalcPoint(point2, tb);
    curlb = (point2[0] - xy[0])*ij[1] - (point2[1] - xy[1])*ij[0];

    if (curla*curlb < 0)
      seedGood = true;
    else
    {
      width = (pt1 - pt0)/100;

      ta = *tsol - width;
      CalcPoint(point, ta);
      curla = (point[0] - xy[0])*ij[1] - (point[1] - xy[1])*ij[0];

      tb = *tsol + width;
      CalcPoint(point2, tb);
      curlb = (point2[0] - xy[0])*ij[1] - (point2[1] - xy[1])*ij[0];
      if (curla*curlb < 0)
        seedGood = true;
    }
  }

  if (!seedGood)
  {
    tb = pt0;
    CalcPoint(point, tb);
    curlb = (point[0] - xy[0])*ij[1] - (point[1] - xy[1])*ij[0];
    int numCheck = isRay ? 45 : 12;
    delta = (pt1 - pt0)/(numCheck - 1);


    bool saved = false;
    double saveta=0.0, savetb=0.0, savecurla=0.0, savecurlb=0.0;
    double bestcheck = -2.0;
    for (i=0; i<numCheck; i++)  // find a decent seed
    {
      ta = tb;
      curla = curlb;

      tb += delta;
      CalcPoint(point, tb);
      curlb = (point[0] - xy[0])*ij[1] - (point[1] - xy[1])*ij[0];

      if (curla*curlb < 0.0)
      {
        // if ray intersection, favor solution that is in right direction...
        double checker = (point[0] - xy[0])*ij[0] + (point[1] - xy[1])*ij[1];

        if (isRay)
        {
          // this dot product should be closer to 1, but might be sampling points sparesly, just make sure it's not
          if ((checker > -0.5) && (checker > bestcheck))
          {
            bestcheck = checker;
            saveta = ta;
            savetb = tb;
            savecurla = curla;
            savecurlb = curlb;
            saved = true;
          }
        }
        else
          break;
      }
    }

    if (i == numCheck)
    {
      if (saved)
      {
        ta = saveta;
        tb = savetb;
        curla = savecurla;
        curlb = savecurlb;
      }
      else
      {
        //bugout(5, "LineIntersect failed, no bracket");
        return 0;
      }
    }
  }

  delta = 0.5*(tb - ta);
  tc = ta + delta;
  while (delta > 1.0e-8*(pt1 - pt0))
  {
    tc = ta + delta;

    CalcPoint(point, tc);
    curlc = (point[0] - xy[0])*ij[1] - (point[1] - xy[1])*ij[0];

    if (curlc == 0.0)
    {
      if (sol)
      {
        sol[0] = point[0];
        sol[1] = point[1];
      }

      if (tsol)
        *tsol = tc;
      if (dist)
        *dist = _hypot(xy[0] - point[0], xy[1] - point[1]);
      return 1;
    }

    if (curla*curlc < 0.0)
    {
      //curlb = curlc; these are actually never used
      //tb = tc;
    }
    else
    {
      curla = curlc;
      ta = tc;
    }

    delta *= 0.5;
  }

  if (sol)
  {
    sol[0] = point[0];
    sol[1] = point[1];
  }

  if (tsol)
    *tsol = tc;

  if (dist)
    *dist = _hypot(xy[0] - point[0], xy[1] - point[1]);

  return 1;
}


// Intersect a Circle with the spline
// Used on Mean Camber Line as part of Edge Width calculations.
// Note this should be called so that there is but one solution in the
// interval [pt0, pt1].
// if there the points corresponding to pt0 and pt1 are both further or both
// closer to xy, the end point that is closest to r away from is selected,
// but zero will be returned in case it matters.

int CCurve::CircIntersect(const double* xy, double r, double* sol, double pt0, double pt1, double* tsol,
                          double* tanv) const
{
  return Hexagon::Blade::circleIntersection(*this, xy, r, sol, pt0, pt1, tsol, tanv);
}

// find the high point on a curve relative to line defined by p and v.

int CCurve::HighPoint(double *v, double *et, double *ep, double pt0, double pt1, bool endOK)
{
  if (Type() == CIRCLE_TYPE || Type() == POINT_TYPE)  // should not be here
    return 0;

  if (pt0 == 0.0 && pt1 == 0.0)
  {
    pt0 = T0();
    pt1 = T1();
  }
  //double q1[2],q2[2];
  //CalcPoint(q1,pt0);
  //bugout(0, _T("10 %f %f H0"), q1[0], q1[1]);
  //CalcPoint(q2,pt1);
  //bugout(0, _T("10 %f %f H1"), q2[0], q2[1]);

  int i;
  double bestt=0.0,t,td,xyz[2],nv[2],maxd = -1.0e20;
  double orig[2] = {0.0, 0.0};  // point doesn't matter

  normalize(nv,v);
  td = (pt1 - pt0)/19.0;
  t = pt0;
  for (i=0; i<20; i++)  // find starting point to look
  {
    CalcPoint(xyz, t);
    double d = projdist(orig, nv, xyz);
    if (d > maxd)
    {
      maxd = d;
      bestt = t;
    }
    t += td;
  }

  double at0 = bestt - td;
  double at1 = bestt + td;
  if (at0 < pt0)
    at0 = pt0;
  if (at1 > pt1)
    at1 = pt1;

  //CalcPoint(q1,at0);
  //bugout(0, _T("10 %f %f Z0"), q1[0], q1[1]);
  //CalcPoint(q2,at1);
  //bugout(0, _T("10 %f %f Z1"), q2[0], q2[1]);
  //bugout(0, _T("8 %f %f %f %f Z1-Z2"), q1[0], q1[1], q2[0], q2[1]);

  return Extreme(nv, et, ep, at0, at1, 0, endOK);
}

// this is used to advance along the mcl, breaking it up into somewhat equal segments
// note that this expects to go from nose to tail.

bool CCurve::AdvanceByLength(double delta, double t0, double t1, double *pt, double *tsol)
{

  double tlen = t1 - t0;
  if (tlen < 1.0e-5)
    return false;

  double p0[2], start[2], p1[2];

  CalcPoint(p0, t0);
  CalcPoint(p1, t1);

  if (dist(p0, p1) < delta)
  {
    pt[0] = p1[0];
    pt[1] = p1[1];
    *tsol = t1;
    return false;
  }

  start[0] = p0[0];
  start[1] = p0[1];

  for (int i=0; i<30; i++)
  {
    *tsol = 0.5*(t0 + t1);
    CalcPoint(pt, *tsol);
    double d = dist(start, pt);

    if (fabs(d - delta) < 1.0e-5)
      break;

    if (d < delta)
      t0 = *tsol;
    else
      t1 = *tsol;
  }

  return true;
}

// implementations here instead of in header file to help avoid warnings

double CCurve::compute_length(double, double)
{
  return 0.0;
}

void CCurve::Extreme(double extr)
{
  m_extreme = extr;
}

double CCurve::TotalLength(int, double, double)
{
  return 0.0;
}

void CCurve::Walk(double &, double, double *, int)
{
  return;
}

double CCurve::arc_inverse(double, double, int &)
{
  return 0.0;
}


void CCurve::CenterOfMass(double *CM, bool)
{
  CM[0] = CM[1] = 0.0;
}

double CCurve::ComputeArea(double, double)
{
  return 0.0;
}

void CCurve::Align(CAlignment *al)
{
  m_palign = al;
}

void CCurve::IsEnglish(bool eng)
{
  m_isEnglish = eng;
}

bool CCurve::IsEnglish()
{
  return m_isEnglish;
}

double CCurve::Extreme() const
{
  return m_extreme;
}

double CCurve::Period()
{
  return m_period;
}

int CCurve::CalcPoints(double* outPoints, const double* inT, ptrdiff_t numberOfT, double* outTangents /*= nullptr*/,
                       double* outCurves /*= nullptr*/) const
{
  for(ptrdiff_t i = 0; i < numberOfT; i++)
  {
    double* xyz = outPoints ? outPoints + 2 * i : nullptr;
    double* tan = outTangents ? outTangents + 2 * i : nullptr;
    double* curv = outCurves ? outCurves + 2 * i : nullptr;
    if(!CalcPoint(xyz, inT[i], tan, curv))
      return 0;
  }
  return 1;
}

void CCurve::FindClosestTValues(double* outT, const double* inPoints, ptrdiff_t numberOfPoints)
{
  Hexagon::Blade::findClosestTValues(*this, outT, inPoints, numberOfPoints);
}
