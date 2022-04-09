#include "stdafx.h"
#include "SmallestCircle.h"

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

CSmallestCircle::CSmallestCircle(int n, double *x, double *y, double *leTarget, double *teTarget)
{
  m_numPoints = n;

  m_X = x;
  m_Y = y;

  m_leTarget[0] = leTarget[0];
  m_leTarget[1] = leTarget[1];

  m_teTarget[0] = teTarget[0];
  m_teTarget[1] = teTarget[1];

  FindCircle();
}

CSmallestCircle::CSmallestCircle(int n, double *x, double *y)  // quicker version.  Use when you only need approximate values
{
  m_circumference = 0.0;
  
  int i, ep1 = 0, ep2 = 0;
  double d = 0.0, dMax = 0.0;
  
  // first pass, compute circumference of points and determine point furthest from first point
  for(i=0; i<n; i++)
  {
    if (i > 0)
      d = sqrt((x[i] - x[i-1])*(x[i] - x[i-1]) + (y[i] - y[i-1])*(y[i] - y[i-1]));
    else
      d = sqrt((x[0] - x[n-1])*(x[0] - x[n-1]) + (y[0] - y[n-1])*(y[0] - y[n-1]));

    m_circumference += d;
    
    d = sqrt((x[0] - x[i])*(x[0] - x[i]) + (y[0] - y[i])*(y[0] - y[i]));
    if (d > dMax)
    {
      dMax = d;
      ep1 = i;
    }
  }

  // second pass, determine point furthest from extreme point found in first pass
  dMax = 0.0;
  for(i=0; i<n; i++)
  {
    d = sqrt((x[ep1] - x[i])*(x[ep1] - x[i]) + (y[ep1] - y[i])*(y[ep1] - y[i]));
    if (d > dMax)
    {
      dMax = d;
      ep2 = i;
    }
  }
  
  // don't need to be completely accurate here, so stop after two passes
  
  m_radius = 0.5*sqrt((x[ep1] - x[ep2])*(x[ep1] - x[ep2]) + (y[ep1] - y[ep2])*(y[ep1] - y[ep2]));
}

CSmallestCircle::~CSmallestCircle()
{
}


void CSmallestCircle::FindCircle()
{
  double dx0, dy0, dx1, dy1;

  //determine a point P with the smallest y value
  double vx, vy, ymin = 1.0e8;
  double P[2];
  P[0] = P[1] = 0.0;
  int i, pi=0;
  for (i = 0; i < m_numPoints; i++)
  {
    if (m_Y[i] < ymin)
    {
      ymin = m_Y[i];
      pi = i;
      P[0] = m_X[i];
      P[1] = m_Y[i];
    }
  }

  double px = P[0];
  double py = P[1];

  // find a point Q such that the angle of the line segment PQ with the x axis is minimal
  double dot_max = -1.0e20, dot;
  double Q[2];
  Q[0] = Q[1] = 0.0;
  int qi = 0;
  for (i = 0; i < m_numPoints; i++)
  {
    if (i == pi)
      continue;

    dx0 = m_X[i] - px;
    dy0 = m_Y[i] - py;

    dot = (dx0 < 0 ? -dx0 : dx0) / sqrt(dx0 * dx0 + dy0 * dy0);
    if (dot > dot_max)
    {
      dot_max = dot;
      Q[0] = m_X[i];
      Q[1] = m_Y[i];
      qi = i;
    }
  }

  double qx = Q[0];
  double qy = Q[1];

  double rx = 0.0, ry = 0.0;
  for (i = 0; i < m_numPoints; i++)
  {
    dot_max = -1.0e20;

    //find R such that the absolute value of the angle PRQ is minimal
    int j;
    double R[2];
    R[0] = R[1] = 0.0;
    int ri = 0;
    for (j = 0; j < m_numPoints; j++)
    {
      if (j == pi || j == qi)
        continue;

      vx = m_X[j];
      vy = m_Y[j];

      dx0 = px - vx;
      dy0 = py - vy;
      
      dx1 = qx - vx;
      dy1 = qy - vy;

      dot = (dx0 * dx1 + dy0 * dy1) / (sqrt(dx0 * dx0 + dy0 * dy0) * sqrt(dx1 * dx1 + dy1 * dy1));
      if (dot > dot_max)
      {				
        dot_max = dot;
        R[0]  = vx;
        R[1]  = vy;
        ri = j;
      }
    }

    rx = R[0];
    ry = R[1];

    //check for case 1 (angle PRQ is obtuse), the circle is determined
    //by two points, P and Q. radius = |(P-Q)/2|, center = (P+Q)/2
    if (dot_max < 0)
    {
      dx0 = px - qx;
      dy0 = py - qy;

      m_center[0] = (px + qx) / 2;
      m_center[1] = (py + qy) / 2;
      m_radius = sqrt(((dx0 * dx0) / 4) + ((dy0 * dy0) / 4));

      // assign P and Q to be LE and TE points, depending on which is closer to target

      double check1 = (px - m_leTarget[0])*(px - m_leTarget[0]) + (py - m_leTarget[1])*(py - m_leTarget[1]);
      double check2 = (qx - m_leTarget[0])*(qx - m_leTarget[0]) + (qy - m_leTarget[1])*(qy - m_leTarget[1]);

      if (check1 < check2)  // P is closer to LE
      {
        m_lePoint[0] = px;
        m_lePoint[1] = py;
        m_tePoint[0] = qx;
        m_tePoint[1] = qy;
      }
      else  // Q is closer to LE
      {
        m_lePoint[0] = qx;
        m_lePoint[1] = qy;
        m_tePoint[0] = px;
        m_tePoint[1] = py;
      }

      return;
    }

    //check if angle RPQ is acute
    dx0 = rx - px;
    dy0 = ry - py;

    dx1 = qx - px;
    dy1 = qy - py;

    dot = (dx0 * dx1 + dy0 * dy1) / (sqrt(dx0 * dx0 + dy0 * dy0) * sqrt(dx1 * dx1 + dy1 * dy1));

    // if angle RPQ is 
    if (dot < 0)
    {
      P[0] = R[0];
      P[1] = R[1];
      px = rx;
      py = ry;
      continue;
    }

    // angle PQR is acute ?
    dx0 = px - qx;
    dy0 = py - qy;

    dx1 = rx - qx;
    dy1 = ry - qy;

    dot = (dx0 * dx1 + dy0 * dy1) / (sqrt(dx0 * dx0 + dy0 * dy0) * sqrt(dx1 * dx1 + dy1 * dy1));

    if (dot < 0)
    {
      Q[0] = R[0];
      Q[1] = R[1];
      qx = rx;
      qy = ry;
      continue;
    }

    //all angles in PQR are acute; quit
    break;
  }

  double mPQx = 0.5*(px + qx);
  double mPQy = 0.5*(py + qy);
  double mQRx = 0.5*(qx + rx);
  double mQRy = 0.5*(qy + ry);

  double numer = -(-mPQy * ry + mPQy * qy + mQRy * ry - mQRy * qy - mPQx * rx + mPQx * qx + mQRx * rx - mQRx *qx);
  double denom =  (-qx * ry + px * ry - px * qy + qy * rx - py * rx + py * qx);

  double t = numer / denom;

  m_center[0] = -t * (qy - py) + mPQx;
  m_center[1] =  t * (qx - px) + mPQy;

  dx0 = m_center[0] - px;
  dy0 = m_center[1] - py;

  m_radius = sqrt(dx0 * dx0 + dy0 * dy0);

  // check P, Q and R to see which are closer to the target to assign LE and TE points

  double check1 = (px - m_leTarget[0])*(px - m_leTarget[0]) + (py - m_leTarget[1])*(py - m_leTarget[1]);
  double check2 = (qx - m_leTarget[0])*(qx - m_leTarget[0]) + (qy - m_leTarget[1])*(qy - m_leTarget[1]);
  double check3 = (rx - m_leTarget[0])*(rx - m_leTarget[0]) + (ry - m_leTarget[1])*(ry - m_leTarget[1]);

  if (check1 < check2 && check1 < check3)  // P is closer to LE
  {
    m_lePoint[0] = px;
    m_lePoint[1] = py;

    check1 = (qx - m_teTarget[0])*(qx - m_teTarget[0]) + (qy - m_teTarget[1])*(qy - m_teTarget[1]);
    check2 = (rx - m_teTarget[0])*(rx - m_teTarget[0]) + (ry - m_teTarget[1])*(ry - m_teTarget[1]);

    if (check1 < check2)
    {
      m_tePoint[0] = qx;
      m_tePoint[1] = qy;
    }
    else
    {
      m_tePoint[0] = rx;
      m_tePoint[1] = ry;
    }
  }
  else if (check2 < check1 && check2 < check3) // Q is closer to LE
  {
    m_lePoint[0] = qx;
    m_lePoint[1] = qy;

    check1 = (px - m_teTarget[0])*(px - m_teTarget[0]) + (py - m_teTarget[1])*(py - m_teTarget[1]);
    check2 = (rx - m_teTarget[0])*(rx - m_teTarget[0]) + (ry - m_teTarget[1])*(ry - m_teTarget[1]);
    if (check1 < check2)
    {
      m_tePoint[0] = px;
      m_tePoint[1] = py;
    }
    else
    {
      m_tePoint[0] = rx;
      m_tePoint[1] = ry;
    }
  }
  else  // R is closer to LE
  {
    m_lePoint[0] = rx;
    m_lePoint[1] = ry;

    check1 = (px - m_teTarget[0])*(px - m_teTarget[0]) + (py - m_teTarget[1])*(py - m_teTarget[1]);
    check2 = (qx - m_teTarget[0])*(qx - m_teTarget[0]) + (qy - m_teTarget[1])*(qy - m_teTarget[1]);
    if (check1 < check2)
    {
      m_tePoint[0] = px;
      m_tePoint[1] = py;
    }
    else
    {
      m_tePoint[0] = qx;
      m_tePoint[1] = qy;
    }
  }

}
