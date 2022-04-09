#include "stdafx.h"
#include "Alignment.h"

#pragma warning(disable : 4714)
#include <Eigen/Dense>

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

void CAlignment::Clear()
{
  int i,j;
  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
      m_mat[i][j] = i == j ? 1.0 : 0.0;
    m_borig[i] = 0.0;
    m_morig[i] = 0.0;
  }
}

CAlignment:: CAlignment()
{
  n = 2;
  Clear();
}

CAlignment::CAlignment(int d)
{
  n = d == 3 ? 3 : 2;
  Clear();
}

CAlignment::~CAlignment()
{}

CAlignment::CAlignment(CAlignment &obj)   // copy constructor
{
  n = obj.n;

  int i,j;
  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
      m_mat[i][j] = obj.m_mat[i][j];
    m_borig[i] = obj.m_borig[i];
    m_morig[i] = obj.m_morig[i];
  }
}

CAlignment& CAlignment::operator=(const CAlignment &obj) // assignment operator
{
  n = obj.n;

  int i,j;
  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
      m_mat[i][j] = obj.m_mat[i][j];
    m_borig[i] = obj.m_borig[i];
    m_morig[i] = obj.m_morig[i];
  }

  return *this;
}

bool CAlignment::Rotate(int z, double ang)
{
  if (z < 0 || z > 2)
    return false;

  int x = (z + 1)%3;
  int y = (z + 2)%3;

  double c = cos(ang);
  double s = sin(ang);

  double rmat[3][3], pmat[3][3];

  rmat[x][x] = c;
  rmat[x][y] = -s;
  rmat[x][z] = 0.0;
  rmat[y][x] = s;
  rmat[y][y] = c;
  rmat[y][z] = 0.0;
  rmat[z][x] = 0.0;
  rmat[z][y] = 0.0;
  rmat[z][z] = 1.0;

  int i,j;

  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      pmat[i][j] = m_mat[i][0]*rmat[0][j] + m_mat[i][1]*rmat[1][j] +
        m_mat[i][2]*rmat[2][j];
    }
  }

  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      m_mat[i][j] = pmat[i][j];

  return true;
}

// Builds rotational matrix where input vector becomes z axis.
// Note it should be reasonable close to begin with.

bool CAlignment::MatrixFromAxis(double *zaxis)
{
  double xaxis[3], yaxis[3];

  if (fabs(zaxis[0]) > fabs(zaxis[1]))
  {
    xaxis[0] = 1.0;
    xaxis[1] = 0.0;
    xaxis[2] = 0.0;

    normalize(zaxis, zaxis, 3);
    cross(zaxis, xaxis, yaxis);
    normalize(yaxis, yaxis, 3);
    cross(yaxis, zaxis, xaxis);
    normalize(xaxis, xaxis, 3);
  }
  else
  {
    yaxis[0] = 0.0;
    yaxis[1] = 1.0;
    yaxis[2] = 0.0;

    normalize(zaxis, zaxis, 3);
    cross(yaxis, zaxis, xaxis);
    normalize(xaxis, xaxis, 3);
    cross(zaxis, xaxis, yaxis);
    normalize(yaxis, yaxis, 3);
  }

  m_mat[0][0] = xaxis[0];
  m_mat[0][1] = xaxis[1];
  m_mat[0][2] = xaxis[2];
  m_mat[1][0] = yaxis[0];
  m_mat[1][1] = yaxis[1];
  m_mat[1][2] = yaxis[2];
  m_mat[2][0] = zaxis[0];
  m_mat[2][1] = zaxis[1];
  m_mat[2][2] = zaxis[2];

  return true;
}


void CAlignment::BestToMeas(double *best, int point, double *meas) const
{
  int i,j;
  double temp[3];

  for(i=0;i<n;i++)
    temp[i] = (point ? best[i] - m_borig[i] : best[i]);

  for(i=0;i<n;i++)
  {
    meas[i] = 0.0;
    for(j=0;j<n;j++)
      meas[i] += (temp[j]*m_mat[j][i]);
  }

  if (point)
    for(i=0;i<n;i++)
      meas[i] += m_morig[i];
}

void CAlignment::MeasToBest(const double *meas, int point, double *best) const
{
  int i,j;
  double temp[3];

  for(i=0;i<n;i++)
    temp[i] = (point ? meas[i] - m_morig[i] : meas[i]);

  for(i=0;i<n;i++)
  {
    best[i] = 0.0;
    for(j=0;j<n;j++)
      best[i] += (temp[j]*m_mat[i][j]);
  }

  if (point)
    for(i=0;i<n;i++)
      best[i] += m_borig[i];
}

void CAlignment::MeasToBestVector2Xd(double* meas, int point, double* best, ptrdiff_t numberOfElements) const
{
  //alwaysAssert(n == 2);

#ifdef _DEBUG
  Eigen::Vector2d testBest;
  MeasToBest(meas, point, testBest.data());
#endif

  Eigen::Map<const Eigen::Matrix2Xd> eigenMeas(meas, 2, numberOfElements);
  Eigen::Map<Eigen::Matrix2Xd> eigenBest(best, 2, numberOfElements);
  Eigen::Matrix2d mat;
  mat << m_mat[0][0], m_mat[0][1], m_mat[1][0], m_mat[1][1];
  if(point)
  {
    Eigen::Map<const Eigen::Vector2d> borig(m_borig);
    Eigen::Map<const Eigen::Vector2d> morig(m_morig);
    // below, the eval() prevents aliasing issues; it's probably unnecessary
    eigenBest = ((mat * (eigenMeas.colwise() - morig)).colwise() + borig).eval();
  }
  else
  {
    eigenBest = (mat * eigenMeas).eval(); // the eval() is to prevent aliasing issues; it's probably unnecessary
  }

#ifdef _DEBUG
  assert((testBest - eigenBest.col(0)).norm() < 1e-10);
#endif
}

/* If need this, change CString to WCHAR
void  CAlignment::BugOut(CString &str)
{
  bugout(0, _T("%s"), str);
  bugout(0, _T("B %9.4f %9.4f %9.4f"), m_borig[0], m_borig[1], m_borig[2]);
  bugout(0, _T("M %9.4f %9.4f %9.4f"), m_morig[0], m_morig[1], m_morig[2]);
  bugout(0, _T("  %9.4f %9.4f %9.4f"), m_mat[0][0], m_mat[0][1], m_mat[0][2]);
  bugout(0, _T("  %9.4f %9.4f %9.4f"), m_mat[1][0], m_mat[1][1], m_mat[1][2]);
  bugout(0, _T("  %9.4f %9.4f %9.4f"), m_mat[2][0], m_mat[2][1], m_mat[2][2]);
}*/

