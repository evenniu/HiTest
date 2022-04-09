#pragma once
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
struct CuppingResults
{
	bool cuppingResultIsFilledOut;

	double camberPointA[2];
	double surfacePointA[2];

	double camberPointB[2];
	double surfacePointB[2];

	double camberPointC[2];
	double surfacePointC[2];

	double circleCenter[2];
	double circleRadius;
};

class DLLEXPORT CSection
{
public: // as part of the big push to the dll, this was changed from public to private. make this  private and make
		// access functions where needed
	int m_leType;              // NORMAL, SQUARE or PARTIAL
	int m_teType;              // NORMAL, SQUARE or PARTIAL
	int m_skewReport;          // if skewed, 0=report in skew plane, 1=report in blade coords
	double m_zValue;           // height of section
	CAlignment* m_skewalign;   // alignment for skewed section
	wchar_t m_name[MAXBUFSZ];  // name of section
	CCurve* m_nomCurve;        // nominal closed curve pointer
	CCurve* m_meaCurve;        // measured closed curve pointer
	CCurve* m_BCCurve;         // measured ball center closed curve pointer
	CCurve* m_nomPart[5];      // nominal curve pieces, CV, CC, LE, TE, MC
	CCurve* m_meaPart[5];      // measured curve pieces, CV, CC, LE, TE, MC
	double m_nomPitch[3];      // x, y and d of nominal pitch circles
	double m_meaPitch[3];      // x, y and d of measured pitch circles
	CuppingResults m_cuppingResults[10];
	double m_leNomWid1[10][2]; // location of point 1 for le nominal width
	double m_leActWid1[10][2]; // location of point 1 for le actual width
	double m_teNomWid1[10][2]; // location of point 1 for te nominal width
	double m_teActWid1[10][2]; // location of point 1 for te actual width
	double m_leNomWid2[10][2]; // location of point 2 for le nominal width
	double m_leActWid2[10][2]; // location of point 2 for le actual width
	double m_teNomWid2[10][2]; // location of point 2 for te nominal width
	double m_teActWid2[10][2]; // location of point 2 for te actual width
	int m_totalPoints;         // total number of raw points
	int* m_partOf;             // point originally part of CVC, CCC, LEC or TEC
	double m_arcRangeCV[2];    // used by fullblade le arc fit, the range[s] of the one or two segments used in the fit;
	double m_arcRangeCC[2];    // used by fullblade le arc fit, the range[s] of the one or two segments used in the fit;

	double m_xy_chord[8][2];//nyc 1,2行表示交点 3，4行表示touch1和touch2的弦切线，5,6:表示旋转后的弦线， 7,8:表示平移之后的弦线
	double m_chdthck[2];//0:LECHDTHCK 1:LECHDTHCK
	double* m_mxpt; // raw x values of curve points
	double* m_mypt; // raw y values of curve points

	double* m_cxpt; // ball center x values of curve points
	double* m_cypt; // ball center y values of curve points
	double* m_czpt; // ball center z values of curve points

	double* m_ival; // i component of measured surface normal
	double* m_jval; // j component of measured surface normal
	double* m_kval; // k component of measured surface normal

	double* m_nxpt; // original nominal x values of curve points
	double* m_nypt; // original nominal y values of curve points

	double* m_nomt; // t parameter for nominal points

	int m_numNomPoints;
	const Hexagon::Blade::MeanCamberCurveParameters2016* nominalMCLParams;
	double* m_nomx;
	double* m_nomy;
	double* m_nomi;
	double* m_nomj;
	double* m_nomk;
	double* m_mtol;
	double* m_ptol;

	int m_tolsegCount; //公差段个数
					   // CStringArray tolsegNamearray; //公差段名称集合
	int* m_start_point;
	int* m_end_point;
	int* m_curvature_coef;
	CMatrix* m_tolSegPoints_start;
	CMatrix* m_tolSegPoints_end;
	double* m_pTol_start;
	double* m_pTol_end;
	double* m_mTol_start;
	double* m_mTol_end;
	double* m_t_start;
	double* m_t_end;
	BOOL m_tolCCW;

	int    m_openareaCount;         //开口区域个数最大3个
	double m_areaStartPoint[3][2]; //开口区域起始点坐标
	double m_areaEndPoint[3][2];   //开口区域终止点坐标
	int    m_areaIndex[3][2];      //开口区域起始-终止点索引
	int m_numBestFits;
	int m_fixedAxis;           //坐标轴固定: 0-X; 1-Y; 2-Z
public:
	bool m_notEnoughToleranced;

	//CCompRecordArray* m_compArray;

	CSection();
	virtual ~CSection();
	double NomT(int i)
	{
		return m_nomt[i];
	}
	int GetNumNomPoints()
	{
		return m_numNomPoints;
	}
	bool GetNomPoint(int i, double* xyk)
	{
		if (i < 0 || i >= m_numNomPoints)
			return false;
		xyk[0] = m_nomx[i];
		xyk[1] = m_nomy[i];
		xyk[2] = m_nomk[i];
		return true;
	}
	bool GetTols(int i, double* mtol, double* ptol)
	{
		if (i < 0 || i >= m_numNomPoints)
			return false;
		*mtol = m_mtol[i];
		*ptol = m_ptol[i];
		return true;
	}
	bool SetTols(int i, double mtol, double ptol)
	{
		if (i < 0 || i >= m_numNomPoints)
			return false;
		m_mtol[i] = mtol;
		m_ptol[i] = ptol;
		return true;
	}
	int GetNumMeaPoints()
	{
		return m_totalPoints;
	}
	bool GetMeaPoint(int i, double* xyz, double* ijk = 0)
	{
		if (i < 0 || i >= m_totalPoints)
			return false;
		xyz[0] = m_mxpt[i];
		xyz[1] = m_mypt[i];
		xyz[2] = m_zValue;
		if (ijk)
		{
			ijk[0] = m_ival[i];
			ijk[1] = m_jval[i];
			ijk[2] = m_kval[i];
		}
		return true;
	}

	void Name(wchar_t* name)
	{
		wcscpy_s(m_name, name);
	}

	wchar_t* Name()
	{
		return m_name;
	}

	void ZValue(const double zvalue)
	{
		m_zValue = zvalue;
	}
	double ZValue()
	{
		return m_zValue;
	}

	void LEType(const int letype)
	{
		m_leType = letype;
	}
	int LEType() const
	{
		return m_leType;
	}
	void TEType(const int tetype)
	{
		m_teType = tetype;
	}
	int TEType() const
	{
		return m_teType;
	}

	void NomCurve(CCurve* nc)
	{
		m_nomCurve = nc;
	}
	CCurve* NomCurve() const
	{
		return m_nomCurve;
	}
	void MeaCurve(CCurve* nc)
	{
		m_meaCurve = nc;
	}
	CCurve* MeaCurve() const
	{
		return m_meaCurve;
	}
	void BCCurve(CCurve* nc)
	{
		m_BCCurve = nc;
	}
	CCurve* BCCurve()
	{
		return m_BCCurve;
	}
	void NomPart(int i, CCurve* curve)
	{
		m_nomPart[i] = curve;
	}
	CCurve* NomPart(int i) const
	{
		return m_nomPart[i];
	}
	void MeaPart(int i, CCurve* curve)
	{
		m_meaPart[i] = curve;
	}
	CCurve* MeaPart(int i) const
	{
		return m_meaPart[i];
	}

	void NomPitch(const double* p)
	{
		m_nomPitch[0] = p[0];
		m_nomPitch[1] = p[1];
		m_nomPitch[2] = p[2];
	}
	const double* NomPitch();
	void MeaPitch(const double* p)
	{
		m_meaPitch[0] = p[0];
		m_meaPitch[1] = p[1];
		m_meaPitch[2] = p[2];
	}
	const double* MeaPitch();



	double FindKValue(double* pt);
	void MakeNomArrays(int numPts, const Hexagon::Blade::MeanCamberCurveParameters2016* mclParams);
	void AddNomXYIJK(int i, double* xyijk);
	void AddTol(int i, double* mptol);
	bool AssignPoints(double* xv, double* yv, int n, int* start, int* end);
	int Chord(int flg, double* lcp, double* tcp, double* lctr, double* tctr, double* ltv, double* ttv, double* m = 0,
		double* w = 0, double* tew = 0, double* zeroPt = 0);
};