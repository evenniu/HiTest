#pragma once
#include "Flavor.h"
#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

class DLLEXPORT CAnalysisSect
{
public:
	bool m_goodVectors;
	bool m_bigGap;
	int m_numPoints;
	int m_inose;
	int m_itail;
	int m_inter1;
	int m_inter2;
	int m_skewReport;
	int m_phantomIndexLE;
	int m_phantomIndexTE;
	wchar_t m_sectName[MAXBUFSZ];
	double* x;
	double* y;
	double* z;
	double* ox;
	double* oy;
	double* oz;
	size_t numberOfBallCenters;
	double* ballCenterX;
	double* ballCenterY;
	double* ballCenterZ;
	double* i;
	double* j;
	double* oi;
	double* oj;
	double m_nose[2];
	double m_tail[2];
	//CAlignment* skewalign;
	CAnalysisSect();
	~CAnalysisSect();
	CAnalysisSect(const CAnalysisSect& obj);               // copy constructor
	CAnalysisSect& operator= (const CAnalysisSect& obj);  // assignment operator
};

class DLLEXPORT CAnalysisCell
{
public:
	int m_show;
	int m_decimals;
	wchar_t m_label[MAXBUFSZ];  // only for extra dimensions
	RECT m_box;
	double m_nom;
	double m_act;
	double m_ltol;
	double m_utol;
	double m_outtol; // for variable tolerances
	CAnalysisCell();
	void EmptyBox()
	{
		m_box.bottom = m_box.left = m_box.right = m_box.top = 0;
	}
	bool PtInRect(int x, int y)
	{
		if (x < m_box.left || x > m_box.right)
			return false;
		if (y < m_box.top || y > m_box.bottom)
			return false;

		return true;
	}
};

class DLLEXPORT CAnalysis
{
public:
	CBladeError* m_error;     // last error holder
	wchar_t** m_tnames;
	wchar_t** m_traces;
	wchar_t** m_sectionNames;      // section names in report
	wchar_t** m_calcLabels;        // array of calculation labels use \n if two lines long
	int m_decAng;
	int m_decMea;

	int m_numSect;
	int m_numCalc;
	int m_numPlat;
	int m_numTraces;

	int m_lotID;                      // lotplot trace field indices and other lotplot stuff
	int m_lotSeq;
	int m_lotTotal;
	int m_lotSize;
	int m_numLotTrans;

	int m_refSect;                    // special section indices and flags
	int m_rootSect;
	int m_tipSect;
	int m_refRow;
	int m_rootRow;
	int m_tipRow;

	int m_numSuspicious;              //  number of sections with suspicious vectors

	bool m_good;                      // calcs okay?
	bool m_refChecked;
	bool m_rootChecked;
	bool m_tipChecked;
	bool m_bowChecked;
	bool m_refchdangchecked;
	bool m_tipLEChecked;
	bool m_rootLEChecked;
	bool m_bowLEChecked;

	double m_bowCorrect;
	double m_twistCorrect;
	double m_dispCorrect;
	double m_refchdangnom;
	double m_refchdangact;
	double m_extraTolerance;      // from round off;

	double m_rootXYZ[3];
	double m_refXYZ[3];
	double m_tipXYZ[3];
	double m_bowXYZ[3];

	double m_rootLEXYZ[3];
	double m_tipLEXYZ[3];
	double m_bowLEXYZ[3];
	double m_rootLENOMXYZ[3];
	double m_tipLENOMXYZ[3];
	double m_bowLENOMXYZ[3];
	double m_bowLEDIR[3];

	double m_refNomCentroid[2];
	double m_refActCentroid[2];

	double m_lastXY[2];     // for ADJ STACKX and STACKY calcs
	double m_lastTXY[2];    // for ADJ STACKT calcs (in case using different fit for some reason.
	double m_lastA;
	double m_lastNomCentroid[2];
	double m_lastActCentroid[2];
	double m_lastChdAngNom;
	double m_lastChdAngAct;

	double m_rootNomXYZ[3];  // for CONVEX and CONCAVE bow calculations
	double m_tipNomXYZ[3];   // for CONVEX and CONCAVE bow calculations
	double m_bowNomXYZ[3];   // for CONVEX and CONCAVE bow calculations


	CMatrix* m_pPlat;                 // Platform Points
	CMatrix* m_pZone;                 // Zone Form Calculations Stuff
	CMatrix* m_pNX;
	CMatrix* m_pNY;
	CMatrix* m_pNI;
	CMatrix* m_pNJ;

	double** m_pZoneNew;              // for new zone form stuff, may have diff number for each section
	double** m_pZoneX;
	double** m_pZoneY;

	double m_probeRad;
	CAnalysisSect *m_sect;
	unsigned short* m_calc;
	CBlade* m_pBlade;
	/*CTolerance* m_pTol;*/
	CFlavor* m_pFlavor;
	
	int* m_pBSect;                    // map from matrix sections into blade document sections

	int** m_pBestFitSection;          // Best fit alignments
	CAnalysisCell** m_pCell;          // Matrix of calculations

	//HWND m_statusBarHWND;  // window to set status messages
	wchar_t m_statusSection[MAXBUFSZ];
	wchar_t m_statusTrimAndOrder[MAXBUFSZ];
	wchar_t m_statusCompensateAndFit[MAXBUFSZ];
	wchar_t m_statusBestFit[MAXBUFSZ];
	wchar_t m_statusCalculating[MAXBUFSZ];
	wchar_t m_statusReadingPoints[MAXBUFSZ];
	wchar_t m_processNomSection[MAXBUFSZ];
	wchar_t m_nomFilePath[MAXBUFSZ];   // needed for variable tolerance
	
	CAnalysis();
	~CAnalysis(void);

	int ReOrder(CAnalysisSect& sect, int types);
	int ReOrder(CAnalysisSect& sect, int* n1, int* n2, int types);
	int FitSplines();
	void Initialize();
	bool CalcAlign(int r, BladeBestFitType typ, int doingBow = 0, int bfind = 0, double* mtols = NULL, double* ptols = NULL);
	bool Locate(int r, double* xy, int doingBow);
	bool FillCells();
	int GetMethod(int c, int ts);
	double inchSize() const
	{
		//return m_pBlade->IsEnglish() ? 1.0 : 25.4;
		return 25.4;
	}
};

