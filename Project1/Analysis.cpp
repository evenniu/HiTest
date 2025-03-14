#include "stdafx.h"
#include "Analysis.h"
#include "Nurb.h"
#include "SUBCURVE.H"
#include "SectionCurve.h"
#include "MeanCamberCurve.h"
#include "ArraySlicing.h"
#include "EigenAbstractCurve.h"

CAnalysis::CAnalysis(): m_error(new CBladeError())
{
	Initialize();
}
int CAnalysis::ReOrder(CAnalysisSect& sect, int types)
{
	bugout(0, L"Reorder:entered");
	sect.m_bigGap = false;
	int s;
	if (types & 4 && types & 16) // le or te partial
	{
		//ErrorStruct es(BE_CANNOTCOMPLETEBOTH);
		//m_error->AddError(&es);
		return 0;
	}
}
int CAnalysis::ReOrder(CAnalysisSect& sect, int* n1, int* n2, int types)
{
	bugout(0,L"Reorder:entered");
	sect.m_bigGap = false;
	int s;
	return 0;
}
int CAnalysis::GetMethod(int c, int ts)
{
	int method = m_pFlavor->m_method[m_calc[c]];
	//int overrideMethod = m_pTol->m_sect[ts]->m_dim[m_calc[c]].m_methodOverride;
	//if (overrideMethod != -1)
	//	method = overrideMethod;

	return method;
}
int CAnalysis::FitSplines()
{
	int i;
	bugout(0, L"FitSplines:Processing analysis file");
	if (!m_pBlade)
	{
		return 0;
	}
	for (i = 0; i < m_pBlade->NumSect(); i++)
		m_pBlade->m_section[i]->ResetCurves();

	double t0[4], t1[4];
	CCurve* lec, * tec, * cvc, * ccc;
	m_numSuspicious = 0;
	m_pBSect = new int[m_numSect];
	m_pBestFitSection = new int* [m_numSect];

	for (i = 0; i < m_numSect; i++)
	{
		m_pBestFitSection[i] = new int[MAXFITS];
		for (int jj = 0; jj < MAXFITS; jj++)
			m_pBestFitSection[i][jj] = -1;
	}
	for (i = 0; i < m_numSect; i++)//实测截面个数
	{
		m_pBSect[i] = i;
		double mtle = -1.0, mtte = -1.0;
		double ler = -1.0, ter = -1.0;

		if (m_sect[i].m_numPoints < 10) // trail trim problem or ???
			continue;


		//int ts = toleranceSectionIndex(m_pTol, m_sect[i].m_sectName);
		    
		int s = 0;
		CCurve* ncp = m_pBlade->m_section[s]->NomCurve();
		CCurve* whole = 0;
		int nomfixdat = myGetProfileInt(L"NominalRemove", 0) == 0 ? FALSE : TRUE;
		int nomtension = myGetProfileInt(L"NominalTension", 0) == 0 ? FALSE : TRUE;
		int numNompoints = m_pBlade->m_section[s]->GetNumNomPoints();
		double nomxyk[3], origmeasxy_start[2], origmeasxy_end[2], Min = 1000;
		origmeasxy_start[0] = m_sect[i].x[0];
		origmeasxy_start[1] = m_sect[i].y[0];
		origmeasxy_end[0] = m_sect[i].x[m_sect[i].m_numPoints - 1];
		origmeasxy_end[1] = m_sect[i].y[m_sect[i].m_numPoints - 1];
		int nomindex1 = 0, nomindex2 = 0;
		double nt_start, nt_end;
		double tmp_t0, tmp_t1, xy[2];
		tmp_t0 = m_pBlade->m_section[s]->NomCurve()->T0();
		tmp_t1 = m_pBlade->m_section[s]->NomCurve()->T1();
		m_pBlade->m_section[s]->NomCurve()->ClosestPoint(origmeasxy_start, xy, &nt_start, NULL, tmp_t0, tmp_t1);
		m_pBlade->m_section[s]->NomCurve()->ClosestPoint(origmeasxy_end, xy, &nt_end, NULL, tmp_t0, tmp_t1);
		if (nt_start > nt_end) // do swap
		{
			int measnum = m_sect[i].m_numPoints;
			double* mx;
			double* my;
			double* mz;
			double* mi;
			double* mj;
			mx = new double[measnum];
			my = new double[measnum];
			mz = new double[measnum];
			mi = new double[measnum];
			mj = new double[measnum];
			for (int j = 0; j < measnum; j++)
			{
				mx[j] = m_sect[i].x[j];
				my[j] = m_sect[i].y[j];
				mz[j] = m_sect[i].z[j];
				mi[j] = m_sect[i].i[j];
				mj[j] = m_sect[i].j[j];
			}
			for (int j = 0; j < measnum; j++)
			{
				int k = (measnum - 1 - j) % measnum;
				m_sect[i].x[k] = mx[j];
				m_sect[i].y[k] = my[j];
				m_sect[i].z[k] = mz[j];
				m_sect[i].i[k] = mi[j];
				m_sect[i].j[k] = mj[j];
			}
			delete[] mx;
			delete[] my;
			delete[] mz;
			delete[] mi;
			delete[] mj;
		}
		whole =new CNurbCurve(m_sect[i].m_numPoints, m_sect[i].x, m_sect[i].y, 0, true, 0, 1, 1, 0.0, 0, 0, 0, 2.0, true);
		if(!whole)
		{
			return 0;
		}
		t0[CVC] = whole->T0();
		t1[CVC] = whole->T1();
		double xxx[2];

		double period = whole->T1() - whole->T0();
		whole->CalcPoint(xxx, whole->T0());
		bugout(0, L"FinalSplines (%f ,%f) T0(%f)", xxx[0], xxx[1], whole->T0());
		whole->CalcPoint(xxx, period / 2.0);
		bugout(0, L"FinalSplines (%f ,%f) P/2(%f)", xxx[0], xxx[1], period / 2.0);
		whole->CalcPoint(xxx, period);
		bugout(0, L"FinalSplines (%f ,%f)  T1(%f)", xxx[0], xxx[1], whole->T1());

		if (!whole)
		{
			ErrorStruct es(BE_COMPENSATIONFAILED, m_sect[i].m_sectName);
			m_error->AddError(&es);
			break;
		}
		period = whole->T1() - whole->T0();
		bugout(0, L"FinalSplineFit: whole Meascurve period(%lf)", period);
		// associate points with curve components

		double d, minNose = 1.0e20, minTail = 1.0e20;
		int start[4], end[4];
		double np[2],dummy[2], tp[2], nbt;
		np[0] = 0;
		np[1] = 1;
		double dddd = whole->ClosestPoint(np, dummy, &nbt, tp, 0.0, 0.0, 400);
		start[LEC] = start[TEC] = -1;
		end[LEC] = end[TEC] = -1;

		m_pBlade->m_section[i]->MeaCurve(whole);

		t0[LEC] = whole->T0();
		t1[LEC] = whole->T1();
		lec = new CSubCurve(whole, t0[LEC], t1[LEC], period);
		lec->Extreme(t0[LEC]);
		t0[TEC] = period / 2;
		t1[TEC] = whole->T1();
		t0[CCC] = period / 2 + period/4;
		t1[CCC] = period / 2 + period/2;

		tec = new CSubCurve(whole, t0[TEC], t1[TEC], period / 2);
		tec->Extreme(t1[LEC]);

		cvc = new CSubCurve(whole, t0[CVC], t1[CVC], period);
		ccc = new CSubCurve(whole, t0[CCC], t1[CCC], period);


		m_pBlade->m_section[i]->MeaCurve(whole);
		m_pBlade->m_section[i]->MeaPart(LEC, lec);
		m_pBlade->m_section[i]->MeaPart(TEC, tec);
		m_pBlade->m_section[i]->MeaPart(CVC, cvc);
		m_pBlade->m_section[i]->MeaPart(CCC, ccc);

		CCurve* mcc = NULL;
		if (mcc)
			delete mcc;
		

		m_pBlade->m_section[i]->AssignPoints(m_sect[i].x, m_sect[i].y, m_sect[i].m_numPoints, start, end);
		


	}
	if (i < m_numSect)
		return 0;

	
	return 1;
}
void CAnalysis::Initialize()
{
	m_numSect = 0;
	m_numCalc = 0;
	m_numPlat = 0;
	m_numTraces = 0;
	//m_calc = NULL;
	m_tnames = NULL;
	m_traces = NULL;
	m_sectionNames = NULL;
	m_calcLabels = NULL;
	m_probeRad = -1.0;
	m_refchdangchecked = false;
	m_refChecked = false;
	m_rootChecked = false;
	m_tipChecked = false;
	m_refSect = -1;
	m_rootSect = -1;
	m_tipSect = -1;

	m_refRow = -1;
	m_rootRow = -1;
	m_tipRow = -1;

	m_decAng = 1;
	m_decMea = 4;

	m_pFlavor = NULL;
	//m_pTol = NULL;
	//m_pPlat = NULL;
	m_pBlade = NULL;
	m_sect = NULL;
	m_pBSect = NULL;
	m_pZone = NULL;
	m_pZoneNew = NULL;
	m_pZoneX = NULL;
	m_pZoneY = NULL;
	m_pNX = NULL;
	m_pNY = NULL;
	m_pNI = NULL;
	m_pNJ = NULL;
	m_pCell = NULL;
	m_pBestFitSection = NULL;

	m_lotID = -1;
	m_lotSeq = -1;
	m_lotTotal = -1;
	m_lotSize = -1;
	m_numLotTrans = 0;
	m_refChecked = m_rootChecked = m_tipChecked = m_bowChecked = false;
	m_refchdangchecked = m_rootLEChecked = m_tipLEChecked = m_bowLEChecked = false;
	m_bowCorrect = m_twistCorrect = m_dispCorrect = 0.0;
	m_refchdangnom = m_refchdangact = 0.0;
	m_extraTolerance = 0.0;
}


// as of October 2018, the following functions are defined in SECTION.CPP
CBestFit* createMeasuredPointsToNominalCurveBestFit(const Hexagon::Blade::SectionCurve& nominalCurves,
	const Eigen::Isometry2d& measuredToNominalTransform,
	const Eigen::Ref<const Eigen::Matrix2Xd>& measuredPoints, const Eigen::Ref<const Eigen::ArrayXb>& isUsedInFit,const CFitParams& fitParams, double* mtols = nullptr, double* ptols = nullptr);
Eigen::Matrix2Xd constructPointMatrix(const double* xValues, const double* yValues, const ptrdiff_t n);

/// <summary>
/// 
/// </summary>
/// <param name="r">section index</param>
/// <param name="typ"></param>
/// <param name="doingBow"></param>
/// <param name="bfind"></param>
/// <param name="mtols"></param>
/// <param name="ptols"></param>
/// <returns></returns>
bool CAnalysis::CalcAlign(int r, BladeBestFitType typ, int doingBow, int bfind, double* mtols, double* ptols)
{
	bugout(0, L"CalcAlign(): enterdd");
	CFitParams fp;
	fp.usenominals = 1;// m_pFlavor->m_usenominals[bfind];
	fp.weightcurve[CVC] = 1;
	fp.weightcurve[CCC] = 1;
	fp.weightcurve[LEC] = 1;
	fp.weightcurve[TEC] = 1;
	fp.rotfit = 0;
	int bs = m_pBSect[r];
	bs = 0;
	int ts = 0;// toleranceSectionIndex(m_pTol, m_pBlade->m_section[bs]->Name()); // index into m_pTol->Sect;
	if (ts < 0)
		return false;
	if (typ == BladeBestFitType::BestFitFullBlade ||
		typ == BladeBestFitType::BestFitFullBladeLELS) // fitting all sections at once.
	{
		if (m_pBestFitSection[r][bfind] >= 0) // calculations already done.
			return true;
		int sec1 = -1, sec2 = -1;
		if (m_pFlavor->m_rootTip[bfind] && m_rootChecked && m_tipChecked && m_rootSect != StackOriginPointGhostSection)
		{
			sec1 = m_tipSect;
			sec2 = m_rootSect;
		}
		CBestFit* wholeFit = NULL;
		if (typ == BladeBestFitType::BestFitFullBlade)
		{
			//if (!m_pBlade->FitBlade(&fp, sec1, sec2)) // Perform the fit
			//	return false;

			//wholeFit = m_pBlade->GetBestFit();
		}
		else
		{
			fp.fitcurve[CVC] = m_pFlavor->m_useCV[bfind] ? 1 : 0;
			fp.fitcurve[CCC] = m_pFlavor->m_useCC[bfind] ? 1 : 0;
			fp.fitcurve[LEC] = 0;
			fp.fitcurve[TEC] = 0;

			//m_pTol->m_defaultSection.GetLERadii(&fp.leoff1, &fp.leoff2);
			//if (!m_pBlade->FitBladeLELS(&fp, this)) // Perform the fit
			//	return false;

			//wholeFit = m_pBlade->GetBestFitLE();
		}

		// make a copy of the best fit.
		fp.algorithm = BestFitAlgorithm::None; // no fit
		fp.fitcurve[CVC] = 1;
		fp.fitcurve[CCC] = 1;
		fp.fitcurve[LEC] = 1;
		fp.fitcurve[TEC] = 1;
		fp.weightcurve[CVC] = 1;
		fp.weightcurve[CCC] = 1;
		fp.weightcurve[LEC] = 1;
		fp.weightcurve[TEC] = 1;
		const Eigen::Matrix2Xd measuredPoints = constructPointMatrix(
			m_pBlade->m_section[bs]->m_mxpt, m_pBlade->m_section[bs]->m_mypt, m_pBlade->m_section[bs]->m_totalPoints);
		const Eigen::Map<const Eigen::ArrayXi> partOf(m_pBlade->m_section[bs]->m_partOf, m_pBlade->m_section[bs]->m_totalPoints);
		m_pBlade->m_section[bs]->m_numBestFits++;
		return true;
	}
	if (typ == BladeBestFitType::BestFitNone) // no fit
	{
		

		fp.algorithm = BestFitAlgorithm::None; // no fit
		fp.fitcurve[CVC] = 1;
		fp.fitcurve[CCC] = 1;
		fp.fitcurve[LEC] = 1;
		fp.fitcurve[TEC] = 1;

		if (m_pBlade->m_section[bs]->FitPoints(fp, m_pBestFitSection[r][bfind], inchSize(), mtols, ptols))
			return true;

		return false;
	}
/*createMeasuredPointsToNominalCurveBestFit(Hexagon::Blade::nominalSectionCurve(m_pBlade->m_section[bs]),
			Hexagon::Blade::toIsometry2d(*wholeFit->GetAlign()), measuredPoints,
			Eigen::ArrayXb::Ones(partOf.size()), fp, mtols, ptols)*/;


		fp.algorithm = BestFitAlgorithm::LeastSquares; // ls fit
		fp.fitcurve[CVC] = 1;
		fp.fitcurve[CCC] = 1;
		fp.fitcurve[LEC] = 1;
		fp.fitcurve[TEC] = 1;
		fp.weightcurve[CVC] = 1;
		fp.weightcurve[CCC] = 1;
		fp.weightcurve[LEC] = 1;
		fp.weightcurve[TEC] = 1;
		fp.lepercent = 5.0;
		fp.tepercent = 95.0;
		fp.fitcurve[CVC] = m_pFlavor->m_useCV[bfind] ? 1 : 0;
		fp.fitcurve[CCC] = m_pFlavor->m_useCC[bfind] ? 1 : 0;
		fp.fitcurve[LEC] = m_pFlavor->m_useLE[bfind] ? 1 : 0;
		fp.fitcurve[TEC] = m_pFlavor->m_useTE[bfind] ? 1 : 0;

		fp.weightcurve[CVC] = m_pFlavor->m_weightCV[bfind];
		fp.weightcurve[CCC] = m_pFlavor->m_weightCC[bfind];
		fp.weightcurve[LEC] = m_pFlavor->m_weightLE[bfind];
		fp.weightcurve[TEC] = m_pFlavor->m_weightTE[bfind];
		if (m_pFlavor->m_noTranslate[bfind])
			fp.tranfit = 1; // no translation
		else
			fp.tranfit = 0; // full translation allowed

		if (m_pFlavor->m_noRotate[bfind])
			fp.rotfit = 1; // no rotation
		else
			fp.rotfit = 0; // full rotation allowed

		if (fp.fitcurve[LEC] && !fp.fitcurve[CVC] && !fp.fitcurve[CCC] && !fp.fitcurve[TEC])
			fp.pivot = 1;
		else if (fp.fitcurve[TEC] && !fp.fitcurve[CVC] && !fp.fitcurve[CCC] && !fp.fitcurve[LEC])
			fp.pivot = 3;
		switch (m_pFlavor->m_Transfit_bf[bfind])
		{
		case 0:
			fp.tranfit = 0;
			break;
		case 1:
			fp.tranfit = 1;
			break;
		case 2:
			fp.tranfit = 4;
			break;
		case 3:
			fp.tranfit = 5;
			break;
		}
	//}
	try
	{
		//if (m_pBlade->m_section[bs]->FitPoint(fp, m_pBestFitSection[r][bfind], inchSize(), mtols, ptols))
		int tmp_index = m_pBestFitSection[r][bfind];
		if (m_pBlade->m_section[bs]->FitPoints(fp,tmp_index, inchSize(), mtols, ptols))
			return true;
	}
	catch (...)
	{
		return false;
	}
	bugout(0, L"CalcAlign: rotfit %d tranfit %d ****", fp.rotfit, fp.tranfit);

	return false;
}

/// <summary>
/// 
/// </summary>
/// <param name="r">section index</param>
/// <param name="xy"></param>
/// <param name="doingBow"></param>
/// <returns></returns>
bool CAnalysis::Locate(int r, double* xy, int doingBow)
{
	bugout(0, L"Locate(): enterdd");
	int bs = m_pBSect[r];
	int ts = 0;// toleranceSectionIndex(m_pTol, m_pBlade->m_section[bs]->Name()); // index into m_pTol->Sect;

	for (int bfind = 0; bfind < MAXFITS; bfind++)
	{
		bool thisFitUsed = bfind == 0 ? true : false;

		double mtols[4] = { -100, -100, -100, -100 };
		double ptols[4] = { 100, 100, 100, 100 };
		for (int c = 0; c < m_numCalc; c++)
		{
			if (m_pFlavor->m_fitToUse[m_calc[c]] == bfind)
			{
				thisFitUsed = true;

				if (ts >= 0)
				{
					int toltype = 0;//m_pTol->m_sect[ts]->m_dim[m_calc[c]].m_type;
					if (toltype >= 0)
					{
						double mtol = -0.05;// m_pTol->m_sect[ts]->m_dim[m_calc[c]].m_mtol;
						double ptol = 0.05;//m_pTol->m_sect[ts]->m_dim[m_calc[c]].m_ptol;

						switch (m_calc[c])
						{
						case CalcLEContour:
						case CalcLEContour2:
							mtols[LEC] = -0.5 * ptol;
							ptols[LEC] = 0.5 * ptol;
							break;
						}
					}
				}
			}
		}
		if (!thisFitUsed) // not used;
			continue;
		if (CalcAlign(r, m_pFlavor->m_fitType[bfind], doingBow, bfind, mtols, ptols))
		{
			bugout(0, L"Locate: after CalcAlign, will cal GetBestFitV1");
			//* thisFit = m_pBlade->m_section[bs]->GetBestFitV1(m_pBestFitSection[r][bfind]);
		}

	}
	return false;
}

bool CAnalysis::Locate(int secid, BladeBestFitType fitType,int fitToMiddleOfZone, int Transfit, bool noRotate, int rotfit, int useNominal)
{
	bugout(0, L"Locate: entered");
	double mtols[4] = { -100, -100, -100, -100 };
	double ptols[4] = { 100, 100, 100, 100 };
	int bfind = 0;
	m_pFlavor->m_Transfit_bf[bfind] = Transfit;
	m_pFlavor->m_noTranslate[bfind] = noRotate;
	m_pFlavor->m_noRotate[bfind] = noRotate;
	m_pFlavor->m_usenominals[bfind] = useNominal;
	m_pFlavor->m_fitToMiddleOfZone[bfind] = fitToMiddleOfZone;

	if (CalcAlign(secid, fitType, 0, bfind, mtols, ptols))
	{
		bugout(0, L"Locate: after CalcAlign, will cal GetBestFitV1");
	#if 1//输出打印调试信息
			CBestFit* thisFit = m_pBlade->m_section[secid]->GetBestFitV1(m_pBestFitSection[secid][bfind]);
			thisFit->ReportFit(m_pFlavor->m_reportFit[bfind]);

			double x, y, ang, np[2], bp[2];
			thisFit->ReturnFit(&x, &y, &ang);
			bugout(0, L"Locate: after CalcAlign %d, will cal GetBestFitV1{%lf,%lf, %lf}",
				m_pBestFitSection[secid][bfind], x, y, ang);
			int numPoints = thisFit->NumPoints();
			//for(int i = 0; i < numPoints; i++)
			//{
			//  np[0] = thisFit->m_noms->m[i][0];
			//  np[1] = thisFit->m_noms->m[i][1];
			//  bp[0] = thisFit->m_infs->m[i][0];
			//  bp[1] = thisFit->m_infs->m[i][1];
			//  bugout(3, L"after CalcAlign  %lf %lf %lf %lf}", bp[0], bp[1], np[0], np[1]);
			//}
	#endif 
	}
	return false;
}

bool CAnalysis::FillCells()
{
	bugout(0,L"Fillcells(): enterd ***");
	m_refNomCentroid[0] = m_refNomCentroid[1] = m_refActCentroid[0] = m_refActCentroid[1] = -1.0e20;
	m_good = true;
	int i, j, bs;
	int doingBow = 0;

	if (m_numSect > 0)
	{
		m_pBSect = new int[m_numSect]; // indices in m_pBlade->m_section[]

		if (m_pFlavor->m_specials[SpecialZoneForm] == 1 || m_pFlavor->m_specials[SpecialZoneForm] == 3)
		{
			m_pZone = new CMatrix(m_numSect, 30);
			m_pNX = new CMatrix(m_numSect, 30);
			m_pNY = new CMatrix(m_numSect, 30);
			m_pNI = new CMatrix(m_numSect, 30);
			m_pNJ = new CMatrix(m_numSect, 30);
		}
		else if (m_pFlavor->m_specials[SpecialZoneForm] > 1)
		{
			m_pZoneNew = new double* [m_numSect];
			m_pZoneX = new double* [m_numSect];
			m_pZoneY = new double* [m_numSect];
		}

		if (!m_pZoneNew && m_pFlavor->m_specials[SpecialSaveXYZFile] && m_pFlavor->m_xyzFormat == 1)
		{
			m_pZoneNew = new double* [m_numSect];
			m_pZoneX = new double* [m_numSect];
			m_pZoneY = new double* [m_numSect];
		}
		// will need to fill m_pZoneNew if form calcs use nom file method

		for (i = 0; i < m_numCalc && !m_pZoneNew; i++)
		{
			if (m_calc[i] == CalcMinForm || m_calc[i] == CalcMaxForm)
			{
				int method = GetMethod(i, 0);
				if (method == MethodFormVariableNomFile || method == MethodFormVariableVarFile)
				{
					m_pZoneNew = new double* [m_numSect];
					m_pZoneX = new double* [m_numSect];
					m_pZoneY = new double* [m_numSect];
					break;
				}
			}
		}
		m_pBestFitSection = new int* [m_numSect];
		for (i = 0; i < m_numSect; i++)
		{
			m_pBestFitSection[i] = new int[MAXFITS];
			for (int jj = 0; jj < MAXFITS; jj++)
				m_pBestFitSection[i][jj] = -1;
		}
		for (int ss = 0; ss < m_numSect; ss++)
		{
			double origin[2], axis[2], xy[2], oxy[2], nxy[2];
			if (Locate(ss, xy, doingBow))
			{
				for (int i = 0; i < m_numCalc; i++)
				{
					//Calculate(j, i, bs, ts);
				}
			}
		}
	}
	return false;
}

CAnalysis::~CAnalysis(void)
{
	if (m_error)
	{
		delete m_error;
		m_error = nullptr;
	}
	if (m_pPlat)
		delete m_pPlat;

	if (m_pBlade)
		delete m_pBlade;

	if (m_sect)
		delete[] m_sect;

	if (m_pBSect)
		delete[] m_pBSect;

	if (m_pZone)
		delete m_pZone;

	if (m_pZoneNew)
	{
		for (int s = 0; s < m_numSect; s++)
			delete[] m_pZoneNew[s];
		delete[] m_pZoneNew;
	}
	if (m_pZoneX)
	{
		for (int s = 0; s < m_numSect; s++)
			delete[] m_pZoneX[s];
		delete[] m_pZoneX;
	}
	if (m_pZoneY)
	{
		for (int s = 0; s < m_numSect; s++)
			delete[] m_pZoneY[s];
		delete[] m_pZoneY;
	}

	if (m_pNX)
		delete m_pNX;
	if (m_pNY)
		delete m_pNY;
	if (m_pNI)
		delete m_pNI;
	if (m_pNJ)
		delete m_pNJ;

	//if (m_pCell)
	//{
	//	for (int s = 0; s < m_numSect; s++)
	//		delete[] m_pCell[s];

	//	delete[] m_pCell;
	//}

	if (m_pBestFitSection)
	{
		for (int s = 0; s < m_numSect; s++)
			delete m_pBestFitSection[s];
		delete[] m_pBestFitSection;
	}

	if (m_numTraces > 0 && m_tnames != NULL && m_traces != NULL)
	{
		int i;
		for (i = 0; i < m_numTraces; i++)
		{
			delete[] m_tnames[i];
			delete[] m_traces[i];
		}
		delete[] m_tnames;
		delete[] m_traces;
	}

	if (m_calc != NULL)
		delete[] m_calc;

	if (m_sectionNames)
	{
		for (int s = 0; s < m_numSect; s++)
			delete[] m_sectionNames[s];
		delete[] m_sectionNames;
	}

	if (m_calcLabels)
	{
		for (int i = 0; i < m_numCalc; i++)
			delete[] m_calcLabels[i];
		delete[] m_calcLabels;
	}
}

CAnalysisSect::CAnalysisSect()
{
	m_goodVectors = true;
	m_bigGap = false;
	m_numPoints = 0;
	numberOfBallCenters = 0;
	m_phantomIndexLE = -1;
	m_phantomIndexTE = -1;
	//skewalign = 0; // don't delete in destructor, owned by Section
	x = 0;
	y = 0;
	z = 0;
	i = 0;
	j = 0;
	ox = 0;
	oy = 0;
	oz = 0;
	oi = 0;
	oj = 0;
	ballCenterX = nullptr;
	ballCenterY = nullptr;
	ballCenterZ = nullptr;
}

CAnalysisSect::~CAnalysisSect()
{
	if (x)
		delete[] x;
	if (y)
		delete[] y;
	if (z)
		delete[] z;
	if (i)
		delete[] i;
	if (j)
		delete[] j;
	if (ox)
		delete[] ox;
	if (oy)
		delete[] oy;
	if (oz)
		delete[] oz;
	if (oi)
		delete[] oi;
	if (oj)
		delete[] oj;
	if (ballCenterX)
		delete[] ballCenterX;
	if (ballCenterY)
		delete[] ballCenterY;
	if (ballCenterZ)
		delete[] ballCenterZ;
}

CAnalysisSect::CAnalysisSect(const CAnalysisSect& obj)
{
	m_numPoints = obj.m_numPoints;
	numberOfBallCenters = obj.numberOfBallCenters;
	m_inose = obj.m_inose;
	m_itail = obj.m_itail;
	m_inter1 = obj.m_inter1;
	m_inter2 = obj.m_inter2;
	m_phantomIndexLE = obj.m_phantomIndexLE;
	m_phantomIndexTE = obj.m_phantomIndexTE;
	wcscpy_s(m_sectName, obj.m_sectName);
	m_nose[0] = obj.m_nose[0];
	m_nose[1] = obj.m_nose[1];
	m_tail[0] = obj.m_tail[0];
	m_tail[1] = obj.m_tail[1];
	if (m_numPoints < 1)
	{
		m_numPoints = 0;
		numberOfBallCenters = 0;
		x = NULL;
		y = NULL;
		z = NULL;
		i = NULL;
		j = NULL;
		ox = NULL;
		oy = NULL;
		oz = NULL;
		oi = NULL;
		oj = NULL;
		ballCenterX = nullptr;
		ballCenterY = nullptr;
		ballCenterZ = nullptr;
	}
	else
	{
		x = new double[m_numPoints];
		y = new double[m_numPoints];
		z = new double[m_numPoints];
		i = new double[m_numPoints];
		j = new double[m_numPoints];
		ox = new double[m_numPoints];
		oy = new double[m_numPoints];
		oz = new double[m_numPoints];
		oi = new double[m_numPoints];
		oj = new double[m_numPoints];
		ballCenterX = new double[numberOfBallCenters];
		ballCenterY = new double[numberOfBallCenters];
		ballCenterZ = new double[numberOfBallCenters];

		for (int ii = 0; ii < m_numPoints; ii++)
		{
			x[ii] = obj.x[ii];
			y[ii] = obj.y[ii];
			z[ii] = obj.z[ii];
			i[ii] = obj.i[ii];
			j[ii] = obj.j[ii];
			ox[ii] = obj.ox[ii];
			oy[ii] = obj.oy[ii];
			oz[ii] = obj.oz[ii];
			oi[ii] = obj.oi[ii];
			oj[ii] = obj.oj[ii];
		}
		for (size_t ii = 0; ii < numberOfBallCenters; ii++)
		{
			ballCenterX[ii] = obj.ballCenterX[ii];
			ballCenterY[ii] = obj.ballCenterY[ii];
			ballCenterZ[ii] = obj.ballCenterZ[ii];
		}
	}
}

CAnalysisSect& CAnalysisSect::operator=(const CAnalysisSect& obj) // assignment operator
{
	if (m_numPoints)
	{
		if (x)
			delete[] x;
		if (y)
			delete[] y;
		if (z)
			delete[] z;
		if (i)
			delete[] i;
		if (j)
			delete[] j;
		if (ox)
			delete[] ox;
		if (oy)
			delete[] oy;
		if (oz)
			delete[] oz;
		if (oi)
			delete[] oi;
		if (oj)
			delete[] oj;
		if (ballCenterX)
			delete[] ballCenterX;
		if (ballCenterY)
			delete[] ballCenterY;
		if (ballCenterZ)
			delete[] ballCenterZ;
	}

	m_numPoints = obj.m_numPoints;
	numberOfBallCenters = obj.numberOfBallCenters;
	m_inose = obj.m_inose;
	m_itail = obj.m_itail;
	m_inter1 = obj.m_inter1;
	m_inter2 = obj.m_inter2;
	m_phantomIndexLE = obj.m_phantomIndexLE;
	m_phantomIndexTE = obj.m_phantomIndexTE;
	wcscpy_s(m_sectName, obj.m_sectName);
	m_nose[0] = obj.m_nose[0];
	m_nose[1] = obj.m_nose[1];
	m_tail[0] = obj.m_tail[0];
	m_tail[1] = obj.m_tail[1];
	if (m_numPoints < 1)
	{
		m_numPoints = 0;
		x = NULL;
		y = NULL;
		z = NULL;
		i = NULL;
		j = NULL;
		ox = NULL;
		oy = NULL;
		oz = NULL;
		oi = NULL;
		oj = NULL;
		ballCenterX = nullptr;
		ballCenterY = nullptr;
		ballCenterZ = nullptr;
	}
	else
	{
		x = new double[m_numPoints];
		y = new double[m_numPoints];
		z = new double[m_numPoints];
		i = new double[m_numPoints];
		j = new double[m_numPoints];
		ox = new double[m_numPoints];
		oy = new double[m_numPoints];
		oz = new double[m_numPoints];
		oi = new double[m_numPoints];
		oj = new double[m_numPoints];
		ballCenterX = new double[numberOfBallCenters];
		ballCenterY = new double[numberOfBallCenters];
		ballCenterZ = new double[numberOfBallCenters];

		for (int ii = 0; ii < m_numPoints; ii++)
		{
			x[ii] = obj.x[ii];
			y[ii] = obj.y[ii];
			z[ii] = obj.z[ii];
			i[ii] = obj.i[ii];
			j[ii] = obj.j[ii];
			ox[ii] = obj.ox[ii];
			oy[ii] = obj.oy[ii];
			oz[ii] = obj.oz[ii];
			oi[ii] = obj.oi[ii];
			oj[ii] = obj.oj[ii];
		}
		for (size_t ii = 0; ii < numberOfBallCenters; ii++)
		{
			ballCenterX[ii] = obj.ballCenterX[ii];
			ballCenterY[ii] = obj.ballCenterY[ii];
			ballCenterZ[ii] = obj.ballCenterZ[ii];
		}
	}

	return *this;
}

CAnalysisCell::CAnalysisCell()
{
	m_label[0] = 0;
	m_show = 0;
	m_decimals = 4;
	m_nom = 0.0;
	m_act = 0.0;
	m_ltol = 0.0;
	m_utol = 0.0;
	m_outtol = 0.0;
	m_box.bottom = m_box.left = m_box.right = m_box.top = 0;
}