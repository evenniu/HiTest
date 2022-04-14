#include "stdafx.h"
#include "Analysis.h"
#include "Nurb.h"
#include "SUBCURVE.H"
#include "SectionCurve.h"
#include "MeanCamberCurve.h"

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
int CAnalysis::FitSplines()
{
	int i;
	bugout(0, L"FitSplines:Processing analysis file");

	for (i = 0; i < m_pBlade->NumSect(); i++)
		m_pBlade->m_section[i]->ResetCurves();

	double t0[4], t1[4];
	CCurve* lec, * tec, * cvc, * ccc;
	m_numSuspicious = 0;
	
	for (i = 0; i < m_numSect; i++)
	{
		double voff = 0.0, uoff = 0.0;
		double mtle = -1.0, mtte = -1.0;
		double ler = -1.0, ter = -1.0;

		if (m_sect[i].m_numPoints < 10) // trail trim problem or ???
			continue;

		wchar_t analName[MAXBUFSZ];
		wcscpy_s(analName, m_sect[i].m_sectName);
		CCurve* whole = 0;
		int newPhantomIndexLE = -1;
		int newPhantomIndexTE = -1;
		whole =new CNurbCurve(m_sect[i].m_numPoints, m_sect[i].x, m_sect[i].y, 0, true, 0, 1, 1, 0.0, 0, 0, 0, 2.0, true);
		if(!whole)
		{
			return 0;
		}
		t0[0] = whole->T0();
		t1[0] = whole->T1();
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
		int noseClosest = 0, tailClosest = 0; // save indices that are closest to ends
		int lePartial = 0, tePartial = 0;
		
		if (newPhantomIndexLE > 0 && m_pBlade->m_section[i]->LEType() == EDGE_PARTIAL)
		{
			lePartial = 1;
			// so these points will be assigned to CV and CC side
			start[LEC] = newPhantomIndexLE + 1;
			end[LEC] = newPhantomIndexLE;
		}

		if (newPhantomIndexTE > 0 && m_pBlade->m_section[i]->TEType() == EDGE_PARTIAL)
		{
			tePartial = 1;
			// so these points will be assigned to CV and CC side
			start[TEC] = newPhantomIndexTE + 1;
			end[TEC] = newPhantomIndexTE;
		}

		//if (!m_pFlavor->m_fromStack && m_pTol->m_sect[ts]->m_leChange > 0.0)
		//	voff = 1.5 * m_pTol->m_sect[ts]->m_leChange; // tie these back to the tol file.

	/*	if (!m_pFlavor->m_fromStack && m_pTol->m_sect[ts]->m_teChange > 0.0)
			uoff = 1.5 * m_pTol->m_sect[ts]->m_teChange;*/

		double voff2 = -1.0;
		double uoff2 = -1.0;

		//if (m_pFlavor->m_fromStack && m_pTol->m_sect[ts]->m_leChange > 0.0)
		//	voff2 = m_pTol->m_sect[ts]->m_leChange; // tie these back to the tol file.

		//if (m_pFlavor->m_fromStack && m_pTol->m_sect[ts]->m_teChange > 0.0)
		//	uoff2 = m_pTol->m_sect[ts]->m_teChange;
		//
		double nvec[2], tvec[2];
		normalize(nvec, m_sect[i].m_nose);
		normalize(tvec, m_sect[i].m_tail);
		double origin[2] = { 0.0, 0.0 };
		for (int q = 0; q < m_sect[i].m_numPoints; q++)
		{
			if (!lePartial)
			{
				d = _hypot(m_sect[i].m_nose[0] - m_sect[i].x[q], m_sect[i].m_nose[1] - m_sect[i].y[q]);
				if (d < minNose)
				{
					minNose = d;
					noseClosest = q;
				}
			}
		}
		m_pBlade->m_section[i]->MeaCurve(whole);
		m_pBlade->m_section[i]->NomCurve(whole);

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
		double mcv1[2], mcc0[2];
		cvc->CalcPoint(mcv1, t1[CVC]);
		ccc->CalcPoint(mcc0, t0[CCC]);


		m_pBlade->m_section[i]->MeaCurve(whole);
		m_pBlade->m_section[i]->MeaPart(LEC, lec);
		m_pBlade->m_section[i]->MeaPart(TEC, tec);
		m_pBlade->m_section[i]->MeaPart(CVC, cvc);
		m_pBlade->m_section[i]->MeaPart(CCC, ccc);

		CCurve* mcc = NULL;
		if (mcc)
			delete mcc;
		Hexagon::Blade::MeanCamberCurveParameters mccParams;
		mccParams.section = m_pBlade->m_section[i];
		mccParams.analysisSection = &(m_sect[i]);
		mccParams.flavor = m_pFlavor;
		mccParams.ler = ler;
		mccParams.ter = ter;
		mccParams.uoff = uoff;
		mccParams.voff = voff;
		mccParams.mtle = mtle;
		mccParams.mtte = mtte;
		const bool isEnglish = true;//ccc->IsEnglish();
		Hexagon::Blade::MeanCamberResult meanCamberResult;
		try
		{
			meanCamberResult = Hexagon::Blade::createMeasuredMeanCamberCurve(mccParams, isEnglish);
		}
		catch (...)
		{
			std::wstring errorSection = L"measured section " + std::wstring(m_sect[i].m_sectName);
			ErrorStruct es(BE_MEANCAMBERFAILED, errorSection.c_str());
			m_error->AddError(&es);
			break;
		}
		mcc = meanCamberResult.meanCamberCurve;
		//if (!mcc->Valid())
		//{
		//	delete mcc;
		//	/* ErrorStruct es(BE_MEANCAMBERFAILED, m_sect[i].m_sectName);
		//	 m_error->AddError(&es);*/
		//	 // return 1;
		//}
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
	//if (typ == BladeBestFitType::BestFitLeastSquares)
	//{
	//	if (m_pBestFitSection[r][bfind] >= 0)
	//		return true;
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
	//	//fp.fitcurve[CVC] = m_pFlavor->m_useCV[bfind] ? 1 : 0;
	//	//fp.fitcurve[CCC] = m_pFlavor->m_useCC[bfind] ? 1 : 0;
	//	//fp.fitcurve[LEC] = m_pFlavor->m_useLE[bfind] ? 1 : 0;
	//	//fp.fitcurve[TEC] = m_pFlavor->m_useTE[bfind] ? 1 : 0;

	//	//fp.weightcurve[CVC] = m_pFlavor->m_weightCV[bfind];
	//	//fp.weightcurve[CCC] = m_pFlavor->m_weightCC[bfind];
	//	//fp.weightcurve[LEC] = m_pFlavor->m_weightLE[bfind];
	//	//fp.weightcurve[TEC] = m_pFlavor->m_weightTE[bfind];
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
	//bugout(0, L"CalcAlign: rotfit %d tranfit %d ****", fp.rotfit, fp.tranfit);

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

		}
		//if (CalcAlign(r, m_pFlavor->m_fitType[bfind], doingBow, bfind, mtols, ptols))
		//{
		//	bugout(0, L"Locate: after CalcAlign, will cal GetBestFitV1");
		//	//CBestFit* thisFit = //m_pBlade->m_section[bs]->GetBestFitV1(m_pBestFitSection[r][bfind]);
		//	//thisFit->ReportFit(m_pFlavor->m_reportFit[bfind]);
		//	if (bfind == 0)
		//	{
		//		double theta;
		//		//thisFit->ReturnFit(&xy[0], &xy[1], &theta);
		//	}
		//}
		//else if (bfind == 0)
		//{
		//	return false;
		//}

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