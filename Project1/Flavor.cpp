#include "Flavor.h"

CFlavor::CFlavor(void)
{
    m_numCalc = 0;
    m_summary = 0;

    m_DPDir[0] = 0;
    m_DPPlusDir[0] = 0;
    m_DPPlusDir1[0] = 0;
    m_outputRegDimDir[0] = 0;
    m_headerFile[0] = 0;
    m_comment[0] = 0;
    m_dstFile[0] = 0;
    m_exeFile[0] = 0;

    int i;
    for (i = 0; i < MAXFITS; i++)
    {
        m_fitType[i] = i == 0 ? BladeBestFitType::BestFitLeastSquares : BladeBestFitType::BestFitNone;
        m_usenominals[i] = 0;
        m_reportFit[i] = 0;
        m_percent[i][0] = 5;
        m_percent[i][1] = 95;
        m_rootTip[i] = FALSE;

        m_useCV[i] = TRUE;
        m_useCC[i] = TRUE;
        m_useLE[i] = TRUE;
        m_useTE[i] = TRUE;

        m_noTranslate[i] = FALSE;
        m_noRotate[i] = FALSE;
        m_Transfit_bf[i] = 0;//FULL
        m_weightCV[i] = 1;
        m_weightCC[i] = 1;
        m_weightLE[i] = 1;
        m_weightTE[i] = 1;

        m_fitToMiddleOfZone[i] = FALSE;
        m_complexEdgeZoneIndex[i] = -1;
        m_chordZoneIndex[i] = -1;
    }

    for (i = 0; i < CalcArraySize; i++)
    {
        m_calcs[i] = CalcUnset;
        m_method[i] = 0;
        m_mult[i] = 1.0;
        m_offset[i] = 0.0;
        m_fitToUse[i] = 0;
        m_head1[i][0] = 0;
        m_head2[i][0] = 0;
        m_RoadTHCKRotate[i] = 0.0;
        m_RoadTHCKoffset[i] = 0;
        m_RoadTHCKRotateXaxis[i] = 0.0; //初始化
        m_useExtreamABS[i] = false;//初始化
    }

    for (i = 0; i < SpecialArraySize; i++)
        m_specials[i] = 0;

    m_specials[SpecialAutoPrint] = 1;
    m_specials[SpecialPlotEdges] = 1;

    m_endMag = 1.0;
    m_sideMag = 1.0;
    m_printScale = -1.0;


    m_fromStack = FALSE;
    m_xyzBestFit = FALSE;
    m_xyzFormat = 0;
    m_avgFindBestFit = TRUE;
    m_avgOutBestFit = FALSE;
    m_calcCMM = FALSE;

    m_stockMinMax = FALSE;
    m_stockPlot = 0;

    m_compensateMethod = 0;
    m_filterSpacing = 3;
    m_resortEnds = 1;
    m_removeSquare = 0;
    m_notTrimTail = 0;
    m_regDime = FALSE;
    // default these to current config
}
CFlavor::~CFlavor(void)
{
}
bool CFlavor::SetNumCalcs(int numCalcs)
{
	if (numCalcs < 0)
		return false;

	m_numCalc = numCalcs;

	return true;
}