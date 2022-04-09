#pragma once
#include "stdafx.h"

class  CFlavor
{
public:
	CFlavor(void);
	~CFlavor(void);
	bool SetNumCalcs(int numCalcs);
	BladeBestFitType m_fitType[MAXFITS];      // best fit method
	int m_percent[MAXFITS][2];                // for two point fits
	int m_numCalc;                            // if calcs are hard coded, the number of calcs.
	BladeCalculations m_calcs[CalcArraySize]; // if calcs are hard coded, the calcs.
	int m_method[CalcArraySize];              // defines which method to use for some calcs.
	BOOL m_useExtreamABS[CalcArraySize];      //Extream是否取绝对值
	int m_specials[SpecialArraySize];         // array of flags to implement special options.
	int m_usenominals[MAXFITS];               // use nominal points rather than measured.

	wchar_t m_DPDir[MAXBUFSZ];                // DataPage Directory
	wchar_t m_DPPlusDir[MAXBUFSZ];            // DataPage Directory
	wchar_t m_DPPlusDir1[MAXBUFSZ];            // QDAS Directory
	wchar_t m_outputRegDimDir[MAXBUFSZ];     // 常规尺寸json Directory
	wchar_t m_headerFile[MAXBUFSZ];           // Header File
	wchar_t m_comment[MAXBUFSZ];              // Note
	wchar_t m_dstFile[MAXBUFSZ];              // Point Distribution file
	wchar_t m_exeFile[MAXBUFSZ];              // Name of an optional executable file to execute after all output files have been created

	wchar_t m_head1[CalcArraySize][MAXBUFSZ]; // user defined calculation column headers.
	wchar_t m_head2[CalcArraySize][MAXBUFSZ]; // user defined calculation column headers.
	double m_mult[CalcArraySize];             // multipliers - usually 1.0 or -1.0
	double m_RoadTHCKRotate[CalcArraySize];           // 旋转角度 
	double m_RoadTHCKRotateXaxis[CalcArraySize];      // 绕X轴旋转角度 -
	double m_RoadTHCKoffset[CalcArraySize];           // 偏移距离 offsets -
	double m_offset[CalcArraySize];           // offsets - usually 0.0
	int m_fitToUse[CalcArraySize];            // which best fit should be applied, if applicable?
	double m_endMag;                          // multiplier for end form plot deviations
	double m_sideMag;                         // multiplier for side form plot deviations
	double m_printScale;                      // scale factor for printing form plots
	UINT m_summary;                           // bit map of summarizations of each caclulation;
	double m_forgeAngle[MAXFITS];             // used for two point forge fit

	BOOL m_rootTip[MAXFITS];                  // use only root and tip sections in full blade fitting
	BOOL m_useCV[MAXFITS];                    // use CV if this is a full blade fitting
	BOOL m_useCC[MAXFITS];                    // use CC if this is a full blade fitting
	BOOL m_useLE[MAXFITS];                    // use LE if this is a full blade fitting
	BOOL m_useTE[MAXFITS];                    // use TE if this is a full blade fitting
	BOOL m_noTranslate[MAXFITS];              // translation turned off for lsq fits
	BOOL m_noRotate[MAXFITS];                 // rotation turned off for lsq fits
	int  m_Transfit_bf[MAXFITS];
	UINT m_weightCV[MAXFITS];                 // weight of CV if this is a full blade fitting
	UINT m_weightCC[MAXFITS];                 // weight of CC if this is a full blade fitting
	UINT m_weightLE[MAXFITS];                 // weight of LE if this is a full blade fitting
	UINT m_weightTE[MAXFITS];                 // weight of TE if this is a full blade fitting
	int  m_reportFit[MAXFITS];                // report versus pivot point or stack point
	int  m_fitToMiddleOfZone[MAXFITS];        // fit to the middle of the tolerance zone, or fit to nominal?
	int m_complexEdgeZoneIndex[MAXFITS];
	int m_chordZoneIndex[MAXFITS];
	BOOL m_fromStack;                         // measure offsets from stack point
	BOOL m_xyzBestFit;                        // XYZ reports use bestfit alignment
	BOOL m_avgFindBestFit;                    // AVG reports to used best fit to find points
	BOOL m_avgOutBestFit;                     // AVG reports to output points in best fit coordinates
	BOOL m_calcCMM;                           // CMM reports organized by calculation rather than section
	UINT m_xyzFormat;                         // XYZ file format; 0=PC-DMIS, 1=TOLERANCED, 2=DEVIATION...
	BOOL m_stockMinMax;                       // whether to show min and max stock points on form plots
	int  m_stockPlot;                         // what to show on form plots

	// values that are also in overall setup

	int  m_compensateMethod;                  // 0 means compute normals, 1 means use normals from scans
	UINT m_filterSpacing;                     // use every nth point when building ball center curve to start compensation
	BOOL m_resortEnds;                        // resort the end points after probe compensation
	BOOL m_removeSquare;                      // for square ends with sharp edges, remove ball center data at corners before compensation
	BOOL m_notTrimTail;                       //不需要trim tail没有过扫描时配置，有劈缝的测量结果如果trimtail会出错
	BOOL m_regDime;                           //常规尺寸
};

