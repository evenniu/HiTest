#pragma once
// blade constands and enumerations

#define UM_SELECTCHANGED WM_USER + 2
#define UM_COPYDATA      WM_USER + 3

#define MAXFITS     10  // if this is made bigger than 10, need to add more property page string ID's
#define MAX_SPECIAL_ZONE_REGIONS 50

#define GOODTYPE   20170428 // a date in numeric format YYYYMMDD is perhaps a good way to track changes in MTH file formats

#define CVC  0  // position of curves in section curve array
#define CCC  1
#define LEC  2
#define TEC  3
#define MCC  4

#define EDGE_NORMAL  0
#define EDGE_SQUARE  1
#define EDGE_PARTIAL 2

#define NOMINAL           0
#define PLUSTOL           1
#define MINUSTOL          2
#define PLUSINTOL         3
#define MINUSINTOL        4
#define PLUSOUTTOL        5
#define MINUSOUTTOL       6

#define BLACK   (COLORREF)0x00000000
#define RED     (COLORREF)0x000000FF
#define GREEN   (COLORREF)0x0000FF00
#define BLUE    (COLORREF)0x00FF0000
#define CYAN    (COLORREF)0x00FFFF00
#define MAGENTA (COLORREF)0x00FF00FF
#define YELLOW  (COLORREF)0x0000FFFF
#define WHITE   (COLORREF)0x00FFFFFF

#define GRAY    (COLORREF)0x00777777
#define LTGRAY  (COLORREF)0x00BBBBBB
#define DKGRAY  (COLORREF)0x00555555

#define SHOWDEVS      0
#define SHOWACTUALS   1
#define SHOWNOMINALS  2

#define CLICKINSPECTED   0
#define CLICKANY         1
#define CLICKWIDTH       2

#define DIM_SKIP    -1   // methods of dimensioning
#define DIM_TOL      0
#define DIM_SPEC     1
#define DIM_NOMTOL   2
#define DIM_NOMSPEC  3

#define SHOW_DEV     0   // what to display
#define SHOW_ACT     1
#define SHOW_NOM     2
#define SHOW_OUT     3

enum BladeParams
{
  ParamLEOffset = 0,        // location to measure LE Thickness
  ParamLEOffset2,           // location to measure LE Thickness
  ParamLEOffset3,           // location to measure LE Thickness
  ParamLEOffset4,           // location to measure LE Thickness
  ParamLEOffset5,           // location to measure LE Thickness
  ParamLEOffset6,           // location to measure LE Thickness
  ParamLEOffset7,           // location to measure LE Thickness
  ParamLEOffset8,           // location to measure LE Thickness
  ParamLEOffset9,           // location to measure LE Thickness
  ParamLEOffset10,          // location to measure LE Thickness
  ParamTEOffset,            // location to measure TE Thickness
  ParamTEOffset2,           // location to measure TE Thickness
  ParamTEOffset3,           // location to measure TE Thickness
  ParamTEOffset4,           // location to measure TE Thickness
  ParamTEOffset5,           // location to measure TE Thickness
  ParamTEOffset6,           // location to measure TE Thickness
  ParamTEOffset7,           // location to measure TE Thickness
  ParamTEOffset8,           // location to measure TE Thickness
  ParamTEOffset9,           // location to measure TE Thickness
  ParamTEOffset10,          // location to measure TE Thickness
  ParamLEChange,            // location to start contours
  ParamTEChange,            // location to end contours
  ParamLEChange2,           // location to start contours
  ParamTEChange2,           // location to end contours
  ParamRootSect,            // name of root section for bow calcs
  ParamTipSect,             // name of tip section for bow calcs
  ParamRefSect,             // name of reference section for relative calcs
  ParamRootAngSect,         // name of root angle section
  ParamUnits,               // units, MM or INCH
  ParamRootAng,             // angle of root section
  ParamSetSpec,             // following tol entries are spec limits
  ParamSetTol,              // following tol entries are +/- tolerances
  ParamLEAngOffset,         // locations to measure LE Angle
  ParamTEAngOffset,         // locations to measure TE Angle
  ParamLEVarTol,            // LE tolerances for variable tolerancing
  ParamTEVarTol,            // TE tolerances for variable tolerancing
  ParamMidVarTol,           // Stack axis tolerances for variable tolerancing
  ParamChordVarTol,         // TE chord tolerances for variable tolerancing
  ParamGageOffsets,         // Gage Angle offset values
  ParamLEPoints,            // location to measure LE Thickness
  ParamLEPoints2,           // location to measure LE Thickness
  ParamLEPoints3,           // location to measure LE Thickness
  ParamLEPoints4,           // location to measure LE Thickness
  ParamLEPoints5,           // location to measure LE Thickness
  ParamLEPoints6,           // location to measure LE Thickness
  ParamLEPoints7,           // location to measure LE Thickness
  ParamLEPoints8,           // location to measure LE Thickness
  ParamLEPoints9,           // location to measure LE Thickness
  ParamLEPoints10,          // location to measure LE Thickness
  ParamTEPoints,            // location to measure TE Thickness
  ParamTEPoints2,           // location to measure TE Thickness
  ParamTEPoints3,           // location to measure TE Thickness
  ParamTEPoints4,           // location to measure TE Thickness
  ParamTEPoints5,           // location to measure TE Thickness
  ParamTEPoints6,           // location to measure TE Thickness
  ParamTEPoints7,           // location to measure TE Thickness
  ParamTEPoints8,           // location to measure TE Thickness
  ParamTEPoints9,           // location to measure TE Thickness
  ParamTEPoints10,          // location to measure TE Thickness
  ParamBowPlaneAng,         // angle of bow plane in XY plane from X axis, for singapore customer
  ParamBowPlaneX,           // x coord of a point on the bow plane
  ParamBowPlaneY,           // x coord of a point on the bow plane
  ParamLEZone,              // multipler for LE INZONE Calc
  ParamTEZone,              // multipler for TE INZONE Calc
  ParamExclude1,            // sides to exclude from best fit 1
  ParamExclude2,            // sides to exclude from best fit 2
  ParamExclude3,            // sides to exclude from best fit 3
  ParamExclude4,            // sides to exclude from best fit 4
  ParamExclude5,            // sides to exclude from best fit 5
  ParamExclude6,            // sides to exclude from best fit 5
  ParamExclude7,            // sides to exclude from best fit 5
  ParamExclude8,            // sides to exclude from best fit 5
  ParamExclude9,            // sides to exclude from best fit 5
  ParamExclude10,           // sides to exclude from best fit 5
  ParamExtremeAng1,         // extreme angle
  ParamExtremeAng2,         // ditto for second calculation
  ParamExtremeAng3,         // ditto for third calculation
  ParamExtremeAng4,         // ditto for fourth calculation
  ParamLEAngle,             // angle to find offset for LE Thickness
  ParamLEAngle2,            // angle to find offset for LE Thickness
  ParamLEAngle3,            // angle to find offset for LE Thickness
  ParamLEAngle4,            // angle to find offset for LE Thickness
  ParamLEAngle5,            // angle to find offset for LE Thickness
  ParamLEAngle6,            // angle to find offset for LE Thickness
  ParamLEAngle7,            // angle to find offset for LE Thickness
  ParamLEAngle8,            // angle to find offset for LE Thickness
  ParamLEAngle9,            // angle to find offset for LE Thickness
  ParamLEAngle10,           // angle to find offset for LE Thickness
  ParamTEAngle,             // angle to find offset for TE Thickness
  ParamTEAngle2,            // angle to find offset for TE Thickness
  ParamTEAngle3,            // angle to find offset for TE Thickness
  ParamTEAngle4,            // angle to find offset for TE Thickness
  ParamTEAngle5,            // angle to find offset for TE Thickness
  ParamTEAngle6,            // angle to find offset for TE Thickness
  ParamTEAngle7,            // angle to find offset for TE Thickness
  ParamTEAngle8,            // angle to find offset for TE Thickness
  ParamTEAngle9,            // angle to find offset for TE Thickness
  ParamTEAngle10,           // angle to find offset for TE Thickness
  ParamMinStockDist,        // Delta for checks for MINSTOCK
  ParamWaveWidth,           // width of window for waviness calculation
  ParamLERadii,             // radii of two arcs to swing for LE arc fits and arc contours
  ParamTERadii,             // radii of two arcs to swing for TE arc fits and arc contours
  ParamCamberOffsets,       // backoff distances used for computing MCL half-spaces (with the 2016 MCL algorithm)
  ParamLEFlatSize,          // amount of edge-flat for computing LE thickness direction
  ParamLEFlatSize2,         // amount of edge-flat for computing LE thickness direction
  ParamLEFlatSize3,         // amount of edge-flat for computing LE thickness direction
  ParamLEFlatSize4,         // amount of edge-flat for computing LE thickness direction
  ParamLEFlatSize5,         // amount of edge-flat for computing LE thickness direction
  ParamLEFlatSize6,         // amount of edge-flat for computing LE thickness direction
  ParamLEFlatSize7,         // amount of edge-flat for computing LE thickness direction
  ParamLEFlatSize8,         // amount of edge-flat for computing LE thickness direction
  ParamLEFlatSize9,         // amount of edge-flat for computing LE thickness direction
  ParamLEFlatSize10,        // amount of edge-flat for computing LE thickness direction
  ParamTEFlatSize,          // amount of edge-flat for computing TE thickness direction
  ParamTEFlatSize2,         // amount of edge-flat for computing TE thickness direction
  ParamTEFlatSize3,         // amount of edge-flat for computing TE thickness direction
  ParamTEFlatSize4,         // amount of edge-flat for computing TE thickness direction
  ParamTEFlatSize5,         // amount of edge-flat for computing TE thickness direction
  ParamTEFlatSize6,         // amount of edge-flat for computing TE thickness direction
  ParamTEFlatSize7,         // amount of edge-flat for computing TE thickness direction
  ParamTEFlatSize8,         // amount of edge-flat for computing TE thickness direction
  ParamTEFlatSize9,         // amount of edge-flat for computing TE thickness direction
  ParamTEFlatSize10,        // amount of edge-flat for computing TE thickness direction
  ParamComplexEdgeZone,     // specification for a complex (variable) tolerance zone near an edge
  ParamComplexEdgeZone2,    // specification for a complex (variable) tolerance zone near an edge
  ParamComplexEdgeZone3,    // specification for a complex (variable) tolerance zone near an edge
  ParamComplexEdgeZone4,    // specification for a complex (variable) tolerance zone near an edge
  ParamComplexEdgeZone5,    // specification for a complex (variable) tolerance zone near an edge
  ParamComplexEdgeZone6,    // specification for a complex (variable) tolerance zone near an edge
  ParamComplexEdgeZone7,    // specification for a complex (variable) tolerance zone near an edge
  ParamComplexEdgeZone8,    // specification for a complex (variable) tolerance zone near an edge
  ParamComplexEdgeZone9,    // specification for a complex (variable) tolerance zone near an edge
  ParamComplexEdgeZone10,   // specification for a complex (variable) tolerance zone near an edge
  ParamChordZone,           // specification for a complex (variable) tolerance zone along a chord
  ParamChordZone2,          // specification for a complex (variable) tolerance zone along a chord
  ParamChordZone3,          // specification for a complex (variable) tolerance zone along a chord
  ParamChordZone4,          // specification for a complex (variable) tolerance zone along a chord
  ParamChordZone5,          // specification for a complex (variable) tolerance zone along a chord
  ParamChordZone6,          // specification for a complex (variable) tolerance zone along a chord
  ParamChordZone7,          // specification for a complex (variable) tolerance zone along a chord
  ParamChordZone8,          // specification for a complex (variable) tolerance zone along a chord
  ParamChordZone9,          // specification for a complex (variable) tolerance zone along a chord
  ParamChordZone10,         // specification for a complex (variable) tolerance zone along a chord
  ParamCuppingLEOffset1,    // LE offset for cupping calculations
  ParamCuppingLEOffset2,    // LE offset for cupping calculations
  ParamCuppingLEOffset3,    // LE offset for cupping calculations
  ParamCuppingLEOffset4,    // LE offset for cupping calculations
  ParamCuppingLEOffset5,    // LE offset for cupping calculations
  ParamCuppingLEOffset6,    // LE offset for cupping calculations
  ParamCuppingLEOffset7,    // LE offset for cupping calculations
  ParamCuppingLEOffset8,    // LE offset for cupping calculations
  ParamCuppingLEOffset9,    // LE offset for cupping calculations
  ParamCuppingLEOffset10,   // LE offset for cupping calculations
  ParamCuppingTEOffset1,    // TE offset for cupping calculations
  ParamCuppingTEOffset2,    // TE offset for cupping calculations
  ParamCuppingTEOffset3,    // TE offset for cupping calculations
  ParamCuppingTEOffset4,    // TE offset for cupping calculations
  ParamCuppingTEOffset5,    // TE offset for cupping calculations
  ParamCuppingTEOffset6,    // TE offset for cupping calculations
  ParamCuppingTEOffset7,    // TE offset for cupping calculations
  ParamCuppingTEOffset8,    // TE offset for cupping calculations
  ParamCuppingTEOffset9,    // TE offset for cupping calculations
  ParamCuppingTEOffset10,   // TE offset for cupping calculations
  ParamCuppingTEOffset1_B,  // TE offset B for cupping calculations
  ParamCuppingTEOffset2_B,  // TE offset B for cupping calculations
  ParamCuppingTEOffset3_B,  // TE offset B for cupping calculations
  ParamCuppingTEOffset4_B,  // TE offset B for cupping calculations
  ParamCuppingTEOffset5_B,  // TE offset B for cupping calculations
  ParamCuppingTEOffset6_B,  // TE offset B for cupping calculations
  ParamCuppingTEOffset7_B,  // TE offset B for cupping calculations
  ParamCuppingTEOffset8_B,  // TE offset B for cupping calculations
  ParamCuppingTEOffset9_B,  // TE offset B for cupping calculations
  ParamCuppingTEOffset10_B, // TE offset B for cupping calculations
  ParamCuppingTEOffset1_C,  // TE offset C for cupping calculations
  ParamCuppingTEOffset2_C,  // TE offset C for cupping calculations
  ParamCuppingTEOffset3_C,  // TE offset C for cupping calculations
  ParamCuppingTEOffset4_C,  // TE offset C for cupping calculations
  ParamCuppingTEOffset5_C,  // TE offset C for cupping calculations
  ParamCuppingTEOffset6_C,  // TE offset C for cupping calculations
  ParamCuppingTEOffset7_C,  // TE offset C for cupping calculations
  ParamCuppingTEOffset8_C,  // TE offset C for cupping calculations
  ParamCuppingTEOffset9_C,  // TE offset C for cupping calculations
  ParamCuppingTEOffset10_C, // TE offset C for cupping calculations
  ParamConcavePositionReferenceAngle, // Angle for computing concave position
  ParamArraySize
};

enum BladeCalculations
{
  CalcUnset = -1,
  CalcBowX = 0,
  CalcBowY,
  CalcBowT,
  CalcLEContour,
  CalcLEContour2,
  CalcTEContour,
  CalcTEContour2,
  CalcCVContour,
  CalcCCContour,
  CalcLEThick,
  CalcLEThick2,
  CalcLEThick3,
  CalcLEThick4,
  CalcLEThick5,
  CalcLEThick6,
  CalcLEThick7,
  CalcLEThick8,
  CalcLEThick9,
  CalcLEThick10,
  CalcTEThick,
  CalcTEThick2,
  CalcTEThick3,
  CalcTEThick4,
  CalcTEThick5,
  CalcTEThick6,
  CalcTEThick7,
  CalcTEThick8,
  CalcTEThick9,
  CalcTEThick10,
  CalcMaxThick,
  CalcMaxThick2,
  CalcChordLength,
  CalcChordAngle,//32
  CalcRelChordAngle,
  CalcAdjChordAngle,
  CalcLEPosition,
  CalcLEPosition2,
  CalcTEPosition,
  CalcTEPosition2,
  CalcCCPosition,
  CalcLERadius,//40
  CalcTERadius,
  CalcStackX,
  CalcStackY,
  CalcStackT,
  CalcStackU,
  CalcStackV,
  CalcTwist,
  CalcRelStackX,
  CalcRelStackY,
  CalcRelStackT,
  CalcRelTwist,
  CalcRelTwist2,
  CalcLEAngle,
  CalcTEAngle,
  CalcMinForm,
  CalcMaxForm,
  CalcCVMin,//57
  CalcCVMax,
  CalcCCMin,
  CalcCCMax,
  CalcLEMin,
  CalcLEMax,
  CalcTEMin,
  CalcTEMax,
  CalcNumOutOfTol,
  CALC_PERCOUTTOL,
  CalcExtreme,
  CalcExtreme2,
  CalcExtreme3,
  CalcExtreme4,
  CalcAdjStackX,
  CalcAdjStackY,
  CalcAdjStackT,
  CalcAdjTwist,
  CalcMinStock,
  CalcMaxStock,
  CalcMinStockLE,
  CalcMaxStockLE,
  CalcMinStockTE,
  CalcMaxStockTE,
  CalcCentroidX,
  CalcCentroidY,
  CalcCentroidT,
  CalcRelCentroidX,
  CalcRelCentroidY,
  CalcRelCentroidT,
  CalcAdjCentroidX,
  CalcAdjCentroidY,
  CalcAdjCentroidT,
  CalcArea,
  CalcCVWaviness,
  CalcCCWaviness,
  CalcLEWaviness,
  CalcTEWaviness,
  CalcWaviness,
  CalcLEBowX,
  CalcLEBowY,
  CalcLEBowU,
  CalcLEBowV,
  CalcContArc,
  CalcContArcCV,
  CalcContArcCC,
  CalcContArc2,
  CalcContArcCV2,
  CalcContArcCC2,
  CalcTipLean,
  CalcBaseLean,
  CalcCupping1,
  CalcCupping2,
  CalcCupping3,
  CalcCupping4,
  CalcCupping5,
  CalcCupping6,
  CalcCupping7,
  CalcCupping8,
  CalcCupping9,
  CalcCupping10,
  CalcRelExtreme,
  CalcRelExtreme2,
  CalcRelExtreme3,
  CalcRelExtreme4,
  CalcRelLECHDTHCK,//TL-LECHDTHCK
  CalcRelLECALTHCK, //LECALTHCK
  CalcRelTECHDTHCK,//TT-TECHDTHCK
  CalcRelTECALTHCK, //TECALTHCK
  CalcReTEADDXAXISTHCK,//b - TECHDAXTHCK
  CalcReTECALAXISTHCK,   //TECALAXTHCK
  CalcReLEADDXAXISTHCK,  // LECHDAXTHCK
  CalcReLECALAXISTHCK,   // LECALAXTHCK
  CalcReTEOCLP,			//A- TETSTNTDIST
  CalcReTESTCALDist,    //TESTCALDIST
  CalcReLETSTNTDIST,    //LETSTNTDIST
  CalcReLEMCLDist,
  CalcReTEMCLDist,
  CalcReLETTHCK,
  CalcReTETTHCK,
  CalcReLERTHCK,
  CalcReLERTHCK2,
  CalcReLERTHCK3,
  CalcReLERTHCK4,
  CalcReTERTHCK,
  CalcReTERTHCK2,
  CalcReTERTHCK3,
  CalcReTERTHCK4,
  CalcReMAXTHKCHDV,
  CalcReLECCTANTHCK,//time:2021-12-1
  CalcReTECCTANTHCK,
  CalcReLECCTANAxTHCK,//136
  CalcReTECCTANAxTHCK,//137
  CalcReTHROATYCCTECV,//138 //time:2022-1-5
  CalcArraySize
};

enum BladeSpecialRequests
{
  SpecialLotPlot=0,
  SpecialZoneForm,
  SpecialDrillCorrection,
  SpecialShowTolerances,
  SpecialUseDMS,
  SpecialSaveASFile,
  SpecialSaveXYZFile,
  SpecialSaveIGESFile,
  SpecialSaveDataPage,
  SpecialSaveCSVFile,
  SpecialSaveCMMFile,
  SpecialAutoPrint,
  SpecialReportOption,
  SpecialAutoPlot,
  SpecialShowWhiskers,
  SpecialPlotEdges,
  SpecialSaveBCDFile,
  SpecialUseHeader,
  SpecialUseCustomerNote,
  SpecialSpecifyColumns,
  SpecialCreateExplodedPlots,
  SpecialMagnify,
  SpecialSaveMCLFile,
  SpecialSaveAVGFile,
  SpecialCreatePointExceptionReport,
  SpecialPrinterScale,
  SpecialSaveBestFitFile,
  SpecialUseMCLToSetLETransition,
  SpecialUseMCLToSetTETransition,
  SpecialUseSideToSetLETransition,
  SpecialUseSideToSetTETransition,
  SpecialEmbedZoneFormInMatrix,
  SpecialSaveDataPagePlus,
  SpecialPostExecutable,
  SpecialSaveDataPagePlus1,//SAVEDATAQDAS
  SpecialSaveRegDimData,//³£¹æ³ß´ç
  SpecialArraySize
};

#define TOTAL_SIZE ParamArraySize + CalcArraySize + SpecialArraySize

enum class BladeBestFitType
{
  BestFitNone = -1,                   // No Fit
  BestFitLeastSquares = 0,            // Least Squares fit
  BestFitConstrainedLeastSquares = 1, // Constrained Least Squares fit
  BestFitTwoPointMCL = 2,             // 5% - 95% fit
  BestFitTwoPointTailMCL = 3,         // 5% - 95% fit - from tail.
  BestFitFullBlade = 4,               // least squares fit performed to all sections
  BestFitLELeastSquares = 5,          // least squares on leading edge, points defined by two arcs.
  BestFitTELeastSquares = 6,          // least squares on trailing edge, points defined by two arcs.
  BestFitMinMax = 7,                  // Min Max Fit
  BestFitTwoPointForge = 8,           // forge fit, translation allowed only along forge axis
  BestFitFullBladeLELS = 9,           // all sections, only the leading edge
  BestFitLEComplexZoneLS = 10,        // least squares on leading edge, with complex zone
  BestFitTEComplexZoneLS = 11,        // least squares on trailing edge, with complex zone
  BestFitLEComplexZoneMM = 12,        // min-max on leading edge, with complex zone
  BestFitTEComplexZoneMM = 13,        // min-max on trailing edge, with complex zone
  BestFitCCChordZoneLS = 14,          // least squares on concave side, with chord zone
  BestFitCVChordZoneLS = 15,          // least squares on convex side, with chord zone
  BestFitCCChordZoneMM = 16,          // min-max on concave side, with chord zone
  BestFitCVChordZoneMM = 17,          // min-max on convex side, with chord zone
};

enum class BestFitAlgorithm
{
  None = 0,                         // no fit
  LeastSquares = 1,                 // least squares fit
  TwoPointsOnMCLFromNose = 2,       // two points on the MCL, measured from the nose
  TwoPointsOnMCLFromTail = 3,       // two points on the MCL, measured from the tail
  MinMax = 4,                       // min max fit
  TwoPointsOnMCLForForgedBlade = 5, // two points on the MCL, designed for forge fitting
};

// the two bow enumerations can be or'ed together
enum MethodBowReportEnds
{
  MethodBowZero =  0, // report root and tip sections as zero
  MethodBowDLB =   1  // P&W, displacement, lean, bow
};
enum MethodBowPoint
{
  MethodBowStack =    0, // use stack point after fit as bow point.
  MethodBowConvex =   2, // use intersection of datum plane (XZ for BOWY or YZ for BOWX) and CVX side
  MethodBowConcave =  4, // use intersection of datum plane (XZ for BOWY or YZ for BOWX) and CCV side
  MethodBowCentroid = 6  // use centroid
};

// the two Width methods can be or'ed together
enum MethodWidthReport
{
  MethodWidthNormal =  0, // for LE and TE widths, which component
  MethodWidthX      =  1,
  MethodWidthY      =  2,
  MethodWidthSquare =  3
};
enum MethodWidthReference
{
  MethodWidthTipArc      = 0, // for LE and TE widths, swing arc from...
  MethodWidthStackArc    = 4,
  MethodWidthStackX      = 8,
  MethodWidthStackY      = 12,
  MethodWidthPoints      = 16,
  MethodWidthEndX        = 20,
  MethodWidthEndY        = 24,
  MethodWidthAngleOffset = 28,
  MethodWidthConvexX     = 32,
  MethodWidthConcaveX    = 36,
  MethodWidthPercentX    = 40,
  MethodWidthOffsetFlat  = 44
};

enum MethodChord
{
  MethodChordTangent     = 0, // P&W M dimension - proj nose to tail distance
  MethodChordNoseTail    = 1, // distance from nose point to tail point
  MethodChordRadiiCenter = 2, // SW: distance from edge centers + radii
  MethodChordCaliper     = 3, // sit cc side, flag, find extremes
  MethodChordXAxis       = 4, // sit cc side, flag, find extremes
  MethodChordBelac       = 5, // special method for Belac
  MethodChordMinCircle   = 6 // MinimumEnclosingCircle
};

enum MethodChordAnglePoints
{
  MethodChordAngleNoseTail = 0, // angle from nose point to tail point
  MethodChordAngleCaliper  = 1, // sit cc side, flag, find extremes
  MethodChordAngleLEDist   = 2, // gage angle, LE offset and gage distance supplied
  MethodChordAngleLETE     = 3, // gage angle, LE and TE offsets supplied relative to the LE and TE.
  MethodChordAngleStack    = 4, // gage angle, LE and TE offsets supplied relative to the stack point
  MethodChordAngleBelac    = 5  // special method for BELAC
};
enum MethodChordAngleSide
{
  MethodChordAngleConcave =  0, // Chord Angle measured on Concave Side.
  MethodChordAngleConvex =   8, // Chord Angle measured on Convex Side.
  MethodChordAngleMCL =     16  // Chord Angle measured on MCL.
};

// the two end position methods can be or'ed together
enum MethodEndPositionNormalTo
{
  MethodEndPositionChord  = 0, // for LEPOS and TEPOS
  MethodEndPositionXAxis  = 1,
  MethodEndPositionYAxis  = 2,
  MethodEndPositionRadial = 3,
  MethodEndPositionEnds   = 4
};
enum MethodEndPositionReport
{
  MethodEndPositionStack  = 0,
  MethodEndPositionBow    = 8
};

enum MethodSidePositionNormalTo
{
  MethodSidePositionAngle = 0, // for CCPOS (and CVPOS if ever added)
  MethodSidePositionMidChord = 1,
};
enum MethodSidePositionReport
{
  MethodSidePositionStack = 0,
};

enum MethodEdgeRadius
{
  MethodRadiusInZone = 0  // use all points with a nominal radius of the end point
};

enum MethodStack
{
  MethodRefStackZero  = 0,
  MethodRefStackDatum = 1
};

// the two end angle methods can be or'ed together
enum MethodEndAngleMeasuredFrom
{
  MethodXPlus       = 0, // for LEANGLE and TEANGLE, which axis?
  MethodXMinus      = 1,
  MethodYPlus       = 2,
  MethodYMinus      = 3, // first five also used for EXTREME
  MethodRotation    = 4,
  MethodSquareChord = 5,
  MethodConvexSide  = 6,
  MethodConcaveSide = 7,
  MethodConcaveConvexMidpoint = 8
};
enum MethodEndAngleAlignment
{
  MethodMeasured    = 0,
  MethodBestFit     = 16
};

enum MethodMinMaxForm // for MINFORM and MAXFORM
{
  MethodFormConstant        = 0, // constant tolerance
  MethodFormGEVariable      = 1, // variable tolerance (GE spec)
  MethodFormVariableNomFile = 2, // variable defined in nominal file
  MethodFormVariableVarFile = 3  // variable tolerance defined in variable tol file 
};

enum MethodContourReporting
{
  MethodContourReportMaxMinusMin = 0, // report (MAX - MIN)
  MethodContourReportDoubleWorstDeviation = 16, // report the larger of 2*abs(MAX) and 2*abs(MIN)
  MethodContourReportDeviationRatio = 32 // report the ratio of the enveloping zone to the original zone
};

enum MethodLeanCalculations // for CalcTipLean and CalcBaseLean
{
  MethodLeanCaliper = 0 // only one calculation method for now
};
enum MethodLeanReporting // for CalcTipLean and CalcBaseLean
{
  MethodLeanAllSections = 0, // report all sections
  MethodLeanRootSectionOnly = 16, // report only the root section
  MethodLeanTipSectionOnly = 32 // report only the tip section
};

enum MethodCuppingReporting // for CalcCupping
{
  MethodCuppingFullConvex = 0,
  MethodCuppingFullConcave = 16,
  MethodCuppingTrailingConvex = 32,
  MethodCuppingTrailingConcave = 48
};

enum MethodThicknessCalculations
{
  MethodThicknessMaxInscribed = 0, // default calculation method
  MethodThicknessTwoPoint     = 1 // new method requested by China in BRK-639
};

// Matching Width Methods, except POINTS and ANGLEOFFSET. The user can define how to determine the LE and TE transition points
enum EndTransitionPointOnMCL
{
  TransitionPointTipArc  = 0, // for LE and TE widths, swing arc from...
  TransitionPointStackArc = 1,
  TransitionPointStackX = 2,
  TransitionPointStackY = 3,
  TransitionPointEndX   = 4,
  TransitionPointEndY   = 5
};
enum EndTransitionPointIntersectDirection
{
  TransitionPointDirectionNormal  = 0,
  TransitionPointDirectionX       = 1,
  TransitionPointDirectionY       = 2,
  TransitionPointDirectionSquare  = 3
};

enum StackOriginSpecialPoints
{
  StackOriginPointGhostSection = -1723
};
