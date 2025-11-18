#include "stdafx.h"
#include "Circle.h"
#include "SubCurve.h"
//#include "BestFit.h"

#include "SmallestCircle.h"
#include "SECTION.h"
#include "SectionCurve.h"
#include "MINMAX.H"
#include "BestFits.h"
#include "ArraySlicing.h"
#include "HermiteCurve.h"
#include "EigenAbstractCurve.h"
#include "TemplateHermiteSpline.h"
#include "ToleranceSection.h"
#pragma warning(push)
#pragma warning(disable : 4267)
#pragma warning(disable : 4244)
#pragma warning(disable : 4100)
#pragma warning(pop)

#include "CurvePolygon.h"
#include "MeanCamberCurve.h"
using namespace Hexagon;

#ifdef _DEBUG
#ifndef DBG_NEW
#define DBG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DBG_NEW
#endif
#endif // _DEBUG
CSection::CSection()
{
    m_name[0] = 0;
    m_numNomPoints = 0;
    m_zValue = 0.0;
    m_leType = m_teType = EDGE_NORMAL;
    m_nomCurve = m_meaCurve = m_BCCurve = NULL;
    //m_minStock = NULL;
    m_tolSegPoints_start = NULL;
    m_tolSegPoints_end = NULL;
    for (int i = 0; i < 5; i++)
    {
        m_nomPart[i] = NULL;
        m_meaPart[i] = NULL;
    }

    for (int i = 0; i < 3; i++)
    {
        m_nomPitch[i] = 0.0;
        m_meaPitch[i] = 0.0;
    }

    m_cxpt = m_cypt = m_czpt = 0;
    m_ival = 0;
    m_jval = 0;
    m_kval = 0;
    m_mxpt = 0;
    m_mypt = 0;
    m_nxpt = 0;
    m_nypt = 0;
    m_nomt = 0;
    m_partOf = 0;
    m_skewalign = 0;
    m_skewReport = 0;
    m_nomx = 0;
    m_nomy = 0;
    m_nomi = 0;
    m_nomj = 0;
    m_nomk = 0;
    m_mtol = 0;
    m_ptol = 0;
    m_tolsegCount = 0;
    m_start_point = 0;
    m_end_point = 0;
    m_curvature_coef = 0;
    m_tolSegPoints_start = 0;
    m_tolSegPoints_end = 0;
    m_pTol_start = 0;
    m_pTol_end = 0;
    m_mTol_start = 0;
    m_mTol_end = 0;
    m_t_start = 0;
    m_t_end = 0;
    nominalMCLParams = nullptr;

    for (int i = 0; i < 10; i++)
    {
        m_leNomWid1[i][0] = m_leNomWid1[i][1] = 1.0e20;
        m_teNomWid1[i][0] = m_teNomWid1[i][1] = 1.0e20;
        m_leNomWid2[i][0] = m_leNomWid2[i][1] = 1.0e20;
        m_teNomWid2[i][0] = m_teNomWid2[i][1] = 1.0e20;
        m_leActWid1[i][0] = m_leActWid1[i][1] = 1.0e20;
        m_teActWid1[i][0] = m_teActWid1[i][1] = 1.0e20;
        m_leActWid2[i][0] = m_leActWid2[i][1] = 1.0e20;
        m_teActWid2[i][0] = m_teActWid2[i][1] = 1.0e20;

        m_cuppingResults[i].cuppingResultIsFilledOut = false;
        m_cuppingResults[i].camberPointA[0] = m_cuppingResults[i].camberPointA[1] = 0.0;
        m_cuppingResults[i].camberPointB[0] = m_cuppingResults[i].camberPointB[1] = 0.0;
        m_cuppingResults[i].camberPointC[0] = m_cuppingResults[i].camberPointC[1] = 0.0;
        m_cuppingResults[i].surfacePointA[0] = m_cuppingResults[i].surfacePointA[1] = 0.0;
        m_cuppingResults[i].surfacePointB[0] = m_cuppingResults[i].surfacePointB[1] = 0.0;
        m_cuppingResults[i].surfacePointA[0] = m_cuppingResults[i].surfacePointA[1] = 0.0;
        m_cuppingResults[i].circleCenter[0] = m_cuppingResults[i].circleCenter[1] = 0.0;
        m_cuppingResults[i].circleRadius = 0.0;
    }

    m_meaPitch[2] = m_nomPitch[2] = -1.0;

   // m_compArray = NULL;
    BOOL dumpComp = 1;//myGetProfileInt(L"SaveCompInfo", 1);
    if (dumpComp)
        //m_compArray = new CCompRecordArray();

    //for (int i = 0; i < 100; i++)
    //    m_bestFits[i] = NULL;
    m_numBestFits = 0;

    m_arcRangeCV[0] = m_arcRangeCV[1] = m_arcRangeCC[0] = m_arcRangeCC[1] = -1.0;
    m_openareaCount = 0;
    for (int i = 0; i < 3; i++)
    {
        m_areaStartPoint[i][0] = 0.0;
        m_areaStartPoint[i][1] = 0.0;
        m_areaEndPoint[i][0] = 0.0;
        m_areaEndPoint[i][1] = 0.0;
        m_areaIndex[i][0] = 0;
        m_areaIndex[i][1] = 0;
    }
    m_fixedAxis = -1;
}
CSection::~CSection()
{
    if (m_nomCurve)
        delete m_nomCurve;

    if (m_meaCurve)
        delete m_meaCurve;

    if (m_BCCurve)
        delete m_BCCurve;

    if (m_skewalign)
        delete m_skewalign;

    if (m_cxpt)
        delete[] m_cxpt;
    if (m_cypt)
        delete[] m_cypt;
    if (m_czpt)
        delete[] m_czpt;

    if (m_ival)
        delete[] m_ival;
    if (m_jval)
        delete[] m_jval;
    if (m_kval)
        delete[] m_kval;

    if (m_mxpt)
        delete[] m_mxpt;
    if (m_mypt)
        delete[] m_mypt;

    if (m_nxpt)
        delete[] m_nxpt;
    if (m_nypt)
        delete[] m_nypt;

    if (m_nomt)
        delete[] m_nomt;

    if (m_partOf)
        delete[] m_partOf;

    if (m_nomx)
        delete[] m_nomx;
    if (m_nomy)
        delete[] m_nomy;
    if (m_nomi)
        delete[] m_nomi;
    if (m_nomj)
        delete[] m_nomj;
    if (m_nomk)
        delete[] m_nomk;
    if (m_mtol)
        delete[] m_mtol;
    if (m_ptol)
        delete[] m_ptol;
    if (nominalMCLParams)
        delete nominalMCLParams;

    int i;
    for (i = 0; i < 5; i++)
    {
        if (m_nomPart[i])
            delete m_nomPart[i];
        if (m_meaPart[i])
            delete m_meaPart[i];
    }

    //for (i = 0; i < m_numBestFits; i++)
    //    if (m_bestFits[i])
    //        delete m_bestFits[i];

    //if (m_compArray->size() > 0)
    //{
    //    for (i = 0; i < (int)m_compArray->size(); i++)
    //        delete m_compArray->at(i);
    //}

   // m_compArray->clear();
    //delete m_compArray;


    if (m_start_point)
        delete[] m_start_point;
    if (m_end_point)
        delete[] m_end_point;
    if (m_curvature_coef)
        delete[] m_curvature_coef;
    if (m_pTol_start)
        delete[] m_pTol_start;
    if (m_pTol_end)
        delete[] m_pTol_end;
    if (m_mTol_start)
        delete[] m_mTol_start;
    if (m_mTol_end)
        delete[] m_mTol_end;
    if (m_tolSegPoints_start)
        delete m_tolSegPoints_start;
    if (m_tolSegPoints_end)
        delete m_tolSegPoints_end;
    if (m_t_start)
        delete[] m_t_start;
    if (m_t_end)
        delete[] m_t_end;

}
void CSection::AddNomXYIJK(int i, double* xyijk)
{
    if (i < 0 || i >= m_numNomPoints)
        return;

    m_nomx[i] = xyijk[0];
    m_nomy[i] = xyijk[1];
    m_nomi[i] = xyijk[2];
    m_nomj[i] = xyijk[3];
    m_nomk[i] = xyijk[4];
}
double CSection::FindKValue(double* pt)
{
    if (m_numNomPoints < 1)
        return 0.0;

    int i, bi = 0, ci = 0;
    double minD = 1.0e20; // closest point that projects onto a nominal segment.
    double minC = 1.0e20; // closest point
    double br = -2000.0, kv = 0.0;
    for (i = 0; i < m_numNomPoints; i++)
    {
        int i2 = (i + 1) % m_numNomPoints;

        double dx = pt[0] - m_nomx[i];
        double dy = pt[1] - m_nomy[i];
        double d = dx * dx + dy * dy;
        if (d < minC) // save closest point in case pt doesn't project onto any segment.
        {
            minC = d;
            ci = i;
        }

        double lp[2], lv[2];
        lp[0] = m_nomx[i];
        lp[1] = m_nomy[i];
        lv[0] = m_nomx[i2] - m_nomx[i];
        lv[1] = m_nomy[i2] - m_nomy[i];
        double seglen = normalize(lv, lv);

        d = ptlinedist(lp, lv, pt);

        if (d < minD) // closest so far...
        {
            // but need to make sure that pt projects onto line segment [i, i+1]

            double pd = projdist(lp, lv, pt);
            double r = pd / seglen;
            if (r >= 0.0 && r <= 1.0)
            {
                minD = d;
                bi = i;
                kv = m_nomk[i] + r * (m_nomk[i2] - m_nomk[i]);
                br = r;
            }
        }
    }

    if (br < -1000.0 || (2 * minC < minD && minC > 0.001)) // point fell in a gap.  just use k value of closest point
    {
        kv = m_nomk[ci];
        // bugout(0, _T("9 %f %f %f %f  FindKValue segment GAP"), pt[0], pt[1], m_nomx[ci], m_nomy[ci]);
        // bugout(0, _T("10 %f %f GAP"), pt[0], pt[1]);
    }
    else
    {
        // bugout(0, _T("9 %f %f %f %f  FindKValue segment r %f"), pt[0], pt[1], m_nomx[bi], m_nomy[bi], br);
    }

    return kv;
}
void CSection::MakeNomArrays(int numPts, const Hexagon::Blade::MeanCamberCurveParameters2016* mclParams)
{
    m_numNomPoints = numPts;

    if (m_numNomPoints < 1)
        return;

    nominalMCLParams = mclParams ? new Blade::MeanCamberCurveParameters2016(*mclParams) : nullptr;

    m_nomx = new double[numPts];
    m_nomy = new double[numPts];
    m_nomi = new double[numPts];
    m_nomj = new double[numPts];
    m_nomk = new double[numPts];
    m_mtol = new double[numPts];
    m_ptol = new double[numPts];
}
void CSection::AddTol(int i, double* mptol)
{
    if (i < 0 || i >= m_numNomPoints)
        return;

    m_mtol[i] = mptol[0];
    m_ptol[i] = mptol[1];
}
Eigen::Matrix2Xd constructPointMatrix(const double* xValues, const double* yValues, const ptrdiff_t n)
{
    Eigen::Matrix2Xd points(2, n);
    points.row(0) = Eigen::Map<const Eigen::RowVectorXd>(xValues, n);
    points.row(1) = Eigen::Map<const Eigen::RowVectorXd>(yValues, n);
    return points;
}
double minMaxObjective(const CBestFit& bf, const CFitParams& fitParams)
{
    return 0;
    /*std::vector<double> objectives;
    for (int side = 0; side < 4; side++)
    {
        if (fitParams.fitcurve[side] && fitParams.weightcurve[side] > 0.0)
        {
            objectives.push_back(std::max(std::abs(bf.m_mindev[side]), std::abs(bf.m_maxdev[side])));
        }
    }
    if (objectives.empty())
    {
        return 0.0;
    }
    return *std::max_element(objectives.begin(), objectives.end());*/
}
void updateBestFit_measuredPointsToNominalCurve(CBestFit* bf, const Hexagon::Blade::SectionCurve& nominalCurves,
                                                const Eigen::Isometry2d& measuredToNominalTransform,
                                                const Eigen::Ref<const Eigen::Matrix2Xd>& measuredPoints,
                                                const Eigen::Ref<const Eigen::ArrayXb>& isUsedInFit,
                                                const CFitParams& fitParams, double* mtols = nullptr,
                                                double* ptols = nullptr)
{
    // figure out the closest nominal points
    const ptrdiff_t N = measuredPoints.cols();
    const Eigen::Matrix2Xd alignedMeasuredPoints = measuredToNominalTransform * measuredPoints;
    Eigen::VectorXd nominalT(N);
    nominalCurves.whole->FindClosestTValues(nominalT.data(), alignedMeasuredPoints.data(), N);
    Eigen::Matrix2Xd nominalPoints(2, N);
    Eigen::Matrix2Xd nominalTangents(2, N);
    nominalCurves.whole->CalcPoints(nominalPoints.data(), nominalT.data(), N, nominalTangents.data());
    Eigen::Matrix2Xd ijk(2, N);
    ijk.row(0) = nominalTangents.row(1);
    ijk.row(1) = -nominalTangents.row(0);
    ijk.colwise().normalize();
    const Eigen::ArrayXd distances =
        (alignedMeasuredPoints - nominalPoints).cwiseProduct(ijk).colwise().sum().transpose();

    // assign the points to a curve, after fitting
    Eigen::ArrayXi fittedBestPartOf(N);
    fittedBestPartOf =
        Hexagon::Blade::tIsInSubcurve_eigen(nominalT, *nominalCurves.leading, nominalCurves.whole->period())
        .select(LEC, fittedBestPartOf);

    // construct the result
    for (int i = 0; i < N; i++)
    {
        bf->PutVal(i, measuredPoints(0, i), measuredPoints(1, i));
        bf->PutNom(i, nominalPoints(0, i), nominalPoints(1, i));
        bf->PutT(i, nominalT[i]);
        bf->PutVec(i, ijk(0, i), ijk(1, i));
        bf->PutInf(i, alignedMeasuredPoints(0, i), alignedMeasuredPoints(1, i));
        bf->m_bestPartOf[i] = fittedBestPartOf[i];
        bf->m_valWasUsedInFit[i] = isUsedInFit[i];
    }

    // add the summaries
    bf->m_totalBad = 0;
    bf->m_totalChecked = 0;
    for (int side = 0; side < 4; side++)
    {
        const Eigen::ArrayXd sideDeviations = Hexagon::Blade::sliceVector(distances, fittedBestPartOf == side);
        if (sideDeviations.size() > 0)
        {
            bf->m_mindev[side] = sideDeviations.minCoeff();
            bf->m_maxdev[side] = sideDeviations.maxCoeff();
            bf->m_meandev[side] = sideDeviations.mean();
        }
        if (sideDeviations.size() > 1)
        {
            bf->m_stddev[side] = std::sqrt((sideDeviations - sideDeviations.mean()).cwiseAbs2().sum() /
                static_cast<double>(sideDeviations.size() - 1));
        }
        if (mtols && ptols)
        {
            bf->m_totalChecked += static_cast<int>(sideDeviations.size());
            bf->m_totalBad += static_cast<int>((sideDeviations < mtols[side] && sideDeviations > ptols[side]).count());
        }
    }
    // add the transformation as well
    bf->m_align = Hexagon::Blade::toCAlignment(measuredToNominalTransform);
    bf->m_fitParams = fitParams;
}

CBestFit* createMeasuredPointsToNominalCurveBestFit(const Hexagon::Blade::SectionCurve& nominalCurves,
                                                    const Eigen::Isometry2d& measuredToNominalTransform,
                                                    const Eigen::Ref<const Eigen::Matrix2Xd>& measuredPoints,
                                                    const Eigen::Ref<const Eigen::ArrayXb>& isUsedInFit,
                                                    const CFitParams& fitParams, double* mtols = nullptr,
                                                    double* ptols = nullptr)
{
    auto result = std::make_unique<CBestFit>(static_cast<int>(measuredPoints.cols()));
    updateBestFit_measuredPointsToNominalCurve(result.get(), nominalCurves, measuredToNominalTransform, measuredPoints,
        isUsedInFit, fitParams, mtols, ptols);
    return result.release();
}
const double infinity = std::numeric_limits<double>::infinity();

struct ChordInformation
{
    Eigen::Vector2d leadingPoint, trailingPoint, leadingCenter, trailingCenter, leadingVector, trailingVector;
};
void setRotationOptions(Hexagon::Blade::FitOptions& options, const CFitParams& fp)
{
    // what kinds of rotation are allowed?
    switch (fp.rotfit)
    {
    case 0: // no limits to rotation
        options.allowRotation = true;
        options.rotationLimits = Eigen::Vector2d(-infinity, infinity);
        break;
    case 1: // no rotation allowed at all
        options.allowRotation = false;
        break;
    case 2: // rotation allowed within limits
        options.allowRotation = true;
        // convert the rotation limits to radians
        options.rotationLimits = M_PI * Eigen::Vector2d(fp.rotMTol, fp.rotPTol) / 180.0;
        break;
    default:
        throw std::logic_error("This should be impossible.");
    }
}


void setTranslationOptions(Hexagon::Blade::FitOptions& options, const CFitParams& fp, const ChordInformation& target,
    const ChordInformation& fitted)
{
    switch (fp.tranfit)
    {
    case 0: // no limits to rotation
    case 4: // no limits to rotation
    case 5: // no limits to rotation
        break;
    case 1: // no translation allowed at all
        options.allowTranslation = false;
        break;
    case 2: // translation allowed within limits
        options.translationXLimits = Eigen::Vector2d(fp.tranMTol[0], fp.tranPTol[0]);
        options.translationYLimits = Eigen::Vector2d(fp.tranMTol[1], fp.tranPTol[1]);
        break;
    case 3: // a pivot point is selected; no translation allowed
    {
        options.allowTranslation = false;
        switch (fp.pivot)
        {
            case 0: // pivot about LE center
                options.targetCurvePivotPoint = target.leadingCenter;
                options.pointsPivotPoint = fitted.leadingCenter;
                break;
            case 1: // pivot about LE nose
                options.targetCurvePivotPoint = target.leadingPoint;
                options.pointsPivotPoint = fitted.leadingPoint;
                break;
            case 2: // pivot about TE center
                options.targetCurvePivotPoint = target.trailingCenter;
                options.pointsPivotPoint = fitted.trailingCenter;
                break;
            case 3: // pivot about TE tail
                options.targetCurvePivotPoint = target.trailingPoint;
                options.pointsPivotPoint = fitted.trailingPoint;
                break;
            default:
                throw std::logic_error("This should be impossible.");
        }
    }
    break;
    default:
        throw std::logic_error("This should be impossible.");
    }
}

Eigen::Isometry2d createGuessTransform(const CFitParams& fp, const ChordInformation& targetChord,
    const ChordInformation& fittedChord, const Hexagon::Blade::SectionCurve& fittedCurves, const int leType,
    const int teType)
{
    const double targetChordLength = (targetChord.leadingPoint - targetChord.trailingPoint).norm();
    const double r = 0.75 * targetChordLength;
    Eigen::Vector2d targetPoint, fittedPoint;
    bugout(0, L"createGuessTransform:pivot(%d) tranfit(%d) entered ****", fp.pivot, fp.tranfit);
    if (fp.pivot == 0 && fp.tranfit == 3) // le center
    {
       // targetCurves.meanCamber->CircIntersect(targetChord.leadingPoint.data(), r, targetPoint.data());
        fittedCurves.meanCamber->CircIntersect(fittedChord.leadingPoint.data(), r, fittedPoint.data());
 /*       return Hexagon::Blade::twoPointBestFit(targetChord.leadingCenter, targetPoint, fittedChord.leadingCenter,
            fittedPoint);*/
    }
    //else if (fp.pivot == 1 && fp.tranfit == 3) // le nose
    //{
    //    targetCurves.meanCamber->CircIntersect(targetChord.leadingPoint.data(), r, targetPoint.data());
    //    fittedCurves.meanCamber->CircIntersect(fittedChord.leadingPoint.data(), r, fittedPoint.data());
    //    return Hexagon::Blade::twoPointBestFit(targetChord.leadingPoint, targetPoint, fittedChord.leadingPoint,
    //        fittedPoint);
    //}
    //else if (fp.pivot == 2 && fp.tranfit == 3) // te center
    //{
    //    targetCurves.meanCamber->CircIntersect(targetChord.trailingPoint.data(), r, targetPoint.data());
    //    fittedCurves.meanCamber->CircIntersect(fittedChord.trailingPoint.data(), r, fittedPoint.data());
    //    return Hexagon::Blade::twoPointBestFit(targetChord.trailingCenter, targetPoint, fittedChord.trailingCenter,
    //        fittedPoint);
    //}
    //else if (fp.pivot == 3 && fp.tranfit == 3) // te tail
    //{
    //    targetCurves.meanCamber->CircIntersect(targetChord.trailingPoint.data(), r, targetPoint.data());
    //    fittedCurves.meanCamber->CircIntersect(fittedChord.trailingPoint.data(), r, fittedPoint.data());
    //    return Hexagon::Blade::twoPointBestFit(targetChord.trailingPoint, targetPoint, fittedChord.trailingPoint,
    //        fittedPoint);
    //}
    //else if (fp.fitcurve[LEC] == 1 && fp.fitcurve[TEC] == 0)
    //{
    //    // LE and no TE: let's start with nose points aligned
    //    if (leType == EDGE_NORMAL || leType == EDGE_SQUARE)
    //    {
    //        targetCurves.whole->CircIntersect(targetChord.leadingPoint.data(), r, targetPoint.data());
    //        fittedCurves.whole->CircIntersect(fittedChord.leadingPoint.data(), r, fittedPoint.data());
    //        return Hexagon::Blade::twoPointBestFit(targetChord.leadingPoint, targetPoint, fittedChord.leadingPoint,
    //            fittedPoint);
    //    }
    //    else // partial edge
    //    {
    //        Eigen::Vector2d target0, target1, fitted0, fitted1;
    //        targetCurves.leading->CalcPoint(target0.data(), targetCurves.leading->T0());
    //        targetCurves.leading->CalcPoint(target1.data(), targetCurves.leading->T1());
    //        fittedCurves.leading->CalcPoint(fitted0.data(), fittedCurves.leading->T0());
    //        fittedCurves.leading->CalcPoint(fitted1.data(), fittedCurves.leading->T1());
    //        return Hexagon::Blade::twoPointBestFit(target0, target1, fitted0, fitted1);
    //    }
    //}
    //else if (fp.fitcurve[LEC] == 0 && fp.fitcurve[TEC] == 1)
    //{
    //    // TE and no LE: let's start with tail points aligned
    //    if (teType == EDGE_NORMAL || teType == EDGE_SQUARE)
    //    {
    //        targetCurves.meanCamber->CircIntersect(targetChord.trailingPoint.data(), r, targetPoint.data());
    //        fittedCurves.meanCamber->CircIntersect(fittedChord.trailingPoint.data(), r, fittedPoint.data());
    //        return Hexagon::Blade::twoPointBestFit(targetChord.trailingPoint, targetPoint, fittedChord.trailingPoint,
    //            fittedPoint);
    //    }
    //    else // partial edge
    //    {
    //        Eigen::Vector2d target0, target1, fitted0, fitted1;
    //        targetCurves.trailing->CalcPoint(target0.data(), targetCurves.trailing->T0());
    //        targetCurves.trailing->CalcPoint(target1.data(), targetCurves.trailing->T1());
    //        fittedCurves.trailing->CalcPoint(fitted0.data(), fittedCurves.trailing->T0());
    //        fittedCurves.trailing->CalcPoint(fitted1.data(), fittedCurves.trailing->T1());
    //        return Hexagon::Blade::twoPointBestFit(target0, target1, fitted0, fitted1);
    //    }
    //}
    //else if (fp.fitcurve[LEC] == 1 && fp.fitcurve[TEC] == 1 && leType == EDGE_PARTIAL)
    //{
    //    targetCurves.meanCamber->CircIntersect(targetChord.trailingPoint.data(), r, targetPoint.data());
    //    fittedCurves.meanCamber->CircIntersect(fittedChord.trailingPoint.data(), r, fittedPoint.data());
    //    return Hexagon::Blade::twoPointBestFit(targetChord.trailingPoint, targetPoint, fittedChord.trailingPoint,
    //        fittedPoint);
    //}
    //else if (fp.fitcurve[LEC] == 1 && fp.fitcurve[TEC] == 1 && teType == EDGE_PARTIAL)
    //{
    //    targetCurves.meanCamber->CircIntersect(targetChord.leadingPoint.data(), r, targetPoint.data());
    //    fittedCurves.meanCamber->CircIntersect(fittedChord.leadingPoint.data(), r, fittedPoint.data());
    //    return Hexagon::Blade::twoPointBestFit(targetChord.leadingPoint, targetPoint, fittedChord.leadingPoint,
    //        fittedPoint);
    //}
    //else if (fp.fitcurve[CCC] == 1 && fp.fitcurve[CVC] == 0 && fp.fitcurve[LEC] == 0 && fp.fitcurve[TEC] == 0)
    //{
    //    // fit only the CCC curve; start computing the endpoints
    //    const Eigen::Vector2d target0 = Hexagon::Blade::evaluate(*targetCurves.concave, targetCurves.concave->t0());
    //    const Eigen::Vector2d target1 = Hexagon::Blade::evaluate(*targetCurves.concave, targetCurves.concave->t1());
    //    const Eigen::Vector2d fitted0 = Hexagon::Blade::evaluate(*fittedCurves.concave, fittedCurves.concave->t0());
    //    const Eigen::Vector2d fitted1 = Hexagon::Blade::evaluate(*fittedCurves.concave, fittedCurves.concave->t1());

    //    // find the midpoints and direction vectors
    //    const Eigen::Vector2d midTarget = 0.5 * (target0 + target1);
    //    const Eigen::Vector2d midFitted = 0.5 * (fitted0 + fitted1);
    //    const Eigen::Vector2d diffTarget = target0 - target1;
    //    const Eigen::Vector2d diffFitted = fitted0 - fitted1;
    //    const Eigen::Vector2d midDirectionTarget = midTarget + diffTarget;
    //    const Eigen::Vector2d midDirectionFitted = midFitted + diffFitted;

    //    // return the two-point fit
    //    return Hexagon::Blade::twoPointBestFit(midTarget, midDirectionTarget, midFitted, midDirectionFitted);
    //}
    //else if (fp.fitcurve[CVC] == 1 && fp.fitcurve[CCC] == 0 && fp.fitcurve[LEC] == 0 && fp.fitcurve[TEC] == 0)
    //{
    //    // fit only the CCC curve; start computing the endpoints
    //    const Eigen::Vector2d target0 = Hexagon::Blade::evaluate(*targetCurves.convex, targetCurves.convex->t0());
    //    const Eigen::Vector2d target1 = Hexagon::Blade::evaluate(*targetCurves.convex, targetCurves.convex->t1());
    //    const Eigen::Vector2d fitted0 = Hexagon::Blade::evaluate(*fittedCurves.convex, fittedCurves.convex->t0());
    //    const Eigen::Vector2d fitted1 = Hexagon::Blade::evaluate(*fittedCurves.convex, fittedCurves.convex->t1());

    //    // find the midpoints and direction vectors
    //    const Eigen::Vector2d midTarget = 0.5 * (target0 + target1);
    //    const Eigen::Vector2d midFitted = 0.5 * (fitted0 + fitted1);
    //    const Eigen::Vector2d diffTarget = target0 - target1;
    //    const Eigen::Vector2d diffFitted = fitted0 - fitted1;
    //    const Eigen::Vector2d midDirectionTarget = midTarget + diffTarget;
    //    const Eigen::Vector2d midDirectionFitted = midFitted + diffFitted;

    //    // return the two-point fit
    //    return Hexagon::Blade::twoPointBestFit(midTarget, midDirectionTarget, midFitted, midDirectionFitted);
    //}

    // no special guess
    return Eigen::Isometry2d::Identity();

}




// create a function that looks for t-values of closest approach to a point
template <class TreeType>
double findNearestTValue(CCurve* curve, const TreeType& tree, ptrdiff_t numberOfPointsInTree, double* treeTValues,
    double* point)
{
    ptrdiff_t index;
    double squaredDistance;
    tree.query(point, 1, &index, &squaredDistance);
    double lowT = (index > 0) ? treeTValues[index - 1] : treeTValues[numberOfPointsInTree - 1] - curve->Period();
    double highT = (index < numberOfPointsInTree - 1) ? treeTValues[index + 1] : treeTValues[0] + curve->Period();
    double t;
    Eigen::Vector2d trash;
    curve->NewClosestPoint(point, trash.data(), &t, nullptr, lowT, highT, 1);
    return t;
}

/// <summary>
/// 
/// </summary>
/// <param name="index"></param>
/// <param name="inchSize"></param>
/// <param name="mtols"></param>
/// <param name="ptols"></param>
/// <returns></returns>
bool CSection::FitPoints(CFitParams& fp,int& index, double inchSize, double* mtols, double* ptols)
{
    const auto sectionCurve =Hexagon::Blade::nominalSectionCurve(this);
    const Eigen::Matrix2Xd measuredPoints = constructPointMatrix(m_mxpt, m_mypt, m_totalPoints);
    bugout(0, L"FitPoints:m_totalPoints(%d) entered ****", m_totalPoints);
    for (const auto side : { LEC, TEC })
    {
        if (fp.weightcurve[side] == 0)
        {
            fp.fitcurve[side] = 0;
        }
    }
    const Eigen::Map<const Eigen::ArrayXi> partOf(m_partOf, m_totalPoints);
    // no fit
    if (fp.algorithm == BestFitAlgorithm::None)
    {
        index = m_numBestFits;
        
        m_bestFits[m_numBestFits] = createMeasuredPointsToNominalCurveBestFit(sectionCurve, Eigen::Isometry2d::Identity(), measuredPoints,
                Eigen::ArrayXb::Ones(m_totalPoints), fp, mtols, ptols);
        m_numBestFits++;
        bugout(0, L"FitPoints:entered * m_numBestFits=%d ***", m_numBestFits);

        return true;
    }

    ChordInformation nominalChordInfo;
    //if (!Chord(0, nominalChordInfo.leadingPoint.data(), nominalChordInfo.trailingPoint.data(),
    //    nominalChordInfo.leadingCenter.data(), nominalChordInfo.trailingCenter.data(),
    //    nominalChordInfo.leadingVector.data(), nominalChordInfo.trailingVector.data()))
    //{
    //    //return false;//¡Ÿ ±◊¢ ÕµÙ
    //}
    ChordInformation measuredChordInfo;
    //if (!Chord(1, measuredChordInfo.leadingPoint.data(), measuredChordInfo.trailingPoint.data(),
    //    measuredChordInfo.leadingCenter.data(), measuredChordInfo.trailingCenter.data(),
    //    measuredChordInfo.leadingVector.data(), measuredChordInfo.trailingVector.data()))
    //{
    //    //return false;//¡Ÿ ±◊¢ ÕµÙ
    //}
    // create an initial guess
    Hexagon::Blade::FitOptions options;

    const Eigen::Isometry2d guessTransform =
        createGuessTransform(fp, nominalChordInfo, measuredChordInfo,
            Hexagon::Blade::measuredSectionCurve(this), LEType(), TEType());
   // Hexagon::Blade::FitOptions options;

    // construct inner and outer tolerance curves if applicable
    const ptrdiff_t numFineSamples = 16384;
    const Eigen::VectorXd fineNominalTValues =
        Eigen::VectorXd::LinSpaced(numFineSamples + 1, NomCurve()->t0(), NomCurve()->t1()).head(numFineSamples);
    const Eigen::Matrix2Xd fineNominalPoints = Hexagon::Blade::evaluate(*NomCurve(), fineNominalTValues);

    if (fp.fitToMiddleOfZone)
    {
        //const Eigen::Matrix2Xd nominalPoints = constructPointMatrix(m_nomx, m_nomy, m_numNomPoints);
        //auto coarseMinusTolerances = createToleranceCurve(*sectionCurve.whole, nominalPoints,
        //    Eigen::Map<const Eigen::ArrayXd>(m_mtol, m_numNomPoints));
        //auto coarsePlusTolerances = createToleranceCurve(*sectionCurve.whole, nominalPoints,
        //    Eigen::Map<const Eigen::ArrayXd>(m_ptol, m_numNomPoints));
        //Eigen::ArrayXd fineMinusTolerances = Hexagon::Blade::evaluate(*coarseMinusTolerances, fineNominalTValues);
        //Eigen::ArrayXd finePlusTolerances = Hexagon::Blade::evaluate(*coarsePlusTolerances, fineNominalTValues);

        //// if there are profile tolerances from the form dialog box, use those instead of the ones from the .NOM file
        //if (fp.profilePTol > fp.profileMTol)
        //{
        //    fineMinusTolerances = Eigen::ArrayXd::Constant(numFineSamples, fp.profileMTol);
        //    finePlusTolerances = Eigen::ArrayXd::Constant(numFineSamples, fp.profilePTol);
        //} // create the tolerance curves themselves
        //if (itMakesSenseToCreateInnerAndOuterToleranceCurves(fineMinusTolerances, finePlusTolerances))
        //{
        //    options.innerTolerance = createToleranceCurve(*sectionCurve.whole, fineNominalPoints, fineMinusTolerances);
        //    options.outerTolerance = createToleranceCurve(*sectionCurve.whole, fineNominalPoints, finePlusTolerances);
        //}
    }

    Eigen::VectorXd weightFittedPoints(m_totalPoints);
    std::vector<Hexagon::Blade::LinearDeviation> linearDeviations;
    std::unique_ptr<const Hexagon::Blade::Curve<2>> reducedCurveToFit;
   // figureOutWeightingAndEndpointConstraints(this);//≤‚ ‘figureOutWeightingAndEndpointConstraints

    options.weightFittedPoints = Eigen::VectorXd::Zero(m_totalPoints);
    for (int m = 0; m < m_totalPoints; m++)
    {
        options.weightFittedPoints[m] = 1;
    }


    const Hexagon::Blade::Curve<2>* curveToFit = NomCurve();
    if (reducedCurveToFit)//always empty
    {
        curveToFit = reducedCurveToFit.get();
    }
    // set the pivot points (may get overwritten later; that's OK)
    if (weightFittedPoints.sum() > 0.0)
    {
        options.pointsPivotPoint = measuredPoints * weightFittedPoints / weightFittedPoints.sum();
        options.targetCurvePivotPoint = guessTransform * options.pointsPivotPoint.head<2>();
    }
    else
    {
        options.pointsPivotPoint = Eigen::Vector2d::Zero();
        options.targetCurvePivotPoint = Eigen::Vector2d::Zero();
    }

    const Eigen::VectorXd distancesToPivot =
        (measuredPoints).colwise().norm().transpose();
    //const double scale = (weightFittedPoints.array() > 0.0).select(distancesToPivot, 0.0).maxCoeff();
    const double scale = 1000;
    options.translationXLimits = Eigen::Vector2d(-scale, scale);
    options.translationYLimits = Eigen::Vector2d(-scale, scale);

    // set the rotation and translation options
    setRotationOptions(options, fp);
    setTranslationOptions(options, fp, nominalChordInfo, measuredChordInfo);
    std::vector<Hexagon::Blade::LinearDeviation> linearDeviations1;

    if (fp.algorithm == BestFitAlgorithm::LeastSquares)
    {
       /* for (int m = 0; m < m_totalPoints; m++)
        {
            options.weightFittedPoints[m] = 1;
        }*/
        auto fitTransform = Hexagon::Blade::computeLeastSquaresBestFit(*curveToFit, measuredPoints, guessTransform, options,
            linearDeviations, inchSize);
        index = m_numBestFits;
        if (fp.tranfit == 4) // To X Axis
        {
            fitTransform(1, 2) = 0; // y offset
        }
        else if (fp.tranfit == 5) // To Y Axis
        {
            fitTransform(0, 2) = 0; // x offset
        }
        m_bestFits[m_numBestFits] = createMeasuredPointsToNominalCurveBestFit(
            sectionCurve, fitTransform, measuredPoints, options.weightFittedPoints.array() > 0.0, fp, mtols, ptols);
        for (int i = 0; i < m_totalPoints; i++)
        {
            m_bestFits[m_numBestFits]->Omega(i, options.weightFittedPoints[i]);
        }
        m_numBestFits++;
    }
    return true;
}
bool CSection::FitPointsV42(CFitParams& fp, int& index, double* mtols, double* ptols)
{
    // bugout(0, L"FitPoints for %s algorithm %d", m_name, fp.algorithm);
    // bugout(0, L"fitcurve %d %d %d %d", fp.fitcurve[0], fp.fitcurve[1], fp.fitcurve[2], fp.fitcurve[3]);
    // bugout(0, L"weightcurve %d %d %d %d", fp.weightcurve[0], fp.weightcurve[1], fp.weightcurve[2], fp.weightcurve[3]);

    CBestFit* bf = new CBestFit(m_totalPoints);

    double nose[2], range[2];

    bool checkNose = false;
    bool checkRange = false;
    // zzz
    if (fp.leoff2 > fp.leoff1 && NomPart(LEC) && NomPart(LEC)->Extreme() != -1.0)
    {
        NomPart(LEC)->CalcPoint(nose, NomPart(LEC)->Extreme());
        checkNose = true;
        // bugout(0, _T("10 %f %f N %s"), nose[0], nose[1], m_name);
        if (fp.fitcurve[CVC] && !fp.fitcurve[CCC])
        {
            double cvp0[2], cvp1[2], lep0[2];

            NomPart(CVC)->CalcPoint(cvp0, NomPart(CVC)->T0());
            NomPart(CVC)->CalcPoint(cvp1, NomPart(CVC)->T1());
            NomPart(LEC)->CalcPoint(lep0, NomPart(LEC)->T0());

            if (dist(cvp0, lep0) < 1.0e-3 || dist(cvp1, lep0) < 1.0e-3) // want to keep the beginning of the le
            {
                range[0] = NomPart(LEC)->T0();
                range[1] = NomPart(LEC)->Extreme();
            }
            else // want to keep end of the le
            {
                range[0] = NomPart(LEC)->Extreme();
                range[1] = NomPart(LEC)->T1();
            }
            checkRange = true;
        }
        else if (fp.fitcurve[CCC] && !fp.fitcurve[CVC])
        {
            double ccp0[2], ccp1[2], lep0[2];

            NomPart(CCC)->CalcPoint(ccp0, NomPart(CCC)->T0());
            NomPart(CCC)->CalcPoint(ccp1, NomPart(CCC)->T1());
            NomPart(LEC)->CalcPoint(lep0, NomPart(LEC)->T0());

            if (dist(ccp0, lep0) < 1.0e-3 || dist(ccp1, lep0) < 1.0e-3) // want to keep the beginning of the le
            {
                range[0] = NomPart(LEC)->T0();
                range[1] = NomPart(LEC)->Extreme();
            }
            else // want to keep end of the le
            {
                range[0] = NomPart(LEC)->Extreme();
                range[1] = NomPart(LEC)->T1();
            }
            checkRange = true;
        }
    }
    // zzz

    // bugout(0, L"Maxiters %d", fp.maxiters);
    int i, s;
    for (i = 0; i < m_totalPoints; i++) // assign meas and nom points to bestfit
    {
        s = m_partOf[i];
        bf->m_bestPartOf[i] = m_partOf[i];

        bf->PutVal(i, m_mxpt[i], m_mypt[i]);
        bf->PutNom(i, m_nxpt[i], m_nypt[i]);
        bf->PutT(i, m_nomt[i]);
        // bugout(0, L"9 %f %f %f %f FitPoints", m_mxpt[i], m_mypt[i], m_nxpt[i], m_nypt[i]);
        // bugout(0, L"7 %f %f NomPoint", m_nxpt[i], m_nypt[i]);
        // bugout(0, L"14 %f %f MeaPoint", m_mxpt[i], m_mypt[i]);

        bf->Omega(i, 0.0); // default is no weight to fit.

        if (s >= 0 && s <= 3)
        {
            if (fp.fitcurve[s] && m_nomPart[s])
            {
                bool usePoint = true;

                if (checkNose)
                {
                    double np[2];
                    np[0] = m_nxpt[i];
                    np[1] = m_nypt[i];
                    double d = dist(np, nose);
                    if (d < fp.leoff1 || d > fp.leoff2)
                        usePoint = false;
                    else if (s == LEC && checkRange) // may need to toss point if on LEC, but on side of nose point that isn't
                                                    // being kept
                    {
                        if (m_nomt[i] < range[0] || m_nomt[i] > range[1])
                            usePoint = false;
                    }
                }

                if (usePoint)
                {
                    // if (checkNose)
                    //  bugout(0, _T("7 %f %f %s"), m_nxpt[i], m_nypt[i], m_name);
                    bf->Omega(i, (double)fp.weightcurve[s]);
                }
            }
        }
    }

    int maxiters = fp.maxiters;
    if (maxiters < 1)
        maxiters = 1;
    if (maxiters > MAXITERS)
        maxiters = MAXITERS;

    if (fp.algorithm == BestFitAlgorithm::None) // no fit
    {
        for (i = 0; i < m_totalPoints; i++) // assign meas and nom points to bestfit
            bf->PutInf(i, m_mxpt[i], m_mypt[i]);

        RefindNomsV42(bf, fp, true, 0, NULL, mtols, ptols); // for form calculations

        maxiters = 0;
    }

    double nlcp[2], ntcp[2], nlctr[2], ntctr[2], nltv[2], nttv[2];
    double mlcp[2], mtcp[2], mlctr[2], mtctr[2], mltv[2], mttv[2];
    if (!Chord(0, nlcp, ntcp, nlctr, ntctr, nltv, nttv) || !Chord(1, mlcp, mtcp, mlctr, mtctr, mltv, mttv))
    {
        delete bf;
        return false;
    }

    double cl = dist(nlcp, ntcp);
    double clm = dist(mlcp, mtcp);

    CCurve* nommc = m_nomPart[MCC];

    bool good = true;

    bool usePivot = false;
    double nomPiv[2], meaPiv[2];
    if (fp.algorithm == BestFitAlgorithm::LeastSquares) // least squares
    {
        if (fp.pivot == 0 && fp.tranfit == 3) // le center
        {
            double r = 0.75 * cl;
            double npt[2], mpt[2];
            nommc->CircIntersect(nlcp, r, npt);
            m_meaPart[MCC]->CircIntersect(mlcp, r, mpt);
            if (bf->TwoPointFit(nlctr, npt, mlctr, mpt))
            {
                usePivot = true;
                nomPiv[0] = nlctr[0];
                nomPiv[1] = nlctr[1];
                meaPiv[0] = mlctr[0];
                meaPiv[1] = mlctr[1];
                RefindNomsV42(bf, fp, false);
            }
            else
                good = false;
        }
        else if (fp.pivot == 1 && fp.tranfit == 3) // le nose
        {
            double r = 0.75 * cl;
            double npt[2], mpt[2];
            nommc->CircIntersect(nlcp, r, npt);
            m_meaPart[MCC]->CircIntersect(mlcp, r, mpt);
            if (bf->TwoPointFit(nlcp, npt, mlcp, mpt))
            {
                usePivot = true;
                nomPiv[0] = nlcp[0];
                nomPiv[1] = nlcp[1];
                meaPiv[0] = mlcp[0];
                meaPiv[1] = mlcp[1];
                RefindNomsV42(bf, fp, false);
            }
            else
                good = false;
        }
        else if (fp.pivot == 2 && fp.tranfit == 3) // te center
        {
            double r = 0.75 * cl;
            double npt[2], mpt[2];
            nommc->CircIntersect(ntcp, r, npt);
            m_meaPart[MCC]->CircIntersect(mtcp, r, mpt);
            if (bf->TwoPointFit(ntctr, npt, mtctr, mpt))
            {
                usePivot = true;
                nomPiv[0] = ntctr[0];
                nomPiv[1] = ntctr[1];
                meaPiv[0] = mtctr[0];
                meaPiv[1] = mtctr[1];
                RefindNomsV42(bf, fp, false);
            }
            else
                good = false;
        }
        else if (fp.pivot == 3 && fp.tranfit == 3) // te tail
        {
            double r = 0.75 * cl;
            double npt[2], mpt[2];
            nommc->CircIntersect(ntcp, r, npt);
            m_meaPart[MCC]->CircIntersect(mtcp, r, mpt);
            if (bf->TwoPointFit(ntcp, npt, mtcp, mpt))
            {
                usePivot = true;
                nomPiv[0] = ntcp[0];
                nomPiv[1] = ntcp[1];
                meaPiv[0] = mtcp[0];
                meaPiv[1] = mtcp[1];
                RefindNomsV42(bf, fp, false);
            }
            else
                good = false;
        }
        else if (fp.fitcurve[LEC] == 1 && fp.fitcurve[TEC] == 0 && fp.fitcurve[CVC] == 0 && fp.fitcurve[CCC] == 0)
        {
            // LE only want to start with nose points aligned

            if (LEType() == EDGE_NORMAL)
            {
                double r = 0.75 * cl;
                double npt[2], mpt[2];
                nommc->CircIntersect(nlcp, r, npt);
                m_meaPart[MCC]->CircIntersect(mlcp, r, mpt);
                if (bf->TwoPointFit(nlcp, npt, mlcp, mpt))
                    RefindNomsV42(bf, fp, false);
                else
                    good = false;
            }
            else // square LE
            {
                double n0[2], n1[2], m0[2], m1[2];
                NomPart(LEC)->CalcPoint(n0, NomPart(LEC)->T0());
                NomPart(LEC)->CalcPoint(n1, NomPart(LEC)->T1());
                MeaPart(LEC)->CalcPoint(m0, MeaPart(LEC)->T0());
                MeaPart(LEC)->CalcPoint(m1, MeaPart(LEC)->T1());
                if (bf->TwoPointFit(n0, n1, m0, m1))
                    RefindNomsV42(bf, fp, false);
                else
                    good = false;
            }
        }
        else if (fp.fitcurve[LEC] == 0 && fp.fitcurve[TEC] == 1 && fp.fitcurve[CVC] == 0 && fp.fitcurve[CCC] == 0)
        {
            // TE only want to start with tail points aligned
            if (TEType() == EDGE_NORMAL)
            {
                double r = 0.75 * cl;
                double npt[2], mpt[2];
                nommc->CircIntersect(ntcp, r, npt);
                m_meaPart[MCC]->CircIntersect(mtcp, r, mpt);
                if (bf->TwoPointFit(ntcp, npt, mtcp, mpt))
                    RefindNomsV42(bf, fp, false);
                else
                    good = false;
            }
            else
            {
                double n0[2], n1[2], m0[2], m1[2];
                NomPart(TEC)->CalcPoint(n0, NomPart(TEC)->T0());
                NomPart(TEC)->CalcPoint(n1, NomPart(TEC)->T1());
                MeaPart(TEC)->CalcPoint(m0, MeaPart(TEC)->T0());
                MeaPart(TEC)->CalcPoint(m1, MeaPart(TEC)->T1());
                if (bf->TwoPointFit(n0, n1, m0, m1))
                    RefindNomsV42(bf, fp, false);
                else
                    good = false;
            }
        }
    }

    // if(fp.algorithm == 2) // guillotine
    //{
    //  // guillotine only operates on CV or CC side
    //  CCurve* gcurve = 0;

    //  fp.fitcurve[LEC] = fp.fitcurve[TEC] = 0;

    //  if(fp.fitcurve[CVC])
    //    gcurve = m_nomPart[CVC];
    //  else if(fp.fitcurve[CCC])
    //    gcurve = m_nomPart[CCC];
    //  else
    //    good = false;

    //  if(good)
    //  {
    //    if(bf->GuillotineFit(gcurve, fp))
    //      RefindNoms(bf, fp, true, 0, NULL, mtols, ptols);
    //    else
    //      good = false;
    //  }

    //  maxiters = 0;
    //}

    if (fp.algorithm == BestFitAlgorithm::TwoPointsOnMCLFromNose) // two points on MCL - from nose
    {
        double r = 0.01 * fp.lepercent * cl;
        double onom[2], pnom[2], omea[2], pmea[2];

        if (!nommc->CircIntersect(nlcp, r, onom))
            int_circ_line(nlcp, r, nlcp, nlctr, onom);
        if (!m_meaPart[MCC]->CircIntersect(mlcp, r, omea))
            int_circ_line(mlcp, r, mlcp, mlctr, omea);

        r = 0.01 * fp.tepercent * clm;
        if (!nommc->CircIntersect(nlcp, r, pnom))
            int_circ_line(nlcp, r, ntctr, ntcp, pnom);
        if (!m_meaPart[MCC]->CircIntersect(mlcp, r, pmea))
        {
            int_circ_line(mlcp, r, mtctr, mtcp, pmea);
        }
        if (bf->TwoPointFit(onom, pnom, omea, pmea))
            RefindNomsV42(bf, fp, true, 0, NULL, mtols, ptols);
        else
            good = false;
        // bugout(0, _T("10 %f %f NL"), onom[0], onom[1]);
        // bugout(0, _T("10 %f %f NT"), pnom[0], pnom[1]);
        double xy[2];
        bf->GetAlign()->MeasToBest(omea, 1, xy);
        // bugout(0, _T("10 %f %f ML"), xy[0], xy[1]);
        bf->GetAlign()->MeasToBest(pmea, 1, xy);
        // bugout(0, _T("10 %f %f MT"), xy[0], xy[1]);

        maxiters = 0;
    }

    if (fp.algorithm == BestFitAlgorithm::TwoPointsOnMCLFromTail) // two points on MCL - from tail
    {
        double r = 0.01 * fp.tepercent * cl;
        double onom[2], pnom[2], omea[2], pmea[2];

        if (!nommc->CircIntersect(ntcp, r, onom))
            int_circ_line(ntcp, r, ntcp, ntctr, onom);
        // bugout(0, _T("10 %lf %lf NT"), onom[0], onom[1]);

        if (!m_meaPart[MCC]->CircIntersect(mtcp, r, omea))
            int_circ_line(mtcp, r, mtcp, mtctr, omea);
        // bugout(0, _T("10 %lf %lf MT"), omea[0], omea[1]);

        r = 0.01 * fp.lepercent * clm;
        if (!nommc->CircIntersect(ntcp, r, pnom))
            int_circ_line(ntcp, r, nlctr, nlcp, pnom);
        // bugout(0, _T("10 %lf %lf NN"), pnom[0], pnom[1]);

        if (!m_meaPart[MCC]->CircIntersect(mtcp, r, pmea))
            int_circ_line(mtcp, r, mlctr, mlcp, pmea);
        // bugout(0, _T("10 %lf %lf MN"), pmea[0], pmea[1]);

        if (bf->TwoPointFit(onom, pnom, omea, pmea))
            RefindNomsV42(bf, fp, true, 0, NULL, mtols, ptols);
        else
            good = false;

        maxiters = 0;
    }
    CMinMax* mmf = NULL;

    if (fp.algorithm == BestFitAlgorithm::MinMax) // min max fit
    {
        fp.tranfit = fp.rotfit = -1; // who knows why this makes it happy.

        if (bf->LeastSquaresFit(fp))
            RefindNomsV42(bf, fp, false);

        for (i = 0; i < m_totalPoints; i++) // assign meas and nom points to bestfit
            bf->PutVec(i, m_ival[i], m_jval[i]);

        mmf = new CMinMax(m_totalPoints);
    }

    if (fp.algorithm == BestFitAlgorithm::TwoPointsOnMCLForForgedBlade) // two points on MCL for forged blade
    {
        double lv[2]; // line vector in forge plane
        double nv[2]; // normal vector of vorge plane
        double ang = fp.forgeAngle * M_PI / 180;
        lv[0] = cos(ang);
        lv[1] = sin(ang);
        nv[0] = lv[1];
        nv[1] = -lv[0];

        double rle = 0.01 * fp.lepercent * cl;
        double rte = 0.01 * fp.tepercent * cl;

        double tmid, orig[2], nple[2], npte[2], mple[2], mpte[2];
        orig[0] = orig[1] = 0.0;

        // find le point on nom mcl
        if (!nommc->CircIntersect(nlcp, rle, nple))
            int_circ_line(nlcp, rle, nlcp, nlctr, nple);
        double srle = _hypot(nple[0], nple[1]);
        // bugout(0, _T("10 %f %f NL"), nple[0], nple[1]);

        // find te point on nom mcl
        if (!nommc->CircIntersect(nlcp, rte, npte))
            int_circ_line(nlcp, rte, ntctr, ntcp, npte);
        double srte = _hypot(npte[0], npte[1]);
        // bugout(0, _T("10 %f %f NT"), npte[0], npte[1]);

        double norig[2], morig[2], ndir[2], mdir[2], nother[2], mother[2];

        ndir[0] = npte[0] - nple[0];
        ndir[1] = npte[1] - nple[1];
        normalize(ndir, ndir);
        int_line_line(orig[0], orig[1], nv[0], nv[1], nple[0], nple[1], ndir[0], ndir[1], &norig[0], &norig[1]);

        double deltaLE = dist(norig, nple);
        double deltaTE = dist(norig, npte);

        m_meaPart[MCC]->ClosestPoint(orig, lv, &tmid); // find approximate location of origin on mcl

        if (!m_meaPart[MCC]->CircIntersect(orig, srle, mple, m_meaPart[MCC]->T0(), tmid))
            int_circ_line(orig, srle, mlctr, mlcp, mple);
        // bugout(0, _T("10 %f %f ML"), mple[0], mple[1]);

        if (!m_meaPart[MCC]->CircIntersect(orig, srte, mpte, tmid, m_meaPart[MCC]->T1()))
            int_circ_line(orig, srte, mtctr, mtcp, mpte);
        // bugout(0, _T("10 %f %f MT"), mpte[0], mpte[1]);

        mdir[0] = mpte[0] - mple[0];
        mdir[1] = mpte[1] - mple[1];
        normalize(mdir, mdir);
        int_line_line(orig[0], orig[1], nv[0], nv[1], mple[0], mple[1], mdir[0], mdir[1], &morig[0], &morig[1]);

        if (deltaLE > deltaTE)
        {
            nother[0] = nple[0];
            nother[1] = nple[1];
            mother[0] = morig[0] - deltaLE * mdir[0];
            mother[1] = morig[1] - deltaLE * mdir[1];
        }
        else
        {
            nother[0] = npte[0];
            nother[1] = npte[1];
            mother[0] = morig[0] + deltaTE * mdir[0];
            mother[1] = morig[1] + deltaTE * mdir[1];
        }

        // bugout(0, _T("10 %f %f N"), norig[0], norig[1]);
        // bugout(0, _T("10 %f %f N"), nother[0], nother[1]);
        // bugout(0, _T("10 %f %f M"), morig[0], morig[1]);
        // bugout(0, _T("10 %f %f M"), mother[0], mother[1]);

        if (bf->TwoPointFit(norig, nother, morig, mother))
            RefindNomsV42(bf, fp, true, 0, NULL, mtols, ptols);
        else
            good = false;

        maxiters = 0;
    }
    if (!good)
    {
        delete bf;
        return false;
    }

    double lsq;

    CMatrix* lastInf = new CMatrix(m_totalPoints, 2);

    int iter;
    for (iter = 0; iter < maxiters; iter++)
    {
        if (iter > 0)
        {
            for (i = 0; i < m_totalPoints; i++)
            {
                lastInf->m[i][0] = bf->m_infs->m[i][0];
                lastInf->m[i][1] = bf->m_infs->m[i][1];
            }
        }
        /* FIX ME LATER

        if (fp.algorithm == 5 && vf)
        {
        vf->TransferIn(bf->m_noms, bf->m_vals, bf->m_ijks, bf->m_omega, &bf->m_align);
        good = vf->vect_fit();
        vf->TransferOut(bf->m_infs, &bf->m_align);
        }
        else
        */
        if (fp.algorithm == BestFitAlgorithm::MinMax && mmf)
        {
            mmf->TransferIn(bf->m_noms, bf->m_vals, bf->m_omega, &bf->m_align);
            good = mmf->MinMaxFit();
            // bugout(0, _T("MinMaxFit iteration %d"), iter);
            mmf->TransferOut(bf->m_infs, &bf->m_align);
        }
        else
        {
            if (usePivot)
                good = bf->LeastSquaresFit(fp, nomPiv, meaPiv);
            else
                good = bf->LeastSquaresFit(fp);
        }

        if (!good)
            break;

        bool breakOut = false;
        if (iter == maxiters - 1) // need to check alignment somehow too
        {
            breakOut = true;
        }

        if (iter > 0)
        {
            // compute max distance distance between vals and infs
            double maxDist = 0.0;
            for (i = 0; i < m_totalPoints; i++)
            {
                double last[2], inf[2];
                inf[0] = bf->m_infs->m[i][0];
                inf[1] = bf->m_infs->m[i][1];
                last[0] = lastInf->m[i][0];
                last[1] = lastInf->m[i][1];
                double d = dist(last, inf);
                if (d > maxDist)
                    maxDist = d;

                // if (iter == 0)
                //  bugout(0, _T("8 %lf %lf %lf %lf %lf"), last[0], last[1], inf[0], inf[1], d);
                // else if (iter == 9)
                //  bugout(0, _T("9 %lf %lf %lf %lf %lf"), last[0], last[1], inf[0], inf[1], d);
            }
            // bugout(0, _T("iter %d maxDist %lf stopDist %lf"), iter, maxDist, fp.stopDist);

            if (maxDist < 0) //42memo
                breakOut = true;
        }
        // test maxDist and break out if it is sufficiently small.

        RefindNomsV42(bf, fp, breakOut, 0, &lsq, mtols, ptols);
        if (breakOut)
            break;
    }

    delete lastInf;

    index = m_numBestFits;
    m_bestFits[m_numBestFits] = bf;
    m_numBestFits++;

    if (mmf)
        delete mmf;

    // bugout(0, _T("Fit Measured Points %s (%d) LSQ = %lf"), m_name, index, lsq);

    /* FIX ME LATER
    if (vf)
    delete vf;
    */

    return true;
}

bool CSection::RefindNomsV42(CBestFit* bf, CFitParams& fp, bool finalTime, int offset, double* lsq, double* mtols, double* ptols, bool useRanges)
{
    if (!bf || !m_meaCurve || !m_meaPart[CVC] || !m_meaPart[CCC] || !m_meaPart[LEC] || !m_meaPart[TEC])
        return false;

    if (!m_nomCurve)
        return false;

    // if (finalTime) bugout(0, _T("FitCurves %d %d %d %d"), fp.fitcurve[0], fp.fitcurve[1], fp.fitcurve[2],
    // fp.fitcurve[3]);

    int s;
    double t0[4], t1[4], nt0[4], nt1[4];
    for (s = 0; s < 4; s++)
    {
        t0[s] = m_meaPart[s]->T0();
        t1[s] = m_meaPart[s]->T1();
        nt0[s] = m_nomPart[s]->T0();
        nt1[s] = m_nomPart[s]->T1();
        bf->m_mindev[s] = 1.0e20;
        bf->m_maxdev[s] = -1.0e20;
        bf->m_meandev[s] = bf->m_stddev[s] = 0.0;
    }

    double period = m_meaCurve->Period();
    double nperiod = m_nomCurve->Period();
    m_meaCurve->Align(bf->GetAlign());
    double sum[4], sumsq[4];
    int num[4];
    int totalBad = 0, totalChecked = 0;
    sum[0] = 0.0;
    sumsq[0] = 0.0;
    num[0] = 0;
    sum[1] = 0.0;
    sumsq[1] = 0.0;
    num[1] = 0;
    sum[2] = 0.0;
    sumsq[2] = 0.0;
    num[2] = 0;
    sum[3] = 0.0;
    sumsq[3] = 0.0;
    num[3] = 0;

    double bnt = 0.0, bmt = 0.0;
    double lastd = 2000.0;
    bool kickstart = false;

    double totalsumsq = 0.0;
    int numChecked = 0;
    for (int j = 0; j < m_totalPoints; j++)
    {
        double best[2], nom[2], nv[2], ntv[2];

        bf->GetInf(j + offset, best);

        // This is much quicker if each point can seed the next point (seed = -1).
        // But this doesn't work for partial sections.
        // Fully checking every point is slow, but don't know where the break is at this point.
        // First hint that something is wrong, is if point does not project onto measured curve.
        // This distance should be zero.

        int mseed = j == 0 ? 400 : -1; // exhaustive search for first point only.
        int nnseed = j == 0 ? 400 : -1;
        int nbseed = j == 0 ? 400 : -1;

        double d = m_meaCurve->ClosestPoint(best, nom, &bmt, 0, 0., 0., mseed);
        if (d > 0.0001)
        {
            mseed = 400; // found wrong point, perform more exhaustive search
            d = m_meaCurve->ClosestPoint(best, nom, &bmt, 0, 0., 0., mseed);
        }

        // this is probably useless, we set below based on where the point projects nominally...
        for (s = 0; s < 4; s++)
            if ((bmt >= t0[s] && bmt < t1[s]) || (bmt + period >= t0[s] && bmt + period < t1[s]))
                break;

        nv[0] = m_ival[j];
        nv[1] = m_jval[j];

        normalize(nv, nv);
        int usenorm = ((fabs(nv[0]) <= 1.0) && (fabs(nv[1]) <= 1.0) && (l2norm(nv) > 0.0));
        d = 2000.0;
        // bugout(0, _T("A bnt %f"), bnt);

        if (usenorm)
        {
            if (kickstart)
                nnseed = 400;

            // d = m_nomCurve->ClosestNominal(best, nv, nom, &bnt, 0, 0., 0., j == 0 ? 400 : -1);
            d = m_nomCurve->ClosestNominal(best, nv, nom, &bnt, ntv, 0., 0., nnseed);
            // if (finalTime) bugout(0, _T("d=%f lastd=%f %s j=%d"), d, lastd, m_name, j);
            if (lastd < 1000.0 && lastd > 1.0e-6 && d / lastd > 50)
                d += 2000.0;
            // bugout(0, _T("A bnt %f (d=%f)"), bnt, d);
        }

        if (d > 1000.0) // not usenorm of ClosestNominal failed
        {
            // d = m_nomCurve->ClosestPoint(best, nom, &bnt, 0, 0., 0., j == 0 ? 400 : -1);
            if (usenorm)
            {
                kickstart = true;
                nbseed = 600;
            }
            d = m_nomCurve->ClosestPoint(best, nom, &bnt, ntv, 0., 0., nbseed);
        }

        lastd = d;

        double dist, ijk[2];
        ijk[0] = best[0] - nom[0];
        ijk[1] = best[1] - nom[1];
        normalize(ijk, ijk);
        normalize(ntv, ntv);

        int nns;
        for (nns = 0; nns < 4; nns++)
            if ((bnt >= nt0[nns] && bnt < nt1[nns]) || (bnt + nperiod >= nt0[nns] && bnt + nperiod < nt1[nns]))
                break;

        if (nns < 4)
            s = nns;

        // if this is a partial section and we are close to the end missing, we should do closest point
        if (m_leType == EDGE_PARTIAL)
        {
            if (bnt > nt0[LEC] - 0.04 * nperiod && bnt < nt1[LEC] + 0.04 * nperiod)
            {
                if (fabs(dot(ijk, ntv)) > 0.05) // more than about 3 degrees off
                    s = LEC;
            }
        }

        if (m_teType == EDGE_PARTIAL)
        {
            if (bnt > nt0[TEC] - 0.04 * nperiod && bnt < nt1[TEC] + 0.04 * nperiod)
            {
                if (fabs(dot(ijk, ntv)) > 0.05) // more than about 3 degrees off
                    s = TEC;
            }
        }

        if (s < 4)
            bf->m_bestPartOf[j + offset] = s;
        else
            s = bf->m_bestPartOf[j + offset]; // shouldn't happen

        if (bnt > m_nomCurve->Period()) // need to get this back in the domain, for next loop
            bnt -= m_nomCurve->Period();
        if (bnt < 0.0)
            bnt += m_nomCurve->Period();

        double oldnom[2], shift[2];
        bf->GetNom(j + offset, oldnom);
        shift[0] = nom[0] - oldnom[0];
        shift[1] = nom[1] - oldnom[1];

        // double sd = dot(shift, shift);
        // if (sd > maxshift)
        // maxshift = sd;
        // sumshift[0] += shift[0];
        // sumshift[1] += shift[1];

        bf->PutNom(j + offset, nom[0], nom[1]);
        bf->PutT(j + offset, bnt);
        // don't change Omega values when doing LE Arc LS fit

        if (fp.leoff1 >= fp.leoff2)
            bf->Omega(j + offset, 0.0); // default is no weight to fit.

        if (fp.algorithm == BestFitAlgorithm::MinMax && s == LEC && fp.fitcurve[LEC] && !fp.fitcurve[TEC] && !fp.fitcurve[CVC] &&
            !fp.fitcurve[CCC]) // LE Only
        {
        }
        else if (fp.algorithm == BestFitAlgorithm::MinMax && s == TEC && fp.fitcurve[TEC] && !fp.fitcurve[LEC] &&
            !fp.fitcurve[CVC] &&
            !fp.fitcurve[CCC]) // TE Only
        {
        }

        if (fp.leoff1 >= fp.leoff2)
            if (s >= 0 && s <= 3)
                if (fp.fitcurve[s])
                    bf->Omega(j + offset, (double)fp.weightcurve[s]);

        /*  Try weighting points based on projected deviation.
        if (!finalTime)  // don't need to do this until bestfit is completed
        continue;
        */

        // for full blade LE arc fit, set weight to zero if nominal parameter isn't in the range
        if (useRanges)
        {
            if (!ParameterInRange(bnt, m_arcRangeCV[0], m_arcRangeCV[1], m_nomCurve->Period()) &&
                !ParameterInRange(bnt, m_arcRangeCC[0], m_arcRangeCC[1], m_nomCurve->Period()))
                bf->Omega(j + offset, 0.0);
            else
            {
                // bugout(0, L"14 %lf %lf %s j=%d", nom[0], nom[1], m_name, j);
            }
        }

        int flip = 1;

        double dummy[2], tp[2];
        m_nomCurve->CalcPoint(dummy, /*m_nomt[j]*/ bnt, tp);

        if (curl(ijk, tp) < 0.0)
            flip *= -1;

        ijk[0] *= flip;
        ijk[1] *= flip;

        dist = (best[0] - nom[0]) * ijk[0] + (best[1] - nom[1]) * ijk[1];
        // bugout(0, _T("dist %f"), dist);

        bf->PutVec(j + offset, ijk[0], ijk[1]);

        if (dist < bf->m_mindev[s]) // for form calculations
            bf->m_mindev[s] = dist;
        if (dist > bf->m_maxdev[s])
            bf->m_maxdev[s] = dist;

        num[s]++;
        sum[s] += dist;
        sumsq[s] += dist * dist;
        totalsumsq += dist * dist; // for debug

        if (mtols && ptols)
        {
            if (dist < mtols[s] || dist > ptols[s])
                totalBad++;
            totalChecked++;
        }

        if (bf->Omega(j + offset) <= 0.0)
            continue;
        numChecked++;
        // if (finalTime) bugout(2, _T("point %3d dist %.4f dist^2 %.8f"), j, dist, dist*dist);
    }

    if (lsq)
        *lsq = totalsumsq;

    bf->m_totalBad = totalBad;
    bf->m_totalChecked = totalChecked;

    // double sss = 0.0;
    for (s = 0; s < 4; s++)
    {
        if (num[s] > 0)
            bf->m_meandev[s] = sum[s] / num[s];

        if (num[s] > 1)
        {
            // bugout(0, _T(" sum sq (%d) %f %s"), s, sumsq[s], m_name);
            // sss += sumsq[s];
            bf->m_stddev[s] = (sumsq[s] - (sum[s] * sum[s]) / (double)num[s]) / (num[s] - 1.0);
            bf->m_stddev[s] = sqrt(bf->m_stddev[s]);
        }
    }
    m_meaCurve->Align(NULL);

    return true;
}

bool CSection::RefindMeasV42(CBestFit* bf, CFitParams& fp, bool finalTime, bool firsttime, int offset, double* lsq)
{
    return false;
}

bool CSection::AssignPoints(double* xv, double* yv, int n, int* /*start*/, int* /*end*/)
{
    // this is for closed curve for P&W AS file read and analysis file
    bugout(0, L"CSection::AssignPoints: m_totalPoints %d", n);

    m_totalPoints = n;
    m_cxpt = new double[m_totalPoints];
    m_cypt = new double[m_totalPoints];
    m_czpt = new double[m_totalPoints];
    m_ival = new double[m_totalPoints];
    m_jval = new double[m_totalPoints];
    m_kval = new double[m_totalPoints];
    m_mxpt = new double[m_totalPoints];
    m_mypt = new double[m_totalPoints];
    m_nxpt = new double[m_totalPoints];
    m_nypt = new double[m_totalPoints];
    m_nomt = new double[m_totalPoints];
    m_partOf = new int[m_totalPoints];
    if (!m_cxpt || !m_cypt || !m_czpt || !m_mxpt || !m_mypt || !m_ival || !m_jval || !m_kval || !m_nxpt || !m_nypt ||
        !m_nomt)
        return false;
    // make sure the ijk vectors are initialized to zero
    Eigen::Map<Eigen::VectorXd>(m_ival, m_totalPoints).setZero();
    Eigen::Map<Eigen::VectorXd>(m_jval, m_totalPoints).setZero();
    Eigen::Map<Eigen::VectorXd>(m_kval, m_totalPoints).setZero();
    CAlignment* bfa = 0;
    double nlcp[2], ntcp[2], nlctr[2], ntctr[2], nltv[2], nttv[2];
    if ((m_leType != EDGE_PARTIAL || m_teType != EDGE_PARTIAL) && // if both ends are partial, don't try to improve
        Chord(0, nlcp, ntcp, nlctr, ntctr, nltv, nttv))
    {
        double mlcp[2], mtcp[2], mlctr[2], mtctr[2], mltv[2], mttv[2];
        if (Chord(1, mlcp, mtcp, mlctr, mtctr, mltv, mttv))
        {
            CCurve* nommc = m_nomPart[MCC];

            double r, cl, nle[2], mle[2], nte[2], mte[2], norig[2], morig[2];

            cl = _hypot(nlcp[0] - ntcp[0], nlcp[1] - ntcp[1]);

            r = 0.05 * cl;
            nommc->CircIntersect(nlcp, r, nle);
            if (dist(nlcp, nle) - r > 1.0e-4)
            {
                double ij[2];
                ij[0] = nle[0] - nlcp[0];
                ij[1] = nle[1] - nlcp[1];
                normalize(ij, ij);
                nle[0] = nlcp[0] + r * ij[0];
                nle[1] = nlcp[1] + r * ij[1];
            }

            m_meaPart[MCC]->CircIntersect(mlcp, r, mle);
            if (dist(mlcp, mle) - r > 1.0e-4)
            {
                double ij[2];
                ij[0] = mle[0] - mlcp[0];
                ij[1] = mle[1] - mlcp[1];
                normalize(ij, ij);
                mle[0] = mlcp[0] + r * ij[0];
                mle[1] = mlcp[1] + r * ij[1];
            }

            nommc->CircIntersect(ntcp, r, nte);
            if (dist(ntcp, nte) - r > 1.0e-4)
            {
                double ij[2];
                ij[0] = nte[0] - ntcp[0];
                ij[1] = nte[1] - ntcp[1];
                normalize(ij, ij);
                nte[0] = ntcp[0] + r * ij[0];
                nte[1] = ntcp[1] + r * ij[1];
            }

            m_meaPart[MCC]->CircIntersect(mtcp, r, mte);
            if (dist(mtcp, mte) - r > 1.0e-4)
            {
                double ij[2];
                ij[0] = mte[0] - mtcp[0];
                ij[1] = mte[1] - mtcp[1];
                normalize(ij, ij);
                mte[0] = mtcp[0] + r * ij[0];
                mte[1] = mtcp[1] + r * ij[1];
            }

            if (m_leType == EDGE_PARTIAL)
            {
                norig[0] = nte[0];
                norig[1] = nte[1];
                morig[0] = mte[0];
                morig[1] = mte[1];

                // modify nle and mle
                nommc->CircIntersect(ntcp, 0.25 * cl, nle);
                m_meaPart[MCC]->CircIntersect(mtcp, 0.25 * cl, mle);
            }
            else if (m_teType == EDGE_PARTIAL)
            {
                norig[0] = nle[0];
                norig[1] = nle[1];
                morig[0] = mle[0];
                morig[1] = mle[1];

                // modify nte and mte
                nommc->CircIntersect(nlcp, 0.25 * cl, nte);
                m_meaPart[MCC]->CircIntersect(mlcp, 0.25 * cl, mte);
            }
            else
            {
                norig[0] = 0.5 * (nle[0] + nte[0]);
                norig[1] = 0.5 * (nle[1] + nte[1]);
                morig[0] = 0.5 * (mle[0] + mte[0]);
                morig[1] = 0.5 * (mle[1] + mte[1]);
            }

            // Note, this isn't equivalent to TwoPointFit, the origins are (usually) set to the midpoints.

            //CBestFit bf(2);

            //bf.PutNom(0, nle[0], nle[1]);
            //bf.PutNom(1, nte[0], nte[1]);
            //bf.PutVal(0, mle[0], mle[1]);
            //bf.PutVal(1, mte[0], mte[1]);
            //bf.Omega();
            //CFitParams fp;
            //fp.tranfit = 0;
            //fp.rotfit = 0;
            //fp.algorithm = BestFitAlgorithm::LeastSquares;

            //bf.LeastSquaresFit(fp, norig, morig);

            //bfa = new CAlignment;
            //bf.Align(bfa);
        }
    }

    Eigen::Vector4d t0, t1;
    /*for (int s = 0; s < 4; s++)
    {
        t0[s] = m_meaPart[s]->T0();
        t1[s] = m_meaPart[s]->T1();
    }*/
    // polygonalize the measured curve
    auto measuredPolygon = Blade::polygonalizeWithT(*m_meaCurve, 2048, 1e-4);
    Eigen::MatrixX2d measuredPoints = std::get<0>(measuredPolygon).transpose();
    Eigen::VectorXd measuredT = std::get<1>(measuredPolygon);

    // polygonalize the nominal curve
    auto nominalPolygon = Blade::polygonalizeWithT(*m_nomCurve, 2048, 1e-4);
    Eigen::MatrixX2d nominalPoints = std::get<0>(nominalPolygon).transpose();
    Eigen::VectorXd nominalT = std::get<1>(nominalPolygon);

    // construct KD trees of the measured and nominal curves
    //typedef nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixX2d, 2> KDTree;
    //KDTree nominalKDtree(2 /*dimensions*/, nominalPoints, 10 /*max leaf*/);
    //KDTree measuredKDtree(2 /*dimensions*/, measuredPoints, 10 /*max leaf*/);

    //// map arrays of the inputs
    //Eigen::Matrix2Xd points(2, m_totalPoints);
    //points.row(0) = Eigen::Map<const Eigen::VectorXd>(xv, m_totalPoints).transpose();
    //points.row(1) = Eigen::Map<const Eigen::VectorXd>(yv, m_totalPoints).transpose();



    for (int i = 0; i < m_totalPoints; i++)
    {
        double xyz[2], bxy[2], tanv[2], nv[2];
        xyz[0] = m_mxpt[i] = xv[i];
        xyz[1] = m_mypt[i] = yv[i];
        int seed = i == 0 ? 400 : -1; // full search for 1st, point otherwise quick search
    }
    for (int i = 0; i < m_totalPoints; i++)
    {
        m_partOf[i] = 0;
    }
    return true;
}
int CSection::Chord(int flg, double* lcp, double* tcp, double* lctr, double* tctr, double* ltv, double* ttv, double* m,
    double* w, double* tew, double* zeroPt)
{
    CCurve* lec, * tec, * cvc, * ccc, * mcc;
    double mclLERatio, mclTERatio;
    if (flg == 0)
    {
        lec = m_nomPart[LEC];
        tec = m_nomPart[TEC];
        cvc = m_nomPart[CVC];
        ccc = m_nomPart[CCC];
        mcc = m_nomPart[MCC];

        mclLERatio = 0.001;// myGetProfileDouble(L"NomMclLEBackoff", 0.001); // NEED TO FIX THESE
        mclTERatio = 0.001;//myGetProfileDouble(L"NomMclTEBackoff", 0.001);
    }
    else
    {
        lec = m_meaPart[LEC];
        tec = m_meaPart[TEC];
        cvc = m_meaPart[CVC];
        ccc = m_meaPart[CCC];
        mcc = m_meaPart[MCC];

        mclLERatio = 0.004;// myGetProfileDouble(L"MeaMclLEBackoff", 0.004); // NEED TO FIX THESE
        mclTERatio = 0.004;//myGetProfileDouble(L"MeaMclTEBackoff", 0.004);
    }
    if (!lec || !tec || !cvc || !ccc || !mcc)
        return 0;
    double p0[2], p1[2], t0[2], t1[2], lep[2], tep[2];

    // determine which end of mean camber goes with which edge

    double deltaLE = mclLERatio * (mcc->T1() - mcc->T0()); // avoid hook on end
    double deltaTE = mclTERatio * (mcc->T1() - mcc->T0()); // avoid hook on end

    mcc->CalcPoint(p0, mcc->T0() + deltaLE, t0);
    mcc->CalcPoint(p0, mcc->T0());
    mcc->CalcPoint(p1, mcc->T1() - deltaTE, t1);
    mcc->CalcPoint(p1, mcc->T1());
    lec->CalcPoint(lep, lec->T0());
    tec->CalcPoint(tep, tec->T0());
    normalize(t0, t0);
    normalize(t1, t1);
    if (_hypot(lep[0] - p0[0], lep[1] - p0[1]) < _hypot(tep[0] - p0[0], tep[1] - p0[1]))
    {
        // t0 and p0 go with le
        lctr[0] = p0[0];
        lctr[1] = p0[1];
        ltv[0] = -t0[0]; // point away
        ltv[1] = -t0[1];
        tctr[0] = p1[0];
        tctr[1] = p1[1];
        ttv[0] = t1[0];
        ttv[1] = t1[1];
    }
    else
    {
        // t0 and p0 go with te
        tctr[0] = p0[0];
        tctr[1] = p0[1];
        ttv[0] = -t0[0];
        ttv[1] = -t0[1];
        lctr[0] = p1[0];
        lctr[1] = p1[1];
        ltv[0] = t1[0];
        ltv[1] = t1[1];
    }

    // make sure vectors point away from center

    if ((lctr[0] - tctr[0]) * ltv[0] + (lctr[1] - tctr[1]) * ltv[1] < -0.5) // changed these from 0.0 to -0.5 for big hook
                                                                           // blades
    {
        ltv[0] *= -1.0;
        ltv[1] *= -1.0;
    }

    if ((tctr[0] - lctr[0]) * ttv[0] + (tctr[1] - lctr[1]) * ttv[1] < -0.5)
    {
        ttv[0] *= -1.0;
        ttv[1] *= -1.0;
    }

    if (lec->Extreme() != -1.0 && tec->Extreme() != -1.0)
    {
        lec->CalcPoint(lcp, lec->Extreme());
        tec->CalcPoint(tcp, tec->Extreme());
    }
    else
    {
        if (lec->Type() == POINT_TYPE || !lec->LineIntersect(lctr, ltv, lcp))
        {
            lcp[0] = lep[0];
            lcp[1] = lep[1];
        }

        if (tec->Type() == POINT_TYPE || !tec->LineIntersect(tctr, ttv, tcp))
        {
            tcp[0] = tep[0];
            tcp[1] = tep[1];
        }
    }

    if (m || w || tew)
    {
        // bugout(0, _T("Chord %s %d"), m_name, flg);
        // bugout(0, _T("10 %f %f Ns%d"), lcp[0], lcp[1], flg);
        // bugout(0, _T("10 %f %f Tl%d"), tcp[0], tcp[1], flg);
        double et, cv[2], ltpt[2], ttpt[2];
        cv[0] = tcp[0] - lcp[0];
        cv[1] = tcp[1] - lcp[1];
        normalize(cv, cv);

        double bv[2];
        bv[0] = -cv[0];
        bv[1] = -cv[1];
        // -1 is a kluge, see curve.cpp

        if (m || w)
        {
            if (!lec->Extreme(bv, &et, ltpt, 0.0, 0.0, -1))
                return 0;
            // bugout(0, _T("10 %f %f Le%d %s"), ltpt[0], ltpt[1], flg, m_name);
        }

        if (m || tew)
        {
            if (!tec->Extreme(cv, &et, ttpt, 0.0, 0.0, -1))
                return 0;
        }

        if (m)
        {
            double projpt[2]; // just for debug
            *m = fabs(projdist(ltpt, cv, ttpt, projpt));
            // bugout(0, _T("10 %f %f Te%d %s dist %f"), ttpt[0], ttpt[1], flg, m_name, *m);
            // bugout(0, _T("10 %f %f PJ"), projpt[0], projpt[1]);
        }

        double origin[2] = { 0.0, 0.0 };
        if (zeroPt != NULL)
        {
            origin[0] = zeroPt[0];
            origin[1] = zeroPt[1];
        }
        else
            origin[0] = origin[1] = 0.0;

        if (w)
        {
            *w = projdist(ltpt, cv, origin);
        }

        if (tew)
        {
            *tew = projdist(ttpt, bv, origin);
        }
    }

    return 1;
}



Eigen::VectorXd lowMagnitudeFinitePeriodicDifference(const Eigen::Ref<const Eigen::ArrayXd>& t,
    const Eigen::Ref<const Eigen::ArrayXd>& f, const double period)
{
    const ptrdiff_t N = t.size();

    // the central difference is one of the estimates
    Eigen::ArrayXd centralDifference(N);
    centralDifference.segment(1, N - 2) = (f.tail(N - 2) - f.head(N - 2)) / (t.tail(N - 2) - t.head(N - 2));
    centralDifference[0] = (f[1] - f[N - 1]) / (t[1] - t[N - 1] + period);
    centralDifference[N - 1] = (f[0] - f[N - 2]) / (t[0] - t[N - 2] + period);

    // we can also look at the forward and backward differences
    Eigen::ArrayXd forwardDifference(N);
    forwardDifference.head(N - 1) = (f.tail(N - 1) - f.head(N - 1)) / (t.tail(N - 1) - t.head(N - 1));
    forwardDifference[N - 1] = (f[0] - f[N - 1]) / (t[0] - t[N - 1]);
    Eigen::ArrayXd reverseDifference(N);
    reverseDifference.tail(N - 1) = (f.tail(N - 1) - f.head(N - 1)) / (t.tail(N - 1) - t.head(N - 1));
    reverseDifference[0] = (f[0] - f[N - 1]) / (t[0] - t[N - 1]);

    // construct the low-magnitude result
    const Eigen::ArrayXd lowMagnitudeOneSidedDifference =
        (forwardDifference.cwiseAbs() < reverseDifference.cwiseAbs()).select(forwardDifference, reverseDifference);
    const Eigen::ArrayXd result = (centralDifference.cwiseAbs() < lowMagnitudeOneSidedDifference.cwiseAbs())
        .select(centralDifference, lowMagnitudeOneSidedDifference);

    // all done
    return result;
}

// this function should only be called when the circle only intersects the curve in one place
double computeCurveCircleIntersectionT(const Hexagon::Blade::Curve<2>& curve,
    const Eigen::Ref<const Eigen::Vector2d>& circleCenter, const double circleRadius)
{
    const Eigen::Vector2d contiguousCenter = circleCenter;

    double intersectionT;
    Eigen::Vector2d intersectionPoint;
    if (circleRadius < 1e-5 || !Hexagon::Blade::circleIntersection(curve, contiguousCenter.data(), circleRadius,
        intersectionPoint.data(), 0.0, 0.0, &intersectionT))
    {
        Hexagon::Blade::closestPoint(curve, contiguousCenter.data(), intersectionPoint.data(), &intersectionT);
    }
    return intersectionT;
}


double computeCurveLineIntersectionT(const Hexagon::Blade::Curve<2>& curve,
    const Eigen::Ref<const Eigen::Vector2d>& linePoint,
    const Eigen::Ref<const Eigen::Vector2d>& lineDirection)
{
    const Eigen::Vector2d contiguousPoint = linePoint;
    const Eigen::Vector2d contiguousDirection = lineDirection.normalized();

    double intersectionT;
    Eigen::Vector2d intersectionPoint;
    if (!Hexagon::Blade::lineIntersection(curve, contiguousPoint.data(), contiguousDirection.data(),
        intersectionPoint.data(), 0.0, 0.0, &intersectionT))
    {
        Hexagon::Blade::closestPoint(curve, contiguousPoint.data(), intersectionPoint.data(), &intersectionT);
    }
    return intersectionT;
}



// this function should only be called when the line only intersects the curve in one place


void CSection::ResetCurves()
{
    int i;

    for (i = 0; i < 5; i++)
    {
        if (m_meaPart[i])
            delete m_meaPart[i];
        m_meaPart[i] = 0;
    }

    if (m_meaCurve)
        delete m_meaCurve;
    m_meaCurve = 0;

    if (m_BCCurve)
        delete m_BCCurve;
    m_BCCurve = 0;

    if (m_cxpt)
        delete[] m_cxpt;
    if (m_cypt)
        delete[] m_cypt;
    if (m_czpt)
        delete[] m_czpt;

    if (m_ival)
        delete[] m_ival;
    if (m_jval)
        delete[] m_jval;
    if (m_kval)
        delete[] m_kval;

    if (m_mxpt)
        delete[] m_mxpt;
    if (m_mypt)
        delete[] m_mypt;

    if (m_nxpt)
        delete[] m_nxpt;
    if (m_nypt)
        delete[] m_nypt;

    if (m_nomt)
        delete[] m_nomt;

    if (m_partOf)
        delete[] m_partOf;

    m_totalPoints = 0;

    m_mxpt = m_mypt = m_nxpt = m_nypt = m_nomt = 0;
    m_ival = m_jval = m_kval = m_cxpt = m_cypt = m_czpt = 0;
    m_partOf = 0;

    m_meaPitch[2] = m_nomPitch[2] = -1.0;
}

CBestFit* CSection::GetBestFit(int index)
{
    CBestFit* bestfit = NULL;
    if (index >= 0 && index < m_numBestFits)
        bestfit = m_bestFits[index];
    bestfit = m_bestFits[0];
    return bestfit;
}

CBestFit* CSection::GetBestFitV1(int index)
{
    CBestFit* bestfit = NULL;
    // if(index >= 0 && index < m_numBestFits)
    bestfit = m_bestFits[0];
    return bestfit;
}



bool isNear_periodic(double t1, double t2, double period)
{
    return std::abs(std::remainder(t2 - t1, period)) < 1e-3;
}
bool isWithinOrNear_periodic(double t, const Hexagon::Blade::Curve<2>& curve)
{
    return isNear_periodic(t, curve.t0(), curve.period()) || isNear_periodic(t, curve.t1(), curve.period()) ||
        Hexagon::Blade::tIsInSubcurve(t, curve, curve.period());
}

std::unique_ptr<const Hexagon::Blade::LinearDeviation> makeLinearDeviationFromTValue(
    const Hexagon::Blade::Curve<2>& nominalCurve, const Hexagon::Blade::Curve<2>& measuredCurve,
    std::function<bool(double)> canMakeLinearDeviationHere, const double nominalTValue, const double measuredTValue)
{
    Eigen::Vector2d nominalPoint, nominalTangent;
    std::tie(nominalPoint, nominalTangent) =
        Hexagon::Blade::evaluateWithDerivative(nominalCurve, Eigen::Map<const Eigen::VectorXd>(&nominalTValue, 1));
    const Eigen::Vector2d measuredPoint =
        Hexagon::Blade::evaluate(measuredCurve, Eigen::Map<const Eigen::VectorXd>(&measuredTValue, 1));
    const Eigen::Vector2d nominalNormal = (Hexagon::Blade::makeRotate90() * nominalTangent).normalized();
    if (canMakeLinearDeviationHere(nominalTValue))
    {
        return std::make_unique<const Hexagon::Blade::LinearDeviation>(nominalPoint, nominalNormal, measuredPoint);
    }
    return nullptr;
}
double findMiddleT(const Hexagon::Blade::Curve<2>& curve, const Eigen::Vector2d& bounds)
{
    const Eigen::Matrix2d endPoints = Hexagon::Blade::evaluate(curve, bounds);
    const Eigen::Vector2d midPoint = endPoints.rowwise().mean();
    const Eigen::Vector2d crossLineDirection =
        (Hexagon::Blade::makeRotate90() * (endPoints.col(1) - endPoints.col(0))).normalized();
    return computeCurveLineIntersectionT(curve, midPoint, crossLineDirection);
}

std::unique_ptr<const Hexagon::Blade::LinearDeviation>
makeLinearDeviationFromMidpoint(const Hexagon::Blade::Curve<2>& nominalCurve,
    const Hexagon::Blade::Curve<2>& measuredCurve,
    std::function<bool(double)> canMakeLinearDeviationHere)
{
    const auto nominalBounds = Hexagon::Blade::parametricBounds(nominalCurve);
    const auto measuredBounds = Hexagon::Blade::parametricBounds(measuredCurve);
    const bool ok0 = canMakeLinearDeviationHere(nominalBounds[0]);
    const bool ok1 = canMakeLinearDeviationHere(nominalBounds[1]);
    if (ok0 == ok1) // use midpoint if both true or if both false
    {
        const double nominalMiddleT = findMiddleT(nominalCurve, nominalBounds);
        const double measuredMiddleT = findMiddleT(measuredCurve, measuredBounds);
        return makeLinearDeviationFromTValue(nominalCurve, measuredCurve, canMakeLinearDeviationHere, nominalMiddleT,
            measuredMiddleT);
    }
    else if (ok0)
    {
        return makeLinearDeviationFromTValue(nominalCurve, measuredCurve, canMakeLinearDeviationHere, nominalBounds[0],
            measuredBounds[0]);
    }
    else if (ok1)
    {
        return makeLinearDeviationFromTValue(nominalCurve, measuredCurve, canMakeLinearDeviationHere, nominalBounds[1],
            measuredBounds[1]);
    }
    return nullptr;
}

std::unique_ptr<const Hexagon::Blade::LinearDeviation>
meanLinearDeviation(std::unique_ptr<const Hexagon::Blade::LinearDeviation> a,
    std::unique_ptr<const Hexagon::Blade::LinearDeviation> b)
{
    if (!a && !b)
    {
        return nullptr;
    }
    if (a && !b)
    {
        return a;
    }
    if (!a && b)
    {
        return b;
    }
    if (a && b)
    {
        Eigen::Matrix2d vectors;
        vectors << a->nominalTangentDirection, b->nominalTangentDirection;
        const Eigen::Vector2d newTangent = vectors.jacobiSvd(Eigen::ComputeFullU).matrixU().col(0);
        return std::make_unique<const Hexagon::Blade::LinearDeviation>(
            0.5 * (a->nominalPoint + b->nominalPoint), newTangent, 0.5 * (a->measuredPoint + b->measuredPoint));
    }
    throw std::logic_error("This should be impossible.");
}

std::tuple<Eigen::VectorXd, std::vector<Hexagon::Blade::LinearDeviation>,
    std::unique_ptr<const Hexagon::Blade::Curve<2>>>
    figureOutWeightingAndEndpointConstraints(const CFitParams& fp, const CSection* section)
{
    // we need linear deviations sometimes to keep things confined
    const auto nominalSectionCurve = Hexagon::Blade::nominalSectionCurve(section);
    const int leType = section->LEType();
    const int teType = section->TEType();
    auto canMakeLinearDeviationHere = [&nominalSectionCurve, leType, teType](double nominalT) -> bool {
        if (isWithinOrNear_periodic(nominalT, *nominalSectionCurve.leading))
        {
            return leType != EDGE_PARTIAL;
        }
        if (isWithinOrNear_periodic(nominalT, *nominalSectionCurve.trailing))
        {
            return teType != EDGE_PARTIAL;
        }
        return true;
    };

    typedef std::unique_ptr<const CSubCurve> PointerType;
    typedef std::vector<Hexagon::Blade::LinearDeviation> Deviations;
    typedef std::tuple<Eigen::VectorXd, Deviations, std::unique_ptr<const Hexagon::Blade::Curve<2>>> ReturnType;

    // make some half curves, because we use them sometimes
    PointerType nominalConcaveHalfCurve = MakeNominalHalfCurve(section, CCC);
    PointerType nominalConvexHalfCurve = MakeNominalHalfCurve(section, CVC);
    PointerType measuredConcaveHalfCurve = MakeMeasuredHalfCurve(section, CCC);
    PointerType measuredConvexHalfCurve = MakeMeasuredHalfCurve(section, CVC);

    // first, figure out the points weighting, and some nominal points and vectors
    Eigen::VectorXd weightFittedPoints = Eigen::VectorXd::Zero(section->m_totalPoints);
    bool offsetsAreUnusable = (fp.leoff1 == 0.0 && fp.leoff2 == 0.0 && fp.teoff1 == 0.0 && fp.teoff2 == 0.0) ||
        fp.leoff1 < 0.0 || fp.leoff2 < 0.0 || fp.teoff1 < 0.0 || fp.teoff2 < 0.0 ||
        fp.leoff1 > fp.leoff2 || fp.teoff1 > fp.teoff2 ||
        (fp.leoff1 == fp.leoff2 && fp.teoff1 == fp.teoff2);
    if (offsetsAreUnusable && !fp.complexEdgeZone && !fp.chordZone)
    {
        // if there are no leading/trailing edge offsets, or any of the offsets are nonsense,
        // then the weighting is purely based on which part of the curve each point is assigned to
        const Eigen::Map<const Eigen::ArrayXi> partOf(section->m_partOf, section->m_totalPoints);
        weightFittedPoints = (partOf == CCC).select(fp.fitcurve[CCC] ? fp.weightcurve[CCC] : 0.0, weightFittedPoints);
        weightFittedPoints = (partOf == CVC).select(fp.fitcurve[CVC] ? fp.weightcurve[CVC] : 0.0, weightFittedPoints);
        weightFittedPoints =
            (partOf == LEC)
            .select(fp.fitcurve[LEC] && leType != EDGE_PARTIAL ? fp.weightcurve[LEC] : 0.0, weightFittedPoints);
        weightFittedPoints =
            (partOf == TEC)
            .select(fp.fitcurve[TEC] && teType != EDGE_PARTIAL ? fp.weightcurve[TEC] : 0.0, weightFittedPoints);

        // Do we need any linear deviations? That depends on which curves are involved in the fit
        // figure out the linear deviations
        if (fp.fitcurve[LEC] || fp.fitcurve[TEC])
        {
            // if the LEC and/or TEC are involved in the fit, we can get away without any linear deviations
            return ReturnType(weightFittedPoints, Deviations{}, nullptr);
        }
        else if (fp.fitcurve[CCC] && fp.fitcurve[CVC])
        {
            alwaysAssert(!fp.fitcurve[LEC] && !fp.fitcurve[TEC]);
            // With the CCC and CVC but not the LEC and not the TEC being fitted,
            // we need linear deviations to avoid the risk that the fit will
            // "slide" from side to side
            auto deviationCCC = makeLinearDeviationFromMidpoint(*nominalConcaveHalfCurve, *measuredConcaveHalfCurve,
                canMakeLinearDeviationHere);
            auto deviationCVC = makeLinearDeviationFromMidpoint(*nominalConvexHalfCurve, *measuredConvexHalfCurve,
                canMakeLinearDeviationHere);
            auto linearDeviation = meanLinearDeviation(std::move(deviationCCC), std::move(deviationCVC));
            alwaysAssert(linearDeviation);
            return ReturnType(weightFittedPoints, Deviations{ *linearDeviation }, nullptr);
        }
        else if (fp.fitcurve[CCC])
        {
            alwaysAssert(!fp.fitcurve[LEC] && !fp.fitcurve[TEC] && !fp.fitcurve[CVC] && fp.fitcurve[CCC]);
            auto linearDeviation = makeLinearDeviationFromMidpoint(*nominalConcaveHalfCurve, *measuredConcaveHalfCurve,
                canMakeLinearDeviationHere);
            alwaysAssert(linearDeviation);
            return ReturnType(weightFittedPoints, Deviations{ *linearDeviation }, std::move(nominalConcaveHalfCurve));
        }
        else if (fp.fitcurve[CVC])
        {
            alwaysAssert(!fp.fitcurve[LEC] && !fp.fitcurve[TEC] && fp.fitcurve[CVC] && !fp.fitcurve[CCC]);
            auto linearDeviation = makeLinearDeviationFromMidpoint(*nominalConvexHalfCurve, *measuredConvexHalfCurve,
                canMakeLinearDeviationHere);
            alwaysAssert(linearDeviation);
            return ReturnType(weightFittedPoints, Deviations{ *linearDeviation }, std::move(nominalConvexHalfCurve));
        }
    }

    // what are the measured points and their t-values?
    const Eigen::Matrix2Xd measuredPoints =
        constructPointMatrix(section->m_mxpt, section->m_mypt, section->m_totalPoints);
    Eigen::VectorXd measuredT(section->m_totalPoints);
    Hexagon::Blade::findClosestTValues(*section->MeaCurve(), measuredT.data(), measuredPoints.data(),
        section->m_totalPoints);

    // at this point, we either have a LEARC, a TEARC, or a complex edge fit
    int edgeCurve = -1;
    Eigen::Vector2d edgeOffset;
    bool edgeIsPartial = false;
    {
        alwaysAssert(fp.teoff1 < fp.teoff2);
        alwaysAssert(fp.teoff1 >= 0.0);
        edgeOffset << fp.teoff1, fp.teoff2;
        edgeCurve = TEC;
        edgeIsPartial = section->TEType() == EDGE_PARTIAL;
    }

    // figure out the tip locations
    const double nominalTipT = section->NomPart(edgeCurve)->Extreme();
    Eigen::Vector2d nominalTip;
    section->NomPart(edgeCurve)->CalcPoint(nominalTip.data(), nominalTipT);
    const double measuredTipT = section->MeaPart(edgeCurve)->Extreme();
    Eigen::Vector2d measuredTip;
    section->MeaPart(edgeCurve)->CalcPoint(measuredTip.data(), measuredTipT);

    // now, do we have an edge-arc fit or a complex fit?
    if (!offsetsAreUnusable && !fp.complexEdgeZone)
    {
        const Eigen::ArrayXd distanceToTip = (measuredPoints.colwise() - measuredTip).colwise().norm().transpose();
        const Eigen::ArrayXb isBetweenArcs = distanceToTip >= edgeOffset[0] && distanceToTip <= edgeOffset[1];

        // if we've gotten to this point, it's because we have a leading-edge-arc or trailing-edge-arc best-fit
        // (and not a complex fit)
        alwaysAssert(fp.leoff1 > 0.0 || fp.leoff2 > 0.0 || fp.teoff1 > 0.0 || fp.teoff2 > 0.0);
        alwaysAssert(fp.leoff1 >= 0.0 && fp.leoff2 >= 0.0 && fp.teoff1 >= 0.0 && fp.teoff2 >= 0.0);
        alwaysAssert(fp.leoff1 < fp.leoff2 || fp.teoff1 < fp.teoff2);
        alwaysAssert(fp.leoff1 <= fp.leoff2 && fp.teoff1 <= fp.teoff2);

        // are we fitting both sides or just one?
        alwaysAssert(fp.fitcurve[CCC] || fp.fitcurve[CVC]);
        if (fp.fitcurve[CCC] && fp.fitcurve[CVC])
        {
            // we're going to fit both sides, so we only need to weight based on distance to the leading edge
            weightFittedPoints = isBetweenArcs.select(1.0, weightFittedPoints);

            const double nominalConcaveT =
                computeCurveCircleIntersectionT(*nominalConcaveHalfCurve, nominalTip, edgeOffset[0]);
            const double nominalConvexT = computeCurveCircleIntersectionT(*nominalConvexHalfCurve, nominalTip, edgeOffset[0]);
            const double measuredConcaveT =
                computeCurveCircleIntersectionT(*measuredConcaveHalfCurve, measuredTip, edgeOffset[0]);
            const double measuredConvexT =
                computeCurveCircleIntersectionT(*measuredConvexHalfCurve, measuredTip, edgeOffset[0]);
            auto deviationCCC = makeLinearDeviationFromTValue(*section->NomPart(CCC), *section->MeaPart(CCC),
                canMakeLinearDeviationHere, nominalConcaveT, measuredConcaveT);
            auto deviationCVC = makeLinearDeviationFromTValue(*section->NomPart(CVC), *section->MeaPart(CVC),
                canMakeLinearDeviationHere, nominalConvexT, measuredConvexT);
            Deviations devs;
            if (!edgeIsPartial && deviationCCC && deviationCVC &&
                (deviationCCC->nominalPoint - deviationCVC->nominalPoint).norm() > edgeOffset[1] - edgeOffset[0])
            {
                devs.push_back(*makeLinearDeviationFromTValue(*section->NomPart(edgeCurve), *section->MeaPart(edgeCurve),
                    canMakeLinearDeviationHere, nominalTipT, measuredTipT));
                devs.push_back(devs.at(0));
                devs.at(1).nominalTangentDirection =
                    Eigen::Vector2d(devs.at(1).nominalTangentDirection[1], -devs.at(1).nominalTangentDirection[0]);
            }
            else
            {
                auto linearDeviation = meanLinearDeviation(std::move(deviationCCC), std::move(deviationCVC));
                if (linearDeviation)
                {
                    devs.push_back(*linearDeviation);
                }
            }
            return ReturnType(weightFittedPoints, devs, nullptr);
        }
        else if (fp.fitcurve[CCC])
        {
            // only use the ones on the concave side that are within the distance range
            alwaysAssert(!fp.fitcurve[CVC]);
            const Eigen::ArrayXb isConcave =
                Hexagon::Blade::tIsInSubcurve_eigen(measuredT, *measuredConcaveHalfCurve, section->MeaCurve()->Period());
            weightFittedPoints = (isBetweenArcs && isConcave).select(1.0, weightFittedPoints);

            const double nominalConcaveT =
                computeCurveCircleIntersectionT(*nominalConcaveHalfCurve, nominalTip, edgeOffset[0]);
            const double measuredConcaveT =
                computeCurveCircleIntersectionT(*measuredConcaveHalfCurve, measuredTip, edgeOffset[0]);
            auto deviationCCC = makeLinearDeviationFromTValue(*section->NomPart(CCC), *section->MeaPart(CCC),
                canMakeLinearDeviationHere, nominalConcaveT, measuredConcaveT);
            Deviations devs;
            if (deviationCCC)
            {
                devs.push_back(*deviationCCC);
            }
            return ReturnType(weightFittedPoints, devs, std::move(nominalConcaveHalfCurve));
        }
        else if (fp.fitcurve[CVC])
        {
            // only use the ones on the convex side that are within the distance range
            alwaysAssert(!fp.fitcurve[CCC]);
            const Eigen::ArrayXb isConvex =
                Hexagon::Blade::tIsInSubcurve_eigen(measuredT, *measuredConvexHalfCurve, section->MeaCurve()->Period());
            weightFittedPoints = (isBetweenArcs && isConvex).select(1.0, weightFittedPoints);

            const double nominalConvexT = computeCurveCircleIntersectionT(*nominalConvexHalfCurve, nominalTip, edgeOffset[0]);
            const double measuredConvexT =
                computeCurveCircleIntersectionT(*measuredConvexHalfCurve, measuredTip, edgeOffset[0]);
            auto deviationCVC = makeLinearDeviationFromTValue(*section->NomPart(CVC), *section->MeaPart(CVC),
                canMakeLinearDeviationHere, nominalConvexT, measuredConvexT);
            Deviations devs;
            if (deviationCVC)
            {
                devs.push_back(*deviationCVC);
            }
            return ReturnType(weightFittedPoints, devs, std::move(nominalConvexHalfCurve));
        }
        else
        {
            throw std::logic_error("This should be impossible; " + (__FILE__ + std::to_string(__LINE__)));
        }
    }
    else
    {
        alwaysAssert(fp.complexEdgeZone);
        const Hexagon::Blade::Curve<2>* const nominalMCL = section->NomPart(MCC);
        const Hexagon::Blade::Curve<2>* const measuredMCL = section->MeaPart(MCC);

        Eigen::Vector2d nominalMCLTRange, measuredMCLTRange;
        nominalMCLTRange[0] = computeCurveCircleIntersectionT(*nominalMCL, nominalTip, edgeOffset[0]);
        nominalMCLTRange[1] = computeCurveCircleIntersectionT(*nominalMCL, nominalTip, edgeOffset[1]);
        measuredMCLTRange[0] = computeCurveCircleIntersectionT(*measuredMCL, measuredTip, edgeOffset[0]);
        measuredMCLTRange[1] = computeCurveCircleIntersectionT(*measuredMCL, measuredTip, edgeOffset[1]);

        Eigen::Matrix2d nominalMCLPoints, nominalMCLTangents;
        std::tie(nominalMCLPoints, nominalMCLTangents) =
            Hexagon::Blade::evaluateWithDerivative(*nominalMCL, nominalMCLTRange);
        nominalMCLTangents.colwise().normalize();
        Eigen::Matrix2d nominalMCLOrthogonals;
        nominalMCLOrthogonals.row(0) = nominalMCLTangents.row(1);
        nominalMCLOrthogonals.row(1) = -nominalMCLTangents.row(0);
        nominalMCLOrthogonals.colwise().normalize();
        const Eigen::Matrix2d measuredMCLPoints = Hexagon::Blade::evaluate(*measuredMCL, measuredMCLTRange);

        const Eigen::Isometry2d alignToNominal = Hexagon::Blade::twoPointBestFit(
            nominalMCLPoints.col(0), nominalMCLPoints.col(1), measuredMCLPoints.col(0), measuredMCLPoints.col(1));

        // now we have the zones in the measured space
        const Eigen::Matrix2d measuredZonePoints = alignToNominal.inverse() * nominalMCLPoints;
        const Eigen::Matrix2d measuredZoneNormals =
            (alignToNominal.inverse().linear() * nominalMCLTangents).colwise().normalized();
        const Eigen::Matrix2d measuredZoneDirections =
            (alignToNominal.inverse().linear() * nominalMCLOrthogonals).colwise().normalized();
        const double normalSign = Hexagon::Blade::sign(nominalMCLTRange[1] - nominalMCLTRange[0]);
        const Eigen::ArrayXb isWithinZone0 =
            ((measuredPoints.colwise() - measuredZonePoints.col(0)).transpose() * measuredZoneNormals.col(0) * normalSign)
            .array() >= 0.0;
        const Eigen::ArrayXb isWithinZone1 =
            ((measuredPoints.colwise() - measuredZonePoints.col(1)).transpose() * measuredZoneNormals.col(1) * normalSign)
            .array() <= 0.0;
        const Eigen::ArrayXb isWithinZone =
            (edgeOffset[0] > 1e-5) ? Eigen::ArrayXb(isWithinZone0 && isWithinZone1) : isWithinZone1;

        // intersect the measured curves and the nominal curves with the zone boundaries
        // are we fitting both sides or just one?
        alwaysAssert(fp.fitcurve[CCC] || fp.fitcurve[CVC]);
        if (fp.fitcurve[CCC] && fp.fitcurve[CVC])
        {
            // we're going to fit both sides, so we only need to weight based on distance to the leading edge
            weightFittedPoints = isWithinZone.select(1.0, weightFittedPoints);

            std::unique_ptr<const Hexagon::Blade::LinearDeviation> deviationCCC, deviationCVC;
            if (edgeOffset[0] > 1e-5)
            {
                const double nominalConcaveT = computeCurveLineIntersectionT(*nominalConcaveHalfCurve, nominalMCLPoints.col(0),
                    nominalMCLOrthogonals.col(0));
                const double nominalConvexT = computeCurveLineIntersectionT(*nominalConvexHalfCurve, nominalMCLPoints.col(0),
                    nominalMCLOrthogonals.col(0));
                const double measuredConcaveT = computeCurveLineIntersectionT(
                    *measuredConcaveHalfCurve, measuredZonePoints.col(0), measuredZoneDirections.col(0));
                const double measuredConvexT = computeCurveLineIntersectionT(
                    *measuredConvexHalfCurve, measuredZonePoints.col(0), measuredZoneDirections.col(0));
                deviationCCC = makeLinearDeviationFromTValue(*section->NomPart(CCC), *section->MeaPart(CCC),
                    canMakeLinearDeviationHere, nominalConcaveT, measuredConcaveT);
                deviationCVC = makeLinearDeviationFromTValue(*section->NomPart(CVC), *section->MeaPart(CVC),
                    canMakeLinearDeviationHere, nominalConvexT, measuredConvexT);
            }
            else
            {
                deviationCCC = makeLinearDeviationFromTValue(*section->NomPart(edgeCurve), *section->MeaPart(edgeCurve),
                    canMakeLinearDeviationHere, nominalTipT, measuredTipT);
                deviationCVC = makeLinearDeviationFromTValue(*section->NomPart(edgeCurve), *section->MeaPart(edgeCurve),
                    canMakeLinearDeviationHere, nominalTipT, measuredTipT);
            }
            Deviations devs;
            if (!edgeIsPartial && deviationCCC && deviationCVC &&
                (deviationCCC->nominalPoint - deviationCVC->nominalPoint).norm() > edgeOffset[1] - edgeOffset[0])
            {
                devs.push_back(*makeLinearDeviationFromTValue(*section->NomPart(edgeCurve), *section->MeaPart(edgeCurve),
                    canMakeLinearDeviationHere, nominalTipT, measuredTipT));
                devs.push_back(devs.at(0));
                devs.at(1).nominalTangentDirection =
                    Eigen::Vector2d(devs.at(1).nominalTangentDirection[1], -devs.at(1).nominalTangentDirection[0]);
            }
            else
            {
                auto linearDeviation = meanLinearDeviation(std::move(deviationCCC), std::move(deviationCVC));
                if (linearDeviation)
                {
                    devs.push_back(*linearDeviation);
                }
            }
            return ReturnType(weightFittedPoints, devs, nullptr);
        }
        else if (fp.fitcurve[CCC])
        {
            // only use the ones on the concave side that are within the distance range
            alwaysAssert(!fp.fitcurve[CVC]);
            const Eigen::ArrayXb isConcave =
                Hexagon::Blade::tIsInSubcurve_eigen(measuredT, *measuredConcaveHalfCurve, section->MeaCurve()->Period());
            weightFittedPoints = (isWithinZone && isConcave).select(1.0, weightFittedPoints);

            std::unique_ptr<const Hexagon::Blade::LinearDeviation> deviationCCC;
            if (edgeOffset[0] > 1e-5)
            {
                const double nominalConcaveT = computeCurveLineIntersectionT(*nominalConcaveHalfCurve, nominalMCLPoints.col(0),
                    nominalMCLOrthogonals.col(0));
                const double measuredConcaveT = computeCurveLineIntersectionT(
                    *measuredConcaveHalfCurve, measuredZonePoints.col(0), measuredZoneDirections.col(0));
                deviationCCC = makeLinearDeviationFromTValue(*section->NomPart(CCC), *section->MeaPart(CCC),
                    canMakeLinearDeviationHere, nominalConcaveT, measuredConcaveT);
            }
            else
            {
                deviationCCC = makeLinearDeviationFromTValue(*section->NomPart(edgeCurve), *section->MeaPart(edgeCurve),
                    canMakeLinearDeviationHere, nominalTipT, measuredTipT);
            }
            Deviations devs;
            if (deviationCCC)
            {
                devs.push_back(*deviationCCC);
            }
            return ReturnType(weightFittedPoints, devs, std::move(nominalConcaveHalfCurve));
        }
        else if (fp.fitcurve[CVC])
        {
            // only use the ones on the convex side that are within the distance range
            alwaysAssert(!fp.fitcurve[CCC]);
            const Eigen::ArrayXb isConvex =
                Hexagon::Blade::tIsInSubcurve_eigen(measuredT, *measuredConvexHalfCurve, section->MeaCurve()->Period());
            weightFittedPoints = (isWithinZone && isConvex).select(1.0, weightFittedPoints);

            std::unique_ptr<const Hexagon::Blade::LinearDeviation> deviationCVC;
            if (edgeOffset[0] > 1e-5)
            {
                const double nominalConvexT = computeCurveLineIntersectionT(*nominalConvexHalfCurve, nominalMCLPoints.col(0),
                    nominalMCLOrthogonals.col(0));
                const double measuredConvexT = computeCurveLineIntersectionT(
                    *measuredConvexHalfCurve, measuredZonePoints.col(0), measuredZoneDirections.col(0));
                deviationCVC = makeLinearDeviationFromTValue(*section->NomPart(CVC), *section->MeaPart(CVC),
                    canMakeLinearDeviationHere, nominalConvexT, measuredConvexT);
            }
            else
            {
                deviationCVC = makeLinearDeviationFromTValue(*section->NomPart(edgeCurve), *section->MeaPart(edgeCurve),
                    canMakeLinearDeviationHere, nominalTipT, measuredTipT);
            }
            Deviations devs;
            if (deviationCVC)
            {
                devs.push_back(*deviationCVC);
            }
            return ReturnType(weightFittedPoints, devs, std::move(nominalConvexHalfCurve));
        }
        else
        {
            throw std::logic_error("This should be impossible; " + (__FILE__ + std::to_string(__LINE__)));
        }
    }
}