#include "stdafx.h"
#include "Circle.h"
#include "SubCurve.h"
#include "SmallestCircle.h"
#include "SECTION.h"
#include "SectionCurve.h"
//#include "BestFits.h"
#include "ArraySlicing.h"
#include "HermiteCurve.h"
#include "EigenAbstractCurve.h"
#include "TemplateHermiteSpline.h"
#include "ToleranceSection.h"
#pragma warning(push)
#pragma warning(disable : 4267)
#pragma warning(disable : 4244)
#pragma warning(disable : 4100)
#include <nanoflann/nanoflann.hpp>
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
}double CSection::FindKValue(double* pt)
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
bool CSection::AssignPoints(double* xv, double* yv, int n, int* /*start*/, int* /*end*/)
{
    // this is for closed curve for P&W AS file read and analysis file
    bugout(0, L"CSection::AssignPoints: m_totalPoints %d", m_totalPoints);

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
    for (int s = 0; s < 4; s++)
    {
        t0[s] = m_meaPart[s]->T0();
        t1[s] = m_meaPart[s]->T1();
    }
    // polygonalize the measured curve
    auto measuredPolygon = Blade::polygonalizeWithT(*m_meaCurve, 2048, 1e-4);
    Eigen::MatrixX2d measuredPoints = std::get<0>(measuredPolygon).transpose();
    Eigen::VectorXd measuredT = std::get<1>(measuredPolygon);
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