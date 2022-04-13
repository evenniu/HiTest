#include "stdafx.h"
#include "CBLADE.h"
#include "Nurb.h"
#include "SubCurve.h"
#include "Analysis.h"
#include <io.h>
#include <share.h>
#include <fcntl.h>
#include "HermiteCurve.h"
#include "MeanCamberCurve.h"
//#include "BestFits.h"
#include "SectionCurve.h"
#include "ArraySlicing.h"
#include "EigenAbstractCurve.h"
#ifdef _DEBUG
#ifndef DBG_NEW
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
#define new DBG_NEW
#endif
#endif  // _DEBUG

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CBlade::CBlade()
{
    m_nomValid = false;
    m_measLoaded = false;
    m_numSections = 0;
    //m_pBestFit = NULL;
    //m_pBestFitLE = NULL;
    m_section = NULL;
}
CBlade::CBlade(wchar_t* mathFileName)
{
    m_nomValid = false;
    m_measLoaded = false;
    m_numSections = 0;
    //m_pBestFit = NULL;
    //m_pBestFitLE = NULL;
    m_section = NULL;
    wcscpy_s(m_mathFileName, mathFileName);

    FILE* fp;
    _wfopen_s(&fp, mathFileName, L"rb");

    if (fp != NULL)
    {
        if (ReadFile(fp))
            m_nomValid = true;

        fclose(fp);
    }
}
CBlade::~CBlade()
{
    if (m_section)
    {
        for (int i = 0; i < m_numSections; i++)
            delete m_section[i];
        delete[]m_section;
    }

    //if (m_pBestFit)
    //    delete m_pBestFit;

    //if (m_pBestFitLE)
    //    delete m_pBestFitLE;
}

bool CBlade::IsValid()
{
    return m_nomValid;
}
bool isHermiteOpenCurve(FILE* fp, long storeSize)
{
    using namespace Hexagon;
    using namespace Blade;

    long curvePosition = ftell(fp);
    if (storeSize >= static_cast<long>(sizeof(wchar_t) * HermiteOpenCurve::fileDescriptorLength))
    {
        wchar_t hermiteDescriptor[HermiteOpenCurve::fileDescriptorLength];
        fread(hermiteDescriptor, sizeof(wchar_t), HermiteOpenCurve::fileDescriptorLength, fp);
        if (wcsncmp(hermiteDescriptor, HermiteOpenCurve::fileDescriptor, HermiteOpenCurve::fileDescriptorLength) == 0)
        {
            fseek(fp, curvePosition, 0);
            return true;
        }
    }
    fseek(fp, curvePosition, 0);
    return false;
}
CCurve* readCurve(FILE* fp, const bool isEnglish)
{
    long storeSize;
    fread(&storeSize, sizeof(long), 1, fp);
    bugout(0, L"readCurve:storeSize(%ld) ****", storeSize);
    if (storeSize > 0)
    {
        if (isHermiteOpenCurve(fp, storeSize))
        {
            return new Hexagon::Blade::HermiteOpenCurve(fp, isEnglish);
        }
        else
        {
            return new CNurbCurve(fp, isEnglish);
        }
    }
    return nullptr;
}
bool CBlade::ReadFile(FILE* fp)
{
    short ft;
    fread(&ft, sizeof(short), 1, fp);
    fread(&m_numSections, sizeof(short), 1, fp);

    fread(&m_english, sizeof(bool), 1, fp);
    bugout(0, L"CBlade ReadFile:num sect %d, m_english=%d", m_numSections, m_english);

    m_section = new CSection * [m_numSections];
    wchar_t tbuf[MAXBUFSZ];
    int i;
    for (i = 0; i < m_numSections; i++)
    {
        CSection* sect = new CSection;
        m_section[i] = sect;
        double d;
        fread(&d, sizeof(double), 1, fp);
        m_section[i]->ZValue(d);

        fread(tbuf, 1, 24, fp);
        m_section[i]->Name(tbuf);        //wcscpy_s(tbuf, MAXBUFSZ, L"A-A");//由于math文件空，测试时时防止出现错误暂赋初值
        bugout(0, L"CBlade::ReadFile name %s", tbuf);

        short s=0;
        fread(&s, sizeof(short), 1, fp);
        m_section[i]->LEType(s);
        //bugout(0, L"letype %d", s);

        fread(&s, sizeof(short), 1, fp);
        m_section[i]->TEType(s);
        //bugout(0, L"tetype %d", s);

        short skewed;
        fread(&skewed, sizeof(short), 1, fp);
        bugout(0, L"CBlade::ReadFile skewed %d", skewed);
        if (skewed)
        {
            double zaxis[3], origin[3];

            if (!fread(zaxis, sizeof(double), 3, fp))
                break;
            if (!fread(origin, sizeof(double), 3, fp))
                break;

            CAlignment* align = new CAlignment(3);
            align->MatrixFromAxis(zaxis);
            align->m_morig[0] = origin[0];
            align->m_morig[1] = origin[1];
            align->m_morig[2] = origin[2];
            align->m_borig[0] = origin[0];
            align->m_borig[1] = origin[1];
            align->m_borig[2] = origin[2];
           // m_section[i]->SkewAlign(align);

            if (skewed > 9)
                m_section[i]->m_skewReport = 1;
        }
        else
        {
            if (!fread(tbuf, 1, 48, fp))   // extra stuff
                break;
        }
        m_section[i]->NomCurve(readCurve(fp, m_english));

        double period, leExtreme, teExtreme, pitch[3];
        fread(&period, sizeof(double), 1, fp);    // could probably get this from whole curve
        bugout(0, L"CBlade ReadFile period %lf ** L 202", period);
        fread(&leExtreme, sizeof(double), 1, fp);
        fread(&teExtreme, sizeof(double), 1, fp);
        fread(pitch, sizeof(double), 3, fp);

        m_section[i]->NomPitch(pitch);
    
        for (int j = 0; j < 4; j++)
        {
            double t0, t1;

            fread(&t0, sizeof(double), 1, fp);
            fread(&t1, sizeof(double), 1, fp);
            bugout(0, L"CBlade ReadFile: j=%d t0=%f, t1=%f", j, t0, t1);
            CSubCurve* sc = new CSubCurve(m_section[i]->NomCurve(), t0, t1, period);
            if (j == LEC)
                sc->Extreme(leExtreme);
            else if (j == TEC)
                sc->Extreme(teExtreme);
            m_section[i]->NomPart(j, sc);
        }
        // read in the mean camber line straight from the file
        bugout(0, L"readCurve:readCurve  ****");
        m_section[i]->NomPart(4, readCurve(fp, m_english));
    }
    // read raw nominal points
    for (i = 0; i < m_numSections; i++)
    {
        int numPts;
        fread(&numPts, sizeof(int), 1, fp);
        auto mclParams = Hexagon::Blade::readMeanCamberCurveParameters2016(fp);
        m_section[i]->MakeNomArrays(numPts, mclParams.get());
        //bugout(0, L"ReadFile: AddTol,  SectionID=%d ****,numpoints=%d", i, numPts);
        for (int j = 0; j < numPts; j++)
        {
            double xyijk[5];
            fread(xyijk, sizeof(double), 5, fp);
            m_section[i]->AddNomXYIJK(j, xyijk);
            double tolerances[2];
            fread(tolerances, sizeof(double), 2, fp);
            m_section[i]->AddTol(j, tolerances);
        }
    }
    /************OpenArea infos*************************************/
    for (int s = 0; s < m_numSections; s++)
    {
        int tmp_areaCount = 0;
        fread(&tmp_areaCount, sizeof(int), 1, fp); //读取开口区域个数
        bugout(0, L"CBLADE::ReadFile: SectionID=%d, areaCount=%d", s, tmp_areaCount);
        m_section[s]->m_openareaCount = tmp_areaCount;
        if (tmp_areaCount > 0)
        {
            for (int i = 0; i < tmp_areaCount; i++)
            {
                double tmp_areaStartpoint[2];
                double tmp_areaEndpoint[2];
                fread(tmp_areaStartpoint, sizeof(double), 2, fp);
                fread(tmp_areaEndpoint, sizeof(double), 2, fp);

                m_section[s]->m_areaStartPoint[i][0] = tmp_areaStartpoint[0];
                m_section[s]->m_areaStartPoint[i][1] = tmp_areaStartpoint[1];
                m_section[s]->m_areaEndPoint[i][0] = tmp_areaEndpoint[0];
                m_section[s]->m_areaEndPoint[i][1] = tmp_areaEndpoint[1];
                int tmp_areaIndex[2];
                fread(tmp_areaIndex, sizeof(int), 2, fp);
                m_section[s]->m_areaIndex[i][0] = tmp_areaIndex[0];
                m_section[s]->m_areaIndex[i][1] = tmp_areaIndex[1];
                bugout(0, L"CBLADE::ReadFile:(%f,%f)(%f %f) -[%d %d] ", tmp_areaStartpoint[0], tmp_areaStartpoint[1], tmp_areaEndpoint[0], tmp_areaEndpoint[1], 
                    tmp_areaIndex[1], tmp_areaIndex[1]);
            }
        }
    }
    return true;
}
// as of November 2017, these are defined in SECTION.CPP
Eigen::Matrix2Xd constructPointMatrix(const double* xValues, const double* yValues, const ptrdiff_t n);

//bool CBlade::FitBlade(const CFitParams* fp, int sec1, int sec2)
//{
//}
