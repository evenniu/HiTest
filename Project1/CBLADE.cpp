#include "stdafx.h"
#include "CBLADE.h"
#include "Nurb.h"
#include "SubCurve.h"
#include "Analysis.h"
#include <io.h>
#include <share.h>
#include <fcntl.h>
#include "HermiteCurve.h"
//#include "MeanCamberCurve.h"
#include "BestFits.h"
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
    //m_section = NULL;
}
CBlade::CBlade(wchar_t* mathFileName)
{
    m_nomValid = false;
    m_measLoaded = false;
    m_numSections = 0;
    //m_pBestFit = NULL;
    //m_pBestFitLE = NULL;
    //m_section = NULL;
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
    //if (m_section)
    //{
    //    for (int i = 0; i < m_numSections; i++)
    //        delete m_section[i];
    //    delete[]m_section;
    //}

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
    return true;
}
//bool CBlade::FitBlade(const CFitParams* fp, int sec1, int sec2)
//{
//}
