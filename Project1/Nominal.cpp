#include "stdafx.h"
#include "Nominal.h"
#include "MeanCamberCurve.h"
#ifdef _DEBUG
#ifndef DBG_NEW
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
#define new DBG_NEW
#endif
#endif  // _DEBUG

CNominal::CNominal(void)
{
    m_sections = NULL;
}

CNominal::~CNominal(void)
{
    if (m_sections)
    {
        int i;
        for (i = 0; i < m_numSect; i++)
            delete m_sections[i];
        delete[] m_sections;
    }
}

bool CNominal::WriteFile(wchar_t* filename)
{
    FILE* fp;
    _wfopen_s(&fp, filename, L"wb");
    if (fp == NULL)
        return false;

#pragma warning(suppress : 4309)
    short ft = static_cast<short>(GOODTYPE);
    fwrite(&ft, sizeof(short), 1, fp);   // file format
    fwrite(&m_numSect, sizeof(short), 1, fp);   // number of sections
    fwrite(&m_english, sizeof(bool), 1, fp);

    wchar_t tbuf[MAXBUFSZ];

    int s;
    for (s = 0; s < m_numSect; s++)
    {
        fwrite(&m_sections[s]->m_zval, sizeof(double), 1, fp);
        wcscpy_s(tbuf, MAXBUFSZ, m_sections[s]->m_name);
        fwrite(tbuf, 1, 24, fp);
        fwrite(&m_sections[s]->m_letype, sizeof(short), 1, fp);
        fwrite(&m_sections[s]->m_tetype, sizeof(short), 1, fp);

        short skew = m_sections[s]->m_skewed;
        if (skew != 0)
            skew += (short)(10 * m_sections[s]->m_skewReport);

        fwrite(&skew, sizeof(short), 1, fp);
        if (m_sections[s]->m_skewed)
        {
            fwrite(m_sections[s]->m_skew, sizeof(double), 3, fp);
            fwrite(m_sections[s]->m_skewOrigin, sizeof(double), 3, fp);
        }
        else
            fwrite(tbuf, 1, 48, fp);

        long storeSize = 0L;
        if (m_sections[s]->m_whole)
            storeSize = m_sections[s]->m_whole->StoreSize();
        fwrite(&storeSize, sizeof(long), 1, fp);
        bugout(0, L"WriteMATHFile: storeSize(%ld) m_period(%f)", storeSize, &m_sections[s]->m_period);

        if (storeSize > 0)
            m_sections[s]->m_whole->Write(fp);
        fwrite(&m_sections[s]->m_period, sizeof(double), 1, fp);
        double extreme = m_sections[s]->m_lec->Extreme();
        fwrite(&extreme, sizeof(double), 1, fp);
        if (m_sections[s]->m_tec != NULL)
        {
            extreme = m_sections[s]->m_tec->Extreme();
        }
        else
        {
            extreme = 0;
        }

        fwrite(&extreme, sizeof(double), 1, fp);
        fwrite(m_sections[s]->m_pitch, sizeof(double), 3, fp);

        int i;
        for (i = 0; i < 4; i++)
        {
            double t0 = m_sections[s]->m_t0[i];
            double t1 = m_sections[s]->m_t1[i];
            bugout(0, L"WriteMATHFile:type=%d, t0(%f) t1(%f)", i, t0, t1);

            fwrite(&t0, sizeof(double), 1, fp);
            fwrite(&t1, sizeof(double), 1, fp);
        }

        storeSize = 0L;
        if (m_sections[s]->m_mc)
            storeSize = m_sections[s]->m_mc->StoreSize();
        fwrite(&storeSize, sizeof(long), 1, fp);
        if (storeSize > 0)
            m_sections[s]->m_mc->Write(fp);

        //fwrite(&m_sections[s]->m_fixedAxis, sizeof(int), 1, fp);
    }

    // save raw points
    for (s = 0; s < m_numSect; s++)
    {
        fwrite(&m_sections[s]->m_npts, sizeof(int), 1, fp);
        Hexagon::Blade::writeMeanCamberCurveParameters2016(m_sections[s]->m_mclParams, fp);
        for (int i = 0; i < m_sections[s]->m_npts; i++)
        {
            fwrite(&m_sections[s]->m_ox[i], sizeof(double), 1, fp);
            fwrite(&m_sections[s]->m_oy[i], sizeof(double), 1, fp);
            double ival = m_sections[s]->m_havek ? m_sections[s]->m_oi[i] : 0.0;
            fwrite(&ival, sizeof(double), 1, fp);
            double jval = m_sections[s]->m_havek ? m_sections[s]->m_oj[i] : 0.0;
            fwrite(&jval, sizeof(double), 1, fp);
            double kval = m_sections[s]->m_havek ? m_sections[s]->m_ok[i] : 0.0;
            fwrite(&kval, sizeof(double), 1, fp);
            double mtol = m_sections[s]->m_haveTol ? m_sections[s]->m_omtol[i] : 0.0;
            fwrite(&mtol, sizeof(double), 1, fp);
            double ptol = m_sections[s]->m_haveTol ? m_sections[s]->m_optol[i] : 0.0;
            fwrite(&ptol, sizeof(double), 1, fp);
        }
    }
    
    /************OpenArea infos*************************************/
    for (s = 0; s < m_numSect; s++)
    {
        int tmp_AreaCount = 0;
        tmp_AreaCount = m_sections[s]->m_areaCount;
        fwrite(&tmp_AreaCount, sizeof(int), 1, fp); //保存开口区域个数
        if (tmp_AreaCount > 0)
        {
            for (int i = 0; i < tmp_AreaCount; i++)
            {
                fwrite(&m_sections[s]->m_areaStartPoint[i], sizeof(double), 2, fp);
                fwrite(&m_sections[s]->m_areaEndPoint[i], sizeof(double), 2, fp);
                fwrite(&m_sections[s]->m_areaIndex[i], sizeof(int), 2, fp);
            }
        }
    }

    fclose(fp);
    return true;
}