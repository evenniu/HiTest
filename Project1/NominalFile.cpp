#include "stdafx.h"
#include "NominalFile.h"

#include <fstream>
#include "Plane.h"
#include "MeanCamberCurve.h"

#include <string>

CNominalFile::CNominalFile(wchar_t* nominalFile, wchar_t* rptPath, wchar_t* mathFile, HWND statusHWND, wchar_t* processNomSection, bool isCreating)
{
	m_nom = NULL;
	m_statusBarHWND = statusHWND;
	wcscpy_s(m_processNomSection, processNomSection);

	wcscpy_s(m_mathFile, mathFile);
	wcscpy_s(m_nomFile, nominalFile);

	//if (!m_error->isOK())
	//	return;
	if (m_mathFile[0] == 0)
	{
		int index = ReverseFind(nominalFile, wchar_t('.'));

		wcscpy_s(m_mathFile, nominalFile);
		if (index > 0)
			m_mathFile[index] = 0;

		wcscat_s(m_mathFile, L".MTH");
	}

	m_readonly = false;
	if (wcscmp(m_mathFile, L"READONLY") == 0)
		m_readonly = true;

	wcscpy_s(mathFile, 200, m_mathFile);

	m_nom = new CNominal();

	m_firstWarning = -1;
	bugout(0, L"new CNominalFile(): will cal Translate ");
	Translate(&m_firstWarning);

	bugout(0, L"CNominalFile:Translate done ! ");
	//fclose(m_fp);
}

CNominalFile::~CNominalFile()
{
}

bool CNominalFile::Translate(int* firstWarning)
{
	bugout(0, L"CNominalFile ::Translate : entered");
	int i;
	wchar_t key[MAXBUFSZ], rest[MAXBUFSZ], buf[MAXBUFSZ];
	key[0] = rest[0] = buf[0] = 0;
	wchar_t tmp_secName[MAXBUFSZ];
	tmp_secName[0] = 0;
	m_nom->m_numSect = 0;
	int ok = 0;
	std::ifstream  fin(m_nomFile , std::ios::in);
	if (!fin)
	{
		return false;
	}
	string line;
	string::size_type idx;
	const char* strline;
	int strlinesize = 0;
	while (getline(fin, line)) // line中不包括每行的换行符
	{
		idx = line.find("UNITS");
		if (idx != string::npos)
		{
			ok = 1;
			line = line.substr(6);
			strline = line.c_str();
			strlinesize = line.length();
			MultiByteToWideChar(CP_ACP, 0, line.c_str(), strlen(strline) + 1, rest, MAXBUFSZ);
			if (wcsncmp(rest, L"MM", 2) == 0)
			{
				m_nom->m_english = false;
			}
			else if (wcsncmp(rest, L"IN", 2) == 0)
			{
				m_nom->m_english = true;
			}
			else
			{
				ok = 0;
			}
			continue;
		}
		if (!ok)
		{
			/*ErrorStruct es(BE_BADNOMUNITS);
			m_error->AddError(&es);*/
			return false;
		}
		
		idx = line.find("NUM_SECT");
		if (idx != string::npos)
		{
			line = line.substr(9);
			strline = line.c_str();
			strlinesize = line.length();
			MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, rest, MAXBUFSZ);
			int tmp_m_numSect = 0;
			int trash_usedToBeVoronoi;
			int ns = swscanf_s(rest, L"%hd %d", &tmp_m_numSect, &trash_usedToBeVoronoi);
			if (ns > 1)
			{
				m_nom->m_numSect = tmp_m_numSect;
				ok = 1;
			}
		}
		if (!ok)
		{
			/*ErrorStruct es(BE_BADNOMUNITS);
			m_error->AddError(&es);*/
			return false;
		}

		m_nom->m_sections = new CNominalSection * [m_nom->m_numSect];
		*firstWarning = -1;
		int s;
		for (s = 0; s < m_nom->m_numSect; s++)
			m_nom->m_sections[s] = NULL;

		for (s = 0; s < m_nom->m_numSect; s++)
		{
			int ns, n;
			m_nom->m_sections[s] = new CNominalSection();
			getline(fin, line );
			
			wchar_t tmp[MAXBUFSZ], tmp2[MAXBUFSZ], tmp3[MAXBUFSZ];

			if (line.find("SECTION") != string::npos)
			{
				
				//line = line.substr(8);
				strline = line.c_str();
				strlinesize = line.length();
				MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, rest, MAXBUFSZ);
				ns = swscanf_s(rest, L"%s %s %d %lf %s %s", key, MAXBUFSZ, tmp, MAXBUFSZ, &n, &m_nom->m_sections[s]->m_zval, tmp2, MAXBUFSZ, tmp3, MAXBUFSZ);
				if(ns<4 || wcscmp(key, L"SECTION") != 0 || n < 10)
				{
					ErrorStruct es(BE_MISSINGORBADSECTION, s + 1);
					return false;
				}
				if (ns > 4)
				{
					if (wcscmp(tmp2, L"SKEWED") == 0)
					{
						m_nom->m_sections[s]->m_skewed = 1;
						m_nom->m_sections[s]->m_skewReport = 0;
					}
					else if (wcscmp(tmp2, L"SKEWED_BLADE") == 0)
					{
						m_nom->m_sections[s]->m_skewed = 1;
						m_nom->m_sections[s]->m_skewReport = 1;
					}
					else if (wcscmp(tmp2, L"+X") == 0)
						m_nom->m_sections[s]->m_convexSide = 1;
					else if (wcscmp(tmp2, L"-X") == 0)
						m_nom->m_sections[s]->m_convexSide = 2;
					else if (wcscmp(tmp2, L"+Y") == 0)
						m_nom->m_sections[s]->m_convexSide = 3;
					else if (wcscmp(tmp2, L"-Y") == 0)
						m_nom->m_sections[s]->m_convexSide = 4;
				}
				if (ns > 5)
				{
					if (wcscmp(tmp3, L"SKEWED") == 0)
					{
						m_nom->m_sections[s]->m_skewed = 1;
						m_nom->m_sections[s]->m_skewReport = 0;
					}
					else if (wcscmp(tmp3, L"SKEWED_BLADE") == 0)
					{
						m_nom->m_sections[s]->m_skewed = 1;
						m_nom->m_sections[s]->m_skewReport = 1;
					}
					else if (wcscmp(tmp3, L"+X") == 0)
						m_nom->m_sections[s]->m_convexSide = 1;
					else if (wcscmp(tmp3, L"-X") == 0)
						m_nom->m_sections[s]->m_convexSide = 2;
					else if (wcscmp(tmp3, L"+Y") == 0)
						m_nom->m_sections[s]->m_convexSide = 3;
					else if (wcscmp(tmp3, L"-Y") == 0)
						m_nom->m_sections[s]->m_convexSide = 4;
				}

				wcscpy_s(m_nom->m_sections[s]->m_name, tmp);
				MakeUpper(m_nom->m_sections[s]->m_name);
				m_nom->m_sections[s]->SetSize(n, 1, 1, m_nom->m_english);  // 1 means has vectors, second 1 means has tolerances

				bugout(0, L"Translate: sect=%d m_npts:%d", s, m_nom->m_sections[s]->m_npts);
			}
			m_nom->m_sections[s]->m_skewed = 0;
			m_nom->m_sections[s]->m_skewReport = 0;
			m_nom->m_sections[s]->m_convexSide = 0;
			/*if (m_statusBarHWND != 0 && m_processNomSection[0] != 0)
			{
				wchar_t str[MAXBUFSZ];
				swprintf_s(str, MAXBUFSZ, m_processNomSection, m_nom->m_sections[s]->m_name, s + 1, m_nom->m_numSect);
				::SetWindowText(m_statusBarHWND, str);
			}
			else if (m_statusBarHWND)
			{
				wchar_t str[MAXBUFSZ];
				swprintf_s(str, MAXBUFSZ, L">>>>> %s", m_nom->m_sections[s]->m_name);
				::SetWindowText(m_statusBarHWND, str);
			}*/

			double loff = 0.0, toff = 0.0;
			if (!getline(fin, line))
			{
				ErrorStruct es(BE_INCOMPLETESECTION, s + 1);
				//m_error->AddError(&es);
				return false;
			}
			m_nom->m_sections[s]->m_noseindex[0] = -1;  // these values are used by nominal file editor
			m_nom->m_sections[s]->m_noseindex[1] = -1;
			m_nom->m_sections[s]->m_tailindex[0] = -1;
			m_nom->m_sections[s]->m_tailindex[1] = -1;
			m_nom->m_sections[s]->m_noseforce = 0;
			m_nom->m_sections[s]->m_tailforce = 0;
			m_nom->m_sections[s]->m_letype = EDGE_NORMAL;
			m_nom->m_sections[s]->m_tetype = EDGE_NORMAL;
			for (int end = 0; end < 2; end++)
			{
				int force = 0, ns2, index[3];
				wchar_t typelabel[20], forcelabel[20];
				strline = line.c_str();
				strlinesize = line.length();
				MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, rest, MAXBUFSZ);
				ns = swscanf_s(rest, L"%s %d %d %s %s", key, MAXBUFSZ, &index[0], &index[1], typelabel, 20, forcelabel, 20);
				if (wcscmp(key, L"NOSE") == 0)
				{
					m_nom->m_sections[s]->m_noseindex[0] = index[0] == 0 ? -99 : index[0] - 1;
					if (force)
						m_nom->m_sections[s]->m_noseforce = 1;
					if (ns > 2)
						m_nom->m_sections[s]->m_noseindex[1] = index[1] - 1;
					if (ns > 3)
					{
						MakeUpper(typelabel);
						if (wcscmp(typelabel, L"SQUARE") == 0)
							m_nom->m_sections[s]->m_letype = EDGE_SQUARE;
						else if (wcscmp(typelabel, L"PARTIAL") == 0)
							m_nom->m_sections[s]->m_letype = EDGE_PARTIAL;
					}
				}
				else if (wcscmp(key, L"TAIL") == 0)
				{
					m_nom->m_sections[s]->m_tailindex[0] = index[0] == 0 ? -99 : index[0] - 1;
					if (force)
						m_nom->m_sections[s]->m_tailforce = 1;
					if (ns > 2)
						m_nom->m_sections[s]->m_tailindex[1] = index[1] - 1;
					if (ns > 3)
					{
						MakeUpper(typelabel);
						if (wcscmp(typelabel, L"SQUARE") == 0)
							m_nom->m_sections[s]->m_tetype = EDGE_SQUARE;
						else if (wcscmp(typelabel, L"PARTIAL") == 0)
							m_nom->m_sections[s]->m_tetype = EDGE_PARTIAL;
					}
				}
				else
					break;

				if (!getline(fin, line))
				{
					ErrorStruct es(BE_INCOMPLETESECTION, s + 1);
					//m_error->AddError(&es);
					return false;
				}
				
				
			}
			int skipread = 1;
			if (line.find("OFFSETS") != string::npos || line.find("CHANGES") != string::npos)
			{
				wchar_t typelabel[20], forcelabel[20];
				strline = line.c_str();
				strlinesize = line.length();
				MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, rest, MAXBUFSZ);
				ns = swscanf_s(rest, L"%s %lf %lf", key, MAXBUFSZ, &loff, &toff);
				if (ns < 3)
				{
					ErrorStruct es(BE_BADOFFSETSINNOM, s + 1);
					//m_error->AddError(&es);
					return false;
				}
				skipread = 0;
			}
			else
				loff = toff = 0.0;
			CMatrix val(m_nom->m_sections[s]->m_npts, 8);
			CMatrix planePts(m_nom->m_sections[s]->m_npts, 3);
			CMatrix planeVcs(m_nom->m_sections[s]->m_npts, 3);
			int havek = 0;
			int haveTol = 0;
			for (i = 0; i < m_nom->m_sections[s]->m_npts; i++)
			{
				if (!skipread)
				{
					if (!getline(fin, line))
					{
						ErrorStruct es(BE_MISSINGPOINTNOM, s + 1);
						//m_error->AddError(&es);
						return false;
					}
				}
				skipread = 0;
				strline = line.c_str();
				strlinesize = line.length();
				MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, rest, MAXBUFSZ);
				Replace(rest, wchar_t(','), wchar_t(' '));
				double d[8];
				ns = swscanf_s(rest, L"%lf %lf %lf %lf %lf %lf %lf %lf", &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7]);
				if (ns < 2 || n == 4 || n == 7)
				{
					ErrorStruct es(BE_BADPOINTINNOM, i, s + 1);
					//m_error->AddError(&es);
					return false;
				}
				if (ns == 2)
					d[2] = m_nom->m_sections[s]->m_zval;

				if (ns == 5)  // tolerances, no vectors
				{
					d[6] = d[3];
					d[7] = d[4];
					d[3] = d[4] = d[5] = 0.0;
				}
				if (ns == 3)  // only points
					d[3] = d[4] = d[5] = d[6] = d[7] = 0.0;

				if (ns == 6)  // vector, no tolerances
					d[6] = d[7] = 0.0;

				val.m[i][0] = d[0];
				val.m[i][1] = d[1];
				val.m[i][2] = d[2];
				val.m[i][3] = d[3];
				val.m[i][4] = d[4];
				val.m[i][5] = d[5];
				val.m[i][6] = d[6];
				val.m[i][7] = d[7];
				planePts.m[i][0] = d[0];
				planePts.m[i][1] = d[1];
				planePts.m[i][2] = d[2];

				planeVcs.m[i][0] = d[3];
				planeVcs.m[i][1] = d[4];
				planeVcs.m[i][2] = d[5];
				if (fabs(d[5]) > 1.0e-4)
					havek = 1;
				if (m_readonly)  // BladeRunner needs to save IJK in file if the vectors are defined
					if ((fabs(d[3]) > 1.0e-4) || (fabs(d[3]) > 1.0e-4))
						havek = 1;

				if (fabs(d[6]) > 1.0e-5 || fabs(d[7]) > 1.0e-5)
					haveTol = 1;
			}// for (i = 0; i < m_nom->m_sections[s]->m_npts; i++)
			m_nom->m_sections[s]->m_havek = havek;
			m_nom->m_sections[s]->m_haveTol = haveTol;



			double nose[2][3], tail[2][3];
			for (int eind = 0; eind < 2; eind++)
			{
				if (m_nom->m_sections[s]->m_noseindex[eind] >= 0)
				{
					nose[eind][0] = val.m[m_nom->m_sections[s]->m_noseindex[eind]][0];
					nose[eind][1] = val.m[m_nom->m_sections[s]->m_noseindex[eind]][1];
					nose[eind][2] = val.m[m_nom->m_sections[s]->m_noseindex[eind]][2];
				}
				if (m_nom->m_sections[s]->m_tailindex[eind] >= 0)
				{
					tail[eind][0] = val.m[m_nom->m_sections[s]->m_tailindex[eind]][0];
					tail[eind][1] = val.m[m_nom->m_sections[s]->m_tailindex[eind]][1];
					tail[eind][2] = val.m[m_nom->m_sections[s]->m_tailindex[eind]][2];
				}
			}
			int useextr = 1 * m_nom->m_sections[s]->m_noseforce + 2 * m_nom->m_sections[s]->m_tailforce;

			m_nom->m_sections[s]->m_change[0] = loff;  // these are used by the editor, make Calculate use some day?
			m_nom->m_sections[s]->m_change[1] = toff;
			if (m_nom->m_sections[s]->m_skewed && !m_readonly)
			{
				double pnt[3], vec[3];
				calc_plane(m_nom->m_sections[s]->m_npts, planePts.m, pnt, vec);

				if (fabs(vec[2]) > 1.0e-6)
				{
					double zprime = dot(pnt, vec, 3) / vec[2];
					pnt[0] = pnt[1] = 0.0;
					pnt[2] = zprime;
				}
				m_nom->m_sections[s]->ApplySkew(vec, pnt);

			}// build up the input data-structure for the 2016 mean camber line algorithm
			auto mclParams = computeMCLParams(s, val);
			bugout(0, L"Translate: loff=%lf, toff=%lf", loff, toff);
			bool vecWarn = false;
			int rv;
			m_errorCode = NS_OK;
			if (m_nom->m_sections[s]->m_noseindex[1] >= 0 || m_nom->m_sections[s]->m_tailindex[1] >= 0)
				rv = m_nom->m_sections[s]->Calculate(val.m, m_nom->m_sections[s]->m_noseindex[0],
					m_nom->m_sections[s]->m_noseindex[1], m_nom->m_sections[s]->m_tailindex[0],
					m_nom->m_sections[s]->m_tailindex[1], m_nom->m_sections[s]->m_letype,
					m_nom->m_sections[s]->m_tetype, &loff, &toff, useextr, mclParams.get());
			else
				rv = m_nom->m_sections[s]->Calculate(val.m, m_nom->m_sections[s]->m_noseindex[0] >= 0 ? nose[0] : 0,
					m_nom->m_sections[s]->m_tailindex[0] >= 0 ? tail[0] : 0, &loff, &toff,
					useextr, &vecWarn, mclParams.get());
			if (!rv)
			{
				m_errorCode = m_nom->m_sections[s]->m_errorCode;

				ErrorStruct es(BE_NOCALCSECTIONNOM, m_nom->m_sections[s]->m_name);
				//m_error->AddError(&es);
				return false;
			} if (vecWarn)
			{
				if (*firstWarning < 0)
					*firstWarning = s;
			}
			bugout(0, L"Translate:s=%d m_npts:%d", s, m_nom->m_sections[s]->m_npts);
		}//for (s = 0; s < m_nom->m_numSect; s++)
		/****************读取变公差配置*****************************/
		for (int s1 = 0; s1 < m_nom->m_numSect; s1++)
		{
			if (!getline(fin, line))
			{
				ErrorStruct es(BE_MISSINGPOINTNOM, s + 1);
				//m_error->AddError(&es);
				return false;
			}
			strline = line.c_str();
			strlinesize = line.length();
			MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, rest, MAXBUFSZ);
			int tolsegCount = 0;
			int tol_ccw = 0, ns;
			ns = swscanf_s(rest, L"%s %s %d", key, MAXBUFSZ, tmp_secName, MAXBUFSZ, &tolsegCount);
			if (wcscmp(key, L"PERIODTOL_NUM") == 0)
			{
				if (tol_ccw != 1)
				{
					tol_ccw = 0;
				}
				if (ns == 2)
				{
					tol_ccw = 1;
				}
				bugout(0, L"Translate: read PERIODTOL_NUM:tolsegCount=%d, ccw=%d", tolsegCount, tol_ccw);
				if (tolsegCount > 0)
				{
					m_nom->m_sections[s1]->MakeTolSegStore(tolsegCount);
					m_nom->m_sections[s1]->m_tolCCW = tol_ccw;
				}
				for (int m = 0; m < tolsegCount; m++)
				{
					wchar_t segname[MAXBUFSZ], tolsegName[MAXBUFSZ];
					int start_poi, end_poi, curv_coef = 0;
					int idx = 0;
					double m_pTol_start;
					double m_pTol_end;
					double m_mTol_start;
					double m_mTol_end;
					if (!getline(fin, line))
					{
						ErrorStruct es(BE_MISSINGPOINTNOM, s + 1);
						//m_error->AddError(&es);
						return false;
					}
					strline = line.c_str();
					strlinesize = line.length();
					MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, rest, MAXBUFSZ);
					// PERIODTOL
					ns = swscanf_s(rest, L"%s %d %s %d %d %lf %lf %lf %lf", segname, MAXBUFSZ, &idx, tolsegName, MAXBUFSZ,
						&start_poi, &end_poi, &m_pTol_start, &m_pTol_end, &m_mTol_start, &m_mTol_end);
					if (wcscmp(segname, L"PERIODTOL") == 0)
					{
						//bugout(0, L"Translate: read PERIODTOL_NUM:tolSegName:%s,idx=%d,(%d~%d) U(%f %f) L(%f %f)\n", tolsegName, idx,
						//       start_poi, end_poi, m_pTol_start, m_pTol_end, m_mTol_start, m_mTol_end);
					}
					m_nom->m_sections[s1]->m_start_point[m] = start_poi;
					m_nom->m_sections[s1]->m_end_point[m] = end_poi;
					m_nom->m_sections[s1]->m_pTol_start[m] = m_pTol_start;
					m_nom->m_sections[s1]->m_pTol_end[m] = m_pTol_end;
					m_nom->m_sections[s1]->m_mTol_start[m] = m_mTol_start;
					m_nom->m_sections[s1]->m_mTol_end[m] = m_mTol_end;
					if (curv_coef < 0)
					{
						curv_coef = 0;
					}
					m_nom->m_sections[s1]->m_curvature_coef[m] = curv_coef;
					int iSize =
						WideCharToMultiByte(CP_ACP, 0, tolsegName, -1, NULL, 0, NULL, NULL); // iSize =wcslen(pwsUnicode)+1=6
					char* pszMultiByte;
					pszMultiByte = (char*)malloc(iSize * sizeof(char));
					WideCharToMultiByte(CP_ACP, 0, tolsegName, -1, pszMultiByte, iSize, NULL, NULL);
					//m_nom->m_sections[s1]->m_tolsegNameVec.push_back(pszMultiByte);
					bugout(0, L"Translate: read from :Validate idx=%d,startIndex:%d,endIndex:%d U(%f %f) L(%f %f)\n", idx,
					   m_nom->m_sections[s1]->m_start_point[m], m_nom->m_sections[s1]->m_end_point[m],
					   m_nom->m_sections[s1]->m_pTol_start[m], m_nom->m_sections[s1]->m_pTol_end[m],
					   m_nom->m_sections[s1]->m_mTol_start[m], m_nom->m_sections[s1]->m_mTol_end[m]);
				}
			}
		}

		/****************读取开口区域*****************************/
		for (int s1 = 0; s1 < m_nom->m_numSect; s1++)
		{
			if (!getline(fin, line))
			{
				ErrorStruct es(BE_MISSINGPOINTNOM, s + 1);
				//m_error->AddError(&es);
				return false;
			}
			strline = line.c_str();
			strlinesize = line.length();
			MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, rest, MAXBUFSZ);
			int areaCount = 0;
			int tol_ccw = 0, ns;
			ns = swscanf_s(rest, L"%s %s %d", key, MAXBUFSZ, tmp_secName, MAXBUFSZ, &areaCount);
			if (wcscmp(key, L"OPENAREA_NUM") == 0)
			{
				bugout(0, L"Translate: read OPENAREA_NUM:secName=%s,AreaCount=%d", tmp_secName, areaCount);
				m_nom->m_sections[s1]->m_areaCount = areaCount;
				for (int m = 0; m < areaCount; m++)
				{
					wchar_t segname[MAXBUFSZ], tolsegName[MAXBUFSZ];

					double m_startX;
					double m_startY;
					double m_endX;
					double m_endY;
					int s_index = 0;
					int e_index = 0;
					if (!getline(fin, line))
					{
						ErrorStruct es(BE_MISSINGPOINTNOM, s1 + 1);
						//m_error->AddError(&es);
						return false;
					}
					strline = line.c_str();
					strlinesize = line.length();
					MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, rest, MAXBUFSZ);
					// PERIODTOL
					ns = swscanf_s(rest, L"%s %lf %lf %lf %lf %d %d", segname, MAXBUFSZ, &m_startX, &m_startY, &m_endX, &m_endY, &s_index, &e_index);
					m_nom->m_sections[s1]->m_areaStartPoint[m][0] = m_startX;
					m_nom->m_sections[s1]->m_areaStartPoint[m][1] = m_startY;
					m_nom->m_sections[s1]->m_areaEndPoint[m][0] = m_endX;
					m_nom->m_sections[s1]->m_areaEndPoint[m][1] = m_endY;
					m_nom->m_sections[s1]->m_areaIndex[m][0] = s_index;
					m_nom->m_sections[s1]->m_areaIndex[m][1] = e_index;
					//bugout(0, L"Translate:  OPENAREA_NUM:s-e(%d-%d)", s_index, e_index);
				}
			}

		}
	}
	if (m_readonly)
		return true;
	bugout(0, L"Translate； WriteFile save MATH File ");

	if (!m_nom->WriteFile(m_mathFile))
	{
		ErrorStruct es(BE_NOFILECREATE);
		//m_error->AddError(&es);
		return false;
	}
	return false;
}

void computeEdge(int edgeType, int index0, int index1, const CMatrix& val, int numberOfPoints,
	Eigen::Ref<Eigen::Vector2d> outEdgePoint, Eigen::Ref<Eigen::Vector2d> outEdgeNormal)
{
	switch (edgeType)
	{
	case EDGE_NORMAL:
	{
		//alwaysAssert(index0 >= 0);
		//alwaysAssert(index1 == -1);
		outEdgePoint = Eigen::Vector2d(val.m[index0][0], val.m[index0][1]);
		outEdgeNormal = Eigen::Vector2d(val.m[index0][3], val.m[index0][4]);
		outEdgeNormal.normalize();
		return;
	}
	case EDGE_SQUARE:
	{
		//alwaysAssert(index0 >= 0);
		//alwaysAssert(index1 >= 0);
		int edgeIndex;
		if (abs(index0 - index1) < (numberOfPoints / 2))
		{
			edgeIndex = (index0 + index1) / 2;
		}
		else
		{
			edgeIndex = ((index0 + index1 + numberOfPoints) / 2) % numberOfPoints;
		}
		outEdgePoint = Eigen::Vector2d(val.m[edgeIndex][0], val.m[edgeIndex][1]);
		outEdgeNormal = Eigen::Vector2d(val.m[edgeIndex][3], val.m[edgeIndex][4]);
		outEdgeNormal.normalize();
		return;
	}
	case EDGE_PARTIAL:
	{
		//alwaysAssert(index0 >= 0);
		//alwaysAssert(index1 >= 0);
		Eigen::Vector2d point0(val.m[index0][0], val.m[index0][1]);
		Eigen::Vector2d point1(val.m[index1][0], val.m[index1][1]);
		outEdgePoint = 0.5 * point0 + 0.5 * point1;
		Eigen::Vector2d difference = point0 - point1;
		outEdgeNormal[0] = difference[1];
		outEdgeNormal[1] = -difference[0];
		outEdgeNormal.normalize();
		return;
	}
	default:
		throw std::logic_error("Unknown edge type: neither normal, SQUARE, nor PARTIAL");
	}
}


std::shared_ptr<Hexagon::Blade::MeanCamberCurveParameters2016> CNominalFile::computeMCLParams(int s,const CMatrix& val) const
{
	//return nullptr;
	/*if (m_ptol == nullptr)
	{
		return nullptr;
	}
	int ts = toleranceSectionIndex(m_ptol->m_tol, m_nom->m_sections[s]->m_name);
	if (ts < 0)
	{
		return nullptr;
	}
	
	CToleranceSection* tolerance = m_ptol->m_tol->m_sect[ts];
	if (!std::isfinite(tolerance->m_camber2016_offsets[0]) && !std::isfinite(tolerance->m_camber2016_offsets[1]))
	{
		return nullptr;
	}
	*/
	CNominalSection* nominal = m_nom->m_sections[s];
	auto mclParams = std::make_shared<Hexagon::Blade::MeanCamberCurveParameters2016>();
	mclParams->wholeCurve = nullptr;

	// look at the nose 
	Eigen::Vector2d nosePoint;
	Eigen::Vector2d noseNormal;
	computeEdge(nominal->m_letype, nominal->m_noseindex[0], nominal->m_noseindex[1], val, m_nom->m_sections[s]->m_npts,nosePoint, noseNormal);

	// look at the tail 
	Eigen::Vector2d tailPoint;
	Eigen::Vector2d tailNormal;
	computeEdge(nominal->m_tetype, nominal->m_tailindex[0], nominal->m_tailindex[1], val, m_nom->m_sections[s]->m_npts,
		tailPoint, tailNormal);

	// in some circumstances (partial edges, in particular) it is possible for the normal
	// to point the wrong way
	Eigen::Vector2d tailToNose = nosePoint - tailPoint;
	if (tailToNose.dot(noseNormal) < 0.0)
	{
		noseNormal *= -1.0;
	}
	if (tailToNose.dot(tailNormal) > 0.0)
	{
		tailNormal *= -1.0;
	}

	// set up the half spaces
	//if (nominal->m_letype != EDGE_PARTIAL)
	//{
	//	double noseBackoffDistance = std::isfinite(tolerance->m_camber2016_offsets[0]) ? tolerance->m_camber2016_offsets[0] : 0.0;
	//	mclParams->noseBackoff = std::make_shared<Hexagon::Blade::CamberBackoff>();
	//	mclParams->noseBackoff->point.reset(new Eigen::Array2d(nosePoint));
	//	mclParams->noseBackoff->normal.reset(new Eigen::Array2d(noseNormal));
	//	mclParams->noseBackoff->backoffDistance = noseBackoffDistance;
	//}
	//if (nominal->m_tetype != EDGE_PARTIAL)
	//{
	//	double tailBackoffDistance = std::isfinite(tolerance->m_camber2016_offsets[1]) ? tolerance->m_camber2016_offsets[1] : 0.0;
	//	mclParams->tailBackoff = std::make_shared<Hexagon::Blade::CamberBackoff>();
	//	mclParams->tailBackoff->point.reset(new Eigen::Array2d(tailPoint));
	//	mclParams->tailBackoff->normal.reset(new Eigen::Array2d(tailNormal));
	//	mclParams->tailBackoff->backoffDistance = tailBackoffDistance;
	//}

	// all done
	return mclParams;
}