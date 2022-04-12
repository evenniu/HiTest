#include "stdafx.h"
#include <fstream>
//#include "Plane.h"
#include "NominalFile.h"
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
	bugout(0, L"new CNominalFile will cal Translate ");
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
			//continue;
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
				//calc_plane(m_nom->m_sections[s]->m_npts, planePts.m, pnt, vec);

				if (fabs(vec[2]) > 1.0e-6)
				{
					double zprime = dot(pnt, vec, 3) / vec[2];
					pnt[0] = pnt[1] = 0.0;
					pnt[2] = zprime;
				}
			}
		}//for (s = 0; s < m_nom->m_numSect; s++)
		for (s = 0; s < m_nom->m_numSect; s++)
		{

		}
	}
	return false;
}
