#include "StdAfx.h"
#include "ANAFILE.h"
#include <cstring>
#include <fstream>
#include <string>
namespace MyFuncClassApplication
{
	int CAnalysisFile::ReadFile(bool& allowOverride)
	{
		wchar_t rptPath[MAXBUFSZ];
		wchar_t flvName[MAXBUFSZ];
		wchar_t m_rest[MAXBUFSZ] = L"E:\\BladeExap\\219Template\\A-A.rpt";
		wcscpy_s(m_fileName, MAXBUFSZ, L"E:\\BladeExap\\219Template\\A-A.rpt");
	
		wcscpy_s(rptPath, MAXBUFSZ, m_fileName);

		int pos = ReverseFind(rptPath, wchar_t('\\'));
		if (pos < 0)
			pos = ReverseFind(rptPath, wchar_t('/'));
		Left(rptPath, pos, rptPath);
		std::ifstream  fin(m_rest, std::ios::in);
		if (!fin)
		{

		}
		string line;
		std::string  x = "";
		std::string  y = "";
		std::string  z = "";
		wchar_t textINfo[MAXBUFSZ];
		char charStrs[MAXBUFSZ];
		char* pCStrKey;
		string::size_type idx;
		const char* strline;
		int strlinesize = 0;
		while (getline(fin, line)) // line中不包括每行的换行符
		{
			idx = line.find("FLAVOR");
			if (idx != string::npos)
			{
				line = line.substr(7);
				strline = line.c_str();
				strlinesize = line.length();
				MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, flvName, MAXBUFSZ);
				//ewcscpy_s(flvName,)
				//swprintf(flvName, MAXBUFSZ,);
				printf("idx=%d", idx);
				continue;
			}
			idx = line.find("NOMINAL");
			if (idx != string::npos)
			{
				line = line.substr(8);
				strline = line.c_str();
				strlinesize = line.length();
				MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, m_nomFileName, MAXBUFSZ);
				continue;
			}
			idx = line.find("TOLERANCE");
			if (idx != string::npos)
			{
				line = line.substr(11);
				strline = line.c_str();
				strlinesize = line.length();
				MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, m_tolFileName, MAXBUFSZ);
				continue;
			}
			idx = line.find("PART");
			if (idx != string::npos)
			{
				line = line.substr(5);
				strline = line.c_str();
				strlinesize = line.length();
				MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, m_partName, MAXBUFSZ);
				continue;
			}
			idx = line.find("DATE");
			if (idx != string::npos)
			{
				line = line.substr(5);
				strline = line.c_str();
				strlinesize = line.length();
				MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, m_rest, MAXBUFSZ);
				ptrdiff_t dateLen = wcslen(m_rest);
				for (int c = 0; c < dateLen; c++)
				{
					wchar_t thisChar = m_rest[c];
					if (thisChar == wchar_t(' ') || iswdigit(thisChar))
						continue;

					// something other than a space or a decimal digit, map to space

					Replace(m_rest, thisChar, wchar_t(' '));
				}
				if (swscanf_s(m_rest, L"%d %d %d %d %d %d", &m_mo, &m_da, &m_yr, &m_hr, &m_mn, &m_sc) < 6)
				{
					return 0;
				}
				if (m_yr < 70)
					m_yr += 2000;
				else if (m_yr < 1900)
					m_yr += 1900;
				continue;
			}

			idx = line.find("RADIUS");
			if (idx != string::npos)
			{
				line = line.substr(7);
				strline = line.c_str();
				strlinesize = line.length();
				MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, m_rest, MAXBUFSZ);
				if (swscanf_s(m_rest, L"%lf", &m_analysis->m_probeRad) < 1)
				{
					return 0;
				}
				continue;
			}
			idx = line.find("NUM_SECT");
			if (idx != string::npos)
			{
				line = line.substr(9);
				strline = line.c_str();
				strlinesize = line.length();
				MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, m_rest, MAXBUFSZ);
				if (swscanf_s(m_rest, L"%d", &m_analysis->m_numSect) < 1)
					m_analysis->m_numSect = 0;

				if (m_analysis->m_numSect < 1)
					continue;
				m_analysis->m_sect = new CAnalysisSect[m_analysis->m_numSect];
				for (int i = 0; i < m_analysis->m_numSect; i++)
				{
					getline(fin, line);
					strline = line.c_str();
					strlinesize = line.length();
					MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, m_rest, MAXBUFSZ);
					wchar_t dummy[80], sectname[80];
					int numberOfElementsReadFromSectionHeading, n, n1, n2;

					numberOfElementsReadFromSectionHeading = swscanf_s(m_rest, L"%s %s %d %d %d", dummy, 80, sectname, 80, &n, &n1, &n2);
					if (numberOfElementsReadFromSectionHeading < 3)
					{
						return 0;
					}
					m_analysis->m_sect[i].m_numPoints = n;
					m_analysis->m_sect[i].numberOfBallCenters = n;
					m_analysis->m_sect[i].x = new double[n];
					m_analysis->m_sect[i].y = new double[n];
					m_analysis->m_sect[i].z = new double[n];
					m_analysis->m_sect[i].ox = new double[n];
					m_analysis->m_sect[i].oy = new double[n];
					m_analysis->m_sect[i].oz = new double[n];
					m_analysis->m_sect[i].i = new double[n];
					m_analysis->m_sect[i].j = new double[n];
					m_analysis->m_sect[i].oi = new double[n];
					m_analysis->m_sect[i].oj = new double[n];
					m_analysis->m_sect[i].ballCenterX = new double[n];
					m_analysis->m_sect[i].ballCenterY = new double[n];
					m_analysis->m_sect[i].ballCenterZ = new double[n];
					wcscpy_s(m_analysis->m_sect[i].m_sectName, sectname);
					MakeUpper(m_analysis->m_sect[i].m_sectName);

					int s;
					//for (s = 0; s < m_analysis->m_pBlade->NumSect(); s++)   // loop thru blade sections
					//	if (wcscmp(m_analysis->m_sect[i].m_sectName, m_analysis->m_pBlade->m_section[s]->Name()) == 0)
					//		break;
				/*	if (s == m_analysis->m_pBlade->NumSect())   // did not find it
					{
						ErrorStruct es(BE_BADSECTION, m_analysis->m_sect[i].m_sectName);
						m_error->AddError(&es);
						while (head)
						{
							tail = head;
							head = head->Next;
							delete tail;
						}
						return 0;
					}
					double znom = m_analysis->m_pBlade->m_section[s]->ZValue();
					double ztol = m_analysis->m_pBlade->IsEnglish() ? 0.5 : 12.0;
					m_analysis->m_sect[i].skewalign = m_analysis->m_pBlade->m_section[s]->SkewAlign();
					m_analysis->m_sect[i].m_skewReport = m_analysis->m_pBlade->m_section[s]->SkewReport();

					m_analysis->UpdateStatus(sectname, m_analysis->m_statusReadingPoints);
		*/

					int j;
					int tossed = 0;
					for (j = 0; j < n; j++)
					{
						getline(fin, line);
						strline = line.c_str();
						strlinesize = line.length();
						MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, m_rest, MAXBUFSZ);
						int ns = swscanf_s(m_rest, L"%lf %lf %lf %lf %lf", &m_analysis->m_sect[i].x[j], &m_analysis->m_sect[i].y[j], &m_analysis->m_sect[i].z[j], &m_analysis->m_sect[i].i[j], &m_analysis->m_sect[i].j[j]);
						if (ns < 3)
						{
							return 0;
						}
						if (ns < 5)
						{
							m_analysis->m_sect[i].m_goodVectors = false;
							m_analysis->m_sect[i].i[j] = m_analysis->m_sect[i].j[j] = 0.0;
						}
						/*	if (m_analysis->m_sect[i].skewalign)
							{
								double xyz[3];
								xyz[0] = m_analysis->m_sect[i].x[j];
								xyz[1] = m_analysis->m_sect[i].y[j];
								xyz[2] = m_analysis->m_sect[i].z[j];
								//m_analysis->m_sect[i].skewalign->MeasToBest(xyz, 1, xyz);  // don't rotate vector (don't have a K anyway??)

								m_analysis->m_sect[i].x[j] = xyz[0];
								m_analysis->m_sect[i].y[j] = xyz[1];
								m_analysis->m_sect[i].z[j] = xyz[2];
							}
							*/
							// throw away wild points, check the z height

					   /*	if (fabs(m_analysis->m_sect[i].z[j] - znom) > ztol)
						   {
							   bugout(0, L"Discarded point with bad Z:\n %s", m_buf);
							   if (numberOfElementsReadFromSectionHeading > 4)
							   {
								   if (j < n1)
									   n1--;
								   else
									   n2--;
							   }

							   tossed++;
							   j--;
							   n--;
						   }
						   */
					}
					if (tossed > n / 5)
					{
						/*ErrorStruct es(BE_ZMISMATCH, m_analysis->m_sect[i].m_sectName);
						m_error->AddError(&es);
						while (head)
						{
							tail = head;
							head = head->Next;
							delete tail;
						}*/
						return 0;
					}
					//m_analysis->UpdateStatus(sectname, m_analysis->m_statusTrimAndOrder);
					int types = 0;
					/*	int letype = m_analysis->m_pBlade->m_section[s]->LEType();
					int tetype = m_analysis->m_pBlade->m_section[s]->TEType();

					if (letype == EDGE_SQUARE)
						types |= 2;
					else if (letype == EDGE_PARTIAL)
						types |= 4;

					if (tetype == EDGE_SQUARE)
						types |= 8;
					else if (tetype == EDGE_PARTIAL)
						types |= 16;
					*/
					if (numberOfElementsReadFromSectionHeading < 5)
					{
						n = m_analysis->ReOrder(m_analysis->m_sect[i], types);
						m_analysis->m_sect[i].m_inter1 = m_analysis->m_sect[i].m_inter2 = 0;
					}
					else
					{
						//bugout(3, "n1 %d n2 %d", n1, n2);
						n = m_analysis->ReOrder(m_analysis->m_sect[i], &n1, &n2, types);
						//bugout(3, "n %d (tossed %d =%d+%d-%d)", n, n1+n2-n,n1,n2,n);
					}
					if (n < 10)
					{
						//nointer++;
						//ErrorStruct es(BE_POINTORDERINGFAILEDONKNOWNSECTION, m_analysis->m_sect[i].m_sectName);  // Does this used to be concatenated to existing error?
						//m_error->AddError(&es);
						continue;
					}
				}
				continue;
			}
		}
		fin.clear();
		fin.close();
		m_pFlavorFile = new FlavorFile(flvName, rptPath);
		m_analysis->m_pFlavor = m_pFlavorFile->m_flav;


	}
	int CAnalysisFile::ProcessAnalysisData()
	{
		return m_analysis->FitSplines();
	}
	CAnalysisFile::CAnalysisFile(void)
	{
		wchar_t rptPath[MAXBUFSZ];
		m_fileName[0] = 0;
		//wcscpy_s(rptPath, m_fileName);
		//m_pFlavorFile = NULL;
		m_analysis = new CAnalysis();
		m_buf[0] = 0;
		m_autoSave = false;
		m_saved = false;

	}
	CAnalysisFile::CAnalysisFile(wchar_t* fn, bool aut, HWND statHWND)
	{
		//wchar_t rptPath[MAXBUFSZ];
		//wcscpy_s(rptPath, m_fileName);
		//m_pFlavorFile = NULL;
		m_analysis = new CAnalysis();

		m_autoSave = false;
		m_saved = false;
	}
	CAnalysisFile::~CAnalysisFile()
	{
		//if (m_pFlavorFile)
		//	delete m_pFlavorFile;
		if (m_analysis)
			delete m_analysis;
	}
}
