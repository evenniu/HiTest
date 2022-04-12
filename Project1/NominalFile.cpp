#include "stdafx.h"
#include <fstream>

#include "NominalFile.h"
#include <string>

CNominalFile::CNominalFile(wchar_t* nominalFile, wchar_t* rptPath, wchar_t* mathFile, HWND statusHWND, wchar_t* processNomSection, bool isCreating)
{
	//m_nom = NULL;
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

	//m_nom = new CNominal();

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
	//m_nom->m_numSect = 0;
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
		idx = line.find("MM");
		if (idx != string::npos)
		{
			ok = 1;
			if (wcsncmp(rest, L"MM", 2) == 0)
			{
			}
			else if (wcsncmp(rest, L"IN", 2) == 0)
			{

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
		ok = 0;
		idx = line.find("NUM_SECT");
		if (idx != string::npos)
		{
			line = line.substr(9);
			strline = line.c_str();
			strlinesize = line.length();
			MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, rest, MAXBUFSZ);
		}
	}
	return false;
}
