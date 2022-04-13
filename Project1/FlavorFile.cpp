#include "FlavorFile.h"
#include "stdafx.h"
#include <fstream>
#include <iostream>
#include <string>


FlavorFile::FlavorFile(wchar_t* flavorFile, wchar_t* rptDir)
{
	
    wcscpy_s(m_flvfileName, MAXBUFSZ, flavorFile);
    wcscpy_s(m_rptPath, MAXBUFSZ, rptDir);
	m_flav = new CFlavor();
	m_fakeCalc = false;
	bool noError = true;
	noError = ReadFile();
}
FlavorFile::FlavorFile()
{
	m_flav = new CFlavor();
	m_fakeCalc = false;
	bool noError = true;
	noError = ReadFile();
}
FlavorFile::~FlavorFile()
{
	if (m_flav)
		delete m_flav;
}
bool FlavorFile::ReadFile()
{
	ifstream myfile(m_flvfileName, std::ios::in);
	if (!myfile)
	{
		myfile.close();
		return false;
	}
	string line;
	const char* strline;
	int strlinesize = 0;
	wchar_t key[MAXBUFSZ], rest[MAXBUFSZ];
	int id, type;
	while (getline(myfile, line)) // line中不包括每行的换行符
	{
		if (line.find("CALCS") != string::npos)
		{
			m_flav->m_numCalc = 10;
			continue;
		}
		if (line.find("FORGED") != string::npos)
		{
			continue;
		}
		if (line.find("FULLBLADE") != string::npos)
		{
			continue;
		}
		if (line.find("WEIGHTS") != string::npos)
		{
			continue;
		}		
		if (line.find("USENOMINALS") != string::npos)
		{
			continue;
		}
		strline = line.c_str();
		strlinesize = line.length();
		MultiByteToWideChar(CP_ACP, 0, strline, strlen(strline) + 1, rest, MAXBUFSZ);
	}

	for (int i = 0; i < m_flav->m_numCalc; i++)
	{
		m_flav->m_calcs[i] = (BladeCalculations)(i+2);
	}
	m_flav->m_compensateMethod = 1;
	m_flav->m_endMag = 1;
	m_flav->m_sideMag = 1;
	swprintf_s(rest, MAXBUFSZ, L"E:\\BladeExap2\\Tenon333");
	wcscpy_s(m_flav->m_DPDir, rest);
	return true;
}
bool FlavorFile::SaveFile()
{
	ofstream fout("test.txt", std::ios::out);
	string line;
	if (!fout) {
		fout.close(); //程序结束前不能忘记关闭以前打开过的文件
		cout << "error opening destination file." << endl;
		return 0;
	}

	fout << "123" << "你是好孩子" << endl;
	fout << "第二次写文件" << endl;
	fout.clear();
	fout.close();
	return false;
}
