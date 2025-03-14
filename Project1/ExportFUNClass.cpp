#include "ExportFUNClass.h"
// 创建对象
extern "C" __declspec(dllexport)  void CreateSimple()
{
	g_pSimple = new FUNClass();
}
extern "C" __declspec(dllexport)  int Add(int x, int y, int& sum)
{
	return	g_pSimple->Add(x, y, sum);
	//return	x+y;
}
extern "C" __declspec(dllexport) int Divide(int a, int b)
{
	return a / b;
}
extern "C" __declspec(dllexport) void TestMatrix(int np)
{
	g_pSimple->TestMatrix(np);
}
extern "C" __declspec(dllexport)  double Multiply(double a, double b)
{
	return g_pSimple->Multiply(a, b);
}
extern "C" __declspec(dllexport)  void SetName(LPCTSTR sName)
{
	//g_pSimple->SetName(sName);
}
/* return a - b */
extern "C" __declspec(dllexport)  double Subtract(int a, int b)
{
	return g_pSimple->Subtract(a, b);
}
extern "C" __declspec(dllexport) void enycode(LPCTSTR plaintext)
{
	//g_pSimple->enycode(plaintext);
}
extern "C" __declspec(dllexport) void Calculate(LPCTSTR plaintext)
{
	g_pSimple->Calculate(plaintext);
}
extern "C" __declspec(dllexport) void Release()
{
	if (NULL != g_pSimple)
	{
		delete g_pSimple;
		g_pSimple = NULL;
	}
}

void GetFunc(int* res, int size)
{
	int numSect = g_pSimple->m_cAnalysisFile->m_analysis->m_numSect;
	for (int s = 0; s < size; s++)
	{
		res[s] = s;
	}
}

void GetSecNum(int* res, int& size)
{
	if (g_pSimple->m_cAnalysisFile)
	{
		size = 5;// g_pSimple->m_cAnalysisFile->m_analysis->m_numSect;
		for (int s = 0; s < size; s++)
		{
			res[s] = s * 2;
		}
	}
}

bool LoadPoint(double** meas, int numpoint)
{
	if (g_pSimple)
	{
		return g_pSimple->ReadData(meas, numpoint, nullptr,0);
	}
	return false;
}

bool CalcBestFit(int BestFitType, int fitToMiddleOfZone, int Transfit, bool noRotate, int rotfit, int useNominal)
{
	if (g_pSimple)
	{
		CAnalysis* analysis = g_pSimple->m_cAnalysisFile->m_analysis;

		//		
		int bfind = 0; 
		double mtols[4] = { -100, -100, -100, -100 };
		double ptols[4] = { 100, 100, 100, 100 };
		BladeBestFitType fitType;
		fitType = (BladeBestFitType)BestFitType;
		
		if (analysis->Locate(0, fitType, fitToMiddleOfZone,Transfit,noRotate,rotfit, useNominal))
		{
#if 0//输出打印调试信息
			CBestFit* thisFit = m_pBlade->m_section[i]->GetBestFitV1(m_pBestFitSection[i][bfind]);
			thisFit->ReportFit(m_pFlavor->m_reportFit[bfind]);

			double x, y, ang, np[2], bp[2];
			thisFit->ReturnFit(&x, &y, &ang);
			bugout(2, L"Locate: after CalcAlign %d, will cal GetBestFitV1{%lf,%lf, %lf}",
				m_pBestFitSection[i][bfind], x, y, ang);
			int numPoints = thisFit->NumPoints();
			//for(int i = 0; i < numPoints; i++)
			//{
			//  np[0] = thisFit->m_noms->m[i][0];
			//  np[1] = thisFit->m_noms->m[i][1];
			//  bp[0] = thisFit->m_infs->m[i][0];
			//  bp[1] = thisFit->m_infs->m[i][1];
			//  bugout(3, L"after CalcAlign  %lf %lf %lf %lf}", bp[0], bp[1], np[0], np[1]);
			//}
#endif 
		}
	}
	return false;
}
