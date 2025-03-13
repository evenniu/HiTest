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
		return g_pSimple->ReadData(meas, numpoint);
	}
	return false;
}
