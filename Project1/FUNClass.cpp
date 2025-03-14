#include "FUNClass.h"
#include <stdexcept>

using namespace std;
namespace MyApp
{
	unsigned char key[] = "he0xsoft93";
	FUNClass::FUNClass(void)
	{
		m_iResult = -1;
		m_sName = L"";
		m_aes = new AES(key);
		m_cAnalysisFile = new CAnalysisFile();
		
	}
	FUNClass::~FUNClass(void)
	{
		if (m_aes)
			delete m_aes;
		if (m_noms)
			delete m_noms;
		if (m_cAnalysisFile)
			delete m_cAnalysisFile;
	}
	int FUNClass::Add(int a, int b, int& sum)
	{
		m_iResult = m_aes->add(a, b);
		sum = m_iResult;
		m_aes->TestEigen();
		return  m_iResult;
	}
	int FUNClass::Add1(int a, int b)
	{
		return a + b;
	}
	int FUNClass::Subtract(int a, int b) {
		return a - b;
	}

	void FUNClass::GetFun(int* res, int size)
	{
		for (int i = 0;i < size; i++)
		{
			res[i] = i * 3;
		}
		// 在这里进行计算，并将结果存储在 data 数组中
	}

	/* return a * b */
	double FUNClass::Multiply(double a, double b)
	{
		return a * b;
	}
	double FUNClass::Divide(double a, double b)
	{
		return a / b;
	}
	void FUNClass::TestMatrix(int np)
	{
		if (np < 0)
			return;
		m_noms = new CMatrix(np, 2);
		for (int i = 0;i < np;i++)
		{
			m_noms->m[i][0] = i;
			m_noms->m[i][1] = i*i;
		}
		//const Eigen::Matrix2Xd measuredPoints = constructPointMatrix(m_mxpt, m_mypt, m_totalPoints);
	}
	void FUNClass::Calculate(wstring plaintext)
	{
		std::string str;
		int nLen = (int)plaintext.length();
		int readFileReturnValue = 0;

		try
		{
			readFileReturnValue = m_cAnalysisFile->ReadFile(g_allowOverride);
			if (readFileReturnValue == false)
			{
				return;
			}

			readFileReturnValue = m_cAnalysisFile->ProcessAnalysisData();//FitSplines
		}
		catch (const std::exception& e)
		{
			//const CString errorMessage = CString(L"Unhandled exception during reading input files: ") + e.what();
		   // AfxMessageBox(errorMessage);
			return;
		}
	}
	bool FUNClass::ReadData(double* measxyzijk[8], int numpoints, double* nomxyzijk[8], int numNomPt)
	{
		if (numpoints <5)
		{
			return false;
		}  
		
		CAnalysis* analysis = m_cAnalysisFile->m_analysis;
		analysis->m_pBlade = new CBlade();
		
		int n = numpoints;
		m_cAnalysisFile->m_pFlavorFile = new FlavorFile();//CFlavorFile
		analysis->m_pFlavor = m_cAnalysisFile->m_pFlavorFile->m_flav;
		analysis->m_numSect = 1;
		int id = 0;
		analysis->m_sect = new CAnalysisSect[m_cAnalysisFile->m_analysis->m_numSect];
		analysis->m_sect[id].m_numPoints = numpoints;
		analysis->m_sect[id].numberOfBallCenters = numpoints;
		analysis->m_sect[id].x = new double[n];
		analysis->m_sect[id].y = new double[n];
		analysis->m_sect[id].z = new double[n];
		analysis->m_sect[id].ox = new double[n];
		analysis->m_sect[id].oy = new double[n];
		analysis->m_sect[id].oz = new double[n];
		analysis->m_sect[id].i = new double[n];
		analysis->m_sect[id].j = new double[n];
		analysis->m_sect[id].oi = new double[n];
		analysis->m_sect[id].oj = new double[n];
		analysis->m_sect[id].ballCenterX = new double[n];
		analysis->m_sect[id].ballCenterY = new double[n];
		analysis->m_sect[id].ballCenterZ = new double[n];
		double* kval = new double[n];
		for (int i = 0; i < numpoints; i++)
		{
			analysis->m_sect[id].x[i] = measxyzijk[i][0];
			analysis->m_sect[id].y[i] = measxyzijk[i][1];
			analysis->m_sect[id].z[i] = measxyzijk[i][2];
			analysis->m_sect[id].ox[i] = measxyzijk[i][0];
			analysis->m_sect[id].oy[i] = measxyzijk[i][1];
			analysis->m_sect[id].oz[i] = measxyzijk[i][2];
			analysis->m_sect[id].i[i] = measxyzijk[i][3];
			analysis->m_sect[id].j[i] = measxyzijk[i][4];
			kval[i] = 0;
		}       
		
		analysis->m_pBlade->ReadNomdata(numpoints, analysis->m_sect[id].x, analysis->m_sect[id].y,kval);
		if (kval)
		{
			delete[] kval;
		}
		if (!analysis->FitSplines())
		{
			return false;
		}
	}
}


