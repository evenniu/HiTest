#pragma once
#include <string>
#include "AES.h"
#include "stdafx.h"
using namespace std;
namespace MyApp
{
	class __declspec(dllexport) FUNClass
	{
	private:
		wstring		m_sName;
		int			m_iResult;
		bool g_allowOverride;
		AES* m_aes;
		CMatrix* m_noms;

	public:
		FUNClass(void);
		~FUNClass(void);
		int Add(int a, int b, int& sum);
		int Add1(int a, int b);
		int Subtract(int a, int b);

		/* return a * b */
		double Multiply(double a, double b);
		double Divide(double a, double b);
		void TestMatrix(int np);
	};
}
