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
		
	}
	FUNClass::~FUNClass(void)
	{
		if (m_aes)
			delete m_aes;
		if (m_noms)
			delete m_noms;
	}
	int FUNClass::Add(int a, int b, int& sum)
	{
		m_iResult = m_aes->add(a, b);
		sum = m_iResult;
		return  m_iResult;
	}
	int FUNClass::Add1(int a, int b)
	{
		return a + b;
	}
	int FUNClass::Subtract(int a, int b) {
		return a - b;
	}

	/* return a * b */
	double FUNClass::Multiply(double a, double b)
	{
		return a * b + 1;
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
	}
}


