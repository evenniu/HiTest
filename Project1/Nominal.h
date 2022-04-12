#pragma once
#include "NominalSection.h"
#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

class  DLLEXPORT CNominal
{
public:
	CNominal(void);
	~CNominal(void);

	short m_numSect;
	bool m_english;
	CNominalSection** m_sections;

	int NumSect()
	{
		return m_numSect;
	}

	bool WriteFile(wchar_t* filename);//Write to math  file
};

