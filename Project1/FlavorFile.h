#pragma once
#include "Flavor.h"

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif
class DLLEXPORT FlavorFile
{
private:
	bool ReadFile();
	wchar_t m_flvfileName[MAXBUFSZ];
public:
	CFlavor* m_flav;
	FlavorFile();
	FlavorFile(wchar_t* flavorFile, wchar_t* rptDir);
	~FlavorFile();
	bool SaveFile();

	bool m_fakeCalc;
};

