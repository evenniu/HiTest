#pragma once
#include "Analysis.h"
#include "FlavorFile.h"

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#define EXTERN_TEMPLATE
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#define EXTERN_TEMPLATE extern
#endif

namespace MyFuncClassApplication
{
	class DLLEXPORT CAnalysisFile
	{
	public:
		CAnalysisFile(void);
		CAnalysisFile(wchar_t* fn, bool aut, HWND statHWND);
		~CAnalysisFile();
		FlavorFile* m_pFlavorFile;//∂¡»°flvŒƒº˛≈‰÷√
		CAnalysis* m_analysis;

		wchar_t m_rest[MAXBUFSZ];
		wchar_t m_buf[MAXBUFSZ];
		wchar_t m_fileName[MAXBUFSZ];  // name of file to open
		wchar_t m_partName[MAXBUFSZ];
		wchar_t m_nomFileName[MAXBUFSZ];
		wchar_t m_tolFileName[MAXBUFSZ];
		wchar_t m_mathFileName[MAXBUFSZ];

		int m_numExtraDimensions;

		bool m_autoSave;
		bool m_good;                      // calcs okay?
		bool m_saved;                     // save performed yet?


		int m_mo, m_da, m_yr, m_hr, m_mn, m_sc; // rpt file date

		int ReadFile(bool& allowOverride);
		int ProcessAnalysisData();
	};
}


