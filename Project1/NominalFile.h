#pragma once
#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif


#include "NominalSection.h"
#include "Nominal.h"


class DLLEXPORT CNominalFile
{
private:
	//CToleranceFile* m_ptol;
	wchar_t m_mathFile[MAXBUFSZ];
	wchar_t m_nomFile[MAXBUFSZ];
	bool m_readonly;
public:
	CNominalFile(wchar_t* nominalFile, wchar_t* rptPath, wchar_t* mathFile, HWND statusHWDN, wchar_t* processNomSection, bool isCreating);
	virtual ~CNominalFile();
	CNominal* m_nom;

	bool m_isCreating;  // used when file is edited or created, to know whether to display file missing error
	int m_errorCode;
	int m_firstWarning;
	HWND m_statusBarHWND;  // window to set status messages
	wchar_t m_processNomSection[MAXBUFSZ];


	bool Translate(int* firstWarning);
	void PreProcessFile();  // convert UA to NOM if necessary

	const wchar_t* GetMathFile()
	{
		return m_mathFile;
	}

	//void UpdateStatus(wchar_t* sectName);
	std::shared_ptr<Hexagon::Blade::MeanCamberCurveParameters2016> computeMCLParams(int s, const CMatrix& val) const;

};

