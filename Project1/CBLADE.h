#pragma once
#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

//class CViewMatrixComponent;
class CAnalysis;


class DLLEXPORT CBlade
{
private:
	short m_numSections;
	bool m_english;
	bool m_nomValid;
	bool m_measLoaded;
	wchar_t m_mathFileName[MAXBUFSZ];
	//CBestFit* m_pBestFit;
	//CBestFit* m_pBestFitLE;
public:

	bool IsValid();
	int NumSect()
	{
		return m_numSections;
	}
	bool IsEnglish()
	{
		return m_english;
	}
	void MeasLoaded(bool loaded)
	{
		m_measLoaded = loaded;
	}
	bool MeasLoaded()
	{
		return m_measLoaded;
	}
	//CBestFit* GetBestFit()
	//{
	//	return m_pBestFit;
	//}

	//CBestFit* GetBestFitLE()
	//{
	//	return m_pBestFitLE;
	//}
	//bool FitBlade(const CFitParams* fp, int sec1, int sec2);
	//bool FitBladeLELS(const CFitParams* fp, CAnalysis* pAnal);

	CBlade();
	CBlade(wchar_t* mathFileName);
	virtual ~CBlade();

	const wchar_t* GetFileName()
	{
		return m_mathFileName;
	}

	bool ReadFile(FILE* fp);
};

