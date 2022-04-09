

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

#define MAXBUFSZ 1000

DLLEXPORT  unsigned int myGetProfileInt(wchar_t *lpszEntry, int nDefault );
DLLEXPORT void myGetProfileString(wchar_t *lpszEntry, wchar_t *returnValue, wchar_t *defaultValue );
DLLEXPORT double myGetProfileDouble(wchar_t *lpszEntry, double dDefault);
DLLEXPORT bool IsCustomerID(wchar_t *id);

DLLEXPORT void TrimLeft(wchar_t *buf);
DLLEXPORT void TrimRight(wchar_t *buf);
DLLEXPORT void MakeUpper(wchar_t *buf);
DLLEXPORT void Replace(wchar_t *buf, wchar_t oldch, wchar_t newch);

DLLEXPORT int ReverseFind(wchar_t *buf, wchar_t ch);
DLLEXPORT int FindOneOf(wchar_t *buf, wchar_t *delims);
DLLEXPORT void SpanExcluding( wchar_t *outbuf, int len, wchar_t *inbuf, wchar_t *delims);
DLLEXPORT void Left(wchar_t *inbuf, int n, wchar_t *outbuf);
DLLEXPORT void Right(wchar_t *inbuf, int n, wchar_t *outbuf);
DLLEXPORT void Mid(wchar_t *inbuf, int s, wchar_t *outbuf, int n=0);
DLLEXPORT int Find(const wchar_t *buf, const wchar_t *substr, int start=0);

