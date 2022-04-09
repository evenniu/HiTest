#include "stdafx.h"
#include <wchar.h>

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

unsigned int myGetProfileInt(wchar_t *lpszEntry, int nDefault)
{
  int rv = nDefault;

  HKEY hk;
	DWORD dwType;
	DWORD dwLength;
	long ival;
	if (RegOpenKeyEx(HKEY_CURRENT_USER, L"Software\\WAI\\Blade\\Blade", 0, KEY_QUERY_VALUE, &hk) ==  ERROR_SUCCESS)
  {
    dwLength = sizeof(ival);
    LSTATUS st = RegQueryValueEx(hk, lpszEntry, NULL, &dwType, (LPBYTE)(&ival), &dwLength);
	  if (st == ERROR_SUCCESS)
	  {
      rv = ival;
	  }
    RegCloseKey(hk);
  }

  return rv;
}

double myGetProfileDouble(wchar_t *lpszEntry, double dDefault)
{
  double rv = dDefault;

  HKEY hk;
	DWORD dwType;
	DWORD dwLength;
	TCHAR buffer[1024];
	if (RegOpenKeyEx(HKEY_CURRENT_USER, L"Software\\WAI\\Blade\\Blade", 0, KEY_QUERY_VALUE, &hk) ==  ERROR_SUCCESS)
  {
    dwLength = sizeof(buffer);
	  LSTATUS st = RegQueryValueEx(hk, lpszEntry, NULL, &dwType, (LPBYTE)buffer, &dwLength);
    if (st == ERROR_SUCCESS)
	  {
	    if (swscanf_s(buffer, L"%lf", &rv) < 1)
        rv = dDefault;
	  }
    RegCloseKey(hk);
  }
  return rv;
}

void myGetProfileString(wchar_t *lpszEntry, wchar_t *returnValue, wchar_t *defaultValue )
{
  wcscpy_s(returnValue, MAXBUFSZ, defaultValue);

  HKEY hk;
	DWORD dwType;
	DWORD dwLength;
  wchar_t buffer[MAXBUFSZ];
	if (RegOpenKeyEx(HKEY_CURRENT_USER, L"Software\\WAI\\Blade\\Blade", 0, KEY_QUERY_VALUE, &hk) ==  ERROR_SUCCESS)
  {
    dwLength = MAXBUFSZ;
	  if (RegQueryValueEx(hk, lpszEntry, NULL, &dwType, (LPBYTE)buffer, &dwLength) == ERROR_SUCCESS)
      wcscpy_s(returnValue, MAXBUFSZ, buffer);

    RegCloseKey(hk);
  }
}

//bool IsCustomerID(wchar_t *id)
//{
//  wchar_t  customer[MAXBUFSZ];
//  myGetProfileString(L"CustomerID", customer, L"");
//  if (wcscmp(customer, id) == 0)
//    return true;
//  return false;
//}


// some string functions, move these to a better location later so usable by all

// inbuf and outbuf can be the same variable
void TrimLeft(wchar_t *buf)
{
  int i;
  for(i=0; buf[i]; i++)
    if (!iswspace(buf[i]))
      break;

  if (i == 0) // nothing to trim
    return;

  int j = i;
  for(i=j; buf[i]; i++)
    buf[i-j] = buf[i];
  buf[i-j] = 0;
}

void TrimRight(wchar_t *buf)
{
  ptrdiff_t l = wcslen(buf);

  ptrdiff_t i;
  for(i=l-1; i>0; i--)
    if (!iswspace(buf[i]))
      break;
  if (i < l-1)
    buf[i+1] = 0;
}

void Left(wchar_t *inbuf, int n, wchar_t *outbuf)
{
  int i;
  for(i=0; i<n; i++)
    outbuf[i] = inbuf[i];
  outbuf[n] = 0;
}

void Mid(wchar_t *inbuf, int s, wchar_t *outbuf, int n)
{
  ptrdiff_t l = wcslen(inbuf);
  
  if (s < 0 || s >= l)
  {
    outbuf[0] = 0;
    return;
  }

  if (n == 0)
    n = static_cast<int>(l - s);


  int i;
  for(i=0; i<n; i++)
    outbuf[i] = inbuf[s+i];
  outbuf[n] = 0;
}

void Right(wchar_t *inbuf, int n, wchar_t *outbuf)
{
  ptrdiff_t l = wcslen(inbuf);
  ptrdiff_t s = l - n;
  if (s < 0)
    s = 0;

  int i;
  for (i=0; i<n; i++)
    outbuf[i] = inbuf[i+s];
  outbuf[n] = 0;
}

extern void Replace(wchar_t *buf, wchar_t oldch, wchar_t newch)
{
  int i;
  for(i=0; buf[i]; i++)
    if (buf[i] == oldch)
      buf[i] = newch;
}

int Find(const wchar_t *buf, const wchar_t *substr, int start)
{
  ptrdiff_t l = wcslen(substr);
  ptrdiff_t max = wcslen(buf) - l + 1;
  int i;
  for (i=start; i<max; i++)
    if (wcsncmp(&buf[i], substr, l) == 0)
      break;

  if (i < max)
    return i;

  return -1;
}

void MakeUpper(wchar_t *buf)
{
  int i;
  for (i=0; buf[i]; i++)
    buf[i] = towupper(buf[i]);
}

int ReverseFind(wchar_t *buf, wchar_t ch)
{
  ptrdiff_t i = wcslen(buf) - 1;

  while(i >= 0)
  {
    if (buf[i] == ch)
      break;
    i--;
  }

  return static_cast<int>(i);
}

int FindOneOf(wchar_t *buf, wchar_t *delims)
{
  if ( buf[0] < 32)
  {
    int i = 0;
    i++;
  }

  ptrdiff_t li = wcslen(buf);
  int i;
  for(i=0; i<li; i++)
  {
    ptrdiff_t ld = wcslen(delims);
    int j;
    for(j=0; j<ld; j++)
      if (delims[j] == buf[i])
        break;

    if (j < ld)
      return i;
  }

  return -1;
}

void SpanExcluding(wchar_t *outbuf, int len, wchar_t *inbuf, wchar_t *delims)
{
  wcscpy_s(outbuf, len, inbuf);
  int index = FindOneOf(inbuf, delims);
  if (index >= 0)  // no delims, copy entire string
    outbuf[index] = 0;
}
