#include "stdafx.h"
#include "BladeError.h"

#ifdef _DEBUG
   #ifndef DBG_NEW
      #define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
      #define new DBG_NEW
   #endif
#endif  // _DEBUG

void ErrorStruct::Init()
{
  m_errorCode = BE_OK;
  m_subError = NULL;
  m_lineNumber = 0;
  m_intVal[0] = m_intVal[1] = m_intVal[2] = 0;
  m_doubleVal[0] = m_doubleVal[1] = m_doubleVal[2] = 0.0;
  m_stringVal[0][0] = m_stringVal[1][0] = m_stringVal[2][0] = 0;
  m_argType = ArgsNoError;
}

CBladeError::CBladeError()
{
  ClearErrors();
}

void CBladeError::ClearErrors()
{
  m_errors.clear();
//   m_lastError = -1;
}

CBladeError::~CBladeError()
{
}

bool CBladeError::LastError(ErrorStruct *err)
{
  if(m_errors.size() > 0)
  {
    return getError(m_errors.size() - 1, err);
  }

  return false;
}

bool CBladeError::isOK()
{
  return m_errors.empty();
}


void CBladeError::AddError(ErrorStruct *err)
{
  m_errors.resize(m_errors.size()+1);

  err->CopyTo(&m_errors.back());
}

DLLEXPORT std::size_t CBladeError::numberOfErrors() const
{
  return m_errors.size();
}

DLLEXPORT bool CBladeError::getError(std::size_t index, ErrorStruct *err) const
{
  if(index < m_errors.size())
  {
    m_errors.at(index).CopyTo(err);
    return true;
  }

  return false;
}


ErrorStruct::ErrorStruct()
{
  Init();
}

void ErrorStruct::CopyTo(ErrorStruct *err) const
{
  err->m_errorCode = m_errorCode;
  err->m_subError = NULL;
  err->m_lineNumber = m_lineNumber;
  err->m_argType = m_argType;

  for (int i=0; i<3; i++)
  {
    err->m_intVal[i] = m_intVal[i];
    err->m_doubleVal[i] = m_doubleVal[i];
    wcscpy_s(err->m_stringVal[i], m_stringVal[i]);
  }

  if (m_subError != NULL)
  {
    m_subError->m_subError = NULL;  // these shouldn't be nested, make sure recursion doesn't go beyond this one copy.
    err->m_subError = new ErrorStruct();
    m_subError->CopyTo(err->m_subError);
  }
}

ErrorStruct::ErrorStruct(BladeErrorType errnum)
{
  Init();
  m_errorCode = errnum;
  m_subError = NULL;
  m_argType = ArgsNone;
}

ErrorStruct::ErrorStruct(BladeErrorType errnum, const wchar_t *s1)
{
  Init();
  m_errorCode = errnum;
  m_subError = NULL;
  m_argType = ArgsS;
  wcscpy_s(m_stringVal[0], s1);
}

ErrorStruct::ErrorStruct(BladeErrorType errnum, int d1)
{
  Init();
  m_errorCode = errnum;
  m_subError = NULL;
  m_argType = ArgsD;
  m_intVal[0] = d1;
}

ErrorStruct::ErrorStruct(BladeErrorType errnum, int d1, int d2)
{
  Init();
  m_errorCode = errnum;
  m_subError = NULL;
  m_argType = ArgsDD;
  m_intVal[0] = d1;
  m_intVal[1] = d2;
}

ErrorStruct::ErrorStruct(BladeErrorType errnum, const wchar_t *s1, const wchar_t *s2)
{
  Init();
  m_errorCode = errnum;
  m_subError = NULL;
  m_argType = ArgsSS;
  wcscpy_s(m_stringVal[0], s1);
  wcscpy_s(m_stringVal[1], s2);
}

ErrorStruct::ErrorStruct(BladeErrorType errnum, int d1, int d2, const wchar_t *s3)
{
  Init();
  m_errorCode = errnum;
  m_subError = NULL;
  m_argType = ArgsDDS;
  m_intVal[0] = d1;
  m_intVal[1] = d2;
  wcscpy_s(m_stringVal[2], s3);
}

ErrorStruct::ErrorStruct(BladeErrorType errnum, const wchar_t *s1, int d2, int d3)
{
  Init();
  m_errorCode = errnum;
  m_subError = NULL;
  m_argType = ArgsSDD;
  wcscpy_s(m_stringVal[0], s1);
  m_intVal[1] = d2;
  m_intVal[2] = d3;
}

ErrorStruct::ErrorStruct(BladeErrorType errnum, ErrorStruct *subError, const wchar_t *s1)
{
  Init();
  m_errorCode = errnum;
  m_argType = ArgsSubS;
  m_subError = new ErrorStruct();
  wcscpy_s(m_stringVal[0], s1);
  subError->CopyTo(m_subError);
}

ErrorStruct::~ErrorStruct()
{
  if (m_subError != NULL)
    delete m_subError;
}
