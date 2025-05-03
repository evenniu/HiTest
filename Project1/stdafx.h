// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once
using namespace std;

#pragma warning(disable:4820)
#pragma warning(disable:4365)
#pragma warning(disable:4710)
#pragma warning(disable:4127)
#pragma warning(disable:4668)
#pragma warning(disable:4350)
#pragma warning(disable:4266)
#pragma warning(disable:4714) // we can ignore warning C4714: function '...' marked as __forceinline not inlined


#include "targetver.h"

#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
// Windows Header Files:
#include <windows.h>
#include <math.h>
#include <cstdlib>
#include <string>
#include <time.h>
using namespace std;
#include <vector>
#include <utility>
#include "StringUtil.h"
#include "MATRIX.H"
#include "MATHUTIL.H"
#include "CURVE.H"
#include "CONSTANTS.H"
#include "BESTFIT.H"
#include "SECTION.H"
#include "BladeError.h"
#include "CBLADE.H"

#include "AlwaysAssert.h"

#define M_PI 3.14159265358979323846
extern void bugout(int level, const WCHAR* fmt, ...);
