#pragma once

#ifdef BLADEMATH_EXPORTS
#define DLLEXPORT __declspec(dllexport)
#else // !defined (BLADEMATH_EXPORTS)
#define DLLEXPORT __declspec(dllimport)
#endif

enum BladeErrorType {
  BE_OK = 0,                     // No error
  BE_HARDCODEDTEXT,              // should fix these, there are some hard coded error text
  BE_BADCALCSFLAVOR,             // "Invalid CALCS line in Flavor file"
  BE_TOOMANYFLAVORCALCS,         // "Too many CALCS in Flavor file"
  BE_INCOMPLETECALCSBLOCK,       // "Incomplete CALCS block in Flavor file"
  BE_CALCSLINE,                  // "Invalid line in CALCS block in Flavor file (%s)"
  BE_BADBESTFITTYPE,             // "Invalid bestfit type in Flavor file"
  BE_BAD2POINTPERCENTS,          // "Invalid 2Point percentages in Flavor file"
  BE_BAD2POINTMODIFIER,          // "Invalid 2Point modifier in Flavor file"
  BE_BADBESTFIT,                 // "Unknown bestfit type in Flavor file"
  BE_BADFLAVORCOMMAND,           // "Unknown command in Flavor file"
  BE_MISSINGCALCMODIFIER,        // "Missing calculation modifier (%s)"
  BE_MISSINGMETHOD,              //  "Missing METHOD"
  BE_BADMETHOD,                  // "Unknown METHOD %s"
  BE_BADMULTIPLIER,              // "Invalid multiplier (%s)"
  BE_BADOFFSET,                  // "Invalid offset (%s)"
  BE_BADUSEBESTFIT,              // "Invalid bestfit index (%s)"
  BE_BADCALCMODIFIER,            // "Invalid calculation modifier (%s)"
  BE_BADFLAVORUSEFIT,            // "BestFit %d is Referenced, but not Defined in Flavor File"
  BE_NOFILECREATE,               // "Could Not Create File"
  BE_BAD_VECTOR_SECTION,         // "Warning:  Section %s contains suspicious vectors, please review!"
  BE_NOOPENNOM,                  // "Could not open converted NOM file (%s)???"
  BE_BADNOMUNITS,                // "UNITS Not Defined in NOM file"
  BE_BADNUMSECT,                 // "Invalid NUM_SECT line"
  BE_MISSINGSECTIONNOM,          // "Missing SECTION (number %d) in NOM file"
  BE_MISSINGORBADSECTION,        // "Missing or Bad SECTION (number %d) in NOM file"
  BE_INCOMPLETESECTION,          // "incomplete SECTION (number %d) in NOM file"
  BE_BADOFFSETSINNOM,            // "Bad CHANGES (number %d) in NOM file"
  BE_MISSINGPOINTNOM,            // "Missing point %d for Section %d in NOM file"
  BE_BADPOINTINNOM,              // "Bad point %d for Section %d in NOM file"
  BE_NOCALCSECTIONNOM,           // "Could not calculate Section %s in NOM file"
  BE_BADNOSEINDICES,             // "Invalid Indices in NOSE Definition (%s)"
  BE_BADTAILINDICES,             // "Invalid Indices in TAIL Definition (%s)"
  BE_FILENOTFOUND,               // "File (%s) Not Found"
  BE_DIFFERENTPART,              // "Part Name is Different (%s vs %s)!"
  BE_NOCREATETOLFILE,            // "Could Not Create Tolerance File (%s)"
  BE_PARSINGERROR,               // "File Parsing Error"
  BE_NOTITLE,                    // "No TITLE Line"
  BE_TOOSHORTGETINT,             // "GetInt [%d,%d] String Too Short"
  BE_BADGETINT,                  // "GetInt [%d,%d] [%s] Is Not An Int"
  BE_TOOSHORTGETDOUBLE,          // "GetDouble [%d,%d] String Too Short"
  BE_BADGETDOUBLE,               // "GetDouble [%d,%d] [%s] Not A Double"
  BE_TOOSHORTGETSTRING,          // "GetString [%d,%d] String Too Short"
  BE_NOUNITS,                    // "No ENTLISH/METRIC INDICATOR Line"
  BE_NOTPLANAR,                  // "PROFILE SECTION TYPE Must Be PLANAR"
  BE_CANNOTCOMPLETEBOTH,         // "Can not Complete Both Ends"
  BE_BADSECTION,                 // "Unknown SECTION NAME (%s)"
  BE_MISSINGLENOM,               // "Missing LE Nominal for section %s"
  BE_MISSINGTENOM,               // "Missing TE Nominal for section %s"
  BE_NOTENOUGHPOINTS,            // "Not enough points in scans for Section %s"
  BE_NOTRIMINTERSECTION,         // "Could Not Find Trim Intersection Point"
  BE_NOFLAVOR,                   // "No FLAVOR Record In File (%s)"
  BE_NOPART,                     // "No PART Record In File (%s)"
  BE_NODATE,                     // "No DATE Record In File (%s)"
  BE_BADDATE,                    // "Bad DATE (%s)"
  BE_NONOMINAL,                  // "No NOMINAL Record In File (%s)"
  BE_NOTOLERANCE,                // "No TOLERANCE Record In File (%s)"
  BE_DIFFNOMTOL,                 // "Different number of sections in Nom (%d) and Tol (%d) files"
  BE_DIFFMATHTOL,                // "Different number of sections in Math (%d) and Tol (%d) files"
  BE_NORADIUS,                   // "No RADIUS Record In File (%s)"
  BE_BADRADIUS,                  // "Bad RADIUS Record (%s)"
  BE_BADTEXTLINE,                // "Bad TEXT line (%s)"
  BE_MISSINGCALC,                // "Missing Calculation"
  BE_MISSINGPLATFORM,            // "Missing Platform Point"
  BE_BADPLATFORM,                // "Bad Platform Point"
  BE_MISSINGSECTIONHEADING,      // "Missing Section Heading"
  BE_BADSECTIONHEADING,          // "Bad Section Heading"
  BE_MISSINGACTUAL,              // "Missing Actual Point"
  BE_BADMEASUREDPOINT,           // "Invalid Point in RPT file, Section %s, Point %d of %d"
  BE_ZMISMATCH,                  // "Z values do not seem to agree between NOM and RPT files for section %s"
  BE_POINTORDERINGFAILED,        // "Point Ordering Failed"
  BE_POINTORDERINGFAILEDONKNOWNSECTION,  // "Point Ordering Failed on Section %s"
  BE_BADKEYWORD,                 // "Invalid keyword (%s)"
  BE_COMPENSATIONFAILED,         // "Probe Compensation Failed for Section %s!!!"
  BE_MEANCAMBERFAILED,           // "Mean Camber Line Failed for %s"
  BE_NOOPENTOL,                  // "Could not open converted Tol file (%s)???"
  BE_FLAVORFILEERROR,            // "Flavor File Error: %s"
  BE_TOLFILEERROR,               // "Tol File Error: %s"
  BE_NOMERROR,                   // "Nom File Error: %s"
  BE_BAD,                        // "Bad %s (%s)"  <-  actually the string resource is just "Bad", will need to build up the string for this error
  BE_TOOMANYNUMSECT,             // "Too many NUM_SECT lines"
  BE_NUMSECTBEFORESECTION,       // "SECTION lines must follow NUM_SECT line"
  BE_TOOMANYSECTIONS,            // "Too many SECTION lines"
  BE_MISSINGSECTION,             // "Missing section name"
  BE_NOLOTPLOTINTE,              // "TE File does not support Lot Plot"
  BE_ASPROBLEMS,                 // "AS File Creation Failed\nProblems encountered building AS file.\n%d Sections Processed, %d Sections skipped."
  BE_NODSTFILE,                  // "Distribution File (%s) not found."
  BE_BADDSTFILE,                 // "Bad Distribution File Format"
  BE_NEEDDPDATABASE,             // "DataPage Database must be specified"
  BE_NOTENOUGHTOLERANCED,        // "Section %s has no points toleranced, but is using nominal point fit."
  BE_TIPARCDOESNOTINTERSECTMCL,  // "The tip-arc does not intersect the mean camber line in section %s; try increasing L/TE_OFFSET."
  BE_NOMINALMEANCAMBERNOTALLOWED,// "NOMINALMEANCAMBER is no longer a valid option; we recommend using the new mean camber line specified with the CAMBER_OFFSETS parameter."
  
  BE_UNKNOWNERROR
};

enum BladeErrorArgTypes {
  ArgsNoError = -1, // No Error
  ArgsNone = 0,     // No args
  ArgsS,            // String
  ArgsD,            // Int
  ArgsDD,           // Int, Int
  ArgsSS,           // String, String
  ArgsDDS,          // Int, Int, String
  ArgsSDD,          // String, int, int
  ArgsSubS          // SubError, String
};


class DLLEXPORT ErrorStruct
{
public:
  void Init();
  ErrorStruct();
  ErrorStruct(BladeErrorType errnum);
  ErrorStruct(BladeErrorType errnum, const wchar_t *s1);
  ErrorStruct(BladeErrorType errnum, int d1);
  ErrorStruct(BladeErrorType errnum, int d1, int d2);
  ErrorStruct(BladeErrorType errnum, const wchar_t *s1, const wchar_t *s2);
  ErrorStruct(BladeErrorType errnum, int d1, int d2, const wchar_t *s3);
  ErrorStruct(BladeErrorType errnum, const wchar_t *s1, int d2, int d3);
  ErrorStruct(BladeErrorType errnum, ErrorStruct *subError, const wchar_t *s1); 
  ~ErrorStruct();

  BladeErrorType m_errorCode;   // what the error is
  ErrorStruct *m_subError;      // to call out an error in a file and perserve the lower level error
  int m_lineNumber;             // if reading a file, what line did it occur on
  int m_intVal[3];              // int values
  double m_doubleVal[3];        // double values (not used yet)
  wchar_t m_stringVal[3][MAXBUFSZ];  // string values
  BladeErrorArgTypes m_argType; // defines the arguments so the error reporter will know how to put together the message
  void CopyTo(ErrorStruct *err) const;
};

// this class will hold an array of errors.  The calling program should create one of these and pass it into the methods that can return errors

class CBladeError
{
  std::vector<ErrorStruct> m_errors;

public:
//   DLLEXPORT int m_lastError;
//   ErrorStruct m_errors[5];  // circular buffer, last 5 errors
  DLLEXPORT void ClearErrors();
  DLLEXPORT void AddError(ErrorStruct *err);
  DLLEXPORT bool LastError(ErrorStruct *err);
  DLLEXPORT std::size_t numberOfErrors() const;
  DLLEXPORT bool getError(std::size_t index, ErrorStruct *err) const;
  DLLEXPORT bool isOK();
  DLLEXPORT CBladeError();
  DLLEXPORT ~CBladeError();
};

