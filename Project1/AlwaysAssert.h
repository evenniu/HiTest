#pragma once

template <class ValueType>
void alwaysAssertCore(const ValueType& value, const std::string& valueString, int lineNumber,
                      const std::string& fileName)
{
  if(!value)
  {
    throw std::logic_error("always-assert failed: '" + valueString + "' at line " + std::to_string(lineNumber) +
                           " of file " + fileName);
  }
}

#define alwaysAssert(value) alwaysAssertCore((value), #value, __LINE__, __FILE__)
