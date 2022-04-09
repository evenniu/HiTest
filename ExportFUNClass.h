#pragma once
#include "windows.h"

#include "FUNClass.h"

using namespace MyApp;
FUNClass* g_pSimple = NULL;

// 创建对象
extern "C" __declspec(dllexport)  void CreateSimple();
extern "C" __declspec(dllexport)  int Add(int x, int y, int& sum);
extern "C" __declspec(dllexport)  int Divide(int a, int b);
extern "C" __declspec(dllexport)  void TestMatrix(int np);
extern "C" __declspec(dllexport)  double Multiply(double a, double b);
extern "C" __declspec(dllexport)  void SetName(LPCTSTR sName);
/* return a - b */
extern "C" __declspec(dllexport)  double Subtract(int a, int b);
extern "C" __declspec(dllexport) void enycode(LPCTSTR plaintext);
extern "C" __declspec(dllexport) void Calculate(LPCTSTR plaintext);
extern "C" __declspec(dllexport) void Release();



