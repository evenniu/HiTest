using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

namespace HiTest
{
    public static class ImportSimleDLL
    {
        const string DLLPath = @"BladeMath.dll";
        [DllImport(DLLPath, EntryPoint = "CreateSimple", CallingConvention = CallingConvention.StdCall, CharSet = CharSet.None)]
        public static extern void CreateSimple();

        [DllImport(DLLPath, EntryPoint = "Add", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern int Add(int x, int z, ref int sum);

        [DllImport(DLLPath, EntryPoint = "Divide", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern double Divide(double a, double b);

        [DllImport(DLLPath, EntryPoint = "Multiply", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern double Multiply(double a, double b);

        [DllImport(DLLPath, EntryPoint = "TestMatrix", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern void TestMatrix(int np);

        [DllImport(DLLPath, EntryPoint = "Calculate", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern void Calculate([MarshalAs(UnmanagedType.LPTStr)] string plaintext);

        [DllImport(DLLPath, EntryPoint = "Release", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern void Release();//释放

        [DllImport(DLLPath, EntryPoint = "GetFunc", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern void GetFunc(int[] resultArray, int arraySize);        
        
        [DllImport(DLLPath, EntryPoint = "GetSecNum", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern void GetSecNum(int[] resultArray, ref int size);        
        [DllImport(DLLPath, EntryPoint = "LoadPoint", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern bool LoadPoint(IntPtr measPoints,  int numpoints);

        [DllImport(DLLPath, EntryPoint = "LoadPoints", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern bool LoadPoints(IntPtr measPoints, int numpoints, IntPtr nomPoints, int numNompoints, int cols);

        [DllImport(DLLPath, EntryPoint = "CalcBestFit", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern bool CalcBestFit(int BestFitType, int fitToMiddleOfZone, int Transfit, bool noRotate, int rotfit, int useNominal);


        [DllImport(DLLPath, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern bool GetDev(out IntPtr outdev, ref double maxdev, ref double mindev, out int size);

        [DllImport(DLLPath, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern void FreeMemory(IntPtr outdevPtr);



        [DllImport(DLLPath, CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern void GetFitResult(int report, ref double x, ref double y, ref double ang);
    }
}
