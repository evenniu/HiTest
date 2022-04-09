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
        const string DLLPath = @"Project1.dll";
        [DllImport(DLLPath, EntryPoint = "CreateSimple", CallingConvention = CallingConvention.StdCall, CharSet = CharSet.None)]
        public static extern void CreateSimple();
        [DllImport(DLLPath, EntryPoint = "Add", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern int Add(int x, int z, ref int sum);

        [DllImport(DLLPath, EntryPoint = "Divide", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern double Divide(double a, double b);

        [DllImport(DLLPath, EntryPoint = "Multiply", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern int Multiply(double a, double b);

        [DllImport(DLLPath, EntryPoint = "TestMatrix", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern void TestMatrix(int np);
        [DllImport(DLLPath, EntryPoint = "Calculate", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern void Calculate([MarshalAs(UnmanagedType.LPTStr)] string plaintext);

        [DllImport(DLLPath, EntryPoint = "Release", CallingConvention = CallingConvention.Cdecl, CharSet = CharSet.None)]
        public static extern void Release();//释放
    }
}
