using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace HiTest
{
    /// <summary>
    /// MainWindow.xaml 的交互逻辑
    /// </summary>
    public partial class MainWindow : Window
    {
        static List<double[]> ReadDataFromFile(string filePath)
        {
            List<double[]> data = new List<double[]>();

            // 检查文件是否存在
            if (!File.Exists(filePath))
            {
                Console.WriteLine("File not found!");
                return data;
            }

            // 逐行读取文件内容
            foreach (string line in File.ReadLines(filePath))
            {
                // 分割每一行的数据
                string[] values = line.Split(new char[] { ' ', ',', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                double[] row = new double[values.Length];

                // 将字符串转换为 double 类型
                for (int i = 0; i < values.Length; i++)
                {
                    if (double.TryParse(values[i], out double value))
                    {
                        row[i] = value;
                    }
                    else
                    {
                        Console.WriteLine($"Invalid value: {values[i]}");
                    }
                }

                // 添加到数据列表
                data.Add(row);
            }

            return data;
        }


        // 将 List<double[]> 转换为指针
        public static IntPtr Convert2DArrayToPointer(List<double[]> array)
        {
            int rows = array.Count;
            int cols = array[0].Length;

            // 创建一个指针数组
            IntPtr[] pointers = new IntPtr[rows];
            for (int i = 0; i < rows; i++)
            {
                pointers[i] = Marshal.AllocHGlobal(cols * sizeof(double));
                Marshal.Copy(array[i], 0, pointers[i], cols);
            }

            // 创建一个指针数组的指针
            IntPtr pointerToPointer = Marshal.AllocHGlobal(rows * IntPtr.Size);
            Marshal.Copy(pointers, 0, pointerToPointer, rows);

            return pointerToPointer;
        }

        // 释放指针数组的内存
        public static void FreePointerToPointer(IntPtr pointerToPointer, int rows)
        {
            // 获取指针数组
            IntPtr[] pointers = new IntPtr[rows];
            Marshal.Copy(pointerToPointer, pointers, 0, rows);

            // 释放每一行的内存
            for (int i = 0; i < rows; i++)
            {
                Marshal.FreeHGlobal(pointers[i]);
            }

            // 释放指针数组的内存
            Marshal.FreeHGlobal(pointerToPointer);
        }
        int sum = 0;
        const int COLUM = 8;
        public MainWindow()
        {
            InitializeComponent();
            ImportSimleDLL.CreateSimple();
            string filePath = "meas.txt";
            List<double[]> Measdata = ReadDataFromFile(filePath);
           
            filePath = "nom.txt";
            List<double[]> Nomdata = ReadDataFromFile(filePath);
            IntPtr ptr_MeasPoints = IntPtr.Zero;
            IntPtr ptr_NomPoints = IntPtr.Zero;
            int rows, cols;
            try
            {
                // 将二维数组转换为一维数组
                rows = Measdata.Count();
                cols = Measdata[0].Length;
                double[] flatArray = new double[rows * cols];
                for (int i = 0; i < rows; i++)
                {
                    for (int j = 0; j < cols; j++)
                    {
                        flatArray[i * cols + j] = Measdata[i][j];
                    }
                }

                // 创建一个指针数组
                IntPtr[] pointers = new IntPtr[rows];
                for (int i = 0; i < rows; i++)
                {
                    pointers[i] = Marshal.AllocHGlobal(cols * sizeof(double));
                    Marshal.Copy(flatArray, i * cols, pointers[i], cols);
                }

                // 创建一个指针数组的指针
                IntPtr pointerToPointer = Marshal.AllocHGlobal(rows * IntPtr.Size);
                Marshal.Copy(pointers, 0, pointerToPointer, rows);
                bool result = ImportSimleDLL.LoadPoint(pointerToPointer, rows);
                // 释放内存
                for (int i = 0; i < rows; i++)
                {
                    Marshal.FreeHGlobal(pointers[i]);
                }
                Marshal.FreeHGlobal(pointerToPointer);
                if (result)
                {
                    ImportSimleDLL.CalcBestFit(0, 0, 0, false, 0, 0);
                }
            }
            catch(Exception ex)
            {

            }finally
            {
                ImportSimleDLL.Release();
            }

           


        }

        private void Win_closed(object sender, EventArgs e)
        {
            ImportSimleDLL.Release();
        }
    }
}
