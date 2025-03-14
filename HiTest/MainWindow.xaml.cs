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
                if(values.Length < 6)
                {
                    continue;
                }
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
        const int COLUM = 8;
        public MainWindow()
        {
            InitializeComponent();
            ImportSimleDLL.CreateSimple();
            string filePath = "./test/meas.xyz";
            List<double[]> measData = ReadDataFromFile(filePath);
           
            filePath = "./test/nom.xyz";
            List<double[]> nomData = ReadDataFromFile(filePath);
            IntPtr ptr_MeasPoints = IntPtr.Zero;
            IntPtr ptr_NomPoints = IntPtr.Zero;
            int rows, cols; 
            IntPtr p_measPoi = IntPtr.Zero;
            IntPtr p_nomPoi  = IntPtr.Zero;
            try
            {
                rows = measData.Count();
               
                p_measPoi= Convert2DArrayToPointer(measData);
                p_nomPoi = Convert2DArrayToPointer(nomData);
                cols = nomData[0].Length;
                bool result = ImportSimleDLL.LoadPoints(p_measPoi, rows,p_nomPoi, nomData.Count,cols);
                
                if (result)
                {
                    ImportSimleDLL.CalcBestFit(0, 0, 0, false, 0, 0);
                }
            }
            catch(Exception ex)
            {
                MessageBox.Show(ex.ToString());
            }
            finally
            {
                FreePointerToPointer(p_measPoi, measData.Count);
                FreePointerToPointer(p_nomPoi, nomData.Count);
                ImportSimleDLL.Release();
            }
        }

        private void Win_closed(object sender, EventArgs e)
        {
            ImportSimleDLL.Release();
        }
    }
}
