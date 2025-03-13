using System;
using System.Collections.Generic;
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
        int sum = 0;

        public MainWindow()
        {
            InitializeComponent();
            ImportSimleDLL.CreateSimple();
            int iResult = ImportSimleDLL.Add(10, 3, ref sum);
            double d2 = ImportSimleDLL.Divide(10.0, 3.0);
            double m1 = ImportSimleDLL.Multiply(10, 10);
            string plainText = "12d33qwe";
            int numpts = 16;
            double[,] points = new double[16, 6]
            {
                {51.333500, -79.210700, 111.010000, -0.725300, -0.661900, 0.189200},
                {51.267300, -79.135800, 111.010000, -0.740500, -0.640700, 0.202200},
                {51.203100, -79.059100, 111.010000, -0.754900, -0.619600, 0.214200},
                {51.140700, -78.981000, 111.010000, -0.768400, -0.598800, 0.225100},
                {51.079800, -78.901600, 111.010000, -0.781200, -0.578300, 0.235300},
                {51.021700, -78.820300, 111.010000, -0.792700, -0.558400, 0.244000},
                {50.965100, -78.737900, 111.010000, -0.803600, -0.539000, 0.251800},
                {50.910400, -78.654100, 111.010000, -0.813700, -0.520300, 0.258700},
                {50.857500, -78.569300, 111.010000, -0.823000, -0.502400, 0.264600},
                {50.806000, -78.483600, 111.010000, -0.831800, -0.485000, 0.269700},
                {50.756800, -78.396500, 111.010000, -0.839700, -0.468800, 0.273700},
                {50.708300, -78.309100, 111.010000, -0.847300, -0.453000, 0.277200},
                {50.663700, -78.219600, 111.010000, -0.853400, -0.439700, 0.278700},
                {50.619400, -78.129900, 111.010000, -0.859500, -0.426700, 0.280100},
                {50.575200, -78.040300, 111.010000, -0.865500, -0.413700, 0.281500},
                {50.531000, -77.950500, 111.010000, -0.871400, -0.400800, 0.282800}
            };

            // 获取二维数组的指针

            // 将二维数组转换为一维数组
            int rows = points.GetLength(0);
            int cols = points.GetLength(1);
            double[] flatArray = new double[rows * cols];
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    flatArray[i * cols + j] = points[i, j];
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
            
            ImportSimleDLL.Release();


        }

        private void Win_closed(object sender, EventArgs e)
        {
            ImportSimleDLL.Release();
        }
    }
}
