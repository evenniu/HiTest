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
            int iResult = ImportSimleDLL.Add(10, 3, ref sum);
            double d2 = ImportSimleDLL.Divide(10.0, 3.0);
            double m1 = ImportSimleDLL.Multiply(10, 10);
            string plainText = "12d33qwe";
            int numpts = 16;
            string filePath = "meas.txt";
            List<double[]> Measdata = ReadDataFromFile(filePath);
           
            filePath = "nom.txt";
            List<double[]> Nomdata = ReadDataFromFile(filePath);
            IntPtr ptr_MeasPoints = IntPtr.Zero;
            IntPtr ptr_NomPoints = IntPtr.Zero;
            int rows, cols;
            try
            {
     
                double[,] points = new double[63, 6]
                {
                    {   30.988247   ,   6.020396    ,   18.014885   ,   0.897296    ,   -0.44143    ,   0   },
                    {   31.097179   ,   6.239536    ,   18.014845   ,   0.876845    ,   -0.480773   ,   0   },
                    {   31.245785   ,   6.484317    ,   18.013525   ,   0.832046    ,   -0.554706   ,   0   },
                    {   31.420549   ,   6.723003    ,   18.01368    ,   0.775016    ,   -0.631941   ,   0   },
                    {   31.627651   ,   6.950818    ,   18.013124   ,   0.702666    ,   -0.711519   ,   0   },
                    {   31.84613    ,   7.145291    ,   18.01388    ,   0.63355 ,   -0.773702   ,   0   },
                    {   32.067757   ,   7.312995    ,   18.013113   ,   0.571805    ,   -0.820389   ,   0   },
                    {   32.29871    ,   7.460814    ,   18.012386   ,   0.508472    ,   -0.861078   ,   0   },
                    {   32.562943   ,   7.602392    ,   18.009466   ,   0.417274    ,   -0.908781   ,   0   },
                    {   32.817093   ,   7.702025    ,   18.010332   ,   0.35126 ,   -0.936278   ,   0   },
                    {   33.046227   ,   7.785419    ,   18.010115   ,   0.310925    ,   -0.950434   ,   0   },
                    {   33.258369   ,   7.846795    ,   18.007092   ,   0.288407    ,   -0.957508   ,   0   },
                    {   33.421432   ,   7.898388    ,   18.009266   ,   0.313653    ,   -0.949538   ,   0   },
                    {   33.607906   ,   7.963253    ,   18.007101   ,   0.318796    ,   -0.947823   ,   0   },
                    {   33.80756    ,   8.027608    ,   18.00544    ,   0.306375    ,   -0.951911   ,   0   },
                    {   33.995724   ,   8.089098    ,   18.005358   ,   0.313629    ,   -0.949546   ,   0   },
                    {   34.184566   ,   8.151443    ,   18.004932   ,   0.313729    ,   -0.949512   ,   0   },
                    {   34.373447   ,   8.213798    ,   18.003645   ,   0.312852    ,   -0.949802   ,   0   },
                    {   34.564529   ,   8.277116    ,   18.000866   ,   0.314323    ,   -0.949316   ,   0   },
                    {   34.757824   ,   8.340645    ,   17.999826   ,   0.310008    ,   -0.950734   ,   0   },
                    {   34.950672   ,   8.40298 ,   17.998703   ,   0.305997    ,   -0.952033   ,   0   },
                    {   35.139469   ,   8.463786    ,   17.996902   ,   0.310404    ,   -0.950605   ,   0   },
                    {   35.327763   ,   8.525486    ,   17.996264   ,   0.309406    ,   -0.95093    ,   0   },
                    {   35.519703   ,   8.588112    ,   17.995132   ,   0.31035 ,   -0.950622   ,   0   },
                    {   35.712097   ,   8.650501    ,   17.993967   ,   0.306863    ,   -0.951754   ,   0   },
                    {   35.888557   ,   8.707065    ,   17.9939 ,   0.320675    ,   -0.947189   ,   0   },
                    {   36.072506   ,   8.772489    ,   17.991741   ,   0.324   ,   -0.946057   ,   0   },
                    {   36.256699   ,   8.83327 ,   17.990072   ,   0.329943    ,   -0.944001   ,   0   },
                    {   36.395153   ,   8.88404 ,   17.988939   ,   0.382451    ,   -0.923976   ,   0   },
                    {   36.511822   ,   8.937014    ,   17.989382   ,   0.444927    ,   -0.895567   ,   0   },
                    {   36.613743   ,   8.991911    ,   17.988281   ,   0.520438    ,   -0.853899   ,   0   },
                    {   36.731232   ,   9.072944    ,   17.986122   ,   0.569428    ,   -0.822041   ,   0   },
                    {   36.849991   ,   9.155322    ,   17.985647   ,   0.610683    ,   -0.791875   ,   0   },
                    {   36.940285   ,   9.232775    ,   17.984505   ,   0.675604    ,   -0.737264   ,   0   },
                    {   37.047001   ,   9.338263    ,   17.982931   ,   0.711119    ,   -0.703071   ,   0   },
                    {   37.141964   ,   9.435429    ,   17.982943   ,   0.752553    ,   -0.658531   ,   0   },
                    {   37.213768   ,   9.5273  ,   17.983379   ,   0.805894    ,   -0.592059   ,   0   },
                    {   37.280952   ,   9.624676    ,   17.982901   ,   0.851405    ,   -0.52451    ,   0   },
                    {   37.332489   ,   9.718536    ,   17.983183   ,   0.897941    ,   -0.440117   ,   0   },
                    {   37.385979   ,   9.842963    ,   17.982449   ,   0.923129    ,   -0.384491   ,   0   },
                    {   37.453552   ,   10.015656   ,   17.978603   ,   0.929731    ,   -0.368241   ,   0   },
                    {   37.519722   ,   10.179262   ,   17.978167   ,   0.937696    ,   -0.347457   ,   0   },
                    {   37.55827    ,   10.289532   ,   17.978212   ,   0.963436    ,   -0.267937   ,   0   },
                    {   37.579906   ,   10.390918   ,   17.978256   ,   0.984552    ,   -0.175093   ,   0   },
                    {   37.596344   ,   10.515795   ,   17.976778   ,   0.995081    ,   -0.099065   ,   0   },
                    {   37.605183   ,   10.647989   ,   17.977451   ,   0.999503    ,   -0.031526   ,   0   },
                    {   37.604969   ,   10.778255   ,   17.977745   ,   0.999372    ,   0.035429    ,   0   },
                    {   37.595959   ,   10.90496    ,   17.978498   ,   0.993843    ,   0.110795    ,   0   },
                    {   37.575306   ,   11.039435   ,   17.978853   ,   0.984515    ,   0.175303    ,   0   },
                    {   37.545849   ,   11.184471   ,   17.97975    ,   0.97491 ,   0.222599    ,   0   },
                    {   37.513771   ,   11.311358   ,   17.981236   ,   0.95608 ,   0.293107    ,   0   },
                    {   37.468044   ,   11.439245   ,   17.982605   ,   0.93677 ,   0.349946    ,   0   },
                    {   37.418243   ,   11.565848   ,   17.983744   ,   0.911219    ,   0.411922    ,   0   },
                    {   37.350422   ,   11.698461   ,   17.983973   ,   0.889491    ,   0.456953    ,   0   },
                    {   37.263977   ,   11.86365    ,   17.983992   ,   0.88389 ,   0.467694    ,   0   },
                    {   37.171043   ,   12.037066   ,   17.983665   ,   0.882515    ,   0.470284    ,   0   },
                    {   37.076607   ,   12.214514   ,   17.983614   ,   0.881489    ,   0.472205    ,   0   },
                    {   36.984669   ,   12.386287   ,   17.983603   ,   0.880392    ,   0.474248    ,   0   },
                    {   36.885208   ,   12.56918    ,   17.982431   ,   0.884068    ,   0.467358    ,   0   },
                    {   36.793602   ,   12.747507   ,   17.982218   ,   0.883765    ,   0.467931    ,   0   },
                    {   36.704472   ,   12.910538   ,   17.981791   ,   0.875925    ,   0.482447    ,   0   },
                    {   36.600647   ,   13.098824   ,   17.98418    ,   0.884106    ,   0.467287    ,   0   },
                    {   36.520866   ,   13.257736   ,   17.985554   ,   0.873506    ,   0.486813    ,   0   }

                };

                // 将二维数组转换为一维数组
                rows = points.GetLength(0);
                cols = points.GetLength(1);
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
                if (result)
                {
                    ImportSimleDLL.CalcBestFit(0, 0, 0, false, 0, 0);
                }
            }
            catch(Exception ex)
            {

            }

            ImportSimleDLL.Release();


        }

        private void Win_closed(object sender, EventArgs e)
        {
            ImportSimleDLL.Release();
        }
    }
}
