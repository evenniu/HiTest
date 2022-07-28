using System;
using System.Collections.Generic;
using System.Linq;
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
            iResult = ImportSimleDLL.Multiply(10, 10);
            string plainText = "12d33qwe";
            
            ImportSimleDLL.TestMatrix(10);
            ImportSimleDLL.Calculate(plainText);
            

        }

        private void Win_closed(object sender, EventArgs e)
        {
            ImportSimleDLL.Release();
        }
    }
}
