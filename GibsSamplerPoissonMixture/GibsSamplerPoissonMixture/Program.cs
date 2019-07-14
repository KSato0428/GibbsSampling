using System;

namespace GibsSamplerPoissonMixture
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Hello World!");
            double[] a = new double[] { 0.3, 0.3, 0.4 };
            double[,] b = new double[,] { { 4.0 }, { 15.0 }, { 29.0 } };
            PoissonMixture pm = new PoissonMixture(a, b);
            int n = 10000;
            var xs = pm.Sampling(n);
            double[][] dxs = new double[n][];
            for(int i = 0; i < n; i++)
            {
                dxs[i] = new double[1];
                dxs[i][0] = xs[i][0];
            }
            /*
            for(int i = 0; i < n; i++)
            {
                Console.WriteLine(xs[i][0]);
            }*/

            int[] chrt = new int[100];
            for (int i = 0; i < 100; i++) chrt[i] = 0;
            for(int i = 0; i < n; i++)
            {
                if (xs[i][0] < 100) chrt[xs[i][0]]++;
            }
            for(int i = 0; i < 100; i++)
            {
                Console.Write(chrt[i].ToString() + "," );
            }

            PMGibbsSampling gs = new PMGibbsSampling(n, a.Length, b.GetLength(1), dxs);
            Console.WriteLine("Begin");
            gs.Sampling();
            Console.WriteLine("End");
            for(int i = 0; i < 200; i++)
            {
                Func<double[], double[,], double> f = (ap, bp) =>
                 {
                     double pi = 0.0;
                     for(int k = 0; k < ap.Length; k++)
                     {
                         double tmp = Math.Log(ap[k]) - bp[k, 0] + i * Math.Log(bp[k, 0]) - gs.LogFact(i);
                         
                         pi += Math.Exp(tmp);
                     }
                     return pi;
                 };
                double ef = gs.ExpectedF(f);
                Console.Write(string.Format("{0:F8}",ef)+ ",");
            }
            Console.Read();
        }
    }
}
