using System;
using System.Collections.Generic;
using System.Text;
using MathNet.Numerics;
using MathNet.Numerics.Distributions;

namespace GibsSamplerPoissonMixture
{
    class PoissonMixture
    {
        double[] a;
        double[,] b;
        int H, M;
        Random rnd;
        public PoissonMixture(double[] a, double[,] b)
        {
            this.rnd = new Random();
            this.a = new double[a.Length];
            this.b = new double[b.GetLength(0), b.GetLength(1)];
            this.H = b.GetLength(0);
            this.M = b.GetLength(1);
            a.CopyTo(this.a, 0);
            for(int i = 0; i < H; i++)
            {
                for (int j = 0; j < M; j++) this.b[i, j] = b[i, j];
            }
            for(int i = 0; i < a.Length; i++)
            {

            }
        }

        public int[][] Sampling(int numb = 1000)
        {
            int[][] ans = new int[numb][];
            for (int i = 0; i < numb; i++) ans[i] = new int[M];
            for(int l = 0; l < numb; l++)
            {
                double r = rnd.NextDouble();
                int id = 0;
                double sumofa = 0.0;
                for(int i = 0; i < H; i++)
                {
                    sumofa += a[i];
                    if (r < sumofa)
                    {
                        id = i;
                        break;
                    }
                }
                for(int m = 0; m < M; m++)
                {
                    ans[l][m] = (new Poisson(b[id, m], rnd)).Sample();
                }
            }
            return ans;
        } 
    }
}
