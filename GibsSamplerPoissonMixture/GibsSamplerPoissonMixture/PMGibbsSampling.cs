using System;
using System.Collections.Generic;
using System.Text;

using MathNet.Numerics;
using MathNet.Numerics.Distributions;

namespace GibsSamplerPoissonMixture
{
    class PMGibbsSampling
    {
        double alpha = 1.0, beta = 1.0, gamma = 1.0;
        int burnIn = 1000;
        int Rgs = 10000;
        bool sampq;

        double[] a;
        double[,] b;
        int[] y;
        int n,H,M;

        public double[][] ra;
        public double[][,] rb;

        double[][] xs;
        //n: sample size, K: mixture components, M: dimension size, xs: data
        public PMGibbsSampling(int n, int H, int M, double[][] xs, int Rgs = 10000)
        {
            sampq = false;
            a = new double[H];
            b = new double[H, M];
            y = new int[n];
            this.Rgs = Rgs;
            this.n = n;
            this.H = H;
            this.M = M;
            ra = new double[this.Rgs][];
            rb = new double[this.Rgs][,];
            for(int i = 0; i < this.Rgs; i++)
            {
                ra[i] = new double[H];
                rb[i] = new double[H, M];
            }
            if (xs.GetLength(0) != n) throw new Exception();
            this.xs = new double[n][];
            for (int i = 0; i < n; i++)
            {
                this.xs[i] = new double[xs[i].Length];
                xs[i].CopyTo(this.xs[i], 0);
            }
            logFactBd = new double[1000];
            logFactBd[0] = logFactBd[1] = 0;
            for (int i = 2; i < 1000; i++)
            {
                logFactBd[i] = logFactBd[i - 1] + Math.Log(i);
            }
        }

        double[] logFactBd;
        public double LogFact(int x)
        {
            if (x < 1000)
            {
                return logFactBd[x];
            } else
            {
                double ans = 0;
                for (int i = 2; i < x+1; i++)
                {
                    ans = ans + Math.Log(i);
                }
                return ans;
            }
        }
        public void Sampling()
        {
            sampq = true;
            double[] d = new double[H];
            double[] dalpha = new double[H];
            double[,] dbeta = new double[H,M], dgamma = new double[H,M];
            Random rnd = new Random();

            //initialize a,b
            double A = 0;
            for(int k = 0; k < H; k++)
            {
                a[k] = rnd.NextDouble();
                A += a[k];
            }
            for (int k = 0; k < H; k++) a[k] /= A;

            for(int k = 0; k < H; k++)
            {
                for(int m = 0; m < M; m++)
                {
                    b[k,m] = 2.0 + rnd.NextDouble();
                }
            }

            for(int cnt = 0; cnt < burnIn + Rgs; cnt++)
            {
                // generate y^n
                for(int i = 0; i < n; i++)
                {
                    double D = 0;
                    for(int j = 0; j < H; j++)
                    {
                        double dij = Math.Log(a[j]);
                        for(int m = 0; m < M; m++)
                        {
                            dij += -b[j, m] + xs[i][m] * Math.Log(b[j, m]) - LogFact((int)xs[i][m]);
                        }
                        d[j] = Math.Exp(dij);
                        D += d[j];
                    }
                    d[0] /= D;
                    for (int j = 1; j < H; j++)
                    {
                        d[j] /= D;
                        d[j] += d[j - 1];
                    }
                    double r = rnd.NextDouble();
                    int yi = 0;
                    for (int j = 0; j < H; j++)
                    {
                        if(r < d[j])
                        {
                            yi = j;
                            break;
                        }
                    }
                    y[i] = yi;
                }

                //generate a
                for(int k = 0; k < H; k++)
                {
                    dalpha[k] = alpha;
                }
                for(int i = 0; i < n; i++)
                {
                    dalpha[y[i]] += 1.0;
                }
                a = (new Dirichlet(dalpha, rnd)).Sample();
                //generate b
                for (int k = 0; k < H; k++)
                {
                    for (int m = 0; m < M; m++)
                    {
                        dbeta[k, m] = beta;
                        dgamma[k, m] = 1 / gamma;
                    }
                }
                for (int i = 0; i < n; i++)
                {
                    for(int m = 0; m < M; m++)
                    {
                        dbeta[y[i], m] += xs[i][m];
                        dgamma[y[i], m] += 1.0;
                    }
                }
                for (int k = 0; k < H; k++)
                {
                    for(int m = 0; m < M; m++)
                    {
                        b[k, m] = (new Gamma(dbeta[k, m], dgamma[k, m], rnd)).Sample();
                    }
                }
                //memorize
                if(cnt >= burnIn)
                {
                    for(int k = 0; k < H; k++)
                    {
                        ra[cnt - burnIn][k] = a[k];
                    }
                    for(int k = 0; k < H; k++)
                    {
                        for(int m = 0; m < M; m++)
                        {
                            rb[cnt - burnIn][k, m] = b[k, m];
                        }
                    }
                }
            }
        }
        public double ExpectedF(Func<double[],double[,], double> f)
        {
            if (!sampq) throw new Exception();
            double sm = 0.0;
            for(int i = 0; i < Rgs; i++)
            {
                sm += f(ra[i], rb[i]);
            }
            return sm/((double)Rgs);
        }
        public void Print()
        {
            for(int i = 0; i < Rgs; i++)
            {
                Console.WriteLine("a");
                for(int j = 0; j < H; j++)
                {
                    Console.Write(ra[i][j].ToString() + ",");
                }
                Console.WriteLine("b");
                for(int j = 0; j < H; j++)
                {
                    for(int k = 0; k < M; k++)
                    {
                        Console.Write(rb[i][j, k].ToString() + ",");
                    }
                }
            }
        }
    }
}
