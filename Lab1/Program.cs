using System;

namespace Lab1
{
    class Program
    {
        static int variant = 23;
        static private int indexinArray = 0;

        public static int GetIndexinArray()
        {
            return indexinArray;
        }

        public static void SetIndexinArray(int value)
        {
            indexinArray = value;
        }


        static void Main(){
            Random rand = new Random();
            float M = getM(index: GetIndexinArray());
            float Rkr = getRkr(index: GetIndexinArray());
                int y_max = (30 - variant) * 10;
                int y_min = (20 - variant) * 10;
                int x1_min = -30;
                int x1_max = 0;
                int x2_min = -15;
                int x2_max = 35;
                int[,] xn = {
                {-1,-1},
                {-1, 1},
                { 1,-1}
                };

                float[][] y = new float[3][];

                for (int i = 0; i < 3; i++)
                {
                y[i] = new float[(int)M];
                }

                for (int i = 0; i < y.Length; i++)
                {
                    for (int j = 0; j < y[i].Length; j++)
                    {
                        y[i][j] = rand.Next(y_min, y_max);

                    }
                }

                float[] avY = aver(y);
                float sigmateta0 =(2 * (2 * M - 2) / (M * (M - 4)));
                double sigmaTeta = Math.Sqrt(sigmateta0);

                float[] Fuv = new float[3];
                double[] teta = new double[3];
                double[] Ruv = new double[3];

                Fuv[0] = (fuv(disp(y)[0], disp(y)[1]));
                Fuv[1] = (fuv(disp(y)[2], disp(y)[0]));
                Fuv[2] = (fuv(disp(y)[2], disp(y)[1]));

                double m1 = (double)(M - 2) / M;
                teta[0] = m1 * Fuv[0];
                teta[1] = m1 * Fuv[1];
                teta[2] = m1 * Fuv[2];

                Ruv[0] = Math.Abs((teta[0] - 1) / sigmaTeta);
                Ruv[1] = Math.Abs((teta[1] - 1) / sigmaTeta);
                Ruv[2] = Math.Abs((teta[2] - 1) / sigmaTeta);

                int b = 0;
                for (int a = 0; a < Ruv.Length; a++)
                {
                    if (Ruv[a] > Rkr)
                    {
                        //Main();
                        //Console.WriteLine("Помилка, повторіть експеримент");
                    }
                    else
                    {
                        b = b + 1;
                    }
                }


            if (b == 3)
            {
                if (indexinArray==0) {
                    Console.WriteLine("Y = b0 + b1*X1 + b2*X2");
                }
                float mx1 = (float)(xn[0, 0] + xn[1, 0] + xn[2, 0]) / 3;
                float mx2 = (float)(xn[0, 1] + xn[1, 1] + xn[2, 1]) / 3;
                float my = (avY[0] + avY[1] + avY[2]) / 3;

                float a1 = xn[0, 0] * xn[0, 0] + xn[1, 0] * xn[1, 0] + xn[2, 0] * xn[2, 0];
                float a2 = xn[0, 0] * xn[0, 1] + xn[1, 0] * xn[1, 1] + xn[2, 0] * xn[2, 1];
                float a3 = xn[0, 1] * xn[0, 1] + xn[1, 1] * xn[1, 1] + xn[2, 1] * xn[2, 1];

                float a11 = (xn[0, 0] * avY[0] + xn[1, 0] * avY[1] + xn[2, 0] * avY[2]) / 3;
                float a22 = (xn[0, 1] * avY[0] + xn[1, 1] * avY[1] + xn[2, 1] * avY[2]) / 3;

                float b0 = discr(my, mx1, mx2, a11, a1, a2, a22, a2, a3) / discr(1, mx1, mx2, mx1, a1, a2, mx2, a2, a3);
                float b1 = discr(1, my, mx2, mx1, a11, a2, mx2, a22, a3) / discr(1, mx1, mx2, mx1, a1, a2, mx2, a2, a3);
                float b2 = discr(1, mx1, my, mx1, a1, a11, mx2, a2, a22) / discr(1, mx1, mx2, mx1, a1, a2, mx2, a2, a3);

                float y_pr1 = b0 + b1 * xn[0, 0] + b2 * xn[0, 1];
                float y_pr2 = b0 + b1 * xn[1, 0] + b2 * xn[1, 1];
                float y_pr3 = b0 + b1 * xn[2, 0] + b2 * xn[2, 1];

                float dx1 = Math.Abs(x1_max - x1_min) / 2;
                float dx2 = Math.Abs(x2_max - x2_min) / 2;
                float x10 = (float)(x1_max + x1_min) / 2;
                float x20 = (float)(x2_max + x2_min) / 2;

                float koef0 = b0 - (b1 * x10 / dx1) - (b2 * x20 / dx2);
                float koef1 = b1 / dx1;
                float koef2 = b2 / dx2;

                float yP1 = koef0 + koef1 * x1_min + koef2 * x2_min;
                float yP2 = koef0 + koef1 * x1_max + koef2 * x2_min;
                float yP3 = koef0 + koef1 * x1_min + koef2 * x2_max;

                Console.WriteLine("Матриця планування для m =" + M);
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < M; j++)
                    {
                        Console.Write(y[i][j] + " ");
                    }
                    Console.WriteLine();
                }
                Console.WriteLine("Експериментальні значення критерію Романовського:");
                for (int i = 0; i < 3; i++)
                {
                    Console.WriteLine(Ruv[i]);
                }
                Console.WriteLine("Натуралізовані коефіцієнти: \na0 =" + koef0 + " a1 =" + koef1 + " a2 =" + koef2);
                Console.WriteLine("У практичний " + y_pr1 + " " + y_pr2 + " " + y_pr3 + "\nУ середній " + avY[0] + " " + avY[1] + " " + avY[2]);
                Console.WriteLine("У практичний норм." + Math.Round(yP1, 4) + " " + Math.Round(yP3, 4) + " " + Math.Round(yP2, 4));
                Console.WriteLine();
                SetIndexinArray(indexinArray + 1);
                if (indexinArray < 5) {
                    Main();
                }
            } else {
                    Main();
            }
        }
        




        public static float[] aver(float[][] y)
        {
            float[] avY = new float[y.Length];
            for (int i = 0; i < y.Length; i++)
            {
                float s = 0;
                for (int j = 0; j < y[i].Length; j++)
                {
                    s += y[i][j];
                }
                avY[i] = (float)s / y[i].Length;
            }
            return avY;
        }

        public static int getM(int index)
        {
            int[] M = { 6, 10, 12, 15, 20};
            int m = M[index];
            return m;
        }
        public static float getRkr(int index)
        {
            float[] Rkr = { 2, 2.29f, 2.39f, 2.49f, 2.62f };
            float rkr = Rkr[index];
            return rkr;
        }


        public static float[] disp(float[][] y)
        {
            float[] disp = new float[y.Length];
            for (int i = 0; i < y.Length; i++)
            {
                float s = 0;
                for (int j = 0; j < y[i].Length; j++)
                {
                    s += (j - aver(y)[i]) * (j - aver(y)[i]);
                }
                disp[i] = (float)(s / y[i].Length);
            }
            return disp;
        }




        public static float fuv(float u, float v)
        {
            if (u >= v)
            {
                return u / v;
            }
            else
            {
                return v / u;
            }
        }




        public static float discr(float x11, float x12, float x13, float x21, float x22, float x23, float x31, float x32, float x33)
        {
            return x11 * x22 * x33 + x12 * x23 * x31
                + x32 * x21 * x13 - x13 * x22 * x31
                - x32 * x23 * x11 - x12 * x21 * x33;
        }
    }
}