using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace rekuperatorProekt.Models
{
    public class Model
    {//Обозначение переменных
        public double B { get; set; }
        public double ALFA { get; set; }
        public double LALFA { get; set; }
        public double VALFA { get; set; }
        public double VL { get; set; }
        public double TRO { get; set; }
        public double BETTAT { get; set; }
        public double CV1 { get; set; }
        public double TV1 { get; set; }
        public double CV2 { get; set; }
        public double TV2 { get; set; }
        public double ETA { get; set; }
        public double TD1 { get; set; }
        public double TD2 { get; set; }
        public double KSIDELTAT { get; set; }
        public double Z1 { get; set; }
        public double Z2 { get; set; }
        public double Z3 { get; set; }
        public double ETAD { get; set; }
        public double ETAV { get; set; }
        public double DDDK { get; set; }
        public double ALFADL1 { get; set; }
        public double KT { get; set; }
        public double DDVK { get; set; }
        public double S { get; set; }
        public double FIBL { get; set; }
        //Обозначение формул
        public double VB1 { get; set; }
        public double VBPR { get; set; }
        public double VD1 { get; set; }
        public double ID1 { get; set; }
        public double IV1 { get; set; }
        public double IV2 { get; set; }
        public double ID2 { get; set; }
        public double DELTAT1 { get; set; }
        public double P { get; set; }
        public double R { get; set; }
        public double DELTAT { get; set; }
        public double ZCUMMA { get; set; }
        public double WOD { get; set; }
        public double WOV { get; set; }
        public double TD { get; set; }
        public double ALFADK { get; set; }
        public double ALFADL { get; set; }
        public double TV { get; set; }
        public double TCT { get; set; }
        public double TVPG1 { get; set; }
        public double DELTATNAP { get; set; }
        public double ALFAV { get; set; }
        public double lSH { get; set; }
        public double ALFAD { get; set; }
        public double K { get; set; }
        public double QV { get; set; }
        public double F { get; set; }
        public double F1 { get; set; }
        public double FRaz { get; set; }


        public void Calculate()
        {
            CalculateBezPodbora();
            Padbor();
            Graph();
        }

        public void CalculateBezPodbora()
        {
            if (B == 0 || ALFA == 0 || LALFA == 0 || VALFA == 0 || VL == 0 || TRO == 0 || BETTAT == 0 || CV1 == 0 || TV1 == 0 ||
                CV2 == 0 || TV2 == 0 || ETA == 0 || TD1 == 0 || TD2 == 0 || KSIDELTAT == 0 || Z1 == 0 || Z1 == 0 || Z3 == 0 || ETAD == 0 ||
                 ETAV == 0 || DDDK == 0 || ALFADL1 == 0 || KT == 0 || DDVK == 0 || FIBL == 0 || S == 0)
            {
                //MessageBox.Show("Это значение не может равняться нулю");
            }
            else
            {
                VB1 = Math.Round(1.1 * B * LALFA, 1);
                VBPR = Math.Round(VB1 - (B * LALFA), 2);
                VD1 = Math.Round(VALFA * B, 2);
                ID1 = Math.Round(1716 + DELTAT * VL);
                IV1 = Math.Round(CV1 * TV1, 1);
                IV2 = Math.Round(CV2 * TV2);
                ID2 = Math.Round((VD1 * ID1 - (1 / ETA) * (VB1 - VBPR) * (IV2 - IV1)) / (VD1 + VBPR));
                DELTAT1 = Math.Round(((TD1 - TV2) - (TD2 - TV1)) / Math.Log((TD1 - TV2) / (TD2 - TV1)));
                P = Math.Round((TV2 - TV1) / (TD1 - TV1), 2);
                R = Math.Round((TD1 - TD2) / (TV2 - TV1), 2);
                DELTAT = Math.Round(DELTAT1 * KSIDELTAT);
                ZCUMMA = Math.Round(Z1 * Z2 * Z3);
                WOD = Math.Round((VD1 + (0.5 * VBPR)) / ((Z1 - 1) * 0.0274 * (Z2 / 2)), 2);
                WOV = Math.Round((VB1 - 0.5 * VBPR) / (0.013 * (Z1 * Z3) / ETAV), 1);
                TD = Math.Round(0.5 * (TD1 + TD2));
                ALFADK = Math.Round((4.617 + 0.00415 * TD) * ((Math.Pow(3.62, 0.8)) / (Math.Pow(0.158, 0.2))), 1);
                ALFADL = Math.Round(ALFADL1 * KT, 1);
                TV = Math.Round(0.5 * (TV1 + TV2));
                TCT = Math.Round(0.5 * (TD + TV));
                TVPG1 = Math.Round(0.5 * (TCT + TV));
                DELTATNAP = Math.Round(TCT - TV);
                ALFAV = Math.Round((0.9 + (0.00048 * TVPG1)) * ((Math.Pow(WOV, 0.2)) / (Math.Pow(DDVK, 0.5))) * (Math.Pow(DELTATNAP, 0.1)), 1);
                lSH = Math.Round(0.9 + 2.3 * (Math.Pow(10, -4)) * TCT, 2);
                ALFAD = Math.Round(ALFADK + ALFADL, 1);
                K = Math.Round(1 / ((1 / ALFAD) + (S / lSH) + (1 / ALFAV)), 1);
                QV = Math.Round((VB1 - (0.75 * VBPR)) * (IV2 - IV1) * 1000, 2);
                F = Math.Round(QV / (K * DELTAT));
                F1 = Math.Round(ZCUMMA * FIBL);
                FRaz = Math.Round(F1 - F);
            }
        }
        //График
        public static double[] ArrayX;
        public static double[] ArrayY;
        public void Graph()
        {
            int razm = Convert.ToInt32(Math.Round(B + 0.5)) + 2;
            ArrayX = new double[razm];
            ArrayY = new double[razm];

            for (int i = 0; i < razm-1; i++)
            {
                if (B<i)
                {
                    ArrayX[i] = B;
                    ArrayY[i] = Math.Round(VALFA * B, 2);
                    ArrayX[i+1] = i;
                    ArrayY[i+1] = Math.Round(VALFA * i, 2);
                }
                else
                {
                    ArrayX[i] = i;
                    ArrayY[i] = Math.Round(VALFA * i, 2);
                }
            }
        }
        // метод подбора параметров для TV2 
        public void Padbor()
        {
            if (FRaz > 0)
            {
                while (TV2 < 1000)
                {
                    TV2++;
                    CalculateBezPodbora();
                    if (FRaz <= 0)
                        return;
                }
            }
            else if (FRaz < 0)
            {
                while (TV2 > 500)
                {
                    TV2--;
                    CalculateBezPodbora();
                    if (FRaz >= 0)
                        return;
                }
            }
        }
        public ResultModel Rachet()
        {
            return new ResultModel
            {
                VB1 = (double)VB1,
                VBPR = (double)VBPR,
                VD1 = (double)VD1,
                ID1 = (double)ID1,
                IV1 = (double)IV1,
                IV2 = (double)IV2,
                ID2 = (double)ID2,
                DELTAT1 = (double)DELTAT1,
                P = (double)P,
                R = (double)R,
                DELTAT = (double)DELTAT,
                ZCUMMA = (double)ZCUMMA,
                WOD = (double)WOD,
                WOV = (double)WOV,
                TD = (double)TD,
                ALFADK = (double)ALFADK,
                ALFADL = (double)ALFADL,
                TV = (double)TV,
                TCT = (double)TCT,
                TVPG1 = (double)TVPG1,
                DELTATNAP = (double)DELTATNAP,
                ALFAV = (double)ALFAV,
                lSH = (double)lSH,
                ALFAD = (double)ALFAD,
                K = (double)K,
                QV = (double)QV,
                F = (double)F,
                F1 = (double)F1,
                FRaz = (double)FRaz,
                TV2 = (double)TV2,
            };
        }
    }
}
