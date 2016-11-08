using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Pearson_Correlation {
    class Cal_Pearson_Correlation {
        // Given 2 equal length vectors, it will return a correlation of the two vectors
        //private double []PrecursorVector { get; set; }
        //private double []FragmentVector { get; set; }
        public static double CalPearsonCorrelation(double[] PrecursorVector, double[] FragmentVector) {
            double cor = new double ();
            double meanPrecursor, meanFragment;
            double sumPrecursor =0, sumFragment =0;
            double a =0, b=0, c=0;
            for (int i = 0; i < PrecursorVector.GetLength(0); i++) {
                sumPrecursor += PrecursorVector[i];
            }
            for (int i = 0; i < FragmentVector.GetLength(0); i++) {
                sumFragment += FragmentVector[i];
            }
            meanPrecursor = sumPrecursor / PrecursorVector.GetLength(0);
            meanFragment = sumFragment / FragmentVector.GetLength(0);
            // Pearson Correlation Equation
            for (int i = 0; i < PrecursorVector.GetLength(0); i++) {
                a += (FragmentVector[i]-meanFragment)*(PrecursorVector[i]-meanPrecursor);
                b += Math.Pow((FragmentVector[i] - meanFragment), 2);
                c += Math.Pow((PrecursorVector[i] - meanPrecursor), 2);
            }
            cor = a / Math.Sqrt(b*c);
            return cor;
        }

        public static List<double[]>HeatMapPearsonCorrelation(string readfilePath, string outpath){
            List<string[]> meanList = FileProcess.ReadCSV(readfilePath);
            List<double[]> pclist = new List<double[]>();
            for (int i = 0; i < meanList.Count; i++) {
                double[] pcline = new double[meanList.Count];
                double[] line = new double[meanList[i].GetLength(0)];
                for (int m = 0; m < meanList[i].GetLength(0); m++) {
                    line[m] = double.Parse(meanList[i][m]);
                }
                for (int j = 0; j < meanList.Count; j++) {
                    double[] lineX = new double[meanList[j].GetLength(0)];
                    for (int n = 0; n < meanList[j].GetLength(0); n++) {
                         lineX[n] = double.Parse(meanList[j][n]);
                    }
                    double pc = Cal_Pearson_Correlation.CalPearsonCorrelation(line,lineX);
                    pcline[j] = pc;
                }
                pclist.Add(pcline);
            }
            return pclist;
        }

        public static List<double[]> HalfHeatMapPearsonCorrelation(List<double[]> pclist) {
            List<double[]> halfpclist = new List<double[]>();
            for (int i = 0; i < pclist.Count; i++) {
                double[] line = new double[i + 1];
                for (int j = 0; j < line.GetLength(0); j++) {
                    line[j] = pclist[i][j];
                }
                halfpclist.Add(line);
            }
            return halfpclist;
        }
    }
}
