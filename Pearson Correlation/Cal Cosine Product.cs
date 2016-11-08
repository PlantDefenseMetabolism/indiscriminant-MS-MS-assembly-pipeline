using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Pearson_Correlation {
    class Cal_Cosine_Product {
        //Given 2 spectra, it calculates Cosine Product(Normalized dot product) as a result
        public static double CalCosineProduct(SpectrumData querySpec, SpectrumData librarySpec) {
            double NDP = 0;
            List<SpectrumData> commonSpectrumList = new List<SpectrumData>();
            commonSpectrumList = GetCommonSpectrum(querySpec, librarySpec);
            double a=0, b=0, c=0;
            // in commonSpectrum list, query(0) is in front, library(1) is afterwards
            for (int i = 0; i < commonSpectrumList[0].peakList.Count; i++) {
                a += CalWeight(commonSpectrumList[1].peakList[i].intensity, commonSpectrumList[1].peakList[i].mz) * CalWeight(commonSpectrumList[0].peakList[i].intensity, commonSpectrumList[0].peakList[i].mz);
                //b += Math.Pow(CalWeight(commonSpectrumList[1].peakList[i].intensity, commonSpectrumList[1].peakList[i].mz),2);
                //c += Math.Pow(CalWeight(commonSpectrumList[0].peakList[i].intensity, commonSpectrumList[0].peakList[i].mz),2);
            }
            for (int i = 0; i < querySpec.peakList.Count; i++) {
                b +=Math.Pow( CalWeight(querySpec.peakList[i].intensity, querySpec.peakList[i].mz), 2);
            }
            for (int i = 0; i < librarySpec.peakList.Count; i++) {
                c += Math.Pow( CalWeight(librarySpec.peakList[i].intensity,librarySpec.peakList[i].mz),2);
            }
            if (b*c==0) {
                NDP = 0;
            }
            else NDP = Math.Round( (a * a) / (b * c),3);
            return NDP;
        }

        private static List<SpectrumData> GetCommonSpectrum(SpectrumData querySpec, SpectrumData librarySpec) {
            //commonSpectrum contains 2 items indicating equal length of common peaks between querySpec and librarySpec
            List<SpectrumData> commonSpectrumList = new List<SpectrumData>();
            SpectrumData commonQuerySpec = new SpectrumData();
            commonQuerySpec.peakList = new List<PeakData>();
            SpectrumData commonLibrarySpec = new SpectrumData();
            commonLibrarySpec.peakList = new List<PeakData>();
            for (int i = 0; i <librarySpec.peakList.Count ; i++) {
                for (int j = 0; j < querySpec.peakList.Count; j++) {
                    //if 2 peaks' mz are within 0.01Da, then they will be treated as common peaks between 2 spectra
                    if (IsEqualMZ ( librarySpec.peakList[i].mz,querySpec.peakList[j].mz)) { //judge whether 2 mz values <=0.01
                        commonQuerySpec.peakList.Add(querySpec.peakList[j]);
                        commonLibrarySpec.peakList.Add(librarySpec.peakList[i]);
                        break;
                    }
                }
            }
            // in commonSpectrum list, query is in front, library is afterwards
            commonSpectrumList.Add(commonQuerySpec);
            commonSpectrumList.Add(commonLibrarySpec);
            if (commonQuerySpec.peakList.Count!=commonLibrarySpec.peakList.Count) {
                Console.WriteLine("Common spectrum peak list is not equal between Query and Library!");
            }
            return commonSpectrumList;
        }

        private static double CalWeight(double peakIntensity, double mass) { // W = [peak intensity]m [mass]n  where m=0.5 n=2
            double w = Math.Pow(peakIntensity, 0.5)*Math.Pow(mass,2);
            return w;
        }

        public static bool IsEqualMZ(double mzA,double mzB) {
            if (Math.Abs(mzA - mzB) <= 0.01) {
                return true;
            }
            else return false;
        }
        public static bool IsEqualRT(double RTA, double RTB) {
            if (Math.Abs(RTA - RTB) <= 3) {
                return true;
            }
            else return false;
        }
       
    }
}
