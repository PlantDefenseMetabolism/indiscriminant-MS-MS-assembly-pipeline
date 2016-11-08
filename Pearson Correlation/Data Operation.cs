using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Pearson_Correlation {
    class Data_Operation {
        /* This method is to fill those 0 values from XCMS peaktable ( after "fill peaks" function from XCMS, some of the intensities are still 0 )
            This is done by generating 2 peaktables with file1  with "fillpeaks" function, file2 not, then sum all the "NA" areas in file1 and take average then assign to zero area.
         //*/
        public static void Fill0Peaks(string fillPeaksDataPath, string noneFillPeaksDataPath, string outputFilePath) {
            List<string[]> fill0ValueList = new List<string[]>();
            List<string> fillPeaksDataList = FileProcess.ReadFileToList(fillPeaksDataPath);
            List<string> noneFillPeaksDataList = FileProcess.ReadFileToList(noneFillPeaksDataPath);
            if (fillPeaksDataList.Count != noneFillPeaksDataList.Count) {
                Console.WriteLine("Fillpeaks and noneFillpeaks data are not match the same length! check the original data!");
                return;
            }
            string[] titleLine_fill = null; string[] titleLine_nonefill = null;
            if (fillPeaksDataList.Count != 0) { // extract tittle row
                titleLine_fill = fillPeaksDataList[0].Split(',');
                titleLine_nonefill = noneFillPeaksDataList[0].Split(',');
            }
            double sumBackgroundNoise = 0, BackgroundNoise = 0;
            int NAcount = 0;
            for (int i = 1; i < noneFillPeaksDataList.Count; i++) {
                string[] line = noneFillPeaksDataList[i].Split(',');
                for (int j = 0; j < line.GetLength(0); j++) {
                    //if (titleLine[j].Replace("\"", "").StartsWith("X") &&!titleLine[j].Replace("\"", "").Contains("X.") && line[j] == "NA") {
                    if (titleLine_nonefill[j].Contains("_") && line[j] == "NA") {
                        for (int k = 0; k < line.GetLength(0); k++) {
                            if (titleLine_nonefill[j].Replace("\"", "") == titleLine_fill[k].Replace("\"", "")) {
                                sumBackgroundNoise += double.Parse(fillPeaksDataList[i].Split(',')[k]);
                                NAcount++;
                                break;
                            }
                        }

                    }
                }
            }
            if (NAcount != 0) {
                BackgroundNoise = sumBackgroundNoise / NAcount;
            }
            else BackgroundNoise = 0;
            for (int i = 0; i < fillPeaksDataList.Count; i++) {
                string[] line = fillPeaksDataList[i].Split(',');
                for (int j = 0; j < line.GetLength(0); j++) {
                    //if (titleLine[j].Replace("\"", "").StartsWith("X") && !titleLine[j].Replace("\"", "").Contains("X.") && line[j] == "0") {
                    if (titleLine_fill[j].Contains("_") && line[j] == "0") {
                        line[j] = BackgroundNoise.ToString();
                    }
                }
                fill0ValueList.Add(line);
            }
            Console.WriteLine("BackgroundNoise: " + BackgroundNoise);
            FileProcess.WritePeakList(fill0ValueList, outputFilePath);
        }


        /* This method is to extract all the precursors' mz values.
            [M+H]+; [M+Na]+; [M+K]+; [M+NH4]+
         //*/
        public static void FindPrecursorMZ(List<SpectrumData> allSpectraList, string outputFilePath) {
            List<double> mzList = new List<double>();
            for (int i = 0; i < allSpectraList.Count; i++) {
                for (int j = 0; j < allSpectraList[i].peakList.Count; j++) {
                    string iso =allSpectraList[i].peakList[j].isotopes;
                    if (!IdMSMS_Deconvolution.IsIsotope(iso)) {
                    //if (true) {
                        mzList.Add(allSpectraList[i].peakList[j].mz);
                    }
                }
            }
            using (StreamWriter sw = new StreamWriter (outputFilePath)) {
                for (int i = 0; i < mzList.Count; i++) {
                    if (i == mzList.Count - 1) {
                        sw.Write(mzList[i]);
                    }
                    else {
                        sw.Write(mzList[i] + ",");
                        //if (i%300==0&&i!=0) {
                        //    sw.WriteLine();
                        //}
                    }
                }
            }
            Console.WriteLine(mzList.Count);
        }
        
        /* This method is to calculate all the neutral loss in each spectrum.
            cross over all peaks!
         //*/
        public static List<double> FindNeutralLoss(SpectrumData spectrum) {
            List<double> neutralLossList = new List<double>();
            List<PeakData> sortedPeakList = spectrum.peakList.OrderByDescending(x=>x.mz).Where(x=>!IdMSMS_Deconvolution.IsIsotope(x.isotopes)).ToList(); //并计算和isotope的NL
            double mzPrecursor = double.Parse(spectrum.group.Split('_')[0]);
            for (int i = 0; i < sortedPeakList.Count; i++) { // add precursor neutral loss
                double neutralLoss = mzPrecursor - sortedPeakList[i].mz;
                if (neutralLoss > 14) {
                    neutralLossList.Add(neutralLoss);
                }     
            }
            for (int i = 0; i < sortedPeakList.Count; i++) { // add peaks neutral loss
                for (int j = i+1; j < sortedPeakList.Count; j++) {
                    double neutralLoss = sortedPeakList[i].mz - sortedPeakList[j].mz;
                    if (neutralLoss>14) {
                        neutralLossList.Add(neutralLoss);
                    }                    
                }
            }
            neutralLossList = neutralLossList.Distinct().ToList();
            return neutralLossList;
        }

        public static bool IsPrecursor(string adduct) {
            if (adduct.Contains("[M+H]+") || adduct.Contains("[M+Na]+") || adduct.Contains("[M+K]+") || adduct.Contains("[M+NH4]+") || adduct.Contains("[M+H+NH3]+")) {
                return true;
            }
            else return false;
        }


    }
}
