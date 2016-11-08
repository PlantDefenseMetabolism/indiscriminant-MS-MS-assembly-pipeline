using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Pearson_Correlation {
    class FileProcess {

        /**
            * Read CSV file to string list
            */
        public static List<string[]> ReadCSV(string csvFilePath) { 
            List<string[]> csvList = new List<string[]>();
            List<string> allFileList = ReadFileToList(csvFilePath);
            for (int j = 0; j < allFileList.Count; j++) {
                string text = allFileList[j];
                var text_array = new List<string[]>();
                var line = new List<string>();
                var field = new StringBuilder();
                //if inside the quotation marks
                bool in_quata = false;
                //if field starts
                bool field_start = true;
                for (int i = 0; i < text.Length; i++) {
                    char ch = text[i];
                    if (in_quata) {
                        //if already inside the quotation marks
                        if (ch == '\"') {
                            //if there are two quotation marks, consider as normal quotation marks
                            if (i < text.Length - 1 && text[i + 1] == '\"') {
                                field.Append('\"');
                                i++;
                            }
                            else
                                //otherwise exit field of quotation marks
                                in_quata = false;
                        }
                        else { //symbles inside the quotation marks（except for quotation marks）consider as normal char

                            field.Append(ch);
                        }
                    }
                    else {
                        switch (ch) {
                            case ',': //new field
                                line.Add(field.ToString());
                                field.Remove(0, field.Length);
                                field_start = true;
                                break;
                            case '\"'://case of quotation marks
                                if (field_start)
                                    in_quata = true;
                                else
                                    field.Append(ch);
                                break;
                            //case '\r': //new record
                            //    if (field.Length > 0 || field_start) {
                            //        line.Add(field.ToString());
                            //        field.Remove(0, field.Length);
                            //    }
                            //    text_array.Add(line.ToArray());
                            //    line.Clear();
                            //    field_start = true;
                            //    //in windows，\r\n  commonly paired，thus skip
                            //    if (i < text.Length - 1 && text[i + 1] == '\n')
                            //        i++;
                            //    break;
                            default:
                                field_start = false;
                                field.Append(ch);
                                break;
                        }
                    }
                }
                //file end
                if (field.Length > 0 || field_start)
                    line.Add(field.ToString());
                //if (line.Count > 0)
                //    text_array.Add(line.ToArray());
                csvList.Add(line.ToArray());
                //string[] line = allFileList[i].Split(',');
                //csvList.Add(line);
            }
            return csvList;
        }

        
        
        public static List<SpectrumData> ReadKnownMSMSSpectraList(string knownMSMSPath) { // read all known compands MSMS data to a list
            List<SpectrumData> allSpectrumList = new List<SpectrumData>();
            List<string> txtList = ReadFileToList(knownMSMSPath);
            for (int i = 0; i < txtList.Count-1; i++) {                
                if (txtList[i].Contains("eV:")) { // indecating sampleName line
                    string[] tittleLine = txtList[i].Split(',');
                    SpectrumData standardSpectrum = new SpectrumData();
                    List<PeakData> spectrumPeakList = new List<PeakData>();
                    for (int j = i+1; j < txtList.Count; j++) {
                        string[] line = txtList[j].Split(default(char[]), StringSplitOptions.RemoveEmptyEntries);
                        if (line.GetLength(0) > 1) {
                            PeakData peak = new PeakData();
                            peak.mz =double.Parse(line[0]);
                            peak.intensity = double.Parse(line[1]);
                            peak.sampleName = txtList[i].Replace(",","_");
                            peak.rt = double.Parse(tittleLine[tittleLine.GetLength(0) - 2].Split(default(char[]), StringSplitOptions.RemoveEmptyEntries)[0]) * 60;
                            spectrumPeakList.Add(peak);
                        }
                        else break;
                    }
                    standardSpectrum.peakList = spectrumPeakList;
                    allSpectrumList.Add(standardSpectrum);
                }
            }
            return allSpectrumList;
        }

        
        private static SpectrumData RemoveLessThan50s(SpectrumData allPeakDataInEachGroup) { // This is to remove all the retention time <50s peaks cause they are calibration peaks
            SpectrumData removed50sSpectrum = new SpectrumData();
            for (int i = 0; i <allPeakDataInEachGroup.peakList.Count; i++) {
                double rt = allPeakDataInEachGroup.peakList[i].rt;
                if (rt<50) {
                    allPeakDataInEachGroup.peakList.RemoveAt(i);
                    i = i - 1;
                }
            }
            removed50sSpectrum = allPeakDataInEachGroup;
            return removed50sSpectrum;
        }

        private static SpectrumData GenerateMergedSpectrum(SpectrumData allPeakDataInEachGroup) { 
            // adMSMS data(20ev-50ev-MS) of each group is treated by first calculating relative intensities among group,
            // then calculating the average intensities in each CID, next taking the max of the value among CIDs
            SpectrumData mergedSpectrum = new SpectrumData();
            mergedSpectrum =CalMaxAverageIntensity(IdMSMS_Deconvolution.CalRelativeIntensity(allPeakDataInEachGroup));
            return mergedSpectrum;
        }

        
        private static SpectrumData CalMaxAverageIntensity(SpectrumData allRelativeIntensity) { // calculate maximum average intensity
            SpectrumData maxAverageSpectrum = new SpectrumData();
            maxAverageSpectrum.peakList = new List<PeakData>();
            List<PeakData> spectrumList = new List<PeakData>();
            var mz = from peaks in allRelativeIntensity.peakList
                     group peaks by peaks.mz;
            foreach (var peaks in mz) {
                PeakData maxPeak = new PeakData();                
                List<PeakData> ev20 = new List<PeakData>(); //store 20ev samples
                List<PeakData> ev30 = new List<PeakData>();
                List<PeakData> ev40 = new List<PeakData>();
                List<PeakData> ev50 = new List<PeakData>();
                List<PeakData> ms = new List<PeakData>(); 
                foreach (var peak in peaks) {
                    if (peak.sampleName.Contains("20eV")) {
                        ev20.Add(peak);
                    }
                    if (peak.sampleName.Contains("30eV")) {
                        ev30.Add(peak);
                    }
                    if (peak.sampleName.Contains("40eV")) {
                        ev40.Add(peak);
                    }
                    if (peak.sampleName.Contains("50eV")) {
                        ev50.Add(peak);
                    }
                    if (peak.sampleName.Contains(".Ms.")) {
                        ms.Add(peak);
                    }
                }
                List<PeakData> mzAverageList = new List<PeakData>(); // store average intensity 
                if (ev20.Count != 0) {
                    ev20[0].intensity = ev20.Average(x => x.intensity);
                    mzAverageList.Add(ev20[0]);
                }
                if (ev30.Count != 0) {
                    ev30[0].intensity = ev30.Average(x => x.intensity);
                    mzAverageList.Add(ev30[0]);
                }
                if (ev40.Count != 0) {
                    ev40[0].intensity = ev40.Average(x => x.intensity);
                    mzAverageList.Add(ev40[0]);
                }
                if (ev50.Count != 0) {
                    ev50[0].intensity = ev50.Average(x => x.intensity);
                    mzAverageList.Add(ev50[0]);
                }
                if (ms.Count != 0) {
                    ms[0].intensity = ms.Average(x => x.intensity);
                    mzAverageList.Add(ms[0]);
                }          

                maxPeak = mzAverageList[0];
                maxPeak.intensity = mzAverageList.Max(x=>x.intensity);
                maxPeak.sampleName = "MergedSpectrum"+maxPeak.pcgroup.ToString();
                spectrumList.Add(maxPeak);
            }
            maxAverageSpectrum.peakList = spectrumList;
            maxAverageSpectrum.group = spectrumList[0].precursor;
            return maxAverageSpectrum;
        }
 

        public static List<string> ReadFileToList(string filePath) { //read file to string list
            using (StreamReader sr = new StreamReader (filePath)) {
                List<string> fileList = new List<string>();
                while (sr.Peek()!=-1) {
                    fileList.Add(sr.ReadLine());
                }
                return fileList;
            }
        }

        public static void WriteNDPListToFile(List<double[]>NDPResultList, string[] queryTitleLine, string[] librarytitleLine,string outputFilePath){ //output NDP list
            using(StreamWriter sw = new StreamWriter(outputFilePath)) {
                sw.Write(",");
                for (int i = 0; i < librarytitleLine.GetLength(0); i++) {
                    sw.Write(librarytitleLine[i] + ",");
                }
                sw.WriteLine();
                for (int i = 0; i < NDPResultList.Count; i++) {
                    sw.Write(queryTitleLine[i]+",");
                    for (int j = 0; j < NDPResultList[i].GetLength(0); j++) {
                        sw.Write(NDPResultList[i][j]+",");
                    }
                    sw.WriteLine();
                }
            }
        }

        public static void WriteNDPResultForSytoscape(List<SpectrumData> allSpecList, string outputPath) {  //output NDP for Sytoscape input
            using (StreamWriter sw = new StreamWriter (outputPath)) {
                for (int i = 0; i < allSpecList.Count-1; i++) {
                    for (int j = i+1; j < allSpecList.Count; j++) {
                        double ndp = Cal_Cosine_Product.CalCosineProduct(allSpecList[i], allSpecList[j]);
                        if (ndp >= 0.01) {
                            List<double> commonNeutralLostList = new List<double>();
                            for (int m = 0; m < allSpecList[i].neutralLossList.Count; m++) {
                                for (int n = 0; n < allSpecList[j].neutralLossList.Count; n++) {
                                    if (Cal_Cosine_Product.IsEqualMZ(allSpecList[i].neutralLossList[m], allSpecList[j].neutralLossList[n])) { // judge whether 2 mz values <=0.01
                                        commonNeutralLostList.Add(Math.Round(allSpecList[i].neutralLossList[m], 2));
                                    }
                                }
                            }
                            commonNeutralLostList = commonNeutralLostList.Distinct().ToList();
                            commonNeutralLostList = commonNeutralLostList.OrderBy(x => x).ToList();
                            sw.Write(allSpecList[i].group.Replace(" ", "") + "," + ndp + "," + allSpecList[j].group.Replace(" ", "") + ",");
                            for (int k = 0; k < commonNeutralLostList.Count; k++) {
                                sw.Write(commonNeutralLostList[k] + ",");
                            }
                            sw.WriteLine();
                        }                                                                   
                    }
                }
            }
        }

        public static void WriteNDPNLHeatmapForSytoscape(string csvPath, string outputPath) { //prepare network data for sytoscape from heapmap by DiffCoEx
            List<string[]> heapmapList = ReadCSV(csvPath);
            using (StreamWriter sw = new StreamWriter(outputPath)) {
                for (int i = 1; i < heapmapList.Count; i++) {
                    for (int j = 1; j < heapmapList[i].GetLength(0); j++) {
                        if (double.Parse(heapmapList[i][j]) >= 0.8) {
                            if (i < j) {//NDP
                                sw.WriteLine(heapmapList[i][0] + " " + "NDP" + " " + heapmapList[0][j] + " " + heapmapList[i][j]);
                            }
                            if (i > j) {//NL
                                sw.WriteLine(heapmapList[i][0] + " " + "NL" + " " + heapmapList[0][j] + " " + heapmapList[i][j]);
                            }
                        }

                        //if (i < j) {//NDP
                        //    sw.WriteLine(heapmapList[i][0] + " (" + "NDP" + ") " + heapmapList[0][j] + " " + "=" + " " + heapmapList[i][j]);
                        //}
                        //if (i > j) {//NL
                        //    sw.WriteLine(heapmapList[i][0] + " (" + "NL" + ") " + heapmapList[0][j] + " " + "=" + " " + heapmapList[i][j]);
                        //}
                    }
                }
            }
            
        }     

        public static void WritePeakList(List<double[]> peakList, string outputFilePath) {
            using (StreamWriter sw = new StreamWriter(outputFilePath)) {
                for (int i = 0; i < peakList.Count; i++) {
                    for (int j = 0; j < peakList[i].GetLength(0); j++) {
                        if (j==peakList[i].GetLength(0)-1) {
                            sw.Write(peakList[i][j]);
                        }
                        else sw.Write(peakList[i][j]+",");
                    }
                    sw.WriteLine();
                }
            }
        }
        public static void WritePeakList(List<List<double>> peakList, string outputFilePath) {
            using (StreamWriter sw = new StreamWriter(outputFilePath)) {
                for (int i = 0; i < peakList.Count; i++) {
                    for (int j = 0; j < peakList[i].Count; j++) {
                        if (j == peakList[i].Count - 1) {
                            sw.Write(peakList[i][j]);
                        }
                        else sw.Write(peakList[i][j] + ",");
                    }
                    sw.WriteLine();
                }
            }
        }
        public static void WritePeakList(List<string[]> peakList, string outputFilePath) {
            using (StreamWriter sw = new StreamWriter(outputFilePath)) {
                for (int i = 0; i < peakList.Count; i++) {
                    for (int j = 0; j < peakList[i].GetLength(0); j++) {
                        if (j == peakList[i].GetLength(0) - 1) {
                            sw.Write(peakList[i][j]);
                        }
                        else sw.Write(peakList[i][j] + ",");
                    }
                    sw.WriteLine();
                }
            }
        }

        public static void WriteSpectrumList(List<SpectrumData> specList, string outputPath) { //output spectra list
            using (StreamWriter sw = new StreamWriter(outputPath)) {
                sw.WriteLine("pvalue" + "," + "mz" + "," + "rt" + "," + "inten" + "," + "iso" + "," + "add" + "," + "pcgroup"+","+"precursor");
                for (int i = 0; i < specList.Count; i++) {
                    for (int j = 0; j < specList[i].peakList.Count; j++) {
                        sw.WriteLine(specList[i].peakList[j].pval + "," + specList[i].peakList[j].mz + "," + specList[i].peakList[j].rt + "," + specList[i].peakList[j].intensity + "," + specList[i].peakList[j].isotopes + "," +
                            specList[i].peakList[j].adduct + "," + specList[i].peakList[j].pcgroup+","+specList[i].peakList[j].precursor);
                    }
                }
            }
        }

        public static void OutputNeutralLoss(List<SpectrumData> specList, string outputPath) {// output results contain 3 columns: "hit number; NL; all hit spectra"
            List<double> allNeutralList = new List<double>();
            for (int i = 0; i < specList.Count; i++) {
                for (int j = 0; j < specList[i].neutralLossList.Count; j++) {
                    allNeutralList.Add(Math.Round(specList[i].neutralLossList[j], 2));
                }
            }
            allNeutralList = allNeutralList.Distinct().ToList();
            using (StreamWriter sw = new StreamWriter(outputPath)) {
                sw.WriteLine("hit number"+","+"target NL"+","+"all hit spectra");
                for (int i = 0; i < allNeutralList.Count; i++) {
                    List<string> hitSpecList = new List<string>();
                    for (int j = 0; j < specList.Count; j++) {
                        for (int k = 0; k < specList[j].neutralLossList.Count; k++) {
                            if (allNeutralList[i] == Math.Round(specList[j].neutralLossList[k], 2)) {
                                hitSpecList.Add(specList[j].group);
                            }
                        }
                    }
                    hitSpecList = hitSpecList.Distinct().ToList();
                    sw.Write(hitSpecList.Count + ",");
                    sw.Write(allNeutralList[i] + ",");
                    for (int j = 0; j < hitSpecList.Count; j++) {
                        sw.Write(hitSpecList[j] + ",");
                    }
                    sw.WriteLine();
                }
            }
        }
        
        public static void OutputSelectedNLto360Spec(List<SpectrumData> specList, string NLpath, string outputPath) {// this function is to output idMSMS spectra which hit to selected NLs for DiffCoEx(0 or 1); col: selected NL, row: idMSMSspectra
            List<string[]> SelectedNLList = ReadCSV(NLpath); //columnwise
            string [] allNL = new string[SelectedNLList.Count];
            for (int i = 0; i < SelectedNLList.Count; i++) {
                allNL[i] = SelectedNLList[i][0];
            }
            List<string[]> hitList = new List<string[]>();
            hitList.Add(allNL);
            for (int i = 0; i <specList.Count ; i++) {
                string[] row = new string[allNL.GetLength(0)+1];
                row[allNL.GetLength(0)] = specList[i].group;
                List<double> specNLList = specList[i].neutralLossList;
                for (int j = 0; j < allNL.GetLength(0); j++) {
                    string nl = allNL[j];
                    row[j] = "0";
                    for (int k = 0; k < specNLList.Count; k++) {
                        if (Math.Abs(double.Parse(nl)-specNLList[k])<=0.01) {
                            row[j] = "1";
                            break;
                        }
                    }
                }
                hitList.Add(row);
                WritePeakList(hitList,outputPath);
            }
        }

        public static void WriteSpectrumQueryForMassBank(List<SpectrumData> specList, string outputPath) { //output Massbank query format
            using (StreamWriter sw = new StreamWriter(outputPath)) {
                for (int i = 0; i < specList.Count; i++) {
                    sw.WriteLine("Name: "+specList[i].group);
                    for (int j = 0; j < specList[i].peakList.Count; j++) {
                        sw.WriteLine(specList[i].peakList[j].mz + " " + specList[i].peakList[j].intensity);
                    }
                    sw.WriteLine();
                }
            }
        }

      

    }
}
