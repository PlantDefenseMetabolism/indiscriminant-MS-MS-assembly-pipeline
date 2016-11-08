using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Pearson_Correlation {
    class IdMSMS_Deconvolution {
        //parameters:
        static double bg = 2.2235; //background value(resulted from fill 0 values in XCMS peaktable)
        public static int weight =3; // times higher than bg noise is considered as a true peak
        public static int minSampleNrForCor =6; //the minimum number of samples within all samples for correlation calculation
        public static double ndpThreshold = 0.9;//if NDP>threshold in same pcgroup, merge sub-idMS/MS to one idMS/MS
        public static double rtWindow = 9;//  rt window to precursor
        //

        public static List<SpectrumData> MergeIdMSMSOfDiffCIDs(string rootDataPath) {// using idMSMS precursor-product method by calculated pvals, this method will merge 20-50ev data by selecting highest intensity(the peaks are different among CIDs)
            List<SpectrumData> mergedSpecList = new List<SpectrumData>();
            List<SpectrumData> mergedCIDList = new List<SpectrumData>();
            List<SpectrumData> AllCIDList = new List<SpectrumData>(); // all CID
            List<SpectrumData> CID20evList = ProcessPeaktableToSpectralist(rootDataPath.Replace("MS", "20ev"));
            AllCIDList.AddRange(CID20evList.Where(x => double.Parse(x.group.Split('_')[1]) > 50).ToList()); //remove rt<50
            CID20evList = null;
            GC.Collect();
            Console.WriteLine("20ev done!");
            List<SpectrumData> CID30evList = ProcessPeaktableToSpectralist(rootDataPath.Replace("MS", "30ev")); 
            AllCIDList.AddRange(CID30evList.Where(x => double.Parse(x.group.Split('_')[1]) > 50).ToList());
            CID30evList = null;
            GC.Collect();
            Console.WriteLine("30ev done!");
            List<SpectrumData> CID40evList = ProcessPeaktableToSpectralist(rootDataPath.Replace("MS", "40ev"));
            AllCIDList.AddRange(CID40evList.Where(x => double.Parse(x.group.Split('_')[1]) > 50).ToList());
            CID40evList = null;
            GC.Collect();
            Console.WriteLine("40ev done!");
            List<SpectrumData> CID50evList = ProcessPeaktableToSpectralist(rootDataPath.Replace("MS", "50ev"));
            AllCIDList.AddRange(CID50evList.Where(x => double.Parse(x.group.Split('_')[1]) > 50).ToList());
            CID50evList = null;
            GC.Collect();
            Console.WriteLine("50ev done!");
            List<SpectrumData> MSList = ProcessPeaktableToSpectralist(rootDataPath); // Note: pcgroup is different among different CIDs
            List<PeakData> allMSPeakList = new List<PeakData>();
            for (int i = 0; i < MSList.Count; i++) {
                allMSPeakList.AddRange(MSList[i].peakList.Where(x => x.rt > 50).ToList()); //already averaged, meaning that each PeakData=one mz line 
            }         
            MSList = null;
            GC.Collect();
            var allSpec = from specs in AllCIDList
                          group specs by specs.peakList[0].precursor;
            foreach (var specs in allSpec) {
                PeakData precursorPeak = allMSPeakList.Where(x => x.mz == double.Parse(specs.Key.Split('_')[0])).ToList()[0]; // find this precursor in MS
                double precursorPeakRT = precursorPeak.rt;
                /**   
                 * Filter presursor: 1,precursor should pass the "noise"threshold; 2, isotope is removed from precursor list（but isotopes can be fragments; 3, precursors will meet the criteria that "diffrt"(retention time difference) and "feature"(mass difference) are zero 
                 */
                List<PeakData> peaksofeachprecursor = new List<PeakData>();
                for (int k = 0; k < specs.ToList().Count; k++) {
                    peaksofeachprecursor.AddRange(specs.ToList()[k].peakList);
                }
                if (precursorPeak.isNoise == true || IsIsotope(precursorPeak.isotopes) || !IsSelfCor(peaksofeachprecursor)) { // determine whether precursor peak contains minimum number of "minSampleNrForCor" or is isotope, if not, not include
                    continue;
                }

                /**   
                 * Filter fragments: 1, fragments should not be "noise" peaks，2, fragment's "sample position" should be the subcollection of precursor's position，3, fragment's m/z less than precursor（The highest intensity of fragment is selected eventually，and then calculate its relativeintensity）  
                 */
                List<string> precursorPeakPosition = precursorPeak.sampleName.Split('_').ToList();
                double precursorMZ = precursorPeak.mz;
                List<PeakData> allPeaksInOneGroup = new List<PeakData>();
                List<PeakData> MList = new List<PeakData>();
                foreach (var spec in specs) {
                    allPeaksInOneGroup.AddRange(spec.peakList.Where(x => (x.isNoise == false && x.diffrt < rtWindow)).ToList());// determine whether fragment peak contains minimum number of "minSampleNrForCor", if not, not include
                }
                for (int i = 0; i < allPeaksInOneGroup.Count; i++) {
                    List<string> fragmentPeakPosition = allPeaksInOneGroup[i].sampleName.Split('_').ToList();
                    double fragmentMZ = allPeaksInOneGroup[i].mz; //until here, the higher mz is not removed, so fragment mz can be higher than precursor mz here                  
                    if (precursorPeakPosition.Intersect(fragmentPeakPosition).Count() < minSampleNrForCor) { // min samples for correlation
                        if (allPeaksInOneGroup[i].isotopes.Contains("[M]")) {//after saving [M]+, check whether it is in the isotope list，if so, add this [M]+
                            MList.Add(allPeaksInOneGroup[i]);
                        }
                        allPeaksInOneGroup.RemoveAt(i);
                        i = i - 1;
                    }
                }
                for (int i = 0; i < MList.Count; i++) {
                    for (int j = 0; j < allPeaksInOneGroup.Count; j++) {
                        if (IsMainPeak(allPeaksInOneGroup[j], MList[i])) {
                            allPeaksInOneGroup.Add(MList[i]);
                            break;
                        }
                    }
                }
                SpectrumData maxSpectrum = new SpectrumData();
                maxSpectrum.peakList = allPeaksInOneGroup;
                maxSpectrum.group = specs.Key;
                List<SpectrumData> sd = new List<SpectrumData>();
                sd.Add(maxSpectrum);
                maxSpectrum = MergeSpectraToOne(sd); // merge 20-50ev
                maxSpectrum.peakList = maxSpectrum.peakList.OrderBy(x => x.mz).ToList();
                if (maxSpectrum.peakList.Count != 0) {
                    mergedCIDList.Add(CalRelativeIntensity(maxSpectrum)); //calculate relative intensity 
                }
            }
            AllCIDList = null;
            GC.Collect();
            Console.WriteLine("Filter done!");
            /**   
             * merging sub-idMS/MS function: merge those idMSMSs which are in the same pcgroup && ndp>thresholdNDP， then remove higher mz than precursor, calculate relative intensity， calculate neutral loss 
             */
            int n = 0;
            var pcgrps = mergedCIDList.GroupBy(x => int.Parse(x.group.Split('_')[2])); //same pcgroup
            foreach (var specs in pcgrps) {
                List<SpectrumData> specList = specs.ToList();
                if (specList.Count > 1) {//check if the pcgroup contains nore than one idMSMS
                    List<SpectrumData> similarList = new List<SpectrumData>();
                    specList = specList.OrderBy(x => double.Parse(x.group.Split('_')[0])).ToList();//ascending order mz, merge from bottom up
                    for (int i = 0; i < specList.Count - 1; i++) {
                        PeakData precursorPeak = allMSPeakList.Where(x => x.mz == double.Parse(specList[i].group.Split('_')[0])).ToList()[0]; // find this precursor in MS
                        if (!Data_Operation.IsPrecursor(precursorPeak.adduct)) { // if the mz is not precursor, merge upwards
                            for (int j = i + 1; j < specList.Count; j++) {
                                double ndp = 0;
                                ndp = Cal_Cosine_Product.CalCosineProduct(specList[i], specList[j]);
                                if (ndp > ndpThreshold) {
                                    List<SpectrumData> candiList = new List<SpectrumData>();
                                    candiList.Add(specList[i]);
                                    candiList.Add(specList[j]);
                                    SpectrumData MeSpec = MergeSpectraToOne(candiList);
                                    specList[j] = MeSpec; //merge a and b (merged result is stored in b and delete a)
                                    specList.RemoveAt(i);
                                    i--;
                                    n++;
                                    break;
                                }
                            }
                        }
                    }
                   
                }
                for (int i = 0; i < specList.Count; i++) {
                    SpectrumData removeSpec = RemoveMZHigherThanPrecursor(specList[i]); //remove higher mz than precursor                    
                    if (removeSpec.peakList.Count != 0) {
                        SpectrumData finalSpec = CalRelativeIntensity(removeSpec);// cal relative intensity
                        finalSpec.peakList = finalSpec.peakList.OrderByDescending(x => x.mz).ToList();
                        finalSpec.neutralLossList = Data_Operation.FindNeutralLoss(finalSpec);
                        mergedSpecList.Add(finalSpec); 
                    }
                }
            }
            Console.WriteLine("n:"+n);
            Console.WriteLine("merge done!");
            return mergedSpecList;
        }

        /**   
          * Merge multiple idMSMS (if mz difference less than 0.01, pick the higher intensity one)
         * This function is used in two cases：1, merge 20-50ev idMSMSs， 2, merge idMSMSs in the same pcgroup with NDP>threshold
          */
        private static SpectrumData MergeSpectraToOne(List<SpectrumData>spectraWithSamePcgrp) { 
            SpectrumData maxSpectrum = new SpectrumData();
            maxSpectrum.peakList = new List<PeakData>();
            maxSpectrum.group = spectraWithSamePcgrp.OrderByDescending(x => double.Parse(x.group.Split('_')[0])).ToList()[0].group; //choose the highest mz as precursor(if this is 20-50ev, group is the same, then no influnce)
            List<PeakData> allPeaksInOneGroup = new List<PeakData>();
            for (int i = 0; i < spectraWithSamePcgrp.Count; i++) {
                allPeaksInOneGroup.AddRange(spectraWithSamePcgrp[i].peakList);
            }
            var allpeaks = from peaks in allPeaksInOneGroup
                           group peaks by peaks.mz;
            foreach (var peaks in allpeaks) {
                PeakData maxPeak = peaks.OrderByDescending(x => x.intensity).ToList()[0]; //choose the highest intensity peak in one mz
                maxSpectrum.peakList.Add(maxPeak);
            }
            for (int m = 0; m < maxSpectrum.peakList.Count; m++) { //remove redundancy(two mz <0.01)
                for (int n = m + 1; n < maxSpectrum.peakList.Count; n++) {
                    if (Cal_Cosine_Product.IsEqualMZ(maxSpectrum.peakList[m].mz, maxSpectrum.peakList[n].mz)) {
                        if (maxSpectrum.peakList[m].intensity >= maxSpectrum.peakList[n].intensity) {
                            maxSpectrum.peakList.RemoveAt(n);
                            n = n - 1;
                        }
                        else {
                            maxSpectrum.peakList[m] = maxSpectrum.peakList[n];
                            n = n - 1;
                        }
                    }
                }
            }
            for (int i = 0; i < maxSpectrum.peakList.Count; i++) { //assign max precursor to all peaks
                maxSpectrum.peakList[i].precursor = maxSpectrum.group;
            }
            maxSpectrum.peakList = maxSpectrum.peakList.OrderBy(x => x.mz).ToList();
            return maxSpectrum;
        }

        private static SpectrumData RemoveMZHigherThanPrecursor(SpectrumData spec) {
            SpectrumData removedSpec = spec;
            List<PeakData> peaks = removedSpec.peakList;
            peaks = peaks.Where(x => x.mz <= double.Parse(removedSpec.group.Split('_')[0])).ToList();
            removedSpec.peakList = peaks;
            return removedSpec;
        }


        /**   
         * After using this function，list contains a <PeakData> List， each index in the list coresponds to one mz （value averaged）  
         */
        public static List<SpectrumData> ProcessPeaktableToSpectralist(string csvFilePath) { //read result CSV file generated by XCMS and CAMERA to a list of spectra
            List<SpectrumData> allSpectrumList = new List<SpectrumData>();
            List<PeakData> allPeakList = new List<PeakData>();
            List<string> CSVList = new List<string>();
            CSVList = FileProcess.ReadFileToList(csvFilePath);
            //CSVList = BigFileReader.GetFileContent(csvFilePath);
            string[] titleLine = null;
            if (CSVList.Count != 0) { // extract tittle row
                titleLine = CSVList[0].Split(',');
            }
            else Console.WriteLine("CSVFile is empty!");
            int diffrtIndex = 0, pvalIndex = 0, featureIndex = 0, mzIndex = 4, rtIndex = 5, 
                isotopesIndex = titleLine.GetLength(0) - 4, adductIndex = titleLine.GetLength(0) - 3, pcgroupIndex = titleLine.GetLength(0) - 2, precursorIndex = titleLine.GetLength(0) - 1;
            //store the index of interested tittles, allocate initial default values.
            for (int j = 0; j < titleLine.GetLength(0); j++) {
                titleLine[j] = titleLine[j].Replace("\"", "");
                switch (titleLine[j]) {
                    case "mz":
                        mzIndex = j;
                        break;
                    case "rt":
                        rtIndex = j;
                        break;
                    case "isotopes":
                        isotopesIndex = j;
                        break;
                    case "adduct":
                        adductIndex = j;
                        break;
                    case "pcgroup":
                        pcgroupIndex = j;
                        break;
                    case "pval":
                        pvalIndex = j;
                        break;
                    case "feature":
                        featureIndex = j;
                        break;
                    case "diffrt":
                        diffrtIndex = j;
                        break;
                    case "precursor":
                        precursorIndex = j;
                        break;
                    default:
                        break;
                }
            }

            for (int i = 1; i < CSVList.Count; i++) { // extract data rows, each time extract 1 row and take title row as a pair.
                string[] line = CSVList[i].Replace("\"", "").Split(',');
                List<PeakData> mzPeaks = new List<PeakData>();
                for (int j = 0; j < line.GetLength(0); j++) {
                    if (titleLine[j].Contains("_") && titleLine[j].Length > 1) { // all the sample names contain "_" char
                        PeakData peak = new PeakData();
                        peak.mz = double.Parse(line[mzIndex]);
                        peak.rt = double.Parse(line[rtIndex]);
                        peak.isotopes = line[isotopesIndex];
                        peak.adduct = line[adductIndex];
                        peak.pcgroup = int.Parse(line[pcgroupIndex]);
                        peak.intensity = double.Parse(line[j]);
                        peak.sampleName = titleLine[j];
                        if (diffrtIndex != 0) {//MS doesnt contain pval, feature, diffrt
                            peak.pval = double.Parse(line[pvalIndex]);
                            peak.feature = double.Parse(line[featureIndex]);
                            peak.diffrt = double.Parse(line[diffrtIndex]);                         
                        }
                        peak.precursor = line[precursorIndex].Replace(" ", "");    
                        if (double.Parse(line[j]) >= weight * bg) { //intensity which is 3 times higher than the average noise will be deemed as real peak! 
                            peak.isNoise = false; 
                        }
                        else peak.isNoise = true;
                        mzPeaks.Add(peak);
                    }                 
                }
                PeakData averPeak = AverageIntensityofMZline(mzPeaks);
                allPeakList.Add(averPeak);// after AverageIntensityofMZline, only one peakdata left
            }
            var groupResult = from p in allPeakList
                              group p by p.precursor;
            foreach (var groups in groupResult) {
                List<PeakData> peakList = new List<PeakData>();
                peakList = groups.ToList<PeakData>();
                SpectrumData spec = new SpectrumData();
                spec.peakList = peakList; //spec contains all original peak data
                spec.group = groups.Key.Replace(" ", "");
                if (peakList.Count != 0) { //0  means no peaks, 1means single peak group, should be removed
                    allSpectrumList.Add(spec); //here merged all samples and CIDs into allSpectrumList
                }
                else continue;
            }
            allSpectrumList = allSpectrumList.OrderBy(x => x.group).ToList();
            return allSpectrumList;
        }

        private static PeakData AverageIntensityofMZline(List<PeakData> mzPeakList) { // this method is to calculate average intensity for each mz line with the one 3 times higher than noise in one of the 20-50ev spectrum(calculated by pvals). to create a list containing each mz as list[i] with position of this mz
            PeakData averPeak = new PeakData();
            double sumIntensity = 0, averageIntensity = 0;
            int n = 0;
            List<string> peakPositionList = new List<string>();
            StringBuilder peakPosition = new StringBuilder();
            for (int i = 0; i < mzPeakList.Count; i++) {
                if (mzPeakList[i].isNoise == false) {
                    sumIntensity += mzPeakList[i].intensity;
                    peakPositionList.Add(mzPeakList[i].sampleName.Split('_')[0]); 
                    n++;
                }
            }
            averPeak = mzPeakList[0]; // choose the first peak in each mz line
            averPeak.isNoise = true;
            if (n != 0) {
                averageIntensity = sumIntensity / n;
                if (n >= minSampleNrForCor) {// minimun samples for correlation analysis
                    averPeak.isNoise = false;
                }
            }
            for (int i = 0; i < peakPositionList.Count; i++) {
                if (i == peakPositionList.Count - 1) {
                    peakPosition.Append(peakPositionList[i]);
                }
                else peakPosition.Append(peakPositionList[i] + "_");
            }
            averPeak.intensity = averageIntensity; 
            averPeak.sampleName = peakPosition.ToString(); //name formats as 53_196_115
            return averPeak;
        }
        private static SpectrumData AverageIntensityofMZline(SpectrumData allPeakDataInEachGroup) { // this method is to calculate average intensity for each mz line with the one 3 times higher than noise in one of the 20-50ev spectrum(calculated by pvals). to create a list containing each mz as list[i] with position of this mz
            SpectrumData averageSpectrum = new SpectrumData();
            averageSpectrum.peakList = new List<PeakData>();
            var mzLines = from peaks in allPeakDataInEachGroup.peakList
                          group peaks by peaks.mz;
            foreach (var mzline in mzLines) {
                double sumIntensity = 0, averageIntensity = 0;
                int n = 0;
                List<string> peakPositionList = new List<string>();
                StringBuilder peakPosition = new StringBuilder();
                foreach (var peak in mzline) {
                    if (peak.isNoise == false) {
                        sumIntensity += peak.intensity;
                        peakPositionList.Add(peak.sampleName.Split('_')[0]); 
                        n++;
                    }
                }
                PeakData averPeak = mzline.ToList<PeakData>()[0]; // choose the first peak in each mz line
                averPeak.isNoise = true;
                if (n != 0) {
                    averageIntensity = sumIntensity / n;
                    if (n >= minSampleNrForCor) {// minimun samples for correlation analysis
                        averPeak.isNoise = false;
                    }
                }
                for (int i = 0; i < peakPositionList.Count; i++) {
                    if (i == peakPositionList.Count - 1) {
                        peakPosition.Append(peakPositionList[i]);
                    }
                    else peakPosition.Append(peakPositionList[i] + "_");
                }
                averPeak.intensity = averageIntensity; 
                averPeak.sampleName = peakPosition.ToString(); //name formats as 53_196_115
                averageSpectrum.peakList.Add(averPeak);
            }
            averageSpectrum.group = allPeakDataInEachGroup.group;
            return averageSpectrum;
        }

        public static SpectrumData CalRelativeIntensity(SpectrumData allPeakDataInEachGroup) { 
            SpectrumData relativeSpectrum = allPeakDataInEachGroup;
            double maxIntensity = relativeSpectrum.peakList.Max(x => x.intensity);
            for (int i = 0; i < relativeSpectrum.peakList.Count; i++) {
                relativeSpectrum.peakList[i].intensity = 100 * relativeSpectrum.peakList[i].intensity / maxIntensity;
            }
            return relativeSpectrum;
        }

        private static SpectrumData RemoveIsotopes(SpectrumData allPeakDataInEachGroup) { // This is to remove all the isotopes in each group [M+1]+[M+2]+[M+3]+[M+4]+ in order to raise the NDP
            SpectrumData removedIsotopesSpectrum = new SpectrumData();
            for (int i = 0; i < allPeakDataInEachGroup.peakList.Count; i++) {
                string isotopes = allPeakDataInEachGroup.peakList[i].isotopes;
                if (IsIsotope(isotopes)) {
                    allPeakDataInEachGroup.peakList.RemoveAt(i);
                    i = i - 1;
                }
            }
            removedIsotopesSpectrum = allPeakDataInEachGroup;
            return removedIsotopesSpectrum;
        }

        public static bool IsIsotope(string isotopes) {
            if (isotopes.Contains("[M+1]") || isotopes.Contains("[M+2]") || isotopes.Contains("[M+3]") || isotopes.Contains("[M+4]")) {
                return true;
            }
            else return false;
        }

        /**   
          * Function as one of the precursor filters: precursors will meet the criteria that "diffrt"(retention time difference) and "feature"(mass difference) are zero 
          */
        private static bool IsSelfCor(List<PeakData>allpeak) {
            int n=0;
            for (int i = 0; i < allpeak.Count; i++) {
                if (allpeak[i].feature==0&allpeak[i].diffrt==0) {
                    n=1;
                    break;
                }
            }
            if (n == 1) {
                return true;
            }
            else return false;
        }

        private static bool IsMainPeak(PeakData iso, PeakData main) {
            bool reslt = false;
            if (IsIsotope(iso.isotopes) && main.isotopes.Contains("[M]")) {
                int n = int.Parse(iso.isotopes.Split('+')[1].Replace("]", ""));
                if (((iso.mz - main.mz) < 1.01 * n) && (iso.mz - main.mz) > 1 * n) {
                    reslt = true;
                }
            }
            return reslt;
        }
    }
}
