using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;

namespace Pearson_Correlation {
    class Program {
        static void Main(string[] args) {
            Stopwatch watch = new Stopwatch();
            watch.Start();
            string rootPath = @"E:\Different tissues\idMSMS different tissue\MIX Results\"; 

            /**
             * Fill 0 values in XCMS peaktable
             */
            //string fillpeaksPath = rootPath + "NormPeakArea.csv";
            //string nonefillpeaksPath = rootPath + "data_minfrac0.1_sn3_nonfill.csv";
            //string outPath = rootPath + "NormPeakArea_fill0.csv";
            //Data_Operation.Fill0Peaks(fillpeaksPath, nonefillpeaksPath, outPath);

            /**
             * export precursors (change parapeters to export all/without isotope/only adducts)
             */
            //string CSVFilePath = rootPath + "NormPeakArea.csv";
            //List<SpectrumData> allSpectraList = IdMSMS_Deconvolution.ProcessPeaktableToSpectralist(CSVFilePath);
            //Data_Operation.FindPrecursorMZ(allSpectraList, outPrecursorPath);

            /** 
             *Assemble idMS/MS with multi-eV  
             */
            string CSVFilePath = rootPath + "AllDeconvolutedSpectra_MS_pvals_9s.csv"; 
            string outputFilePath = rootPath + "idmsms_3_6_0.9_minCor_rmv50_50ev.csv";
            string NLOutputFilePath = rootPath+"NL_motif.csv";
            string NDPmatrixPath = rootPath + "NDPmatrix.csv";
            List<SpectrumData> allMergedList = IdMSMS_Deconvolution.MergeIdMSMSOfDiffCIDs(CSVFilePath);          
            FileProcess.WriteSpectrumList(allMergedList, outputFilePath); //output idMSMS spectra
            FileProcess.OutputNeutralLoss(allMergedList, NLOutputFilePath); //output NL_motif
            string commonNLPath = rootPath + "Common NL table.csv";
            string outNLHitPath = rootPath + "NLmatrix.csv";
            FileProcess.OutputSelectedNLto360Spec(allMergedList, commonNLPath, outNLHitPath); //output binary NL hit table
            //List<SpectrumData> allSpectraList = FileProcess.ReadSpectralCSVFile(CSVFilePath); // use if there is spectra list to read
            //List<SpectrumData> allKnownSpectraList = FileProcess.ReadKnownMSMSSpectraList(knownMSMSPath); // use if there is known spectra list
            List<SpectrumData> allKnownSpectraList = allMergedList;
            List<double[]> NDPResultList = new List<double[]>();
            string[] queryTitleLine = new string[allMergedList.Count]; // add QueryTitle line to result
            for (int i = 0; i < allMergedList.Count; i++) {
                queryTitleLine[i] = allMergedList[i].group.ToString();
            }
            string[] libraryTitleLine = new string[allKnownSpectraList.Count]; // add Library title to result
            for (int i = 0; i < allKnownSpectraList.Count; i++) {
                libraryTitleLine[i] = allKnownSpectraList[i].peakList[0].sampleName;
            }
            for (int i = 0; i < allMergedList.Count; i++) {
                double[] line = new double[allMergedList.Count]; //
                for (int j = 0; j < allMergedList.Count; j++) { //
                    line[j] = Cal_Cosine_Product.CalCosineProduct(allMergedList[i], allMergedList[j]); //
                }
                NDPResultList.Add(line);
            }
            FileProcess.WriteNDPListToFile(NDPResultList, queryTitleLine, queryTitleLine, NDPmatrixPath); // output NDPmatrix

            watch.Stop();
            Console.WriteLine("weight: " + IdMSMS_Deconvolution.weight+" minSamp:"+IdMSMS_Deconvolution.minSampleNrForCor+" ndpThreh:"+IdMSMS_Deconvolution.ndpThreshold);
            Console.WriteLine("Deconvoluted spectra number: " + allMergedList.Count);
            Console.WriteLine("calculation time: " + watch.Elapsed);
            Console.WriteLine("OK");
            Console.Read();
        }
    }
}
