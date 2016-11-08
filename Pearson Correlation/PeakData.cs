using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Pearson_Correlation {
    class PeakData {
        public double mz { get; set; }
        public double rt { get; set; }
        public double intensity { get; set; }
        public string isotopes { get; set; }
        public string adduct { get; set; }
        public int pcgroup { get; set; }
        public string sampleName { get; set; }
        public double pval { get; set; }
        public double feature { get; set; }
        public bool isNoise { get; set; }
        public double diffrt { get; set; }
        public string precursor { get; set; }
    }
}
