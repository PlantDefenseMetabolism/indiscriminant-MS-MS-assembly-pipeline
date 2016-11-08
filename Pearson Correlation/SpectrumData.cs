using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Pearson_Correlation {
    class SpectrumData {
        public List<PeakData> peakList { get; set; }
        public string group { get; set; }
        public List<double> neutralLossList { get; set; }
    }
}
