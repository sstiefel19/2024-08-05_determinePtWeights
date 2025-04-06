#include "TPython.h"


void RunPythonCode() {

    TPython::Exec(R"""(
from Ostap.Histos import histo
# Create a histogram
histo = ROOT.TH1F("histo", "Example Histogram", 100, 0, 10)
for i in range(1000):
    histo.Fill(ROOT.gRandom.Gaus(5, 1))

# Convert to TF1
f = histo.tf1(interpolation="linear")  # Specify interpolation type

# Evaluate the function at any x value
value = f.Eval(6.5)
print(f"Interpolated value at x=6.5: {value}");
)""");
}

void CheckPythonVersion() {
    TPython::Exec("import sys; print(sys.executable)");
    TPython::Exec("import sys; print(sys.version)");
}


void RunPythonFromROOT() {
    // TPython::Exec("print('Hello from Python!')");
    // RunPythonCode();
    CheckPythonVersion();
}

// void EvaluatePythonExpression() {
//     double result = TPython::Eval("3.14 * 2");
//     std::cout << "Result of Python expression: " << result << std::endl;
// }
