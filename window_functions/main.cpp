#include "window_functions.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

void SaveToFile(const std::vector<double> &window, const char *fileName)
{
    std::ofstream file(fileName);
    if(!file.is_open()) {
        return;
    }

    for(int i = 0; i < window.size(); i++) {
        file << window[i] << std::endl;
    }
}

int main(int argc, char **argv)
{
    std::vector<double> window(33);
    window.assign(window.size(), 1);

    //winHamming_64f(window.data(), 0.54, 0.46, window.size());
    //winNormalize_64f(window.data(), window.size());
    winTriangle_64f(&window[0], window.size());

    //winKaiser_64f(window.data(), 3, window.size());
    //winNormalize_64f(window.data(), window.size());
    //double enbw = winMetricENBW(window.data(), window.size());
    //double bw3 = winMetricBandwidth(window.data(), window.size(), 3);
    //double bw6 = winMetricBandwidth(window.data(), window.size(), 6);

    window.assign(window.size(), 1);
    winSine_64f(&window[0], window.size());
    SaveToFile(window, "sine.csv");

    window.assign(window.size(), 1);
    winHanning_64f(&window[0], window.size());
    SaveToFile(window, "hanning.csv");

    window.assign(window.size(), 1);
    winNuttall_64f(&window[0], window.size());
    SaveToFile(window, "nuttall.csv");

    window.resize(1024);
    window.assign(window.size(), 1);
    //winFlattop_64f(&window[0], window.size());
    //winDolphChebyshev_64f(&window[0], 4.0, window.size());
    winKaiser_64f(&window[0], 2.0, window.size());

    double enbw = winMetricENBW(window.data(), window.size());
    double bw3dB = winMetricBandwidth(window.data(), window.size(), 3);
    double bw6dB = winMetricBandwidth(window.data(), window.size(), 6);
    double scallopLoss = winMetricScallopLoss(&window[0], window.size());

    return 0;
}
