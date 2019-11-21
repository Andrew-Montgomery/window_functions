// Copyright (c) 2019 Andrew Montgomery

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef WINDOW_FUNCTIONS_H
#define WINDOW_FUNCTIONS_H

double sincNormalized(double x);

// Standard sinc filter kernel
// fc = frequency cutoff between [0.0, 0.5]
// offset = fractional sample offset, set to 0 for standard kernel
void winSincFilter_64f(double *dst, int len, double fc, double offset);

void winTriangle_64f(double *dst, int len);
void winSine_64f(double *dst, int len);

// Generalized cosine window
void winCosineSum_64f(double *dst, int len, int K, double *coeffs);

// Hann
// a0 = 0.5, a1 = 0.5
void winHanning_64f(double *dst, int len);
// a0 = 0.54, a1 = 0.46
void winHamming_64f(double *dst, int len);
// a0 = 0.42, a1 = 0.5, a2 = 0.08
void winBlackman_64f(double *dst, int len);
void winNuttall_64f(double *dst, int len);
void winFlattop_64f(double *dst, int len);
// sigma <= 0.5, smaller values equal larger bandwidths
void winGaussian_64f(double *dst, double sigma, int len);
// sigma = temporal width
void winConfinedGaussian_64f(double *dst, double sigma, int len);
// Sidelobe level = -20*alpha, runtime O(len^2)
void winDolphChebyshev_64f(double *dst, double alpha, int len);
// Kaiser-Bessel, typical alpha between 2 and 4
void winKaiser_64f(double *dst, double alpha, int len);

void winLanczos_64f(double *dst, int len);

// t = time constant for exponential decay where t is specified in # of samples
// exponential function decays approx 8.69dB per time constant
// sym = generate a symmetrical window, if true, len must be odd
// Non symmetrical windows is often called the response window
void winExponential_64f(double *dst, double t, int len, bool sym);

// Normalize so the sum of the window is equal to the length of the window
void winNormalize_64f(double *srcDst, int len);
// Scale so the max element is 1
void winNormalizeOne_64f(double *srcDst, int len);

// Get the noise bandwidth of a window function
double winMetricENBW(const double *win, int len);
// Get n-dB banwidth, commonly 3 or 6
double winMetricBandwidth(const double *win, int len, double n);
double winMetricScallopLoss(const double *win, int len);

#endif // WINDOW_FUNCTIONS_H
