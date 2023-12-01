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

#include "window_functions.h"

#include <math.h>
#include <stdlib.h>
#include <assert.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#ifndef M_2PI
#define M_2PI (2.0 * M_PI)
#endif

double sincNormalized(double x) 
{
    if(x == 0.0) return 1.0;
    return sin(M_PI*x) / (M_PI*x);
}

void winSincFilter_64f(double *dst, int len, double fc, double offset)
{
    double n2 = ((double)len - 1.0) / 2.0;
    for(int i = 0; i < len; i++) {
        double x = i - n2 + offset;
        dst[i] = sincNormalized(2*fc*x);
    }
}

void winTriangle_64f(double *dst, int len)
{
    double d = double(len - 1) / 2.0;
    for(int i = 0; i < len; i++) {
        dst[i] *= 1.0 - fabs((i - d) / d);
    }
}

void winSine_64f(double *dst, int len)
{
    for(int i = 0; i < len; i++) {
        dst[i] *= sin(M_PI*i / (len-1));
    }
}

void winCosineSum_64f(double *dst, int len, int K, double *coeffs)
{
    for(int i = 0; i < len; i++) {
        double sum = 0.0;
        for(int k = 0; k < K; k++) {
            sum += pow(-1.0, k) * coeffs[k] * cos((M_2PI * k * i) / (len -1));
        }
        dst[i] *= sum;
    }
}

void winHanning_64f(double *dst, int len)
{
    double a[2] = { 0.5, 0.5 };
    winCosineSum_64f(dst, len, 2, a);
}

void winHamming_64f(double *dst, int len)
{
    double a[2] = { 0.54, 0.46 };
    winCosineSum_64f(dst, len, 2, a);
}

void winBlackman_64f(double *dst, int len)
{
    double a[3] = { 0.42, 0.5, 0.08 };
    winCosineSum_64f(dst, len, 3, a);
}

void winNuttall_64f(double *dst, int len)
{
    double a[4] = { 0.355768, 0.487396, 0.144232, 0.012604 };
    winCosineSum_64f(dst, len, 4, a);
}

void winFlattop_64f(double *dst, int len)
{
    double a[5] = { 1.0, 1.93, 1.29, 0.388, 0.028 };
    winCosineSum_64f(dst, len, 5, a);
}

void winGaussian_64f(double *dst, double sigma, int len)
{
    double a = double(len - 1.0) / 2.0;
    for(int i = 0; i < len; i++) {
        double b = (i - a) / (sigma * a);
        b *= b;
        dst[i] *= exp(-0.5 * b);
    }
}

static double confinedGaussianHelper(int x, double sigma, int len)
{
    double e = ((x - ((len-1)/2.0))/(2.0*sigma));
    return exp(-(e*e));
}

void winConfinedGaussian_64f(double *dst, double sigma, int len)
{
#define G(x) confinedGaussianHelper(x, sigma, len)
    sigma *= len;

    double c = G(-0.5);

    for(int i = 0; i < len; i++) {
        double num = G(i+len) + G(i-len);
        double denom = G(-0.5+len) + G(-0.5-len);

        dst[i] *= G(i) - (c * num / denom);
    }
#undef G
}

// Custom acosh()
static double m_acosh(double x)
{
    if(x < 1.0) {
        return (M_PI / 2.0) - atan(x / sqrt(1.0 - x*x));
    } else {
        return log(x + sqrt(x*x - 1.0));
    }
}

void winDolphChebyshev_64f(double *dst, double alpha, int len)
{
    int N = len;
    if(len & 1) {
        N = len - 1;
    } 

    double b = cosh((1.0/N) * m_acosh(pow(10.0, alpha)));

    for(int i = 0; i < N; i++) {
        double a = fabs(b * cos((M_PI * i) / N));
        double cheby = 0.0;

        if(a > 1) {
            cheby = cosh(N * m_acosh(a));
        } else {
            cheby = cos(N * acos(a));
        }

        // Alternating series
        dst[i] *= ((i & 1) ? - 1.0 : 1.0) * cheby;
    }

    // Inv DFT, only interested in real result
    // Running time, O(N^2)
    // Potential optimizations, eliminate cos() or replace with FFT
    double *w = (double*)malloc(len * sizeof(double)); 
    double sum = 0.0;
    for(int i = 0; i < N; i++) {
        sum = 0.0;
        for(int j = 0; j < N; j++) {
            sum += dst[j]*cos(M_2PI*i*j / N);
        }
        double t = sum / N;
        w[i] = sum / N;
    }
    for(int i = 0; i < N; i++) {
        dst[i] = w[i];
    }
    free(w);

    // Scale and make symettric if odd length
    dst[0] /= 2.0;
    if(len & 1) {
        dst[N] = dst[0];
    }
}

int factorial(int n)
{
    assert(n <= 12);
    if(n == 0) return 1.0;
    int f = n--;
    while(n > 0) f *= n--;
    return f;
}

// Modified bessel function
static double mbf(double x)
{
    double sum = 0.0;
    for(int k = 0; k < 13; k++) {
        double n = pow(x/2, k) / factorial(k);
        n *= n;
        sum += n;
    }

    return sum;
}

// Values generated between [0, N/2]
// For even lengths, copy and flip to lower portion of array
// For odd lengths, make symmetrical around 
void winKaiser_64f(double *dst, double alpha, int len)
{
    int offset = (len)/2;

    for(int i = 0; i < (len+1)/2; i++) {
        double a = (2.0*i) / (len - 1);
        dst[i + offset] *= mbf(M_PI*alpha*sqrt(1.0 - (a*a))) / mbf(M_PI*alpha);
    }

    for(int i = 0; i < (len)/2; i++) {
        dst[i] *= dst[len - i - 1];
    }
}

void winLanczos_64f(double *dst, int len)
{
    for(int i = 0; i < len; i++) {
        dst[i] *= sincNormalized(((2.0*i) / (len-1)) - 1);
    }
}

void winExponential_64f(double *dst, double t, int len, bool sym)
{
    if(sym && ((len & 1) == 0)) {
        assert(false);
    }

    int n2 = (len - 1) / 2;
    if(!sym) n2 = 0;
    double tinv = 1.0 / t;
    for(int i = 0; i < len; i++) {
        dst[i] *= exp(-abs(i - n2) * tinv);
    }
}

void winNormalize_64f(double *srcDst, int len)
{
    double sum = 0.0, inverseSum = 0.0;
    for(int i = 0; i < len; i++) {
        sum += srcDst[i];
    }
    sum /= len;
    inverseSum = 1.0 / sum;

    for(int i = 0; i < len; i++) {
        srcDst[i] *= inverseSum;
    }
}

void winNormalizeOne_64f(double *srcDst, int len)
{
    if(!srcDst || len < 1) return;

    double peak = srcDst[0];
    for(int i = 1; i < len; i++) {
        if(srcDst[i] > peak) {
            peak = srcDst[i];
        }
    }
    double scale = 1.0 / peak;
    for(int i = 0; i < len; i++) {
        srcDst[i] *= scale;
    }
}

double winMetricENBW(const double *win, int len)
{
    double sqrsum = 0.0, sum = 0.0;

    for(int i = 0; i < len; i++) {
        sum += win[i];
        sqrsum += win[i] * win[i];
    }

    return len * (sqrsum / (sum * sum));
}

// Single point DFT
// Frequency specified between 0.0 and 0.5
// Returns full scale dB value
static double spdft(const double *src, double f, int len)
{
    double re = 0.0, im = 0.0;
    for(int i = 0; i < len; i++) {
        re += src[i] * cos(M_2PI * f * i);
        im += src[i] * sin(M_2PI * f * i);
    }
    re /= len;
    im /= len;
    return 10.0 * log10(re*re + im*im);
}

double winMetricBandwidth(const double *win, int len, double n)
{
    double peak = spdft(win, 0.0, len);
    double meas = 0.0;
    double f = 0.5 / len;
    double step = 0.001 / len;
    do {
        f += step;
        meas = spdft(win, f, len);
    } while((meas + n) > peak);

    return 2.0 * f * len;
}

double winMetricScallopLoss(const double *win, int len)
{
    double peak = spdft(win, 0.0, len);
    double midBin = spdft(win, 0.5 / len, len);

    return peak - midBin;
}
