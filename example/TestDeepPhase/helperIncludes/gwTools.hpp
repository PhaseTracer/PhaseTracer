#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <interpolation.h> //alglib

namespace gwTools {

inline double 
get_h() 
{
    return 0.678;
}

inline double 
get_H0() 
{
    return (100 * get_h()) / (3.086e19);
}

inline double 
convertSpectralToOmega(const double& S, double f) 
{
    double H0 = get_H0();
    double h = get_h();
    return 4 * M_PI*M_PI / (3 * H0*H0) * f*f*f * S * h*h;
}

inline double 
LISA_PSD(double f) 
{
    double Poms = 3.6e-41;
    double PaccNum = (2*M_PI*f) * (2*M_PI*f) * (2*M_PI*f) * (2*M_PI*f);
    double Pacc = 1.44e-48 / PaccNum * (1 + (0.4e-3/f)*(0.4e-3/f));

    double L = 2.5e9;
    double c = 2.99792458e8;
    double fs = c/(2*M_PI*L);

    double freqRatio = 3.*f/4.*fs;

    return 40./3. * (Poms + 4 * Pacc) * (1 + freqRatio*freqRatio);
}

inline double
LISA_omegahsq(double f) 
{
    double Sn = LISA_PSD(f);
    return convertSpectralToOmega(Sn, f);
}

inline double 
unresolvedGalacticBinary_PSD(double f)
{
    double fm = f * 1e3;
    double A = 8e-38;
    double f_ref = 1000;
    double a = 0.138;
    double b = -0.221;
    double c = 0.521;
    double d = 1.680;
    double fk = 1.13;

    double Sn;

    // above 0.01 Hz this is dominated by LISA noise
    if (fm > 10.0) 
    {
        Sn = 0;
    } else {
        Sn = A * pow(fm, 7./3.) * exp(- pow(fm/f_ref, a) - b*fm*sin(c*fm)) * (1 + tanh(d*(fk-fm)));
    }

    return Sn;
}

inline double
unresolvedGalacticBinary_omegahsq(double f)
{
    double Sn = unresolvedGalacticBinary_PSD(f);
    return convertSpectralToOmega(Sn, f);
}

inline double
extragalacticBinary_omegahsq(double f) 
{
    double amp = 8.9e-10;
    double f_ref = 25.0;
    double h = get_h();
    return amp * pow(f/f_ref, 2./3.) * h*h;
}

inline std::pair<double, double>
getBounds(const std::vector<double>& x)
{
    double xMin = x.front();
    double xMax = x.back();
    if(xMin > xMax) 
    {
        std::swap(xMin, xMax);
    }
    return {xMin, xMax};
}

inline double
trapezoidalIntegrate(const std::vector<double>& xVals, const std::vector<double>& yVals)
{
    if (xVals.size() != yVals.size()) 
    {
        throw std::invalid_argument("xVals and yVals must have the same size.");
    }
    if (xVals.size() < 2) 
    {
        throw std::invalid_argument("At least two points are required for integration.");
    }

    double integral = 0.0;
    for (size_t i = 0; i < xVals.size() - 1; ++i) 
    {
        double dx = xVals[i + 1] - xVals[i];
        double avgY = (yVals[i] + yVals[i + 1]) * 0.5;
        integral += dx * avgY;
    }
    return integral;
}

inline double
getSNR(const std::vector<double>& frequencyValues, const std::vector<double>& amplitudeValues, const double& Tyears = 4.0)
{
    std::vector<double> integrand;
    for(size_t i = 0; i < frequencyValues.size(); i++)
    {
        double fq = frequencyValues[i];
        double num = amplitudeValues[i];
        double noiseLISA = LISA_omegahsq(fq);
        double noiseUGB = unresolvedGalacticBinary_omegahsq(fq);
        double noiseEGB = extragalacticBinary_omegahsq(fq);
        double den = noiseLISA + noiseUGB + noiseEGB;
        double temp = num/den;
        integrand.push_back(temp*temp);
    }

    double Tobs = Tyears * 365.25*24*60*60;

    double result = trapezoidalIntegrate(frequencyValues, integrand);

    return sqrt(2 * Tobs * result);
}

inline alglib::spline1dinterpolant
createCubicSpline(const std::vector<double>& xVals, const std::vector<double>& yVals)
{
    alglib::real_1d_array xArray, yArray;
    xArray.setcontent(xVals.size(), xVals.data());
    yArray.setcontent(yVals.size(), yVals.data());

    alglib::spline1dinterpolant outputSpline;
    alglib::spline1dbuildcubic(xArray, yArray, outputSpline);

    return outputSpline;
}

inline std::vector<double>
convertToLog10(const std::vector<double>& X) 
{
    std::vector<double> logX(X.size());
    for(size_t i = 0; i < X.size(); i++) {logX[i] = std::log10(X[i]);}
    return logX;
}

inline double
noiseWeightedInnerProduct(
    const std::vector<double>& theoreticalFrequencyValues,
    const std::vector<double>& theoreticalAmplitudeValues,
    const std::vector<double>& fiducialFrequencyValues,
    const std::vector<double>& fiducialAmplitudeValues,
    const double& Tyear = 4.0,
    const int& N = 1000
)
{
    const std::vector<double> logTheoreticalFrequency = convertToLog10(theoreticalFrequencyValues);
    const std::vector<double> logTheoreticalAmplitude = convertToLog10(theoreticalAmplitudeValues);
    const std::vector<double> logFiducialFrequency = convertToLog10(fiducialFrequencyValues);
    const std::vector<double> logFiducialAmplitude = convertToLog10(fiducialAmplitudeValues);

    auto theoreticalSpline = createCubicSpline(logTheoreticalFrequency, logTheoreticalAmplitude);
    auto fiducialSpline = createCubicSpline(logFiducialFrequency, logFiducialAmplitude);

    auto fiducialBounds = getBounds(logFiducialFrequency);
    auto theoreticalBounds = getBounds(logTheoreticalFrequency);

    const double minLogFreq = std::max(fiducialBounds.first, theoreticalBounds.first);
    const double maxLogFreq = std::min(fiducialBounds.second, theoreticalBounds.second);
    
    const double logFreqInterval = (maxLogFreq - minLogFreq)/(N-1);

    std::vector<double> logFreqDomain, integrand;
    for (int i = 0; i < N; ++i)
    {
        double logFreq = minLogFreq + i * logFreqInterval;
        double freq = std::pow(10, logFreq);
        logFreqDomain.push_back(logFreq);

        double logTheoreticalAmplitude = alglib::spline1dcalc(theoreticalSpline, logFreq);
        double logFiducialAmplitude = alglib::spline1dcalc(fiducialSpline, logFreq);
        double logNoise = std::log10(LISA_omegahsq(freq) + unresolvedGalacticBinary_omegahsq(freq) + extragalacticBinary_omegahsq(freq));

        double fraction = logTheoreticalAmplitude + logFiducialAmplitude - 2 * logNoise;
        integrand.push_back( std::pow(10,fraction) * std::log(10) * freq);
    }

    const auto result = trapezoidalIntegrate(logFreqDomain, integrand);
    double Tobs = Tyear * 365.25*24*60*60;

    return 2 * Tobs * result;
}

inline double
getLogLikelihood(
    const std::vector<double>& theoreticalFrequencyValues,
    const std::vector<double>& theoreticalAmplitudeValues,
    const std::vector<double>& fiducialFrequencyValues,
    const std::vector<double>& fiducialAmplitudeValues,
    const double& Tyear = 4.0,
    const int& N = 1000,
    const bool includeFiducialSNR = false
)
{   
    double snrFiducial = 0.0;
    if(includeFiducialSNR)
    {
        snrFiducial = getSNR(fiducialFrequencyValues, fiducialAmplitudeValues, Tyear);
    }
    double snrTheory = getSNR(theoreticalFrequencyValues, theoreticalAmplitudeValues, Tyear);
    double innerProduct = noiseWeightedInnerProduct(theoreticalFrequencyValues, theoreticalAmplitudeValues, fiducialFrequencyValues, fiducialAmplitudeValues, Tyear, N);

    double logLikelihood = -0.5 * snrTheory*snrTheory - 0.5*snrFiducial*snrFiducial + innerProduct;
    return logLikelihood;
}

} // namespace gwTools