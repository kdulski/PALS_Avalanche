#ifndef PALS_avalanche_deconvolution_H
#define PALS_avalanche_deconvolution_H

#include <vector>
#include <numeric>
#include <string>
#include <functional>
#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <cmath>
#define  ARMA_DONT_USE_WRAPPER
#define  ARMA_USE_LAPACK
#include <armadillo>

#include "PALS_avalanche_libTools.h"

using namespace arma;

class DCV
{
	public:
        DCV();
        double FindingOptimalAreaParameter( std::vector< double > Sig );
        double FindingOptimalSigmaParameter( std::vector< double > Sig, std::string ShapeOfComponent );
		double PreIterate( std::vector< double > Sig, double Step, double Lambda, double PreviousChi, std::string ShapeOfComponent );
		double Iterate( vec Intensities, double Step, double Lambda, double PreviousChi );
		double IterateExp( vec Intensities, double Step, double Lambda, double PreviousChi );
        
        void CopyParameter( std::vector< double > LifetimeGridOut, std::vector< DiscreteFitResult > LifetimesFromDiscreteOut, std::vector< DiscreteFitResult > IntensitiesFromDiscreteOut, std::vector< double > ValuesOut, mat GradMaxOut, double BackgroundOut, double SDBackgroundOut, double oPsLFLimitOut, unsigned NumberOfPreIterationsOut, unsigned NumberOfIterationsOut, double ChiDiffBetweenIterationsOut, std::vector< double > Intensities_afterIterOut );
        double AreaParameter;
        std::vector< double > SSig;
        unsigned iterator;
        double ChiDiffBetweenIterations;
        std::vector<double> Intensities_afterIter;
        std::vector< double > LifetimeGrid;
        std::vector< DiscreteFitResult > LifetimesFromDiscrete;
        std::vector< DiscreteFitResult > IntensitiesFromDiscrete;
        std::vector< double > Values;
        mat GradMax;
        double Background;
        double SDBackground;
        double oPsLFLimit;
        unsigned NumberOfPreIterations;
        unsigned NumberOfIterations;
        
        void CopyIteratorAndChiDiff( unsigned iteratorOut, double ChiDiffBetweenIterationsOut );
        unsigned GetInterator();
        std::vector< double > GetSigmas();
        std::vector< double > GetIntensities();
        void ClearIteratorAndChiSQ( int Mode );
};

double GaussDistr( double X, double Mean, double Sigma );
double LogGaussDistr( double X, double Mean, double Sigma );
double FindPeak( std::vector< double > Argument, std::vector< double > Values, unsigned Lowerfit, unsigned Higherfit );

#endif
