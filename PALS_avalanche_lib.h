#ifndef PALS_avalanche_lib_H
#define PALS_avalanche_lib_H

#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <functional>
#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <cmath>
#include <sys/stat.h>
#include <fstream>
#include <sstream>
#include <TString.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TFile.h>
#include <TCanvas.h>
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <TH1F.h>
#include <TLegend.h>
#include <TF1.h>
#include <TF2.h>
#include <TF12.h>
#include <TGraphErrors.h>
#include "TFile.h"
#define  ARMA_DONT_USE_WRAPPER
#define  ARMA_USE_LAPACK
#include <armadillo>
#include <time.h>

using namespace arma;

class LifetimeComponent
{
	public:
		float Lifetime;
		float Intensity;
		std::string Type;

		LifetimeComponent( float lifetime, float intensity, std::string type );	
		LifetimeComponent();
};

class DiscreteFitResult
{
	public:
		double Parameter;
		double Uncertainity;
		std::string Type;
		
		DiscreteFitResult();
		DiscreteFitResult( double Par, double Unc );
};


class Fit
{
	public:
		std::string Path;
		std::string PathWithDate;
		std::string TypeOfData;
		std::string TypeOfDataExtended;
		double BinWidth;
		unsigned LastBinMinValue;
		unsigned FirstBinMinValue;
		unsigned BackgroundEstimationNumbOfPoints;
		std::string SideOfBackgroundEstimation;
		unsigned ShiftForBackgroundEstimation;
		double EndOfFitValue;
		double StartOfFitValue;
		double FirstBinCenter;
		std::string FixGauss;
		std::vector<LifetimeComponent> Lifetimes;
		std::vector<LifetimeComponent> LifetimesNotFixed;
		std::vector<LifetimeComponent> Resolution;
		std::vector< double > Residue;

		std::vector< double > Times;
		std::vector< double > Arguments;
		std::vector< double > Values;
		unsigned BinMax;
		unsigned LastBin;
		unsigned FirstBin;
		unsigned Range_From;
		unsigned Range_To;
		double Background;
		double SDBackground;
		unsigned NmbOfBins;
		unsigned NmbrOfIterations;
		float VarLvl;
		std::string ShapeOfComponent;
		int TypeOfContinousFit;

		std::vector< DiscreteFitResult > LifetimesFromDiscrete;
		std::vector< DiscreteFitResult > IntensitiesFromDiscrete;		
		std::vector< DiscreteFitResult > ResolutionFromDiscrete;
		std::vector< DiscreteFitResult > OffsetsFromDiscrete;
		std::vector< DiscreteFitResult > FractionsFromDiscrete;
		double AreaFromDiscrete;

		int DecoOption;
		double Scaler;
		double FracMinLF;
		double FracMaxLF;
		vec LFGrid;
		mat GradMax;
		int MaxShiftForDeco;
		unsigned LinFilterRange;
		double EndOfFitMultiplicity;
		std::vector< double > LifetimeGrid;
		std::vector< double > Intensities_afterIter;
		std::vector< std::vector < double > > GradientMatrix;
		std::vector<double> SSig;
		unsigned iterator;
		double ChiDiffBetweenIterations;
		double FindMaxLF();
		double FindMinLF();
		double SigmasDefaultFraction;
		double AreaParameter;
		double StepForPreIteration;
		unsigned NumberOfPreIterations;
		double StepForIteration;
		unsigned NumberOfIterations;

		Fit();
		Fit( std::string path, std::string pathForDetails );
		Fit( std::string path, std::string nameOfHistogram, std::string pathForDetails );
		void LinearFilter( unsigned FilterRange );
		void RangeBackgroundData();
		int Discrete();
		int Discrete_old();
		int Discrete_exp();
		void SortLifetimesByType();
		void Deconvolution();
		double FindingOptimalAreaParameter( std::vector< double > Sig );
		double PreIterateGauss( std::vector< double > Sig, double Step, double Lambda, double PreviousChi );
		double PreIterateLogGauss( std::vector< double > Sig, double Step, double Lambda, double PreviousChi );
		double PreIterateMixed( std::vector< double > Sig, double Step, double Lambda, double PreviousChi );
		double Iterate( vec Intensities, double Step, double Lambda, double PreviousChi );
		double IterateExp( vec Intensities, double Step, double Lambda, double PreviousChi );
		double IterateExp2( vec Intensities, double Step, double Lambda, double PreviousChi );
};

Double_t SetIntensityParameter( std::vector<LifetimeComponent> Parameters, unsigned NumberOfParameter );
Double_t GetIntensityParameter( std::vector<Double_t> Parameters, unsigned NumberOfParameter );
Double_t GetIntensityParameterNew( Double_t *P, int type, unsigned startIndex, unsigned NumberOfParameter );
Double_t GetIntensityParameterErrorNew( Double_t *P, unsigned NumberOfParameter );

Double_t DiscreteFitFunctionNoPs( Double_t *A, Double_t *P );
Double_t DiscreteFitFunctionNoPs_old( Double_t *A, Double_t *P );
Double_t DiscreteFitFunctionNoPs_exp_1( Double_t *A, Double_t *P );
Double_t DiscreteFitFunctionNoPs_exp_2( Double_t *A, Double_t *P );
Double_t DiscreteFitFunctionNoPs_exp_3( Double_t *A, Double_t *P );
Double_t DiscreteFitFunctionPs( Double_t *A, Double_t *P );
Double_t DiscreteFitFunctionPs_old( Double_t *A, Double_t *P );
Double_t DiscreteFitFunctionPs_exp_1( Double_t *A, Double_t *P );
Double_t DiscreteFitFunctionPs_exp_2( Double_t *A, Double_t *P );
Double_t DiscreteFitFunctionPs_exp_3( Double_t *A, Double_t *P );

Double_t ExMoGa( Double_t A, Double_t P1, Double_t P2, Double_t P3, Double_t P4 );
Double_t ExMoGa_2( Double_t A, Double_t P1, Double_t P2, Double_t P3, Double_t P4 );
double GaussDistr( double X, double Mean, double Sigma );
double LogGaussDistr( double X, double Mean, double Sigma );
TH1F* FillHistogram( std::string Name, int BinNumber, double Start, double Stop, std::vector<double>& Data );
std::string NumberToChar( double Number, int Precision);
bool PathCheck( std::string Path );
bool FileCheck(const std::string& NameOfFile);
bool CheckTheLine( char FirstCharacter, char TestCharacter, unsigned NumberOfLine );

void TestIntensities();
const std::string GetTime();
double FindPeak( std::vector< double > Argument, std::vector< double > Values, unsigned Lowerfit, unsigned Higherfit );

#endif
