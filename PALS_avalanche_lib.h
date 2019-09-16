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

#include "PALS_avalanche_libTools.h"
#include "PALS_avalanche_fileTools.h"
//#include "PALS_avalanche_fitFunctions.h"
#include "PALS_avalanche_resultsSaver.h"
#include "PALS_avalanche_deconvolution.h"

#define oPsLFLimit 0.7f
#define pPsLFLimit 0.2f

//using namespace arma;

class Fit
{
	public:
		FileTools fileTools;
        DCV DecoTools;
		//ResultsSaver Results;
		std::string TypeOfFit;
        
		std::string Path;
		std::string PathWithDate;
		std::string TypeOfData;
        int ROOTFileTestValue;
		std::string ROOTDirectory;
		std::string ROOTHistogram;
		std::string NameOfTheFileForEXCEL;
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
		double LastBinCenter;
		std::string FixGauss;
        int BackgroundVarLvl;
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

		unsigned NotFixedIterator;
		unsigned FixedIterator;
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
		Fit( std::string path, std::string pathForDetails, int ROOTFileTest, std::string FitType );
		int RenormalizeComponentsIntensities();
		void SortLifetimesByType();
		void GetHistogramToVector( std::string path, int ROOTFileTest );
		void LinearFilter( unsigned FilterRange );
		void RangeBackgroundData();
		int Discrete();
		void Deconvolution();
        
};

TH1F* FillHistogram( std::string Name, int BinNumber, double Start, double Stop, std::vector<double>& Data );
const std::string GetTime();
#endif
