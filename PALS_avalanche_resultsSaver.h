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

#include <PALS_avalanche_lib.h>
#include <PALS_avalanche_fileTools.h>
#include <PALS_avalanche_fitFunctions.h>

class ResultsSaver
{
	public:
		ResultsSaver();
		ResultsSaver( std::string TypeOfFittedModel );
		std::string TypeOfFit;
		
		FileTools fileTools;
		FitFunction fitFunction;
		void saveDiscreteFit( TH1F *histogram, TF1 *DiscreteFit, const char *RootFile, std::string PathOfFile, std::string ResultsPath, std::string PathOfFileWithDate );
		std::vector< std::vector< DiscreteFitResult > > saveDiscreteFitResultsToTXTandExcel( std::string Path, std::string Prefix, std::string PathWithDate, 
								    TF1* Discrete, double Background, double SDBackground, Double_t* ResolutionsFromFit, Double_t* ResolutionsFromFitErrors, 
								    int pPsIndex, double pPsIntensity, Double_t* FreeParameters, Double_t* FreeParametersErrors, 
								    double FixedIntensities, double FreeIntensitiesTopPs, double FixedFixedIntensity, std::vector<LifetimeComponent> Lifetimes, 
								    double oPsLifetimeCutoff, double pPsLifetimeCutoff, std::string PathWithDateForExcel );
		void saveResiduals( TH1F* histogram, double MinArgument, double MaxArgument, unsigned MinBin, unsigned MaxBin, TF1* Discrete, std::string FileName, std::string Path, std::string PathWithDate );
		void saveLFvsIntensities( std::string FileName, std::string Path, std::string PathWithDate, std::vector< DiscreteFitResult > LifetimesFromDiscrete, std::vector< DiscreteFitResult > IntensitiesFromDiscrete );
}

#endif
