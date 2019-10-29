#ifndef PALS_avalanche_fitFunctions_H
#define PALS_avalanche_fitFunctions_H

#include <iostream>
#include <vector>
#include <numeric>
#include <string>
#include <functional>
#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <cmath>
#include <TROOT.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TF1.h>
#include <TF2.h>
#include <TF12.h>
#include <TGraphErrors.h>

#include "PALS_avalanche_libTools.h"
#include "PALS_avalanche_lib.h" 
#include "PALS_avalanche_fileTools.h"

//class LifetimeComponent;;

class FitFunction
{
	typedef Double_t ( FitFunction::*SelectedFunctionType )( Double_t *A, Double_t *P );
  
	public:
		FitFunction();
		FitFunction( std::string Approach, int pPsOption, int NumberOfResolutionComp, double HistoStart, double HistoEnd, unsigned NumberOfParameters, double NumberOfPointsForDrawing );
		FileTools fileTools;
		
		//double oPsLFLimit;
		
		//Double_t *SelectedFunction;
		SelectedFunctionType SelectedFunction;
		
		TH1F* histogramToFit;
		TF1 *Discrete;
		void generateFitFunction( double HistoStart, double HistoEnd, unsigned NumberOfParameters, double NumberOfPointsForDrawing );
		void generateParameter( unsigned NumberOfParameter, double InitValue, const char *NameOfParamter, double LowerLimit, double HigherLimit, std::string FixingOption );
		void generateResolutionParameter( unsigned NumberOfFirstParameter, std::vector<LifetimeComponent> Resolution, std::string FixGaussOption, double MaxPos );
		void generateInitLifetimeParameter( unsigned NumberOfFirstParameter, std::string TypeOfFit, std::vector<LifetimeComponent> Lifetimes, std::vector<LifetimeComponent> LifetimesNotFixed );
		void generateIterLifetimeParameter( unsigned NumberOfFirstParameter, std::string TypeOfFit, std::vector<LifetimeComponent> Lifetimes, std::vector<LifetimeComponent> LifetimesNotFixed, /*TF1 *DiscreteFitted, */unsigned Iteration, double VarLvl );
		TF1 * getFitFunction();
		void Fit( TH1F* histogram );
		
		void generateParameterOutside( unsigned NumberOfParameter, double InitValue, const char * NameOfParamter, double LowerLimit, double HigherLimit, std::string FixingOption, TF1* DiscreteFunction );
		void generateInitLifetimeParameterBRR( unsigned NumberOfFirstParameter, std::string TypeOfFit, std::vector<LifetimeComponent> Lifetimes, std::vector<LifetimeComponent> LifetimesNotFixed, TF1* DiscreteFunction );
		void generateIterLifetimeParameterBRR( unsigned NumberOfFirstParameter, std::string TypeOfFit, std::vector<LifetimeComponent> Lifetimes, std::vector<LifetimeComponent> LifetimesNotFixed, /*TF1 *DiscreteFitted, */unsigned Iteration, double VarLvl, TF1* DiscreteFunction  );
};
Double_t DiscreteFitFunctionNoPsSeparateComp( Double_t *A, Double_t *P );

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

Double_t DiscreteFitFunctionNoPs_forTests( Double_t *A, Double_t *P );

Double_t ExMoGa( Double_t A, Double_t P1, Double_t P2, Double_t P3, Double_t P4 );
Double_t ExMoGa_2( Double_t A, Double_t P1, Double_t P2, Double_t P3, Double_t P4 );

Double_t SetIntensityParameter( std::vector<LifetimeComponent> Parameters, unsigned NumberOfParameter );
Double_t GetIntensityParameter( std::vector<Double_t> Parameters, unsigned NumberOfParameter );
Double_t GetIntensityParameterNew( Double_t *P, int type, unsigned startIndex, unsigned NumberOfParameter );
Double_t GetIntensityParameterErrorNew( Double_t *P, unsigned NumberOfParameter );
#endif
