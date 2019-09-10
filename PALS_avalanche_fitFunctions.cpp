#include "PALS_avalanche_fitFunctions.h"

FitFunction::FitFunction( std::string Approach, int pPsOption, int NumberOfResolutionComp, double HistoStart, double HistoEnd, unsigned NumberOfParameters, double NumberOfPointsForDrawing )
{
	if( Approach == "" )
	{
		if( pPsOption == 0 ) // pPs option equal to zero means that there is pPs
		{
			//SelectedFunction = & FitFunction::DiscreteFitFunctionNoPs;
			Discrete = new TF1( "Discrete", DiscreteFitFunctionNoPs, HistoStart, HistoEnd, NumberOfParameters );
			Discrete -> SetNpx( NumberOfPointsForDrawing );
		}
		else
		{
			//SelectedFunction = & FitFunction::DiscreteFitFunctionPs;
			Discrete = new TF1( "Discrete", DiscreteFitFunctionPs, HistoStart, HistoEnd, NumberOfParameters );
			Discrete -> SetNpx( NumberOfPointsForDrawing );          
		}
	}
	else if( Approach == "exp" )
	{
		if( pPsOption == 0 ) // pPs option equal to zero means that there is no pPs
		{
			if( NumberOfResolutionComp == 1 )
			{
				//SelectedFunction = & FitFunction::DiscreteFitFunctionNoPs_exp_1;        
				Discrete = new TF1( "Discrete", DiscreteFitFunctionNoPs_exp_1, HistoStart, HistoEnd, NumberOfParameters );
				Discrete -> SetNpx( NumberOfPointsForDrawing );
			}
			else if( NumberOfResolutionComp == 2 )
			{
				//SelectedFunction = & FitFunction::DiscreteFitFunctionNoPs_exp_2;
				Discrete = new TF1( "Discrete", DiscreteFitFunctionNoPs_exp_2, HistoStart, HistoEnd, NumberOfParameters );
				Discrete -> SetNpx( NumberOfPointsForDrawing );
			}
			else if( NumberOfResolutionComp == 3 )
			{
				//SelectedFunction = & FitFunction::DiscreteFitFunctionNoPs_exp_3;  
				Discrete = new TF1( "Discrete", DiscreteFitFunctionNoPs_exp_3, HistoStart, HistoEnd, NumberOfParameters );
				Discrete -> SetNpx( NumberOfPointsForDrawing );
			}
			else 
			{
				SelectedFunction = NULL;
			}
		}
		else
		{
			if( NumberOfResolutionComp == 1 )
			{
				//SelectedFunction = & FitFunction::DiscreteFitFunctionPs_exp_1;    
				Discrete = new TF1( "Discrete", DiscreteFitFunctionPs_exp_1, HistoStart, HistoEnd, NumberOfParameters );
				Discrete -> SetNpx( NumberOfPointsForDrawing );
			}
			else if( NumberOfResolutionComp == 2 )
			{
				//SelectedFunction = & FitFunction::DiscreteFitFunctionPs_exp_2; 
				Discrete = new TF1( "Discrete", DiscreteFitFunctionPs_exp_2, HistoStart, HistoEnd, NumberOfParameters );
				Discrete -> SetNpx( NumberOfPointsForDrawing );
			}
			else if( NumberOfResolutionComp == 3 )
			{
				//SelectedFunction = & FitFunction::DiscreteFitFunctionPs_exp_3;
				Discrete = new TF1( "Discrete", DiscreteFitFunctionPs_exp_3, HistoStart, HistoEnd, NumberOfParameters );
				Discrete -> SetNpx( NumberOfPointsForDrawing );
			}
			else 
			{
				SelectedFunction = NULL;
			}
		}
	}
	else
	{
		if( pPsOption >= 0 ) // pPs option equal t zero means that there is pPs
		{
			//SelectedFunction = & FitFunction::DiscreteFitFunctionNoPs_old;
			Discrete = new TF1( "Discrete", DiscreteFitFunctionNoPs_old, HistoStart, HistoEnd, NumberOfParameters );
			Discrete -> SetNpx( NumberOfPointsForDrawing );
		}
		else
		{
			//SelectedFunction = & FitFunction::DiscreteFitFunctionPs_old;
			Discrete = new TF1( "Discrete", DiscreteFitFunctionPs_old, HistoStart, HistoEnd, NumberOfParameters );
			Discrete -> SetNpx( NumberOfPointsForDrawing );
		}        
	}
}

void FitFunction::generateParameter( unsigned NumberOfParameter, double InitValue, const char * NameOfParamter, double LowerLimit, double HigherLimit, std::string FixingOption )
{
	//std::cout << NumberOfParameter << " " << InitValue << " " << NameOfParamter << " " << LowerLimit << " " << HigherLimit << " " << FixingOption << std::endl;
 	Discrete -> SetParameter( NumberOfParameter, InitValue );
	Discrete -> SetParName( NumberOfParameter, NameOfParamter );
	if( FixingOption == "NoFixing" )
	    Discrete -> SetParLimits( NumberOfParameter, LowerLimit, HigherLimit );
	else if( FixingOption == "Fix" )
	    Discrete -> FixParameter( NumberOfParameter, InitValue );
}

void FitFunction::generateParameterOutside( unsigned NumberOfParameter, double InitValue, const char * NameOfParamter, double LowerLimit, double HigherLimit, std::string FixingOption, TF1* DiscreteFunction )
{
 	DiscreteFunction -> SetParameter( NumberOfParameter, InitValue );
	DiscreteFunction -> SetParName( NumberOfParameter, NameOfParamter );
	if( FixingOption == "NoFixing" )
	    DiscreteFunction -> SetParLimits( NumberOfParameter, LowerLimit, HigherLimit );
	else if( FixingOption == "Fix" )
	    DiscreteFunction -> FixParameter( NumberOfParameter, InitValue );
}

void FitFunction::generateResolutionParameter( unsigned NumberOfFirstParameter, std::vector<LifetimeComponent> Resolution, std::string FixGaussOption, double MaxPos )
{
	for( unsigned i=0; i<Resolution.size(); i++ )
	{
		if( FixGaussOption == "yes" )
		{
			generateParameter( NumberOfFirstParameter + 3*i, Resolution[i].Lifetime, ("Sigma for " + fileTools.NumberToChar( i+1, 0 ) + " Gauss").c_str(), 0., 0., "Fix" );
			generateParameter( NumberOfFirstParameter + 1 + 3*i, SetIntensityParameter( Resolution, i+1 ), ("Fraction parameter for " + fileTools.NumberToChar( i+1, 0 ) + " Gauss").c_str(), 0., 0., "Fix" );
			generateParameter( NumberOfFirstParameter + 2 + 3*i, MaxPos + i*0.33, ("Offset for " + fileTools.NumberToChar( i+1, 0 ) + " Gauss").c_str(), MaxPos - 0.5, MaxPos + i, "NoFixing" );
		}
		else
		{
			generateParameter( NumberOfFirstParameter + 3*i, Resolution[i].Lifetime, ("Sigma for " + fileTools.NumberToChar( i+1, 0 ) + " Gauss").c_str(), 0.01, 5, "NoFixing" );
			if( i == 0 )
				generateParameter( NumberOfFirstParameter + 1 + 3*i, SetIntensityParameter( Resolution, i+1 ), ("Fraction parameter for " + fileTools.NumberToChar( i+1, 0 ) + " Gauss").c_str(), 0., 1., "NoFixing" );
			else
				generateParameter( NumberOfFirstParameter + 1 + 3*i, SetIntensityParameter( Resolution, i+1 ), ("Fraction parameter for " + fileTools.NumberToChar( i+1, 0 ) + " Gauss").c_str(), 0., 3.14159/2, "NoFixing" );                
			generateParameter( NumberOfFirstParameter + 2 + 3*i, MaxPos + i*0.33, ("Offset for " + fileTools.NumberToChar( i+1, 0 ) + " Gauss").c_str(), MaxPos - 0.5, MaxPos + i + 1, "NoFixing" );
		}
	}
}

void FitFunction::generateInitLifetimeParameter( unsigned NumberOfFirstParameter, std::string TypeOfFit, std::vector<LifetimeComponent> Lifetimes, std::vector<LifetimeComponent> LifetimesNotFixed )
{
	unsigned NotFixedIterator = 0, FixedIterator = 0;
	for( unsigned j = 0; j < Lifetimes.size(); j++ )
	{
		if( Lifetimes[j].Type == "f" || Lifetimes[j].Type == "pf" || Lifetimes[j].Type == "lpf" || Lifetimes[j].Type == "lf" || Lifetimes[j].Type == "lff" )
		{
			generateParameter( NumberOfFirstParameter + 2*j, Lifetimes[j].Lifetime, ("Lifetime for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0., 0., "Fix" );
			generateParameter( NumberOfFirstParameter + 1 + 2*j, Lifetimes[j].Intensity, ("Intensity for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0., 0., "Fix" );
			FixedIterator++;
		}
		else if( Lifetimes[j].Type == "ps" )
		{
			generateParameter( NumberOfFirstParameter + 2*j, Lifetimes[j].Lifetime, ("Lifetime for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0., 0., "Fix" );
			generateParameter( NumberOfFirstParameter + 1 + 2*j, Lifetimes[j].Intensity, ("Intensity for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0., 0., "Fix" );
		}
		else
		{
			NotFixedIterator++;
			generateParameter( NumberOfFirstParameter + 2*j, Lifetimes[j].Lifetime, ("Lifetime for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0.08, 142., "Fix" );
			if( TypeOfFit == "" )
			{
				if( NotFixedIterator == 1 )
					generateParameter( NumberOfFirstParameter + 1 + 2*j, SetIntensityParameter( LifetimesNotFixed, NotFixedIterator ), ("Intensity for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0., 1., "NoFixing" );
				else
					generateParameter( NumberOfFirstParameter + 1 + 2*j, SetIntensityParameter( LifetimesNotFixed, NotFixedIterator ), ("Intensity for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0., 3.14159/2, "NoFixing" );
			}
			else
			{
				generateParameter( NumberOfFirstParameter + 1 + 2*j, Lifetimes[j].Intensity, ("Intensity for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0., 1., "NoFixing" );
			}
		}
	}
	generateParameter( NumberOfFirstParameter + 2 + 2*(Lifetimes.size()-1), FixedIterator, "Number of somehow fixed Components not ps type", 0., 0., "Fix" );
}

void FitFunction::generateIterLifetimeParameter( unsigned NumberOfFirstParameter, std::string TypeOfFit, std::vector<LifetimeComponent> Lifetimes, std::vector<LifetimeComponent> LifetimesNotFixed, /*TF1 *DiscreteFitted, */unsigned Iteration, double VarLvl )
{
	unsigned NotFixedIterator = 0, FixedIterator = 0;
	for( unsigned j = 0; j < Lifetimes.size(); j++ )
	{	
		if( Lifetimes[j].Type == "f" )
		{
			generateParameter( NumberOfFirstParameter + 2*j, Lifetimes[j].Lifetime, ("Lifetime for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0., 0., "Fix" );
			generateParameter( NumberOfFirstParameter + 1 + 2*j, Lifetimes[j].Intensity, ("Intensity for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0., 0., "Fix" );
			FixedIterator++;
		}
		else if( Lifetimes[j].Type == "pf" )
		{
			generateParameter( NumberOfFirstParameter + 2*j, Discrete -> GetParameter(NumberOfFirstParameter + 2*j), ("Lifetime for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 
					    Discrete -> GetParameter(NumberOfFirstParameter + 2*j) - VarLvl*Discrete -> GetParameter(NumberOfFirstParameter + 2*j), 
					    Discrete -> GetParameter(NumberOfFirstParameter + 2*j) + VarLvl*Discrete -> GetParameter(NumberOfFirstParameter + 2*j), "NoFixing" );
			generateParameter( NumberOfFirstParameter + 1 + 2*j, Discrete -> GetParameter(NumberOfFirstParameter + 1 + 2*j), ("Intensity for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0., 1, "NoFixing" );
			FixedIterator++;
		}
		else if( Lifetimes[j].Type == "lpf" )
		{
			generateParameter( NumberOfFirstParameter + 2*j, Discrete -> GetParameter(NumberOfFirstParameter + 2*j), ("Lifetime for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 
					    Discrete -> GetParameter(NumberOfFirstParameter + 2*j) - VarLvl*Discrete -> GetParameter(NumberOfFirstParameter + 2*j), 
					    Discrete -> GetParameter(NumberOfFirstParameter + 2*j) + VarLvl*Discrete -> GetParameter(NumberOfFirstParameter + 2*j), "NoFixing" );
			generateParameter( NumberOfFirstParameter + 1 + 2*j, Discrete -> GetParameter(NumberOfFirstParameter + 1 + 2*j), ("Intensity for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 
					    Discrete -> GetParameter(NumberOfFirstParameter + 1 + 2*j) - VarLvl*Discrete -> GetParameter(NumberOfFirstParameter + 1 + 2*j), 
					    Discrete -> GetParameter(NumberOfFirstParameter + 1 + 2*j) + VarLvl*Discrete -> GetParameter(NumberOfFirstParameter + 1 + 2*j), "NoFixing" );
			FixedIterator++;
		}
		else if( Lifetimes[j].Type == "lf" )
		{
			generateParameter( NumberOfFirstParameter + 2*j, Lifetimes[j].Lifetime, ("Lifetime for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0., 0., "Fix" );
			generateParameter( NumberOfFirstParameter + 1 + 2*j, Discrete -> GetParameter(NumberOfFirstParameter + 1 + 2*j), ("Intensity for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0., 1, "NoFixing" );
			FixedIterator++;
		}
		else if( Lifetimes[j].Type == "lff" )
		{
			generateParameter( NumberOfFirstParameter + 2*j, Lifetimes[j].Lifetime, ("Lifetime for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0., 0., "Fix" );
			generateParameter( NumberOfFirstParameter + 1 + 2*j, Discrete -> GetParameter(NumberOfFirstParameter + 1 + 2*j), ("Intensity for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 
					    Discrete -> GetParameter(NumberOfFirstParameter + 1 + 2*j) - VarLvl*Discrete -> GetParameter(NumberOfFirstParameter + 1 + 2*j), 
					    Discrete -> GetParameter(NumberOfFirstParameter + 1 + 2*j) + VarLvl*Discrete -> GetParameter(NumberOfFirstParameter + 1 + 2*j), "NoFixing" );
			FixedIterator++;
		}
		else if( Lifetimes[j].Type == "ps" )
		{
			generateParameter( NumberOfFirstParameter + 2*j, Discrete -> GetParameter(NumberOfFirstParameter + 2*j), ("Lifetime for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 
					    Discrete -> GetParameter(NumberOfFirstParameter + 2*j) - VarLvl*Discrete -> GetParameter(NumberOfFirstParameter + 2*j), 
					    Discrete -> GetParameter(NumberOfFirstParameter + 2*j) + VarLvl*Discrete -> GetParameter(NumberOfFirstParameter + 2*j), "NoFixing" );
			generateParameter( NumberOfFirstParameter + 1 + 2*j, Lifetimes[j].Intensity, ("Intensity for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0., 0., "Fix" );
		}
		else
		{
			NotFixedIterator++;	
			generateParameter( NumberOfFirstParameter + 2*j, Discrete -> GetParameter(NumberOfFirstParameter + 2*j), ("Lifetime for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0.08, 142., "NoFixing" );

			if( TypeOfFit == "" )
			{
				if( NotFixedIterator == 1 )
					generateParameter( NumberOfFirstParameter + 1 + 2*j, Discrete -> GetParameter(NumberOfFirstParameter + 1 + 2*j), ("Intensity for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0., 1., "NoFixing" );
				else
					generateParameter( NumberOfFirstParameter + 1 + 2*j, Discrete -> GetParameter(NumberOfFirstParameter + 1 + 2*j), ("Intensity for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0, 3.14159/2, "NoFixing" );
			}
			else
			{
				generateParameter( NumberOfFirstParameter + 1 + 2*j, Discrete -> GetParameter(NumberOfFirstParameter + 1 + 2*j), ("Intensity for " + fileTools.NumberToChar( j+1, 0 ) + " Component").c_str(), 0., 1., "NoFixing" );
			}
		}
	}
	generateParameter( NumberOfFirstParameter + 2 + 2*(Lifetimes.size()-1), FixedIterator, "Number of somehow fixed Components not ps type", 0., 0., "Fix" );
}

TF1 * FitFunction::getFitFunction()
{
	Discrete -> Update();
	return Discrete;
}

void FitFunction::Fit( TH1F* histogram )
{
	histogram -> Fit(Discrete,"RMW");
}

Double_t DiscreteFitFunctionNoPs( Double_t *A, Double_t *P )	//First parameter - nmbr of comp, Second - nmb of Gauss
{
	Double_t sum = 0., FixedIntensity = 0.;
	for( unsigned i = 0; i < P[1]; i++ )
	{
		for( unsigned j = 0; j < P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]; j++ )
		{
			sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
			FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
		}
		if( FixedIntensity >=1 )
		{
			return sum + P[0];
		}
		for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]; j < P[2]; j++ )
		{
			sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, 
						GetIntensityParameterNew( P, 0, 5 + 3*(int)P[1] + 2*P[6 + 3*(int)P[1] + 2*((int)P[2]-1)], (j-P[6 + 3*(int)P[1] + 2*((int)P[2]-1)])+1 )*(1 - FixedIntensity )/* Intensity */ );
		}	
		FixedIntensity = 0.;
	}
	return sum + P[0];
}

Double_t DiscreteFitFunctionNoPs_old( Double_t *A, Double_t *P )	//First parameter - nmbr of comp, Second - nmb of Gauss
{
	Double_t sum = 0., FixedIntensity = 0., FreeIntensity = 0.;	
	for( unsigned j = P[6 + 3*(int)P[1] + 2*(int)P[2]]; j < P[2]; j++ )
	{	  
		FreeIntensity += P[5 + 3*(int)P[1] + 2*j];
	}
	for( unsigned i = 0; i < P[1]; i++ )
	{
		for( unsigned j = 0; j < P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]; j++ )
		{
			sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
			FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
		}
		for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]; j < P[2]; j++ )
		{
			if( FixedIntensity < 1 )
				sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, 
																	  P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/FreeIntensity/* Intensity */ );
		}	
		FixedIntensity = 0.;
	}
	return sum + P[0];
}

Double_t DiscreteFitFunctionNoPs_exp_1( Double_t *A, Double_t *P )	//First parameter - nmbr of comp, Second - nmb of Gauss
{
	Double_t sum = 0., FixedIntensity = 0., FreeIntensity = 0.;	
	for( unsigned j = P[6 + 3*(int)P[1] + 2*(int)P[2]]; j < P[2]; j++ )
	{	  
		FreeIntensity += P[5 + 3*(int)P[1] + 2*j];
	}
	for( unsigned j = 0; j < P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]; j++ )
	{
		sum += P[3]/* Area */ * 1/* Intensity */ *ExMoGa( A[0], P[4]/* Sigma */, P[6]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
		FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
	}
	for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]; j < P[2]; j++ )
	{
		if( FixedIntensity < 1 )
			sum += P[3]/* Area */ * 1/* Intensity */ *ExMoGa( A[0], P[4]/* Sigma */, P[6]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/FreeIntensity/* Intensity */ );
	}	
	FixedIntensity = 0.;
	return sum + P[0];
}

Double_t DiscreteFitFunctionNoPs_exp_2( Double_t *A, Double_t *P )	//First parameter - nmbr of comp, Second - nmb of Gauss
{
	Double_t sum = 0., FixedIntensity = 0., FreeIntensity = 0.;	
	for( unsigned j = P[6 + 3*(int)P[1] + 2*(int)P[2]]; j < P[2]; j++ )
	{	  
		FreeIntensity += P[5 + 3*(int)P[1] + 2*j];
	}
	for( unsigned i = 0; i < P[1]; i++ )
	{
		for( unsigned j = 0; j < P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]; j++ )
		{
			if( i )
			{
				sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
				FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
			}
			else
			{
				sum += P[3]/* Area */ * P[5]/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
				FixedIntensity += P[5 + 3*(int)P[1] + 2*j];			  
			}
		}
		if( FixedIntensity < 1 )
		{
			for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]; j < P[2]; j++ )
			{
				if( i )
				{
					sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/FreeIntensity/* Intensity */ );
				}
				else
				{
					sum += P[3]/* Area */ * P[5]/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/FreeIntensity/* Intensity */ );			  
				}
			}
		}
		FixedIntensity = 0.;
	}
	return sum + P[0];
}

Double_t DiscreteFitFunctionNoPs_exp_3( Double_t *A, Double_t *P )	//First parameter - nmbr of comp, Second - nmb of Gauss
{
	Double_t sum = 0., FixedIntensity = 0., FreeIntensity = 0.;	
	for( unsigned j = P[6 + 3*(int)P[1] + 2*(int)P[2]]; j < P[2]; j++ )
	{	  
		FreeIntensity += P[5 + 3*(int)P[1] + 2*j];
	}
	for( unsigned i = 0; i < P[1]; i++ )
	{
		for( unsigned j = 0; j < P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]; j++ )
		{
			if( i == 1 )
			{
				sum += P[3]/* Area */ * P[5]*sin(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
				FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
			}
			else if( i == 2 )
			{
				sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
				FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
			}
			else
			{
				sum += P[3]/* Area */ * P[5]*cos(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
				FixedIntensity += P[5 + 3*(int)P[1] + 2*j];			  
			}
		}
		if( FixedIntensity < 1 )
		{
			for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]; j < P[2]; j++ )
			{
				if( i == 1 )
				{
					sum += P[3]/* Area */ * P[5]*sin(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, 
																		P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/FreeIntensity/* Intensity */ );
				}
				else if( i == 2 )
				{
					sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/FreeIntensity/* Intensity */ );
				}
				else
				{
					sum += P[3]/* Area */ * P[5]*cos(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, 
																		P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/FreeIntensity/* Intensity */ );			  
				}
			}
		}
		FixedIntensity = 0.;
	}
	return sum + P[0];
}

Double_t DiscreteFitFunctionPs( Double_t *A, Double_t *P )	//First parameter - nmbr of comp, Second - nmb of Gauss
{
	Double_t sum = 0., FixedIntensity = 0., FreepPsIntensity = 0.;	
	for( unsigned j = P[6 + 3*(int)P[1] + 2*(int)P[2]] + 4; j < P[2]; j++ )
	{
		if( P[4 + 3*(int)P[1] + 2*j] > 0.7 )
			FreepPsIntensity += P[5 + 3*(int)P[1]]*GetIntensityParameterNew( P, 2, 5 + 3*(int)P[1] + 2*P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 2, (j-P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]) );
	}
	for( unsigned i = 0; i < P[1]; i++ )
	{
		for( unsigned j = 1; j < P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j++ )
		{
			if( P[4 + 3*(int)P[1] + 2*j] > 0.7 )
			{
				sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *(ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + 
												      ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
				FixedIntensity += P[5 + 3*(int)P[1] + 2*j]*(1 + P[5 + 3*(int)P[1]]);
			}
			else
			{
				sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
				FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
			}
		}
		if( FixedIntensity < 1 )
		{
			for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j < P[2]; j++ )
			{
				sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, 
							GetIntensityParameterNew( P, 2, 5 + 3*(int)P[1] + 2*P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 2, (j-P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]) )*(1 - FixedIntensity )/(1 + FreepPsIntensity)/* Intensity */ );
			}
			sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, (1 - FixedIntensity )*FreepPsIntensity/(1 +
			FreepPsIntensity)/* Intensity */ );
		}
		FixedIntensity = 0.;
	}
	return sum + P[0];
}

Double_t DiscreteFitFunctionPs_old( Double_t *A, Double_t *P )	//First parameter - nmbr of comp, Second - nmb of Gauss
{
	Double_t sum = 0., FixedIntensity = 0., FreepPsIntensity = 0., FreeIntensity = 0.;	
	for( unsigned j = P[6 + 3*(int)P[1] + 2*(int)P[2]] + 4; j < P[2]; j++ )
	{
		if( P[4 + 3*(int)P[1] + 2*j] > 0.7 )
			FreepPsIntensity += P[5 + 3*(int)P[1]]*P[5 + 3*(int)P[1] + 2*j];
		FreeIntensity += P[5 + 3*(int)P[1] + 2*j];
	}
	for( unsigned i = 0; i < P[1]; i++ )
	{
		//std::cin >> FreeIntensity;
		for( unsigned j = 1; j < P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j++ )
		{
			if( P[4 + 3*(int)P[1] + 2*j] > 0.7 )
			{
				sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *(ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + 
												      ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
				FixedIntensity += P[5 + 3*(int)P[1] + 2*j]*(1 + P[5 + 3*(int)P[1]]);
			}
			else
			{
				sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
				FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
			}
		}
		if( FixedIntensity < 1 )
		{
			for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j < P[2]; j++ )
			{
				sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*
																				      (1 - FixedIntensity )/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
			}
			sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, (1 - FixedIntensity )*
																					FreepPsIntensity/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
		}
		FixedIntensity = 0.;
	}
	return sum + P[0];
}

Double_t DiscreteFitFunctionPs_exp_1( Double_t *A, Double_t *P )	//First parameter - nmbr of comp, Second - nmb of Gauss
{
	Double_t sum = 0., FixedIntensity = 0., FreepPsIntensity = 0., FreeIntensity = 0.;	
	for( unsigned j = P[6 + 3*(int)P[1] + 2*(int)P[2]] + 4; j < P[2]; j++ )
	{
		if( P[4 + 3*(int)P[1] + 2*j] > 0.7 )
			FreepPsIntensity += P[5 + 3*(int)P[1]]*P[5 + 3*(int)P[1] + 2*j];
		FreeIntensity += P[5 + 3*(int)P[1] + 2*j];
	}
	for( unsigned j = 1; j < P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j++ )
	{
		if( P[4 + 3*(int)P[1] + 2*j] > 0.7 )
		{
			sum += P[3]/* Area */ * 1/* Intensity */ *(ExMoGa( A[0], P[4]/* Sigma */, P[6]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + 
							  ExMoGa( A[0], P[4]/* Sigma */, P[6]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
			FixedIntensity += P[5 + 3*(int)P[1] + 2*j]*(1 + P[5 + 3*(int)P[1]]);
		}
		else
		{
			sum += P[3]/* Area */ * 1/* Intensity */ *ExMoGa( A[0], P[4]/* Sigma */, P[6]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
			FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
		}
	}
	if( FixedIntensity < 1 )
	{
		for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j < P[2]; j++ )
		{
			sum += P[3]/* Area */ * 1/* Intensity */ *ExMoGa( A[0], P[4]/* Sigma */, P[6]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
		}
	}
	return sum + P[0];
}

Double_t DiscreteFitFunctionPs_exp_2( Double_t *A, Double_t *P )	//First parameter - nmbr of comp, Second - nmb of Gauss
{
	Double_t sum = 0., FixedIntensity = 0., FreepPsIntensity = 0., FreeIntensity = 0.;	
	for( unsigned j = P[6 + 3*(int)P[1] + 2*(int)P[2]] + 4; j < P[2]; j++ )
	{
		if( P[4 + 3*(int)P[1] + 2*j] > 0.7 )
			FreepPsIntensity += P[5 + 3*(int)P[1]]*P[5 + 3*(int)P[1] + 2*j];
		FreeIntensity += P[5 + 3*(int)P[1] + 2*j];
	}
	for( unsigned i = 0; i < P[1]; i++ )
	{
		for( unsigned j = 1; j < P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j++ )
		{
			if( P[4 + 3*(int)P[1] + 2*j] > 0.7 )
			{
				if( i )
				{
					sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *(ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + 
										ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
					FixedIntensity += P[5 + 3*(int)P[1] + 2*j]*(1 + P[5 + 3*(int)P[1]]);				  
				}
				else
				{
					sum += P[3]/* Area */ * P[5]/* Intensity */ *(ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + 
										ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
					FixedIntensity += P[5 + 3*(int)P[1] + 2*j]*(1 + P[5 + 3*(int)P[1]]);
				}
			}
			else
			{
				if( i )
				{
					sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
					FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
				}
				else
				{
					sum += P[3]/* Area */ * P[5]/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
					FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
				}
			}
		}
		if( FixedIntensity < 1 )
		{
			for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j < P[2]; j++ )
			{
				if( i )
				{
					sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, 
											  P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );	  
				}
				else
				{
					sum += P[3]/* Area */ * P[5]/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, 
											  P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
				}
			}
		}
		if( i )
		{
			if( FixedIntensity < 1 )
				sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, (1 - FixedIntensity )*FreepPsIntensity/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
					FixedIntensity = 0.;		  
		}
		else
		{
			if( FixedIntensity < 1 )
				sum += P[3]/* Area */ * P[5]/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, (1 - FixedIntensity )*FreepPsIntensity/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
					FixedIntensity = 0.;		  
		}
	}
	return sum + P[0];
}

Double_t DiscreteFitFunctionPs_exp_3( Double_t *A, Double_t *P )	//First parameter - nmbr of comp, Second - nmb of Gauss
{
	Double_t sum = 0., FixedIntensity = 0., FreepPsIntensity = 0., FreeIntensity = 0.;	
	for( unsigned j = P[6 + 3*(int)P[1] + 2*(int)P[2]] + 4; j < P[2]; j++ )
	{
		if( P[4 + 3*(int)P[1] + 2*j] > 0.7 )
			FreepPsIntensity += P[5 + 3*(int)P[1]]*P[5 + 3*(int)P[1] + 2*j];
		FreeIntensity += P[5 + 3*(int)P[1] + 2*j];
	}
	for( unsigned i = 0; i < P[1]; i++ )
	{
		for( unsigned j = 1; j < P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j++ )
		{
			if( P[4 + 3*(int)P[1] + 2*j] > 0.7 )
			{
				if( i == 1 )
				{
					sum += P[3]/* Area */ * P[5]*sin(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *(ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + 
													    ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
					FixedIntensity += P[5 + 3*(int)P[1] + 2*j]*(1 + P[5 + 3*(int)P[1]]);				  
				}
				else if( i == 2 )
				{
					sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *(ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + 
										  ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
					FixedIntensity += P[5 + 3*(int)P[1] + 2*j]*(1 + P[5 + 3*(int)P[1]]);
				}
				else
				{
					sum += P[3]/* Area */ * P[5]*cos(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *(ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + 
													    ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
					FixedIntensity += P[5 + 3*(int)P[1] + 2*j]*(1 + P[5 + 3*(int)P[1]]);
				}
			}
			else
			{
				if( i == 1 )
				{
					sum += P[3]/* Area */ * P[5]*sin(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
					FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
				}
				else if( i == 2 )
				{
					sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
					FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
				}
				else
				{
					sum += P[3]/* Area */ * P[5]*cos(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
					FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
				}
			}
		}
		if( FixedIntensity < 1 )
		{
			for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j < P[2]; j++ )
			{
				switch( i )
				{
					case( 1 ):
						sum += P[3]/* Area */ * P[5]*sin(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, 
																P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
						break;
					case( 2 ):
						sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */,
												  P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
						break;
					case( 0 ):
						sum += P[3]/* Area */ * P[5]*cos(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, 
																P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
						break;
				}
			}
		}
		switch( i )
		{
			case( 1 ):
				if( FixedIntensity < 1 )
					sum += P[3]/* Area */ * P[5]*sin(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, 
						      P[4 + 3*(int)P[1]]/* Lifetime */, (1 - FixedIntensity )*FreepPsIntensity/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
				FixedIntensity = 0.;		  
				break;
			case( 2 ):
				if( FixedIntensity < 1 )
					sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, 
							P[4 + 3*(int)P[1]]/* Lifetime */, (1 - FixedIntensity )*FreepPsIntensity/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
				FixedIntensity = 0.;
				break;
			case( 0 ):
				if( FixedIntensity < 1 )
					sum += P[3]/* Area */ * P[5]*cos(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, 
							P[4 + 3*(int)P[1]]/* Lifetime */, (1 - FixedIntensity )*FreepPsIntensity/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
				FixedIntensity = 0.;
				break;
		}
	}
	return sum + P[0];
}

Double_t ExMoGa( Double_t A, Double_t P1, Double_t P2, Double_t P3, Double_t P4 )			//A - arguments, p - parameters
{
	Double_t Y = P4 / (2*P3) * exp( P1*P1 / (2*P3*P3) - ( A - P2 ) / P3 ) * ( erf( ( A - P1*P1 / P3 - P2 ) / (sqrt(2)*P1) ) - erf( ( - P1*P1 / P3 - P2 ) / (sqrt(2)*P1) ) );
	return Y;
}

Double_t ExMoGa_2( Double_t A, Double_t P1, Double_t P2, Double_t P3, Double_t P4 )			//A - arguments, p - parameters, alternative version of function
{
	Double_t Y = P4 / (2*P3) * exp( P1*P1 / (2*P3*P3) - ( A - P2 ) / P3 ) * ( 1 - erf( ( P2 + P1*P1 / P3 - A ) / (sqrt(2)*P1) ) );
	return Y;
}

double GaussDistr( double X, double Mean, double Sigma )
{
	return 1/sqrt(2*M_PI)/Sigma*exp(-pow(X-Mean, 2)/2/pow(Sigma, 2));
}

double LogGaussDistr( double X, double Mean, double Sigma )
{
	return 1/sqrt(2*M_PI)/Sigma/X*exp(-pow( log(X) - Mean, 2)/2/pow(Sigma, 2));
}

Double_t SetIntensityParameter( std::vector<LifetimeComponent> Parameters, unsigned NumberOfParameter ) // Inverse order from case number 4, because ordering is growing with complexity of formula
{
	double intensity = 0.1;
	unsigned dimension = Parameters.size();
	switch(dimension)
	{
		case 0:
			intensity = 1.;
			break;
		case 1:
			intensity = 1.;
			break;
		case 2:
			if( NumberOfParameter == 1 )
				intensity = Parameters[0].Intensity;
			else
				intensity = Parameters[1].Intensity;
			break;
		case 3:
			if( NumberOfParameter == 1 )
				intensity = 1 - Parameters[2].Intensity;		// Radius
			else if( NumberOfParameter == 2 )
				intensity = atan( Parameters[1].Intensity/Parameters[0].Intensity ); // Phi1 in Radians !!
			else
				intensity = Parameters[2].Intensity;
			break;
		case 4:
			if( NumberOfParameter == 1 )
				intensity = 1 - Parameters[3].Intensity;	// Radius
			else if( NumberOfParameter == 3 )
				intensity = atan( Parameters[2].Intensity/Parameters[1].Intensity ); // Phi1 in Radians !!
			else if( NumberOfParameter == 2 )
				intensity = atan( Parameters[1].Intensity*sqrt( pow( Parameters[2].Intensity/Parameters[1].Intensity, 2 ) + 1)/( Parameters[0].Intensity ) ); // Phi2 in Radians !!
			else
				intensity = Parameters[3].Intensity;
			break;
		case 5:
			if( NumberOfParameter == 1 )
				intensity = 1 - Parameters[4].Intensity;	// Radius
			else if( NumberOfParameter == 4 )
				intensity = atan( Parameters[3].Intensity/Parameters[2].Intensity ); // Phi3 in Radians !!
			else if( NumberOfParameter == 3 )
				intensity = atan( Parameters[2].Intensity/( Parameters[1].Intensity * cos( atan( Parameters[3].Intensity/Parameters[2].Intensity ) ) ) ); // Phi2 in Radians !!
			else if( NumberOfParameter == 2 )
				intensity = atan( Parameters[1].Intensity/( Parameters[0].Intensity * cos( atan( Parameters[2].Intensity/Parameters[1].Intensity ) ) ) ); // Phi1 in Radians !!
			else
				intensity = Parameters[4].Intensity;
			break;
		case 6:
			if( NumberOfParameter == 1 )
				intensity = 1 - Parameters[5].Intensity;	// Radius
			else if( NumberOfParameter == 5 )
				intensity = atan( Parameters[4].Intensity/Parameters[3].Intensity ); // Phi4 in Radians !!
			else if( NumberOfParameter == 4 )
				intensity = atan( Parameters[3].Intensity/( Parameters[2].Intensity * cos( atan( Parameters[4].Intensity/Parameters[3].Intensity ) ) ) ); // Phi3 in Radians !!
			else if( NumberOfParameter == 3 )
				intensity = atan( Parameters[2].Intensity/( Parameters[1].Intensity * cos( atan( Parameters[3].Intensity/Parameters[2].Intensity ) ) ) ); // Phi2 in Radians !!
			else if( NumberOfParameter == 2 )
				intensity = atan( Parameters[1].Intensity/( Parameters[0].Intensity * cos( atan( Parameters[2].Intensity/Parameters[1].Intensity ) ) ) ); // Phi1 in Radians !!
			else
				intensity = Parameters[5].Intensity;
			break;
		case 7:
			if( NumberOfParameter == 1 )
				intensity = 1 - Parameters[6].Intensity;	// Radius
			else if( NumberOfParameter == 6 )
				intensity = atan( Parameters[5].Intensity/Parameters[4].Intensity ); // Phi5 in Radians !!
			else if( NumberOfParameter == 5 )
				intensity = atan( Parameters[4].Intensity/( Parameters[3].Intensity * cos( atan( Parameters[5].Intensity/Parameters[4].Intensity ) ) ) ); // Phi4 in Radians !!
			else if( NumberOfParameter == 4 )
				intensity = atan( Parameters[3].Intensity/( Parameters[2].Intensity * cos( atan( Parameters[4].Intensity/Parameters[3].Intensity ) ) ) ); // Phi3 in Radians !!
			else if( NumberOfParameter == 3 )
				intensity = atan( Parameters[2].Intensity/( Parameters[1].Intensity * cos( atan( Parameters[3].Intensity/Parameters[2].Intensity ) ) ) ); // Phi2 in Radians !!
			else if( NumberOfParameter == 2 )
				intensity = atan( Parameters[1].Intensity/( Parameters[0].Intensity * cos( atan( Parameters[2].Intensity/Parameters[1].Intensity ) ) ) ); // Phi1 in Radians !!
			else
				intensity = Parameters[6].Intensity;
			break;
		case 8:
			if( NumberOfParameter == 1 )
				intensity = 1 - Parameters[7].Intensity;	// Radius
			else if( NumberOfParameter == 7 )
				intensity = atan( Parameters[6].Intensity/Parameters[5].Intensity ); // Phi1 in Radians !!
			else if( NumberOfParameter == 6 )
				intensity = atan( Parameters[5].Intensity/( Parameters[4].Intensity * cos( atan( Parameters[6].Intensity/Parameters[5].Intensity ) ) ) ); // Phi2 in Radians !!
			else if( NumberOfParameter == 5 )
				intensity = atan( Parameters[4].Intensity/( Parameters[3].Intensity * cos( atan( Parameters[5].Intensity/Parameters[4].Intensity ) ) ) ); // Phi3 in Radians !!
			else if( NumberOfParameter == 4 )
				intensity = atan( Parameters[3].Intensity/( Parameters[2].Intensity * cos( atan( Parameters[4].Intensity/Parameters[3].Intensity ) ) ) ); // Phi4 in Radians !!
			else if( NumberOfParameter == 3 )
				intensity = atan( Parameters[2].Intensity/( Parameters[1].Intensity * cos( atan( Parameters[3].Intensity/Parameters[2].Intensity ) ) ) ); // Phi5 in Radians !!
			else if( NumberOfParameter == 2 )
				intensity = atan( Parameters[1].Intensity/( Parameters[0].Intensity * cos( atan( Parameters[2].Intensity/Parameters[1].Intensity ) ) ) ); // Phi6 in Radians !!
			else
				intensity = Parameters[7].Intensity;
			break;
		case 9:
			if( NumberOfParameter == 1 )
				intensity = 1 - Parameters[8].Intensity;	// Radius
			else if( NumberOfParameter == 8 )
				intensity = atan( Parameters[7].Intensity/Parameters[6].Intensity ); // Phi1 in Radians !!
			else if( NumberOfParameter == 7 )
				intensity = atan( Parameters[6].Intensity/( Parameters[5].Intensity * cos( atan( Parameters[7].Intensity/Parameters[6].Intensity ) ) ) ); // Phi2 in Radians !!
			else if( NumberOfParameter == 6 )
				intensity = atan( Parameters[5].Intensity/( Parameters[4].Intensity * cos( atan( Parameters[6].Intensity/Parameters[5].Intensity ) ) ) ); // Phi3 in Radians !!
			else if( NumberOfParameter == 5 )
				intensity = atan( Parameters[4].Intensity/( Parameters[3].Intensity * cos( atan( Parameters[5].Intensity/Parameters[4].Intensity ) ) ) ); // Phi4 in Radians !!
			else if( NumberOfParameter == 4 )
				intensity = atan( Parameters[3].Intensity/( Parameters[2].Intensity * cos( atan( Parameters[4].Intensity/Parameters[3].Intensity ) ) ) ); // Phi5 in Radians !!
			else if( NumberOfParameter == 3 )
				intensity = atan( Parameters[2].Intensity/( Parameters[1].Intensity * cos( atan( Parameters[3].Intensity/Parameters[2].Intensity ) ) ) ); // Phi6 in Radians !!
			else if( NumberOfParameter == 2 )
				intensity = atan( Parameters[1].Intensity/( Parameters[0].Intensity * cos( atan( Parameters[2].Intensity/Parameters[1].Intensity ) ) ) ); // Phi7 in Radians !!
			else
				intensity = Parameters[8].Intensity;
			break;
		case 10:
			if( NumberOfParameter == 1 )
				intensity = 1 - Parameters[9].Intensity;	// Radius
			else if( NumberOfParameter == 9 )
				intensity = atan( Parameters[8].Intensity/Parameters[7].Intensity ); // Phi1 in Radians !!
			else if( NumberOfParameter == 8 )
				intensity = atan( Parameters[7].Intensity/( Parameters[6].Intensity * cos( atan( Parameters[8].Intensity/Parameters[7].Intensity ) ) ) ); // Phi2 in Radians !!
			else if( NumberOfParameter == 7 )
				intensity = atan( Parameters[6].Intensity/( Parameters[5].Intensity * cos( atan( Parameters[7].Intensity/Parameters[6].Intensity ) ) ) ); // Phi3 in Radians !!
			else if( NumberOfParameter == 6 )
				intensity = atan( Parameters[5].Intensity/( Parameters[4].Intensity * cos( atan( Parameters[6].Intensity/Parameters[5].Intensity ) ) ) ); // Phi4 in Radians !!
			else if( NumberOfParameter == 5 )
				intensity = atan( Parameters[4].Intensity/( Parameters[3].Intensity * cos( atan( Parameters[5].Intensity/Parameters[4].Intensity ) ) ) ); // Phi5 in Radians !!
			else if( NumberOfParameter == 4 )
				intensity = atan( Parameters[3].Intensity/( Parameters[2].Intensity * cos( atan( Parameters[4].Intensity/Parameters[3].Intensity ) ) ) ); // Phi6 in Radians !!
			else if( NumberOfParameter == 3 )
				intensity = atan( Parameters[2].Intensity/( Parameters[1].Intensity * cos( atan( Parameters[3].Intensity/Parameters[2].Intensity ) ) ) ); // Phi7 in Radians !!
			else if( NumberOfParameter == 2 )
				intensity = atan( Parameters[1].Intensity/( Parameters[0].Intensity * cos( atan( Parameters[2].Intensity/Parameters[1].Intensity ) ) ) ); // Phi8 in Radians !!
			else
				intensity = Parameters[8].Intensity;
			break;
	}
	return (Double_t)intensity;
}

Double_t GetIntensityParameter( std::vector<Double_t> Parameters, unsigned NumberOfParameter )
{
	double intensity = 0.1;
	unsigned dimension = Parameters.size();
	double Amplitude = 1.; //
	for( unsigned i=dimension-2; i>0; i-- )
	{
		Amplitude = Amplitude*sin( Parameters[i] ) + cos( Parameters[i] );
	}
	switch(dimension)
	{
		case 0:
			intensity = 1.;
			break;
		case 1:
			intensity = 1.;
			break;
		case 2:
			if( NumberOfParameter == 2 )
				intensity = 1 - Parameters[0];
			else
				intensity = Parameters[0];
			break;
		case 3:
			if( NumberOfParameter == 1 )
				intensity = Parameters[0]*cos(Parameters[1]) /Amplitude;		// r cos(Phi) / Amplitude
			else if( NumberOfParameter == 2 )
				intensity = Parameters[0]*sin(Parameters[1]) /Amplitude; // r sin(Phi) / Amplitude
			else
				intensity = 1 - Parameters[0];				// 1 - r
			break;
		case 4:
			if( NumberOfParameter == 1 )
				intensity = Parameters[0]*cos(Parameters[1]) /Amplitude;  // r cos(Phi1) / Amplitude
			else if( NumberOfParameter == 2 )
				intensity = Parameters[0]*sin(Parameters[1])*cos(Parameters[2]) /Amplitude; // r sin(Phi1)*sin(Phi2) / Amplitude
			else if( NumberOfParameter == 3 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2]) /Amplitude;
			else
				intensity = 1 - Parameters[0];
			break;
		case 5:
			if( NumberOfParameter == 1 )
                                intensity = Parameters[0]*cos(Parameters[1]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = Parameters[0]*sin(Parameters[1])*cos(Parameters[2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*cos(Parameters[3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3]) /Amplitude;
			else
				intensity = 1 - Parameters[0];
			break;
		case 6:
			if( NumberOfParameter == 1 )
                                intensity = Parameters[0]*cos(Parameters[1]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = Parameters[0]*sin(Parameters[1])*cos(Parameters[2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*cos(Parameters[3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*cos(Parameters[4]) /Amplitude;
			else if( NumberOfParameter == 5 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*sin(Parameters[4]) /Amplitude;
			else
				intensity = 1 - Parameters[0];
			break;
		case 7:
			if( NumberOfParameter == 1 )
                                intensity = Parameters[0]*cos(Parameters[1]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = Parameters[0]*sin(Parameters[1])*cos(Parameters[2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*cos(Parameters[3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*cos(Parameters[4]) /Amplitude;
			else if( NumberOfParameter == 5 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*sin(Parameters[4])*cos(Parameters[5]) /Amplitude;
			else if( NumberOfParameter == 6 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*sin(Parameters[4])*sin(Parameters[5]) /Amplitude;
			else
				intensity = 1 - Parameters[0];
			break;
		case 8:
			if( NumberOfParameter == 1 )
                                intensity = Parameters[0]*cos(Parameters[1]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = Parameters[0]*sin(Parameters[1])*cos(Parameters[2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*cos(Parameters[3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*cos(Parameters[4]) /Amplitude;
			else if( NumberOfParameter == 5 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*sin(Parameters[4])*cos(Parameters[5]) /Amplitude;
			else if( NumberOfParameter == 6 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*sin(Parameters[4])*sin(Parameters[5])*cos(Parameters[6]) /Amplitude;
			else if( NumberOfParameter == 7 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*sin(Parameters[4])*sin(Parameters[5])*sin(Parameters[6]) /Amplitude;
			else
				intensity = 1 - Parameters[0];
			break;
		case 9:
			if( NumberOfParameter == 1 )
                                intensity = Parameters[0]*cos(Parameters[1]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = Parameters[0]*sin(Parameters[1])*cos(Parameters[2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*cos(Parameters[3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*cos(Parameters[4]) /Amplitude;
			else if( NumberOfParameter == 5 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*sin(Parameters[4])*cos(Parameters[5]) /Amplitude;
			else if( NumberOfParameter == 6 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*sin(Parameters[4])*sin(Parameters[5])*cos(Parameters[6]) /Amplitude;
			else if( NumberOfParameter == 7 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*sin(Parameters[4])*sin(Parameters[5])*sin(Parameters[6])*cos(Parameters[7]) /Amplitude;
			else if( NumberOfParameter == 8 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*sin(Parameters[4])*sin(Parameters[5])*sin(Parameters[6])*sin(Parameters[7]) /Amplitude;
			else
				intensity = 1 - Parameters[0];
			break;
		case 10:
			if( NumberOfParameter == 1 )
                                intensity = Parameters[0]*cos(Parameters[1]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = Parameters[0]*sin(Parameters[1])*cos(Parameters[2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*cos(Parameters[3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*cos(Parameters[4]) /Amplitude;
			else if( NumberOfParameter == 5 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*sin(Parameters[4])*cos(Parameters[5]) /Amplitude;
			else if( NumberOfParameter == 6 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*sin(Parameters[4])*sin(Parameters[5])*cos(Parameters[6]) /Amplitude;
			else if( NumberOfParameter == 7 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*sin(Parameters[4])*sin(Parameters[5])*sin(Parameters[6])*cos(Parameters[7]) /Amplitude;
			else if( NumberOfParameter == 8 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*sin(Parameters[4])*sin(Parameters[5])*sin(Parameters[6])*sin(Parameters[7])*cos(Parameters[8]) /Amplitude;
			else if( NumberOfParameter == 9 )
				intensity = Parameters[0]*sin(Parameters[1])*sin(Parameters[2])*sin(Parameters[3])*sin(Parameters[4])*sin(Parameters[5])*sin(Parameters[6])*sin(Parameters[7])*sin(Parameters[8]) /Amplitude;
			else
				intensity = 1 - Parameters[0];
			break;
	}
	return (Double_t)intensity;
}

Double_t GetIntensityParameterNew( Double_t *P, int Type, unsigned startIndex, unsigned NumberOfParameter )
{
  	unsigned dimension = 0; 
	int type = Type;
	switch(type)
	{
		case 0:
			dimension = P[2] - P[6 + 3*(int)P[1] + 2*((int)P[2]-1)];			//Free Lifetime components parameters without p-Ps
			break;
		case 1:
			dimension = P[1];								//Resolution components (type is one because between different components there are three paramateres difference) And default difference between parameters in lifetimes is two -> two plus 1 is three (For sake not having negative parameters and to unify function)
			break;	
		case 2:
			dimension = P[2] - P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] - 1;			//Free Lifetime components parameters with p-Ps
			type = 0;
			break;
		case 3:											
			dimension = P[1];
			type = 0;
			break;
	}
	double intensity = 0.1;
	double Amplitude = 1.; //
	for( int i=dimension-2; i>0; i-- )
	{
		Amplitude = Amplitude*sin( P[startIndex + (2 + type)*i] ) + cos( P[startIndex + (2 + type)*i] );
	}
	switch(dimension)
	{
		case 0:
			intensity = 1.;
			break;
		case 1:
			intensity = 1.;
			break;
		case 2:
			if( NumberOfParameter == 1 )
				intensity = P[startIndex];
			else
				intensity = 1 - P[startIndex];
			break;
		case 3:
			if( NumberOfParameter == 1 )
				intensity = P[startIndex]*cos(P[startIndex + (2 + type)]) /Amplitude;		// r cos(Phi) / Amplitude
			else if( NumberOfParameter == 2 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)]) /Amplitude; // r sin(Phi) / Amplitude
			else
				intensity = 1 - P[startIndex];				// 1 - r
			break;
		case 4:
			if( NumberOfParameter == 1 )
				intensity = P[startIndex]*cos(P[startIndex + (2 + type)]) /Amplitude;  // r cos(Phi1) / Amplitude
			else if( NumberOfParameter == 2 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*cos(P[startIndex + (2 + type)*2]) /Amplitude; // r sin(Phi1)*sin(Phi2) / Amplitude
			else if( NumberOfParameter == 3 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2]) /Amplitude;
			else
				intensity = 1 - P[startIndex];
			break;
		case 5:
			if( NumberOfParameter == 1 )
                                intensity = P[startIndex]*cos(P[startIndex + (2 + type)]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*cos(P[startIndex + (2 + type)*2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*cos(P[startIndex + (2 + type)*3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3]) /Amplitude;
			else
				intensity = 1 - P[startIndex];
			break;
		case 6:
			if( NumberOfParameter == 1 )
                                intensity = P[startIndex]*cos(P[startIndex + (2 + type)]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*cos(P[startIndex + (2 + type)*2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*cos(P[startIndex + (2 + type)*3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*cos(P[startIndex + (2 + type)*4]) /Amplitude;
			else if( NumberOfParameter == 5 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4]) /Amplitude;
			else
				intensity = 1 - P[startIndex + (2 + type)];
			break;
		case 7:
			if( NumberOfParameter == 1 )
                                intensity = P[startIndex]*cos(P[startIndex + (2 + type)]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*cos(P[startIndex + (2 + type)*2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*cos(P[startIndex + (2 + type)*3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*cos(P[startIndex + (2 + type)*4]) /Amplitude;
			else if( NumberOfParameter == 5 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*cos(P[startIndex + (2 + type)*5]) /Amplitude;
			else if( NumberOfParameter == 6 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5]) /Amplitude;
			else
				intensity = 1 - P[startIndex];
			break;
		case 8:
			if( NumberOfParameter == 1 )
                                intensity = P[startIndex]*cos(P[startIndex + (2 + type)]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*cos(P[startIndex + (2 + type)*2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*cos(P[startIndex + (2 + type)*3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*cos(P[startIndex + (2 + type)*4]) /Amplitude;
			else if( NumberOfParameter == 5 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*cos(P[startIndex + (2 + type)*5]) /Amplitude;
			else if( NumberOfParameter == 6 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*cos(P[startIndex + (2 + type)*6]) /Amplitude;
			else if( NumberOfParameter == 7 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*sin(P[startIndex + (2 + type)*6]) /Amplitude;
			else
				intensity = 1 - P[startIndex];
			break;
		case 9:
			if( NumberOfParameter == 1 )
                                intensity = P[startIndex]*cos(P[startIndex + (2 + type)]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*cos(P[startIndex + (2 + type)*2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*cos(P[startIndex + (2 + type)*3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*cos(P[startIndex + (2 + type)*4]) /Amplitude;
			else if( NumberOfParameter == 5 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*cos(P[startIndex + (2 + type)*5]) /Amplitude;
			else if( NumberOfParameter == 6 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*cos(P[startIndex + (2 + type)*6]) /Amplitude;
			else if( NumberOfParameter == 7 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*sin(P[startIndex + (2 + type)*6])*cos(P[startIndex + (2 + type)*7]) /Amplitude;
			else if( NumberOfParameter == 8 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*sin(P[startIndex + (2 + type)*6])*sin(P[startIndex + (2 + type)*7]) /Amplitude;
			else
				intensity = 1 - P[startIndex];
			break;
		case 10:
			if( NumberOfParameter == 1 )
                                intensity = P[startIndex]*cos(P[startIndex + (2 + type)]) /Amplitude;  // r cos(Phi1) / Amplitude
                        else if( NumberOfParameter == 2 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*cos(P[startIndex + (2 + type)*2]) /Amplitude; // r sin(Phi1)*cos(Phi2) / Amplitude and so on ...
                        else if( NumberOfParameter == 3 )
                                intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*cos(P[startIndex + (2 + type)*3]) /Amplitude;
			else if( NumberOfParameter == 4 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*cos(P[startIndex + (2 + type)*4]) /Amplitude;
			else if( NumberOfParameter == 5 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*cos(P[startIndex + (2 + type)*5]) /Amplitude;
			else if( NumberOfParameter == 6 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*cos(P[startIndex + (2 + type)*6]) /Amplitude;
			else if( NumberOfParameter == 7 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*sin(P[startIndex + (2 + type)*6])*cos(P[startIndex + (2 + type)*7]) /Amplitude;
			else if( NumberOfParameter == 8 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*sin(P[startIndex + (2 + type)*6])*sin(P[startIndex + (2 + type)*7])*cos(P[startIndex + (2 + type)*8]) /Amplitude;
			else if( NumberOfParameter == 9 )
				intensity = P[startIndex]*sin(P[startIndex + (2 + type)])*sin(P[startIndex + (2 + type)*2])*sin(P[startIndex + (2 + type)*3])*sin(P[startIndex + (2 + type)*4])*sin(P[startIndex + (2 + type)*5])*sin(P[startIndex + (2 + type)*6])*sin(P[startIndex + (2 + type)*7])*sin(P[startIndex + (2 + type)*8]) /Amplitude;
			else
				intensity = 1 - P[startIndex];
			break;
	}
	return (Double_t)intensity;
}

Double_t GetIntensityParameterErrorNew( Double_t *P, unsigned NumberOfParameter )
{
  	unsigned dimension = P[0];
	double intensity = 0.1;
	double Amplitude = 1.; //
	for( int i=dimension-2; i>0; i-- )
	{
		Amplitude = Amplitude*sin( P[2*i+1] ) + cos( P[2*i+1] );
	}
	switch(dimension)
	{
		case 0:
			intensity = 1.;
			break;
		case 1:
			intensity = P[2];
			break;
		case 2:
			intensity = P[2];
			break;
		case 3:
			if( NumberOfParameter == 1 )
			{
				intensity = sqrt( pow(cos(P[3])*P[2]/Amplitude, 2) + pow( P[1]*P[4]/pow(Amplitude, 2), 2) );
			}
			else if( NumberOfParameter == 2 )
				intensity = sqrt( pow(sin(P[3])*P[2]/Amplitude, 2) + pow( P[1]*P[4]/pow(Amplitude, 2), 2) );
			else
				intensity = P[2];				// 
			break;
		case 4:
			if( NumberOfParameter == 1 )
				intensity = sqrt( pow(cos(P[3])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*(sin(P[5]) + cos(P[5]))/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*cos(P[3])*(sin(P[5]) - cos(P[5]))/pow(Amplitude, 2), 2) );
			else if( NumberOfParameter == 2 )
				intensity = sqrt( pow(sin(P[3])*cos(P[5])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*cos(P[5])/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*( sin(P[3]) * pow( sin(P[5]),2 ) + cos(P[3]) * pow( cos(P[5]),2 ) + sin(P[5])*cos(P[3]) )/pow(Amplitude, 2), 2) );
			else if( NumberOfParameter == 3 )
				intensity = sqrt( pow(sin(P[3])*sin(P[5])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*sin(P[5])/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*( sin(P[3]) * pow( sin(P[5]),2 ) + sin(P[3]) * pow( cos(P[5]),2 ) + cos(P[5])*cos(P[3]) )/pow(Amplitude, 2), 2) );
			else
				intensity = P[2];
			break;
		case 5:
			if( NumberOfParameter == 1 )
                                intensity = sqrt( pow(cos(P[3])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*(cos(P[5]) + sin(P[5])*( sin(P[7]) + cos(P[7]) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*cos(P[3])*(sin(P[5]) - cos(P[5])*(sin(P[7]) + cos(P[7]) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[8]*sin(P[3])*cos(P[3])*sin(P[5])*(sin(P[7])-cos(P[7]))/pow(Amplitude, 2), 2) );
                        else if( NumberOfParameter == 2 )
                                intensity = sqrt( pow(sin(P[3])*cos(P[5])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*cos(P[5])/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*(sin(P[3])*cos(P[7]) + cos(P[3])*sin(P[5]) + sin(P[3])*sin(P[7]) )/pow(Amplitude, 2), 2) + pow( P[1]*P[8]*pow(sin(P[3]),2)*cos(P[5])*sin(P[5])*(sin(P[7])-cos(P[7]))/pow(Amplitude, 2), 2) );
                        else if( NumberOfParameter == 3 )
                                intensity = sqrt( pow(sin(P[3])*sin(P[5])*cos(P[7])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*sin(P[5])*cos(P[7])/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*cos(P[7])*(cos(P[3])*cos(P[5]) + sin(P[3]) )/pow(Amplitude, 2), 2) + pow( P[1]*P[8]*sin(P[3])*sin(P[5])*( cos(P[3])*sin(P[7]) + sin(P[3])*(sin(P[5]) + cos(P[5])*sin(P[7]) ) )/pow(Amplitude, 2), 2) );
			else if( NumberOfParameter == 4 )
				intensity = sqrt( pow(sin(P[3])*sin(P[5])*sin(P[7])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*sin(P[5])*sin(P[7])/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*sin(P[7])*(cos(P[3])*cos(P[5]) + sin(P[3]) )/pow(Amplitude, 2), 2) + pow( P[1]*P[8]*sin(P[3])*sin(P[5])*( cos(P[3])*cos(P[7]) + sin(P[3])*(sin(P[5]) + cos(P[5])*cos(P[7]) ) )/pow(Amplitude, 2), 2) );
			else
				intensity = P[2];
			break;
		case 6:
			if( NumberOfParameter == 1 )
                                intensity = sqrt( pow(cos(P[3])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*(cos(P[5]) + sin(P[5])*( cos(P[7]) + sin(P[7])*( cos(P[9]) + sin(P[9]) ) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*cos(P[3])*(- sin(P[5]) + cos(P[5])*( cos(P[7]) + sin(P[7])*(sin(P[9]) + cos(P[9]) ) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[8]*sin(P[3])*cos(P[3])*sin(P[5])*( - sin(P[7])+ cos(P[7])*( cos(P[9]) + sin(P[9]) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[10]*sin(P[3])*cos(P[3])*sin(P[5])*sin(P[7])*( cos(P[9]) - sin(P[9]) )/pow(Amplitude, 2), 2) );
                        else if( NumberOfParameter == 2 )
                                intensity = sqrt( pow(sin(P[3])*cos(P[5])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*cos(P[5])/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*( sin(P[3])*cos(P[7]) + cos(P[3])*sin(P[5]) + sin(P[3])*sin(P[7])*(cos(P[9]) + sin(P[9]) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[8]*pow(sin(P[3]),2)*cos(P[3])*sin(P[5])*( - sin(P[7]) + cos(P[7])*(cos(P[9]) + sin(P[9]) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[10]*pow(sin(P[3]),2)*cos(P[5])*sin(P[5])*sin(P[7])*(cos(P[9]) - sin(P[9]) )/pow(Amplitude, 2), 2) );
                        else if( NumberOfParameter == 3 )
                                intensity = sqrt( pow(sin(P[3])*sin(P[5])*cos(P[7])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*sin(P[5])*cos(P[7])/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*cos(P[7])*(cos(P[3])*cos(P[5]) + sin(P[3]) )/pow(Amplitude, 2), 2) + pow( P[1]*P[8]*sin(P[3])*sin(P[5])*( sin(P[3])*sin(P[5])*cos(P[9]) + cos(P[3])*sin(P[7]) + sin(P[3])*( cos(P[5])*sin(P[7]) + sin(P[5])*sin(P[9]) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[10]*pow(sin(P[3]),2)*pow(sin(P[5]),2)*cos(P[5])*sin(P[7])*( cos(P[9]) - sin(P[9]) )/pow(Amplitude, 2), 2) );
			else if( NumberOfParameter == 4 )
				intensity = sqrt( pow(sin(P[3])*sin(P[5])*sin(P[7])*cos(P[9])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*sin(P[5])*sin(P[7])*cos(P[9])/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*sin(P[7])*cos(P[9])*(cos(P[3])*cos(P[5]) + sin(P[3]) )/pow(Amplitude, 2), 2) + pow( P[1]*P[8]*sin(P[3])*sin(P[5])*cos(P[9])*( cos(P[3])*cos(P[7]) + sin(P[3])*( sin(P[5]) + cos(P[5])*cos(P[7]) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[10]*sin(P[3])*sin(P[5])*sin(P[7])*( cos(P[3])*sin(P[9]) + sin(P[3])*( cos(P[5])*sin(P[9]) + sin(P[5])*( sin(P[7]) + cos(P[7])*sin(P[9]) ) ) )/pow(Amplitude, 2), 2) );
			else if( NumberOfParameter == 5 )
				intensity = sqrt( pow(sin(P[3])*sin(P[5])*sin(P[7])*sin(P[9])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*sin(P[5])*sin(P[7])*sin(P[9])/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*sin(P[7])*sin(P[9])*(cos(P[3])*cos(P[5]) + sin(P[3]) )/pow(Amplitude, 2), 2) + pow( P[1]*P[8]*sin(P[3])*sin(P[5])*sin(P[9])*( cos(P[3])*cos(P[7]) + sin(P[3])*( sin(P[5]) + cos(P[5])*cos(P[7]) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[10]*sin(P[3])*sin(P[5])*sin(P[7])*( cos(P[3])*cos(P[9]) + sin(P[3])*( cos(P[5])*cos(P[9]) + sin(P[5])*( sin(P[7]) + cos(P[7])*cos(P[9]) ) ) )/pow(Amplitude, 2), 2) );
			else
				intensity = P[2];
			break;
		case 7:
			/*if( NumberOfParameter == 1 )
                                intensity = sqrt( pow(cos(P[3])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*(cos(P[5]) + sin(P[5])*( cos(P[7]) + sin(P[7])*( cos(P[9]) + sin(P[9])*( cos(P[11]) + sin(P[11]) ) ) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*cos(P[3])*(- sin(P[5]) + cos(P[5])*( cos(P[7]) + sin(P[7])*( cos(P[9]) + sin(P[9])*( cos(P[11]) + sin(P[11]) ) ) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[8]*sin(P[3])*cos(P[3])*sin(P[5])*( - sin(P[7])+ cos(P[7])*( cos(P[9]) + sin(P[9])*( cos(P[11]) + sin(P[11]) ) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[10]*sin(P[3])*cos(P[3])*sin(P[5])*sin(P[7])*( - sin(P[9]) + cos(P[9])*( cos(P[11]) + sin(P[11]) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[12]*sin(P[3])*cos(P[3])*sin(P[5])*sin(P[7])*sin(P[9])( - sin(P[11]) + cos(P[11]) )/pow(Amplitude, 2), 2) );
                        else if( NumberOfParameter == 2 )
                                intensity = sqrt( pow(sin(P[3])*cos(P[5])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*cos(P[5])/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*( sin(P[3])*cos(P[7]) + cos(P[3])*sin(P[5]) + sin(P[3])*sin(P[7])*(cos(P[9]) + sin(P[9])*( cos(P[11]) + sin(P[11]) ) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[8]*pow(sin(P[3]),2)*sin(P[5])*cos(P[5])*( - sin(P[7]) + cos(P[7])*(cos(P[9]) + sin(P[9])*( cos(P[11]) + sin(P[11]) ) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[10]*pow(sin(P[3]),2)*cos(P[5])*sin(P[5])*sin(P[7])*( - sin(P[9]) + cos(P[9])*( cos(P[11]) + sin(P[11]) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[12]*pow(sin(P[3]),2)*cos(P[5])*sin(P[5])*sin(P[7])*sin(P[9])( - sin(P[11]) + cos(P[11]) )/pow(Amplitude, 2), 2) );
                        else if( NumberOfParameter == 3 )
                                intensity = sqrt( pow(sin(P[3])*sin(P[5])*cos(P[7])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*sin(P[5])*cos(P[7])/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*cos(P[7])*(cos(P[3])*cos(P[5]) + sin(P[3]) )/pow(Amplitude, 2), 2) + pow( P[1]*P[8]*sin(P[3])*sin(P[5])*( sin(P[3])*sin(P[5])*cos(P[9]) + cos(P[3])*sin(P[7]) + sin(P[3])*( cos(P[5])*sin(P[7]) + sin(P[5])*sin(P[9])*( cos(P[11]) + sin(P[11]) ) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[10]*pow(sin(P[3]),2)*pow(sin(P[5]),2)*cos(P[7])*sin(P[7])*( - sin(P[9]) + cos(P[9])*( cos(P[11]) + sin(P[11]) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[12]*pow(sin(P[3]),2)*pow(sin(P[5]),2)*sin(P[7])*cos(P[7])*sin(P[9])( - sin(P[11]) + cos(P[11]) )/pow(Amplitude, 2), 2) );
			else if( NumberOfParameter == 4 )
				intensity = sqrt( pow(sin(P[3])*sin(P[5])*sin(P[7])*cos(P[9])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*sin(P[5])*sin(P[7])*cos(P[9])/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*sin(P[5])*cos(P[9])*(cos(P[3])*cos(P[5]) + sin(P[3]) )/pow(Amplitude, 2), 2) + pow( P[1]*P[8]*sin(P[3])*sin(P[5])*cos(P[9])*( cos(P[3])*cos(P[7]) + sin(P[3])*( sin(P[5]) + cos(P[5])*cos(P[7]) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[10]*sin(P[3])*sin(P[5])*sin(P[7])*( sin(P[3])*sin(P[5])*sin(P[7])*cos(P[11]) + cos(P[3])*sin(P[9]) + sin(P[3])*( cos(P[5])*sin(P[9]) + sin(P[5])*( sin(P[7])*cos(P[11]) + cos(P[7])*sin(P[9]) ) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[12]*pow(sin(P[3]),2)*pow(sin(P[5]),2)*pow(sin(P[7]),2)*sin(P[9])*cos(P[9])( - sin(P[11]) + cos(P[11]) )/pow(Amplitude, 2), 2) );
			else if( NumberOfParameter == 5 )
				intensity = sqrt( pow(sin(P[3])*sin(P[5])*sin(P[7])*sin(P[9])*cos(P[11])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*sin(P[5])*sin(P[7])*sin(P[9])*cos(P[11])/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*sin(P[7])*sin(P[9])*cos(P[11])*(cos(P[3])*cos(P[5]) + sin(P[3]) )/pow(Amplitude, 2), 2) + pow( P[1]*P[8]*sin(P[3])*sin(P[5])*cos(P[11])*( cos(P[3])*cos(P[7]) + sin(P[3])*sin(P[9])*( sin(P[5]) + cos(P[5])*cos(P[7]) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[10]*sin(P[3])*sin(P[5])*sin(P[7])*cos(P[11])*( cos(P[3])*cos(P[9]) + sin(P[3])*( cos(P[5])*cos(P[9]) + sin(P[5])*( sin(P[7]) + cos(P[7])*cos(P[9]) ) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[12]*sin(P[3])*sin(P[5])*sin(P[7])*sin(P[9])*( cos(P[3])*sin(P[11]) + sin(P[3])*( cos(P[5])*sin(P[11]) + sin(P[5])*( cos(P[7])*sin(P[11]) + sin(P[7])*( cos(P[9])*sin(P[11]) + sin(P[9]) ) ) ) )/pow(Amplitude, 2), 2) );
			else if( NumberOfParameter == 6 )
				intensity = sqrt( pow(sin(P[3])*sin(P[5])*sin(P[7])*sin(P[9])*sin(P[11])*P[2]/Amplitude, 2) + pow( P[1]*P[4]*sin(P[5])*sin(P[7])*sin(P[9])*sin(P[11])/pow(Amplitude, 2), 2) + pow( P[1]*P[6]*sin(P[3])*sin(P[7])*sin(P[9])*sin(P[11])*(cos(P[3])*cos(P[5]) + sin(P[3]) )/pow(Amplitude, 2), 2) + pow( P[1]*P[8]*sin(P[3])*sin(P[5])*( cos(P[3])*cos(P[7]) + sin(P[3])*sin(P[9])*sin(P[11])*( sin(P[5]) + cos(P[5])*cos(P[7]) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[10]*sin(P[3])*sin(P[5])*sin(P[7])*( cos(P[3])*cos(P[9]) + sin(P[3])*sin(P[11])*( cos(P[5])*cos(P[9]) + sin(P[5])*( sin(P[7]) + cos(P[7])*cos(P[9]) ) ) )/pow(Amplitude, 2), 2) + pow( P[1]*P[12]*sin(P[3])*sin(P[5])*sin(P[7])*sin(P[9])*( cos(P[3])*cos(P[11]) + sin(P[3])*( cos(P[5])*cos(P[11]) + sin(P[5])*( cos(P[7])*cos(P[11]) + sin(P[7])*( cos(P[9])*cos(P[11]) + sin(P[9]) ) ) ) )/pow(Amplitude, 2), 2) );
			else*/ // Need proper look
				intensity = P[2];
			break;
		case 8:
			if( NumberOfParameter == 1 )
                                intensity = P[2];
                        else if( NumberOfParameter == 2 )
                                intensity = P[4];
                        else if( NumberOfParameter == 3 )
                                intensity = P[6];
			else if( NumberOfParameter == 4 )
				intensity = P[8];
			else if( NumberOfParameter == 5 )
				intensity = P[10];
			else if( NumberOfParameter == 6 )
				intensity = P[12];
			else if( NumberOfParameter == 7 )
				intensity = P[14];
			else
				intensity = P[2];
			break;
		case 9:
			if( NumberOfParameter == 1 )
                                intensity = P[2];
                        else if( NumberOfParameter == 2 )
                                intensity = P[4];
                        else if( NumberOfParameter == 3 )
                                intensity = P[6];
			else if( NumberOfParameter == 4 )
				intensity = P[8];
			else if( NumberOfParameter == 5 )
				intensity = P[10];
			else if( NumberOfParameter == 6 )
				intensity = P[12];
			else if( NumberOfParameter == 7 )
				intensity = P[14];
			else if( NumberOfParameter == 8 )
				intensity = P[16];
			else
				intensity = P[2];
			break;
		case 10:
			if( NumberOfParameter == 1 )
                                intensity = P[2];
                        else if( NumberOfParameter == 2 )
                                intensity = P[4];
                        else if( NumberOfParameter == 3 )
                                intensity = P[6];
			else if( NumberOfParameter == 4 )
				intensity = P[8];
			else if( NumberOfParameter == 5 )
				intensity = P[10];
			else if( NumberOfParameter == 6 )
				intensity = P[12];
			else if( NumberOfParameter == 7 )
				intensity = P[14];
			else if( NumberOfParameter == 8 )
				intensity = P[16];
			else if( NumberOfParameter == 9 )
				intensity = P[18];
			else
				intensity = P[2];
			break;
	}
	return (Double_t)intensity;
}


Double_t DiscreteFitFunctionNoPs_forTests( Double_t *A, Double_t *P )	//First parameter - nmbr of comp, Second - nmb of Gauss
{
	Double_t sum = 0., FixedIntensity = 0.;	
	for( unsigned i = 0; i < P[1]; i++ )
	{
		sum += P[2+i];
	}
	return sum + P[0];
}
