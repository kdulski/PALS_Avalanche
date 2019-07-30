#include "PALS_avalanche_resultsSaver.h"

ResultsSaver::ResultsSaver( std::string TypeOfFittedModel )
{
    TypeOfFit = TypeOfFittedModel;
}

void ResultsSaver::saveDiscreteFit( TH1F *histogram, TF1 *DiscreteFit, const char *RootFile, std::string PathOfFile, std::string ResultsPath, std::string PathOfFileWithDate )
{
	TCanvas *c1 = new TCanvas( "c1", "", 710, 500 );

	c1 -> SetFillColor( 0 ); //Option for histogram
	c1 -> SetFrameBorderMode( 0 );
	c1 -> SetBorderSize( 2 );
	c1 -> SetFrameLineWidth( 2 );
	c1 -> SetLeftMargin(0.12);
	c1 -> SetRightMargin(0.08);
	c1 -> SetBorderMode(0);

	gStyle -> SetOptStat(10);
	gStyle -> SetStatX(0.92);
	gStyle -> SetStatY(0.92);
	gStyle -> SetStatW(0.2);
	gStyle -> SetStatH(0.15);
	gStyle -> SetStatBorderSize(2);
	gStyle->SetOptFit(0000);			//Option for fitting parameters showed in statisticbox
    
	histogram -> Draw("");
	DiscreteFit -> Draw("C same");
	
	TFile *resultsFile = new TFile( RootFile, "update" );
	testFile -> mkdir( PathOfFile.c_str() );
	testFile -> cd( PathOfFile.c_str() );
		    
	if( not fileTools.PathCheck( ResultsPath ) )
	{
		std::cout << "Creating directory for results" << std::endl;
		boost::filesystem::path dir( ResultsPath );
		boost::filesystem::create_directories( dir );
	}
	
	c1->SaveAs( ( ResultsPath + "/" + PathOfFileWithDate + ".png" ).c_str() );
	c1->SetLogy();
	c1->SaveAs( ( ResultsPath + "/" + PathOfFileWithDate + "_logscale"+".png" ).c_str() );
	c1 -> Write( PathOfFileWithDate.c_str() );
	
	delete resultsFile;
	delete c1;
}

std::vector< std::vector< DiscreteFitResult > > ResultsSaver::saveDiscreteFitResultsToTXT( std::string Path, std::string Prefix, std::string PathWithDate, 
                                                            TF1* Discrete, double Background, double SDBackground, Double_t* ResolutionsFromFit, Double_t* ResolutionsFromFitErrors, 
                                                            int pPsIndex, double pPsIntensity, Double_t* FreeParameters, Double_t* FreeParametersErrors, 
                                                            double FixedIntensities, double FreeIntensitiesTopPs, double FixedFixedIntensity, std::vector<LifetimeComponent> Lifetimes, 
                                                            double oPsLifetimeCutoff, double pPsLifetimeCutoff, std::string PathForExcel )
{
	std::vector< std::vector< DiscreteFitResult > > ResultsFromDiscrete;
	std::vector< DiscreteFitResult > LifetimesFromDiscrete;
	std::vector< DiscreteFitResult > IntensitiesFromDiscrete;		
	std::vector< DiscreteFitResult > ResolutionFromDiscrete;
	std::vector< DiscreteFitResult > OffsetsFromDiscrete;
	std::vector< DiscreteFitResult > FractionsFromDiscrete;
	int NotFixedIterator = 0;
	std::ofstream res;
	res.open( Res_path + "/" + Prefix + PathWithDate );
	res << "Results from fitting" << std::endl;
	res << "Parameter name\t \t \t \t \t" << "Value\t \t \t" << "Error" << std::endl;
	res << std::endl;
	res << std::endl;
	res << "Background\t \t \t \t \t" << NumberToChar( Background, 5 ) << "\t \t" << SDBackground << std::endl;
	unsigned ResSize = (unsigned)ResolutionsFromFit[1];
	for( unsigned i = 0; i < ResSize; i++ )
	{
		res << ("Sigma for " + NumberToChar( i+1, 0 ) + " Gauss [ns]").c_str() 
			  << "\t \t \t \t " << NumberToChar( Discrete -> GetParameter( 4 + 3*i ), 5) << "\t \t" << Discrete -> GetParError( 4 + 3*i ) << std::endl;
		ResolutionFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 4 + 3*i ), Discrete -> GetParError( 4 + 3*i )  ) );
		res << ("FWHM for " + NumberToChar( i+1, 0 ) + " Gauss [ns]").c_str() 
			  << "\t \t \t \t " << NumberToChar( 2.355*Discrete -> GetParameter( 4 + 3*i ), 5 ) << "\t \t" << 2.355*Discrete -> GetParError( 4 + 3*i ) << std::endl;
		res << ("Fraction for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() 
			  << "\t \t \t \t " << NumberToChar( GetIntensityParameterNew( ResolutionsFromFit, 3, 3, i+1 ),5 ) << "\t \t" << GetIntensityParameterErrorNew( ResolutionsFromFitErrors, i+1 ) << std::endl;
		FractionsFromDiscrete.push_back( DiscreteFitResult( GetIntensityParameterNew( ResolutionsFromFit, 3, 3, i+1 ), 0 ) );
		res << ("Offset for " + NumberToChar( i+1, 0 ) + " Gauss [ns]").c_str() 
			  << "\t \t \t \t " << NumberToChar( Discrete -> GetParameter( 6 + 3*i ), 5 ) << "\t \t" << Discrete -> GetParError( 6 + 3*i ) << std::endl;
		OffsetsFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 6 + 3*i ), Discrete -> GetParError( 6 + 3*i ) ) );
		res << std::endl;
	}
	res << std::endl;
	if( pPsIndex + 1 )
	{
		unsigned j = pPsIndex;
		res << "Lifetime for p-PS Component\t \t \t " << NumberToChar( Discrete -> GetParameter( 4 + 3*ResSize + 2*j ), 5) << "\t \t" << Discrete -> GetParError( 4 + 3*ResSize + 2*j ) << std::endl;
		LifetimesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 4 + 3*ResSize + 2*j ), Discrete -> GetParError( 4 + 3*ResSize + 2*j ) ) );
		res << "Intensity for p-PS Component \t \t \t " <<  NumberToChar( pPsIntensity, 5) << std::endl;
		res << "Intensity for p-PS Component in percent \t " <<  NumberToChar( pPsIntensity * 100, 5) << std::endl;
		res << "Intensity for p-PS Component in percent no fix" << "\t " <<  NumberToChar( pPsIntensity * 100/(1-FixedFixedIntensity), 5 ) << std::endl;	
		res << "Intensity for p-PS Component in percent of o-Ps\t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*ResSize ) * 100, 5)  << std::endl;
		IntensitiesFromDiscrete.push_back( DiscreteFitResult( pPsIntensity, 0 ) );
		res << std::endl;
	}
	for( unsigned j = 0; j < Lifetimes.size(); j++ )
	{
		if( j < FixedIterator + pPsIndex + 1 && (int)j != pPsIndex )
		{
			res << ("Lifetime for " + NumberToChar( j+1, 0 ) + " Component [ns]").c_str()
				  << "\t \t \t " << NumberToChar( Discrete -> GetParameter( 4 + 3*ResSize + 2*j ), 5) << "\t \t" << Discrete -> GetParError( 4 + 3*ResSize + 2*j ) << std::endl;
			LifetimesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 4 + 3*ResSize + 2*j ), Discrete -> GetParError( 4 + 3*ResSize + 2*j ) ) );
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component").c_str() 
				  << "\t \t \t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*ResSize + 2*j ), 5) << "\t \t" << 
					        Discrete -> GetParError( 5 + 3*ResSize + 2*j )  << std::endl;
			IntensitiesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 5 + 3*ResSize + 2*j ), Discrete -> GetParError( 5 + 3*ResSize + 2*j ) ) );
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component in percent").c_str()
		                  << "\t \t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*ResSize + 2*j ) * 100, 5 ) << "\t \t" <<
		                               Discrete -> GetParError( 5 + 3*ResSize + 2*j ) * 100 << std::endl;
			if( Lifetimes[j].Type != "f" )
			{
				  res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component in percent no fix").c_str()
		                  << "\t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*ResSize + 2*j ) * 100/(1-FixedFixedIntensity), 5 ) << "\t \t" <<
		                               Discrete -> GetParError( 5 + 3*ResSize + 2*j ) * 100/(1-FixedFixedIntensity) << std::endl;
			}
			else
			{
				LifetimesFromDiscrete[ LifetimesFromDiscrete.size() - 1 ].Type = "f";
			}
		        std::cout << "Intensity in percent: \t" << Discrete -> GetParameter( 5 + 3*ResSize + 2*j ) * 100 << std::endl;
			res << std::endl;
		}
		else if( (int)j != pPsIndex )
		{
			NotFixedIterator++;
			res << ("Lifetime for " + NumberToChar( j+1, 0 ) + " Component [ns]").c_str()
				  << "\t \t \t " << NumberToChar( Discrete -> GetParameter( 4 + 3*ResSize + 2*j ), 5) << "\t \t" << Discrete -> GetParError( 4 + 3*ResSize + 2*j ) << std::endl;
			LifetimesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 4 + 3*ResSize + 2*j ), Discrete -> GetParError( 4 + 3*ResSize + 2*j ) ) );
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component").c_str() 
				  << "\t \t \t " <<  NumberToChar( GetIntensityParameterNew( FreeParameters, 3, 3, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*ResSize ) ) , 5) << "\t \t" << GetIntensityParameterErrorNew( FreeParametersErrors, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*ResSize ) ) << std::endl;	  
			IntensitiesFromDiscrete.push_back( DiscreteFitResult( GetIntensityParameterNew( FreeParameters, 3, 3, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*ResSize ) ), GetIntensityParameterErrorNew( FreeParametersErrors, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*ResSize ) ) ) );
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component in percent").c_str()
		                  << "\t \t " <<  NumberToChar( GetIntensityParameterNew( FreeParameters, 3, 3, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*ResSize ) ) * 100, 5 ) << "\t \t" << GetIntensityParameterErrorNew( FreeParametersErrors, NotFixedIterator ) *(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*ResSize ) ) * 100 << std::endl;
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component in percent no fix").c_str()
		                  << "\t " <<  NumberToChar( GetIntensityParameterNew( FreeParameters, 3, 3, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*ResSize ) ) * 100/(1-FixedFixedIntensity), 5 ) << "\t \t" << GetIntensityParameterErrorNew( FreeParametersErrors, NotFixedIterator ) *(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*ResSize ) ) * 100/(1-FixedFixedIntensity) << std::endl;	  
		        std::cout << "Intensity in percent: \t" << GetIntensityParameterNew( FreeParameters, 3, 3, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*ResSize ) ) * 100 << std::endl;
			res << std::endl;
		}
	}
	
//---------------------------------------Sorting the lifetime components, so the first component has the lowest lifetime, second, the second lowest, ...
//---------------------------------------Update 30.07.2019 -> I know that I have soritng algorithm implemented, which in gemneral I could use for soritng lifetime components
//---------------------------------------but this method is quicker and I do not have time to rewrite this part, which is working correctly
	
	std::vector<unsigned> Order;
	double minTemp = 200, minPrevious = 0;
	double minTempOld = minTemp;
	unsigned minTempIndex;
	bool minSearch = true;
	while( minSearch )
	{
		for( unsigned i=0; i<LifetimesFromDiscrete.size(); i++ )
		{
			if( LifetimesFromDiscrete[i].Parameter < minTemp && LifetimesFromDiscrete[i].Parameter > minPrevious )
			{
				minTemp = LifetimesFromDiscrete[i].Parameter;
				minTempIndex = i;
			}
		}
		Order.push_back( minTempIndex );
		minPrevious = minTemp;
		minTemp = minTempOld;
		if( Order.size() == LifetimesFromDiscrete.size() )
			minSearch = false;
	}
//---------------------------------------
	
	double DirectAnnihilation3GFraction = 372;
	double MaximaloPSLifetime = 142;
	double oPs_pPs_Fraction = 3;
	double Frac2GParameterNew = 0, Frac3GParameterNew = 0, DeltaParameter = 0;
	if( LifetimesFromDiscrete.size() == 2 )
	{
		Frac2GParameterNew = ( DirectAnnihilation3GFraction - 1 ) / DirectAnnihilation3GFraction;
		Frac3GParameterNew = 1 / DirectAnnihilation3GFraction;
	}
	else
	{
		for( unsigned i=1; i<LifetimesFromDiscrete.size(); i++ )
		{
			if( LifetimesFromDiscrete[ Order[i] ].Parameter < oPsLifetimeCutoff && LifetimesFromDiscrete[ Order[i] ].Parameter > pPsLifetimeCutoff )
			{
				Frac2GParameterNew += ( DirectAnnihilation3GFraction - 1 )*IntensitiesFromDiscrete[ Order[i] ].Parameter / DirectAnnihilation3GFraction;
				Frac3GParameterNew += IntensitiesFromDiscrete[ Order[i] ].Parameter / DirectAnnihilation3GFraction;
			}
			else if( LifetimesFromDiscrete[ Order[i] ].Parameter > oPsLifetimeCutoff )
			{
				Frac2GParameterNew += IntensitiesFromDiscrete[ Order[i] ].Parameter  / oPs_pPs_Fraction;
				Frac2GParameterNew += IntensitiesFromDiscrete[ Order[i] ].Parameter * ( MaximaloPSLifetime - LifetimesFromDiscrete[ Order[i] ].Parameter ) / MaximaloPSLifetime;
				Frac3GParameterNew += IntensitiesFromDiscrete[ Order[i] ].Parameter * LifetimesFromDiscrete[ Order[i] ].Parameter / MaximaloPSLifetime;
				DeltaParameter += (1-(oPs_pPs_Fraction + 1)*IntensitiesFromDiscrete[ Order[i] ].Parameter / oPs_pPs_Fraction)/DirectAnnihilation3GFraction + IntensitiesFromDiscrete[ Order[i] ].Parameter * LifetimesFromDiscrete[ Order[i] ].Parameter / MaximaloPSLifetime;
			}
		}
	}
	res << "Lifetime cutoffs for para- and ortho-PS " << pPsLifetimeCutoff << " " << oPsLifetimeCutoff << std::endl;
	res << "Below parameters can vary depending on the cutoffs chosen for analysis, cuttofs are given in PALS_avalanche_lib_h header" << std::endl;
	res << "Frac2GParameterNew \t" << Frac2GParameterNew << "\t" << "Number of annihilations into 2 gamma quanta based on the oPs components and direct annihilation" << std::endl;	
	res << "Frac3GParameternew \t" << Frac3GParameterNew << "\t" << "Number of annihilations into 3 gamma quanta based on the oPs components and direct annihilation" << std::endl;	
	res << "DeltaParameter \t" << DeltaParameter << "\t" << "Parameter taken from Acta Physica Polonica B48 (2017) 1577" << std::endl;				
	
	double sum = 0.;
	unsigned sumIT = 0;
	for( int i = 0; i < histogram->GetNbinsX(); i++ )
	{
		if( histogram->GetBinCenter(i) >= Arguments[ Range_From ] && histogram->GetBinCenter(i) <= Arguments[ Range_To ] )
		{
			sum += histogram->GetBinContent(i);
			sumIT++;
		}	
	}
	sum /= sumIT;
	res << std::endl;	
	res << "Statistical parameters of the fit, calculated from the histogram" << std::endl;
	std::cout << "Mean of Histogram \t" << sum << std::endl;
	std::cout << "Degrees of freedom of data \t" << sumIT - 1 << std::endl;
	std::cout << "Degrees of freedom of fit \t" << Discrete -> GetNDF() << std::endl;
	double Diff = 0.;
	for( int i = 0; i < histogram->GetNbinsX(); i++ )
	{
		if( histogram->GetBinCenter(i) >= Arguments[ Range_From ] && histogram->GetBinCenter(i) <= Arguments[ Range_To ] )
		{		
			if( histogram->GetBinContent(i) > 0 )			
				Diff += (sum - histogram->GetBinContent(i) )*( sum - histogram->GetBinContent(i) )/histogram->GetBinContent(i);
			else
				Diff += (sum - 1 )*( sum - 1 );
		}
	}
	std::cout << "Diff between mean of histogram and values \t" << Diff << std::endl;

	res << "ChiSquared" << "\t \t \t \t \t" << NumberToChar( Discrete -> GetChisquare(), 4 ) << std::endl;
	res << "Probability of the fit " << "\t \t \t \t" << NumberToChar(1 - Discrete -> GetProb(), 4) << std::endl;
	res << "Rsquared" << "\t \t \t \t \t" << NumberToChar(1 - (Discrete -> GetChisquare())/Diff, 4 ) << std::endl;
	res << "Adjusted Rsquared" << "\t \t \t \t" << NumberToChar(1 - ((sumIT-1)/Discrete -> GetNDF())*(Discrete -> GetChisquare())/Diff, 4 ) << std::endl;
	res << "ChiSquared/DegreesOfFreedom" << "\t \t \t" << NumberToChar( Discrete -> GetChisquare()/Discrete -> GetNDF(), 4 ) << std::endl;
	res.close();
	
	std::vector< DiscreteFitResult > StatisticalParams;
	StatisticalParams.push_back( DiscreteFitResult, Frac2GParameterNew, 0 );
	StatisticalParams.push_back( DiscreteFitResult, Frac3GParameterNew, 0 );
	StatisticalParams.push_back( DiscreteFitResult, DeltaParameter, 0 );
	
	StatisticalParams.push_back( DiscreteFitResult, Discrete -> GetChisquare(), 0 );
	StatisticalParams.push_back( DiscreteFitResult, 1 - Discrete -> GetProb(), 0 );
	StatisticalParams.push_back( DiscreteFitResult, 1 - (Discrete -> GetChisquare())/Diff, 0 );
	StatisticalParams.push_back( DiscreteFitResult, 1 - ((sumIT-1)/Discrete -> GetNDF())*(Discrete -> GetChisquare())/Diff, 0 );
	StatisticalParams.push_back( DiscreteFitResult, Discrete -> GetChisquare()/Discrete -> GetNDF(), 0 );
	
//---------------------------------------Writing to csv	
	bool IsFileExists =  FileCheck( Res_path + "/" + "Discrete_Fit_" + PathForExcel + ".csv" );
	std::ofstream res_csv;
	if( IsFileExists )
		res_csv.open( Res_path + "/" + "Discrete_Fit_" + PathForExcel + ".csv", std::ofstream::app );
	else
		res_csv.open( Res_path + "/" + "Discrete_Fit_" + PathForExcel + ".csv" );
	res_csv << PathWithDate;	
	res_csv << ",Frac2GParameternew,Frac3GParameternew,DeltaParameter,ChiSquared,Probability_of_the_fit,Rsquared,Adjusted_Rsquared,ChiSquared/DegreesOfFreedom,";
	res_csv << "Background,Standard_deviation_of_Background,";
	for( unsigned i = 0; i < Resolution.size(); i++ )
	{
		res_csv << "Sigma_for_" << NumberToChar( i+1, 0 ) << "_Gauss_[ns],";
		res_csv << "Error_of_Sigma_for_" << NumberToChar( i+1, 0 ) << "_Gauss_[ns],";
		res_csv << "FWHM_for_" << NumberToChar( i+1, 0 ) << "_Gauss_[ns],";
		res_csv << "Error_of_FWHM_for_" << NumberToChar( i+1, 0 ) << "_Gauss_[ns],";
		res_csv << "Fraction_for_" << NumberToChar( i+1, 0 ) << "_Gauss_[ns],";
		res_csv << "Error_of_Fraction_for_" << NumberToChar( i+1, 0 ) << "_Gauss_[ns],";
		res_csv << "Offset_for_" << NumberToChar( i+1, 0 ) << "_Gauss_[ns],";
		res_csv << "Error_of_Offset_for_" << NumberToChar( i+1, 0 ) << "_Gauss_[ns],";
	}
	res_csv << ",";
	if( pPsIndex + 1 )
	{
		res_csv << "Lifetime_for_p-PS_Component,";
		res_csv << "Intensity_for_p-PS_Component,";
		res_csv << "Intensity_for_p-PS_Component_in_percent,";
		res_csv << "Intensity_for_p-PS_Component_in_percent_no_fix,";
		res_csv << "Intensity_for_p-PS_Component_in_percent_of_o-Ps,";
	}
	
	double FixedIntensity = 0;
		
	for( unsigned j = 0; j < Lifetimes.size(); j++ )
	{
		if( LifetimesFromDiscrete[ Order[j] ].Type == "f" )
		{
			res_csv << "Lifetime_for_" << NumberToChar( j+1, 0 ) << "_Component_[ns],";
			res_csv << "Error_of_Lifetime_for_" << NumberToChar( j+1, 0 )<< "_Component_[ns],";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component,";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent,";
			FixedIntensity += IntensitiesFromDiscrete[ Order[j] ].Parameter;
		}
		else
		{
			res_csv << "Lifetime_for_" << NumberToChar( j+1, 0 )<< "_Component_[ns],";
			res_csv << "Error_of_Lifetime_for_" << NumberToChar( j+1, 0 )<< "_Component_[ns],";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component,";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent,";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent_no_fix,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent_no_fix,";
		}
	}
	res_csv << "\n";
	res_csv << ",";
	res_csv << NumberToChar( Frac2GParameterNew, 5 ) << "," << NumberToChar( Frac3GParameterNew, 5 ) << "," << NumberToChar( DeltaParameter, 5 ) << ",";
	res_csv << NumberToChar( Discrete -> GetChisquare(), 4 ) << "," << NumberToChar(1 - Discrete -> GetProb(), 4) << "," << NumberToChar(1 - (Discrete -> GetChisquare())/Diff, 4 ) << "," ;
	res_csv << NumberToChar(1 - ((sumIT-1)/Discrete -> GetNDF())*(Discrete -> GetChisquare())/Diff, 4 ) << "," << NumberToChar(  Discrete -> GetChisquare()/Discrete -> GetNDF(), 4 ) << ",";
	res_csv << NumberToChar( Background, 5 ) << "," << SDBackground << ",";
	for( unsigned i = 0; i < Resolution.size(); i++ )
	{
		res_csv << NumberToChar( Discrete -> GetParameter( 4 + 3*i ), 5) << ",";
		res_csv << Discrete -> GetParError( 4 + 3*i ) << ",";
		res_csv << NumberToChar( 2.355*Discrete -> GetParameter( 4 + 3*i ), 5 ) << ",";
		res_csv << 2.355*Discrete -> GetParError( 4 + 3*i ) << ",";
		res_csv << NumberToChar( GetIntensityParameterNew( ResolutionsFromFit, 3, 3, i+1 ),5 ) << ",";
		res_csv << GetIntensityParameterErrorNew( ResolutionsFromFitErrors, i+1 ) << ",";
		res_csv << NumberToChar( Discrete -> GetParameter( 6 + 3*i ), 5 ) << ",";
		res_csv << Discrete -> GetParError( 6 + 3*i ) << ",";
	}
	res_csv << ",";
	if( pPsIndex + 1 )
	{
		unsigned j = pPsIndex;
		res_csv << NumberToChar( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), 5) << "\t \t" << Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) << ",";
		res_csv << NumberToChar( pPsIntensity, 5) << ",";
		res_csv << NumberToChar( pPsIntensity * 100, 5) << ",";
		res_csv << NumberToChar( pPsIntensity * 100/(1-FixedFixedIntensity), 5 ) << ",";
		res_csv << NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() ) * 100, 5) << ",";
	}
	NotFixedIterator = 0;
	for( unsigned j = 0; j < Lifetimes.size(); j++ )
	{
		if( LifetimesFromDiscrete[ Order[j] ].Type == "f" )
		{
			res_csv << NumberToChar( LifetimesFromDiscrete[ Order[j] ].Parameter, 5) << ",";
			res_csv << NumberToChar( LifetimesFromDiscrete[ Order[j] ].Uncertainity, 5) << ",";
			res_csv << NumberToChar( IntensitiesFromDiscrete[ Order[j] ].Parameter, 5) << ",";
			res_csv << NumberToChar( IntensitiesFromDiscrete[ Order[j] ].Uncertainity, 5) << ",";
			res_csv << NumberToChar( 100*IntensitiesFromDiscrete[ Order[j] ].Parameter, 5) << ",";
			res_csv << NumberToChar( 100*IntensitiesFromDiscrete[ Order[j] ].Uncertainity, 5) << ",";		  
		}
		else
		{
			res_csv << NumberToChar( LifetimesFromDiscrete[ Order[j] ].Parameter, 5) << ",";
			res_csv << NumberToChar( LifetimesFromDiscrete[ Order[j] ].Uncertainity, 5) << ",";
			res_csv << NumberToChar( IntensitiesFromDiscrete[ Order[j] ].Parameter, 5) << ",";
			res_csv << NumberToChar( IntensitiesFromDiscrete[ Order[j] ].Uncertainity, 5) << ",";
			res_csv << NumberToChar( 100*IntensitiesFromDiscrete[ Order[j] ].Parameter, 5) << ",";
			res_csv << NumberToChar( 100*IntensitiesFromDiscrete[ Order[j] ].Uncertainity, 5) << ",";
			res_csv << NumberToChar( 100*IntensitiesFromDiscrete[ Order[j] ].Parameter /(1-FixedIntensity), 5) << ",";
			res_csv << NumberToChar( 100*IntensitiesFromDiscrete[ Order[j] ].Uncertainity /(1-FixedIntensity), 5) << ",";			
		}
	}
	res_csv << "\n";
	res_csv.close();
	
	ResultsFromDiscrete.push_back( LifetimesFromDiscrete );
	ResultsFromDiscrete.push_back( IntensitiesFromDiscrete );
	ResultsFromDiscrete.push_back( ResolutionFromDiscrete );
	ResultsFromDiscrete.push_back( OffsetsFromDiscrete );
	ResultsFromDiscrete.push_back( FractionsFromDiscrete );
	ResultsFromDiscrete.push_back( StatisticalParams );
	
	return ResultsFromDiscrete;
}

void ResultsSaver::saveResiduals( TH1F* histogram, double MinArgument, double MaxArgument, unsigned MinBin, unsigned MaxBin, TF1* Discrete, std::string FileName, std::string Path, std::string PathWithDate )
{
	TCanvas *c2 = new TCanvas( "c2", "", 710, 500 );
        c2 -> SetFillColor( 0 ); //Option for histogram
        c2 -> SetFrameBorderMode( 0 );
        c2 -> SetBorderSize( 2 );
        c2 -> SetFrameLineWidth( 2 );
	c2 -> SetLeftMargin(0.12);
	c2 -> SetRightMargin(0.08);
	c2 -> SetBorderMode(0);

	gStyle -> SetOptStat(10);
	gStyle -> SetStatX(0.92);
	gStyle -> SetStatY(0.92);
	gStyle -> SetStatW(0.2);
	gStyle -> SetStatH(0.15);
	gStyle -> SetStatBorderSize(2);
	Double_t Residuals_root[ MaxBin - MinBin ], ResArg_root[ MaxBin - MinBin ];
	for( int i = 0; i < histogram->GetNbinsX(); i++ )
	{
		if( histogram->GetBinCenter(i) > MinArgument && histogram->GetBinCenter(i) < MaxArgument && i >= (int)MinBin )
		{
			if( isnan(histogram->GetBinContent(i) - Discrete->Eval( histogram->GetBinCenter(i) ) )/sqrt( histogram->GetBinContent(i) ) || isinf(histogram->GetBinContent(i) - Discrete->Eval( histogram->GetBinCenter(i) ) )/sqrt( histogram->GetBinContent(i) ) )
				Residuals_root[i-MinBin] = 0;
			else
				Residuals_root[i-MinBin] = (histogram->GetBinContent(i) - Discrete->Eval( histogram->GetBinCenter(i) ) )/sqrt( histogram->GetBinContent(i) );
			ResArg_root[i-MinBin] = histogram->GetBinCenter(i);
		}
	}

	TFile *residuals = new TFile( FileName, "update" );
	residuals -> mkdir( Path.c_str() ); 
	residuals -> cd( Path.c_str() );  
        TGraph *residualsHisto = new TGraph( MaxBin - MinBin, ResArg_root, Residuals_root );
        residualsHisto -> SetMarkerStyle( 21 );
        residualsHisto -> SetMarkerSize( 0.5 );
        residualsHisto -> SetTitle( "Residuals of the fit" );
        residualsHisto -> Draw( "APL" );
        residualsHisto -> GetXaxis() -> SetTitle( "Time difference [ns]" );
        residualsHisto -> GetXaxis() -> SetTitle( "Residuals" );
        residualsHisto -> Write( PathWithDate.c_str() );
        delete residualsHisto;

	residuals->Close();
	delete c2;
}

void ResultsSaver::saveLFvsIntensities( std::string FileName, std::string Path, std::string PathWithDate, std::vector< DiscreteFitResult > LifetimesFromDiscrete, std::vector< DiscreteFitResult > IntensitiesFromDiscrete )
{
	TFile *comparison = new TFile( FileName, "update" );
	comparison -> mkdir( Path.c_str() ); 
	comparison -> cd( Path.c_str() ); 

	TGraphErrors *LF_vs_In = new TGraphErrors();
	LF_vs_In -> GetXaxis() -> SetTitle( "Mean Lifetime [ns]" );
        LF_vs_In -> GetYaxis() -> SetTitle( "Intensity" );
	for( unsigned i=0; i<LifetimesFromDiscrete.size(); i++ )
	{
		LF_vs_In -> SetPoint( i, LifetimesFromDiscrete[i].Parameter, IntensitiesFromDiscrete[i].Parameter );
		LF_vs_In -> SetPointError( i, LifetimesFromDiscrete[i].Uncertainity, IntensitiesFromDiscrete[i].Uncertainity );
	}
	LF_vs_In -> SetMarkerStyle(20);
	LF_vs_In -> Draw("AP");	
	LF_vs_In -> Write( PathWithDate.c_str() );
	
	delete LF_vs_In;
	comparison->Close();  
}
