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
                                                            double FixedIntensities, double FreeIntensitiesTopPs, double FixedFixedIntensity, std::vector<LifetimeComponent> Lifetimes )
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
	ResultsFromDiscrete.push_back( LifetimesFromDiscrete );
    ResultsFromDiscrete.push_back( IntensitiesFromDiscrete );
    ResultsFromDiscrete.push_back( ResolutionFromDiscrete );
    ResultsFromDiscrete.push_back( OffsetsFromDiscrete );
    ResultsFromDiscrete.push_back( FractionsFromDiscrete );    
	return ResultsFromDiscrete;
}

void saveDiscreteFitResultsToEXCEL()
{
    
}
