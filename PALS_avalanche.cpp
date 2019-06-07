#include "PALS_avalanche_lib.h"

	
int main( int argc, char* argv[] ) 		//First argument - File with data, second - File with FitDetails(Path to it), Third is optional if one would like to check other approaches to fitting (1 -> Gauss simplex intensities old type, other - experimental approach with different fucntion body for diferent type of fittign, works only for up to three Gaussian components)
{
        if( argc > 5 )                                                 //Argument check
	{
                std::cout << "Program is not using any more external arguments, so do not be officious!" << std::endl;
		return 0;
        }                                                               //Argument check							//Start counting time - control purpose
       
        std::string TimesPath = argv[1];                                    //Calculated times
	int DecoOption;
	int ROOTFileTest = 0;
	if( TimesPath.length() > 5 )
	{
		std::string TimesPathExtension;
		TimesPathExtension.insert(0, TimesPath, (TimesPath.length() - 5), TimesPath.length());
		if( TimesPathExtension == ".root" )
			ROOTFileTest = 1;
		else if( argc > 4 )
		{
			std::cout << "Program is not using any more external arguments, so do not be officious!" << std::endl;
			return 0;
		}
	}
        if( argc == 2 )
	{
		if( ROOTFileTest )
		{
			std::cout << "You need to provide the name of histogram for fit from root file" << std::endl;
			std::cout << "Without it program do not know what histogram is for fitting" << std::endl;
			return 0;
		}
		Fit fit( TimesPath, "FitDetails" );                                                      //Getting times from file to object
		fit.RangeBackgroundData();
		DecoOption = fit.Discrete();
		if( DecoOption )
	      		fit.Deconvolution();
	}
	else 
	{
		std::string FitDetailsPath;
		std::string HistoDetails;
		if( argc == 3 || (argc == 4 && ROOTFileTest) )
		{
			if( ROOTFileTest )
			{
				if( argc == 3 )
				{
					HistoDetails = argv[2];
					Fit fit( TimesPath, HistoDetails, "FitDetails" );                                                      //Getting times from file to object
					fit.RangeBackgroundData();
					DecoOption = fit.Discrete();
					if( DecoOption )
						fit.Deconvolution();
				}
				else
				{
					FitDetailsPath = argv[3];
					HistoDetails = argv[2];
					Fit fit( TimesPath, HistoDetails, FitDetailsPath );                                                      //Getting times from file to object
					fit.RangeBackgroundData();
					DecoOption = fit.Discrete();
					if( DecoOption )
						fit.Deconvolution();
				}
			}
			else
			{
				FitDetailsPath = argv[2];
				Fit fit( TimesPath, FitDetailsPath );                                                      //Getting times from file to object
				fit.RangeBackgroundData();
				DecoOption = fit.Discrete();
				if( DecoOption )
					fit.Deconvolution();
			}
		}
		else
		{
			int FitOption;
			if( ROOTFileTest )
				FitOption = atoi(argv[4]);
			else
				FitOption = atoi(argv[3]);
			if( FitOption == 1 )
			{
				if( ROOTFileTest )
				{
					FitDetailsPath = argv[3];
					HistoDetails = argv[2];
					Fit fit( TimesPath, HistoDetails, FitDetailsPath );                                                      //Getting times from file to object
					fit.RangeBackgroundData();
					DecoOption = fit.Discrete_old();
					if( DecoOption )
						fit.Deconvolution();
				}
				else
				{
					FitDetailsPath = argv[2];
					Fit fit( TimesPath, FitDetailsPath );                                                      //Getting times from file to object
					fit.RangeBackgroundData();
					DecoOption = fit.Discrete_old();
					if( DecoOption )
						fit.Deconvolution();				  
				}
			}
			else
			{
				std::cout << std::endl;
				std::cout << std::endl;
				std::cout << "Remember that his experimental approach works for now only for max 3 resolution components!!" << std::endl;
				std::cout << std::endl;
				std::cout << std::endl;
				if( ROOTFileTest )
				{
					FitDetailsPath = argv[3];
					HistoDetails = argv[2];
					Fit fit( TimesPath, HistoDetails, FitDetailsPath );                                                      //Getting times from file to object
					fit.RangeBackgroundData();
					DecoOption = fit.Discrete_exp();
					if( DecoOption )
						fit.Deconvolution();
				}
				else
				{
					FitDetailsPath = argv[2];
					Fit fit( TimesPath, FitDetailsPath );                                                      //Getting times from file to object
					fit.RangeBackgroundData();
					DecoOption = fit.Discrete_exp();
					if( DecoOption )
						fit.Deconvolution();
				}
			}
		}
	}
	return 0;
}


