#include "PALS_avalanche_lib.h"

	
int main( int argc, char* argv[] ) 		//First argument - File with data, second - File with FitDetails(Path to it), third for more experienced users -> no argument Full Simplex, "old" old type, no simplex for Lifetime components, "exp" simplex for resolutino but predefined
{
    if( argc > 4 )                                                 //Argument check
	{
        std::cout << "Program is not using any more external arguments, so do not be officious!" << std::endl;
		return 0;
    }                                                               //Argument check							//Start counting time - control purpose
       
    std::string TimesPath = argv[1];                                    //Calculated times
	int DeconvolutionOption;
	int ROOTFileTest = 0;
	if( TimesPath.length() > 5 )
	{
		std::string TimesPathExtension;
		TimesPathExtension.insert(0, TimesPath, (TimesPath.length() - 5), TimesPath.length());
		if( TimesPathExtension == ".root" )
			ROOTFileTest = 1;
	}
	
    if( argc == 2 )
	{
		Fit fit( TimesPath, "FitDetails", ROOTFileTest );                                                      //Getting times from file to object
		fit.RangeBackgroundData();
		DeconvolutionOption = fit.Discrete();
		if( DeconvolutionOption )
	      		fit.Deconvolution();
	}
	else
	{
		std::string FitDetailsPath;
		std::string HistoDetails;
        FitDetailsPath = argv[2];
        std::string FitType = "";
        if( argc > 3 )
        {
            FitType = argv[3];
        }
        Fit fit( TimesPath, FitDetailsPath, ROOTFileTest, FitType );                                                      //Getting times from file to object
        fit.RangeBackgroundData();
		DeconvolutionOption = fit.Discrete();
		if( DeconvolutionOption )
	      		fit.Deconvolution();
	}
	return 0;
}


