#include "PALS_avalanche_lib.h"
#include <iomanip> 

LifetimeComponent::LifetimeComponent()			//Initializing default values for Lifetime Component 
{
	Type = "";
	Lifetime = 0.1;
	Intensity = 0.0001;
}

LifetimeComponent::LifetimeComponent( float lifetime, float intensity, std::string type )
{							//Initializing values for Lifetime Component given by user 
	Lifetime = lifetime;
	Intensity = intensity;
	Type = type;
}

DiscreteFitResult::DiscreteFitResult()			//Initializing default values for Fitted parameter 
{
	Parameter = 1;
	Uncertainity = 0.01;
	Type = "nf";	
}

DiscreteFitResult::DiscreteFitResult( double Par, double Unc )
{							//Initializing values for Fitted parameter 
	Parameter = Par;
	Uncertainity = Unc;
	Type = "nf";
}

Fit::Fit()
{
	std::cout << "Generating fit object" << std::endl;
}

Fit::Fit( std::string path, std::string pathForDetails )
{							//Getting Times and fitting details from file provided by user and FitDetails file
//-----------------------------------------------Clearing the variables used in the process of analysis of PAL spectrum
//-----------------------------------------------(Because reasons and to be sure that it is blank)
	Path = path;
	Times.clear();
	Arguments.clear();
	Values.clear();
//----------Details about histogram and the way of fitting
	Range_From = 0;
	Range_To = 0;
	Background = 0.0;
	SDBackground = 0.0;
	NmbOfBins = 0;
	NmbrOfIterations = 0;
	VarLvl = 0;
//----------Vectors used in deconvolution procedure which use the results from discrete fit
	LifetimesFromDiscrete.clear();
	ResolutionFromDiscrete.clear();
	OffsetsFromDiscrete.clear();
	FractionsFromDiscrete.clear();
	IntensitiesFromDiscrete.clear();
	LifetimeGrid.clear();
	GradientMatrix.clear();
	Residue.clear();
//----------Numbers transferred in the process of deconvolution - continous fitting
	iterator = 0;
	ChiDiffBetweenIterations = 0.;
//-----------------------------------------------
	char line[120];					      	//For lines with description and header
	char question;						//For user handling with overexaggeration with parameters
	std::string option = "", option2 = "", option3 = "";  	//Sometimes there are three arguments in line to get -> (Lifetime Intensity Type)
	std::cout << "-------------Processing fit details-------------" << std::endl;
	std::ifstream fitDetails;				//Stream for reading file
	std::cout << "-------------Checking the files-------------" << std::endl;
//-----------------------------------------------Checking if file with data provided by user and FitDetails file exists
//-----------------------------------------------(If one of them is not existing in provided path it end the work of these procedure)
//-----------------------------------------------(and in consequense the Resolution and Lifetimes vectors are empty -> rest of procedures do not work)	
	if( FileCheck( pathForDetails ) )
	{
		std::cout << "FitDetails file visible" << std::endl;
		if( FileCheck( path ) )
		{
			std::cout << path << " - file with data visible" << std::endl;
			std::cout << "Perfect, we can work like that" << std::endl;
		}
		else
		{
			std::cout << "!!! Missing file with data !!!" << std::endl;
			std::cout << "Fix this or the program have no chance to work" << std::endl;
			return;
		}
	}
	else
	{
		std::cout << "!!! Missing FitDetails file !!!" << std::endl;
		std::cout << "Fix this or the program have no chance to work" << std::endl;
		return;
	}
        fitDetails.open( pathForDetails );			// Opening the file with details of the analysis -> "FitDetails" - name of the file (Check)
//-----------------------------------------------
	std::cout << "-------------Reading the FitDetails-------------" << std::endl;
	for( unsigned i=0; i<3; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); //120 bo 120 to dwie kopy -> kopać można nogami > nogi są dwie || It is carefully taken value, chosen by optimizing one complicated chosing procedure
	}
	if( CheckTheLine( line[0], 'T', 3 ) == 0 )		//Test if Third line (the 3 in arguments is the number of line in order to write in terminal number of line that is bad) starts with 'T'
		return;
	fitDetails >> option;		//Getting the type of data
	if( option != "histogram" && option != "times" )	//Checking if type of data is provided correctly, one of the recognizable type - histogram or times
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Type of data is not one of correct types - 'histogram' or 'times' " << std::endl;
		return;
	}
	TypeOfData = option;
	for( unsigned i=0; i<3; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'W', 6 ) == 0 )		//Similar test as above line test for 6th line
		return;
	fitDetails >> option;		//Getting the Width of the Bin
	if( stod( option ) <= 0 || stod( option ) >= 142 )	//Checking if BinWidth is provided correctly - is should be given in ns, and maximum lifetime is 142 ns, BinWidth should not be higher than it or smaller than zero
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Bin width is not reasonable - less or equal zero or greater than 142 ns " << std::endl;
		return;	  
	}
	BinWidth = stod( option );
	for( unsigned i=0; i<3; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'I', 9 ) == 0 )		//Similar test as above line test for 9th line
		return;
	fitDetails >> option;		//Getting the extended type of data
	if( option != "oscilloscope" && option != "digitizer" && option != "other" )	//Checking if type of times is provided correctly, one of the recognizable type - oscilloscope or digitizer or other
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Type of times is not one of correct types - 'oscilloscope' or 'digitizer' or 'other' " << std::endl;
		return;
	}
	TypeOfDataExtended = option;
	for( unsigned i=0; i<4; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'L', 13 ) == 0 )		//Similar test as above line test for 13th line
		return;
	fitDetails >> option;		//Getting the value for searching last bin
	if( stoi( option ) < 0 )				//Checking if LastBin minimal value is provided correctly - is should be greater than zero from definition
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "LastBin minimal value is less than zero - is the histogram negative?" << std::endl;
		return;	  
	}
	else if( stoi( option ) > 100 )				//LastBin minimal value corresponds to the level of Background - greater than 100 means that statistics is very high, not usually spotted, this check is for user to confirm that user knows what is doing
	{
		std::cout << "Are You sure You have such high Background (More than 100)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line 14 for proper value of iterations" << std::endl;
			return;
		}
	}
	LastBinMinValue = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'F', 15 ) == 0 )		//Similar test as above line test for 15th line
		return;
	fitDetails >> option;		//Getting the value for searching first bin
	if( stoi( option ) < 0 )				//Similar as for LastBin minimal value
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "FirstBin minimal value is less than zero - is the histogram negative?" << std::endl;
		return;	  
	}
	else if( stoi( option ) > 100 )
	{
		std::cout << "Are You sure You have such high Background (More than 100)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line 16 for proper value of iterations" << std::endl;
			return;
		}
	}
	FirstBinMinValue = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'B', 17 ) == 0 )		//Similar test as above line test for 17th line
		return;
	fitDetails >> option;		//Getting the number of points to calculate background to
	if( stoi( option ) < 0 )				//Checking if Number of points for background calculations are given correctly - background can not be estimated from negative number of points
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Background can not be estimated from negative number of points" << std::endl;
		return;	  
	}
	else if( stoi( option ) > 1000 )			//Checking if user knows what is doing - background follows Poison distribution (in general), and it is correctly estimated if the number of pointsis in the order of 100, more than that is just the exaggeration
	{
		std::cout << "Are You sure You need to calculate background of that many of points(More than 1000)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line 18 for proper value of iterations" << std::endl;
			return;
		}
	}
	BackgroundEstimationNumbOfPoints = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'B', 19 ) == 0 )		//Similar test as above line test for 19th line
		return;	
	fitDetails >> option;		//Getting the side from which background is estimated
	if( option != "right" && option != "left" )		//Checking if side for background is provided correctly - it should be from one of the sides - left or right
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Side of histogram to estimate background should be one of the following - 'right' or 'left' " << std::endl;
		return;
	}
	SideOfBackgroundEstimation = option;
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'L', 21 ) == 0 )		//Similar test as above line test for 21th line
		return;	
	fitDetails >> option;		//Getting the number of points which is leaved from background calculation
	if( stoi( option ) < 0 )				//Checking if Shift for background is given properly - shift can not be lower than zero
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Program can not leave negative number of bins" << std::endl;
		return;
	}
	else if( stoi( option ) > 1000 )			//Checking if user knows what is doing - shifting more than 1000 points is the sign of bad shape of histogram - user must confirm that value is all right
	{
		std::cout << "Are You sure You need to leave that many of points(More than 1000)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line 22 for proper value of iterations" << std::endl;
			return;
		}
	}
	ShiftForBackgroundEstimation = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'E', 23 ) == 0 )		//Similar test as above line test for 23rd line
		return;		
	fitDetails >> option;		//Getting the End of the Fit value
	EndOfFitValue = stod( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'S', 25 ) == 0 )		//Similar test as above line test for 25th line
		return;	
	fitDetails >> option;		//Getting the Start of the Fit value
	if( stod( option ) < 0 || stod( option ) > 1 )		//Checking if StartBin for fit is given properly - it is given by fraction of maximum of histogram -> it can not be smaller than 0 (negative value) or grater than zero (value greater than maximum)
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Program can not start at negative point or higher than maximum" << std::endl;
		return;
	}	
	StartOfFitValue = stod( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'F', 27 ) == 0 )		//Similar test as above line test for 27th line
		return;	
	fitDetails >> option;		//Getting the Center for First Bin	
	FirstBinCenter = stod( option ) - BinWidth/2;
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'D', 29 ) == 0 )		//Similar test as above line test for 27th line
		return;	
	fitDetails >> option;		//Getting the Option for fixing resolution component	
	if( option != "yes" && option != "no" )		//Checking if Option for fixing is provided correctly - yes or no
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "FixGauss Option can be just one of the following - 'yes' or 'no' " << std::endl;
		return;
	}
	FixGauss = option;
	for( unsigned i=0; i<4; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'S', 33 ) == 0 )		//Similar test as above line test for 29th line
		return;	
//-----------------------------------------------
	fitDetails >> option >> option2;
	unsigned MaxIterator = 0;
	while( option != "-------" )
	{
		MaxIterator++;
		if( MaxIterator > 9 )				//If this loop works too long, there is a chance that FitDetails file is corrupted. In other hand, program is not prepared for fitting more than 10 Resolution components
		{
			std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
			std::cout << "Too much Resolution components (>9) or bad style of writing" << std::endl;
			std::cout << "It should be like: [Number] [tabulator between them] [Number]" << std::endl;
			std::cout << "and after writing all components ended by new line with -------	----------" << std::endl;
			return;
		}
		if( stod(option) > 0 && stod(option2) > 0 )	//Check if Resolution components are provided correctly
			Resolution.push_back( LifetimeComponent( stod(option), stod(option2), "gauss" ) );
		else
		{
			std::cout << "Intensity or Sigma can not be 0!!!!" << std::endl;
			std::cout << "Error while reading component!!!!" << std::endl;
		}
		fitDetails >> option >> option2;
	}
	MaxIterator = 0;
//-----------------------------------------------
	for( unsigned i=0; i<7; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'L', 40 + Resolution.size() ) == 0 )		//Similar test as above line test for 37th line
		return;	
	fitDetails >> option >> option2 >> option3;
	while( option != "-------" )
	{
		MaxIterator++;
		if( MaxIterator > 9 )				//Similar as for Resolution components
		{
			std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
			std::cout << "Too much Lifetime components (>9) or bad style of writing" << std::endl;
			std::cout << "It should be like: [Number] [tabulator between them] [Number] [tabulator between them] [Type]" << std::endl;
			std::cout << "and after writing all components ended by new line with -------	----- ----" << std::endl;
			return;
		}
		if( stod(option) > 0 && stod(option2) > 0 )
			Lifetimes.push_back( LifetimeComponent( stod(option), stod(option2), option3 ) );
		else
		{
			std::cout << "Intensity or Lifetime can not be 0!!!!" << std::endl;
			std::cout << "Error while reading component!!!!" << std::endl;
		}
		fitDetails >> option >> option2 >> option3;
	}
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'N', 42 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 43rd line
		return;		
	fitDetails >> option;
	if( stoi( option ) < 0 )				//Checking if Number of iteration are given properly, it can not be smaller than zero
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Number of iterations should be greater or equal zero" << std::endl;
		return;
	}
	else if( stoi( option ) > 20 )				//More iterations -> Longer fitting. Sometimes after 20 iterations parameters are not changing more, so it is check if user understands 
	{
		std::cout << "Are You sure You need that much iterations (More than 20)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line 44 for proper value of iterations" << std::endl;
			return;
		}
	}	  
	NmbrOfIterations = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'H', 44 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 45th line
		return;	
	fitDetails >> option;
	if( stod( option ) < 0 )				//Variability level should not be smaller than zero - it corrupts the boundaries of the fit parameters
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Level of variance should be greater or equal zero" << std::endl;
		return;
	}
	else if( stod( option ) > 100 )				//Variability level should not be greater than 1 -> It would provide a chance of fitting negative, unreal parameters
	{
		std::cout << "Are You sure You need that much variability (More than 100%)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line 46 for proper value of variability" << std::endl;
			return;
		}
	}
	VarLvl = stod( option )/100;
	for( unsigned i=0; i<5; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'D', 49 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 50th line
		return;
	fitDetails >> option;
	if( option != "yes" && option != "no" )		//Checking if option for deconvolution is defined one - yes or no
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Deconvolution Option can be just one of the following - 'yes' or 'no' " << std::endl;
		return;
	}
	else if( option == "yes" )
		DecoOption = 1;
	else
		DecoOption = 0;
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'H', 51 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 52th line
		return;
	fitDetails >> option;
	if( stod( option ) < 0 )				//Denisty parameter characterizes the distances between points on Lifetime Grid used in continous fitting - distance should not be smaller than 0
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Lifetime Grid density parameter should be greater than zero" << std::endl;
		return;
	}
	Scaler = stod( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'L', 53 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 54th line
		return;
	fitDetails >> option;
	if( stod( option ) < 0 )				//Minimal point of Lifetime Grid should be greater than zero, because the lifetime of positronium can not be smaller than zero
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Start point should not be negative fraction of Min Lifetime" << std::endl;
		return;
	}
	else if( stod( option ) > 100 )				//Minimal point of Lifetime Grid should be also smaller than minimal value from discrete fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Start point should not be greater than Min Lifetime" << std::endl;
		return;
	}
	FracMinLF = stod( option )/100.;
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'L', 55 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 55th line
		return;
	fitDetails >> option;
	if( stod( option ) < 1 )				//Maximal point should not be smaller than maximal value of lifetime fitted from discrete procedure
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "End point should not be lower than Max Lifetime" << std::endl;
		return;
	}
	FracMaxLF = stod( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'F', 57 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 58th line
		return;
	fitDetails >> option;
	MaxShiftForDeco = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'F', 59 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 58th line
		return;
	fitDetails >> option;
	if( stod( option ) < 0 )				//Number of point for linear smmothing of histogram, should be greater than zero -> zero if histogram is perfect
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Number of point for linear smoothing should be greater than zero" << std::endl;
		return;
	}
	LinFilterRange = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'R', 61 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 60th line
		return;
	fitDetails >> option;
	if( stod( option ) <= 1 )				//Multiplicity of Max Lifeitme for continous fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "End of Range for fitting should not be lower than maximal Lifetime from discrete fitting" << std::endl;
		return;
	}
	EndOfFitMultiplicity = stod( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'D', 63 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 62th line
		return;
	fitDetails >> option;
	if( stod( option ) <= 0 )				//Fraction of lifetime set to be default sigmas for lifetime distribution in continous fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Fraction of Lifetime for sigmas should not be less or equal than zero" << std::endl;
		return;
	}
	SigmasDefaultFraction = stod( option )/100.;
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'S', 65 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 64th line
		return;
	fitDetails >> option;
	if( stod( option ) <= 0 )				//Step that controls how big steps will be done in preiteration procedure - bigger step quicker fitting, smaller more precise results
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Start value of step can not be lower or equal 0" << std::endl;
		return;
	}
	StepForPreIteration = stod( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'N', 67 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 66th line
		return;
	fitDetails >> option;
	if( stod( option ) < 0 )				//Number of iteration to check the condition to end preiteration -> higher value -> precise estimation -> longer fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Number of preiterations can not be lower than 0" << std::endl;
		return;
	}
	NumberOfPreIterations = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'S', 69 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 68th line
		return;
	fitDetails >> option;
	if( stod( option ) <= 0 )				//Step that controls how big steps will be done in iteration procedure - bigger step quicker fitting, smaller more precise results
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Start value of step can not be lower or equal 0" << std::endl;
		return;
	}
	StepForIteration = stod( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'N', 71 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 70th line
		return;
	fitDetails >> option;
	if( stod( option ) < 0 )				//Number of iteration to check the condition to end iteration -> higher value -> precise estimation -> longer fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Number of iterations can not be lower than 0" << std::endl;
		return;
	}
	NumberOfIterations = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'S', 73 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 70th line
		return;
	fitDetails >> option;
	if( option != "Gaussian" && option != "LogGaussian" && option != "Mixed" )				//Number of iteration to check the condition to end iteration -> higher value -> precise estimation -> longer fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Shape must be one of the mentioned - Gaussian, LogGaussian or Mixed" << std::endl;
		return;
	}
	ShapeOfComponent = option;
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'D', 75 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 50th line
		return;
	fitDetails >> option;
	if( option != "yes" && option != "no" )		//Checking if option for deconvolution is defined one - yes or no
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Type of Continous Fit Option can be just one of the following - 'yes' or 'no' " << std::endl;
		return;
	}
	else if( option == "yes" )
		TypeOfContinousFit = 1;
	else
		TypeOfContinousFit = 0;
	std::cout << "Number of Resolution Components " << Resolution.size() << std::endl;
	std::cout << "Number of Lifetimes Components " << Lifetimes.size() << std::endl;
//-----------------------------------------------It is better to have intensites that sum up to 1 -> If it is not summing up, it can cause possible error in fitting
	std::cout << "-------------Checking the correctness of parameters (normalization)-------------" << std::endl;
	double TotalIntensity = 0., FixedIntensity = 0., TotalFreeIntensity = 0;		
//-----------------------------------------------Setting three counters for summing the intensities
//-----------------------------------------------TotalIntensity -> Intensity of free components and partially fixed (not "f" and "ps" type)
//-----------------------------------------------FixedIntensity -> Intensity of fixed ("f" type)
//-----------------------------------------------TotalFreeIntensity -> Intensity of only free components ("nf" type) (and also resolution components)
	if( Resolution.size()*Lifetimes.size() == 0 ) 			//If there was error in reading the details of resolution or lifetime components, the program should remove those uncertain components
	{								//If there was no good resolution or lifetimes components program will stop the work
		std::cout << std::endl;
		std::cout << "Try in the future fill properly parameters line" << std::endl;
		std::cout << "Both, number of Resolution and Lifetime components" << std::endl;
		std::cout << "have to be greater than 0" << std::endl;
		std::cout << "For now the space for the details for them are empty" << std::endl;
		std::cout << "Or they are badly provided not like they should be" << std::endl;
		std::cout << std::endl;
		return;
	}
	for( unsigned i=0; i<Resolution.size(); i++ )			//Counting total intensities of resolution components that were read from FitDetails file -> intensities of resolution components will be then
	{								//given intensity that was read divided by TotalIntenisty -> In that way sum of intenisities should be equal to 1 (normalized)
		TotalIntensity += Resolution[i].Intensity;
	}
	if( fabs(TotalIntensity - 1) > 0.01 )				//Storing the double is ... not fully reliable -> accepting intensities if TotalIntenisty is more or less 1 ( 0.99 - 1.01 )
	{
		std::cout << "-------------Renormalization Resolution intensities-------------" << std::endl;
		std::cout << "Try in the future write intensities such, that they sum up to 1" << std::endl;
		std::cout << "It would spare the time of renormalization intensities," << std::endl;
		std::cout << "and minimize the chance of observing bug in the code" << std::endl;
		std::cout << std::endl;
		for( unsigned i=0; i<Resolution.size(); i++ )
		{
			if( TotalIntensity == 0 )			//This part is here to prevent some calculations problems
			{
				std::cout << "Total intensity is equal to 0" << std::endl;
				std::cout << "Setting all resolution intensities as 1/number of components" << std::endl;
				std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
				Resolution[i].Intensity = 1/Resolution.size();
			}
			Resolution[i].Intensity = Resolution[i].Intensity/TotalIntensity;
		}
	}
	TotalIntensity = 0.;						//Normalization of intensities for Lifetime Components, similar as for Resolution components
	for( unsigned i=0; i<Lifetimes.size(); i++ )
	{
		if( Lifetimes[i].Type == "nf" )
		{
			LifetimesNotFixed.push_back( Lifetimes[i] );
			TotalFreeIntensity += Lifetimes[i].Intensity;
		}
		if( Lifetimes[i].Type == "f" )
			FixedIntensity += Lifetimes[i].Intensity;
		else if( Lifetimes[i].Type != "ps" )
			TotalIntensity += Lifetimes[i].Intensity;
	}	
	if( fabs(TotalIntensity + FixedIntensity - 1) > 0.01 )		//The same as for resolution 
	{
		std::cout << std::endl;
		std::cout << "-------------Renormalization Lifetime intensities-------------" << std::endl;
		std::cout << "Try in the future write intensities such, that they sum up to 1" << std::endl;
		std::cout << "It would spare the time of renormalization intensities," << std::endl;
		std::cout << "and minimize the chance of observing bug in the code" << std::endl;
		std::cout << std::endl;
		for( unsigned i=0; i<Lifetimes.size(); i++ )
		{
			if( (FixedIntensity >= 1 || ( TotalIntensity == 0 && FixedIntensity > 0 ) ) && Lifetimes[i].Type != "ps" )
			{						//Intensities for fixed components are treated differently than for free components
				Lifetimes[i].Intensity = Lifetimes[i].Intensity/(FixedIntensity + TotalIntensity);
			}
			else if( Lifetimes[i].Type != "ps" )
			{
				if( TotalIntensity == 0 && FixedIntensity == 0 )	//This part is here to prevent some calculations problems
				{
					std::cout << "Total intensity is equal to 0" << std::endl;
					std::cout << "Setting all lifetime intensities as 1/number of components" << std::endl;
					std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
					Lifetimes[i].Intensity = 1/Lifetimes.size();
				}
				else if( Lifetimes[i].Type != "f" )
					Lifetimes[i].Intensity = Lifetimes[i].Intensity*(1-FixedIntensity)/TotalIntensity;
			}
		}
	}
	for( unsigned i=0; i<LifetimesNotFixed.size(); i++ )
	{
		LifetimesNotFixed[i].Intensity = LifetimesNotFixed[i].Intensity/TotalFreeIntensity;
	}
//-----------------------------------------------Sorting Lifetime components such: ps first -> then somehow fixed components -> free components 
//-----------------------------------------------(Because knowledge of amount of fixed components is used in body of fitting function and it is lowering the time)
	SortLifetimesByType();
//-----------------------------------------------
	std::cout << "-------------Getting data-------------" << std::endl;
	std::ifstream data;		//Stream for path to data
	std::string count = "";		//For oscilloscope
	std::string time = "";		//Proper string
	std::string risetime1 = "";	//For oscilloscope
	std::string risetime2 = "";	//For oscilloscope
        data.open( ( Path ).c_str() );

	if( TypeOfData == "histogram" )
	{
	      while( data >> time )	//other
	      {
		      if( time != "inf" && time != "nan"){
			  Times.push_back( stod(time) );
			  NmbOfBins++;
		      }
	      }
	}
	else
	{  
		if( TypeOfDataExtended == "oscilloscope" )
		{
	     		 while( data >> count >> time >> risetime1 >> risetime2 )
	      		{
		      		if( time != "inf" && time != "nan"){
			 		 Times.push_back( stod(time) );
		    		  }		
	     		 }
		}
		else			//digitizer
		{
			while( data >> time )
             		 {
                     		 if( time != "inf" && time != "nan"){
                         		 Times.push_back( stod(time) );
                     		 }
             		 }
		}
	}
	std::cout << "Number of Bins " << NmbOfBins << std::endl;
	if( NmbOfBins )						//Checking if provided by user details about the fit is correct
	{							//NmbOfBins must be greater than BackgroundEstimationNumbOfPoints and ShiftForBackgroundEstimation
		if( NmbOfBins < BackgroundEstimationNumbOfPoints + ShiftForBackgroundEstimation )
		{
			std::cout << "Number Of Bins can not be smaller than number of points for background calculations" << std::endl;
			std::cout << "or Number fo points left in background calculations or their sum" << std::endl;
			Resolution.clear();
			Lifetimes.clear();
		}
		else if( NmbOfBins*BinWidth < EndOfFitValue )
		{
			std::cout << "Number of Points * BinWidth " << NmbOfBins*BinWidth << std::endl;
			std::cout << "End of the range of the fit should not be greater than maximal value of histogram " << std::endl;
			Resolution.clear();
			Lifetimes.clear();
		}
	}
	data.close();
}

Fit::Fit( std::string path, std::string nameOfHistogram, std::string pathForDetails )
{							//Getting Times and fitting details from file provided by user and FitDetails file
//-----------------------------------------------Clearing the variables used in the process of analysis of PAL spectrum
//-----------------------------------------------(Because reasons and to be sure that it is blank)
	Path = path;
	Times.clear();
	Arguments.clear();
	Values.clear();
//----------Details about histogram and the way of fitting
	Range_From = 0;
	Range_To = 0;
	Background = 0.0;
	SDBackground = 0.0;
	NmbOfBins = 0;
	NmbrOfIterations = 0;
	VarLvl = 0;
//----------Vectors used in deconvolution procedure which use the results from discrete fit
	LifetimesFromDiscrete.clear();
	ResolutionFromDiscrete.clear();
	OffsetsFromDiscrete.clear();
	FractionsFromDiscrete.clear();
	IntensitiesFromDiscrete.clear();
	LifetimeGrid.clear();
	GradientMatrix.clear();
	Residue.clear();
//----------Numbers transferred in the process of deconvolution - continous fitting
	iterator = 0;
	ChiDiffBetweenIterations = 0.;
//-----------------------------------------------
	char line[120];					      	//For lines with description and header
	char question;						//For user handling with overexaggeration with parameters
	std::string option = "", option2 = "", option3 = "";  	//Sometimes there are three arguments in line to get -> (Lifetime Intensity Type)
	std::cout << "-------------Processing fit details-------------" << std::endl;
	std::ifstream fitDetails;				//Stream for reading file
	std::cout << "-------------Checking the files-------------" << std::endl;
//-----------------------------------------------Checking if file with data provided by user and FitDetails file exists
//-----------------------------------------------(If one of them is not existing in provided path it end the work of these procedure)
//-----------------------------------------------(and in consequense the Resolution and Lifetimes vectors are empty -> rest of procedures do not work)	
	if( FileCheck( pathForDetails ) )
	{
		std::cout << "FitDetails file visible" << std::endl;
		if( FileCheck( path ) )
		{
			std::cout << path << " - file with data visible" << std::endl;
			std::cout << "Perfect, we can work like that" << std::endl;
		}
		else
		{
			std::cout << "!!! Missing file with data !!!" << std::endl;
			std::cout << "Fix this or the program have no chance to work" << std::endl;
			return;
		}
	}
	else
	{
		std::cout << "!!! Missing FitDetails file !!!" << std::endl;
		std::cout << "Fix this or the program have no chance to work" << std::endl;
		return;
	}
        fitDetails.open( pathForDetails );			// Opening the file with details of the analysis -> "FitDetails" - name of the file (Check)
//-----------------------------------------------
	std::cout << "-------------Reading the FitDetails-------------" << std::endl;
	for( unsigned i=0; i<3; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); //120 bo 120 to dwie kopy -> kopać można nogami > nogi są dwie || It is carefully taken value, chosen by optimizing one complicated chosing procedure
	}
	if( CheckTheLine( line[0], 'T', 3 ) == 0 )		//Test if Third line (the 3 in arguments is the number of line in order to write in terminal number of line that is bad) starts with 'T'
		return;
	fitDetails >> option;		//Getting the type of data
	if( option != "histogram" && option != "times" )	//Checking if type of data is provided correctly, one of the recognizable type - histogram or times
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Type of data is not one of correct types - 'histogram' or 'times' " << std::endl;
		return;
	}
	TypeOfData = option;
	for( unsigned i=0; i<3; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'W', 6 ) == 0 )		//Similar test as above line test for 6th line
		return;
	fitDetails >> option;		//Getting the Width of the Bin
	if( stod( option ) <= 0 || stod( option ) >= 142 )	//Checking if BinWidth is provided correctly - is should be given in ns, and maximum lifetime is 142 ns, BinWidth should not be higher than it or smaller than zero
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Bin width is not reasonable - less or equal zero or greater than 142 ns " << std::endl;
		return;	  
	}
	BinWidth = stod( option );
	for( unsigned i=0; i<3; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'I', 9 ) == 0 )		//Similar test as above line test for 9th line
		return;
	fitDetails >> option;		//Getting the extended type of data
	if( option != "oscilloscope" && option != "digitizer" && option != "other" )	//Checking if type of times is provided correctly, one of the recognizable type - oscilloscope or digitizer or other
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Type of times is not one of correct types - 'oscilloscope' or 'digitizer' or 'other' " << std::endl;
		return;
	}
	TypeOfDataExtended = option;
	for( unsigned i=0; i<4; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'L', 13 ) == 0 )		//Similar test as above line test for 13th line
		return;
	fitDetails >> option;		//Getting the value for searching last bin
	if( stoi( option ) < 0 )				//Checking if LastBin minimal value is provided correctly - is should be greater than zero from definition
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "LastBin minimal value is less than zero - is the histogram negative?" << std::endl;
		return;	  
	}
	else if( stoi( option ) > 100 )				//LastBin minimal value corresponds to the level of Background - greater than 100 means that statistics is very high, not usually spotted, this check is for user to confirm that user knows what is doing
	{
		std::cout << "Are You sure You have such high Background (More than 100)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line 14 for proper value of iterations" << std::endl;
			return;
		}
	}
	LastBinMinValue = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'F', 15 ) == 0 )		//Similar test as above line test for 15th line
		return;
	fitDetails >> option;		//Getting the value for searching first bin
	if( stoi( option ) < 0 )				//Similar as for LastBin minimal value
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "FirstBin minimal value is less than zero - is the histogram negative?" << std::endl;
		return;	  
	}
	else if( stoi( option ) > 100 )
	{
		std::cout << "Are You sure You have such high Background (More than 100)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line 16 for proper value of iterations" << std::endl;
			return;
		}
	}
	FirstBinMinValue = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'B', 17 ) == 0 )		//Similar test as above line test for 17th line
		return;
	fitDetails >> option;		//Getting the number of points to calculate background to
	if( stoi( option ) < 0 )				//Checking if Number of points for background calculations are given correctly - background can not be estimated from negative number of points
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Background can not be estimated from negative number of points" << std::endl;
		return;	  
	}
	else if( stoi( option ) > 1000 )			//Checking if user knows what is doing - background follows Poison distribution (in general), and it is correctly estimated if the number of pointsis in the order of 100, more than that is just the exaggeration
	{
		std::cout << "Are You sure You need to calculate background of that many of points(More than 1000)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line 18 for proper value of iterations" << std::endl;
			return;
		}
	}
	BackgroundEstimationNumbOfPoints = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'B', 19 ) == 0 )		//Similar test as above line test for 19th line
		return;	
	fitDetails >> option;		//Getting the side from which background is estimated
	if( option != "right" && option != "left" )		//Checking if side for background is provided correctly - it should be from one of the sides - left or right
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Side of histogram to estimate background should be one of the following - 'right' or 'left' " << std::endl;
		return;
	}
	SideOfBackgroundEstimation = option;
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'L', 21 ) == 0 )		//Similar test as above line test for 21th line
		return;	
	fitDetails >> option;		//Getting the number of points which is leaved from background calculation
	if( stoi( option ) < 0 )				//Checking if Shift for background is given properly - shift can not be lower than zero
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Program can not leave negative number of bins" << std::endl;
		return;
	}
	else if( stoi( option ) > 1000 )			//Checking if user knows what is doing - shifting more than 1000 points is the sign of bad shape of histogram - user must confirm that value is all right
	{
		std::cout << "Are You sure You need to leave that many of points(More than 1000)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line 22 for proper value of iterations" << std::endl;
			return;
		}
	}
	ShiftForBackgroundEstimation = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'E', 23 ) == 0 )		//Similar test as above line test for 23rd line
		return;		
	fitDetails >> option;		//Getting the End of the Fit value
	EndOfFitValue = stod( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'S', 25 ) == 0 )		//Similar test as above line test for 25th line
		return;	
	fitDetails >> option;		//Getting the Start of the Fit value
	if( stod( option ) < 0 || stod( option ) > 1 )		//Checking if StartBin for fit is given properly - it is given by fraction of maximum of histogram -> it can not be smaller than 0 (negative value) or grater than zero (value greater than maximum)
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Program can not start at negative point or higher than maximum" << std::endl;
		return;
	}	
	StartOfFitValue = stod( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'F', 27 ) == 0 )		//Similar test as above line test for 27th line
		return;	
	fitDetails >> option;		//Getting the Center for First Bin	
	FirstBinCenter = stod( option ) - BinWidth/2;
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'D', 29 ) == 0 )		//Similar test as above line test for 27th line
		return;	
	fitDetails >> option;		//Getting the Option for fixing resolution component	
	if( option != "yes" && option != "no" )		//Checking if Option for fixing is provided correctly - yes or no
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "FixGauss Option can be just one of the following - 'yes' or 'no' " << std::endl;
		return;
	}
	FixGauss = option;
	for( unsigned i=0; i<4; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'S', 33 ) == 0 )		//Similar test as above line test for 29th line
		return;	
//-----------------------------------------------
	fitDetails >> option >> option2;
	unsigned MaxIterator = 0;
	while( option != "-------" )
	{
		MaxIterator++;
		if( MaxIterator > 9 )				//If this loop works too long, there is a chance that FitDetails file is corrupted. In other hand, program is not prepared for fitting more than 10 Resolution components
		{
			std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
			std::cout << "Too much Resolution components (>9) or bad style of writing" << std::endl;
			std::cout << "It should be like: [Number] [tabulator between them] [Number]" << std::endl;
			std::cout << "and after writing all components ended by new line with -------	----------" << std::endl;
			return;
		}
		if( stod(option) > 0 && stod(option2) > 0 )	//Check if Resolution components are provided correctly
			Resolution.push_back( LifetimeComponent( stod(option), stod(option2), "gauss" ) );
		else
		{
			std::cout << "Intensity or Sigma can not be 0!!!!" << std::endl;
			std::cout << "Error while reading component!!!!" << std::endl;
		}
		fitDetails >> option >> option2;
	}
	MaxIterator = 0;
//-----------------------------------------------
	for( unsigned i=0; i<7; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'L', 40 + Resolution.size() ) == 0 )		//Similar test as above line test for 37th line
		return;	
	fitDetails >> option >> option2 >> option3;
	while( option != "-------" )
	{
		MaxIterator++;
		if( MaxIterator > 9 )				//Similar as for Resolution components
		{
			std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
			std::cout << "Too much Lifetime components (>9) or bad style of writing" << std::endl;
			std::cout << "It should be like: [Number] [tabulator between them] [Number] [tabulator between them] [Type]" << std::endl;
			std::cout << "and after writing all components ended by new line with -------	----- ----" << std::endl;
			return;
		}
		if( stod(option) > 0 && stod(option2) > 0 )
			Lifetimes.push_back( LifetimeComponent( stod(option), stod(option2), option3 ) );
		else
		{
			std::cout << "Intensity or Lifetime can not be 0!!!!" << std::endl;
			std::cout << "Error while reading component!!!!" << std::endl;
		}
		fitDetails >> option >> option2 >> option3;
	}
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'N', 42 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 43rd line
		return;		
	fitDetails >> option;
	if( stoi( option ) < 0 )				//Checking if Number of iteration are given properly, it can not be smaller than zero
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Number of iterations should be greater or equal zero" << std::endl;
		return;
	}
	else if( stoi( option ) > 20 )				//More iterations -> Longer fitting. Sometimes after 20 iterations parameters are not changing more, so it is check if user understands 
	{
		std::cout << "Are You sure You need that much iterations (More than 20)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line 44 for proper value of iterations" << std::endl;
			return;
		}
	}	  
	NmbrOfIterations = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'H', 44 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 45th line
		return;	
	fitDetails >> option;
	if( stod( option ) < 0 )				//Variability level should not be smaller than zero - it corrupts the boundaries of the fit parameters
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Level of variance should be greater or equal zero" << std::endl;
		return;
	}
	else if( stod( option ) > 100 )				//Variability level should not be greater than 1 -> It would provide a chance of fitting negative, unreal parameters
	{
		std::cout << "Are You sure You need that much variability (More than 100%)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line 46 for proper value of variability" << std::endl;
			return;
		}
	}
	VarLvl = stod( option )/100;
	for( unsigned i=0; i<5; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'D', 49 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 50th line
		return;
	fitDetails >> option;
	if( option != "yes" && option != "no" )		//Checking if option for deconvolution is defined one - yes or no
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Deconvolution Option can be just one of the following - 'yes' or 'no' " << std::endl;
		return;
	}
	else if( option == "yes" )
		DecoOption = 1;
	else
		DecoOption = 0;
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'H', 51 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 52th line
		return;
	fitDetails >> option;
	if( stod( option ) < 0 )				//Denisty parameter characterizes the distances between points on Lifetime Grid used in continous fitting - distance should not be smaller than 0
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Lifetime Grid density parameter should be greater than zero" << std::endl;
		return;
	}
	Scaler = stod( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'L', 53 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 54th line
		return;
	fitDetails >> option;
	if( stod( option ) < 0 )				//Minimal point of Lifetime Grid should be greater than zero, because the lifetime of positronium can not be smaller than zero
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Start point should not be negative fraction of Min Lifetime" << std::endl;
		return;
	}
	else if( stod( option ) > 100 )				//Minimal point of Lifetime Grid should be also smaller than minimal value from discrete fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Start point should not be greater than Min Lifetime" << std::endl;
		return;
	}
	FracMinLF = stod( option )/100.;
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'L', 55 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 55th line
		return;
	fitDetails >> option;
	if( stod( option ) < 1 )				//Maximal point should not be smaller than maximal value of lifetime fitted from discrete procedure
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "End point should not be lower than Max Lifetime" << std::endl;
		return;
	}
	FracMaxLF = stod( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'F', 57 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 58th line
		return;
	fitDetails >> option;
	MaxShiftForDeco = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'F', 59 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 58th line
		return;
	fitDetails >> option;
	if( stod( option ) < 0 )				//Number of point for linear smmothing of histogram, should be greater than zero -> zero if histogram is perfect
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Number of point for linear smoothing should be greater than zero" << std::endl;
		return;
	}
	LinFilterRange = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'R', 61 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 60th line
		return;
	fitDetails >> option;
	if( stod( option ) <= 1 )				//Multiplicity of Max Lifeitme for continous fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "End of Range for fitting should not be lower than maximal Lifetime from discrete fitting" << std::endl;
		return;
	}
	EndOfFitMultiplicity = stod( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'D', 63 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 62th line
		return;
	fitDetails >> option;
	if( stod( option ) <= 0 )				//Fraction of lifetime set to be default sigmas for lifetime distribution in continous fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Fraction of Lifetime for sigmas should not be less or equal than zero" << std::endl;
		return;
	}
	SigmasDefaultFraction = stod( option )/100.;
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'S', 65 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 64th line
		return;
	fitDetails >> option;
	if( stod( option ) <= 0 )				//Step that controls how big steps will be done in preiteration procedure - bigger step quicker fitting, smaller more precise results
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Start value of step can not be lower or equal 0" << std::endl;
		return;
	}
	StepForPreIteration = stod( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'N', 67 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 66th line
		return;
	fitDetails >> option;
	if( stod( option ) < 0 )				//Number of iteration to check the condition to end preiteration -> higher value -> precise estimation -> longer fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Number of preiterations can not be lower than 0" << std::endl;
		return;
	}
	NumberOfPreIterations = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'S', 69 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 68th line
		return;
	fitDetails >> option;
	if( stod( option ) <= 0 )				//Step that controls how big steps will be done in iteration procedure - bigger step quicker fitting, smaller more precise results
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Start value of step can not be lower or equal 0" << std::endl;
		return;
	}
	StepForIteration = stod( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'N', 71 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 70th line
		return;
	fitDetails >> option;
	if( stod( option ) < 0 )				//Number of iteration to check the condition to end iteration -> higher value -> precise estimation -> longer fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Number of iterations can not be lower than 0" << std::endl;
		return;
	}
	NumberOfIterations = stoi( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'S', 73 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 70th line
		return;
	fitDetails >> option;
	if( option != "Gaussian" && option != "LogGaussian" && option != "Mixed" )				//Number of iteration to check the condition to end iteration -> higher value -> precise estimation -> longer fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Shape must be one of the mentioned - Gaussian, LogGaussian or Mixed" << std::endl;
		return;
	}
	ShapeOfComponent = option;
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line,120 ); 
	}
	if( CheckTheLine( line[0], 'D', 75 + Resolution.size() + Lifetimes.size() ) == 0 )		//Similar test as above line test for 50th line
		return;
	fitDetails >> option;
	if( option != "yes" && option != "no" )		//Checking if option for deconvolution is defined one - yes or no
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Type of Continous Fit Option can be just one of the following - 'yes' or 'no' " << std::endl;
		return;
	}
	else if( option == "yes" )
		TypeOfContinousFit = 1;
	else
		TypeOfContinousFit = 0;
	std::cout << "Number of Resolution Components " << Resolution.size() << std::endl;
	std::cout << "Number of Lifetimes Components " << Lifetimes.size() << std::endl;
//-----------------------------------------------It is better to have intensites that sum up to 1 -> If it is not summing up, it can cause possible error in fitting
	std::cout << "-------------Checking the correctness of parameters (normalization)-------------" << std::endl;
	double TotalIntensity = 0., FixedIntensity = 0., TotalFreeIntensity = 0;		
//-----------------------------------------------Setting three counters for summing the intensities
//-----------------------------------------------TotalIntensity -> Intensity of free components and partially fixed (not "f" and "ps" type)
//-----------------------------------------------FixedIntensity -> Intensity of fixed ("f" type)
//-----------------------------------------------TotalFreeIntensity -> Intensity of only free components ("nf" type) (and also resolution components)
	if( Resolution.size()*Lifetimes.size() == 0 ) 			//If there was error in reading the details of resolution or lifetime components, the program should remove those uncertain components
	{								//If there was no good resolution or lifetimes components program will stop the work
		std::cout << std::endl;
		std::cout << "Try in the future fill properly parameters line" << std::endl;
		std::cout << "Both, number of Resolution and Lifetime components" << std::endl;
		std::cout << "have to be greater than 0" << std::endl;
		std::cout << "For now the space for the details for them are empty" << std::endl;
		std::cout << "Or they are badly provided not like they should be" << std::endl;
		std::cout << std::endl;
		return;
	}
	for( unsigned i=0; i<Resolution.size(); i++ )			//Counting total intensities of resolution components that were read from FitDetails file -> intensities of resolution components will be then
	{								//given intensity that was read divided by TotalIntenisty -> In that way sum of intenisities should be equal to 1 (normalized)
		TotalIntensity += Resolution[i].Intensity;
	}
	if( fabs(TotalIntensity - 1) > 0.01 )				//Storing the double is ... not fully reliable -> accepting intensities if TotalIntenisty is more or less 1 ( 0.99 - 1.01 )
	{
		std::cout << "-------------Renormalization Resolution intensities-------------" << std::endl;
		std::cout << "Try in the future write intensities such, that they sum up to 1" << std::endl;
		std::cout << "It would spare the time of renormalization intensities," << std::endl;
		std::cout << "and minimize the chance of observing bug in the code" << std::endl;
		std::cout << std::endl;
		for( unsigned i=0; i<Resolution.size(); i++ )
		{
			if( TotalIntensity == 0 )			//This part is here to prevent some calculations problems
			{
				std::cout << "Total intensity is equal to 0" << std::endl;
				std::cout << "Setting all resolution intensities as 1/number of components" << std::endl;
				std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
				Resolution[i].Intensity = 1/Resolution.size();
			}
			Resolution[i].Intensity = Resolution[i].Intensity/TotalIntensity;
		}
	}
	TotalIntensity = 0.;						//Normalization of intensities for Lifetime Components, similar as for Resolution components
	for( unsigned i=0; i<Lifetimes.size(); i++ )
	{
		if( Lifetimes[i].Type == "nf" )
		{
			LifetimesNotFixed.push_back( Lifetimes[i] );
			TotalFreeIntensity += Lifetimes[i].Intensity;
		}
		if( Lifetimes[i].Type == "f" )
			FixedIntensity += Lifetimes[i].Intensity;
		else if( Lifetimes[i].Type != "ps" )
			TotalIntensity += Lifetimes[i].Intensity;
	}	
	if( fabs(TotalIntensity + FixedIntensity - 1) > 0.01 )		//The same as for resolution 
	{
		std::cout << std::endl;
		std::cout << "-------------Renormalization Lifetime intensities-------------" << std::endl;
		std::cout << "Try in the future write intensities such, that they sum up to 1" << std::endl;
		std::cout << "It would spare the time of renormalization intensities," << std::endl;
		std::cout << "and minimize the chance of observing bug in the code" << std::endl;
		std::cout << std::endl;
		for( unsigned i=0; i<Lifetimes.size(); i++ )
		{
			if( (FixedIntensity >= 1 || ( TotalIntensity == 0 && FixedIntensity > 0 ) ) && Lifetimes[i].Type != "ps" )
			{						//Intensities for fixed components are treated differently than for free components
				Lifetimes[i].Intensity = Lifetimes[i].Intensity/(FixedIntensity + TotalIntensity);
			}
			else if( Lifetimes[i].Type != "ps" )
			{
				if( TotalIntensity == 0 && FixedIntensity == 0 )	//This part is here to prevent some calculations problems
				{
					std::cout << "Total intensity is equal to 0" << std::endl;
					std::cout << "Setting all lifetime intensities as 1/number of components" << std::endl;
					std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
					Lifetimes[i].Intensity = 1/Lifetimes.size();
				}
				else if( Lifetimes[i].Type != "f" )
					Lifetimes[i].Intensity = Lifetimes[i].Intensity*(1-FixedIntensity)/TotalIntensity;
			}
		}
	}
	for( unsigned i=0; i<LifetimesNotFixed.size(); i++ )
	{
		LifetimesNotFixed[i].Intensity = LifetimesNotFixed[i].Intensity/TotalFreeIntensity;
	}
//-----------------------------------------------Sorting Lifetime components such: ps first -> then somehow fixed components -> free components 
//-----------------------------------------------(Because knowledge of amount of fixed components is used in body of fitting function and it is lowering the time)
	SortLifetimesByType();
//-----------------------------------------------
	std::cout << "-------------Getting data-------------" << std::endl;
	std::ifstream data;		//Stream for path to data
	std::string count = "";		//For oscilloscope
	std::string time = "";		//Proper string
	std::string risetime1 = "";	//For oscilloscope
	std::string risetime2 = "";	//For oscilloscope

	TFile* file1 = new TFile( path.c_str(), "READ" );
	TDirectory *dir;
	file1->GetObject("EventCategorizer subtask 0 stats",dir);
	TH1D* histo1 = (TH1D*) dir->Get( nameOfHistogram.c_str() );
	
	int RebinNumber = 1;
	std::cout << "Rebinning of the histogram. Type 1 if You do not want to rebin, type integer number greater than 1 if You want to rebin" << std::endl;
	std::cout << "Meaning of the value that should be typed -> type 10 if You want to merge every 10 bins." << std::endl;
	std::cin >> RebinNumber;
	histo1 -> Rebin(RebinNumber);
	
	for( int i=0; i<histo1->GetXaxis()->GetNbins(); i++ )
        {
		Times.push_back( histo1->GetBinContent(i) );
		NmbOfBins++;
        }
        BinWidth = histo1->GetBinCenter(3) - histo1->GetBinCenter(2);
	FirstBinCenter = histo1->GetBinCenter(0);
	
	file1->Close();
	std::cout << "Number of Bins " << NmbOfBins << std::endl;
	if( NmbOfBins )						//Checking if provided by user details about the fit is correct
	{							//NmbOfBins must be greater than BackgroundEstimationNumbOfPoints and ShiftForBackgroundEstimation
		if( NmbOfBins < BackgroundEstimationNumbOfPoints + ShiftForBackgroundEstimation )
		{
			std::cout << "Number Of Bins can not be smaller than number of points for background calculations" << std::endl;
			std::cout << "or Number fo points left in background calculations or their sum" << std::endl;
			Resolution.clear();
			Lifetimes.clear();
		}
		else if( NmbOfBins*BinWidth < EndOfFitValue )
		{
			std::cout << "Number of Points * BinWidth " << NmbOfBins*BinWidth << std::endl;
			std::cout << "End of the range of the fit should not be greater than maximal value of histogram " << std::endl;
			Resolution.clear();
			Lifetimes.clear();
		}
	}
	data.close();
}

void Fit::RangeBackgroundData()
{
	if( Resolution.size()*Lifetimes.size() == 0 )		//Test if the reading FitDetails file and data file past properly
	{
		return;
	}
        std::cout<<"-------------Preparing data-------------"<<std::endl;
	TH1F* histogram;		//Histogram to fit
	if( ! NmbOfBins )  // If the data is from oscilloscope or digitizer NmbOfBins is equal to 0, if one would like to change range of histogram change te second, third, and fourth argument in Fill Histogram below
	{
		histogram = FillHistogram( "Time differences", 120/BinWidth, FirstBinCenter, 120 + FirstBinCenter, Times );
	}
	else			//If the data type is histogram NmbOfBins > 0
	{
		histogram = FillHistogram( "Time differences", NmbOfBins, FirstBinCenter, NmbOfBins*BinWidth + FirstBinCenter, Times );
		for(unsigned i=0; i<NmbOfBins; i++)		//Setting the histogram such that each input line in source file is the content of bin
		{
			histogram->SetBinContent(i,Times[i]);
		}	  
	}	
//--------------------------------------- Looking For the first and the last nonzero bin of the histogram 
	BinMax = histogram -> GetMaximumBin();					//bin with maximal content -> As the start point for looking for range of fit
	LastBin = histogram -> FindLastBinAbove( LastBinMinValue, 1 ); 		//last bin above LastBinMinValue on 1 -> (x axis) - End of range to analyze histogram
	FirstBin = histogram -> FindFirstBinAbove( FirstBinMinValue, 1 );	//first bin above FirstBinMinValue on 1 -> (x axis) - Start of range to analyze histogram
	std::cout << "MaxBin: \t" << BinMax << std::endl;
	std::cout << "LastBin: \t" << LastBin  << std::endl;
	std::cout << "FirstBin \t" << FirstBin << std::endl;
//--------------------------------------- Pushing Histogram data into two vectors for analysis of background and finding range of the fit
	for(unsigned i=0; i<LastBin+1; i++)
	{
		Values.push_back( histogram->GetBinContent(i) );		//BinContent - Y
		Arguments.push_back( histogram->GetBinCenter(i) );		//Argument - BinCenter - X
	}
	delete histogram;							//Histogram is not longer needed - having vectors to analyze
	double mean = 0;							//Mean for calculation of background
	double sd = 0;								//Standard deviation of the background
//--------------------------------------- Calculating Background as mean of bins from the boundaries of the histogram (+-shift in loop to avoid bins from the ends of histogram which can not follow the right distribution) 
	if( SideOfBackgroundEstimation == "right" )		//Right side -> lower values of arguments (From Last bin to Last bin - Number of points for estimation)
	{
		for(unsigned j=1; j<=BackgroundEstimationNumbOfPoints; j++)
			mean+=Values[LastBin - ShiftForBackgroundEstimation - j];
	}
	else							//Left side -> lower values of arguments (From First bin to Number of points for estimation)
	{
		for(unsigned j=1; j<=BackgroundEstimationNumbOfPoints; j++)
			mean+=Values[FirstBin + ShiftForBackgroundEstimation + j];
	}
	Background = mean/BackgroundEstimationNumbOfPoints;	//Background - mean of the points at the edges of histogram
	std::cout << "Background: \t" << Background << std::endl;
	if( SideOfBackgroundEstimation == "right" )		//The same as above
	{
		for(unsigned k=1; k<=BackgroundEstimationNumbOfPoints; k++)
			sd+= pow(Values[LastBin - ShiftForBackgroundEstimation - k] - Background,2.0);
	}
	else
	{
		for(unsigned k=1; k<=BackgroundEstimationNumbOfPoints; k++)
			sd+= pow(Values[FirstBin + ShiftForBackgroundEstimation + k] - Background,2.0);
	}
	SDBackground = sqrt(sd/BackgroundEstimationNumbOfPoints/(BackgroundEstimationNumbOfPoints - 1) );	//Estimator for standard deviation 
	std::cout << "Background Uncertainity: \t" << SDBackground << std::endl;
//--------------------------------------- Setting range of the fit, Range_From - start of the fit, Range_To - end of the fit. By default Range_From is the first bin than is more than some threshold, calculated as the fraction of maximum bin. End of the fit value given from FitDetails
	short pointToFitCheck = 1, pointFromFitCheck = 1;		//If both checks are set to 0, loop for finding FitRange is ended. First check is for End of fit, second for beginning
	int iterator = 0;						//Iterator that counts the shift from the bin with maximal content for searching for the last ans first bin for fitting
	while( pointToFitCheck )			//In general first loop should find the bin for the start of the fit range - because it is closer to maximum(starting point)
	{
		iterator++;
		if( pointFromFitCheck && (Values[BinMax - iterator] < StartOfFitValue*(Values[BinMax] - Background)) ) 		//-Background to normalize the height of bin contents
		{
			pointFromFitCheck = 0;			//First bin for fitting is found, first check done
			Range_From = BinMax - iterator;
		}
		if( Arguments[BinMax + iterator] >= EndOfFitValue )
		{
			pointToFitCheck = 0;			//Last bin for fitting is found, second check done
			Range_To = BinMax + iterator;
		}
	}	
        std::cout << "Max bin in [Argument] [Value]\t \t \t \t \t \t" << Arguments[BinMax] << "   " << Values[BinMax] << std::endl;
	std::cout << "[Beginning of the range - Bin] [End of the range - Bin] \t \t" << Range_From << "   " << Range_To << std::endl;
	std::cout << "[Beginning of the range - Argument] [End of the range - Argument] \t" << Arguments[Range_From] << "   " << Arguments[Range_To] << std::endl;
	std::cout << "[Beginning of the range - Value] [End of the range - Value] \t \t" << Values[Range_From] << "   " << Values[Range_To] << std::endl;
//-------------------------------------- Cutting the path such that only the name of the file is in it (cutting up to last slash position), cuuting also the extension of the file
	std::string dot = ".";
	std::string slash = "/";
	std::size_t dotPlace = Path.rfind(dot);
	std::size_t slashPlace = Path.rfind(slash);
	Path.replace( dotPlace, Path.length() - dotPlace, "_RF_" + NumberToChar(StartOfFitValue,2) + "_RT_" + NumberToChar(EndOfFitValue,2) );	//RF - Range_From, RT - Range_To, to distinguish differet ranges of the fit
	if( slashPlace != std::string::npos )
	{
		Path.erase( Path.begin(), Path.begin() + slashPlace + 1);
	}
	PathWithDate = Path + "_Date:" + GetTime();
}

int Fit::Discrete()
{
  	if( Resolution.size()*Lifetimes.size() == 0 )				//Test if the reading FitDetails file and data file past properly
	{
		return 0;
	}
        std::cout<<"-------------Fitting data-------------"<<std::endl;
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

	std::cout << "[Argument for LastBin] [Argument for FirstBin] \t" << Arguments[LastBin] << "   " << Arguments[FirstBin] << std::endl;

//---------------------------------------Range of drawing histogram with fit on it. first two lines is default, next to lines is for setting by user, if strict range is needed

	float DrawFrom = Arguments[FirstBin];
	float DrawTo = Arguments[LastBin];

//---------------------------------------

	TH1F* histogram = FillHistogram( "Time differences", (DrawTo - DrawFrom)*1/BinWidth, DrawFrom, DrawTo, Times );
	if( NmbOfBins )
	{
		for(unsigned i=FirstBin; i<LastBin+1; i++)
		{
			histogram->SetBinContent(i-FirstBin+1,Times[i]);
		}
	}	

	histogram -> GetXaxis() -> SetTitle( "Time difference [ns]" );
	histogram -> GetYaxis() -> SetTitle( "Counts" );
	histogram -> GetYaxis()-> SetTitleOffset(1.3);
	histogram -> GetXaxis() -> SetTitleFont(62);
	histogram -> GetXaxis() -> SetLabelFont(62);
	histogram -> GetYaxis() -> SetTitleFont(62);
	histogram -> GetYaxis() -> SetLabelFont(62);
	histogram -> SetTitle( ("Time differences_ " + Path).c_str() );
	histogram -> Draw( "" );
	histogram -> SetLineWidth( 2 );
	histogram -> SetLineColor( 1 );


//---------------------------------------Setting the fit function
	int lastbin = LastBin; 
	int firstbin = FirstBin;
	while( histogram -> GetBinCenter( lastbin ) > Arguments[ Range_To ] )
	{
		lastbin--;
	}
	while( histogram -> GetBinCenter( firstbin ) < Arguments[ Range_From ] )
	{
		firstbin++;
	}
	std::cout << "[LastBin - new!, of the DrawRange] [End of Range] - test if they are in proper order (LastBin < EndOfRange)\t" << lastbin << "   " << Range_To << std::endl;
	double area = 0.013*(histogram -> Integral( firstbin, lastbin ) - Background*( Arguments[ lastbin ] - Arguments[ firstbin ] ) );
		    //0.013 value... do not know how this value appeared, but after changing it to other value area, area is worse estimated..Probably this is somehow the result of DifferenceOfArgumentsInFittingFunction/BinWidthOfHistogram???!!!
	std::cout << "Area \t \t \t \t \t \t " << area << std::endl;
	area = fabs(area);
	std::cout << "Range - Arguments [Start] [End] \t \t" << Arguments[ Range_From ] << "   " << Arguments[ Range_To ] << std::endl;
	std::cout << "Range - Values [Start] [End] \t \t \t" << Values[ Range_From ] << "   " << Values[ Range_To ] << std::endl;
	
	int pPsIndex = -1;	
	for( unsigned i=0; i<Lifetimes.size(); i++ )
	{
		if( Lifetimes[i].Type == "ps" )
			pPsIndex = 0;
	}
	TF1 *Discrete;
	if( pPsIndex + 1 )
		Discrete = new TF1( "Discrete", DiscreteFitFunctionPs, Arguments[ Range_From ], Arguments[ Range_To ], 4 + Lifetimes.size()*2 + 3*Resolution.size() + 1 );
	else
		Discrete = new TF1( "Discrete", DiscreteFitFunctionNoPs, Arguments[ Range_From ], Arguments[ Range_To ], 4 + Lifetimes.size()*2 + 3*Resolution.size() + 1 );
	Discrete -> SetNpx(Range_To - Range_From);		//For proper drawing of fit
	/*Nmbr of arguments - 
	1x Background
	1x nmbr of gauss
	1x nmbr of lf comp
	1x total area
	Lifetimes.size() - 1x Intensity, 1x Lifetime
	Resolution.size() - 1x Sigma, 1x Offset, 1x Fraction
	1x fixed components
	*/
 	Discrete -> SetParameter( 0, Background );		//Background
	Discrete -> SetParName( 0, "Background" );
	Discrete -> SetParLimits( 0, Background - 3*SDBackground, Background + 3*SDBackground );
	Discrete -> FixParameter( 1, Resolution.size() );			//Nmbr of Gauss comp
	Discrete -> SetParName( 1, "Number of Gauss resolution components" );
	Discrete -> FixParameter( 2, Lifetimes.size() );			//Nmbr of Lifetime comp
	Discrete -> SetParName( 2, "Number of Lifetimes components" );
	Discrete -> SetParameter( 3, area );			//Tot area of components
	Discrete -> SetParName( 3, "Total area of components" );

	for( unsigned i = 0; i < Resolution.size(); i++ )
	{
		if( FixGauss == "yes" )
		{
			Discrete -> FixParameter( 4 + 3*i, Resolution[i].Lifetime );
			Discrete -> SetParName( 4 + 3*i, ("Sigma for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
			
			Discrete -> FixParameter( 5 + 3*i, SetIntensityParameter( Resolution, i+1 ) );
			Discrete -> SetParName( 5 + 3*i, ("Fraction parameter for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
			
			Discrete -> FixParameter( 6 + 3*i, Arguments[ BinMax ] + i*0.33 );
			Discrete -> SetParName( 6 + 3*i, ("Offset for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
		}
		else
		{
			Discrete -> SetParameter( 4 + 3*i, Resolution[i].Lifetime );
			Discrete -> SetParName( 4 + 3*i, ("Sigma for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
			Discrete -> SetParLimits( 4 + 3*i, 0.01, 5 );
			
			Discrete -> SetParameter( 5 + 3*i, SetIntensityParameter( Resolution, i+1 ) );
			Discrete -> SetParName( 5 + 3*i, ("Fraction parameter for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
			if( i == 0 )
				Discrete -> SetParLimits( 5 + 3*i, 0, 1 );
			else
				Discrete -> SetParLimits( 5 + 3*i, 0, 3.14159/2 );
			
			Discrete -> SetParameter( 6 + 3*i, Arguments[ BinMax ] + i*0.33 );
			Discrete -> SetParName( 6 + 3*i, ("Offset for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
			Discrete -> SetParLimits( 6 + 3*i, Arguments[ BinMax ] - 0.5, Arguments[ BinMax ] + 1 );
		}
	}
	unsigned NotFixedIterator = 0, FixedIterator = 0;
	for( unsigned j = 0; j < Lifetimes.size(); j++ )
	{
		Discrete -> SetParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
 		Discrete -> SetParName( 4 + 3*Resolution.size() + 2*j, ("Lifetime for " + NumberToChar( j+1, 0 ) + " Component").c_str() );
	
		Discrete -> SetParName( 5 + 3*Resolution.size() + 2*j, ("Intensity for " + NumberToChar( j+1, 0 ) + " Component").c_str() );
		
		if( Lifetimes[j].Type == "f" || Lifetimes[j].Type == "pf" || Lifetimes[j].Type == "lpf" || Lifetimes[j].Type == "lf" || Lifetimes[j].Type == "lff" )
		{
			Discrete -> FixParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
			Discrete -> FixParameter( 5 + 3*Resolution.size() + 2*j, Lifetimes[j].Intensity );
			FixedIterator++;
		}
		else if( Lifetimes[j].Type == "ps" )
		{
			Discrete -> FixParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
			Discrete -> FixParameter( 5 + 3*Resolution.size() + 2*j, Lifetimes[j].Intensity );
			pPsIndex = j;
		}
		else
		{
			NotFixedIterator++;
			Discrete -> SetParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
			Discrete -> SetParLimits( 4 + 3*Resolution.size() + 2*j, 0.1, 142 );

			Discrete -> SetParameter( 5 + 3*Resolution.size() + 2*j, SetIntensityParameter( LifetimesNotFixed, NotFixedIterator ) );
			if( NotFixedIterator == 1 )
				Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, 0, 1 );
			else
				Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, 0, 3.14159/2 );		
		}
	}
	Discrete -> FixParameter( 6 + 3*Resolution.size() + 2*(Lifetimes.size()-1), FixedIterator );
 	Discrete -> SetParName( 6 + 3*Resolution.size() + 2*(Lifetimes.size()-1), "Number of somehow fixed Components not ps type" );
	
        histogram -> Fit(Discrete,"RM");
	std::cout << std::endl;
	std::cout <<"Do not look on intensities of completely free parameters and resolution intensities"<< std::endl;
	std::cout <<"They are normalized in the body of fitting function and at the end when saving results"<< std::endl;
	std::cout << std::endl;
//---------------------------------------
//---------------------------------------Iterations	
	for( unsigned iter = 0; iter < NmbrOfIterations; iter++ )
	{
		Discrete -> SetParameter( 0, Discrete -> GetParameter(0) );		//Background
		Discrete -> SetParLimits( 0, Background - 3*SDBackground, Background + 3*SDBackground );
		Discrete -> FixParameter( 1, Resolution.size() );			//Nmbr of Gauss comp
		Discrete -> FixParameter( 2, Lifetimes.size() );			//Nmbr of Lifetime comp
		Discrete -> SetParameter( 3, Discrete -> GetParameter(3) );

		for( unsigned i = 0; i < Resolution.size(); i++ )
		{	
			if( FixGauss == "yes" )
			{
				Discrete -> FixParameter( 4 + 3*i, Discrete -> GetParameter(4 + 3*i) );

				Discrete -> FixParameter( 5 + 3*i, Discrete -> GetParameter(5 + 3*i) );

				Discrete -> FixParameter( 6 + 3*i, Discrete -> GetParameter(6 + 3*i) );
			}
			else
			{
				Discrete -> SetParameter( 4 + 3*i, Discrete -> GetParameter(4 + 3*i) );
				Discrete -> SetParLimits( 4 + 3*i, 0.01, 5 );

				Discrete -> SetParameter( 5 + 3*i, Discrete -> GetParameter(5 + 3*i) );

				Discrete -> SetParameter( 6 + 3*i, Discrete -> GetParameter(6 + 3*i) );
			}
		}
		NotFixedIterator = 0;
		for( unsigned j = 0; j < Lifetimes.size(); j++ )
		{
			Discrete -> SetParameter( 4 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) );

			Discrete -> SetParameter( 5 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) );
			if( Lifetimes[j].Type == "f" )
			{
				Discrete -> FixParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
				Discrete -> FixParameter( 5 + 3*Resolution.size() + 2*j, Lifetimes[j].Intensity );
			}
			else if( Lifetimes[j].Type == "pf" )
			{
				Discrete -> SetParLimits( 4 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) - VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j), Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) + VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) );
				Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, 0, 1 );			
			}
			else if( Lifetimes[j].Type == "lpf" )
			{
				Discrete -> SetParLimits( 4 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) - VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j), Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) + VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) );
				Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) - VarLvl*Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j), Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) + VarLvl*Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) );
			}
			else if( Lifetimes[j].Type == "lf" )
			{
				Discrete -> FixParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
				Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, 0, 1 );
			}
			else if( Lifetimes[j].Type == "lff" )
			{
				Discrete -> FixParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
				Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) - VarLvl*Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j), Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) + VarLvl*Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) );
			}
			else if( Lifetimes[j].Type == "ps" )
			{
				Discrete -> SetParLimits( 4 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) - VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j), Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) + VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) );
				Discrete -> FixParameter( 5 + 3*Resolution.size() + 2*j, Lifetimes[j].Intensity );		  
			}
			else
			{
				NotFixedIterator++;
				Discrete -> SetParLimits( 4 + 3*Resolution.size() + 2*j, 0.01, 142 );
				if( NotFixedIterator == 1 )
			  		Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, 0, 1 );
			  	else
			  		Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, 0, 3.14159/2 );
			}
		}
		Discrete -> FixParameter( 6 + 3*Resolution.size() + 2*(Lifetimes.size()-1), FixedIterator );
		histogram -> Fit(Discrete,"RM");
		std::cout << std::endl;
		std::cout <<"Do not look on intensities of completely free parameters and resolution intensities"<< std::endl;
		std::cout <<"They are normalized in the body of fitting function and at the end when saving results"<< std::endl;
		std::cout << std::endl;
	}
//---------------------------------------

//---------------------------------------Drawing and saving histograms and results

        Discrete->Draw("same");
	TFile *testFile = new TFile( "Results_from_fit.root", "update" );
        testFile -> mkdir( Path.c_str() );
        testFile -> cd( Path.c_str() );
	std::string Res_path = "Results";

	if( not PathCheck( Res_path ) )
        {
                std::cout << "Creating directory for results" << std::endl;
                boost::filesystem::path dir( Res_path );
                boost::filesystem::create_directories( dir );
        }
	
        c1->SaveAs( ( Res_path + "/" + PathWithDate + ".png" ).c_str() );
        c1->SetLogy();
        c1->SaveAs( ( Res_path + "/" + PathWithDate + "_logscale"+".png" ).c_str() );
	c1 -> Write( PathWithDate.c_str() );

        AreaFromDiscrete = Discrete -> GetParameter( 3 );
	Background = Discrete -> GetParameter( 0 );
	Double_t ResolutionsFromFit[22];
	ResolutionsFromFit[0] = 1;
	ResolutionsFromFit[1] = Resolution.size();
	Double_t ResolutionsFromFitErrors[22];
	ResolutionsFromFitErrors[0] = Resolution.size();
	double FreeIntensities = 0., FixedIntensities = 0., pPsIntensity = 0., FreeIntensitiesTopPs = 0., FixedFixedIntensity = 0.;
	for( unsigned i = 0; i < Resolution.size(); i++ )
	{
		ResolutionsFromFit[2*(i+1)] = Discrete -> GetParameter( 4 + 3*i );
		ResolutionsFromFit[2*(i+1)+1] = Discrete -> GetParameter( 5 + 3*i );
		ResolutionsFromFitErrors[2*i+1] = Discrete -> GetParameter( 5 + 3*i );
		ResolutionsFromFitErrors[2*i+2] = Discrete -> GetParError( 5 + 3*i );
	}
	Double_t FreeParameters[22];
	FreeParameters[0] = 1;
	FreeParameters[1] = NotFixedIterator;
	unsigned AfterFitIterator = 0;
	Double_t FreeParametersErrors[22];
	FreeParametersErrors[0] = NotFixedIterator;
	for( unsigned i = 0; i < Lifetimes.size(); i++ )
	{
		if( Lifetimes[i].Type == "nf" )
		{
			AfterFitIterator++;
			FreeParameters[2*AfterFitIterator] = Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*i );
			FreeParameters[2*AfterFitIterator+1] = Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
			FreeParametersErrors[2*AfterFitIterator-1] = Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
			FreeParametersErrors[2*AfterFitIterator] = Discrete -> GetParError( 5 + 3*Resolution.size() + 2*i );
		}
	}
	AfterFitIterator = 0;
	for( unsigned i = 0; i < Lifetimes.size(); i++ )
	{
		if( Lifetimes[i].Type != "nf" && Lifetimes[i].Type != "ps" )
		{
			if( Lifetimes[i].Type == "f" )
				FixedFixedIntensity += Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
			FixedIntensities += Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
			if( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*i ) > 0.7 && pPsIndex >=0 )
				pPsIntensity += Discrete -> GetParameter( 5 + 3*Resolution.size() )*Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
		}
		else if( Lifetimes[i].Type != "ps" )
		{
			AfterFitIterator++;
			FreeIntensities += GetIntensityParameterNew( FreeParameters, 3, 3, AfterFitIterator );
			if( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*i ) > 0.7 && pPsIndex >=0 )
				FreeIntensitiesTopPs += GetIntensityParameterNew( FreeParameters, 3, 3, AfterFitIterator );
		}
	}
	FixedIntensities += pPsIntensity;
	NotFixedIterator = 0;
	if( pPsIndex >=0 )
	{
		pPsIntensity += Discrete -> GetParameter( 5 + 3*Resolution.size() )*FreeIntensitiesTopPs*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) );
	}
	
	std::ofstream res;
	res.open( Res_path + "/" + "Discrete_Fit_" + PathWithDate );
	res << "Results from fitting" << std::endl;
	res << "Parameter name\t \t \t \t \t" << "Value\t \t \t" << "Error" << std::endl;
	res << std::endl;
	res << std::endl;
	res << "Background\t \t \t \t \t" << NumberToChar( Background, 5 ) << "\t \t" << SDBackground << std::endl;
	for( unsigned i = 0; i < Resolution.size(); i++ )
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
		res << "Lifetime for p-PS Component\t \t \t " << NumberToChar( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), 5) << "\t \t" << Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) << std::endl;
		LifetimesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) ) );
		res << "Intensity for p-PS Component \t \t \t " <<  NumberToChar( pPsIntensity, 5) << std::endl;
		res << "Intensity for p-PS Component in percent \t " <<  NumberToChar( pPsIntensity * 100, 5) << std::endl;
		res << "Intensity for p-PS Component in percent no fix" << "\t " <<  NumberToChar( pPsIntensity * 100/(1-FixedFixedIntensity), 5 ) << std::endl;	
		res << "Intensity for p-PS Component in percent of o-Ps\t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() ) * 100, 5)  << std::endl;
		IntensitiesFromDiscrete.push_back( DiscreteFitResult( pPsIntensity, 0 ) );
		res << std::endl;
	}
	for( unsigned j = 0; j < Lifetimes.size(); j++ )
	{
		if( j < FixedIterator + pPsIndex + 1 && (int)j != pPsIndex )
		{
			res << ("Lifetime for " + NumberToChar( j+1, 0 ) + " Component [ns]").c_str()
				  << "\t \t \t " << NumberToChar( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), 5) << "\t \t" << Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) << std::endl;
			LifetimesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) ) );
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component").c_str() 
				  << "\t \t \t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ), 5) << "\t \t" << 
					        Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j )  << std::endl;
			IntensitiesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ), Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) ) );
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component in percent").c_str()
		                  << "\t \t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ) * 100, 5 ) << "\t \t" <<
		                               Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) * 100 << std::endl;
			if( Lifetimes[j].Type != "f" )
			{
				  res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component in percent no fix").c_str()
		                  << "\t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ) * 100/(1-FixedFixedIntensity), 5 ) << "\t \t" <<
		                               Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) * 100/(1-FixedFixedIntensity) << std::endl;
			}
			else
			{
				LifetimesFromDiscrete[ LifetimesFromDiscrete.size() - 1 ].Type = "f";
			}
		        std::cout << "Intensity in percent: \t" << Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ) * 100 << std::endl;
			res << std::endl;
		}
		else if( (int)j != pPsIndex )
		{
			NotFixedIterator++;
			res << ("Lifetime for " + NumberToChar( j+1, 0 ) + " Component [ns]").c_str()
				  << "\t \t \t " << NumberToChar( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), 5) << "\t \t" << Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) << std::endl;
			LifetimesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) ) );
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component").c_str() 
				  << "\t \t \t " <<  NumberToChar( GetIntensityParameterNew( FreeParameters, 3, 3, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) , 5) << "\t \t" << GetIntensityParameterErrorNew( FreeParametersErrors, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) << std::endl;	  
			IntensitiesFromDiscrete.push_back( DiscreteFitResult( GetIntensityParameterNew( FreeParameters, 3, 3, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ), GetIntensityParameterErrorNew( FreeParametersErrors, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) ) );
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component in percent").c_str()
		                  << "\t \t " <<  NumberToChar( GetIntensityParameterNew( FreeParameters, 3, 3, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100, 5 ) << "\t \t" << GetIntensityParameterErrorNew( FreeParametersErrors, NotFixedIterator ) *(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100 << std::endl;
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component in percent no fix").c_str()
		                  << "\t " <<  NumberToChar( GetIntensityParameterNew( FreeParameters, 3, 3, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/(1-FixedFixedIntensity), 5 ) << "\t \t" << GetIntensityParameterErrorNew( FreeParametersErrors, NotFixedIterator ) *(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/(1-FixedFixedIntensity) << std::endl;	  
		        std::cout << "Intensity in percent: \t" << GetIntensityParameterNew( FreeParameters, 3, 3, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100 << std::endl;
			res << std::endl;
		}
	}

//---------------------------------------Sorting the lifetime components, so the first component has the lowest lifetime, second, the second lowest, ...
	
	std::vector<unsigned> Order;
	double minTemp = 200, minPrevious = 0;
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
		minTemp = 200;
		if( Order.size() == LifetimesFromDiscrete.size() )
			minSearch = false;
	}
//---------------------------------------
	
	double Frac2GParameter = 0, Frac2GParameterNew = 0, Frac3GParameter = 0, Frac3GParameterNew = 0, Frac3GParameterNew2 = 0, DeltaParameter = 0;
	if( LifetimesFromDiscrete.size() == 2 )
	{
		Frac2GParameter = 371/372;
		Frac2GParameterNew = 371/372;
		Frac3GParameter = 1/372;
		Frac3GParameterNew = 1/372;
		Frac3GParameterNew2 = 1/372;
	}
	else
	{
		for( unsigned i=1; i<LifetimesFromDiscrete.size(); i++ )
		{
			if( LifetimesFromDiscrete[ Order[i] ].Parameter < 0.8 )
			{
				Frac2GParameter += 371*IntensitiesFromDiscrete[ Order[i] ].Parameter / 372;
				Frac2GParameterNew += 371*IntensitiesFromDiscrete[ Order[i] ].Parameter / 372;
				Frac3GParameter += IntensitiesFromDiscrete[ Order[i] ].Parameter / 372;
				Frac3GParameterNew += IntensitiesFromDiscrete[ Order[i] ].Parameter / 372;
			}
			else
			{
				Frac2GParameter += IntensitiesFromDiscrete[ Order[i] ].Parameter  / 3;
				Frac2GParameter += IntensitiesFromDiscrete[ Order[i] ].Parameter * 142 / ( 142 - LifetimesFromDiscrete[ Order[i] ].Parameter );
				Frac2GParameterNew += IntensitiesFromDiscrete[ Order[i] ].Parameter  / 3;
				Frac2GParameterNew += IntensitiesFromDiscrete[ Order[i] ].Parameter * ( 142 - LifetimesFromDiscrete[ Order[i] ].Parameter ) / 142;
				Frac3GParameter += IntensitiesFromDiscrete[ Order[i] ].Parameter * LifetimesFromDiscrete[ Order[i] ].Parameter / ( 142 - LifetimesFromDiscrete[ Order[i] ].Parameter );
				Frac3GParameterNew2 += IntensitiesFromDiscrete[ Order[i] ].Parameter * LifetimesFromDiscrete[ Order[i] ].Parameter / ( 142 - LifetimesFromDiscrete[ Order[i] ].Parameter );
				Frac3GParameterNew += IntensitiesFromDiscrete[ Order[i] ].Parameter * LifetimesFromDiscrete[ Order[i] ].Parameter / 142;
				DeltaParameter += (1-4*IntensitiesFromDiscrete[ Order[i] ].Parameter/3)/372 + IntensitiesFromDiscrete[ Order[i] ].Parameter * LifetimesFromDiscrete[ Order[i] ].Parameter / 142;
			}
		}
	}	
	res << "Frac2GParameter \t" << Frac2GParameter << std::endl;
	res << "Frac2GParameterNew \t" << Frac2GParameterNew << std::endl;	
	res << "Frac3GParameter \t" << Frac3GParameter << std::endl;
	res << "Frac3GParameterCorr \t" << Frac3GParameterNew2 << std::endl;
	res << "Frac3GParameternew \t" << Frac3GParameterNew << std::endl;
	res << "DeltaParameter \t" << DeltaParameter << std::endl;			
	
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

	res << std::endl;
	res << "ChiSquared" << "\t \t \t \t \t" << NumberToChar( Discrete -> GetChisquare(), 4 ) << std::endl;
	res << "Probability of the fit " << "\t \t \t \t" << NumberToChar(1 - Discrete -> GetProb(), 4) << std::endl;
	res << "Rsquared" << "\t \t \t \t \t" << NumberToChar(1 - (Discrete -> GetChisquare())/Diff, 4 ) << std::endl;
	res << "Adjusted Rsquared" << "\t \t \t \t" << NumberToChar(1 - ((sumIT-1)/Discrete -> GetNDF())*(Discrete -> GetChisquare())/Diff, 4 ) << std::endl;
	res << "ChiSquared/DegreesOfFreedom" << "\t \t \t" << NumberToChar(  Discrete -> GetChisquare()/Discrete -> GetNDF(), 4 ) << std::endl;
	res.close();
	
//---------------------------------------Writing to csv	
	bool IsFileExists =  FileCheck( Res_path + "/" + "Discrete_Fit_" + Path + ".csv" );
	std::ofstream res_csv;
	if( IsFileExists )
		res_csv.open( Res_path + "/" + "Discrete_Fit_" + Path + ".csv", std::ofstream::app );
	else
		res_csv.open( Res_path + "/" + "Discrete_Fit_" + Path + ".csv" );
	res_csv << PathWithDate;
	res_csv << ",Frac2GParameter,Frac3GParameter,Frac3GParameterCorr,Frac3GParameternew,DeltaParameter,ChiSquared,Probability_of_the_fit,Rsquared,Adjusted_Rsquared,ChiSquared/DegreesOfFreedom,";
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
		/*if( j < FixedIterator + pPsIndex + 1 && (int)j != pPsIndex )
		{
			res_csv << "Lifetime_for_" << NumberToChar( j+1, 0 )<< "_Component_[ns],";
			res_csv << "Error_of_Lifetime_for_" << NumberToChar( j+1, 0 )<< "_Component_[ns],";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component,";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent,";
			if( Lifetimes[j].Type != "f" )
			{
				res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent_no_fix,";
				res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent_no_fix,";			  
			}
		}
		else if( (int)j != pPsIndex )
		{
			res_csv << "Lifetime_for_" << NumberToChar( j+1, 0 )<< "_Component_[ns],";
			res_csv << "Error_of_Lifetime_for_" << NumberToChar( j+1, 0 )<< "_Component_[ns],";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component,";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent,";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent_no_fix,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent_no_fix,";
		}*/
	}
	res_csv << "\n";
	res_csv << ",";
	res_csv << NumberToChar( Frac2GParameter, 5 ) << "," << NumberToChar( Frac3GParameter, 5 ) << NumberToChar( Frac3GParameterNew2, 5 ) << "," << "," << NumberToChar( Frac3GParameterNew, 5 ) << "," << NumberToChar( DeltaParameter, 5 ) << ",";
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
		/*if( j < FixedIterator + pPsIndex + 1 && (int)j != pPsIndex )
		{
			res_csv << NumberToChar( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), 5) << ",";
			res_csv << Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) << ",";
			res_csv << NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ), 5) << ",";
			res_csv << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) << ",";
			res_csv << NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ) * 100, 5 ) << ",";
			res_csv << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) * 100 << ",";
			if( Lifetimes[j].Type != "f" )
			{
				res_csv << NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ) * 100/(1-FixedFixedIntensity), 5 ) << ",";
				res_csv << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) * 100/(1-FixedFixedIntensity) << ",";			  
			}
		}
		else if( (int)j != pPsIndex )
		{
			NotFixedIterator++;
			res_csv << NumberToChar( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), 5) << ",";
			res_csv << Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) << ",";
			res_csv << NumberToChar( GetIntensityParameterNew( FreeParameters, 3, 3, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) , 5) << ",";
			res_csv << GetIntensityParameterErrorNew( FreeParametersErrors, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) << ",";
			res_csv << NumberToChar( GetIntensityParameterNew( FreeParameters, 3, 3, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100, 5 ) << ",";
			res_csv << GetIntensityParameterErrorNew( FreeParametersErrors, NotFixedIterator ) *(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100 << ",";
			res_csv << NumberToChar( GetIntensityParameterNew( FreeParameters, 3, 3, NotFixedIterator )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/(1-FixedFixedIntensity), 5 ) << ",";
			res_csv << GetIntensityParameterErrorNew( FreeParametersErrors, NotFixedIterator ) *(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/(1-FixedFixedIntensity) << ",";
		}*/
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
	testFile->Close();
//---------------------------------------

//---------------------------------------Residuals	
	delete c1;
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

	Double_t Residuals_root[ Range_To - Range_From ], ResArg_root[ Range_To - Range_From ];
	for( int i = 0; i < histogram->GetNbinsX(); i++ )
	{
		if( histogram->GetBinCenter(i) > Arguments[ Range_From ] && histogram->GetBinCenter(i) < Arguments[ Range_To ] && i >= (int)Range_From )
		{
			if( isnan(histogram->GetBinContent(i) - Discrete->Eval( histogram->GetBinCenter(i) ) )/sqrt( histogram->GetBinContent(i) ) || isinf(histogram->GetBinContent(i) - Discrete->Eval( histogram->GetBinCenter(i) ) )/sqrt( histogram->GetBinContent(i) ) )
				Residuals_root[i-Range_From] = 0;
			else
				Residuals_root[i-Range_From] = (histogram->GetBinContent(i) - Discrete->Eval( histogram->GetBinCenter(i) ) )/sqrt( histogram->GetBinContent(i) );
			ResArg_root[i-Range_From] = histogram->GetBinCenter(i);
		}
	}

	TFile *residuals = new TFile( "Residuals.root", "update" );
	residuals -> mkdir( Path.c_str() ); 
	residuals -> cd( Path.c_str() );  
        TGraph *residualsHisto = new TGraph( Range_To - Range_From, ResArg_root, Residuals_root );
        residualsHisto -> SetMarkerStyle( 21 );
        residualsHisto -> SetMarkerSize( 0.5 );
        residualsHisto -> SetTitle( "Residuals of the fit" );
        residualsHisto -> Draw( "APL" );
        residualsHisto -> GetXaxis() -> SetTitle( "Time difference [ns]" );
        residualsHisto -> GetYaxis() -> SetTitle( "Residuals" );
        residualsHisto -> Write( PathWithDate.c_str() );
        delete residualsHisto;

	residuals->Close();
//---------------------------------------
	delete histogram;
	delete Discrete;
	delete c2;
//---------------------------------------Saving graphical view of discret results
	
	TFile *comparison = new TFile( "LF_vs_In.root", "update" );
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
	//LF_vs_In -> SetLineColor(0);
	LF_vs_In -> Draw("AP");	
	LF_vs_In -> Write( PathWithDate.c_str() );
	
	delete LF_vs_In;
	comparison->Close();
//---------------------------------------	
	return DecoOption;
}

int Fit::Discrete_old()
{
  	if( Resolution.size()*Lifetimes.size() == 0 )			//Test if the reading FitDetails file and data file past properly
	{
		return 0;
	}
        std::cout<<"-------------Fitting data-------------"<<std::endl;
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

	std::cout << "[Argument for LastBin] [Argument for FirstBin] \t" << Arguments[LastBin] << "   " << Arguments[FirstBin] << std::endl;

//---------------------------------------Range of drawing histogram with fit on it. first two lines is default, next to lines is for setting by user, if strict range is needed

	float DrawFrom = Arguments[FirstBin];
	float DrawTo = Arguments[LastBin];

//---------------------------------------

	TH1F* histogram = FillHistogram( "Time differences", (DrawTo - DrawFrom)*1/BinWidth, DrawFrom, DrawTo, Times );
	if( NmbOfBins )
	{
		for(unsigned i=FirstBin; i<LastBin+1; i++)
		{
			histogram->SetBinContent(i-FirstBin+1,Times[i]);
		}
	}	

	histogram -> GetXaxis() -> SetTitle( "Time difference [ns]" );
	histogram -> GetYaxis() -> SetTitle( "Counts" );
	histogram -> GetYaxis()-> SetTitleOffset(1.3);
	histogram -> GetXaxis() -> SetTitleFont(62);
	histogram -> GetXaxis() -> SetLabelFont(62);
	histogram -> GetYaxis() -> SetTitleFont(62);
	histogram -> GetYaxis() -> SetLabelFont(62);
	histogram -> SetTitle( ("Time differences_ " + Path).c_str() );
	histogram -> Draw( "" );
	histogram -> SetLineWidth( 2 );
	histogram -> SetLineColor( 1 );

//---------------------------------------Setting the fit function
	int lastbin = LastBin; 
	int firstbin = FirstBin;
	while( histogram -> GetBinCenter( lastbin ) > Arguments[ Range_To ] )
	{
		lastbin--;
	}
	while( histogram -> GetBinCenter( firstbin ) < Arguments[ Range_From ] )
	{
		firstbin++;
	}

	std::cout << "[LastBin - new!, of the DrawRange] [End of Range] - test if they are in proper order (LastBin < EndOfRange)\t" << lastbin << "   " << Range_To << std::endl;
	double area = 0.013*(histogram -> Integral( firstbin, lastbin ) - Background*( Arguments[ lastbin ] - Arguments[ firstbin ] ) );
	
	std::cout << "Area \t \t \t \t \t \t " << area << std::endl;
	area = fabs(area);
	
	std::cout << "Range - Arguments [Start] [End] \t \t" << Arguments[ Range_From ] << "   " << Arguments[ Range_To ] << std::endl;
	std::cout << "Range - Values [Start] [End] \t \t \t" << Values[ Range_From ] << "   " << Values[ Range_To ] << std::endl;
	
	std::cout <<"-------------Old fit way-------------"<< std::endl;
	
	int pPsIndex = -1;	
	for( unsigned i=0; i<Lifetimes.size(); i++ )
	{
		if( Lifetimes[i].Type == "ps" )
			pPsIndex = 0;
	}
	TF1 *Discrete;
	if( pPsIndex + 1 )
		Discrete = new TF1( "Discrete", DiscreteFitFunctionPs_old, Arguments[ Range_From ], Arguments[ Range_To ], 4 + Lifetimes.size()*2 + 3*Resolution.size() + 1 );
	else
		Discrete = new TF1( "Discrete", DiscreteFitFunctionNoPs_old, Arguments[ Range_From ], Arguments[ Range_To ], 4 + Lifetimes.size()*2 + 3*Resolution.size() + 1 );
	
	
	Discrete -> SetNpx(Range_To - Range_From);
	/*Nmbr of arguments - 
	1x Background
	1x nmbr of gauss
	1x nmbr of lf comp
	1x total area
	Lifetimes.size() - 1x Intensity, 1x Lifetime
	Resolution.size() - 1x Sigma, 1x Offset, 1x Fraction
	1x fixed components
	*/
 	Discrete -> SetParameter( 0, Background );		//Background
	Discrete -> SetParName( 0, "Background" );
	Discrete -> SetParLimits( 0, Background - 3*SDBackground, Background + 3*SDBackground );
	Discrete -> FixParameter( 1, Resolution.size() );			//Nmbr of Gauss comp
	Discrete -> SetParName( 1, "Number of Gauss resolution components" );
	Discrete -> FixParameter( 2, Lifetimes.size() );			//Nmbr of Lifetime comp
	Discrete -> SetParName( 2, "Number of Lifetimes components" );
	Discrete -> SetParameter( 3, area );			//Tot area of components
	Discrete -> SetParName( 3, "Total area of components" );

	for( unsigned i = 0; i < Resolution.size(); i++ )
	{
		if( FixGauss == "yes" )
		{
			Discrete -> FixParameter( 4 + 3*i, Resolution[i].Lifetime );
			Discrete -> SetParName( 4 + 3*i, ("Sigma for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
			
			Discrete -> FixParameter( 5 + 3*i, SetIntensityParameter( Resolution, i+1 ) );
			Discrete -> SetParName( 5 + 3*i, ("Fraction parameter for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
			
			Discrete -> FixParameter( 6 + 3*i, Arguments[ BinMax ] + i*0.33 );
			Discrete -> SetParName( 6 + 3*i, ("Offset for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
		}
		else
		{
			Discrete -> SetParameter( 4 + 3*i, Resolution[i].Lifetime );
			Discrete -> SetParName( 4 + 3*i, ("Sigma for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
			Discrete -> SetParLimits( 4 + 3*i, 0.01, 5 );
			
			Discrete -> SetParameter( 5 + 3*i, SetIntensityParameter( Resolution, i+1 ) );
			Discrete -> SetParName( 5 + 3*i, ("Fraction parameter for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
			if( i == 0 )
				Discrete -> SetParLimits( 5 + 3*i, 0, 1 );
			else
				Discrete -> SetParLimits( 5 + 3*i, 0, 3.14159/2 );
			
			Discrete -> SetParameter( 6 + 3*i, Arguments[ BinMax ] + i*0.33 );
			Discrete -> SetParName( 6 + 3*i, ("Offset for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
			Discrete -> SetParLimits( 6 + 3*i, Arguments[ BinMax ] - 0.5, Arguments[ BinMax ] + 1 );
		}
	}
	unsigned NotFixedIterator = 0, FixedIterator = 0;
	for( unsigned j = 0; j < Lifetimes.size(); j++ )
	{
		Discrete -> SetParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
 		Discrete -> SetParName( 4 + 3*Resolution.size() + 2*j, ("Lifetime for " + NumberToChar( j+1, 0 ) + " Component").c_str() );

/*		
ExMoGa( Double_t A, Double_t P1, Double_t P2, Double_t P3, Double_t P4 )			//A - arguments, p - parameters
P4 / (2*P3) * exp( P1*P1 / (2*P3*P3) - ( A - P2 ) / P3 ) * ( erf( 1 / (sqrt(2)*P1) * ( A - P1*P1 / P3 - P2 ) ) - erf( 1 / (sqrt(2)*P1) * ( - P1*P1 / P3 - P2 ) ) );*/
		
		Discrete -> SetParName( 5 + 3*Resolution.size() + 2*j, ("Intensity for " + NumberToChar( j+1, 0 ) + " Component").c_str() );
		
		if( Lifetimes[j].Type == "f" || Lifetimes[j].Type == "pf" || Lifetimes[j].Type == "lpf" || Lifetimes[j].Type == "lf" || Lifetimes[j].Type == "lff" )
		{
			Discrete -> FixParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
			Discrete -> FixParameter( 5 + 3*Resolution.size() + 2*j, Lifetimes[j].Intensity );
			FixedIterator++;
		}
		else if( Lifetimes[j].Type == "ps" )
		{
			Discrete -> FixParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
			Discrete -> FixParameter( 5 + 3*Resolution.size() + 2*j, Lifetimes[j].Intensity );
			pPsIndex = j;
		}
		else
		{
			NotFixedIterator++;
			Discrete -> SetParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
			Discrete -> SetParLimits( 4 + 3*Resolution.size() + 2*j, 0.1, 142 );
			Discrete -> SetParameter( 5 + 3*Resolution.size() + 2*j, Lifetimes[j].Intensity );
			Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, 0, 1 );
		}
	}
	Discrete -> FixParameter( 6 + 3*Resolution.size() + 2*(Lifetimes.size()-1), FixedIterator );
 	Discrete -> SetParName( 6 + 3*Resolution.size() + 2*(Lifetimes.size()-1), "Number of somehow fixed Components not ps type" );
        histogram -> Fit(Discrete,"RM");
	std::cout << std::endl;
	std::cout <<"Do not look on intensities of completely free parameters and resolution intensities"<< std::endl;
	std::cout <<"They are normalized in the body of fitting function and at the end when saving results"<< std::endl;
	std::cout << std::endl;
//---------------------------------------
//---------------------------------------Iterations	
	for( unsigned iter = 0; iter < NmbrOfIterations; iter++ )
	{
		  
		Discrete -> SetParameter( 0, Discrete -> GetParameter(0) );		//Background
		Discrete -> SetParLimits( 0, Background - 3*SDBackground, Background + 3*SDBackground );
		Discrete -> FixParameter( 1, Resolution.size() );			//Nmbr of Gauss comp
		Discrete -> FixParameter( 2, Lifetimes.size() );			//Nmbr of Lifetime comp
		Discrete -> SetParameter( 3, Discrete -> GetParameter(3) );

		for( unsigned i = 0; i < Resolution.size(); i++ )
		{	
			if( FixGauss == "yes" )
			{
				Discrete -> FixParameter( 4 + 3*i, Discrete -> GetParameter(4 + 3*i) );

				Discrete -> FixParameter( 5 + 3*i, Discrete -> GetParameter(5 + 3*i) );

				Discrete -> FixParameter( 6 + 3*i, Discrete -> GetParameter(6 + 3*i) );
			}
			else
			{
				Discrete -> SetParameter( 4 + 3*i, Discrete -> GetParameter(4 + 3*i) );
				Discrete -> SetParLimits( 4 + 3*i, 0.01, 5 );

				Discrete -> SetParameter( 5 + 3*i, Discrete -> GetParameter(5 + 3*i) );

				Discrete -> SetParameter( 6 + 3*i, Discrete -> GetParameter(6 + 3*i) );
			}
		}
		NotFixedIterator = 0;
		for( unsigned j = 0; j < Lifetimes.size(); j++ )
		{
			Discrete -> SetParameter( 4 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) );

			Discrete -> SetParameter( 5 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) );
			if( Lifetimes[j].Type == "f" )
			{
				Discrete -> FixParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
				Discrete -> FixParameter( 5 + 3*Resolution.size() + 2*j, Lifetimes[j].Intensity );
			}
			else if( Lifetimes[j].Type == "pf" )
			{
				Discrete -> SetParLimits( 4 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) - VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j), Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) + VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) );
				Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, 0, 1 );			
			}
			else if( Lifetimes[j].Type == "lpf" )
			{
				Discrete -> SetParLimits( 4 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) - VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j), Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) + VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) );
				Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) - VarLvl*Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j), Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) + VarLvl*Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) );
			}
			else if( Lifetimes[j].Type == "lf" )
			{
				Discrete -> FixParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
				Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, 0, 1 );
			}
			else if( Lifetimes[j].Type == "lff" )
			{
				Discrete -> FixParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
				Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) - VarLvl*Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j), Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) + VarLvl*Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) );
			}
			else if( Lifetimes[j].Type == "ps" )
			{
				Discrete -> SetParLimits( 4 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) - VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j), Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) + VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) );
				Discrete -> FixParameter( 5 + 3*Resolution.size() + 2*j, Lifetimes[j].Intensity );		  
			}
			else
			{
				NotFixedIterator++;
				Discrete -> SetParLimits( 4 + 3*Resolution.size() + 2*j, 0.01, 142 );
			  	Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, 0, 1 );
			}
		}
		Discrete -> FixParameter( 6 + 3*Resolution.size() + 2*(Lifetimes.size()-1), FixedIterator );
		
		histogram -> Fit(Discrete,"RM");
		std::cout << std::endl;
		std::cout <<"Do not look on intensities of completely free parameters and resolution intensities"<< std::endl;
		std::cout <<"They are normalized in the body of fitting function and at the end when saving results"<< std::endl;
		std::cout << std::endl;
	}
//---------------------------------------

//---------------------------------------Drawing and saving histograms and results
        Discrete->Draw("C same");
	TFile *testFile = new TFile( "Results_from_fit.root", "update" );
        testFile -> mkdir( Path.c_str() );
        testFile -> cd( Path.c_str() );
	std::string Res_path = "Results";
		
	if( not PathCheck( Res_path ) )
        {
                std::cout << "Creating directory for results" << std::endl;
                boost::filesystem::path dir( Res_path );
                boost::filesystem::create_directories( dir );
        }
        
        c1->SaveAs( ( Res_path + "/" + PathWithDate + ".png" ).c_str() );
        c1->SetLogy();
        c1->SaveAs( ( Res_path + "/" + PathWithDate + "_logscale"+".png" ).c_str() );
	c1 -> Write( PathWithDate.c_str() );

	AreaFromDiscrete = Discrete -> GetParameter( 3 );
	Background = Discrete -> GetParameter( 0 );
	Double_t ResolutionsFromFit[22];
	ResolutionsFromFit[0] = 1;
	ResolutionsFromFit[1] = Resolution.size();
	Double_t ResolutionsFromFitErrors[22];
	ResolutionsFromFitErrors[0] = Resolution.size();
	double FreeIntensities = 0., FixedIntensities = 0., pPsIntensity = 0., FreeIntensitiesTopPs = 0., FixedFixedIntensity = 0.;
	for( unsigned i = 0; i < Resolution.size(); i++ )
	{
		ResolutionsFromFit[2*(i+1)] = Discrete -> GetParameter( 4 + 3*i );
		ResolutionsFromFit[2*(i+1)+1] = Discrete -> GetParameter( 5 + 3*i );
		ResolutionsFromFitErrors[2*i+1] = Discrete -> GetParameter( 5 + 3*i );
		ResolutionsFromFitErrors[2*i+2] = Discrete -> GetParError( 5 + 3*i );
	}
	for( unsigned i = 0; i < Lifetimes.size(); i++ )
	{
		if( Lifetimes[i].Type != "nf" && Lifetimes[i].Type != "ps" )
		{
			if( Lifetimes[i].Type == "f" )
				FixedFixedIntensity += Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
			FixedIntensities += Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
			if( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*i ) > 0.7 && pPsIndex >=0 )
				pPsIntensity += Discrete -> GetParameter( 5 + 3*Resolution.size() )*Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
		}
		else if( Lifetimes[i].Type != "ps" )
		{
			FreeIntensities += Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
			if( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*i ) > 0.7 && pPsIndex >=0 )
				FreeIntensitiesTopPs += Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
		}
	}
	FixedIntensities += pPsIntensity;
	if( pPsIndex >=0 )
	{
		pPsIntensity += Discrete -> GetParameter( 5 + 3*Resolution.size() )*FreeIntensitiesTopPs*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) );
	}
	
	std::ofstream res;
	res.open( Res_path + "/" + "Discrete_Fit_" + PathWithDate );
	res << "Results from fitting" << std::endl;
	res << "Parameter name\t \t \t \t \t" << "Value\t \t \t" << "Error" << std::endl;
	res << std::endl;
	res << std::endl;
	res << "Background\t \t \t \t \t" << NumberToChar( Background, 5 ) << "\t \t" << SDBackground << std::endl;
	for( unsigned i = 0; i < Resolution.size(); i++ )
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
		res << "Lifetime for p-PS Component\t \t \t " << NumberToChar( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), 5) << "\t \t" << Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) << std::endl;
		LifetimesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) ) );
		res << "Intensity for p-PS Component \t \t \t " <<  NumberToChar( pPsIntensity, 5) << std::endl;
		res << "Intensity for p-PS Component in percent \t " <<  NumberToChar( pPsIntensity * 100, 5) << std::endl;
		res << "Intensity for p-PS Component in percent no fix" << "\t " <<  NumberToChar( pPsIntensity * 100/(1-FixedFixedIntensity), 5 ) << std::endl;	
		res << "Intensity for p-PS Component in percent of o-Ps\t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() ) * 100, 5)  << std::endl;
		IntensitiesFromDiscrete.push_back( DiscreteFitResult( pPsIntensity, 0 ) );
		res << std::endl;
	}
	for( unsigned j = 0; j < Lifetimes.size(); j++ )
	{
		if( j < FixedIterator + pPsIndex + 1 && (int)j != pPsIndex )
		{
			res << ("Lifetime for " + NumberToChar( j+1, 0 ) + " Component [ns]").c_str()
				  << "\t \t \t " << NumberToChar( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), 5) << "\t \t" << Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) << std::endl;
			LifetimesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) ) );
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component").c_str() 
				  << "\t \t \t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ), 5) << "\t \t" << 
					        Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j )  << std::endl;
			IntensitiesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ), Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) ) );
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component in percent").c_str()
		                  << "\t \t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ) * 100, 5 ) << "\t \t" <<
		                               Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) * 100 << std::endl;
			if( Lifetimes[j].Type != "f" )
			{
				  res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component in percent no fix").c_str()
		                  << "\t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ) * 100/(1-FixedFixedIntensity), 5 ) << "\t \t" <<
		                               Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) * 100/(1-FixedFixedIntensity) << std::endl;
			}
			else
			{
				LifetimesFromDiscrete[ LifetimesFromDiscrete.size() - 1 ].Type = "f";
			}
		        std::cout << "Intensity in percent: \t" << Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ) * 100 << std::endl;
			res << std::endl;
		}
		else if( (int)j != pPsIndex )
		{
			res << ("Lifetime for " + NumberToChar( j+1, 0 ) + " Component [ns]").c_str()
				  << "\t \t \t " << NumberToChar( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), 5) << "\t \t" << Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) << std::endl;
			LifetimesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) ) );
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component").c_str() 
				  << "\t \t \t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) )/FreeIntensities , 5) << "\t \t" << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) )/FreeIntensities << std::endl;
			IntensitiesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) )/FreeIntensities, Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) )/FreeIntensities ) );
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component in percent").c_str()
		                  << "\t \t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/FreeIntensities, 5 ) << "\t \t" << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) *(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/FreeIntensities << std::endl;
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component in percent no fix").c_str()
		                  << "\t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/(1-FixedFixedIntensity)/FreeIntensities, 5 ) << "\t \t" << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/(1-FixedFixedIntensity)/FreeIntensities << std::endl;	  
		        std::cout << "Intensity in percent: \t" << Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/FreeIntensities << std::endl;
			res << std::endl;
		}
	}

//---------------------------------------Sorting the lifetime components, so the first component has the lowest lifetime, second, the second lowest, ...
	
	std::vector<unsigned> Order;
	double minTemp = 200, minPrevious = 0;
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
		minTemp = 200;
		if( Order.size() == LifetimesFromDiscrete.size() )
			minSearch = false;
	}
//---------------------------------------	
	
	double Frac2GParameter = 0, Frac2GParameterNew = 0, Frac3GParameter = 0, Frac3GParameterNew = 0, Frac3GParameterNew2 = 0, DeltaParameter = 0;
	if( LifetimesFromDiscrete.size() == 2 )
	{
		Frac2GParameter = 371/372;
		Frac2GParameterNew = 371/372;
		Frac3GParameter = 1/372;
		Frac3GParameterNew = 1/372;
		Frac3GParameterNew2 = 1/372;
	}
	else
	{
		for( unsigned i=1; i<LifetimesFromDiscrete.size(); i++ )
		{
			if( LifetimesFromDiscrete[ Order[i] ].Parameter < 0.8 )
			{
				Frac2GParameter += 371*IntensitiesFromDiscrete[ Order[i] ].Parameter / 372;
				Frac2GParameterNew += 371*IntensitiesFromDiscrete[ Order[i] ].Parameter / 372;
				Frac3GParameter += IntensitiesFromDiscrete[ Order[i] ].Parameter / 372;
				Frac3GParameterNew += IntensitiesFromDiscrete[ Order[i] ].Parameter / 372;
			}
			else
			{
				Frac2GParameter += IntensitiesFromDiscrete[ Order[i] ].Parameter  / 3;
				Frac2GParameter += IntensitiesFromDiscrete[ Order[i] ].Parameter * 142 / ( 142 - LifetimesFromDiscrete[ Order[i] ].Parameter );
				Frac2GParameterNew += IntensitiesFromDiscrete[ Order[i] ].Parameter  / 3;
				Frac2GParameterNew += IntensitiesFromDiscrete[ Order[i] ].Parameter * ( 142 - LifetimesFromDiscrete[ Order[i] ].Parameter ) / 142;
				Frac3GParameter += IntensitiesFromDiscrete[ Order[i] ].Parameter * LifetimesFromDiscrete[ Order[i] ].Parameter / ( 142 - LifetimesFromDiscrete[ Order[i] ].Parameter );
				Frac3GParameterNew2 += IntensitiesFromDiscrete[ Order[i] ].Parameter * LifetimesFromDiscrete[ Order[i] ].Parameter / ( 142 - LifetimesFromDiscrete[ Order[i] ].Parameter );
				Frac3GParameterNew += IntensitiesFromDiscrete[ Order[i] ].Parameter * LifetimesFromDiscrete[ Order[i] ].Parameter / 142;
				DeltaParameter += (1-4*IntensitiesFromDiscrete[ Order[i] ].Parameter/3)/372 + IntensitiesFromDiscrete[ Order[i] ].Parameter * LifetimesFromDiscrete[ Order[i] ].Parameter / 142;
			}
		}
	}	
	res << "Frac2GParameter \t" << Frac2GParameter << std::endl;
	res << "Frac2GParameterNew \t" << Frac2GParameterNew << std::endl;	
	res << "Frac3GParameter \t" << Frac3GParameter << std::endl;
	res << "Frac3GParameterCorr \t" << Frac3GParameterNew2 << std::endl;
	res << "Frac3GParameternew \t" << Frac3GParameterNew << std::endl;
	res << "DeltaParameter \t" << DeltaParameter << std::endl;	
	
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

	res << std::endl;
	res << "ChiSquared" << "\t \t \t \t \t" << NumberToChar( Discrete -> GetChisquare(), 4 ) << std::endl;
	res << "Probability of the fit " << "\t \t \t \t" << NumberToChar(1 - Discrete -> GetProb(), 4) << std::endl;
	res << "Rsquared" << "\t \t \t \t \t" << NumberToChar(1 - (Discrete -> GetChisquare())/Diff, 4 ) << std::endl;
	res << "Adjusted Rsquared" << "\t \t \t \t" << NumberToChar(1 - ((sumIT-1)/Discrete -> GetNDF())*(Discrete -> GetChisquare())/Diff, 4 ) << std::endl;
	res << "ChiSquared/DegreesOfFreedom" << "\t \t \t" << NumberToChar(  Discrete -> GetChisquare()/Discrete -> GetNDF(), 4 ) << std::endl;
	res.close();
	
//---------------------------------------Writing to csv	
	bool IsFileExists =  FileCheck( Res_path + "/" + "Discrete_Fit_" + Path + ".csv" );
	std::ofstream res_csv;
	if( IsFileExists )
		res_csv.open( Res_path + "/" + "Discrete_Fit_" + Path + ".csv", std::ofstream::app );
	else
		res_csv.open( Res_path + "/" + "Discrete_Fit_" + Path + ".csv" );
	res_csv << PathWithDate;	
	res_csv << ",Frac2GParameter,Frac3GParameter,Frac3GParameterCorr,Frac3GParameternew,DeltaParameter,ChiSquared,Probability_of_the_fit,Rsquared,Adjusted_Rsquared,ChiSquared/DegreesOfFreedom,";
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
		/*if( j < FixedIterator + pPsIndex + 1 && (int)j != pPsIndex )
		{
			res_csv << "Lifetime_for_" << NumberToChar( j+1, 0 )<< "_Component_[ns],";
			res_csv << "Error_of_Lifetime_for_" << NumberToChar( j+1, 0 )<< "_Component_[ns],";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component,";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent,";
			if( Lifetimes[j].Type != "f" )
			{
				res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent_no_fix,";
				res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent_no_fix,";			  
			}
		}
		else if( (int)j != pPsIndex )
		{
			res_csv << "Lifetime_for_" << NumberToChar( j+1, 0 )<< "_Component_[ns],";
			res_csv << "Error_of_Lifetime_for_" << NumberToChar( j+1, 0 )<< "_Component_[ns],";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component,";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent,";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent_no_fix,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent_no_fix,";
		}*/
	}
	res_csv << "\n";
	res_csv << NumberToChar( Frac2GParameter, 5 ) << "," << NumberToChar( Frac3GParameter, 5 ) << NumberToChar( Frac3GParameterNew2, 5 ) << "," << "," << NumberToChar( Frac3GParameterNew, 5 ) << "," << NumberToChar( DeltaParameter, 5 ) << ",";
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
		/*if( j < FixedIterator + pPsIndex + 1 && (int)j != pPsIndex )
		{
			res_csv << NumberToChar( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), 5) << ",";
			res_csv << Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) << ",";
			res_csv << NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ), 5) << ",";
			res_csv << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) << ",";
			res_csv << NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ) * 100, 5 ) << ",";
			res_csv << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) * 100 << ",";
			if( Lifetimes[j].Type != "f" )
			{
				res_csv << NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ) * 100/(1-FixedFixedIntensity), 5 ) << ",";
				res_csv << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) * 100/(1-FixedFixedIntensity) << ",";			  
			}
		}
		else if( (int)j != pPsIndex )
		{
			NotFixedIterator++;
			res_csv << NumberToChar( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), 5) << ",";
			res_csv << Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) << ",";
			res_csv << NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) )/FreeIntensities , 5) << ",";
			res_csv << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) )/FreeIntensities << ",";
			res_csv << NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/FreeIntensities, 5 ) << ",";
			res_csv << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) *(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/FreeIntensities << ",";
			res_csv << NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/(1-FixedFixedIntensity)/FreeIntensities, 5 ) << ",";
			res_csv << Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/FreeIntensities << ",";
		}*/
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
	return 0;
//---------------------------------------
	
//---------------------------------------Residuals	
	delete c1;
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
	Double_t Residuals_root[ Range_To - Range_From ], ResArg_root[ Range_To - Range_From ];
	for( int i = 0; i < histogram->GetNbinsX(); i++ )
	{
		if( histogram->GetBinCenter(i) > Arguments[ Range_From ] && histogram->GetBinCenter(i) < Arguments[ Range_To ] && i >= (int)Range_From )
		{
			if( isnan(histogram->GetBinContent(i) - Discrete->Eval( histogram->GetBinCenter(i) ) )/sqrt( histogram->GetBinContent(i) ) || isinf(histogram->GetBinContent(i) - Discrete->Eval( histogram->GetBinCenter(i) ) )/sqrt( histogram->GetBinContent(i) ) )
				Residuals_root[i-Range_From] = 0;
			else
				Residuals_root[i-Range_From] = (histogram->GetBinContent(i) - Discrete->Eval( histogram->GetBinCenter(i) ) )/sqrt( histogram->GetBinContent(i) );
			ResArg_root[i-Range_From] = histogram->GetBinCenter(i);
		}
	}

	TFile *residuals = new TFile( "Residuals.root", "update" );
	residuals -> mkdir( Path.c_str() ); 
	residuals -> cd( Path.c_str() );  
        TGraph *residualsHisto = new TGraph( Range_To - Range_From, ResArg_root, Residuals_root );
        residualsHisto -> SetMarkerStyle( 21 );
        residualsHisto -> SetMarkerSize( 0.5 );
        residualsHisto -> SetTitle( "Residuals of the fit" );
        residualsHisto -> Draw( "APL" );
        residualsHisto -> GetXaxis() -> SetTitle( "Time difference [ns]" );
        residualsHisto -> GetXaxis() -> SetTitle( "Residuals" );
        residualsHisto -> Write( PathWithDate.c_str() );
        delete residualsHisto;

	residuals->Close();
//---------------------------------------	
	delete histogram;
	delete Discrete;
	delete c2;
	testFile->Close();
//---------------------------------------Saving graphical view of discret results
	
	TFile *comparison = new TFile( "LF_vs_In.root", "update" );
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
	//LF_vs_In -> SetLineColor(0);
	LF_vs_In -> Draw("AP");	
	LF_vs_In -> Write( PathWithDate.c_str() );
	
	delete LF_vs_In;
	comparison->Close();
//---------------------------------------
	return DecoOption;
}

int Fit::Discrete_exp()		//Only for 3 resolution components - quicker because for each resolution size, there is dedacted function that speeds up determination of the parameters
{
  	if( Resolution.size()*Lifetimes.size() == 0 )			//Test if the reading FitDetails file and data file past properly
	{
		return 0;
	}
        std::cout<<"-------------Fitting data-------------"<<std::endl;
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

	std::cout << "[Argument for LastBin] [Argument for FirstBin] \t" << Arguments[LastBin] << "   " << Arguments[FirstBin] << std::endl;

//---------------------------------------Range of drawing histogram with fit on it. first two lines is default, next to lines is for setting by user, if strict range is needed

	float DrawFrom = Arguments[FirstBin];
	float DrawTo = Arguments[LastBin];

//---------------------------------------

	TH1F* histogram = FillHistogram( "Time differences", (DrawTo - DrawFrom)*1/BinWidth, DrawFrom, DrawTo, Times );
	if( NmbOfBins )
	{
		for(unsigned i=FirstBin; i<LastBin+1; i++)
		{
			histogram->SetBinContent(i-FirstBin+1,Times[i]);
		}
	}	

	histogram -> GetXaxis() -> SetTitle( "Time difference [ns]" );
	histogram -> GetYaxis() -> SetTitle( "Counts" );
	histogram -> GetYaxis()-> SetTitleOffset(1.3);
	histogram -> GetXaxis() -> SetTitleFont(62);
	histogram -> GetXaxis() -> SetLabelFont(62);
	histogram -> GetYaxis() -> SetTitleFont(62);
	histogram -> GetYaxis() -> SetLabelFont(62);
	histogram -> SetTitle( ("Time differences_ " + Path).c_str() );
	histogram -> Draw( "" );
	histogram -> SetLineWidth( 2 );
	histogram -> SetLineColor( 1 );

//---------------------------------------Setting the fit function
	int lastbin = LastBin; 
	int firstbin = FirstBin;
	while( histogram -> GetBinCenter( lastbin ) > Arguments[ Range_To ] )
	{
		lastbin--;
	}
	while( histogram -> GetBinCenter( firstbin ) < Arguments[ Range_From ] )
	{
		firstbin++;
	}

	std::cout << "[LastBin - new!, of the DrawRange] [End of Range] - test if they are in proper order (LastBin < EndOfRange)\t" << lastbin << "   " << Range_To << std::endl;
	double area = 0.013*(histogram -> Integral( firstbin, lastbin ) - Background*( Arguments[ lastbin ] - Arguments[ firstbin ] ) );
	
	std::cout << "Area \t \t \t \t \t \t " << area << std::endl;
	area = fabs(area);
	
	std::cout << "Range - Arguments [Start] [End] \t \t" << Arguments[ Range_From ] << "   " << Arguments[ Range_To ] << std::endl;
	std::cout << "Range - Values [Start] [End] \t \t \t" << Values[ Range_From ] << "   " << Values[ Range_To ] << std::endl;
	
	std::cout <<"-------------Old fit way-------------"<< std::endl;
	
	int pPsIndex = -1;	
	for( unsigned i=0; i<Lifetimes.size(); i++ )
	{
		if( Lifetimes[i].Type == "ps" )
			pPsIndex = 0;
	}
	TF1 *Discrete;
	if( pPsIndex + 1 )
	{
		if( Resolution.size() == 1 )
			Discrete = new TF1( "Discrete", DiscreteFitFunctionPs_exp_1, Arguments[ Range_From ], Arguments[ Range_To ], 4 + Lifetimes.size()*2 + 3*Resolution.size() + 1 );
		else if( Resolution.size() == 2 )
			Discrete = new TF1( "Discrete", DiscreteFitFunctionPs_exp_2, Arguments[ Range_From ], Arguments[ Range_To ], 4 + Lifetimes.size()*2 + 3*Resolution.size() + 1 );
		else
			Discrete = new TF1( "Discrete", DiscreteFitFunctionPs_exp_3, Arguments[ Range_From ], Arguments[ Range_To ], 4 + Lifetimes.size()*2 + 3*Resolution.size() + 1 );
	}
	else
	{
		if( Resolution.size() == 1 )
			Discrete = new TF1( "Discrete", DiscreteFitFunctionNoPs_exp_1, Arguments[ Range_From ], Arguments[ Range_To ], 4 + Lifetimes.size()*2 + 3*Resolution.size() + 1 );
		else if( Resolution.size() == 2 )
			Discrete = new TF1( "Discrete", DiscreteFitFunctionNoPs_exp_2, Arguments[ Range_From ], Arguments[ Range_To ], 4 + Lifetimes.size()*2 + 3*Resolution.size() + 1 );
		else
			Discrete = new TF1( "Discrete", DiscreteFitFunctionNoPs_exp_3, Arguments[ Range_From ], Arguments[ Range_To ], 4 + Lifetimes.size()*2 + 3*Resolution.size() + 1 );
	}
	
	
	Discrete -> SetNpx(Range_To - Range_From);
	/*Nmbr of arguments - 
	1x Background
	1x nmbr of gauss
	1x nmbr of lf comp
	1x total area
	Lifetimes.size() - 1x Intensity, 1x Lifetime
	Resolution.size() - 1x Sigma, 1x Offset, 1x Fraction
	1x fixed components
	*/
 	Discrete -> SetParameter( 0, Background );		//Background
	Discrete -> SetParName( 0, "Background" );
	Discrete -> SetParLimits( 0, Background - 3*SDBackground, Background + 3*SDBackground );
	Discrete -> FixParameter( 1, Resolution.size() );			//Nmbr of Gauss comp
	Discrete -> SetParName( 1, "Number of Gauss resolution components" );
	Discrete -> FixParameter( 2, Lifetimes.size() );			//Nmbr of Lifetime comp
	Discrete -> SetParName( 2, "Number of Lifetimes components" );
	Discrete -> SetParameter( 3, area );			//Tot area of components
	Discrete -> SetParName( 3, "Total area of components" );

	for( unsigned i = 0; i < Resolution.size(); i++ )
	{
		if( FixGauss == "yes" )
		{
			Discrete -> FixParameter( 4 + 3*i, Resolution[i].Lifetime );
			Discrete -> SetParName( 4 + 3*i, ("Sigma for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
			
			Discrete -> FixParameter( 5 + 3*i, SetIntensityParameter( Resolution, i+1 ) );
			Discrete -> SetParName( 5 + 3*i, ("Fraction parameter for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
			
			Discrete -> FixParameter( 6 + 3*i, Arguments[ BinMax ] + i*0.33 );
			Discrete -> SetParName( 6 + 3*i, ("Offset for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
		}
		else
		{
			Discrete -> SetParameter( 4 + 3*i, Resolution[i].Lifetime );
			Discrete -> SetParName( 4 + 3*i, ("Sigma for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
			Discrete -> SetParLimits( 4 + 3*i, 0.01, 5 );
			
			Discrete -> SetParameter( 5 + 3*i, SetIntensityParameter( Resolution, i+1 ) );
			Discrete -> SetParName( 5 + 3*i, ("Fraction parameter for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
			if( i == 0 )
				Discrete -> SetParLimits( 5 + 3*i, 0, 1 );
			else
				Discrete -> SetParLimits( 5 + 3*i, 0, 3.14159/2 );
			
			Discrete -> SetParameter( 6 + 3*i, Arguments[ BinMax ] + i*0.33 );
			Discrete -> SetParName( 6 + 3*i, ("Offset for " + NumberToChar( i+1, 0 ) + " Gauss").c_str() );
			Discrete -> SetParLimits( 6 + 3*i, Arguments[ BinMax ] - 0.5, Arguments[ BinMax ] + 1 );
		}
	}
	unsigned NotFixedIterator = 0, FixedIterator = 0;
	for( unsigned j = 0; j < Lifetimes.size(); j++ )
	{
		Discrete -> SetParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
 		Discrete -> SetParName( 4 + 3*Resolution.size() + 2*j, ("Lifetime for " + NumberToChar( j+1, 0 ) + " Component").c_str() );

/*		
ExMoGa( Double_t A, Double_t P1, Double_t P2, Double_t P3, Double_t P4 )			//A - arguments, p - parameters
P4 / (2*P3) * exp( P1*P1 / (2*P3*P3) - ( A - P2 ) / P3 ) * ( erf( 1 / (sqrt(2)*P1) * ( A - P1*P1 / P3 - P2 ) ) - erf( 1 / (sqrt(2)*P1) * ( - P1*P1 / P3 - P2 ) ) );*/
		
		Discrete -> SetParName( 5 + 3*Resolution.size() + 2*j, ("Intensity for " + NumberToChar( j+1, 0 ) + " Component").c_str() );
		
		if( Lifetimes[j].Type == "f" || Lifetimes[j].Type == "pf" || Lifetimes[j].Type == "lpf" || Lifetimes[j].Type == "lf" || Lifetimes[j].Type == "lff" )
		{
			Discrete -> FixParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
			Discrete -> FixParameter( 5 + 3*Resolution.size() + 2*j, Lifetimes[j].Intensity );
			FixedIterator++;
		}
		else if( Lifetimes[j].Type == "ps" )
		{
			Discrete -> FixParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
			Discrete -> FixParameter( 5 + 3*Resolution.size() + 2*j, Lifetimes[j].Intensity );
			pPsIndex = j;
		}
		else
		{
			NotFixedIterator++;
			Discrete -> SetParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
			Discrete -> SetParLimits( 4 + 3*Resolution.size() + 2*j, 0.1, 142 );
			Discrete -> SetParameter( 5 + 3*Resolution.size() + 2*j, Lifetimes[j].Intensity );
			Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, 0, 1 );
		}
	}
	Discrete -> FixParameter( 6 + 3*Resolution.size() + 2*(Lifetimes.size()-1), FixedIterator );
 	Discrete -> SetParName( 6 + 3*Resolution.size() + 2*(Lifetimes.size()-1), "Number of somehow fixed Components not ps type" );
        histogram -> Fit(Discrete,"RM");
	std::cout << std::endl;
	std::cout <<"Do not look on intensities of completely free parameters and resolution intensities"<< std::endl;
	std::cout <<"They are normalized in the body of fitting function and at the end when saving results"<< std::endl;
	std::cout << std::endl;
//---------------------------------------
//---------------------------------------Iterations	
	for( unsigned iter = 0; iter < NmbrOfIterations; iter++ )
	{
		  
		Discrete -> SetParameter( 0, Discrete -> GetParameter(0) );		//Background
		Discrete -> SetParLimits( 0, Background - 3*SDBackground, Background + 3*SDBackground );
		Discrete -> FixParameter( 1, Resolution.size() );			//Nmbr of Gauss comp
		Discrete -> FixParameter( 2, Lifetimes.size() );			//Nmbr of Lifetime comp
		Discrete -> SetParameter( 3, Discrete -> GetParameter(3) );

		for( unsigned i = 0; i < Resolution.size(); i++ )
		{	
			if( FixGauss == "yes" )
			{
				Discrete -> FixParameter( 4 + 3*i, Discrete -> GetParameter(4 + 3*i) );

				Discrete -> FixParameter( 5 + 3*i, Discrete -> GetParameter(5 + 3*i) );

				Discrete -> FixParameter( 6 + 3*i, Discrete -> GetParameter(6 + 3*i) );
			}
			else
			{
				Discrete -> SetParameter( 4 + 3*i, Discrete -> GetParameter(4 + 3*i) );
				Discrete -> SetParLimits( 4 + 3*i, 0.01, 5 );

				Discrete -> SetParameter( 5 + 3*i, Discrete -> GetParameter(5 + 3*i) );

				Discrete -> SetParameter( 6 + 3*i, Discrete -> GetParameter(6 + 3*i) );
			}
		}
		NotFixedIterator = 0;
		for( unsigned j = 0; j < Lifetimes.size(); j++ )
		{
			Discrete -> SetParameter( 4 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) );

			Discrete -> SetParameter( 5 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) );
			if( Lifetimes[j].Type == "f" )
			{
				Discrete -> FixParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
				Discrete -> FixParameter( 5 + 3*Resolution.size() + 2*j, Lifetimes[j].Intensity );
			}
			else if( Lifetimes[j].Type == "pf" )
			{
				Discrete -> SetParLimits( 4 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) - VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j), Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) + VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) );
				Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, 0, 1 );			
			}
			else if( Lifetimes[j].Type == "lpf" )
			{
				Discrete -> SetParLimits( 4 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) - VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j), Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) + VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) );
				Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) - VarLvl*Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j), Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) + VarLvl*Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) );
			}
			else if( Lifetimes[j].Type == "lf" )
			{
				Discrete -> FixParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
				Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, 0, 1 );
			}
			else if( Lifetimes[j].Type == "lff" )
			{
				Discrete -> FixParameter( 4 + 3*Resolution.size() + 2*j, Lifetimes[j].Lifetime );
				Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) - VarLvl*Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j), Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) + VarLvl*Discrete -> GetParameter(5 + 3*Resolution.size() + 2*j) );
			}
			else if( Lifetimes[j].Type == "ps" )
			{
				Discrete -> SetParLimits( 4 + 3*Resolution.size() + 2*j, Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) - VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j), Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) + VarLvl*Discrete -> GetParameter(4 + 3*Resolution.size() + 2*j) );
				Discrete -> FixParameter( 5 + 3*Resolution.size() + 2*j, Lifetimes[j].Intensity );		  
			}
			else
			{
				NotFixedIterator++;
				Discrete -> SetParLimits( 4 + 3*Resolution.size() + 2*j, 0.01, 142 );
			  	Discrete -> SetParLimits( 5 + 3*Resolution.size() + 2*j, 0, 1 );
			}
		}
		Discrete -> FixParameter( 6 + 3*Resolution.size() + 2*(Lifetimes.size()-1), FixedIterator );
		
		histogram -> Fit(Discrete,"RM");
		std::cout << std::endl;
		std::cout <<"Do not look on intensities of completely free parameters and resolution intensities"<< std::endl;
		std::cout <<"They are normalized in the body of fitting function and at the end when saving results"<< std::endl;
		std::cout << std::endl;
	}
//---------------------------------------

//---------------------------------------Drawing and saving histograms and results
        Discrete->Draw("C same");
	TFile *testFile = new TFile( "Results_from_fit.root", "update" );
        testFile -> mkdir( Path.c_str() );
        testFile -> cd( Path.c_str() );
	std::string Res_path = "Results";
		
	if( not PathCheck( Res_path ) )
        {
                std::cout << "Creating directory for results" << std::endl;
                boost::filesystem::path dir( Res_path );
                boost::filesystem::create_directories( dir );
        }
	
        c1->SaveAs( ( Res_path + "/" + PathWithDate + ".png" ).c_str() );
        c1->SetLogy();
        c1->SaveAs( ( Res_path + "/" + PathWithDate + "_logscale"+".png" ).c_str() );
	c1 -> Write( PathWithDate.c_str() );

	AreaFromDiscrete = Discrete -> GetParameter( 3 );
	Background = Discrete -> GetParameter( 0 );
	Double_t ResolutionsFromFit[22];
	ResolutionsFromFit[0] = 1;
	ResolutionsFromFit[1] = Resolution.size();
	Double_t ResolutionsFromFitErrors[22];
	ResolutionsFromFitErrors[0] = Resolution.size();
	double FreeIntensities = 0., FixedIntensities = 0., pPsIntensity = 0., FreeIntensitiesTopPs = 0., FixedFixedIntensity = 0.;
	for( unsigned i = 0; i < Resolution.size(); i++ )
	{
		ResolutionsFromFit[2*(i+1)] = Discrete -> GetParameter( 4 + 3*i );
		ResolutionsFromFit[2*(i+1)+1] = Discrete -> GetParameter( 5 + 3*i );
		ResolutionsFromFitErrors[2*i+1] = Discrete -> GetParameter( 5 + 3*i );
		ResolutionsFromFitErrors[2*i+2] = Discrete -> GetParError( 5 + 3*i );
	}
	for( unsigned i = 0; i < Lifetimes.size(); i++ )
	{
		if( Lifetimes[i].Type != "nf" && Lifetimes[i].Type != "ps" )
		{
			if( Lifetimes[i].Type == "f" )
				FixedFixedIntensity += Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
			FixedIntensities += Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
			if( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*i ) > 0.7 && pPsIndex >=0 )
				pPsIntensity += Discrete -> GetParameter( 5 + 3*Resolution.size() )*Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
		}
		else if( Lifetimes[i].Type != "ps" )
		{
			FreeIntensities += Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
			if( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*i ) > 0.7 && pPsIndex >=0 )
				FreeIntensitiesTopPs += Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
		}
	}
	FixedIntensities += pPsIntensity;
	if( pPsIndex >=0 )
	{
		pPsIntensity += Discrete -> GetParameter( 5 + 3*Resolution.size() )*FreeIntensitiesTopPs*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) );
	}
	
	std::ofstream res;
	res.open( Res_path + "/" + "Discrete_Fit_" + PathWithDate );
	res << "Results from fitting" << std::endl;
	res << "Parameter name\t \t \t \t \t" << "Value\t \t \t" << "Error" << std::endl;
	res << std::endl;
	res << std::endl;
	res << "Background\t \t \t \t \t" << NumberToChar( Background, 5 ) << "\t \t" << SDBackground << std::endl;
	for( unsigned i = 0; i < Resolution.size(); i++ )
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
		res << "Lifetime for p-PS Component\t \t \t " << NumberToChar( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), 5) << "\t \t" << Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) << std::endl;
		LifetimesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) ) );
		res << "Intensity for p-PS Component \t \t \t " <<  NumberToChar( pPsIntensity, 5) << std::endl;
		res << "Intensity for p-PS Component in percent \t " <<  NumberToChar( pPsIntensity * 100, 5) << std::endl;
		res << "Intensity for p-PS Component in percent no fix" << "\t " <<  NumberToChar( pPsIntensity * 100/(1-FixedFixedIntensity), 5 ) << std::endl;	
		res << "Intensity for p-PS Component in percent of o-Ps\t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() ) * 100, 5)  << std::endl;
		IntensitiesFromDiscrete.push_back( DiscreteFitResult( pPsIntensity, 0 ) );
		res << std::endl;
	}
	for( unsigned j = 0; j < Lifetimes.size(); j++ )
	{
		if( j < FixedIterator + pPsIndex + 1 && (int)j != pPsIndex )
		{
			res << ("Lifetime for " + NumberToChar( j+1, 0 ) + " Component [ns]").c_str()
				  << "\t \t \t " << NumberToChar( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), 5) << "\t \t" << Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) << std::endl;
			LifetimesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) ) );
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component").c_str() 
				  << "\t \t \t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ), 5) << "\t \t" << 
					        Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j )  << std::endl;
			IntensitiesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ), Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) ) );
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component in percent").c_str()
		                  << "\t \t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ) * 100, 5 )<< "\t \t" <<
		                               Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) * 100 << std::endl;
			if( Lifetimes[j].Type != "f" )
			{
				  res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component in percent no fix").c_str()
		                  << "\t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ) * 100/(1-FixedFixedIntensity), 5 )<< "\t \t" <<
		                               Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) * 100/(1-FixedFixedIntensity) << std::endl;
			}
			else
			{
				LifetimesFromDiscrete[ LifetimesFromDiscrete.size() - 1 ].Type = "f";
			}
		        std::cout << "Intensity in percent: \t" << Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ) * 100 << std::endl;
			res << std::endl;
		}
		else if( (int)j != pPsIndex )
		{
			res << ("Lifetime for " + NumberToChar( j+1, 0 ) + " Component [ns]").c_str()
				  << "\t \t \t " << NumberToChar( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), 5) << "\t \t" << Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) << std::endl;
			LifetimesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) ) );
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component").c_str() 
				  << "\t \t \t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) )/FreeIntensities , 5) << "\t \t" << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) )/FreeIntensities << std::endl;
			IntensitiesFromDiscrete.push_back( DiscreteFitResult( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) )/FreeIntensities, Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) )/FreeIntensities ) );
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component in percent").c_str()
		                  << "\t \t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/FreeIntensities, 5 ) << "\t \t" << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) *(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/FreeIntensities<< std::endl;
			res << ("Intensity for " + NumberToChar( j+1, 0 ) + " Component in percent no fix").c_str()
		                  << "\t " <<  NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/(1-FixedFixedIntensity)/FreeIntensities, 5 ) << "\t \t" << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/(1-FixedFixedIntensity)/FreeIntensities << std::endl;	  
		        std::cout << "Intensity in percent: \t" << Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/FreeIntensities << std::endl;
			res << std::endl;
		}
	}

//---------------------------------------Sorting the lifetime components, so the first component has the lowest lifetime, second, the second lowest, ...
	
	std::vector<unsigned> Order;
	double minTemp = 200, minPrevious = 0;
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
		minTemp = 200;
		if( Order.size() == LifetimesFromDiscrete.size() )
			minSearch = false;
	}
//---------------------------------------	

	double Frac2GParameter = 0, Frac2GParameterNew = 0, Frac3GParameter = 0, Frac3GParameterNew = 0, Frac3GParameterNew2 = 0, DeltaParameter = 0;
	if( LifetimesFromDiscrete.size() == 2 )
	{
		Frac2GParameter = 371/372;
		Frac2GParameterNew = 371/372;
		Frac3GParameter = 1/372;
		Frac3GParameterNew = 1/372;
		Frac3GParameterNew2 = 1/372;
	}
	else
	{
		for( unsigned i=1; i<LifetimesFromDiscrete.size(); i++ )
		{
			if( LifetimesFromDiscrete[ Order[i] ].Parameter < 0.8 )
			{
				Frac2GParameter += 371*IntensitiesFromDiscrete[ Order[i] ].Parameter / 372;
				Frac2GParameterNew += 371*IntensitiesFromDiscrete[ Order[i] ].Parameter / 372;
				Frac3GParameter += IntensitiesFromDiscrete[ Order[i] ].Parameter / 372;
				Frac3GParameterNew += IntensitiesFromDiscrete[ Order[i] ].Parameter / 372;
			}
			else
			{
				Frac2GParameter += IntensitiesFromDiscrete[ Order[i] ].Parameter  / 3;
				Frac2GParameter += IntensitiesFromDiscrete[ Order[i] ].Parameter * 142 / ( 142 - LifetimesFromDiscrete[ Order[i] ].Parameter );
				Frac2GParameterNew += IntensitiesFromDiscrete[ Order[i] ].Parameter  / 3;
				Frac2GParameterNew += IntensitiesFromDiscrete[ Order[i] ].Parameter * ( 142 - LifetimesFromDiscrete[ Order[i] ].Parameter ) / 142;
				Frac3GParameter += IntensitiesFromDiscrete[ Order[i] ].Parameter * LifetimesFromDiscrete[ Order[i] ].Parameter / ( 142 - LifetimesFromDiscrete[ Order[i] ].Parameter );
				Frac3GParameterNew2 += IntensitiesFromDiscrete[ Order[i] ].Parameter * LifetimesFromDiscrete[ Order[i] ].Parameter / ( 142 - LifetimesFromDiscrete[ Order[i] ].Parameter );
				Frac3GParameterNew += IntensitiesFromDiscrete[ Order[i] ].Parameter * LifetimesFromDiscrete[ Order[i] ].Parameter / 142;
				DeltaParameter += (1-4*IntensitiesFromDiscrete[ Order[i] ].Parameter/3)/372 + IntensitiesFromDiscrete[ Order[i] ].Parameter * LifetimesFromDiscrete[ Order[i] ].Parameter / 142;
			}
		}
	}	
	res << "Frac2GParameter \t" << Frac2GParameter << std::endl;
	res << "Frac2GParameterNew \t" << Frac2GParameterNew << std::endl;	
	res << "Frac3GParameter \t" << Frac3GParameter << std::endl;
	res << "Frac3GParameterCorr \t" << Frac3GParameterNew2 << std::endl;
	res << "Frac3GParameternew \t" << Frac3GParameterNew << std::endl;
	res << "DeltaParameter \t" << DeltaParameter << std::endl;	
	
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

	res << std::endl;
	res << "ChiSquared" << "\t \t \t \t \t" << NumberToChar( Discrete -> GetChisquare(), 4 ) << std::endl;
	res << "Probability of the fit " << "\t \t \t \t" << NumberToChar(1 - Discrete -> GetProb(), 4) << std::endl;
	res << "Rsquared" << "\t \t \t \t \t" << NumberToChar(1 - (Discrete -> GetChisquare())/Diff, 4 ) << std::endl;
	res << "Adjusted Rsquared" << "\t \t \t \t" << NumberToChar(1 - ((sumIT-1)/Discrete -> GetNDF())*(Discrete -> GetChisquare())/Diff, 4 ) << std::endl;
	res << "ChiSquared/DegreesOfFreedom" << "\t \t \t" << NumberToChar(  Discrete -> GetChisquare()/Discrete -> GetNDF(), 4 ) << std::endl;
	res.close();
	
//---------------------------------------Writing to csv	
	bool IsFileExists =  FileCheck( Res_path + "/" + "Discrete_Fit_" + Path + ".csv" );
	std::ofstream res_csv;
	if( IsFileExists )
		res_csv.open( Res_path + "/" + "Discrete_Fit_" + Path + ".csv", std::ofstream::app );
	else
		res_csv.open( Res_path + "/" + "Discrete_Fit_" + Path + ".csv" );
	res_csv << PathWithDate;	
	res_csv << ",Frac2GParameter,Frac3GParameter,Frac3GParameterCorr,Frac3GParameternew,DeltaParameter,ChiSquared,Probability_of_the_fit,Rsquared,Adjusted_Rsquared,ChiSquared/DegreesOfFreedom,";
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
		/*if( j < FixedIterator + pPsIndex + 1 && (int)j != pPsIndex )
		{
			res_csv << "Lifetime_for_" << NumberToChar( j+1, 0 )<< "_Component_[ns],";
			res_csv << "Error_of_Lifetime_for_" << NumberToChar( j+1, 0 )<< "_Component_[ns],";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component,";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent,";
			if( Lifetimes[j].Type != "f" )
			{
				res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent_no_fix,";
				res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent_no_fix,";			  
			}
		}
		else if( (int)j != pPsIndex )
		{
			res_csv << "Lifetime_for_" << NumberToChar( j+1, 0 )<< "_Component_[ns],";
			res_csv << "Error_of_Lifetime_for_" << NumberToChar( j+1, 0 )<< "_Component_[ns],";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component,";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent,";
			res_csv << "Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent_no_fix,";
			res_csv << "Error_of_Intensity_for_" << NumberToChar( j+1, 0 )<< "_Component_in_percent_no_fix,";
		}*/
	}
	res_csv << "\n";
	res_csv << NumberToChar( Frac2GParameter, 5 ) << "," << NumberToChar( Frac3GParameter, 5 ) << NumberToChar( Frac3GParameterNew2, 5 ) << "," << "," << NumberToChar( Frac3GParameterNew, 5 ) << "," << NumberToChar( DeltaParameter, 5 ) << ",";
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
		/*if( j < FixedIterator + pPsIndex + 1 && (int)j != pPsIndex )
		{
			res_csv << NumberToChar( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), 5) << ",";
			res_csv << Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) << ",";
			res_csv << NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ), 5) << ",";
			res_csv << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) << ",";
			res_csv << NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ) * 100, 5 ) << ",";
			res_csv << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) * 100 << ",";
			if( Lifetimes[j].Type != "f" )
			{
				res_csv << NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j ) * 100/(1-FixedFixedIntensity), 5 ) << ",";
				res_csv << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) * 100/(1-FixedFixedIntensity) << ",";			  
			}
		}
		else if( (int)j != pPsIndex )
		{
			NotFixedIterator++;
			res_csv << NumberToChar( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*j ), 5) << ",";
			res_csv << Discrete -> GetParError( 4 + 3*Resolution.size() + 2*j ) << ",";
			res_csv << NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) )/FreeIntensities , 5) << ",";
			res_csv << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) )/FreeIntensities << ",";
			res_csv << NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/FreeIntensities, 5 ) << ",";
			res_csv << Discrete -> GetParError( 5 + 3*Resolution.size() + 2*j ) *(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/FreeIntensities << ",";
			res_csv << NumberToChar( Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/(1-FixedFixedIntensity)/FreeIntensities, 5 ) << ",";
			res_csv << Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*j )*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) ) * 100/FreeIntensities << ",";
		}*/
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
	return 0;
//---------------------------------------
	
//---------------------------------------Residuals	
	delete c1;
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
	Double_t Residuals_root[ Range_To - Range_From ], ResArg_root[ Range_To - Range_From ];
	for( int i = 0; i < histogram->GetNbinsX(); i++ )
	{
		if( histogram->GetBinCenter(i) > Arguments[ Range_From ] && histogram->GetBinCenter(i) < Arguments[ Range_To ] && i >= (int)Range_From )
		{
			if( isnan(histogram->GetBinContent(i) - Discrete->Eval( histogram->GetBinCenter(i) ) )/sqrt( histogram->GetBinContent(i) ) || isinf(histogram->GetBinContent(i) - Discrete->Eval( histogram->GetBinCenter(i) ) )/sqrt( histogram->GetBinContent(i) ) )
				Residuals_root[i-Range_From] = 0;
			else
				Residuals_root[i-Range_From] = (histogram->GetBinContent(i) - Discrete->Eval( histogram->GetBinCenter(i) ) )/sqrt( histogram->GetBinContent(i) );
			ResArg_root[i-Range_From] = histogram->GetBinCenter(i);
		}
	}

	TFile *residuals = new TFile( "Residuals.root", "update" );
	residuals -> mkdir( Path.c_str() ); 
	residuals -> cd( Path.c_str() );  
        TGraph *residualsHisto = new TGraph( Range_To - Range_From, ResArg_root, Residuals_root );
        residualsHisto -> SetMarkerStyle( 21 );
        residualsHisto -> SetMarkerSize( 0.5 );
        residualsHisto -> SetTitle( "Residuals of the fit" );
        residualsHisto -> Draw( "APL" );
        residualsHisto -> GetXaxis() -> SetTitle( "Time difference [ns]" );
        residualsHisto -> GetXaxis() -> SetTitle( "Residuals" );
        residualsHisto -> Write( PathWithDate.c_str() );
        delete residualsHisto;

	residuals->Close();
//---------------------------------------	
	delete histogram;
	delete Discrete;
	delete c2;
	testFile->Close();
//---------------------------------------Saving graphical view of discret results
	
	TFile *comparison = new TFile( "LF_vs_In.root", "update" );
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
	//LF_vs_In -> SetLineColor(0);
	LF_vs_In -> Draw("AP");	
	LF_vs_In -> Write( PathWithDate.c_str() );
	
	delete LF_vs_In;
	comparison->Close();	
//---------------------------------------
	return DecoOption;
}

void Fit::Deconvolution()
{
	if( Resolution.size()*Lifetimes.size() == 0 )				//Test if the reading FitDetails file and data file past properly
	{
		return;
	}
	
	std::cout << "-------------Initializing LifetimeGrid-------------" << std::endl;
	double MaxLF = FindMaxLF();		//Looking for minimal and maximal value of lifetimes from Discrete procedure of fitting
	double MinLF = FindMinLF();		//in order to describe whole LifetimeGrid as the range between FracMinLF*MinLF and FracMaxLF*MaxLF
	
	int StartPoint = 0;			//Iterator that characterizes start point of LifetimeGrid
	int cond = 1;				//Condition that control process of finding LifetimeGrid - is zero if whole LifetimeGrid is as user wanted
	unsigned SizeOfLFGrid = 1;		//Variable that stores the information about size of LifetimeGrid
	
	while( cond )
	{
		if( exp((double)( StartPoint )/Scaler)-1 > FracMinLF*MinLF )		//If Value of StartPoint will be such, that it will be greater than lower side of range of LifetimeGrid, then beginning of LifetimeGrid is found
		{
			if( exp((double)( StartPoint + SizeOfLFGrid )/Scaler)-1 < FracMaxLF*MaxLF )	//Looking next for the end of LifetimeGrid as the minmal point that fulfills condition  in if()
				SizeOfLFGrid++;
			else
				cond = 0;
		}
		else
			StartPoint++;
	}
//-----------------------------------------------
	using namespace arma;			//Creating vectors for continous fitting procedure - It uses Armadillo software -> arma is the namespace for them
	vec LFGrid_(SizeOfLFGrid);		//Definiton of mock vector that will be next written to Fit class vector without "_"
	vec Inten(SizeOfLFGrid);		//Vector of Intensities of each Lifetime in LifetimeGrid
	Double_t LifetimeGrid_root[SizeOfLFGrid], Intensities_root[SizeOfLFGrid];	//This Vectors are introduced in order to plot LifetimeGrid and it intesities in root Graph
	for( unsigned i=0; i<SizeOfLFGrid; i++ )
	{
		LifetimeGrid.push_back( exp((double)(i+StartPoint)/Scaler)-1 );	
		LFGrid_(i) = exp((double)(i+StartPoint)/Scaler)-1 ;
	}
	LFGrid = LFGrid_;
//-----------------------------------------------
	std::cout << "-------------Smoothing histogram-------------" << std::endl;
	LinearFilter( LinFilterRange );				//Smoothing the histogram more - in order to obtain better distribution of lifetimes
//-----------------------------------------------	
	std::vector < double > TempArg, TempVal;		//Temporary vectors for trimming Arguments and Values vectors, only to real range of fitting
	auto maxValueIterator = std::distance(Values.begin(), std::max_element(Values.begin()+Range_From, Values.end() - ShiftForBackgroundEstimation)); 	//Finding maxmimum of the distribution to fit only from the maximum to the Range_To (To avoid resolution uncertainity)
	
	for( unsigned i=0; i<Arguments.size() - ShiftForBackgroundEstimation; i++ )
	{
		if( i >= maxValueIterator + MaxShiftForDeco && Arguments[i] <= EndOfFitMultiplicity*MaxLF  )
		{	
			TempArg.push_back( Arguments[i] );
			TempVal.push_back( Values[i] );		
		}
	}
	Values = TempVal;			//Setting Arguments and Values vectors as trimmed vectors only to interesting range, where user wants to fit
	Arguments = TempArg;
//-----------------------------------------------	
	std::cout << "-------------Initializing gradient base for continous fitting-------------" << std::endl;
	double sum = 0.;			//sum is the temporary variable for calculations of gradient matrix components
	mat GradMax_( LifetimeGrid.size(), Arguments.size() );		//As in LFGrid_ situation, creating temporary matrix for Jacobian matrix of fitting function (It should have )
	Double_t par[4];			//Parameter container created such, that function defined to fit in discrete mode is compatible with continous function
						//0 -> X of histogram; 1 -> Sigma of resolution component; 2 -> Offset of resolution component; 3 -> Lifetime
	for( unsigned i=0; i<LifetimeGrid.size(); i++ )			//Gradient calculated for all Lifetimes in LifetiemsGrid
	{
		Intensities_afterIter.push_back( LifetimeGrid[i] );	//First filling of container for Intensities (Y of Lifetime distribution)
		Inten(i) = LifetimeGrid[i];				//First filling of container for Lifetimes (X of Lifetime distribution)
		for( unsigned j=0; j<Arguments.size(); j++ )		//Jacobian calculated for all Arguments in histogram (cutted only up to Fit Range) - Jacobian dimension -> LifetimeGrid x Arguments
		{
			sum = 0;					//sum over all Resolution components
			for( unsigned k=0; k<ResolutionFromDiscrete.size(); k++ )
			{
				par[0] = Arguments[j];
				par[1] = ResolutionFromDiscrete[k].Parameter;
				par[2] = OffsetsFromDiscrete[k].Parameter;
				par[3] = LifetimeGrid[i];
				//sum += AreaFromDiscrete*FractionsFromDiscrete[k].Parameter*ExMoGa( par[0], par[1], par[2], par[3], 1 );		//Derrivative of Exponentially modified Gaussian over Intensity is just ExMoGa with intensity 1
				sum += FractionsFromDiscrete[k].Parameter*ExMoGa( par[0], par[1], par[2], par[3], 1 );		//Derrivative of Exponentially modified Gaussian over Intensity is just ExMoGa with intensity 1
				if( isnan(sum) )			//Sometimes values are close to zero -> error while calculating returning nan
					sum = 0;
			}
			GradMax_( i,j ) = sum;
		}
	}
	GradMax = GradMax_;
						//Preparing first step of continous fit -> Assuming that Lifetime distribution is a sum of Gaussian/LogGaussians components -> Finding optimal sigmas
	SSig.clear();				//Clearing container for sigmas, that will vary in the first step of continous fitting
	for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )
	{
		SSig.push_back( LifetimesFromDiscrete[j].Parameter*SigmasDefaultFraction );		//First values of sigmas 
	}
	ChiDiffBetweenIterations = 0;
	
	AreaParameter = 1;
	AreaParameter = FindingOptimalAreaParameter( SSig );
	if( ShapeOfComponent == "Gaussian" )
		std::cout << "PreIteration status " <<  PreIterateGauss( SSig, StepForPreIteration, 100, 0. ) << std::endl;	//Finding optimal sigmas -> Measure of goddnes of fit -> ChiSq_pre
	else if( ShapeOfComponent == "LogGaussian" )								//1 argument -> Container for sigmas; 2 -> step in changes of sigmas; 3 -> Lambda which controls the impact of Gauss-Newton method (small lambda) and max damping method (big lambda); 4 -> Previous Chi
		std::cout << "PreIteration status " <<  PreIterateLogGauss( SSig, StepForPreIteration, 100, 0. ) << std::endl;
	else
		std::cout << "PreIteration status " <<  PreIterateMixed( SSig, StepForPreIteration, 100, 0. ) << std::endl;
	  
	ChiDiffBetweenIterations = 0.;

	iterator = 0;			//Setting global iterator that is counting number of iteration/preiteration to 0, before starting second part of continous fitting
	double ChiSq = -0.01;			//Second step is done until condition on Chisquared is fulfilled -> Then the Iterate procedure is returning notzero value which stops this loop
	while( ChiSq <= 0 )
	{
		for( unsigned p=0; p<LifetimeGrid.size(); p++ )
			Inten(p) = Intensities_afterIter[p];			//Copying intensities from preiterations to new vector used in Iterate procedure
		ChiDiffBetweenIterations = 0.;
		if( TypeOfContinousFit )
			ChiSq = IterateExp( Inten, StepForIteration, -ChiSq, 0. );
		else
			ChiSq = Iterate( Inten, StepForIteration, -ChiSq, 0. );
	}
	ChiSq = -0.01;
	for( unsigned p=0; p<LifetimeGrid.size(); p++ )
	{
		if( LifetimeGrid[p] < 5*MaxLF && LifetimeGrid[p] > 0.5*MinLF )
			Inten(p) = Intensities_afterIter[p];			//Copying intensities from preiterations to new vector used in Iterate procedure
		else
			Inten(p) = 0;
	}
	while( ChiSq <= 0 )
	{
		for( unsigned p=0; p<LifetimeGrid.size(); p++ )
			Inten(p) = Intensities_afterIter[p];			//Copying intensities from preiterations to new vector used in Iterate procedure*/
		/*for( unsigned p=0; p<LifetimeGrid.size(); p++ )
		{
			if( LifetimeGrid[p] < 5*MaxLF && LifetimeGrid[p] > 0.5*MinLF )
				Inten(p) = Intensities_afterIter[p];			//Copying intensities from preiterations to new vector used in Iterate procedure
			else
				Inten(p) = 0;
		}*/
		ChiDiffBetweenIterations = 0.;
		if( TypeOfContinousFit )
			ChiSq = IterateExp2( Inten, StepForIteration, -ChiSq, 0. );
		else
			ChiSq = Iterate( Inten, StepForIteration, -ChiSq, 0. );
	}
	for( unsigned k=0; k<SizeOfLFGrid; k++ )
	{
		LifetimeGrid_root[k] = LifetimeGrid[k];
	//	if( LifetimeGrid[k] < 5*MaxLF && LifetimeGrid[k] > 0.5*MinLF )
	//	{
			Intensities_root[k] = Intensities_afterIter[k];
			Inten(k) = Intensities_afterIter[k];
	/*	}
		else
		{
			Intensities_root[k] = 0;
			Inten(k) = 0;		  
		}*/
	}		
	TFile *testFile = new TFile( "LFdistr.root", "update" );
        testFile -> mkdir( Path.c_str() );
        testFile -> cd( Path.c_str() );
        
        TGraph *distr = new TGraph( LifetimeGrid.size(), LifetimeGrid_root, Intensities_root );

        distr -> SetMarkerStyle( 21 );
        distr -> SetMarkerSize( 0.5 );
        distr -> SetTitle( "LF distribution" );
        distr -> Draw( "APL" );

        distr -> GetXaxis() -> SetTitle( "Lifetime [ns]" );

        distr -> Write( PathWithDate.c_str() );

        delete distr;
	testFile->Close();

        std::cout<<"-------------Saving resulting distributions option on-------------"<<std::endl;
        std::ofstream results;
        results.open( ( "Results/Lifetime_distribution_" + PathWithDate).c_str() );
        for( unsigned ite = 0; ite < SizeOfLFGrid; ite++ )
        {
                results << LifetimeGrid_root[ite] << "   " << Intensities_root[ite] <<std::endl;
        }
        results.close();

	delete testFile;
	
	TFile *testFile2 = new TFile( "Results_from_fit.root", "update" );
        testFile2 -> mkdir( ("Continous_fit_"+Path).c_str() );
        testFile2 -> cd( ("Continous_fit_"+Path).c_str() );
	
	vec ModelEstimation =  GradMax.t() * Inten;		//Later plus Background
	Double_t Values_root[ Values.size() ], Arguments_root[ Values.size() ], Model_root[ Values.size() ];
	for( unsigned l=0; l<Values.size(); l++ )
	{
		Values_root[l] = Values[l];
		Arguments_root[l] = Arguments[l];
		Model_root[l] = ModelEstimation(l) + Background;
	}
	TGraph *histo = new TGraph( Values.size(), Arguments_root, Values_root );
	TGraph *histoFit = new TGraph( Values.size(), Arguments_root, Model_root );
	
	TCanvas *c1 = new TCanvas( "c1", "", 710, 500 );
        c1 -> SetFillColor( 0 ); //Option for histogram
        c1 -> SetFrameBorderMode( 0 );
        c1 -> SetBorderSize( 2 );
        c1 -> SetFrameLineWidth( 2 );
	c1 -> SetLeftMargin(0.12);
	c1 -> SetRightMargin(0.08);
	c1 -> SetBorderMode(0);
	
	histo -> SetMarkerStyle( 21 );
        histo -> SetMarkerSize( 0.5 );
        histo -> SetTitle( "Histogram from data" );
	histo -> GetXaxis() -> SetTitle( "Time difference [ns]" );
        histo -> Draw( "AP" );
	
	histoFit -> SetMarkerStyle( 20 );
	histoFit -> SetMarkerColor( kRed );
        histoFit -> SetMarkerSize( 0.2 );
	histoFit -> SetLineColor( kRed );
	histoFit -> SetLineWidth( 2 );
        histoFit -> SetTitle( "Histogram from fit" );
        histoFit -> Draw( "same" );

	c1 -> SetLogy();
	histo -> Write( ("Raw_Histo_"+PathWithDate).c_str() );
	histoFit -> Write( ("Fit_Histo_"+PathWithDate).c_str() );
	
	c1 -> Write( PathWithDate.c_str() );
	delete c1;
	delete histo;
	delete histoFit;	
	delete testFile2;
	
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
	Double_t Residuals_root[ Values.size() ];
	
	for( unsigned i = 0; i < Values.size(); i++ )
	{
		Residuals_root[i] = (Values_root[i] - Model_root[i])/sqrt( Values_root[i] );
	}
	TFile *residuals = new TFile( "Residuals.root", "update" );
	residuals -> mkdir( ("Continous_fit_"+Path).c_str() ); 
	residuals -> cd( ("Continous_fit_"+Path).c_str() );  
        TGraph *residualsHisto = new TGraph( Values.size(), Arguments_root, Residuals_root );
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

void Fit::LinearFilter( unsigned FilterRange )
{
	double sum = 0.;
	std::vector<double> tempValues;
	for( unsigned j=0; j<FilterRange; j++ )
	{
		tempValues.push_back( Values[j] );
	}
	for( unsigned i=FilterRange; i<Values.size() - FilterRange; i++ )
	{	
		sum = 0.;
		for( unsigned k=0; k<=2*FilterRange + 1; k++ )			//From 0 to 2*FilterRange + 1, because filter should be symmetric around middle bin -> FilterRange points on the left + FilterRange points on the right + middle Point
		{
			sum += Values[i - FilterRange + k]/(2*FilterRange + 1);
		}
		tempValues.push_back(sum);
	}
	for( unsigned l=Values.size()-FilterRange; l<Values.size(); l++ )
	{
		tempValues.push_back(Values[l]);
	}
	Values = tempValues;
}

double Fit::FindingOptimalAreaParameter( std::vector< double > Sig )
{
	std::cout << "Area Parameter Finding" << std::endl;
	
	double sum = 0;
	vec Intensities( LifetimeGrid.size() );				//Vector for intensities calculated under given distribution with given sigmas
	for( unsigned i=0; i<LifetimeGrid.size(); i++ )			//Calculating intensities over whole LifetimeGrid
	{
		sum = 0;
		for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )	//Each Lifetime is separate Gaussian-like distribution
		{
			sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, Sig[j]);
		}
		Intensities(i) = sum;		//Intensities for calculations in this procedure
	}
	vec Model( Values.size() );		//Vector that characterizes estimation of histogram with sigmas given in the beginning of the function
	Model =  GradMax.t() * Intensities;	//GradMax is jus a matrix of Exponentially modified Gaussians for all bins in histogram and for all lifetimes from LifetimeGrid -> Multiplication by Intensities calculated following Gaussian-like distribution -> Model that is estimated to histogram
	vec Diff( Values.size() );		//For calculations of Gradient for sigmas
	double ChiSq = 0;			//ChiSquared value for sigmas given in the beginning
	double MeanDiff = 0.;
	for( unsigned i=0; i<Values.size(); i++ )
	{		
	  	Diff(i) =  Values[i] - Background - Model(i);		//Values -> Content of histogram;   Model + Background -> Fit to data
		ChiSq += pow( Diff(i), 2 ) / Values[i];
		MeanDiff = (Values[i] - Background)/(Model(i)*Values.size());
	}
	std::cout << "Mean Difference between histogram and model " << MeanDiff << std::endl;
	std::cout << "ChiSq: " << ChiSq << std::endl;
	
	AreaParameter = MeanDiff;
	unsigned NmbOfPointToMinimize = 30;
	std::vector<double> ChiSqVector, AreaParameters, ChiSqDiffVector;
	for( unsigned k=0; k<NmbOfPointToMinimize; k++ )
	{
		AreaParameters.push_back( AreaParameter - ((int)NmbOfPointToMinimize/2 - (int)k)*AreaParameter/50 );
		for( unsigned i=0; i<LifetimeGrid.size(); i++ )			//Calculating intensities over whole LifetimeGrid
		{
			sum = 0;
			for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )	//Each Lifetime is separate Gaussian-like distribution
			{
				sum += ( AreaParameter - ((int)NmbOfPointToMinimize/2 - (int)k)*AreaParameter/50 )*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, Sig[j]);
			}
			Intensities(i) = sum;		//Intensities for calculations in this procedure
		}
		Model =  GradMax.t() * Intensities;
		double ChiSq2 = 0;
		for( unsigned i=0; i<Values.size(); i++ )
		{		
			Diff(i) =  Values[i] - Background - Model(i);		//Values -> Content of histogram;   Model + Background -> Fit to data
			ChiSq2 += pow( Diff(i), 2 ) / Values[i];
		}
		ChiSqVector.push_back(ChiSq2);
		if( k>0 )
			ChiSqDiffVector.push_back( ChiSqVector[k] - ChiSqVector[k-1] );
	}
	unsigned MinElement = std::distance( std::begin(ChiSqVector), std::min_element( std::begin(ChiSqVector), std::end(ChiSqVector) ) );
	std::cout << "Area Parameter Minimal Element " << MinElement << std::endl;
	unsigned LeftSide = (MinElement > 5 ? MinElement - 5 : 0);
	unsigned RightSide = (MinElement < ChiSqVector.size() - 5 ? MinElement + 5 : ChiSqVector.size() - 1);
	double OptimalAreaParameter = FindPeak( AreaParameters, ChiSqDiffVector, LeftSide, RightSide ) + AreaParameter/100;
	std::cout << "Area Parameter Minimum From ChiSq " << OptimalAreaParameter << std::endl;
	return OptimalAreaParameter;
}

double Fit::PreIterateGauss( std::vector< double > Sig, double Step, double Lambda, double PreviousChi )
{
	std::cout << "PreIteration " << ++iterator << std::endl;
	if( iterator%NumberOfPreIterations == 0)			//Every ten iterations condition is checked
	{
		std::cout << "PreIteration ChiSq " << ChiDiffBetweenIterations << std::endl;
		if( ChiDiffBetweenIterations < PreviousChi*0.001 )		//If this condition is fulfilled last Sigmas are saved to container SSig, and PreIteration are ended
		{
			for( unsigned i=0; i<Sig.size(); i++ )
			{
				SSig[i] = Sig[i];
			}
			return 1;
		}
		ChiDiffBetweenIterations = 0.;
	}
	double sum = 0;
	std::vector< double > sums;		//For gradient calculation

	//Initializing vector elements to 0 -> Maybe in future will be introduced some funky method to do that
	for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )
	{
		sums.push_back(0);
	}
	mat Gradient_d( Sig.size(), LifetimeGrid.size() );		//Gradient over sigmas
	vec Intensities( LifetimeGrid.size() );				//Vector for intensities calculated under given distribution with given sigmas
	for( unsigned i=0; i<LifetimeGrid.size(); i++ )			//Calculating intensities over whole LifetimeGrid
	{
		sum = 0;
		for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )	//Each Lifetime is separate Gaussian-like distribution
		{
			sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, Sig[j]);
			//Value following Gaussian distribution
			sums[j] = AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, Sig[j]) * ( 2*pow( LifetimeGrid[i] - LifetimesFromDiscrete[j].Parameter ,2)/pow(Sig[j],3) - 1/Sig[j] );
			//Derivative of Gaussian distribution
			Gradient_d(j,i) = sums[j];
			//Gradient is vector of all derivatives
		}
		Intensities(i) = sum;		//Intensities for calculations in this procedure
		Intensities_afterIter[i] = sum;
	}
	vec Model( Values.size() );		//Vector that characterizes estimation of histogram with sigmas given in the beginning of the function
	Model =  GradMax.t() * Intensities;	//GradMax is jus a matrix of Exponentially modified Gaussians for all bins in histogram and for all lifetimes from LifetimeGrid -> Multiplication by Intensities calculated following Gaussian-like distribution -> Model that is estimated to histogram
	vec Diff( Values.size() );		//For calculations of Gradient for sigmas
	double ChiSq = 0;			//ChiSquared value for sigmas given in the beginning
	for( unsigned i=0; i<Values.size(); i++ )
	{		
	  	Diff(i) =  Values[i] - Background - Model(i);		//Values -> Content of histogram;   Model + Background -> Fit to data
		ChiSq += pow( Diff(i), 2 ) / Values[i];
	}
	std::cout << "ChiSq: " << ChiSq << std::endl;
	
	vec Gradient = Gradient_d*GradMax*Diff;			//Gradient definition in the Levenberg-Marquardt method 
	mat Hessian = Gradient_d *GradMax*GradMax.t() * Gradient_d.t();		//Hessian definition in the Levenberg-Marquardt method 
	mat HessToInv = Hessian + Lambda*diagmat(Hessian);			//Estimation of the real hessian
	vec Delta = solve(HessToInv, Gradient);					//Solving equation: Delta (Sigma) * Real Hessian = Gradient to obtain Delta (Sigma)
	vec IntensitiesDelta(Intensities.n_elem);				//Vector containing new Intensities shifted by values in Delta 
	std::vector< double > SigDelta;						//Vector for new sigmas
	for( unsigned j=0; j<Sig.size(); j++ )
	{
		SigDelta.push_back( (Sig[j] + Step*Delta[j] > 0.005 ? Sig[j] + Step*Delta[j] : 0.005) );		//Assuming that minila sigma is 10 ps - somehow reasonable value
		std::cout << LifetimesFromDiscrete[j].Parameter << " " << Sig[j] << " " << SigDelta[j] << std::endl;
	}
	for( unsigned i=0; i<LifetimeGrid.size(); i++ )		//Creating new values of intensities with new sigmas
	{
		sum = 0;
		for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )
		{
			sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, SigDelta[j]);
		}
		IntensitiesDelta(i) = sum;
	}
	vec Model2 =  GradMax.t() * IntensitiesDelta;		//New model with new sigmas
	double ChiSq2 = 0;					//ChiSquared value for sigmas corrected by Delta
	for( unsigned i=0; i<Values.size(); i++ )
		ChiSq2 += pow( Values[i] - Background - Model2(i), 2 ) / Values[i];
	std::cout << ChiSq2 << std::endl;
	if( ChiSq2 <= ChiSq )					//If new sigmas are better
	{
	  	ChiDiffBetweenIterations += ChiSq - ChiSq2;
		return PreIterateGauss( SigDelta, Step, 0.1*Lambda, ChiSq2 );	//Lambda lowered -> Gauss-Newton approach stronger -> New sigmas are promising
	}
	else
	{
		return PreIterateGauss( Sig, Step, 10*Lambda, ChiSq );		//Lambda increased -> Maximal damping method stronger -> New sigma is not good enough
	}
}

double Fit::PreIterateLogGauss( std::vector< double > Sig, double Step, double Lambda, double PreviousChi )
{
	std::cout << "PreIteration " << ++iterator << std::endl;
	if( iterator%NumberOfPreIterations == 0)			//Every ten iterations condition is checked
	{
		std::cout << "PreIteration ChiSq " << ChiDiffBetweenIterations << std::endl;
		if( ChiDiffBetweenIterations < PreviousChi*0.001 )		//If this condition is fulfilled last Sigmas are saved to container SSig, and PreIteration are ended
		{
			for( unsigned i=0; i<Sig.size(); i++ )
			{
				SSig[i] = Sig[i];
			}
			return 1;
		}
		ChiDiffBetweenIterations = 0.;
	}
	double sum = 0;
	std::vector< double > sums;		//For gradient calculation

	//Initializing vector elements to 0 -> Maybe in future will be introduced some funky method to do that
	for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )
	{
		sums.push_back(0);
	}
	mat Gradient_d( Sig.size(), LifetimeGrid.size() );		//Gradient over sigmas
	vec Intensities( LifetimeGrid.size() );				//Vector for intensities calculated under given distribution with given sigmas
	for( unsigned i=0; i<LifetimeGrid.size(); i++ )			//Calculating intensities over whole LifetimeGrid
	{
		sum = 0;
		for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )	//Each Lifetime is separate Gaussian-like distribution
		{
			sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(Sig[j],2)/2, Sig[j]);
			sums[j] = AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(Sig[j],2)/2, Sig[j] ) * ( 2* pow( log(LifetimeGrid[i]) - log( LifetimesFromDiscrete[j].Parameter ) - pow(Sig[j],2)/2
			,2)/pow(Sig[j],3) - 1/Sig[j] );
			Gradient_d(j,i) = sums[j];
		}
		Intensities(i) = sum;		//Intensities for calculations in this procedure
		Intensities_afterIter[i] = sum;
	}
	vec Model( Values.size() );		//Vector that characterizes estimation of histogram with sigmas given in the beginning of the function
	Model =  GradMax.t() * Intensities;	//GradMax is jus a matrix of Exponentially modified Gaussians for all bins in histogram and for all lifetimes from LifetimeGrid -> Multiplication by Intensities calculated following Gaussian-like distribution -> Model that is estimated to histogram
	vec Diff( Values.size() );		//For calculations of Gradient for sigmas
	double ChiSq = 0;			//ChiSquared value for sigmas given in the beginning
	for( unsigned i=0; i<Values.size(); i++ )
	{		
	  	Diff(i) =  Values[i] - Background - Model(i);		//Values -> Content of histogram;   Model + Background -> Fit to data
		ChiSq += pow( Diff(i), 2 ) / Values[i];
	}
	std::cout << "ChiSq: " << ChiSq << std::endl;
	
	vec Gradient = Gradient_d*GradMax*Diff;			//Gradient definition in the Levenberg-Marquardt method 
	mat Hessian = Gradient_d *GradMax*GradMax.t() * Gradient_d.t();		//Hessian definition in the Levenberg-Marquardt method 
	mat HessToInv = Hessian + Lambda*diagmat(Hessian);			//Estimation of the real hessian
	vec Delta = solve(HessToInv, Gradient);					//Solving equation: Delta (Sigma) * Real Hessian = Gradient to obtain Delta (Sigma)
	vec IntensitiesDelta(Intensities.n_elem);				//Vector containing new Intensities shifted by values in Delta 
	std::vector< double > SigDelta;						//Vector for new sigmas
	for( unsigned j=0; j<Sig.size(); j++ )
	{
		SigDelta.push_back( (Sig[j] + Step*Delta[j] > 0.005 ? Sig[j] + Step*Delta[j] : 0.005) );		//Assuming that minila sigma is 10 ps - somehow reasonable value
		std::cout << LifetimesFromDiscrete[j].Parameter << " " << Sig[j] << " " << SigDelta[j] << std::endl;
	}
	for( unsigned i=0; i<LifetimeGrid.size(); i++ )		//Creating new values of intensities with new sigmas
	{
		sum = 0;
		for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )
		{
			sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(SigDelta[j],2)/2, SigDelta[j] );
		}
		IntensitiesDelta(i) = sum;
	}
	vec Model2 =  GradMax.t() * IntensitiesDelta;		//New model with new sigmas
	double ChiSq2 = 0;					//ChiSquared value for sigmas corrected by Delta
	for( unsigned i=0; i<Values.size(); i++ )
		ChiSq2 += pow( Values[i] - Background - Model2(i), 2 ) / Values[i];
	std::cout << ChiSq2 << std::endl;
	if( ChiSq2 <= ChiSq )					//If new sigmas are better
	{
	  	ChiDiffBetweenIterations += ChiSq - ChiSq2;
		return PreIterateLogGauss( SigDelta, Step, 0.1*Lambda, ChiSq2 );	//Lambda lowered -> Gauss-Newton approach stronger -> New sigmas are promising
	}
	else
	{
		return PreIterateLogGauss( Sig, Step, 10*Lambda, ChiSq );		//Lambda increased -> Maximal damping method stronger -> New sigma is not good enough
	}
}

double Fit::PreIterateMixed( std::vector< double > Sig, double Step, double Lambda, double PreviousChi )
{
	std::cout << "PreIteration " << ++iterator << std::endl;
	if( iterator%NumberOfPreIterations == 0)			//Every ten iterations condition is checked
	{
		std::cout << "PreIteration ChiSq " << ChiDiffBetweenIterations << std::endl;
		if( ChiDiffBetweenIterations < PreviousChi*0.001 )		//If this condition is fulfilled last Sigmas are saved to container SSig, and PreIteration are ended
		{
			for( unsigned i=0; i<Sig.size(); i++ )
			{
				SSig[i] = Sig[i];
			}
			return 1;
		}
		ChiDiffBetweenIterations = 0.;
	}
	double sum = 0;
	std::vector< double > sums;		//For gradient calculation

	//Initializing vector elements to 0 -> Maybe in future will be introduced some funky method to do that
	for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )
	{
		sums.push_back(0);
	}
	mat Gradient_d( Sig.size(), LifetimeGrid.size() );		//Gradient over sigmas
	vec Intensities( LifetimeGrid.size() );				//Vector for intensities calculated under given distribution with given sigmas
	for( unsigned i=0; i<LifetimeGrid.size(); i++ )			//Calculating intensities over whole LifetimeGrid
	{
		sum = 0;
		for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )	//Each Lifetime is separate Gaussian-like distribution
		{
			if( LifetimesFromDiscrete[j].Parameter < 0.5 )
			{
			      sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, Sig[j]);
			      //Value following Gaussian distribution
			      sums[j] = AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, Sig[j]) * ( 2*pow( LifetimeGrid[i] - LifetimesFromDiscrete[j].Parameter ,2)/pow(Sig[j],3) - 1/Sig[j] );
			      //Derivative of Gaussian distribution
			      Gradient_d(j,i) = sums[j];
			      //Gradient is vector of all derivatives
			}
			else		//Option for assuming logGaussian distributions
			{
			      sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(Sig[j],2)/2, Sig[j]);
			      sums[j] = AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(Sig[j],2)/2, Sig[j] ) * ( 2* pow( log(LifetimeGrid[i]) - log( LifetimesFromDiscrete[j].Parameter ) - pow(Sig[j],2)/2
			      ,2)/pow(Sig[j],3) - 1/Sig[j] );
			      Gradient_d(j,i) = sums[j];
			}
		}
		Intensities(i) = sum;		//Intensities for calculations in this procedure
		Intensities_afterIter[i] = sum;
	}
	vec Model( Values.size() );		//Vector that characterizes estimation of histogram with sigmas given in the beginning of the function
	Model =  GradMax.t() * Intensities;	//GradMax is jus a matrix of Exponentially modified Gaussians for all bins in histogram and for all lifetimes from LifetimeGrid -> Multiplication by Intensities calculated following Gaussian-like distribution -> Model that is estimated to histogram
	vec Diff( Values.size() );		//For calculations of Gradient for sigmas
	double ChiSq = 0;			//ChiSquared value for sigmas given in the beginning
	for( unsigned i=0; i<Values.size(); i++ )
	{		
	  	Diff(i) =  Values[i] - Background - Model(i);		//Values -> Content of histogram;   Model + Background -> Fit to data
		ChiSq += pow( Diff(i), 2 ) / Values[i];
	}
	std::cout << "ChiSq: " << ChiSq << std::endl;
	
	vec Gradient = Gradient_d*GradMax*Diff;			//Gradient definition in the Levenberg-Marquardt method 
	mat Hessian = Gradient_d *GradMax*GradMax.t() * Gradient_d.t();		//Hessian definition in the Levenberg-Marquardt method 
	mat HessToInv = Hessian + Lambda*diagmat(Hessian);			//Estimation of the real hessian
	vec Delta = solve(HessToInv, Gradient);					//Solving equation: Delta (Sigma) * Real Hessian = Gradient to obtain Delta (Sigma)
	vec IntensitiesDelta(Intensities.n_elem);				//Vector containing new Intensities shifted by values in Delta 
	std::vector< double > SigDelta;						//Vector for new sigmas
	for( unsigned j=0; j<Sig.size(); j++ )
	{
		SigDelta.push_back( (Sig[j] + Step*Delta[j] > 0.005 ? Sig[j] + Step*Delta[j] : 0.005) );		//Assuming that minila sigma is 10 ps - somehow reasonable value
		std::cout << LifetimesFromDiscrete[j].Parameter << " " << Sig[j] << " " << SigDelta[j] << std::endl;
	}
	for( unsigned i=0; i<LifetimeGrid.size(); i++ )		//Creating new values of intensities with new sigmas
	{
		sum = 0;
		for( unsigned j=0; j<LifetimesFromDiscrete.size(); j++ )
		{
			if( LifetimesFromDiscrete[j].Parameter < 0.5 )
				sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*GaussDistr( LifetimeGrid[i], LifetimesFromDiscrete[j].Parameter, SigDelta[j]);
			else
				sum += AreaParameter*IntensitiesFromDiscrete[j].Parameter*LogGaussDistr( LifetimeGrid[i], log( LifetimesFromDiscrete[j].Parameter ) - pow(SigDelta[j],2)/2, SigDelta[j] );
		}
		IntensitiesDelta(i) = sum;
	}
	vec Model2 =  GradMax.t() * IntensitiesDelta;		//New model with new sigmas
	double ChiSq2 = 0;					//ChiSquared value for sigmas corrected by Delta
	for( unsigned i=0; i<Values.size(); i++ )
		ChiSq2 += pow( Values[i] - Background - Model2(i), 2 ) / Values[i];
	std::cout << ChiSq2 << std::endl;
	if( ChiSq2 <= ChiSq )					//If new sigmas are better
	{
	  	ChiDiffBetweenIterations += ChiSq - ChiSq2;
		return PreIterateMixed( SigDelta, Step, 0.1*Lambda, ChiSq2 );	//Lambda lowered -> Gauss-Newton approach stronger -> New sigmas are promising
	}
	else
	{
		return PreIterateMixed( Sig, Step, 10*Lambda, ChiSq );		//Lambda increased -> Maximal damping method stronger -> New sigma is not good enough
	}
}

double Fit::IterateExp( vec Intensities, double Step, double Lambda, double PreviousChi )		//Similar procedure as in PreIteration but fitting intesities not sigmas
{
	std::cout << "Iteration: " << ++iterator << std::endl;
	vec Model =  GradMax.t() * Intensities;
	vec Diff( Values.size() );
	double ChiSq = 0;
	double ModelDiff = 0.;
	for( unsigned i=0; i<Values.size(); i++ )
	{
	  	Diff(i) =  Values[i] - Background - Model(i);
		ModelDiff += ( Values[i] - Model(i) )/Values.size();
	}
	std::cout << ModelDiff << " " << Background << " " << SDBackground << std::endl;
	double BGCorr = 0.;
	if( ModelDiff < Background - 2*SDBackground )
		BGCorr = Background + 2*SDBackground;
	else if( ModelDiff > Background + 2*SDBackground )
		BGCorr = Background - 2*SDBackground;
	else
		BGCorr = 2*Background - ModelDiff;
	for( unsigned i=0; i<Values.size(); i++ )
	{
	  	Diff(i) =  Values[i] - BGCorr - Model(i);
		ChiSq += pow( Diff(i), 2 ) / Values[i];
	}
	std::cout << "ChiSq: " << ChiSq << std::endl;
	if( iterator%NumberOfIterations == 0 )
	{
		std::cout << "ChiSq Difference: " << ChiDiffBetweenIterations << std::endl;
		if( ChiDiffBetweenIterations < PreviousChi*0.001 )
		{
			for( unsigned p=0; p<Intensities.n_elem; p++)
			{
				Intensities_afterIter.push_back(Intensities(p));
			}
			if( ChiSq > 2000 )
			{
				double Lama;
				std::cout << "You can give new Lambda value if You want" << std::endl;
				std::cin >> Lama;
				for( unsigned i=0; i<LifetimeGrid.size(); i++ )	
					Intensities_afterIter[i] = Intensities(i);
				return -Lama;
			}
			else
			{
				for( unsigned i=0; i<LifetimeGrid.size(); i++ )	
					Intensities_afterIter[i] = Intensities(i);
				return 1.;	    		
			}
		}
		else
		{
			for( unsigned i=0; i<LifetimeGrid.size(); i++ )	
				Intensities_afterIter[i] = Intensities(i);
			return -0.01;
		}
	}
	vec Gradient = GradMax*Diff;
	mat Hessian = GradMax * GradMax.t();
	mat HessToInv = Hessian + Lambda*diagmat(Hessian);
	vec Delta = solve(HessToInv, Gradient);
	vec IntensitiesDelta(Intensities.n_elem);
	for( unsigned j=0; j<Intensities.n_elem; j++ )
	{
		if( Intensities(j) + Step*Delta(j) < 0)
			IntensitiesDelta(j) = 0;
		else	
			IntensitiesDelta(j) = Intensities(j) + Step*Delta(j);
	}
	vec Model2 =  GradMax.t() * IntensitiesDelta;
	double ChiSq2 = 0;
	ModelDiff = 0.;
	for( unsigned i=0; i<Values.size(); i++ )
	{
	  	Diff(i) =  Values[i] - Background - Model2(i);
		ModelDiff += ( Values[i] - Model2(i) )/Values.size();
	}
	std::cout << ModelDiff << " " << Background << " " << SDBackground << std::endl;
	BGCorr = 0.;
	if( ModelDiff < Background - 2*SDBackground )
		BGCorr = Background + 2*SDBackground;
	else if( ModelDiff > Background + 2*SDBackground )
		BGCorr = Background - 2*SDBackground;
	else
		BGCorr = 2*Background - ModelDiff;
	for( unsigned i=0; i<Values.size(); i++ )
		ChiSq2 += pow( Values[i] - BGCorr - Model2(i), 2 ) / Values[i];
	if( ChiSq2 <= ChiSq )
	{
	  	ChiDiffBetweenIterations += ChiSq - ChiSq2;
		return IterateExp( IntensitiesDelta, Step, 0.1*Lambda, ChiSq2 );
	}
	else
	{
		return IterateExp( Intensities, Step, 10*Lambda, ChiSq );
	}
}

double Fit::IterateExp2( vec Intensities, double Step, double Lambda, double PreviousChi )		//Similar procedure as in PreIteration but fitting intesities not sigmas
{
	std::cout << "Iteration: " << ++iterator << std::endl;
	vec Model =  GradMax.t() * Intensities;
	vec Diff( Values.size() );
	double ChiSq = 0;
	double ModelDiff = 0.;
	for( unsigned i=0; i<Values.size(); i++ )
	{
	  	Diff(i) =  Values[i] - Background - Model(i);
		ModelDiff += ( Values[i] - Model(i) )/Values.size();
	}
	std::cout << ModelDiff << " " << Background << " " << SDBackground << std::endl;
	double BGCorr = 0.;
	if( ModelDiff < Background - SDBackground )
		BGCorr = Background + SDBackground;
	else if( ModelDiff > Background + SDBackground )
		BGCorr = Background - SDBackground;
	else
		//BGCorr = 2*Background - ModelDiff;
		BGCorr = Background;
	for( unsigned i=0; i<Values.size(); i++ )
	{
	  	Diff(i) =  Values[i] - BGCorr - Model(i);
		ChiSq += pow( Diff(i), 2 ) / Values[i];
	}
	std::cout << "ChiSq: " << ChiSq << std::endl;
	if( iterator%NumberOfIterations == 0 )
	{
		std::cout << "ChiSq Difference: " << ChiDiffBetweenIterations << std::endl;
		if( ChiDiffBetweenIterations < PreviousChi*0.001 )
		{
			for( unsigned p=0; p<Intensities.n_elem; p++)
			{
				Intensities_afterIter.push_back(Intensities(p));
			}
			if( ChiSq > 2000 )
			{
				double Lama;
				std::cout << "You can give new Lambda value if You want" << std::endl;
				std::cin >> Lama;
				for( unsigned i=0; i<LifetimeGrid.size(); i++ )	
					Intensities_afterIter[i] = Intensities(i);
				return -Lama;
			}
			else
			{
				for( unsigned i=0; i<LifetimeGrid.size(); i++ )	
					Intensities_afterIter[i] = Intensities(i);
				return 1.;	    		
			}
		}
		else
		{
			for( unsigned i=0; i<LifetimeGrid.size(); i++ )	
				Intensities_afterIter[i] = Intensities(i);
			return -0.01;
		}
	}
	vec Gradient = GradMax*Diff;
	mat Hessian = GradMax * GradMax.t();
	mat HessToInv = Hessian + Lambda*diagmat(Hessian);
	vec Delta = solve(HessToInv, Gradient);
	vec IntensitiesDelta(Intensities.n_elem);
	for( unsigned j=0; j<Intensities.n_elem; j++ )
	{
		if( Intensities(j) + Step*Delta(j) < 0)
			IntensitiesDelta(j) = 0;
		else	
			IntensitiesDelta(j) = Intensities(j) + Step*Delta(j);
	}
	vec Model2 =  GradMax.t() * IntensitiesDelta;
	double ChiSq2 = 0;
	ModelDiff = 0.;
	for( unsigned i=0; i<Values.size(); i++ )
	{
	  	Diff(i) =  Values[i] - Background - Model2(i);
		ModelDiff += ( Values[i] - Model2(i) )/Values.size();
	}
	std::cout << ModelDiff << " " << Background << " " << SDBackground << std::endl;
	BGCorr = 0.;
	if( ModelDiff < Background - SDBackground )
		BGCorr = Background + SDBackground;
	else if( ModelDiff > Background + SDBackground )
		BGCorr = Background - SDBackground;
	else
		//BGCorr = 2*Background - ModelDiff;
		BGCorr = Background;
	for( unsigned i=0; i<Values.size(); i++ )
		ChiSq2 += pow( Values[i] - BGCorr - Model2(i), 2 ) / Values[i];
	if( ChiSq2 <= ChiSq )
	{
	  	ChiDiffBetweenIterations += ChiSq - ChiSq2;
		return IterateExp( IntensitiesDelta, Step, 0.1*Lambda, ChiSq2 );
	}
	else
	{
		return IterateExp( Intensities, Step, 10*Lambda, ChiSq );
	}
}

double Fit::Iterate( vec Intensities, double Step, double Lambda, double PreviousChi )		//Similar procedure as in PreIteration but fitting intesities not sigmas
{
	std::cout << "Iteration: " << ++iterator << std::endl;
	vec Model =  GradMax.t() * Intensities;
	vec Diff( Values.size() );
	double ChiSq = 0;
	for( unsigned i=0; i<Values.size(); i++ )
	{
	  	Diff(i) =  Values[i] - Background - Model(i);
		ChiSq += pow( Diff(i), 2 ) / Values[i];
	}
	std::cout << "ChiSq: " << ChiSq << std::endl;
	if( iterator%NumberOfIterations == 0 )
	{
		std::cout << "ChiSq Difference: " << ChiDiffBetweenIterations << std::endl;
		if( ChiDiffBetweenIterations < PreviousChi*0.001 )
		{
			for( unsigned p=0; p<Intensities.n_elem; p++)
			{
				Intensities_afterIter.push_back(Intensities(p));
			}
			if( ChiSq > 2000 )
			{
				double Lama;
				std::cout << "You can give new Lambda value if You want" << std::endl;
				std::cin >> Lama;
				for( unsigned i=0; i<LifetimeGrid.size(); i++ )	
					Intensities_afterIter[i] = Intensities(i);
				return -Lama;
			}
			else
			{
				for( unsigned i=0; i<LifetimeGrid.size(); i++ )	
					Intensities_afterIter[i] = Intensities(i);
				return 1.;	    		
			}
		}
		else
		{
			for( unsigned i=0; i<LifetimeGrid.size(); i++ )	
				Intensities_afterIter[i] = Intensities(i);
			return -0.01;
		}
	}
	vec Gradient = GradMax*Diff;
	mat Hessian = GradMax * GradMax.t();
	mat HessToInv = Hessian + Lambda*diagmat(Hessian);
	vec Delta = solve(HessToInv, Gradient);
	vec IntensitiesDelta(Intensities.n_elem);
	for( unsigned j=0; j<Intensities.n_elem; j++ )
	{
		if( Intensities(j) + Step*Delta(j) < 0)
			IntensitiesDelta(j) = 0;
		else	
			IntensitiesDelta(j) = Intensities(j) + Step*Delta(j);
	}
	vec Model2 =  GradMax.t() * IntensitiesDelta;
	double ChiSq2 = 0;
	for( unsigned i=0; i<Values.size(); i++ )
		ChiSq2 += pow( Values[i] - Background - Model2(i), 2 ) / Values[i];
	if( ChiSq2 <= ChiSq )
	{
	  	ChiDiffBetweenIterations += ChiSq - ChiSq2;
		return Iterate( IntensitiesDelta, Step, 0.1*Lambda, ChiSq2 );
	}
	else
	{
		return Iterate( Intensities, Step, 10*Lambda, ChiSq );
	}
}

void TestIntensities()
{
	std::vector<LifetimeComponent> Parameters;
	Parameters.push_back( LifetimeComponent( 0.1, 0.1, "gauss" ) );
	Parameters.push_back( LifetimeComponent( 0.1, 0.25, "gauss" ) );
	Parameters.push_back( LifetimeComponent( 0.1, 0.4, "gauss" ) );
	Parameters.push_back( LifetimeComponent( 0.1, 0.1, "gauss" ) );
	Parameters.push_back( LifetimeComponent( 0.1, 0.15, "gauss" ) );
	std::vector<Double_t> Parameters2;
	Double_t Hmm1Parameter = SetIntensityParameter( Parameters, 1 );
	Double_t Hmm2Parameter = SetIntensityParameter( Parameters, 2 );
	Double_t Hmm3Parameter = SetIntensityParameter( Parameters, 3 );
	Double_t Hmm4Parameter = SetIntensityParameter( Parameters, 4 );
	Double_t Hmm5Parameter = SetIntensityParameter( Parameters, 5 );
	Parameters2.push_back(Hmm1Parameter);
	Parameters2.push_back(Hmm2Parameter);
	Parameters2.push_back(Hmm3Parameter);
	Parameters2.push_back(Hmm4Parameter);
	Parameters2.push_back(Hmm5Parameter);
	std::cout << Hmm1Parameter << " " << Hmm2Parameter << " " << Hmm3Parameter << " " << Hmm4Parameter << " " << Hmm5Parameter << std::endl;
	std::cout << GetIntensityParameter( Parameters2, 1 ) << " " << GetIntensityParameter( Parameters2, 2 ) << " " << GetIntensityParameter( Parameters2, 3 ) << " " << GetIntensityParameter( Parameters2, 4 ) << " " << GetIntensityParameter( Parameters2, 5 )<< std::endl;
	std::cout << GetIntensityParameter( Parameters2, 1 ) + GetIntensityParameter( Parameters2, 2 ) + GetIntensityParameter( Parameters2, 3 ) + GetIntensityParameter( Parameters2, 4 ) << std::endl;
	Parameters2[1] -= 0.5;
	std::cout << GetIntensityParameter( Parameters2, 1 ) << " " << GetIntensityParameter( Parameters2, 2 ) << " " << GetIntensityParameter( Parameters2, 3 ) << " " << GetIntensityParameter( Parameters2, 4 ) << " " << GetIntensityParameter( Parameters2, 5 )<< std::endl;
	std::cout << GetIntensityParameter( Parameters2, 1 ) + GetIntensityParameter( Parameters2, 2 ) + GetIntensityParameter( Parameters2, 3 ) + GetIntensityParameter( Parameters2, 4 ) << std::endl;
	Parameters2[0] -= 0.1;
	std::cout << GetIntensityParameter( Parameters2, 1 ) << " " << GetIntensityParameter( Parameters2, 2 ) << " " << GetIntensityParameter( Parameters2, 3 ) << " " << GetIntensityParameter( Parameters2, 4 ) << " " << GetIntensityParameter( Parameters2, 5 )<< std::endl;
	std::cout << GetIntensityParameter( Parameters2, 1 ) + GetIntensityParameter( Parameters2, 2 ) + GetIntensityParameter( Parameters2, 3 ) + GetIntensityParameter( Parameters2, 4 ) << std::endl;
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

Double_t DiscreteFitFunction( Double_t *A, Double_t *P )	//First parameter - nmbr of comp, Second - nmb of Gauss
{
	Double_t sum = 0.;
	Double_t FixedFraction = 0.;
	Double_t pPsFixedFraction = 0.;
	Double_t pPsFreeFraction = 0.;
	Double_t pPsFraction = 0.;
	Double_t FreeFraction = 0.;
	std::vector<Double_t> ResolutionIntensities;
	std::vector<unsigned> FreepPsIndexes;
	unsigned FreeIterator = 0, pPsIndex = -1;	
	for( unsigned i = 0; i < P[1]; i++ )
	{
		ResolutionIntensities.push_back( P[5 + 3*i] );
	}
	for( unsigned j = 0; j < P[2]; j++ )
	{
		if( P[4 + 3*(int)P[1] + 2*j] < 0 && P[5 + 3*(int)P[1] + 2*j] > 0 )
		{
			FixedFraction += P[5 + 3*(int)P[1] + 2*j];
			if( P[4 + 3*(int)P[1] + 2*j] < -0.7 )
				pPsFixedFraction += P[5 + 3*(int)P[1] + 2*j];
		}
		else if( P[5 + 3*(int)P[1] + 2*j] > 0 )
		{
			if( P[4 + 3*(int)P[1] + 2*j] > 0.7 )
				FreepPsIndexes.push_back(j);
			FreeFraction += P[5 + 3*(int)P[1] + 2*j];
		}
		else
		{
			pPsFraction -= P[5 + 3*(int)P[1] + 2*j];
			j = pPsIndex;
		}
	}
	if( pPsIndex + 1 )
	{
		for( unsigned j = 0; j < FreepPsIndexes.size(); j++ )
		{
			pPsFreeFraction += pPsFraction*P[5 + 3*(int)P[1] + 2*j];
		}
	}
	Double_t FreeIntensity = 0.;
	for( unsigned i = 0; i < P[1]; i++ )
	{
		//std::cin >> FreeIntensity;
		for( unsigned j = 0; j < P[2]; j++ )
		{
			if( P[4 + 3*(int)P[1] + 2*j] < 0. )
			{
				sum += P[3]/* Area */ * GetIntensityParameter( ResolutionIntensities, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, -P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
			}
			else if( P[5 + 3*(int)P[1] + 2*j] >= 0. )
			{
				FreeIterator++;
				FreeIntensity = (1-FixedFraction-pPsFixedFraction*pPsFraction)*(1-pPsFreeFraction)*P[5 + 3*(int)P[1] + 2*j]/FreeFraction;
				sum += P[3]/* Area */ * GetIntensityParameter( ResolutionIntensities, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, FreeIntensity/* Intensity */ );
			}
			
		}		
	}
	if( pPsIndex + 1 )
	{
		for( unsigned i = 0; i < P[1]; i++ )
		{
			sum += P[3]/* Area */ * GetIntensityParameter( ResolutionIntensities, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*pPsIndex]/* Lifetime */, pPsFixedFraction+(1-FixedFraction-pPsFixedFraction*pPsFraction)*pPsFreeFraction/* Intensity */ );
		}
	}
	return sum + P[0];
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
		for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]; j < P[2]; j++ )
		{
			sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, GetIntensityParameterNew( P, 0, 5 + 3*(int)P[1] + 2*P[6 + 3*(int)P[1] + 2*((int)P[2]-1)], (j-P[6 + 3*(int)P[1] + 2*((int)P[2]-1)])+1 )*(1 - FixedIntensity )/* Intensity */ );
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
			sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/FreeIntensity/* Intensity */ );
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
		for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]; j < P[2]; j++ )
		{
			if( i == 1 )
			{
				sum += P[3]/* Area */ * P[5]*sin(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/FreeIntensity/* Intensity */ );
			}
			else if( i == 2 )
			{
				sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/FreeIntensity/* Intensity */ );
			}
			else
			{
				sum += P[3]/* Area */ * P[5]*cos(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/FreeIntensity/* Intensity */ );			  
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
				sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *(ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
				FixedIntensity += P[5 + 3*(int)P[1] + 2*j]*(1 + P[5 + 3*(int)P[1]]);
			}
			else
			{
				sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
				FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
			}
		}
		for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j < P[2]; j++ )
		{
			sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, GetIntensityParameterNew( P, 2, 5 + 3*(int)P[1] + 2*P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 2, (j-P[6 + 3*(int)P[1] + 2*((int)P[2]-1)]) )*(1 - FixedIntensity )/(1 + FreepPsIntensity)/* Intensity */ );
		}
		sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, (1 - FixedIntensity )*FreepPsIntensity/(1 + FreepPsIntensity)/* Intensity */ );
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
				sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *(ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
				FixedIntensity += P[5 + 3*(int)P[1] + 2*j]*(1 + P[5 + 3*(int)P[1]]);
			}
			else
			{
				sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
				FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
			}
		}
		for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j < P[2]; j++ )
		{
			sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
		}
		sum += P[3]/* Area */ * GetIntensityParameterNew( P, 1, 5, i+1 )/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, (1 - FixedIntensity )*FreepPsIntensity/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
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
			sum += P[3]/* Area */ * 1/* Intensity */ *(ExMoGa( A[0], P[4]/* Sigma */, P[6]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + ExMoGa( A[0], P[4]/* Sigma */, P[6]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
			FixedIntensity += P[5 + 3*(int)P[1] + 2*j]*(1 + P[5 + 3*(int)P[1]]);
		}
		else
		{
			sum += P[3]/* Area */ * 1/* Intensity */ *ExMoGa( A[0], P[4]/* Sigma */, P[6]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ );
			FixedIntensity += P[5 + 3*(int)P[1] + 2*j];
		}
	}
	for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j < P[2]; j++ )
	{
		sum += P[3]/* Area */ * 1/* Intensity */ *ExMoGa( A[0], P[4]/* Sigma */, P[6]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
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
					sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *(ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
					FixedIntensity += P[5 + 3*(int)P[1] + 2*j]*(1 + P[5 + 3*(int)P[1]]);				  
				}
				else
				{
					sum += P[3]/* Area */ * P[5]/* Intensity */ *(ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
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
		for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j < P[2]; j++ )
		{
			if( i )
			{
				sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );	  
			}
			else
			{
				sum += P[3]/* Area */ * P[5]/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
			}
		}
		if( i )
		{
			sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, (1 - FixedIntensity )*FreepPsIntensity/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
			FixedIntensity = 0.;		  
		}
		else
		{
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
					sum += P[3]/* Area */ * P[5]*sin(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *(ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
					FixedIntensity += P[5 + 3*(int)P[1] + 2*j]*(1 + P[5 + 3*(int)P[1]]);				  
				}
				else if( i == 2 )
				{
					sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *(ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
					FixedIntensity += P[5 + 3*(int)P[1] + 2*j]*(1 + P[5 + 3*(int)P[1]]);
				}
				else
				{
					sum += P[3]/* Area */ * P[5]*cos(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *(ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]/* Intensity */ ) + ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*P[5 + 3*(int)P[1]]/* Intensity */ ));
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
		for( unsigned j = P[6 + 3*(int)P[1] + 2*((int)P[2]-1)] + 1; j < P[2]; j++ )
		{
			if( i == 1 )
			{
				sum += P[3]/* Area */ * P[5]*sin(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );	  
			}
			else if( i == 2 )
			{
				sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
			}
			else
			{
				sum += P[3]/* Area */ * P[5]*cos(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1] + 2*j]/* Lifetime */, P[5 + 3*(int)P[1] + 2*j]*(1 - FixedIntensity )/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
			}
		}
		if( i == 1 )
		{
			sum += P[3]/* Area */ * P[5]*sin(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, (1 - FixedIntensity )*FreepPsIntensity/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
			FixedIntensity = 0.;		  
		}
		else if( i == 2 )
		{
			sum += P[3]/* Area */ * (1-P[5])/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, (1 - FixedIntensity )*FreepPsIntensity/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
			FixedIntensity = 0.;
		}
		else
		{
			sum += P[3]/* Area */ * P[5]*cos(P[8])/(cos(P[8]) + sin(P[8]))/* Intensity */ *ExMoGa( A[0], P[4 + 3*i]/* Sigma */, P[6 + 3*i]/* Offset */, P[4 + 3*(int)P[1]]/* Lifetime */, (1 - FixedIntensity )*FreepPsIntensity/(1 + FreepPsIntensity)/FreeIntensity/* Intensity */ );
			FixedIntensity = 0.;		  
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

double Fit::FindMaxLF()
{
	double max = 0.1;
	for( unsigned i=0; i<LifetimesFromDiscrete.size(); i++ )
	{
		if( LifetimesFromDiscrete[i].Parameter > max )
			max = LifetimesFromDiscrete[i].Parameter;
	}
	return max;
}

double Fit::FindMinLF()
{
 	double min = 10000.;
	for( unsigned i=0; i<LifetimesFromDiscrete.size(); i++ )
	{
		if( LifetimesFromDiscrete[i].Parameter < min )
			min = LifetimesFromDiscrete[i].Parameter;
	}
	return min; 
}

Double_t Filter( std::vector<double> Vec, unsigned Pos, unsigned Range )
{
	Double_t Res = 0;
	if( Pos < Range )
		return Res;
	for( unsigned i=Pos-Range; i<Pos+Range+1; i++ )
	{
		Res += Vec[i];
	}
	Res = Res/(2*Range+1);
	if( Range == 0 )
		Res = Vec[Pos];
	return Res;
}

double GaussDistr( double X, double Mean, double Sigma )
{
	return 1/sqrt(2*M_PI)/Sigma*exp(-pow(X-Mean, 2)/2/pow(Sigma, 2));
}

double LogGaussDistr( double X, double Mean, double Sigma )
{
	return 1/sqrt(2*M_PI)/Sigma/X*exp(-pow( log(X) - Mean, 2)/2/pow(Sigma, 2));
}

TH1F* FillHistogram( std::string Name, int BinNumber, double Start, double Stop, std::vector<double>& Data )
{
        TH1F *histogram = new TH1F( Name.c_str(), Name.c_str(), BinNumber, Start, Stop );
        for( unsigned i = 0; i < Data.size(); i++ )
                histogram -> Fill( Data[i], 1 );
        return histogram;
}

std::string NumberToChar( double Number, int Precision)		//Converting Numbers to char
{
        std::ostringstream conv;
        conv << std::fixed << std::setprecision( Precision );
        conv << Number;
        return conv.str();
}

bool PathCheck( std::string Path )
{
	struct stat buffer;
	return stat ( Path.c_str(), &buffer ) == 0 && S_ISDIR(buffer.st_mode);
}

bool FileCheck(const std::string& NameOfFile) 
{
	struct stat buffer;   
	return (stat (NameOfFile.c_str(), &buffer) == 0); 
}

bool CheckTheLine( char FirstCharacter, char TestCharacter, unsigned NumberOfLine )
{
	if( FirstCharacter != FirstCharacter )
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << NumberOfLine << " line is not starting with character = " << TestCharacter << "..."<< std::endl;
		return 0;
	} 
	return 1;
}

void Fit::SortLifetimesByType()
{
	std::vector<unsigned> indexesFixed;
	int indexpPs = -1;
	for( unsigned i=0; i<Lifetimes.size(); i++ )
	{
		if( Lifetimes[i].Type != "nf" && Lifetimes[i].Type != "ps" )
			indexesFixed.push_back(i);
		else if( Lifetimes[i].Type == "ps" )
			indexpPs = i;
	}
	std::vector<LifetimeComponent> NewLifetimes;
	if( indexpPs >= 0 )
		NewLifetimes.push_back( Lifetimes[indexpPs] );
	if( indexesFixed.size() > 0 )
	{
		for( unsigned j=0; j<indexesFixed.size(); j++)
		{
			NewLifetimes.push_back( Lifetimes[indexesFixed[j]] );
		}
	}
	for( unsigned i=0; i<Lifetimes.size(); i++ )
	{
		if( Lifetimes[i].Type == "nf"  )
			NewLifetimes.push_back( Lifetimes[i] );
	}
	Lifetimes = NewLifetimes;
}

const std::string GetTime()
{
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	strftime( buf, sizeof(buf), "%Y-%m-%d_Time:%X", &tstruct );
	return buf;
}

double FindPeak( std::vector< double > Argument, std::vector< double > Values, unsigned Lowerfit, unsigned Higherfit )
{
	double XX=0, XY=0, YY=0, X=0, Y=0;
	unsigned size = Higherfit - Lowerfit;
        for( unsigned k = Lowerfit; k < Higherfit; k++ )
        {
                XX += Argument[k]*Argument[k];
                XY += Argument[k]*Values[k];
                YY += Values[k]*Values[k];
                X += Argument[k];
                Y += Values[k];
        }
        double a = (size*XY-X*Y) / (size*XX - X*X);
        double b = (Y-a*X)/size;
	double peak;
	if(a)
	    peak = -b / a;
	else 
	    peak = 0;
	return peak;
}
