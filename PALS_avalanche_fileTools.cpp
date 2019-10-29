#include "PALS_avalanche_fileTools.h"

FileTools::FileTools()
{
	FitDetailsFile = "FitDetails";
}

FileTools::FileTools( std::string fitDetailsFile )
{
	FitDetailsFile = fitDetailsFile;
}

std::string FileTools::NumberToChar( double Number, int Precision)
{
	std::ostringstream conv;
	conv << std::fixed << std::setprecision( Precision );
	conv << Number;
	return conv.str();    
}

bool FileTools::PathCheck( std::string Path )
{
	struct stat buffer;
	return stat ( Path.c_str(), &buffer ) == 0 && S_ISDIR(buffer.st_mode);
}

bool FileTools::FileCheck(const std::string& NameOfFile) 
{
	struct stat buffer;   
	return (stat (NameOfFile.c_str(), &buffer) == 0); 
}

bool FileTools::CheckTheLine( char FirstCharacter, char TestCharacter, unsigned NumberOfLine )
{
	if( FirstCharacter != FirstCharacter )
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << NumberOfLine << " line is not starting with character = " << TestCharacter << "..."<< std::endl;
		return 0;
	} 
	return 1;
}

std::vector<int> FileTools::GetComponentsMultiplicities( std::string SizesToExtract )
{
	std::stringstream StringToSplit( SizesToExtract );
	std::string tempExtraction;
	std::vector<int> Sizes;
	while (std::getline( StringToSplit, tempExtraction, ' ' )) 
	{
		Sizes.push_back( stoi(tempExtraction) );
	}
	return Sizes;
}

std::vector<std::string> FileTools::GetComponentDetails( std::string DetailsToExtract )
{
	std::stringstream StringToSplit( DetailsToExtract );
	std::string tempExtraction;
	std::vector<std::string> Details;
	while(std::getline( StringToSplit, tempExtraction, ' ' ) ) 
	{
		Details.push_back( tempExtraction );
	}
	return Details;
}

std::vector<std::string> FileTools::GetFitDetails()
{
	std::vector<std::string> FitDetailsData;
	FitDetailsData.push_back("");       //First element as a test of correctness of the FitDetails file
	FitDetailsData.push_back("");       //Second element as a indicator how muh there is Resolution and Lifetime components
	int lineWidth = 120;                 //120 bo 120 to dwie kopy -> kopać można nogami > nogi są dwie || It is carefully taken value, chosen by optimizing one complicated chosing procedure
	char line[120];					      	//For lines with description and header
	char question;						//For user handling with overexaggeration with parameters
	int lineNumber = 0;                //Indicator which line is processed
	std::string option = "", option2 = "", option3 = "";  	//Sometimes there are three arguments in line to get -> (Lifetime Intensity Type)
	std::cout << "-------------Processing fit details-------------" << std::endl;
	std::ifstream fitDetails;				//Stream for reading file
	fitDetails.open( FitDetailsFile );			// Opening the file with details of the analysis -> "FitDetails" - name of the file (Check)
//-----------------------------------------------
	std::cout << "-------------Reading the FitDetails-------------" << std::endl;
	for( unsigned i=0; i<3; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'T', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );  
	fitDetails >> option;		//Getting the type of data

	if( option != "histogram" && option != "times" && option != "ROOT" )	//Checking if type of data is provided correctly, one of the recognizable type - histogram or times
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Type of data is not one of correct types - 'histogram', 'times' or 'ROOT' " << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Type of Data ";
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'N', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	for( unsigned i=0; i<lineWidth; i++ )
	{
		line[i] = '\n';
	}
	char lineTemp[lineWidth];
	for( unsigned i=0; i<lineWidth; i++ )
	{
		lineTemp[i] = '\n';
	}
	fitDetails.getline( lineTemp, lineWidth );
	std::string optionTemp = "";
	int test000 = 1;
	while( test000 )
	{
		optionTemp = optionTemp + lineTemp[test000-1];
		//std::cout << line[test000-1] << std::endl;
		test000++;
		if( lineTemp[test000-1] == '\n' )
			test000 = 0;
	}
	optionTemp.erase( optionTemp.end()-1, optionTemp.end() );
	std::cout << "Histogram name provided: " << optionTemp << std::endl;
	FitDetailsData.push_back( optionTemp );
 	for( unsigned i=0; i<1; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'N', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	FitDetailsData.push_back( option );
 	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	
	if( CheckTheLine( line[0], 'N', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<3; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'W', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;		//Getting the Width of the Bin
	if( stod( option ) <= 0 || stod( option ) >= 142 )	//Checking if BinWidth is provided correctly - is should be given in ns, and maximum lifetime is 142 ns, BinWidth should not be higher than it or smaller than zero
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Bin width is not reasonable - less or equal zero or greater than 142 ns " << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Bin Width ";
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<3; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'I', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;		//Getting the extended type of data
	if( option != "oscilloscope" && option != "digitizer" && option != "other" )	//Checking if type of times is provided correctly, one of the recognizable type - oscilloscope or digitizer or other
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Type of times is not one of correct types - 'oscilloscope' or 'digitizer' or 'other' " << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Type of Times ";
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<4; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'L', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;		//Getting the value for searching last bin
	if( stoi( option ) < 0 )				//Checking if LastBin minimal value is provided correctly - is should be greater than zero from definition
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "LastBin minimal value is less than zero - is the histogram negative?" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Last Bin Min Value ";
	}
	else if( stoi( option ) > 100 )				//LastBin minimal value corresponds to the level of Background - greater than 100 means that statistics is very high, not usually spotted, this check is for user to confirm that user knows what is doing
	{
		std::cout << "Are You sure You have such high Background (More than 100)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line for proper value of iterations" << std::endl;
			FitDetailsData[0] = FitDetailsData[0] + " --- bad Last Bin Min Value ";
		}
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'F', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;		//Getting the value for searching first bin
	if( stoi( option ) < 0 )				//Similar as for LastBin minimal value
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "FirstBin minimal value is less than zero - is the histogram negative?" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad First Bin Min Value ";
	}
	else if( stoi( option ) > 100 )        //Similar as previous argument
	{
		std::cout << "Are You sure You have such high Background (More than 100)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line for proper value of iterations" << std::endl;
			FitDetailsData[0] = FitDetailsData[0] + " --- bad First Bin Min Value ";
		}
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'B', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;		//Getting the number of points to calculate background to
	if( stoi( option ) < 0 )				//Checking if Number of points for background calculations are given correctly - background can not be estimated from negative number of points
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Background can not be estimated from negative number of points" << std::endl;
			FitDetailsData[0] = FitDetailsData[0] + " --- bad Number of Points in Background ";  
	}
	else if( stoi( option ) > 1000 )			//Checking if user knows what is doing - background follows Poison distribution (in general), and it is correctly estimated if the number of pointsis in the order of 100, more than that is just the exaggeration
	{
		std::cout << "Are You sure You need to calculate background of that many of points(More than 1000)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line for proper value of iterations" << std::endl;
			FitDetailsData[0] = FitDetailsData[0] + " --- bad Number of Points in Background ";
		}
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'B', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;		//Getting the side from which background is estimated
	if( option != "right" && option != "left" )		//Checking if side for background is provided correctly - it should be from one of the sides - left or right
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Side of histogram to estimate background should be one of the following - 'right' or 'left' " << std::endl;
        FitDetailsData[0] = FitDetailsData[0] + " --- bad Side for Background ";
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'L', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;		//Getting the number of points which is leaved from background calculation
	if( stoi( option ) < 0 )				//Checking if Shift for background is given properly - shift can not be lower than zero
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Program can not leave negative number of bins" << std::endl;
        FitDetailsData[0] = FitDetailsData[0] + " --- bad Shift for Background ";
	}
	else if( stoi( option ) > 1000 )			//Checking if user knows what is doing - shifting more than 1000 points is the sign of bad shape of histogram - user must confirm that value is all right
	{
		std::cout << "Are You sure You need to leave that many of points(More than 1000)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line for proper value of iterations" << std::endl;
			FitDetailsData[0] = FitDetailsData[0] + " --- bad Shift for Background ";
		}
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'E', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;		//Getting the End of the Fit value
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'S', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;		//Getting the Start of the Fit value
	if( stod( option ) < 0 || stod( option ) > 1 )		//Checking if StartBin for fit is given properly - it is given by fraction of maximum of histogram -> it can not be smaller than 0 (negative value) or grater than zero (value greater than maximum)
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Program can not start at negative point or higher than maximum" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Start of the Fit ";
	}	
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'F', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;		//Getting the Center for First Bin	
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'L', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;		//Getting the Center for Last Bin	
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'D', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;		//Getting the Option for fixing resolution component	
	if( option != "yes" && option != "no" )		//Checking if Option for fixing is provided correctly - yes or no
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "FixGauss Option can be just one of the following - 'yes' or 'no' " << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Option for Fixing Resolution";
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<4; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'S', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
//-----------------------------------------------
	fitDetails >> option >> option2;
	lineNumber++;
	unsigned MaxIterator = 0;
	while( option != "-------" )
	{
		MaxIterator++;
		lineNumber++;
		if( MaxIterator > 9 )				//If this loop works too long, there is a chance that FitDetails file is corrupted. In other hand, program is not prepared for fitting more than 10 Resolution components
		{
			std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
			std::cout << "Too much Resolution components (>9) or bad style of writing" << std::endl;
			std::cout << "It should be like: [Number] [tabulator between them] [Number]" << std::endl;
			std::cout << "and after writing all components ended by new line with -------	----------" << std::endl;
			FitDetailsData[0] = FitDetailsData[0] + " --- bad Resolution Component nr " + NumberToChar( MaxIterator, 0 );
		}
		if( stod(option) > 0 && stod(option2) > 0 )	//Check if Resolution components are provided correctly
			FitDetailsData.push_back( option + " " + option2 );
		else
		{
			std::cout << "Intensity or Sigma can not be 0!!!!" << std::endl;
			std::cout << "Error while reading component!!!!" << std::endl;
			FitDetailsData[0] = FitDetailsData[0] + " --- bad Resolution Component nr " + NumberToChar( MaxIterator, 0 );
		}
		fitDetails >> option >> option2;
	}
	FitDetailsData[1] = FitDetailsData[1] + NumberToChar( MaxIterator, 0 );
	MaxIterator = 0;
//-----------------------------------------------
	for( unsigned i=0; i<7; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'L', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option >> option2 >> option3;
	lineNumber++;
	while( option != "-------" )
	{
		MaxIterator++;
		lineNumber++;
		if( MaxIterator > 9 )				//Similar as for Resolution components
		{
			std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
			std::cout << "Too much Lifetime components (>9) or bad style of writing" << std::endl;
			std::cout << "It should be like: [Number] [tabulator between them] [Number] [tabulator between them] [Type]" << std::endl;
			std::cout << "and after writing all components ended by new line with -------	----- ----" << std::endl;
			FitDetailsData[0] = FitDetailsData[0] + " --- bad Lifetime Component nr " + NumberToChar( MaxIterator, 0 );
		}
		if( stod(option) > 0 && stod(option2) > 0 )
			FitDetailsData.push_back( option + " " + option2 + " " + option3 );
		else
		{
			std::cout << "Intensity or Lifetime can not be 0!!!!" << std::endl;
			std::cout << "Error while reading component!!!!" << std::endl;
			FitDetailsData[0] = FitDetailsData[0] + " --- bad Lifetime Component nr " + NumberToChar( MaxIterator, 0 );
		}
		fitDetails >> option >> option2 >> option3;
	}
	FitDetailsData[1] = FitDetailsData[1] + " " + NumberToChar( MaxIterator, 0 );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'N', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	if( stoi( option ) < 0 )				//Checking if Number of iteration are given properly, it can not be smaller than zero
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Number of iterations should be greater or equal zero" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Number of Iterations ";
	}
	else if( stoi( option ) > 20 )				//More iterations -> Longer fitting. Sometimes after 20 iterations parameters are not changing more, so it is check if user understands 
	{
		std::cout << "Are You sure You need that much iterations (More than 20)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line for proper value of iterations" << std::endl;
			FitDetailsData[0] = FitDetailsData[0] + " --- bad Number of Iterations ";
		}
	}	  
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'H', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	if( stod( option ) < 0 )				//Variability level should not be smaller than zero - it corrupts the boundaries of the fit parameters
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Level of variance should be greater or equal zero" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Variation Level ";
	}
	else if( stod( option ) > 100 )				//Variability level should not be greater than 1 -> It would provide a chance of fitting negative, unreal parameters
	{
		std::cout << "Are You sure You need that much variability (More than 100%)? (y/n)" << std::endl;
		std::cin >> question;
		if( question == 'n' )
		{
			std::cout << "So try changing line for proper value of variability" << std::endl;
			FitDetailsData[0] = FitDetailsData[0] + " --- bad Variation Level ";
		}
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<5; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'D', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	if( option != "yes" && option != "no" )		//Checking if option for deconvolution is defined one - yes or no
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Deconvolution Option can be just one of the following - 'yes' or 'no' " << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Deconovolution Option ";
	}
	else if( option == "yes" )
		FitDetailsData.push_back( NumberToChar( 1, 0 ) );
	else
		FitDetailsData.push_back( NumberToChar( 0, 0 ) );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'H', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	if( stod( option ) < 0 )				//Denisty parameter characterizes the distances between points on Lifetime Grid used in continous fitting - distance should not be smaller than 0
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Lifetime Grid density parameter should be greater than zero" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Lifetime Grid Density Parameter ";
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'L', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	if( stod( option ) < 0 )				//Minimal point of Lifetime Grid should be greater than zero, because the lifetime of positronium can not be smaller than zero
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Start point should not be negative fraction of Min Lifetime" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Min Lifetime on Lifetime Grid valuer ";
	}
	else if( stod( option ) > 100 )				//Minimal point of Lifetime Grid should be also smaller than minimal value from discrete fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Start point should not be greater than Min Lifetime" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Min Lifetime on Lifetime Grid value ";
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'L', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	if( stod( option ) < 1 )				//Maximal point should not be smaller than maximal value of lifetime fitted from discrete procedure
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "End point should not be lower than Max Lifetime" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Max Lifetime on Lifetime Grid value ";
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'F', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'F', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	if( stod( option ) < 0 )				//Number of point for linear smmothing of histogram, should be greater than zero -> zero if histogram is perfect
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Number of point for linear smoothing should be greater than zero" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Number of Points for Linear Smoothing ";
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'R', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	if( stod( option ) <= 1 )				//Multiplicity of Max Lifeitme for continous fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "End of Range for fitting should not be lower than maximal Lifetime from discrete fitting" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad End of the Fit Range for Deconovolution ";
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'D', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	if( stod( option ) <= 0 )				//Fraction of lifetime set to be default sigmas for lifetime distribution in continous fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Fraction of Lifetime for sigmas should not be less or equal than zero" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Fraction of Lifetime for Smearing of Lifetime component ";
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'S', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	if( stod( option ) <= 0 )				//Step that controls how big steps will be done in preiteration procedure - bigger step quicker fitting, smaller more precise results
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Start value of step can not be lower or equal 0" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Step from the First Part of Deconvolution ";
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'N', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	if( stod( option ) < 0 )				//Number of iteration to check the condition to end preiteration -> higher value -> precise estimation -> longer fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Number of preiterations can not be lower than 0" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Number of Iterations from the First Part of Deconvolution ";
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'S', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	if( stod( option ) <= 0 )				//Step that controls how big steps will be done in iteration procedure - bigger step quicker fitting, smaller more precise results
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Start value of step can not be lower or equal 0" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Step from the Second Part of Deconvolution ";
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'N', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	if( stod( option ) < 0 )				//Number of iteration to check the condition to end iteration -> higher value -> precise estimation -> longer fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Number of iterations can not be lower than 0" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Number of Iterations from the Second Part of Deconvolution ";
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'S', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	if( option != "Gaussian" && option != "LogGaussian" && option != "Mixed" )				//Number of iteration to check the condition to end iteration -> higher value -> precise estimation -> longer fitting
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Shape must be one of the mentioned - Gaussian, LogGaussian or Mixed" << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Shape of the Component ";
	}
	FitDetailsData.push_back( option );
	for( unsigned i=0; i<2; i++ ) //Leaving header
	{
		fitDetails.getline( line, lineWidth );
		lineNumber++;
	}
	if( CheckTheLine( line[0], 'D', lineNumber ) == 0 )
		FitDetailsData[0] = FitDetailsData[0] + " --- bad line nr " + NumberToChar( lineNumber, 0 );
	fitDetails >> option;
	if( option != "yes" && option != "no" )		//Checking if option for deconvolution is defined one - yes or no
	{
		std::cout << "Check if there is everything all right in the FitDetails file" << std::endl;
		std::cout << "Type of Continous Fit Option can be just one of the following - 'yes' or 'no' " << std::endl;
		FitDetailsData[0] = FitDetailsData[0] + " --- bad Type of Continous FIt ";
	}
	else if( option == "yes" )
		FitDetailsData.push_back( NumberToChar( 1, 0 ) );
	else
		FitDetailsData.push_back( NumberToChar( 0, 0 ) );
    
	return FitDetailsData;
}
