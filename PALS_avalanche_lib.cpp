#include "PALS_avalanche_lib.h"

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

Fit::Fit( std::string path, std::string pathForDetails, int ROOTFileTest, std::string FitType )
{							//Getting Times and fitting details from file provided by user and FitDetails file
	TypeOfFit = FitType;
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
	fileTools(pathForDetails);
	if( fileTools.FileCheck( pathForDetails ) )
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
	std::vector<std::string> DataAsString = fileTools.GetFitDetails();
	if( DataAsString[0] != "" ) // Zero element is the error indicator
	{
		std::cout << DataAsString[0] << std::endl;
		return;
	}
	unsigned lineNumber = 2;
	TypeOfData = DataAsString[lineNumber]; // Second element is the multiplicity of Resolution and Lifetime components
	lineNumber++;
	ROOTDirectory = DataAsString[lineNumber];
	lineNumber++;
	ROOTHistogram = DataAsString[lineNumber];
	lineNumber++;
	NameOfTheFileForEXCEL = DataAsString[lineNumber];
	lineNumber++;
	BinWidth = stod( DataAsString[lineNumber] );
	lineNumber++;
	TypeOfDataExtended = DataAsString[lineNumber];
	lineNumber++;
	LastBinMinValue = stoi( DataAsString[lineNumber] );
	lineNumber++;
	FirstBinMinValue = stoi( DataAsString[lineNumber] );
	lineNumber++;
	BackgroundEstimationNumbOfPoints = stoi( DataAsString[lineNumber] );
	lineNumber++;
	SideOfBackgroundEstimation = DataAsString[lineNumber];
	lineNumber++;
	ShiftForBackgroundEstimation = stoi( DataAsString[lineNumber] );
	lineNumber++;
	EndOfFitValue = stod( DataAsString[lineNumber] );
	lineNumber++;
	StartOfFitValue = stod( DataAsString[lineNumber] );
	lineNumber++;
	FirstBinCenter = stod( DataAsString[lineNumber] );// - BinWidth/2;
	lineNumber++;
	LastBinCenter = stod( DataAsString[lineNumber] );// - BinWidth/2;
	lineNumber++;
	FixGauss = DataAsString[lineNumber];
	lineNumber++;
    
	std::vector<int> ComponentsMultiplicites = fileTools.GetComponentsMultiplicities( DataAsString[1] );
	int ComponentsShift = ComponentsMultiplicites[0] + ComponentsMultiplicites[1];
	if( ComponentsMultiplicites[0]*ComponentsMultiplicites[1] == 0 )
	{
		std::cout << "Nor Resolution or Lifetime Components" << std::endl;
		return;
	}
	std::vector<double> ComponentsDetails;
	for( int i=0; i<ComponentsMultiplicites[0]; i++ )
	{
		ComponentsDetails = GetComponentDetails( DataAsString[17+i] ); 
		Resolution.push_back( LifetimeComponent( stod(ComponentsDetails[0]), stod(ComponentsDetails[1]), "gauss" ) );
	}
	for( int i=0; i<ComponentsMultiplicites[1]; i++ )
	{
		ComponentsDetails = GetComponentDetails( DataAsString[17+ComponentsMultiplicites[0]+i] ); 
		Lifetimes.push_back( LifetimeComponent( stod(ComponentsDetails[0]), stod(ComponentsDetails[1]), ComponentsDetails[2]) );
	}
	
	NmbrOfIterations = stoi( DataAsString[lineNumber+ComponentsShift] );
	lineNumber++;
	VarLvl = stod( DataAsString[lineNumber+ComponentsShift] )/100;
	lineNumber++;
	DecoOption = stoi( DataAsString[lineNumber+ComponentsShift] );
	lineNumber++;
	Scaler = stod( DataAsString[lineNumber+ComponentsShift] );
	lineNumber++;
	FracMinLF = stod( DataAsString[lineNumber+ComponentsShift] )/100.;
	lineNumber++;
	FracMaxLF = stod( DataAsString[lineNumber+ComponentsShift] );
	lineNumber++;
	MaxShiftForDeco = stoi( DataAsString[lineNumber+ComponentsShift] );
	lineNumber++;
	LinFilterRange = stoi( DataAsString[lineNumber+ComponentsShift] );
	lineNumber++;
	EndOfFitMultiplicity = stod( DataAsString[lineNumber+ComponentsShift] );
	lineNumber++;
	SigmasDefaultFraction = stod( DataAsString[lineNumber+ComponentsShift] )/100.;
	lineNumber++;
	StepForPreIteration = stod( DataAsString[lineNumber+ComponentsShift] );
	lineNumber++;
	NumberOfPreIterations = stoi( DataAsString[lineNumber+ComponentsShift] );
	lineNumber++;
	StepForIteration = stod( DataAsString[lineNumber+ComponentsShift] );
	lineNumber++;
	NumberOfIterations = stoi( DataAsString[lineNumber+ComponentsShift] );
	lineNumber++;
	ShapeOfComponent = DataAsString[lineNumber+ComponentsShift];
	lineNumber++;
	TypeOfContinousFit = stoi( DataAsString[lineNumber+ComponentsShift] );
	lineNumber++;
	
	int ErrorTest = RenormalizeComponentsIntensities();
	if( ErrorTest )
	{
		Resolution.clear();
		return;
	}
	SortLifetimesByType();
	
	GetHistogramToVector( path, ROOTFileTest );
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
		histogram = FillHistogram( "Time differences", (LastBinCenter - FirstBinCenter)/BinWidth, FirstBinCenter, LastBinCenter, Times );
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
	if( NameOfTheFileForEXCEL == "no" )
		NameOfTheFileForEXCEL = Path;
}

int Fit::Discrete()
{
  	if( Resolution.size()*Lifetimes.size() == 0 )				//Test if the reading FitDetails file and data file past properly
	{
		return 0;
	}
	std::cout<<"-------------Fitting data-------------"<<std::endl;
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
	histogram -> GetYaxis() -> SetTitleOffset(1.3);
	histogram -> GetXaxis() -> SetTitleFont(62);
	histogram -> GetXaxis() -> SetLabelFont(62);
	histogram -> GetYaxis() -> SetTitleFont(62);
	histogram -> GetYaxis() -> SetLabelFont(62);
	histogram -> SetTitle( ("Time differences_ " + Path).c_str() );
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
            //do not had time to study this, as it is not very crucial, because it is fitted to histogram
	std::cout << "Area \t \t \t \t \t \t " << area << std::endl;
	area = fabs(area);
	std::cout << "Range - Arguments [Start] [End] \t \t" << Arguments[ Range_From ] << "   " << Arguments[ Range_To ] << std::endl;
	std::cout << "Range - Values [Start] [End] \t \t \t" << Values[ Range_From ] << "   " << Values[ Range_To ] << std::endl;
	
	int pPsIndex = -1;	
	for( unsigned i=0; i<Lifetimes.size(); i++ )
	{
		if( Lifetimes[i].Type == "ps" )
			pPsIndex = i;
	}

	FitFunction fitTools( TypeOfFit, pPsIndex + 1, Resolution.size(), oPsLFLimit );
	
	fitTools.generateFitFunction( Arguments[ Range_From ], Arguments[ Range_To ], 4 + Lifetimes.size()*2 + 3*Resolution.size() + 1, (Range_To - Range_From)/BinWidth );
	fitTools.generateParameter( 0, Background, "Number of Gauss resolution components", Background - 3*SDBackground, Background + 3*SDBackground, "NoFixing" );
	fitTools.generateParameter( 1, Resolution.size(), "Number of Gauss resolution components", 0., 0., "Fix" );
	fitTools.generateParameter( 2, Lifetimes.size(), "Number of Lifetimes components", 0., 0., "Fix" );
	fitTools.generateParameter( 3, area, "Total area of components", 0., 0., "NoLimits" );
	fitTools.generateResolutionParameter( 4, Resolution, FixGauss, Arguments[ BinMax ] );
	fitTools.generateInitLifetimeParameter( 4 + 3*Resolution.size(), TypeOfFit, Lifetimes, LifetimesNotFixed );
	TF1 *Discrete = FitTools.getFitFunction();
    
	histogram -> Fit(Discrete,"RM");
	std::cout << std::endl;
	std::cout <<"Do not look on intensities of completely free parameters and resolution intensities"<< std::endl;
	std::cout <<"They are normalized in the body of fitting function and at the end when saving results"<< std::endl;
	std::cout << std::endl;
    
	for( unsigned iteration = 0; iteration < NmbrOfIterations; iter++ )
	{
		fitTools.generateIterLifetimeParameter( 4 + 3*Resolution.size(), TypeOfFit, Lifetimes, LifetimesNotFixed, Discrete, iter, VarLvl );
		Discrete = FitTools.getFitFunction();
		
		histogram -> Fit(Discrete,"RM");
		std::cout << std::endl;
		std::cout <<"Do not look on intensities of completely free parameters and resolution intensities"<< std::endl;
		std::cout <<"They are normalized in the body of fitting function and at the end when saving results"<< std::endl;
		std::cout << std::endl;
	}
    
	ResultsSaver Results( TypeOfFit );
	Results.saveDiscreteFit( histogram, Discrete, "Results_from_fit.root", Path, "Results", PathWithDate );
	
	//AreaFromDiscrete = Discrete -> GetParameter( 3 );	If neede in future can be uncommented !?|?!
	Background = Discrete -> GetParameter( 0 );
	Double_t ResolutionsFromFit[22]; //Maximal number of Resolution component is 10 -> 20 parameters + tempfirst + numberofresolutioncomponents
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
			if( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*i ) > oPsLFLimit && pPsIndex >=0 )
				pPsIntensity += Discrete -> GetParameter( 5 + 3*Resolution.size() )*Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
		}
		else if( Lifetimes[i].Type != "ps" )
		{
			AfterFitIterator++;
			if( TypeOfFit == "" )
			{
				FreeIntensities += GetIntensityParameterNew( FreeParameters, 3, 3, AfterFitIterator );
				if( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*i ) > oPsLFLimit && pPsIndex >=0 )
					FreeIntensitiesTopPs += GetIntensityParameterNew( FreeParameters, 3, 3, AfterFitIterator );
			}
			else
			{
				FreeIntensities += Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );
				if( Discrete -> GetParameter( 4 + 3*Resolution.size() + 2*i ) > oPsLFLimit && pPsIndex >=0 )
					FreeIntensitiesTopPs += Discrete -> GetParameter( 5 + 3*Resolution.size() + 2*i );                
			}
		}
	}
	FixedIntensities += pPsIntensity;
	NotFixedIterator = 0;
	if( pPsIndex >=0 )
	{
		pPsIntensity += Discrete -> GetParameter( 5 + 3*Resolution.size() )*FreeIntensitiesTopPs*(1-FixedIntensities)/(1 + FreeIntensitiesTopPs*Discrete -> GetParameter( 5 + 3*Resolution.size() ) );
	}
	
	std::vector< std::vector< DiscreteFitResult > > ResultsDiscrete = Results.saveDiscreteFitResultsToTXTandExcel( Path, "Discrete_Fit_", PathWithDate, 
								Discrete, Background, SDBackground, ResolutionsFromFit, ResolutionsFromFitErrors, 
								pPsIndex, pPsIntensity, FreeParameters, FreeParametersErrors, 
								FixedIntensities, FreeIntensitiesTopPs, FixedFixedIntensity, Lifetimes, 
								oPsLFLimit, pPsLFLimit, NameOfTheFileForEXCEL);
	Results.saveResiduals( histogram, Arguments[ Range_From ], Arguments[ Range_To ], Range_From, Range_To, Discrete, "Residuals.root", Path, PathWithDate );
	Results.saveLFvsIntensities( "LF_vs_In.root", Path, PathWithDate, ResultsDiscrete[0], ResultsDiscrete[1] );
	
	return DecoOption;
}

int Fit::RenormalizeComponentsIntensities()
{
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
		return 1;
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
	return 0;
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

void Fit::GetHistogramToVector( std::string path, int ROOTFileTest )
{
	std::cout << "-------------Getting data-------------" << std::endl;
	std::ifstream data;		//Stream for path to data
	std::string count = "";		//For oscilloscope
	std::string time = "";		//Proper string
	std::string risetime1 = "";	//For oscilloscope
	std::string risetime2 = "";	//For oscilloscope
	
	if( ROOTFileTest ) //Indicator that there is ROOT file given to the program
	{
		TFile* file1 = new TFile( path.c_str(), "READ" );
		TDirectory *dir;
		if( ROOTDirectory != "none" )
			file1->GetObject(ROOTDirectory,dir);
		
		TH1D* histo1 = (TH1D*) dir->Get( ROOTHistogram.c_str() );
		
		double BinWidthOfHisto = histo1->GetBinCenter(2) - histo1->GetBinCenter(1);     //Leaving Underflow bin -> 0
		if( BinWidth < BinWidthOfHisto )
		{
			std::cout << "Cannot get smaller bins from the big once. Rebinning impossible, bad BinWidth in FitDetails or no histogram readed" << std::cout;
			std::cout << "Bin from FitDetails " << BinWidth << " and the Bin of the histogram " << BinWidthOfHisto << std::endl;
			return;
		}
		else
		{
			double RebinNumber = BinWidth/BinWidthOfHisto;
			double decimals = RebinNumber % 1;
			if( decimals > 0.01 )
			{
				std::cout << "BinWidth given by the user in FitDetails is not the multiplicity of BinWidth of the histogram in ROOT file" << std::cout;
				std::cout << "Bin from FitDetails " << BinWidth << " and the Bin of the histogram " << BinWidthOfHisto << std::endl;
				return;               
			}
			std::cout << "Rebinning of the histogram by " << (int)RebinNumber << std::endl;
			histo1 -> Rebin( (int)RebinNumber );
		}
		
		for( int i=0; i<histo1->GetXaxis()->GetNbins(); i++ )
		{
			Times.push_back( histo1->GetBinContent(i) );
			NmbOfBins++;
		}

		//FirstBinCenter = histo1->GetBinCenter(0);
		
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
	else
	{
		data.open( ( Path ).c_str() );

		if( TypeOfData == "histogram" )
		{
			while( data >> time )	//other
			{
				if( time != "inf" && time != "nan")
				{
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
					if( time != "inf" && time != "nan")
					{
						Times.push_back( stod(time) );
					}		
				}
			}
			else			//digitizer
			{
				while( data >> time )
				{
					if( time != "inf" && time != "nan")
					{
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
}
